using CSV, DataFrames, Random, Statistics, Plots, XLSX
include("percolation_functions.jl")
##

# En este script se ejecuta una percolación en la que se parte con los países ricos (higer income) cooperando. Se deja, como en otros casos, el parametro alpha
# libre. Cuando un animal visita vairos paises, su peso, 1, se divide uniformemente 
# entre todos los paises que visita. Esto significa que se le asigna a cada país la misma 
# probabilidad de llevarse el total del beneficio por cazarlo.
# cuando se agrega por especies, esta simetría se rompe, ya que aunque haya patrones caracteristicos,
# la trayectoria de un individuo es unica y puede visitar distintos paises.

# Podemos hacer una percolación basada en diversidad, en la que se priorice la proteccion de los
# mas diversos. Para eso utilizamos el indice de Hill. POdemos pintar la curva de conservacion para 
# diferentes q. Con q=0, es el numero de especies, independientemente de su abundancia. Con q=1, es el
# numero efectivo de especies de Shannon. Con q=2, es el numero efectivo de especies de Simpson. 
# A mayor q, mas peso se le da a las especies mas abundantes. 

# ^qD = (\sum p_i^q)^{1/(1-q)}; q =/= 1
# ^1D = exp(-\sum p_i log(p_i)); q = 1
##
# open the full_data_inds.csv.gz file

# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))
eez_to_int = Dict(zip(eez_codes.EEZ, eez_codes.Int))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

# read the eez_to_iso3 file
eez_to_iso3_data = CSV.read("data/eez_to_iso3.csv", DataFrame)
eez_to_iso3 = Dict(zip(eez_to_iso3_data.Country, eez_to_iso3_data.ISO_3digit))
# add High Seas
eez_to_iso3["-1"] = "-1"

mkpath("percolation/Game")

##
# Read the data, consisting on the newid, the species, the eez and the time spent in the eez    
agg_data = CSV.read("data/agg_data.csv", DataFrame)

id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs];
##





##
# Write three functions to make the plots p1, pcost and p3

function plot_protected(prot_number, prot_species, n, n_species, title; xlabel = "EEZs cooperating", ylabel = "Fraction protected", savename = "none", title_location = :center)
    p1 = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(p1, cumsum(prot_number)./ n, label="individuals", color = "black")
    plot!(p1, cumsum(prot_species)./n_species, label="species (50%)", color = "red")
    title!(p1, title, title_location = title_location)
    if savename != "none"
        savefig(p1, "percolation/figures/$(savename).pdf")
        savefig(p1, "percolation/figures/$(savename).png")
    end
    plot!(p1)
    return p1
end

function plot_protection_cost(prot_cost, title; xlabel = "EEZs cooperating", ylabel = "Spents", savename = "none", title_location = :center)
    pcost = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(pcost, prot_cost, label="Per EEZ", color = "black")
    plot!(pcost, cumsum(prot_cost), label="Cumulative", color = "blue")
    title!(pcost, title, title_location = title_location)
    if savename != "none"
        savefig(pcost, "percolation/figures/$(savename).pdf")
        savefig(pcost, "percolation/figures/$(savename).png")
    end
    plot!(pcost)
    return pcost
end

function plot_protection_cost_per_individual(prot_cost, prot_number, prot_species, title; xlabel = "EEZs cooperating", ylabel = "Cost/individuals", savename = "none", title_location = :center)
    p3 = plot(xlabel = xlabel, ylabel = ylabel)
    p3twinx = twinx(p3)
    plot!(p3twinx, ylabel = "Cost/species")
    plot!(p3, cumsum(prot_cost[2:end])./(cumsum(prot_number[2:end])), label="individuals", color = "black")
    plot!(p3twinx, cumsum(prot_cost[2:end])./(cumsum(prot_species[2:end])), label=false, color = "red")
    plot!(p3, [-1], [0], label="species", color = "red")
    title!(p3, title, title_location = title_location)
    # plot!(p3, margin = 5) 
    if savename != "none"
        savefig(p3, "percolation/figures/$(savename).pdf")
        savefig(p3, "percolation/figures/$(savename).png")
    end
    plot!(p3)
    return p3
end


function plot_number_ids(data, prot_eez; xlabel = "EEZs cooperating", ylabel = "Number of individuals", label = "Individuals", color = "black", title="", title_location = ":center")
    ids_at_eez = [sum(data[:, :EEZ] .== eez) for eez in prot_eez]
    p4 = plot(ids_at_eez, xlabel = xlabel, ylabel = ylabel, label = label, color = color, title=title, title_location = title_location)
    return p4
end
###

##

operation = /

##

rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
EEZ_neis_dict = compute_neighbors(agg_data)


## 

# Diversity functions:
function hill_number(p::Vector{Float64}; q::Float64=0., normalization_threshold::Float64=10^-6)
    @assert abs(sum(p) - 1 < normalization_threshold) "probabilities are not normalized"
    @assert isinteger(q) "q must be a float type integer"
    p = p[p .> 0.]
    q == 1. ? (h = exp(-sum(p .* log.(p)))) : (h = sum(p.^q)^(1/(1-q)))
    return h
end

function compute_species_dist(df)
    # df is a dataframe with the columns newid, Species, EEZ
    # restricted to the eez of interest
    id_sp_pairs = unique(df[:, [:newid, :Species]])
    sp_size = combine(groupby(id_sp_pairs, :Species), nrow)
    sizes = Vector{Float64}(sp_size.nrow)
    sizes = sizes ./ sum(sizes)
    # sort the sizes
    sizes = sort(sizes, rev=true)
    return sizes
end


function compute_div(df; q=0.)
    # df is a dataframe with the columns newid, Species, EEZ
    # restricted to the eez of interest
    sizes = compute_species_dist(df)
    df_div = hill_number(sizes, q=q)
    return df_div
end

## 

function run_diversity(data; start_protecting = [-1], α=1., q=1.)

    # 
    ids = unique(data[:, :newid])
    eezs = unique(data[:, :EEZ])
    iterated_eezs = setdiff(eezs, start_protecting)
    Neez = length(iterated_eezs)
    data_div = data # make a copy of the data to compute the diversity of the protected areas


    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, Neez+1)
    prot_cost   = zeros(Float64, Neez+1)
    prot_eezs   = copy(start_protecting)
    global_div = zeros(Float64, Neez+1)
    eez_div    = zeros(Float64, Neez+1)

    
    # Protect the first EEZs
    data = data[data[:, :EEZ] .∈ (iterated_eezs, ), :] # protect the EEZ
    new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
    new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
    if length(new_protected_ids) > 0
        # add the new protected individuals to the list
        prot_number[1] = length(new_protected_ids)
        prot_times[ids .∈ (new_protected_ids,)] .= 0 # add the time at which they were protected
    end

    unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    unprotected_eezs = copy(iterated_eezs)
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])

    # Compute the initial diversity of the protected areas
    initial_global_div = compute_div(data_div[data_div.EEZ .∈ [start_protecting], :], q=q)
    mean_initial_div = mean([compute_div(data_div[data_div.EEZ .== eez, :], q=q) for eez in start_protecting])
    global_div[1] = initial_global_div
    eez_div[1] = mean_initial_div
    # 1- compute the diversity of each unprotected EEZ
    eezs_diversity = zeros(Float64, length(unprotected_eezs))
    eez_sizes = Vector{Vector{Float64}}(undef, length(unprotected_eezs))
    # Initialize a vector of vector{Float64} that potentially will be of different lengths    
    for (i, eez) in enumerate(unprotected_eezs)
        s = compute_species_dist(unique_pairs[unique_pairs.EEZ .== eez, :])
        eez_sizes[i] = s 
        eezs_diversity[i] = hill_number(s, q=q)
    end
    
    for ii in 2:Neez+1
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])

        # 2- find the EEZ with the highest diversity, and protect it. If there is a tie, protect the one with the highest diversity with q+1
        max_diversity = maximum(eezs_diversity)
        arg = findall(eezs_diversity .== max_diversity)
        to_protect = unprotected_eezs[arg]

        q0 = q
        while length(to_protect) > 1
            # if the distribution is exactly the same for all the eez, just pick one randomly
            ss = eez_sizes[arg] 
            if all([ss[1] == s for s in ss[2:end]])
                to_protect = [to_protect[rand(1:end)]]
                break
            end
            q0 += 1
            println("There is a tie, breaking it with q=$(q0)")
            eezs_diversity_tiebraker = zeros(Float64, length(to_protect))
            for (i, size) in enumerate(ss)
                eezs_diversity_tiebraker[i] = hill_number(size, q=q)
            end
            max_diversity = maximum(eezs_diversity_tiebraker)
            arg = findall(eezs_diversity_tiebraker .== max_diversity)
            to_protect = to_protect[arg]
        end
        protect_eez = to_protect[1]

        # 3- compute the id weight only for the individuals that are present in protect_eez
        visiting_ids = unique(unique_pairs[unique_pairs.EEZ .== protect_eez, :newid])
        eez_id_weights = compute_id_weight(unique_pairs[unique_pairs.newid .∈ (visiting_ids, ), :], α=α)

        # 4- compute the cooperation cost of the EEZ
        cost = sum(eez_id_weights)

        # 5- Protect that EEZ
        push!(prot_eezs, protect_eez)
        prot_cost[ii] = cost
        data = data[data[:, :EEZ] .!= protect_eez, :] # remove protect_eez from the data
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[ii] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= ii-1 # add the time at which they were protected
        end

        # drop the diversity of the protected eez
        eez_div[ii] = eezs_diversity[unprotected_eezs .== protect_eez][1]
        eezs_diversity = eezs_diversity[unprotected_eezs .!= protect_eez]
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        # println(length(unique(data[:, :EEZ])))

        # 6- Compute the diversity of the protected areas
        global_div[ii] = compute_div(data_div[.!(data_div.EEZ .∈ [unprotected_ids, ]), :], q=q)
    end
    return prot_times, prot_number, prot_cost, prot_eezs, global_div, eez_div
end

# %%


q = 0.
α = 0.
##
l = @layout [ [a{0.49w} b{0.49w}
               c{0.49w} d{0.49w}] e{0.05w} ]
# get a list of colors for the plots using a gradient from red to blue
rich = rich[1:end-2]
Nqs = 6.
reds = cgrad([RGBA(1,0,0,1), RGBA(1,0,0,0.2)], Int(Nqs), categorical = true, rev=true)
blacks  = cgrad([RGBA(0,0,0,1), RGBA(0,0,0,0.2)], Int(Nqs), categorical = true, rev=true)
blues = cgrad([RGBA(0,0,1,1.), RGBA(0,0,1,0.2)], Int(Nqs), categorical = true, rev=true)
pqs = plot( label = false,
            dpi=300, 
            size=(400, 400), 
            xlabel="EEZs cooperating", 
            ylabel="Fraction protected",
            xticks=(1:20:(length(eezs)- length(rich)),length(rich):20:length(eezs))
            )
cost_plot = plot(label=false, 
                dpi=300, 
                size=(400, 400), 
                xlabel="EEZs cooperating", 
                ylabel="Cost",
                xticks=(1:20:(length(eezs)- length(rich)),length(rich):20:length(eezs))
                ) 
        
cost_species_plot = plot(
    label=false, 
    dpi=300, 
    size=(400, 400), 
    xlabel="EEZs cooperating", 
    ylabel="Cost/especies", 
    yguidefontcolor=:red,
    xticks=(1:20:(length(eezs)- length(rich)),length(rich):20:length(eezs))
    )         

cost_ind_plot = twinx(cost_species_plot)
plot!(cost_ind_plot, ylabel="Cost/individuals", yguidefontcolor=:black)

div_plot = plot(
    label=false, 
    dpi=300, 
    size=(400, 400), 
    xlabel="EEZs cooperating", 
    ylabel="Protected areas diversity (H(q))",
    yguidefontcolor= :black,
    xticks=(1:20:(length(eezs)- length(rich)),length(rich):20:length(eezs))
)

eez_div_plot = twinx(div_plot)
plot!(eez_div_plot, ylabel="EEZ diversity (H(q))",
    yguidefontcolor= :blue)

for (i,q) in enumerate(0.:Nqs-1)
    println("________________")
    println("q = $q")
    println("________________")
    protected_times, protected_number, prot_cost, prot_eezs, global_diversity, eez_div = run_diversity(agg_data; start_protecting = rich, α=α, q=q)
    prot_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)


    plot!(pqs, cumsum(protected_number)./N, 
        label=false, 
        color = blacks[i], 
        lw=1
        )

    plot!(pqs, cumsum(prot_species_number)./N_species,
        label=false, 
        color = reds[i], 
        lw=1
        )

    plot!(cost_plot, cumsum(prot_cost),
        label=false, 
        color = blacks[i], 
        lw=1,
        lty=:dot
        )

    plot!(cost_plot, prot_cost,
        label=false, 
        color = blues[i], 
        lw=1
        )
    
    plot!(cost_species_plot, cumsum(prot_cost)./cumsum(prot_species_number),
        label=false, 
        color = reds[i], 
        lw=1,
        lty=:dot
        )
    plot!(cost_ind_plot, cumsum(prot_cost)./cumsum(protected_number),
        label=false, 
        color = blacks[i], 
        lw=1,
        lty=:dot
        )

    plot!(div_plot, global_diversity,
        label=false, 
        color = blacks[i], 
        lw=1,
        lty=:dot
        )

    plot!(eez_div_plot, eez_div,
        label=false, 
        color = blues[i], 
        lw=1,
        lty=:dot
        )


end

plot!(pqs, [[NaN], [NaN]],
    c = [:red :black],
    label=["Species" "Individuals"],)
plot!(cost_plot, [[NaN], [NaN]],
    c = [:blue :black],
    label=["EEZ" "Cumulated"],)
plot!(cost_species_plot, [[NaN], [NaN]],
    c = [:red :black],
    label=["Species" "Individuals"],)


# set legend title
xx = range(0,1,100)
zz = Matrix{Float64}(undef, 100, 3)
zz[:,:] .= NaN
# zz[:, 1] .= xx

cbar = heatmap(zz,
    xticks=false, 
    yticks=false,#((collect(range(0, Nqs-1, 6)) .+ 0.5) .*100 ./ Nqs, range(0, Nqs-1, 6)), 
    size = (50,300), 
    legend=false, 
    fc=reds, 
    # lims=(0,1),
    framestyle=:none, 
    ylabel="q",
    yguidefontrotation=-90,
    );
# cbar
for ytck in range(0, Nqs-1, 6)
    annotate!(cbar,
                -0.3,
                (ytck+0.5).*100 ./ Nqs,
                (ytck, 8))
end

for (ic, col) in enumerate([reds, blues, blacks])
    zz = Matrix{Float64}(undef, 100, 3)
    zz[:, :] .= NaN
    zz[:, ic] .= xx 
    heatmap!(cbar,
                zz,
                # xticks=false, 
                # yticks=(0.5/(Nqs-1):0.2:1),#, range(0, Int(Nqs)-1, 6)), 
                # ratio=20, 
                legend=false, 
                fc=col, 
                # lims=(0,1),
                framestyle=:box, 
                # ylabel="q",
                # yguidefontrotation=-90,
                );
end
cbar

plot_q = plot!(pqs, cost_plot, cost_species_plot, div_plot, cbar, layout=l, size=(1000, 800), dpi=300, bottom_margin=20Plots.PlotMeasures.px, left_margin=20Plots.PlotMeasures.px)



# plot(pqs)

# %%


q = 0.
α = 0.

protected_times, protected_number, prot_cost, prot_eezs = run_diversity(agg_data; start_protecting = rich, α=α, q=q)
prot_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)


p1 = plot_protected(protected_number, prot_species_number, N, N_species, "Protected individuals and species", savename = "none", title_location = :left)
pcost = plot_protection_cost(prot_cost, "Protection cost", savename = "none", title_location = :left)
p3 = plot_protection_cost_per_individual(prot_cost, protected_number, prot_species_number, "Protection cost per individual and species", title_location = :left)
p4 = plot_number_ids(agg_data, prot_eezs, title="Number of individuals per EEZ", title_location = :left)
q+=1
plot(p1, pcost, p3, p4, layout = (2, 2), size = (1000, 1000), )

