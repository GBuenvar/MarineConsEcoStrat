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

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

# Read the economic data
economic_data = DataFrame(XLSX.readtable("data/CLASS.xlsx", "List of economies"))
eez_to_int = Dict(zip(eez_codes.EEZ, eez_codes.Int))

# read the eez_to_iso3 file
eez_to_iso3_data = CSV.read("data/eez_to_iso3.csv", DataFrame)
eez_to_iso3 = Dict(zip(eez_to_iso3_data.Country, eez_to_iso3_data.ISO_3digit))
# add High Seas
eez_to_iso3["-1"] = "-1"

mkpath("percolation/Game")

##
# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs]
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

function run_diversity(data; start_protecting = [-1], α=1., q=1)
    ids = unique(data[:, :newid])
    eezs = unique(data[:, :EEZ])
    iterated_eezs = setdiff(eezs, start_protecting)
    Neez = length(iterated_eezs)

    data_div = copy(data) # make a copy of the data to compute the diversity of the protected areas


    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, Neez+1)
    prot_cost   = zeros(Float64, Neez+1)
    prot_eezs   = copy(start_protecting)
    global_diversity = zeros(Int64, Neez+1)

    # Compute the initial diversity of the protected areas
    unique_pairs = unique(data_div[:, ["newid", "Species", "EEZ"]])
    sizes = Vector{Float64}(unique_pairs[unique_pairs.EEZ .∈ (start_protecting, ), :Species])
    sizes = sizes ./ sum(sizes)
    initial_div = hill_number(sizes, q=q)
    global_diversity[1] = initial_div
    
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

    # Compute the potential cost of each EEZ being the last to be protected
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    for ii in 2:Neez+1
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])



        # 1- compute the diversity of each EEZ
        # println("computing diversity for q = ", q, " at time ", ii, " with ", length(unprotected_eezs), " EEZs left.")
        eezs_diversity = zeros(Float64, length(unprotected_eezs))
        for (i, eez) in enumerate(unprotected_eezs)
            sizes = Vector{Float64}(unique_pairs[unique_pairs.EEZ .== eez, :Species])
            sizes = sizes ./ sum(sizes)
            eezs_diversity[i] = hill_number(sizes, q=q)
        end
        # 2- find the EEZ with the highest diversity, and protect it. If there is a tie, protect the one with the highest diversity with q+1
        max_diversity = maximum(eezs_diversity)
        arg = findall(eezs_diversity .== max_diversity)
        to_protect = unprotected_eezs[arg]
        while length(to_protect) > 1
            # if the distribution is exactly the same for all the eez, just pick one randomly
            sizes = [Vector{Float64}(unique_pairs[unique_pairs.EEZ .== eez, :Species]) for eez in to_protect]
            sizes = [s./ sum(s) for s in sizes]
            if all([all(sizes[1] .== s) for s in sizes[2:end]])
                to_protect = [to_protect[rand(1:end)]]
                break
            end
            println("There is a tie, breaking it with q=$(q+1)")
            q += 1
            eezs_diversity = zeros(Float64, length(to_protect))
            for (i, size) in enumerate(sizes)
                eezs_diversity[i] = hill_number(size, q=q)
            end
            max_diversity = maximum(eezs_diversity)
            arg = findall(eezs_diversity .== max_diversity)
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
        # println("EEZ: ", protect_eez, " time: ", ii, " # protected: ", length(new_protected_ids), " # unprotected: ", length(new_unprotected_ids), " cost: ", cost)
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        # println(length(unique(data[:, :EEZ])))
    end
    return prot_times, prot_number, prot_cost, prot_eezs
end

# %%


q = 0.
α = 0.
##
l = @layout [ [a{0.49w} b{0.49w}
               c{0.49w} d{0.49w}] e{0.1w} ]
# get a list of colors for the plots using a gradient from red to blue
Nqs = 11.
reds = cgrad([:red, :white], Int(Nqs), categorical = true, rev=true)
blacks  = cgrad([:black, :white], Int(Nqs), categorical = true, rev=true)
blues = cgrad([:blue, :white], Int(Nqs), categorical = true, rev=true)
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
    xticks=(1:20:(length(eezs)- length(rich)),length(rich):20:length(eezs))
    )         
for (i,q) in enumerate(0.:Nqs-1)
    protected_times, protected_number, prot_cost, prot_eezs = run_diversity(agg_data; start_protecting = rich, α=α, q=q)
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
    plot!(cost_species_plot, cumsum(prot_cost)./cumsum(protected_number),
        label=false, 
        color = blacks[i], 
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
zz = zero(xx)' .+ xx
cbar = heatmap(xx, xx, zz, 
                xticks=false, 
                yticks=(0.5/(Nqs-1):0.2:1),#, range(0, Int(Nqs)-1, 6)), 
                ratio=20, 
                legend=false, 
                fc=reds, lims=(0,1),
                framestyle=:box, 
                ylabel="q",
                yguidefontrotation=-90,
                );
# plot!(cbar, title="q")

plot_q = plot!(pqs, cost_plot, cost_species_plot, cost_species_plot, cbar, layout=l, size=(1000, 800), dpi=300, bottom_margin=20Plots.PlotMeasures.px, left_margin=20Plots.PlotMeasures.px)



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

