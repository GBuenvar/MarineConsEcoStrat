using CSV, DataFrames, Random, Statistics, Plots, XLSX, ArgParse
include("percolation_functions.jl")
##

s = ArgParseSettings()
@add_arg_table! s begin
    "--Nqs", "-q"
        help = "Number of q values to plot, from 0 to Nqs-1"
        default = 10.
        arg_type = Float64
    "--alpha", "-a"
        help = "Alpha parameter for the percolation"
        default = 1.
        arg_type = Float64
    "--saving", "-s"
        help = "Save the outputs"
        default = true
        arg_type = Bool

end

p = parse_args(ARGS, s)

N_qs = p["Nqs"]
α = p["alpha"]
saving = p["saving"]
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
# read the economic data
economic_data = DataFrame(XLSX.readtable("data/CLASS.xlsx", "List of economies"))

mkpath("percolation/Diversity")

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

# Initial conditions
rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
EEZ_neis_dict = compute_neighbors(agg_data)
rich = rich[1:end-2] # remove the high seas and the antarctic
q = 0.



## 

# Plot initialization


l = @layout [ [a{0.49w} b{0.49w}
               c{0.49w} d{0.49w}] e{0.05w} ]
reds = cgrad([RGBA(1,0,0,1), RGBA(1,0,0,0.2)], Int(N_qs), categorical = true, rev=true)
blacks  = cgrad([RGBA(0,0,0,1), RGBA(0,0,0,0.2)], Int(N_qs), categorical = true, rev=true)
blues = cgrad([RGBA(0,0,1,1.), RGBA(0,0,1,0.2)], Int(N_qs), categorical = true, rev=true)


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
plot!(cost_ind_plot, 
        ylabel="Cost/individuals", 
        yguidefontcolor=:black)

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
plot!(eez_div_plot, 
    ylabel="EEZ diversity (H(q))",
    yguidefontcolor= :blue)


# iterate over the q values and plot the results
for (i,q) in enumerate(0.:N_qs-1)
    println("________________")
    println("q = $q")
    println("________________")
    protected_times, protected_number, prot_cost, prot_eezs, global_diversity, eez_div = run_diversity(agg_data; start_protecting = rich, α=α, q=q)
    prot_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

    # save the data
    if saving        
        mkpath("percolation/Diversity/q$(Int(q))")
        CSV.write("percolation/Diversity/q$(Int(q))/protected_times.csv.gz", DataFrame(protected_times=protected_times))
        CSV.write("percolation/Diversity/q$(Int(q))/protected_number.csv.gz", DataFrame(protected_number=protected_number))
        CSV.write("percolation/Diversity/q$(Int(q))/protected_species_times.csv.gz", DataFrame(prot_species_times=protected_species_times))
        CSV.write("percolation/Diversity/q$(Int(q))/protected_species_number.csv.gz", DataFrame(prot_species_number=prot_species_number))
    end 



    # Add the plots
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


# Add legends manually to the plots
plot!(pqs, [[NaN], [NaN]],
    c = [:red :black],
    label=["Species" "Individuals"],)
plot!(cost_plot, [[NaN], [NaN]],
    c = [:blue :black],
    label=["EEZ" "Cumulated"],)
plot!(cost_species_plot, [[NaN], [NaN]],
    c = [:red :black],
    label=["Species" "Individuals"],)


# Create the custom colorbar
xx = range(0,1,100)
ZZ = Matrix{Float64}(undef, 100, 3)
ZZ[:,:] .= NaN

cbar = heatmap(ZZ,
    xticks=false, 
    yticks=false,#((collect(range(0, N_qs-1, 6)) .+ 0.5) .*100 ./ N_qs, range(0, N_qs-1, 6)), 
    size = (50,300), 
    legend=false, 
    fc=reds, 
    # lims=(0,1),
    framestyle=:none, 
    ylabel="q",
    yguidefontrotation=-90,
    );
# cbar
for ytck in range(0, N_qs-1, 6)
    annotate!(cbar,
                -0.3,
                (ytck+0.5).*100 ./ N_qs,
                (ytck, 8))
end

for (ic, col) in enumerate([reds, blues, blacks])
    zz = Matrix{Float64}(undef, 100, 3)
    zz[:, :] .= NaN
    zz[:, ic] .= xx 
    heatmap!(cbar,
                zz,
                # xticks=false, 
                # yticks=(0.5/(N_qs-1):0.2:1),#, range(0, Int(N_qs)-1, 6)), 
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

if saving
    savefig(plot_q, "percolation/figures/diversity_q_$(Int(N_qs+1)).png")
    savefig(plot_q, "percolation/figures/diversity_q_$(Int(N_qs+1)).pdf")
end
