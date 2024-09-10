using CSV, DataFrames, Random, Statistics, Plots, XLSX, ArgParse
include("percolation_functions.jl")
##

# En este caso miramos cuales son los individuos que son más fáciles de proteger, es decir,
# los que visitan menos EEZs, P. De las EEZs visitadas por los individuos P que sean más fáciles de proteger,
# se protegen primero las que más individuos P visitan. Cada vez que se protege una EEZ, se recalculan los individuos
# que son más fáciles de proteger, y se repite el proceso. El output es el numero de individuos protegidos en cada paso
# y el numero de EEZ protegidas antes de que cada individuo se proteja.
s = ArgParseSettings()
@add_arg_table! s begin
    "--rich_protect", "-r"
        help = "Start protecting the rich countries"
        default = true
        arg_type = Bool
end

p = parse_args(ARGS, s)
rich_protect = p["rich_protect"]
suffix = rich_protect ? "_rich" : ""


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

mkpath("percolation/RankDegreeAlpha")

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
# Write three functions to make the plots p1, pcost and p3

function plot_protected(prot_number, prot_species, n, n_species, title; xlabel = "EEZs protecting", ylabel = "Fraction protected", savename = "none")
    p1 = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(p1, cumsum(prot_number)./ n, label="individuals", color = "black")
    plot!(p1, cumsum(prot_species)./n_species, label="species (50% of individuals)", color = "red")
    title!(p1, title)
    if savename != "none"
        savefig(p1, "percolation/figures/$(savename).pdf")
        savefig(p1, "percolation/figures/$(savename).png")
    end
    plot!(p1)
    return p1
end

function plot_protection_cost(prot_cost, title; xlabel = "EEZs protecting", ylabel = "Spents", savename = "none")
    pcost = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(pcost, prot_cost, label="Per EEZ", color = "black")
    plot!(pcost, cumsum(prot_cost), label="Cumulative", color = "blue")
    title!(pcost, title)
    if savename != "none"
        savefig(pcost, "percolation/figures/$(savename).pdf")
        savefig(pcost, "percolation/figures/$(savename).png")
    end
    plot!(pcost)
    return pcost
end

function plot_protection_cost_per_individual(prot_cost, prot_number, prot_species, n, n_species, title; xlabel = "# EEZs protecting", ylabel = "Cost/individuals", savename = "none")
    p3 = plot(xlabel = xlabel, ylabel = ylabel)
    p3twinx = twinx(p3)
    plot!(p3twinx, ylabel = "Cost/species")
    plot!(p3, cumsum(prot_cost[2:end])./(cumsum(prot_number[2:end])), label="individuals", color = "black")
    plot!(p3twinx, cumsum(prot_cost[2:end])./(cumsum(prot_species[2:end])), label=false, color = "red")
    plot!(p3, [-1], [0], label="species", color = "red")
    title!(p3, title)
    # plot!(p3, margin = 5) 
    if savename != "none"
        savefig(p3, "percolation/figures/$(savename).pdf")
        savefig(p3, "percolation/figures/$(savename).png")
    end
    plot!(p3)
    return p3
end

## 

rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
rich = rich_protect ? rich : [0,8]

##

α = 1
print("Higher payoff first (a=$α)")
protected_times, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="higher")
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species are protected with $(protected_number[1]) individuals") 


p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Higher payoff first (a=$α)", savename = "Higher_protected_a_$(α)$suffix")
pcost = plot_protection_cost(protected_cost, "Higher payoff first (a=$α)", savename = "Higher_cost_a_$(α)$suffix")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Higher payoff first (a=$α)", savename = "Higher_protected_cost_a_$(α)$suffix")

# show the plots
plot(p1, pcost, p3, layout = (1, 3), size = (1200, 400))

# save the results
# save protected times and number, and species times and number
CSV.write("percolation/RankDegreeAlpha/protected_times_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/RankDegreeAlpha/protected_number_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/RankDegreeAlpha/protected_species_times_desc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_times=protected_species_times))
CSV.write("percolation/RankDegreeAlpha/protected_species_number_desc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_number=protected_species_number))
CSV.write("percolation/RankDegreeAlpha/protected_eezs_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_eezs=protected_eezs))
CSV.write("percolation/RankDegreeAlpha/protected_cost_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_cost=protected_cost))
println("files saved at percolation/RankDegreeAlpha")


##

α = 0
print("Higher payoff first (a=$α)")
protected_times, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="higher")
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species are protected with $(protected_number[1]) individuals")

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Higher payoff first (a=$α)", savename = "Higher_protected_a_$(α)$suffix")
pcost = plot_protection_cost(protected_cost, "Higher payoff first (a=$α)", savename = "Higher_cost_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Higher payoff first (a=$α)", savename = "Higher_protected_cost_a_$(α)$suffix")

# show the plots
plot(p1, pcost, p3, layout = (1, 3), size = (1200, 400))


# save the results
# save protected times and number, and species times and number
CSV.write("percolation/RankDegreeAlpha/protected_times_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/RankDegreeAlpha/protected_number_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/RankDegreeAlpha/protected_species_times_desc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_times=protected_species_times))
CSV.write("percolation/RankDegreeAlpha/protected_species_number_desc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_number=protected_species_number))
CSV.write("percolation/RankDegreeAlpha/protected_eezs_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_eezs=protected_eezs))
CSV.write("percolation/RankDegreeAlpha/protected_cost_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_cost=protected_cost))
println("files saved at percolation/RankDegreeAlpha")


##

# now do the same for α=-1

α = -1
print("Higher payoff first (a=$α)")
protected_times, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="higher")
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species are protected with $(protected_number[1]) individuals") 

##

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Higher payoff first (a=$α)", savename = "Higher_protected_a_$(α)$suffix")
pcost = plot_protection_cost(protected_cost, "Higher payoff first (a=$α)", savename = "Higher_cost_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Higher payoff first (a=$α)", savename = "Higher_protected_cost_a_$(α)$suffix")

# show the plots
plot(p1, pcost, p3, layout = (1, 3), size = (1200, 400))
# save the results
# save protected times and number, and species times and number
CSV.write("percolation/RankDegreeAlpha/protected_times_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/RankDegreeAlpha/protected_number_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/RankDegreeAlpha/protected_species_times_desc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_times=protected_species_times))
CSV.write("percolation/RankDegreeAlpha/protected_species_number_desc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_number=protected_species_number))
CSV.write("percolation/RankDegreeAlpha/protected_eezs_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_eezs=protected_eezs))
CSV.write("percolation/RankDegreeAlpha/protected_cost_desc_alpha_$α$suffix.csv.gz", DataFrame(protected_cost=protected_cost))
println("files saved at percolation/RankDegreeAlpha")





# repeat for order = "lower"
##
α = 1
print("Lower payoff first (a=$α)")
protected_times, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="lower")
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species are protected with $(protected_number[1]) individuals") 

##
p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Lower payoff first (a=$α)", savename = "Lower_protected_a_$(α)$suffix")
pcost = plot_protection_cost(protected_cost, "Lower payoff first (a=$α)", savename = "Lower_cost_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Lower payoff first (a=$α)", savename = "Lower_protected_cost_a_$(α)$suffix")

# show the plots
plot(p1, pcost, p3, layout = (1, 3), size = (1200, 400))

# save the results
# save protected times and number, and species times and number
CSV.write("percolation/RankDegreeAlpha/protected_times_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/RankDegreeAlpha/protected_number_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/RankDegreeAlpha/protected_species_times_asc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_times=protected_species_times))
CSV.write("percolation/RankDegreeAlpha/protected_species_number_asc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_number=protected_species_number))
CSV.write("percolation/RankDegreeAlpha/protected_eezs_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_eezs=protected_eezs))
CSV.write("percolation/RankDegreeAlpha/protected_cost_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_cost=protected_cost))
println("files saved at percolation/RankDegreeAlpha")

##
α =0
print("Lower payoff first (a=$α)")
protected_times, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="lower")
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species are protected with $(protected_number[1]) individuals")

##
p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Lower payoff first (a=$α)", savename = "Lower_protected_a_$(α)$suffix")
pcost = plot_protection_cost(protected_cost, "Lower payoff first (a=$α)", savename = "Lower_cost_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Lower payoff first (a=$α)", savename = "Lower_protected_cost_a_$(α)$suffix")

# show the plots
plot(p1, pcost, p3, layout = (1, 3), size = (1200, 400))
# save the results
# save protected times and number, and species times and number
CSV.write("percolation/RankDegreeAlpha/protected_times_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/RankDegreeAlpha/protected_number_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/RankDegreeAlpha/protected_species_times_asc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_times=protected_species_times))
CSV.write("percolation/RankDegreeAlpha/protected_species_number_asc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_number=protected_species_number))
CSV.write("percolation/RankDegreeAlpha/protected_eezs_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_eezs=protected_eezs))
CSV.write("percolation/RankDegreeAlpha/protected_cost_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_cost=protected_cost))
println("files saved at percolation/RankDegreeAlpha")

##
α = -1
print("Lower payoff first (a=$α)")
protected_times, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="lower")
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species are protected with $(protected_number[1]) individuals") 
##

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Lower payoff first (a=$α)", savename = "Lower_protected_a_$(α)$suffix")
pcost = plot_protection_cost(protected_cost, "Lower payoff first (a=$α)", savename = "Lower_cost_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Lower payoff first (a=$α)", savename = "Lower_protected_cost_a_$(α)$suffix")

# show the plots
plot(p1, pcost, p3, layout = (1, 3), size = (1200, 400))
# save the results
# save protected times and number, and species times and number
CSV.write("percolation/RankDegreeAlpha/protected_times_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/RankDegreeAlpha/protected_number_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/RankDegreeAlpha/protected_species_times_asc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_times=protected_species_times))
CSV.write("percolation/RankDegreeAlpha/protected_species_number_asc_alpha_$α$suffix.csv.gz", DataFrame(prot_species_number=protected_species_number))
CSV.write("percolation/RankDegreeAlpha/protected_eezs_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_eezs=protected_eezs))
CSV.write("percolation/RankDegreeAlpha/protected_cost_asc_alpha_$α$suffix.csv.gz", DataFrame(protected_cost=protected_cost))
println("files saved at percolation/RankDegreeAlpha")

