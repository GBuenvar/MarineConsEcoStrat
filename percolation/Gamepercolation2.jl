using CSV, DataFrames, Random, Statistics, Plots, XLSX
include("percolation_functions.jl")
##

# En este caso miramos cuales son los individuos que son más fáciles de proteger, es decir,
# los que visitan menos EEZs, P. De las EEZs visitadas por los individuos P que sean más fáciles de proteger,
# se protegen primero las que más individuos P visitan. Cada vez que se protege una EEZ, se recalculan los individuos
# que son más fáciles de proteger, y se repite el proceso. El output es el numero de individuos protegidos en cada paso
# y el numero de EEZ protegidas antes de que cada individuo se proteja.

##
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
α = 1
print("Incentives decrease with neighbors (a=$α)")
@time protected_times, protected_number, protected_cost, protected_eezs = run_game_incentives(agg_data, start_protecting = rich, α=α, operation = operation)
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(length(rich)-2) of $(length(rich)+length(poor)-2) countries collaborate, $(protected_species_number[1]) species ($(round(100 * protected_species_number[1]/sum(protected_species_number), digits=1))%) and $(protected_number[1]) individuals ($(round(100 * protected_number[1]/sum(protected_number), digits=1))%) are protected") 

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "A", title_location = (-0.1,1.1))
pcost = plot_protection_cost(protected_cost, "B", title_location = (-0.1,1.1))
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, "C", title_location = (-0.1,1.1))

eez_ids = [sum(agg_data[:, :EEZ] .== eez) for eez in protected_eezs if !(eez ∈ rich)]
p4 = plot_number_ids(agg_data, eez_ids, title="D", title_location = (-0.1,1.1))
# eez_neighbors = [length(EEZ_neis_dict[eez]) for eez in protected_eezs if !(eez ∈ rich)]
# plot!(twinx(p4), eez_neighbors, label = "Neighbors", color = :blue, ylabel = "Number of neighbors")
# p4 = plot(eez_neighbors, color=:blue, label = "Neighbors", ylabel = "Number of neighbors")

p_title = plot(title = "Incentives decrease with neighbors (a=$α)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
ptot = plot(p_title, p1, pcost, p3, p4, layout = @layout([A{0.01h}; [B C]; [D E]]), size = (600, 600))
savefig(ptot, "percolation/figures/Game_Incentives_decrease_with_neighbors.pdf")
savefig(ptot, "percolation/figures/Game_Incentives_decrease_with_neighbors.png")
plot(ptot)

##

α = 0
print("Same incentive always (a=$α)")
@time protected_times, protected_number, protected_cost, protected_eezs = run_game_incentives(agg_data, start_protecting = rich, α=α, operation = operation)
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species ($(round(100 * protected_species_number[1]/sum(protected_species_number), digits=1))%) and $(protected_number[1]) individuals ($(round(100 * protected_number[1]/sum(protected_number), digits=1))%) are protected") 

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "A", title_location = (-0.1,1.1))
pcost = plot_protection_cost(protected_cost, "B", title_location = (-0.1,1.1))
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, "C", title_location = (-0.1,1.1))

ids_at_eez = [sum(agg_data[:, :EEZ] .== eez) for eez in protected_eezs if !(eez ∈ rich)]
p4 = plot(ids_at_eez, xlabel = "EEZs cooperating", ylabel = "Size of the EEZ", label = "#Individuals", color = "black", title="D", title_location = (-0.1,1.1))


# Show the plot
p_title = plot(title = "Same incentive always (a=$α)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
ptot = plot(p_title, p1, pcost, p3, p4, layout = @layout([A{0.01h}; [B C]; [D E]]), size = (600, 600))
savefig(ptot, "percolation/figures/Game_Same_incentive_always.pdf")
savefig(ptot, "percolation/figures/Game_Same_incentive_always.png")
plot(ptot)
##

α = -1
print("Incentives increase with neighbors (a=$α)")
@time protected_times, protected_number, protected_cost, protected_eezs = run_game_incentives(agg_data, start_protecting = rich, α=α, operation = operation)
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species ($(round(100 * protected_species_number[1]/sum(protected_species_number), digits=1))%) and $(protected_number[1]) individuals ($(round(100 * protected_number[1]/sum(protected_number), digits=1))%) are protected") 


p1 = plot_protected(protected_number, protected_species_number, N, N_species, "A", title_location = (-0.1,1.1))
pcost = plot_protection_cost(protected_cost, "B", title_location = (-0.1,1.1))
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, "C", title_location = (-0.1,1.1))

ids_at_eez = [sum(agg_data[:, :EEZ] .== eez) for eez in protected_eezs if !(eez ∈ rich)]
p4 = plot(ids_at_eez, xlabel = "EEZs cooperating", ylabel = "Number of individuals", label = "Individuals", color = "black", title="D", title_location = (-0.1,1.1))


p_title = plot(title = "Incentives increase with neighbors (a=$α)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
ptot = plot(p_title, p1, pcost, p3, p4, layout = @layout([A{0.01h}; [B C]; [D E]]), size = (600, 600))
savefig(ptot, "percolation/figures/Game_Incentives_increase_with_neighbors.pdf")
savefig(ptot, "percolation/figures/Game_Incentives_increase_with_neighbors.png")
plot(ptot)

##
""" With substraction operation:
Incentives decrease with neighbors (a=1)  8.332424 seconds (8.84 M allocations: 25.834 GiB, 11.46% gc time)
Total cost: 5467.335109335109
Initially 55 of 177 countries collaborate, 72 species (64.9%) and 8379 individuals (66.4%) are protected

Same incentive always (a=0)  1.133673 seconds (1.72 M allocations: 2.159 GiB, 14.06% gc time, 17.14% compilation time)
Total cost: 7856.0
Initially 72 species (64.9%) and 8379 individuals (66.4%) are protected

Incentives increase with neighbors (a=-1)  1.556473 seconds (2.88 M allocations: 4.219 GiB, 14.92% gc time)
Total cost: 17578.0
Initially 72 species (64.9%) and 8379 individuals (66.4%) are protected

"""


""" With division operation:
Incentives decrease with neighbors (a=1)  8.332424 seconds (8.84 M allocations: 25.834 GiB, 11.46% gc time)
Total cost: 5467.335109335109
Initially 55 of 177 countries collaborate, 72 species (64.9%) and 8379 individuals (66.4%) are protected

Same incentive always (a=0)  1.133673 seconds (1.72 M allocations: 2.159 GiB, 14.06% gc time, 17.14% compilation time)
Total cost: 7856.0
Initially 72 species (64.9%) and 8379 individuals (66.4%) are protected

Incentives increase with neighbors (a=-1)  1.556473 seconds (2.88 M allocations: 4.219 GiB, 14.92% gc time)
Total cost: 17578.0
Initially 72 species (64.9%) and 8379 individuals (66.4%) are protected
"""
