##
using CSV, DataFrames, Random, Statistics, Plots
include("percolation_functions.jl")
##

# In this code, the goal is to examine how many individuals visit each EEZ. 
# Based on that, protection is applied in ascending/descending order. 
# Each time protection is applied, the number of protected individuals is recorded. 
# Similar to the random approach, the output includes the number of protected 
# individuals at each step and the number of EEZs protected before each individual is protected.


##
# open the full_data_inds.csv.gz file

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

mkpath("percolation/Strat1")


##
# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

##
protected_times, protected_number = @time sorted_percolation(agg_data, "desc")
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("percolation/Strat1/protected_times_asc.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat1/protected_number_asc.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat1/protected_species_times_asc.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat1/protected_species_number_asc.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat1")

p1 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Descending order")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p1, "percolation/figures/Strat1_descending.png")
savefig(p1, "percolation/figures/Strat1_descending.pdf")
plot(p1)
println("figures saved at percolation/figures")
##

protected_times, protected_number = @time sorted_percolation(agg_data, "asc")
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("percolation/Strat1/protected_times_desc.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat1/protected_number_desc.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat1/protected_species_times_desc.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat1/protected_species_number_desc.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat1")

p2 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Ascending order")
plot!(p2, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p2, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p2, "percolation/figures/Strat1_ascending.png")
savefig(p2, "percolation/figures/Strat1_ascending.pdf")
plot(p2)
println("figures saved at percolation/figures")
