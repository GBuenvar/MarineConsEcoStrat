##
using CSV, DataFrames, Random, Statistics, Plots
include("percolation_functions.jl")
##

# In this case, we are looking at which individuals are easier to protect, i.e.,
# those that visit fewer EEZs, P. From the EEZs visited by the individuals P that are easier to protect,
# the ones visited by the most individuals P are protected first. Each time an EEZ is protected, the individuals
# that are easier to protect are recalculated, and the process is repeated. The output is the number of individuals protected at each step
# and the number of EEZs protected before each individual is protected.

##
# open the full_data_inds.csv.gz file

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))
mkpath("percolation/Strat2")

##
# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))
##


##

protected_times, protected_number = easier_ind_protect(agg_data)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
CSV.write("percolation/Strat2/protected_times.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat2/protected_number.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat2/protected_species_times.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat2/protected_species_number.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat2")



p1 = plot(xlabel = "cooperating EEZs", ylabel = "Fraction of protected")
title!("Easier individuals first")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
savefig(p1, "percolation/figures/Strat2.png")
savefig(p1, "percolation/figures/Strat2.pdf")
