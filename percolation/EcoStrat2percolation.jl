using CSV, DataFrames, Random, Statistics, Plots, XLSX
include("percolation_functions.jl")
##

# In this case, we are looking at which individuals are easier to protect, i.e.,
# those that visit fewer EEZs, P. From the EEZs visited by the individuals P that are easier to protect,
# the ones visited by the most individuals P are protected first. Each time an EEZ is protected, the individuals
# that are easier to protect are recalculated, and the process is repeated. The output is the number of individuals protected at each step
# and the number of EEZs protected before each individual is protected.

# Initially, Antarctica, the High Seas, and the EEZs of rich countries are protected (or left out)


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



mkpath("percolation/Strat2Eco")

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

rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)

protected_times, protected_number = easier_ind_protect(agg_data, start_protecting = rich)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
CSV.write("percolation/Strat2Eco/protected_times.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat2Eco/protected_number.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat2Eco/protected_species_times.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat2Eco/protected_species_number.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat2Eco")



p1 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Easier individuals first")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
savefig(p1, "percolation/figures/Strat2Eco.png")
savefig(p1, "percolation/figures/Strat2Eco.pdf")
