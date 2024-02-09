using CSV, DataFrames, Random, Statistics, Plots, XLSX
include("percolation_functions.jl")
##

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

mkpath("percolation/Strat1Eco")


##
# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))


##
eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs]


##
rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)




""" ASCENDING """
##
protected_times, protected_number = @time sorted_percolation(agg_data, "asc"; start_protecting=rich)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("percolation/Strat1Eco/protectedtimes_asc.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat1Eco/protected_number_asc.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat1Eco/protected_species_times_asc.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat1Eco/protected_species_number_asc.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat1Eco")

p1 = plot(xlabel = "EEZs cooperating", ylabel = "Fraction of protected")
title!("Ascending order")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p1, "percolation/figures/Strat1Eco_ascending.png")
savefig(p1, "percolation/figures/Strat1Eco_ascending.pdf")
plot(p1)
# println("figures saved at percolation/figures")



##
"""DESCENDING"""

 
protected_times, protected_number = @time sorted_percolation(agg_data, "desc"; start_protecting=rich)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("percolation/Strat1Eco/protected_times_desc.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat1Eco/protected_number_desc.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat1Eco/protected_species_times_desc.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat1Eco/protected_species_number_desc.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat1Eco")

p2 = plot(xlabel = "EEZs cooperating", ylabel = "Fraction of protected")
title!("Descending order")
plot!(p2, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p2, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p2, "percolation/figures/Strat1Eco_descending.png")
savefig(p2, "percolation/figures/Strat1Eco_descending.pdf")
plot(p2)
# println("figures saved at percolation/figures")

