using CSV, DataFrames, Random, Statistics, Plots, XLSX, ArgParse
include("percolation_functions.jl")
##

# In this case, we are looking at which individuals are easier to protect, i.e.,
# those that visit fewer EEZs, P. From the EEZs visited by the individuals P that are easier to protect,
# the ones visited by the most individuals P are protected first. Each time an EEZ is protected, the individuals
# that are easier to protect are recalculated, and the process is repeated. The output is the number of individuals protected at each step
# and the number of EEZs protected before each individual is protected.

# Initially, Antarctica, the High Seas, and the EEZs of rich countries are protected (or left out)

s = ArgParseSettings()
@add_arg_table! s begin
    "--rich_protect", "-r"
        help = "Start protecting the rich countries"
        default = false
        arg_type = Bool
end

p = parse_args(ARGS, s)
rich_protect = p["rich_protect"]
suffix = rich_protect ? "_rich" : ""


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

mkpath("percolation/IndCenter")

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

rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
rich = rich_protect ? rich : [0,8]

protected_times, protected_number = easier_ind_protect(agg_data, start_protecting = rich)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
CSV.write("percolation/IndCenter/protected_times$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/IndCenter/protected_number$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/IndCenter/protected_species_times$suffix.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/IndCenter/protected_species_number$suffix.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/IndCenter")



p1 = plot(xlabel = "cooperating EEZs", ylabel = "Fraction of protected")
title!("Easier individuals first")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
savefig(p1, "percolation/figures/IndCenter$suffix.png")
savefig(p1, "percolation/figures/IndCenter$suffix.pdf")
