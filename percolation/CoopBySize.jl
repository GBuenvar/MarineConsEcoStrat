using CSV, DataFrames, Random, Statistics, Plots, XLSX
include("percolation_functions.jl")
##

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

folder = "percolation/SizeCoop$suffix"
mkpath(folder)

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




""" ASCENDING """
##
protected_times, protected_number = @time sorted_percolation(agg_data, "asc"; start_protecting=rich)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("$folder/protectedtimes_asc$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("$folder/protected_number_asc$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("$folder/protected_species_times_asc$suffix.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("$folder/protected_species_number_asc$suffix.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at $folder")

p1 = plot(xlabel = "cooperating EEZs", ylabel = "Fraction of protected")
title!("Ascending order")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p1, "percolation/figures/SizeCoop_ascending$suffix.png")
savefig(p1, "percolation/figures/SizeCoop_ascending$suffix.pdf")
plot(p1)
# println("figures saved at percolation/figures")



##
"""DESCENDING"""

 
protected_times, protected_number = @time sorted_percolation(agg_data, "desc"; start_protecting=rich)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("$folder/protected_times_desc$suffix.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("$folder/protected_number_desc$suffix.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("$folder/protected_species_times_desc$suffix.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("$folder/protected_species_number_desc$suffix.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at $folder")

p2 = plot(xlabel = "cooperating EEZs", ylabel = "Fraction of protected")
title!("Descending order")
plot!(p2, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p2, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p2, "percolation/figures/SizeCoop_descending$suffix.png")
savefig(p2, "percolation/figures/SizeCoop_descending$suffix.pdf")
plot(p2)
# println("figures saved at percolation/figures")

