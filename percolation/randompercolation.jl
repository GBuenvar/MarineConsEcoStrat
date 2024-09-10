##
using CSV, DataFrames, Random, Statistics, Plots, ArgParse
include("percolation_functions.jl")
##
s = ArgParseSettings()
@add_arg_table! s begin
    "--nrep", "-n"
        help = "Number of repetitions"
        default = 10
        arg_type = Int64
    "--seed", "-s"
        help = "Seed for the random number generator"
        default = 1234
        arg_type = Int64
    "--threshold", "-t"
        help = "Threshold for the number of individuals of a species that need to be protected"
        default = 0.5
        arg_type = Float64
    "--rich_protect", "-r"
        help = "Start protecting the rich countries"
        default = true
        arg_type = Bool
end

p = parse_args(ARGS, s)

n = p["nrep"]
seed = p["seed"]
threshold = p["threshold"]
rich_protect = p["rich_protect"]

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

mkpath("percolation/random")

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

function plot_trajectories_median(prot_number, median_prot_number, n; xlabel = "# EEZs protected", ylabel = "Fraction of protected individuals", legend = false, dpi = 300, namesave = "none")
    p1 = plot(xlabel = xlabel, ylabel = ylabel, legend = legend, dpi = dpi)
    nn = size(prot_number)[1]
    for i in 1:nn
        plot!(p1, cumsum(prot_number[i, :]) ./ n)
    end
    plot!(p1, median_prot_number[1, :] ./ n, lw = 1.5, color=:black)
    plot!(p1)
    if namesave != "none"
        savefig(p1, namesave*".png")
        savefig(p1, namesave*".pdf")
    end
    return p1
end

##
##
start_protecting = [0, 8]
suffix = ""
if rich_protect
    rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
    start_protecting = rich
    println("Rich countries: ", rich)
    suffix = "_rich"
end

println("compuning random percolation...")
protected_times, protected_number = @time random_perc(agg_data, n, start_protecting = start_protecting, verbose=true)
median_protected_number = @time median_protected(protected_number)



##
p_inds = plot_trajectories_median(protected_number, median_protected_number, N; xlabel = "protecting EEZs", ylabel = "Fraction of protected individuals", legend = false, dpi = 300, namesave = "percolation/figures/random_protected_ids")
# save outputs in compressed files
CSV.write("percolation/random/protected_times$suffix.csv.gz", DataFrame(protected_times, :auto))
CSV.write("percolation/random/protected_number$suffix.csv.gz", DataFrame(protected_number, :auto))
CSV.write("percolation/random/median_protected_number$suffix.csv.gz", DataFrame(median_protected_number, :auto))

##

println("computing random percolation for species...")
species_id, protected_species_number, protected_species_times = @time protected_species_random(protected_number, protected_times, id_to_species_int, newids)
median_protected_species_number = @time median_protected(protected_species_number)

CSV.write("percolation/random/protected_species_number$suffix.csv.gz", DataFrame(protected_species_number, :auto))
CSV.write("percolation/random/protected_species_times$suffix.csv.gz", DataFrame(protected_species_times, :auto))
CSV.write("percolation/random/median_protected_species_number$suffix.csv.gz", DataFrame(median_protected_species_number, :auto))


##
p_species = plot_trajectories_median(protected_species_number, median_protected_species_number, N_species; xlabel = "# EEZs protected", ylabel = "Fraction of protected species", legend = false, dpi = 300, namesave = "percolation/figures/random_protected_species")
