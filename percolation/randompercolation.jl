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
end

p = parse_args(ARGS, s)

n = p["nrep"]
seed = p["seed"]
threshold = p["threshold"]



# open the full_data_inds.csv.gz file

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))
mkpath("percolation/random")

##
# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
unique_pairs = unique(agg_data[:, ["newid", "Species", "EEZ"]])
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))
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
protected_times, protected_number = @time random_perc(agg_data, n)
median_protected_number = @time median_protected(protected_number)


##
p_inds = plot_trajectories_median(protected_number, median_protected_number, N; xlabel = "# EEZs protected", ylabel = "Fraction of protected individuals", legend = false, dpi = 300, namesave = "percolation/figures/random_protected_ids")
# save outputs in compressed files
CSV.write("percolation/random/protected_times.csv.gz", DataFrame(protected_times, :auto))
CSV.write("percolation/random/protected_number.csv.gz", DataFrame(protected_number, :auto))
CSV.write("percolation/random/median_protected_number.csv.gz", DataFrame(median_protected_number, :auto))

##
species_id, protected_species_number, protected_species_times = @time protected_species_random(protected_number, protected_times, id_to_species_int, newids)
median_protected_species_number = @time median_protected(protected_species_number)



##
p_species = plot_trajectories_median(protected_species_number, median_protected_species_number, N_species; xlabel = "# EEZs protected", ylabel = "Fraction of protected species", legend = false, dpi = 300, namesave = "percolation/figures/random_protected_species")
# save outputs in compressed files
CSV.write("percolation/random/protected_species_number.csv.gz", DataFrame(protected_species_number, :auto))
CSV.write("percolation/random/protected_species_times.csv.gz", DataFrame(protected_species_times, :auto))
CSV.write("percolation/random/median_protected_species_number.csv.gz", DataFrame(median_protected_species_number, :auto))