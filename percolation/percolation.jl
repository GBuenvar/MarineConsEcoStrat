##
using CSV, DataFrames, Random, Statistics, Plots, ArgParse
##
s = ArgParseSettings()
@add_arg_table! s begin
    "--nrep", "n"
        help = "Number of repetitions"
        default = 100
        arg_type = Int64
    "--seed", "s"
        help = "Seed for the random number generator"
        default = 1234
        arg_type = Int64
    "--threshold", "t"
        help = "Threshold for the number of individuals of a species that need to be protected"
        default = 0.5
        arg_type = Float64
    "--data", "d"
        help = "Path to the data file"
        default = "data/full_data_inds.csv.gz"
        arg_type = String
    "--eez", "e"
        help = "Path to the eez file"
        default = "data/eez_to_int.csv"
        arg_type = String
    "--species", "sp"
        help = "Path to the species file"
        default = "data/species_to_int.csv"
        arg_type = String
    "--output", "o"
        help = "Path to the output folder"
        default = "percolation"
        arg_type = String
end



# open the full_data_inds.csv.gz file

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))


##

# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
unique_pairs = unique(agg_data[:, ["newid", "Species", "EEZ"]])
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
# First mode: random list of EEZs, 

n = 100
N = length(newids)
##
function random_perc_1_rep(data; rng = MersenneTwister(1234), newids = unique(data[:, :newid]), eezlist = unique(data[:, :EEZ]))

    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(newids))
    prot_number = zeros(Int64, length(eezlist))

    # High Seas is protected from the beginning
    data = data[data[:, :EEZ] .!= -1, :] # protect High Seas
    new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
    new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
    if length(new_protected_ids) > 0
        # add the new protected individuals to the list
        prot_number[1] = length(new_protected_ids)
        prot_times[newids .∈ (new_protected_ids,)] .= 0 # add the time at which they were protected
    end
    random_eezlist = shuffle(rng, eezlist[eezlist .!= -1])
    @views for (t, eez) in enumerate(random_eezlist)
        data = data[data[:, :EEZ] .!= eez, :]    # protect a new EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[t+1] = length(new_protected_ids)
            prot_times[newids .∈ (new_protected_ids,)] .= t # add the time at which they were protected
            end
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        
    end
    return prot_times, prot_number
end


##

# create a function that calls to the previous function n times, saves the results as Vector{Vector{Int64}} and returns the median of the number of protected individuals at each time

function random_perc(data, n; seed = 1234)
    rng = MersenneTwister(seed)
    newids = unique(data[:, :newid])
    eezlist = unique(data[:, :EEZ])
    prot_times = zeros(Int64, (n, length(newids)))
    prot_number = zeros(Int64, (n, length(eezlist)))
    for i in 1:n
        prot_times[i, :], prot_number[i, :] = random_perc_1_rep(data; rng = rng, newids = newids, eezlist = eezlist)
    end   
    return prot_times, prot_number
end

protected_times, protected_number = @time random_perc(agg_data, n)

##
function median_protected(prot_number)
    # compute the cumulative sum of each row
    cum_prot = cumsum(prot_number, dims = 2)
    # compute the median of the number of protected individuals at each time
    median_prot = median(cum_prot, dims = 1)
    return median_prot
end

median_protected_number = @time median_protected(protected_number)
##

# save outputs in compressed files
CSV.write("percolation/protected_times.csv.gz", DataFrame(protected_times, :auto))
CSV.write("percolation/protected_number.csv.gz", DataFrame(protected_number, :auto))
CSV.write("percolation/median_protected_number.csv.gz", DataFrame(median_protected_number, :auto))


##
p1 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected individuals", legend = false)
for i in 1:n
    plot!(p1, cumsum(protected_number[i, :]) ./ N)
end
plot!(p1, median_protected_number[1, :] ./ N, lw = 1.5, color=:black)
plot!(p1)
savefig(p1, "percolation/figures/random_protected_ids.pdf")

##


function protected_species(prot_number, prot_times, dict_id_species, newids; threshold = 0.5)
    species = [dict_id_species[id] for id in newids]
    unique_species = unique(species)
    threshold_species = [Int64(round(threshold * sum(species .== sp))) for sp in unique_species]
    nn = size(prot_number)[1]
    prot_species_number = zeros(Int64, size(prot_number))
    prot_species_times  = zeros(Int64, (nn, length(unique_species)))
    for ii in 1:nn # for each simulation
        for (sp_idx, sp) in enumerate(unique_species) # for each species 
            sp_times = prot_times[ii, species .== sp] 
            sp_threshold = threshold_species[sp_idx]
            n_prot_sp = 0
            t_prot = 1
            while (n_prot_sp < sp_threshold) || (n_prot_sp == 0)
                n_prot_sp = sum(sp_times .<= t_prot)
                t_prot += 1
            end
            prot_species_times[ii, sp_idx] = t_prot
            prot_species_number[ii, t_prot] += 1
        end
    end
    return species, prot_species_number, prot_species_times
end

species_id, protected_species_number, protected_species_times = @time protected_species(protected_number, protected_times, id_to_species_int, newids)
median_protected_species_number = @time median_protected(protected_species_number)

plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected species", legend = false, dpi = 300)
for i in 1:n
    plot!(cumsum(protected_species_number[i, :]) ./ length(unique(species_id)))
end
plot!(median_protected_species_number[1, :] ./ length(unique(species_id)), lw = 1.5, color=:black)
savefig("percolation/figures/random_protected_species.pdf")

##

# save outputs in compressed files
CSV.write("percolation/protected_species_number.csv.gz", DataFrame(protected_species_number, :auto))
CSV.write("percolation/protected_species_times.csv.gz", DataFrame(protected_species_times, :auto))
CSV.write("percolation/median_protected_species_number.csv.gz", DataFrame(median_protected_species_number, :auto))