##
using CSV, DataFrames, Random, Statistics, Plots
##
trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))
species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

# aggregate the data by individual, species and EEZ, and compute the total time spent in each EEZ.
# we are not interested in the rest of fields
agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
unique_pairs = unique(agg_data[:, ["newid", "Species", "EEZ"]])

# First mode: random list of EEZs, 

function random_perc_1_rep(data, rng)
    eezlist = unique(data[:, :EEZ])
    random_eezlist = shuffle(rng, eezlist)
    unprotected_ids = unique(data[:, :newid])
    protected_ids = Vector{Int64}(undef, 0)
    protencted_numer = Vector{Int64}(undef, 0)
    protected_times = Vector{Int64}(undef, 0)
    for (t, eez) in enumerate(random_eezlist)
        new_data = data[data[:, :EEZ] .!= eez, :]    # protect a new EEZ
        new_unprotected_ids = unique(new_data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            append!(protected_ids, new_protected_ids)
            append!(protencted_numer, length(new_protected_ids))
            new_protected_times = ones(length(new_protected_ids)) .* t # add the time at which they were protected
            append!(protected_times, new_protected_times)
            end
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    end
    return protected_ids, protected_times, protencted_numer
end
n = 1000

# create a function that calls to the previous function n times, saves the results as Vector{Vector{Int64}} and returns the median of the number of protected individuals at each time

function random_perc(data, n)
    rng = MersenneTwister(1234)
    protected_ids = Vector{Vector{Int64}}(undef, n)
    protected_times = Vector{Vector{Int64}}(undef, n)
    protected_number = Vector{Vector{Int64}}(undef, n)
    for i in 1:n
        protected_ids[i], protected_times[i], protected_number[i] = random_perc_1_rep(data, rng)
    end
    return protected_ids, protected_times, protected_number
end

protected_ids, protected_times, protected_number = @time random_perc(agg_data, n)



# compute the median of the number of protected individuals at each time.
# for this, create a vector for each possible time values present on every simulation.
# Then, for each time value, compute the median of the number of protected individuals across all simulations
##
function median_protected(prot_times)
    flatten_times = vcat(prot_times...)
    max_time = maximum(flatten_times)
    min_time = minimum(flatten_times)
    # create a vector with all the possible time values
    all_t = collect(min_time:max_time)
    # create a vector to store the median of the number of protected individuals at each time
    median_prot = Vector{Float64}(undef, length(all_t))
    cum_prot = length.(findall.(==(min_time), prot_times))
    median_prot[1] = median(cum_prot)
    @views for (i, t) in enumerate(all_t[2:end])
        # get the indices of the times that are equal to t
        indices = findall.(==(t), prot_times)
        # get the number of protected individuals at those times
        cum_prot .+= length.(indices)
        # println(cum_protected)
        median_prot[i+1] = median(cum_prot)
        # break
    end
    return all_t, median_prot, cum_prot

end

all_times, median_protected_numbers, cum_protected = @time median_protected(protected_times)

# plot the median of the number of protected individuals at each time
plot(xlabel = "Time", ylabel = "Number of protected individuals", legend = false)
for i in 1:n
    plot!(unique(protected_times[i]), cumsum(protected_number[i]))
end
plot!(all_times, median_protected_numbers, xlabel = "Time", ylabel = "Number of protected individuals", legend = false, color=:black)
plot!()





# plot how the number of protected individuals increases with ti
