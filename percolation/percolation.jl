##
using CSV, DataFrames, Random
##

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

# First mode: random list of EEZs, 

function random_perc_1_rep(data)
    eezlist = unique(data[:, :EEZ])
    random_eezlist = shuffle(eezlist)
    unprotected_ids = unique(data[:, :newid])
    protected_ids = []
    protected_times = []
    for (t, eez) in enumerate(random_eezlist)
        new_data = data[data[:, :EEZ] .!= eez, :]    # protect a new EEZ
        new_unprotected_ids = unique(new_data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        push!(protected_ids, new_protected_ids) # add the new protected individuals to the list
        new_protected_times = ones(length(new_protected_ids)) * t # add the time at which they were protected
        push!(protected_times, new_protected_times)
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    end
    return protected_ids, protected_times
end
n = 1000

function random_perc(data, n)
    unique_ids = unique(data[:, :newid])
    protected_ids = zeros(Int, (length(unique_ids), n))
    protected_times = zeros(Int, (length(unique_ids), n))
    for i in 1:n
        new_protected_ids, new_protected_times = random_perc_1_rep(data)
        
        
    end
    return protected_ids, protected_times
end

protected_ids, protected_times = random_perc(agg_data, 3)

##


