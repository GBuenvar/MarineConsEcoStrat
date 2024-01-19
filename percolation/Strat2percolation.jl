##
using CSV, DataFrames, Random, Statistics, Plots
##

# En este caso miramos cuales son los individuos que son más fáciles de proteger, es decir,
# los que visitan menos EEZs, P. De las EEZs visitadas por los individuos P que sean más fáciles de proteger,
# se protegen primero las que más individuos P visitan. Cada vez que se protege una EEZ, se recalculan los individuos
# que son más fáciles de proteger, y se repite el proceso. El output es el numero de individuos protegidos en cada paso
# y el numero de EEZ protegidas antes de que cada individuo se proteja.

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

function easier_ind_protect(data)
    ids = unique(data[:, :newid])
    Neez = length(unique(data[:, :EEZ]))
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, length(unique(data[:, "EEZ"])))

    # High Seas is protected from the beginning
    data = data[data[:, :EEZ] .!= -1, :] # protect High Seas
    new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
    new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
    if length(new_protected_ids) > 0
        # add the new protected individuals to the list
        prot_number[1] = length(new_protected_ids)
        prot_times[newids .∈ (new_protected_ids,)] .= 0 # add the time at which they were protected
    end
    unprotected_ids = new_unprotected_ids # update the list of unprotected individuals

    for i in 2:Neez
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        # 1- Find those individuals that are easier to protect
        ids_eez_count = combine(groupby(unique_pairs, :newid), nrow)
        easy_ids = ids_eez_count[ids_eez_count[:, :nrow] .== minimum(ids_eez_count[:, :nrow]), :newid]

        # 2- Find the EEZs that are visited by those individuals, 
        # identify the EEZ that is visited by the most individuals that are easier to protect
        unique_pairs_easy = unique_pairs[unique_pairs[:, :newid] .∈ (easy_ids, ), :]
        easy_ids_eezs = combine(groupby(unique_pairs_easy, :EEZ), nrow)
        protect_eez = easy_ids_eezs[easy_ids_eezs[:, :nrow] .== maximum(easy_ids_eezs[:, :nrow]), :EEZ][1]

        # 3- Protect that EEZ
        data = data[data[:, :EEZ] .!= protect_eez, :] # protect the EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[i] = length(new_protected_ids)
            prot_times[newids .∈ (new_protected_ids,)] .= i-1 # add the time at which they were protected
        end
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    end
    return prot_times, prot_number
end

function protected_species(prot_number, prot_times, dict_id_species, newids; threshold = 0.5)
    species = [dict_id_species[id] for id in newids]
    unique_species = unique(species)
    threshold_species = [Int64(floor(threshold * sum(species .== sp))) for sp in unique_species]

    prot_species_number = zeros(Int64, size(prot_number))
    prot_species_times  = zeros(Int64, length(unique_species))
    for (sp_idx, sp) in enumerate(unique_species) # for each species 
        sp_times = prot_times[species .== sp] 
        sp_threshold = threshold_species[sp_idx]
        n_prot_sp = 0
        t_prot = 0
        n_prot_sp = sum(sp_times .<= t_prot)
        while (n_prot_sp < sp_threshold) || (n_prot_sp == 0)
            t_prot += 1
            n_prot_sp = sum(sp_times .<= t_prot)
        end
        prot_species_times[sp_idx] = t_prot
        prot_species_number[t_prot+1] += 1
    end
    return prot_species_number, prot_species_times
end





##

protected_times, protected_number = easier_ind_protect(agg_data)
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
CSV.write("percolation/Strat2/protected_times.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat2/protected_number.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat2/protected_species_times.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat2/protected_species_number.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat2")



p1 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Easier individuals first")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
savefig(p1, "percolation/figures/Strat2.png")
savefig(p1, "percolation/figures/Strat2.pdf")
