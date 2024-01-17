##
using CSV, DataFrames, Random, Statistics, Plots, ArgParse
##

# En este lo que se trata es de mirar cuántos individuaos visitan cada EEZ, 
# en funcion de eso, proteger en orden ascendente/descendente. Cada vez que se protege,
# se mira cuántos individuos se han protegido. Igual que en la aleatoria, el output es
# el numero de individuos protegidos en cada paso y el numero de EEZ protegidas antes de
# que cada individuo se proteja.

s = ArgParseSettings()
@add_arg_table! s begin
    "--order", "-o"
        help = "Order in which the EEZs are protected. asc for ascending, desc for descending"
        default = "asc"
        arg_type = String
end

p = parse_args(ARGS, s)
order = p["order"]
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
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

function count_and_sort(df, ord)
    df = combine(groupby(df, :EEZ), nrow)
    rename!(df, Dict(:nrow => :n_ids))
    sort!(df, :n_ids, rev = (ord == "desc"))
    df = vcat(df[df[:, :EEZ] .== -1, :], df[df[:, :EEZ] .!= -1, :])
    return df
end


function sorted_percolation(data, ord)
    ids = unique(data[:, :newid])

    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    sorted_EEZs = count_and_sort(unique_pairs, "asc")
    eezlist = sorted_EEZs[:, :EEZ]
    starting_list = sorted_EEZs[:, :EEZ]
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
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
    unprotected_ids = new_unprotected_ids
    
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    sorted_EEZs = count_and_sort(unique_pairs, ord)
    eezlist = sorted_EEZs[:, :EEZ]
    final_list = [-1]
    # iterate over the rest of EEZs, updating the eezlist at each step
    tt = 1
    while length(eezlist) > 0
        tt += 1
        eez = eezlist[1]
        data = data[data[:, :EEZ] .!= eez, :]    # protect a new EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[tt] = length(new_protected_ids)
            prot_times[newids .∈ (new_protected_ids,)] .= tt # add the time at which they were protected
        end
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        sorted_EEZs = count_and_sort(unique_pairs, ord)
        eezlist = sorted_EEZs[:, :EEZ]
        final_list = vcat(final_list, eez)
    end
    return prot_times, prot_number, starting_list, final_list
end



function protected_species(prot_number, prot_times, dict_id_species, newids; threshold = 0.5)
    species = [dict_id_species[id] for id in newids]
    unique_species = unique(species)
    threshold_species = [Int64(round(threshold * sum(species .== sp))) for sp in unique_species]
    prot_species_number = zeros(Int64, size(prot_number))
    prot_species_times  = zeros(Int64, length(unique_species))
    for (sp_idx, sp) in enumerate(unique_species) # for each species 
        sp_times = prot_times[species .== sp] 
        sp_threshold = threshold_species[sp_idx]
        n_prot_sp = 0
        t_prot = 1
        n_prot_sp = sum(sp_times .<= t_prot)
        while (n_prot_sp < sp_threshold) || (n_prot_sp == 0)
            t_prot += 1
            n_prot_sp = sum(sp_times .<= t_prot)
        end
        prot_species_times[sp_idx] = t_prot
        prot_species_number[t_prot] += 1
    end
    return species, prot_species_number, prot_species_times
end

protected_times, protected_number, starting_list, final_list = sorted_percolation(agg_data, "desc")
species_id, prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)


p1 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Descending order")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species", color = "red")
savefig(p1, "percolation/figures/Strat1_descending.png")
savefig(p1, "percolation/figures/Strat1_descending.pdf")
##

protected_times, protected_number, starting_list, final_list = sorted_percolation(agg_data, "asc")
species_id, prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)


p2 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Ascending order")
plot!(p2, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p2, cumsum(prot_species_number)./N_species, label="species", color = "red")
savefig(p2, "percolation/figures/Strat1_ascending.png")
savefig(p2, "percolation/figures/Strat1_ascending.pdf")