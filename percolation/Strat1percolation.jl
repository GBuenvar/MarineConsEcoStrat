##
using CSV, DataFrames, Random, Statistics, Plots
##

# En este lo que se trata es de mirar cuántos individuaos visitan cada EEZ, 
# en funcion de eso, proteger en orden ascendente/descendente. Cada vez que se protege,
# se mira cuántos individuos se han protegido. Igual que en la aleatoria, el output es
# el numero de individuos protegidos en cada paso y el numero de EEZ protegidas antes de
# que cada individuo se proteja.


##
# open the full_data_inds.csv.gz file

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

mkpath("percolation/Strat1")


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


function sorted_percolation(data, ord; start_protecting=[-1])
    ids = unique(data[:, :newid])

    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    sorted_EEZs = count_and_sort(unique_pairs, ord)
    eezlist = sorted_EEZs[:, :EEZ]
    iterated_eezs = setdiff(eezlist, start_protecting)

    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, length(iterated_eezs)+1)


    # Protect the first EEZs
    data = data[data[:, :EEZ] .∈ (iterated_eezs, ), :] # protect the EEZ
    new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
    new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
    if length(new_protected_ids) > 0
        # add the new protected individuals to the list
        prot_number[1] = length(new_protected_ids)
        prot_times[ids .∈ (new_protected_ids,)] .= 0 # add the time at which they were protected
    end
    unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    println(length(new_protected_ids))
    # iterate over the rest of EEZs, updating the eezlist at each step
    for (tt, eez) in enumerate(iterated_eezs)
        tt_idx = tt + 1
        data = data[data[:, :EEZ] .!= eez, :]    # protect a new EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        println(length(new_protected_ids))
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[tt_idx] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= tt # add the time at which they were protected
        end
        println("EEZ: ", eez, " time: ", tt, " # protected: ", length(new_protected_ids))
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

protected_times, protected_number = @time sorted_percolation(agg_data, "desc")
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("percolation/Strat1/protected_times_asc.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat1/protected_number_asc.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat1/protected_species_times_asc.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat1/protected_species_number_asc.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat1")

p1 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Descending order")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p1, "percolation/figures/Strat1_descending.png")
savefig(p1, "percolation/figures/Strat1_descending.pdf")
plot(p1)
println("figures saved at percolation/figures")
##

protected_times, protected_number = @time sorted_percolation(agg_data, "asc")
prot_species_number, prot_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)

# save protected times and number, and species times and number
CSV.write("percolation/Strat1/protected_times_desc.csv.gz", DataFrame(protected_times=protected_times))
CSV.write("percolation/Strat1/protected_number_desc.csv.gz", DataFrame(protected_number=protected_number))
CSV.write("percolation/Strat1/protected_species_times_desc.csv.gz", DataFrame(prot_species_times=prot_species_times))
CSV.write("percolation/Strat1/protected_species_number_desc.csv.gz", DataFrame(prot_species_number=prot_species_number))
println("files saved at percolation/Strat1")

p2 = plot(xlabel = "# EEZs protected", ylabel = "Fraction of protected")
title!("Ascending order")
plot!(p2, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p2, cumsum(prot_species_number)./N_species, label="species (50% of individuals)", color = "red")
savefig(p2, "percolation/figures/Strat1_ascending.png")
savefig(p2, "percolation/figures/Strat1_ascending.pdf")
plot(p2)
println("figures saved at percolation/figures")
