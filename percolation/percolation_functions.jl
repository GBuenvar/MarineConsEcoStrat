using CSV, DataFrames, Random, Statistics, Plots, XLSX


function Rich_Poor_lists(eezlist, iso3_eez_list, income_data)
    Rich = Vector{Int64}(undef, 0)
    for (eez, iso3) in zip(eezlist, iso3_eez_list)
        if in(iso3, income_data[:, :Code])
            income = income_data[income_data[:, :Code] .== iso3, "Income group"][1]
            if (!ismissing(income)) && (income == "High income")
                push!(Rich, eez)
            end
        end
    end
    # Add High Seas and Antarctica
    push!(Rich, eez_to_int["-1"])
    push!(Rich, eez_to_int["Antarctica"])
    # every other EEZ is Poor
    Poor = setdiff(eezs, Rich)
    return Rich, Poor
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

function compute_neighbors(data)
    pairs = data[:, ["newid", "EEZ"]]
    unique_pairs = unique(pairs)
    eez_neighbors = Dict()
    for (i, eez) in enumerate(unique(unique_pairs[:, :EEZ]))
        visit_eez = unique_pairs[unique_pairs[:, :EEZ] .== eez, :newid]
        neis = unique_pairs[unique_pairs[:, :newid] .∈ (visit_eez, ), :EEZ]
        eez_neis = setdiff(unique(neis), [eez])
        eez_neighbors[eez] = eez_neis
    end
    return eez_neighbors
end


#### Stategy 1

# En este lo que se trata es de mirar cuántos individuaos visitan cada EEZ, 
# en funcion de eso, proteger en orden ascendente/descendente. Cada vez que se protege,
# se mira cuántos individuos se han protegido. Igual que en la aleatoria, el output es
# el numero de individuos protegidos en cada paso y el numero de EEZ protegidas antes de
# que cada individuo se proteja.



function count_and_sort(df, ord)
    df = combine(groupby(df, :EEZ), nrow)
    rename!(df, Dict(:nrow => :n_ids))
    sort!(df, :n_ids, rev = (ord == "desc"))
    df = vcat(df[df[:, :EEZ] .== -1, :], df[df[:, :EEZ] .!= -1, :])
    return df
end

function sorted_percolation(data, ord; start_protecting=[-1], verbose=false)
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
    # iterate over the rest of EEZs, updating the eezlist at each step
    for (tt, eez) in enumerate(iterated_eezs)
        tt_idx = tt + 1
        data = data[data[:, :EEZ] .!= eez, :]    # protect a new EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[tt_idx] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= tt # add the time at which they were protected
        end
        if verbose
            println("EEZ: ", eez, " time: ", tt, " # protected: ", length(new_protected_ids))
        end
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    end
    return prot_times, prot_number
end



#### Strategy 2

function easier_ind_protect(data; start_protecting = [-1])
    ids = unique(data[:, :newid])
    eezs = unique(data[:, :EEZ])
    iterated_eezs = setdiff(eezs, start_protecting)
    Neez = length(iterated_eezs)
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, Neez)

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



#### Game functions

function compute_id_weight(pairs; α=1)
    ids = unique(pairs[:, :newid])
    Neez_ids = zeros(Float64, length(ids))
    for (i, id) in enumerate(ids)
        Neez_ids[i] = length(unique(pairs[pairs[:, :newid] .== id, :EEZ]))
    end
    return Neez_ids.^(-α)
end

function compute_final_payoff(pairs)
    eezs = unique(pairs[:, :EEZ])
    final_po = [sum(pairs[:, :EEZ] .== eez) for eez in eezs]
    return final_po
end



function compute_eez_payoff(pairs, id_weights)
    ids = unique(pairs[:, :newid])
    eezs = unique(pairs[:, :EEZ])
    eezs_payoffs = zeros(Float64, length(eezs))
    for (ii, eez) in enumerate(eezs)
        visiting_ids = unique(pairs[pairs[:, :EEZ] .== eez, :newid])
        eez_po = sum(id_weights[ids .∈ (visiting_ids, )])
        eezs_payoffs[ii] = eez_po
    end
    return eezs_payoffs
end

function run_game(data; start_protecting = [-1], α=1, order="higher")
    ids = unique(data[:, :newid])
    eezs = unique(data[:, :EEZ])
    iterated_eezs = setdiff(eezs, start_protecting)
    Neez = length(iterated_eezs)

    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, Neez+1)
    prot_cost   = zeros(Float64, Neez+1)
    prot_eezs   = copy(start_protecting)

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

    unprotected_eezs = copy(iterated_eezs)
    for ii in 2:Neez+1
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])
        # 1- compute the ids weights
        id_weights = compute_id_weight(unique_pairs; α=α)
        # 2- compute the payoff of each EEZ
        eezs_payoffs = compute_eez_payoff(unique_pairs, id_weights)
        # 3- Find the EEZ with the highest payoff
        if order == "higher"
            arg = argmax(eezs_payoffs)
        elseif order == "lower"
            arg = argmin(eezs_payoffs)
        end
        cost = eezs_payoffs[arg]
        protect_eez = unprotected_eezs[arg]
        # 4- Protect that EEZ
        push!(prot_eezs, protect_eez)
        prot_cost[ii] = cost
        data = data[data[:, :EEZ] .!= protect_eez, :] # protect the EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[ii] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= ii-1 # add the time at which they were protected
        end
        # println("EEZ: ", protect_eez, " time: ", ii, " # protected: ", length(new_protected_ids), " # unprotected: ", length(new_unprotected_ids), " cost: ", cost)
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        # println(length(unique(data[:, :EEZ])))
    end
    return prot_times, prot_number, prot_cost, prot_eezs
end



function run_game_incentives(data; start_protecting = [-1], α=1, operation = /)
    ids = unique(data[:, :newid])
    eezs = unique(data[:, :EEZ])
    iterated_eezs = setdiff(eezs, start_protecting)
    Neez = length(iterated_eezs)

    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(ids))
    prot_number = zeros(Int64, Neez+1)
    prot_cost   = zeros(Float64, Neez+1)
    prot_eezs   = copy(start_protecting)

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

    unprotected_eezs = copy(iterated_eezs)

    # Compute the potential cost of each EEZ being the last to be protected
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
    last_cooperating_payoffs = compute_final_payoff(unique_pairs)
    @assert sum(last_cooperating_payoffs) == size(unique_pairs)[1] "The sum of the last cooperating payoffs is not equal to the number of unprotected individuals"
    for ii in 2:Neez+1
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])
        # 1- compute the ids weights
        id_weights = compute_id_weight(unique_pairs; α=α)
        # 2- compute the payoff of each EEZ
        eezs_payoffs = compute_eez_payoff(unique_pairs, id_weights)
        # 3- Compute the hurry of each EEZ to cooperate
        eezs_hurry = operation.(eezs_payoffs, last_cooperating_payoffs)
        # 4- find the EEZ with the highest hurry, and protect it. If there is a tie, protect the one with the highest payoff

        max_hurry = maximum(eezs_hurry)
        arg = findall(eezs_hurry .== max_hurry)
        if length(arg) > 1
                arg = arg[argmax(eezs_payoffs[arg])]
        end
        arg = arg[1]

        cost = eezs_payoffs[arg]
        protect_eez = unprotected_eezs[arg]

        # 5- update the last cooperating payoffs removing the protected EEZ
        deleteat!(last_cooperating_payoffs, arg)

        # 4- Protect that EEZ
        push!(prot_eezs, protect_eez)
        prot_cost[ii] = cost
        data = data[data[:, :EEZ] .!= protect_eez, :] # protect the EEZ
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[ii] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= ii-1 # add the time at which they were protected
        end
        # println("EEZ: ", protect_eez, " time: ", ii, " # protected: ", length(new_protected_ids), " # unprotected: ", length(new_unprotected_ids), " cost: ", cost)
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        # println(length(unique(data[:, :EEZ])))
    end
    return prot_times, prot_number, prot_cost, prot_eezs
end


## Random percolation

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
    unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
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


function median_protected(prot_number)
    # compute the cumulative sum of each row
    cum_prot = cumsum(prot_number, dims = 2)
    # compute the median of the number of protected individuals at each time
    median_prot = median(cum_prot, dims = 1)
    return median_prot
end


function protected_species_random(prot_number, prot_times, dict_id_species, newids; threshold = 0.5)
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
