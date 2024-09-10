using CSV, DataFrames, Random, Statistics, Plots, XLSX


"""
    Rich_Poor_lists(eezlist, iso3_eez_list, income_data)

Compute the list of rich and poor EEZs based on the income data. 
The income data is a DataFrame with the columns "Code" and "Income group".
The function returns two lists, one with the rich EEZs and the other with the poor EEZs.
"""
function Rich_Poor_lists(eezlist, iso3_eez_list, income_data)::Tuple{Vector{Int64}, Vector{Int64}}
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


"""
    protected_species(prot_number, prot_times, dict_id_species, newids; threshold = 0.5)

Compute the number of individuals and species protected at each time step.
The function returns two vectors, one with the number of individuals protected at each time step, and the other with the time at which each species is protected.


# Arguments
- prot_number::Vector{Int64}: Number of individuals protected at each time step. The length
    of the vector is the number of time steps
- prot_times::Vector{Int64}: Time at which each individual is protected. The length of the vector
    is the number of individuals. If the individual is not protected, the time is set 
    to the length of the prot_number vector plus one.
- dict_id_species::Dict{Int64, String}: Dictionary with the species of each individual
- newids::Vector{Int64}: List of the individuals
- threshold::Float64: Fraction of individuals of a species that need to be protected

# Returns
- prot_species_number::Vector{Int64}: Number of species protected at each time step. 
    The length of the vector is the number of time steps
- prot_species_times::Vector{Int64}: Time at which each species is protected. 
    If the species is not protected, the time is set to the length of the 
    prot_number vector plus one.
    """
function protected_species(
    prot_number,
    prot_times,
    dict_id_species,
    newids;
    threshold = 0.5
    )
    species = [dict_id_species[id] for id in newids]
    unique_species = unique(species)
    threshold_species = [
        Int64(floor(threshold * sum(species .== sp)))
        for sp in unique_species
            ]
    t_non_protected = length(prot_number) + 1
    
    # Initialize the variables with no species protected
    prot_species_number = zeros(Int64, size(prot_number))
    prot_species_times  = fill(
        maximum(prot_number),
        length(unique_species)
        )
    for (sp_idx, sp) in enumerate(unique_species) # for each species 
        sp_times = prot_times[species .== sp] 
        sp_threshold = threshold_species[sp_idx]
        unique_sp_times = sort(unique(sp_times))
        sp_times_count = cumsum([count(==(i), sp_times) for i in unique_sp_times])
        t_prot = unique_sp_times[findfirst(>=(sp_threshold), sp_times_count)]
        if (t_prot != nothing) && (t_prot < t_non_protected)
            prot_species_times[sp_idx] = t_prot
            prot_species_number[t_prot+1] += 1
        end

    end
    return prot_species_number, prot_species_times
end

"""
    compute_neighbors(data)

Compute the neighbors of each EEZ. The function returns 
a dictionary with the EEZs as keys and the neighbors as values.

# Arguments
- data::DataFrame: DataFrame with the columns "newid" and "EEZ"

# Returns
- eez_neighbors::Dict{Int64, Set{Int64}}: Dictionary with the EEZs as keys and the neighbors as values
"""
function compute_neighbors(data::DataFrame)
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
    return df
end

function sorted_percolation(data, ord; start_protecting=[0, 8], verbose=false)
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

function easier_ind_protect(data; start_protecting = [0, 8])
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

function run_game(data; start_protecting = [0,8], α=1, order="higher")
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



function run_game_incentives(data; start_protecting = [0,8], α=1, operation = /)
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

function random_perc_1_rep(data; start_protecting = [0, 8], rng = MersenneTwister(1234), newids = unique(data[:, :newid]), eezlist = unique(data[:, :EEZ]))

    unprotected_ids  = unique(data[:, :newid])
    prot_times  = zeros(Int64, length(newids))
    prot_number = zeros(Int64, length(eezlist) - length(start_protecting) + 1)

    # Protect the first EEZs
    data = data[.!(data[:, :EEZ] .∈ (start_protecting, )), :] # protect the EEZ
    new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
    new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
    if length(new_protected_ids) > 0
        # add the new protected individuals to the list
        prot_number[1] = length(new_protected_ids)
        prot_times[newids .∈ (new_protected_ids,)] .= 0 # add the time at which they were protected
    end
    unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    random_eezlist = shuffle(rng, eezlist[eezlist .∉ (start_protecting, )])
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
function random_perc(data, n; start_protecting = [0, 8], seed = 1234, verbose=false)
    rng = MersenneTwister(seed)
    newids = unique(data[:, :newid])
    eezlist = unique(data[:, :EEZ])
    prot_times = zeros(Int64, (n, length(newids)))
    prot_number = zeros(Int64, (n, length(eezlist) - length(start_protecting) + 1))
    for i in 1:n
        verbose ? println("run $i") : nothing
        prot_times[i, :], prot_number[i, :] = random_perc_1_rep(data; start_protecting= start_protecting, rng = rng, newids = newids, eezlist = eezlist)
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
    t_non_protected = size(prot_number, 2) + 1
    threshold_species = [Int64(round(threshold * sum(species .== sp))) for sp in unique_species]
    nn = size(prot_number)[1]
    prot_species_number = zeros(Int64, size(prot_number))
    prot_species_times  = zeros(Int64, (nn, length(unique_species)))
    for ii in 1:nn # for each simulation
        for (sp_idx, sp) in enumerate(unique_species) # for each species 
            sp_times = prot_times[ii, species .== sp] 
            sp_threshold = threshold_species[sp_idx]
            unique_sp_times = sort(unique(sp_times))
            sp_times_count = cumsum([count(==(i), sp_times) for i in unique_sp_times])
            n_prot_sp = 0
            t_prot = unique_sp_times[findfirst(>=(sp_threshold), sp_times_count)]
            if (t_prot != nothing) && (t_prot < t_non_protected)
                prot_species_times[ii, sp_idx] = t_prot
                prot_species_number[ii, t_prot+1] += 1
            end
        end
    end
    return species, prot_species_number, prot_species_times
end



function run_diversity(data; start_protecting = [0, 8], α=1., q=1.)

    ids = unique(data[:, :newid])
    eezs = unique(data[:, :EEZ])
    iterated_eezs = setdiff(eezs, start_protecting)
    Neez = length(iterated_eezs)
    data_div = data # make a copy of the data to compute the diversity of the protected areas


    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = fill(Neez+1, length(ids)) 
    prot_number = zeros(Int64, Neez+1)
    prot_cost   = zeros(Float64, Neez+1)
    prot_eezs   = copy(start_protecting)
    global_div = zeros(Float64, Neez+1)
    eez_div    = zeros(Float64, Neez+1)

    
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
    unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])

    # Compute the initial diversity of the protected areas
    initial_global_div = compute_div(data_div[data_div.EEZ .∈ [start_protecting], :], q=q)
    mean_initial_div = mean([compute_div(data_div[data_div.EEZ .== eez, :], q=q) for eez in start_protecting])
    global_div[1] = initial_global_div
    eez_div[1] = mean_initial_div
    # 1- compute the diversity of each unprotected EEZ
    eezs_diversity = zeros(Float64, length(unprotected_eezs))
    eez_sizes = Vector{Vector{Float64}}(undef, length(unprotected_eezs))
    # Initialize a vector of vector{Float64} that potentially will be of different lengths    
    for (i, eez) in enumerate(unprotected_eezs)
        s = compute_species_dist(unique_pairs[unique_pairs.EEZ .== eez, :])
        eez_sizes[i] = s 
        eezs_diversity[i] = hill_number(s, q=q)
    end
    
    for ii in 2:Neez+1
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])

        # 2- find the EEZ with the highest diversity, and protect it. If there is a tie, protect the one with the highest diversity with q+1
        max_diversity = maximum(eezs_diversity)
        arg = findall(eezs_diversity .== max_diversity)
        to_protect = unprotected_eezs[arg]

        q0 = q
        while length(to_protect) > 1
            # if the distribution is exactly the same for all the eez, just pick one randomly
            ss = eez_sizes[arg] 
            if all([ss[1] == s for s in ss[2:end]])
                to_protect = [to_protect[rand(1:end)]]
                break
            end
            q0 += 1
            q0 > (q+10.) ? println("There is a hard tie, breaking it with q=$(q0)") : nothing
            eezs_diversity_tiebraker = zeros(Float64, length(to_protect))
            for (i, size) in enumerate(ss)
                eezs_diversity_tiebraker[i] = hill_number(size, q=q)
            end
            max_diversity = maximum(eezs_diversity_tiebraker)
            arg = findall(eezs_diversity_tiebraker .== max_diversity)
            to_protect = to_protect[arg]
        end
        protect_eez = to_protect[1]

        # 3- compute the id weight only for the individuals that are present in protect_eez
        visiting_ids = unique(unique_pairs[unique_pairs.EEZ .== protect_eez, :newid])
        eez_id_weights = compute_id_weight(unique_pairs[unique_pairs.newid .∈ (visiting_ids, ), :], α=α)

        # 4- compute the cooperation cost of the EEZ
        cost = sum(eez_id_weights)

        # 5- Protect that EEZ
        push!(prot_eezs, protect_eez)
        prot_cost[ii] = cost
        data = data[data[:, :EEZ] .!= protect_eez, :] # remove protect_eez from the data
        new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
        if length(new_protected_ids) > 0
            # add the new protected individuals to the list
            prot_number[ii] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= ii-1 # add the time at which they were protected
        end

        # drop the diversity of the protected eez
        eez_div[ii] = eezs_diversity[unprotected_eezs .== protect_eez][1]
        eezs_diversity = eezs_diversity[unprotected_eezs .!= protect_eez]
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
        # println(length(unique(data[:, :EEZ])))

        # 6- Compute the diversity of the protected areas
        global_div[ii] = compute_div(data_div[.!(data_div.EEZ .∈ [unprotected_ids, ]), :], q=q)
    end
    return prot_times, prot_number, prot_cost, prot_eezs, global_div, eez_div
end



# Diversity functions:
function hill_number(p::Vector{Float64}; q::Float64=0., normalization_threshold::Float64=10^-6)
    @assert abs(sum(p) - 1 < normalization_threshold) "probabilities are not normalized"
    @assert isinteger(q) "q must be a float type integer"
    p = p[p .> 0.]
    q == 1. ? (h = exp(-sum(p .* log.(p)))) : (h = sum(p.^q)^(1/(1-q)))
    return h
end

function compute_species_dist(df)
    # df is a dataframe with the columns newid, Species, EEZ
    # restricted to the eez of interest
    id_sp_pairs = unique(df[:, [:newid, :Species]])
    sp_size = combine(groupby(id_sp_pairs, :Species), nrow)
    sizes = Vector{Float64}(sp_size.nrow)
    sizes = sizes ./ sum(sizes)
    # sort the sizes
    sizes = sort(sizes, rev=true)
    return sizes
end


function compute_div(df; q=0.)
    # df is a dataframe with the columns newid, Species, EEZ
    # restricted to the eez of interest
    sizes = compute_species_dist(df)
    df_div = hill_number(sizes, q=q)
    return df_div
end


## Ranking based on the degree of the nodes in the projected network


function projected_graph(
    bipartite_df,
    projected_partition,
    complementary_partition,
    weight_column=nothing
    )

projected_df = DataFrame(
    node_A = Int[],
    node_B = Int[],
    weight = Float64[]
)
nodes = unique(bipartite_df[:, projected_partition])
@views for node in nodes
    # get the zones visited by the individual
    node_EEZs = bipartite_df[
        bipartite_df[:, projected_partition] .== node,
        complementary_partition
        ]
    # get the subset of the bipartite graph with the zones visited by the individual
    neighbouring_subset = bipartite_df[
        bipartite_df[:, complementary_partition] .∈ (node_EEZs,),
        :]
    
    if isnothing(weight_column)
        # count the number of times each individual visits a zone in the subset
        neighbors_times = combine(
            groupby(neighbouring_subset, :newid),
            nrow
            )
        rename!(neighbors_times, :nrow => :weight)
    else
        neighbors_times = combine(
            groupby(neighbouring_subset, :newid),
            weight_column => sum
            )
        rename!(neighbors_times, weight_column*"_sum" => :weight)
    end
    # remove the self node
    neighbors_times = neighbors_times[neighbors_times[:, :newid] .!= node, :]
    # add the node to the projected graph
    append!(
        projected_df, 
        DataFrame(
            node_A = fill(node, size(neighbors_times, 1)),
            node_B = neighbors_times[:, :newid],
            weight = neighbors_times[:, :weight]
        ))
end
return projected_df

end

function degree_ranking(
    projected_df::DataFrame,
    node_N=nothing
    )
    # calculate the degree of each node
    # by summing all the weights of the same node
    degree_df = combine(
        groupby(projected_df, :node_A),
        :weight => sum
    )
    sort!(degree_df, :node_A, rev=true)
    rename!(degree_df, :weight_sum => :degree)
    if !isnothing(node_N)
        degree_df[!, :degree] = degree_df[!, :degree] ./ node_N
    end
    # sort the degree_df
    sort!(degree_df, :degree, rev=true)

    return degree_df
end

##
"""
ranked_ids_remove_eezs(
    data::DataFrame,
    start_protecting::Vector{Int} = [0,8],
    include_eez_resistant::Bool=false,
    weight_column=nothing
    )

"""
function ranked_ids_remove_eezs(
    data::DataFrame,
    start_protecting::Vector{Int} = [0,8],
    include_eez_resistant::Bool=false,
    weight_column=nothing
    )

    ids::Vector{Int64} = unique(data[:, :newid])
    eezs::Vector{Int64} = unique(data[:, :EEZ])
    iterated_eezs::Vector{Int64} = setdiff(eezs, start_protecting)
    Neez::Int = length(iterated_eezs)

    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = fill(Neez+1, length(ids)) 
    prot_number = zeros(Int64, Neez+1)
    prot_eezs   = copy(start_protecting)
    eez_prot_times = fill(length(eezs)+1, length(eezs))
    eez_prot_times[eezs .∈ (start_protecting,)] .= 0

    # Protect the first EEZs
    data = data[data[:, :EEZ] .∈ (iterated_eezs, ), :]
    new_unprotected_ids = unique(data[:, :newid]) # get the list of the individuals that are still not protected
    new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected
    if length(new_protected_ids) > 0
        # add the new protected individuals to the list
        prot_number[1] = length(new_protected_ids)
        prot_times[ids .∈ (new_protected_ids,)] .= 0 # add the time at which they were protected
    end
    unprotected_ids = new_unprotected_ids # update the list of unprotected individuals
    unprotected_eezs = copy(iterated_eezs)
    ii = 2 # start in 2 becouse 1 is the start_protecting
    steps = 0
    while size(data, 1) > 0
        
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])

        # 1- compute the degree of each node and rank them
        if include_eez_resistant
            proj_df = projected_graph(data, :newid, :EEZ, weight_column)
            sorted_N = sort(
                combine(groupby(data, :newid), nrow),
                :newid
                )
            # remove the isolated nodes from the list 
            sorted_N = sorted_N[sorted_N.newid .∈ (proj_df.node_A, ), :][:, :nrow]
            higher_rank = degree_ranking(
                proj_df,
                sorted_N
            )
        else
            higher_rank = degree_ranking(
                projected_graph(data, :newid, :EEZ, weight_column),
                )

        end

        # if the projected network is no connected, simply return the degree_df
        # as the degree in the bipartite network (that is equivalent to rank by number of visited zones)
        if size(higher_rank, 1) == 0
            higher_rank = sort(
                combine(groupby(data, :newid), nrow), 
                :newid
                )
            rename!(higher_rank, :newid => :node_A)
        end
        higher_rank = higher_rank[1, :node_A]

        # 2- get the zones visited by the individual
        protect_eez = data[
            data[:, :newid] .== higher_rank,
            :EEZ
            ]

        # 3- protect the zones visited by the individual.
        # 3.1.- store the new cooperating eezs and remove them
        n_prot = length(protect_eez)
        prot_eezs = vcat(prot_eezs, protect_eez)
        eez_prot_times[eezs .∈ (protect_eez, )] .= steps
        data = data[data[:, :EEZ] .∉ (protect_eez, ), :]        
        
        # 3.2.- store the protected individuals
        # add the protected individual after all 
        # the new EEZs are cooperating
        new_unprotected_ids = unique(data[:, :newid])
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids)
        # index_add = ii + n_prot - 1 # if protected after all protect
        index_add = ii  # if protected before all protect
        prot_number[index_add] = length(new_protected_ids)
        prot_times[ids .∈ (new_protected_ids,)] .= index_add - 1

        # 4- update the variables of the simulation and continue
        unprotected_ids = new_unprotected_ids
        ii += n_prot
        steps += 1
    end
    return prot_times, prot_number, prot_eezs, eez_prot_times
end
