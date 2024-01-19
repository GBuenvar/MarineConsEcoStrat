using CSV, DataFrames, Random, Statistics, Plots, XLSX
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

# Read the economic data
economic_data = DataFrame(XLSX.readtable("data/CLASS.xlsx", "List of economies"))
eez_to_int = Dict(zip(eez_codes.EEZ, eez_codes.Int))

# read the eez_to_iso3 file
eez_to_iso3_data = CSV.read("data/eez_to_iso3.csv", DataFrame)
eez_to_iso3 = Dict(zip(eez_to_iso3_data.Country, eez_to_iso3_data.ISO_3digit))
# add High Seas
eez_to_iso3["-1"] = "-1"

α = 1



mkpath("percolation/Game")

##
# Since I am only interested in some specific fields of the data, I will create a new dataframe with only those fields

agg_data = combine(groupby(trajectories, [:newid, :Species, :EEZ]), "timestay (1/30days)" => sum)
@assert sum(agg_data[:, end])/30 == length(unique(agg_data[:, :newid])) "Total time (/30) in days is not equal to the number of individuals"
id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs]


##
function Rich_Poor_lists(eezlist, iso3_eez_list, income_data)
    Rich = Vector{Int64}(undef, 0)
    for (eez, iso3) in zip(eezlist, iso3_eez_list)
        if in(iso3, income_data[:, :Code])
            income = income_data[income_data[:, :Code] .== iso3, "Income group"][1]
            if (!ismissing(income)) && ((income == "High income") || (income == "Upper middle income"))
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



function compute_id_weight(pairs; α=1)
    ids = unique(pairs[:, :newid])
    Neez_ids = zeros(Float64, length(ids))
    for (i, id) in enumerate(ids)
        Neez_ids[i] = length(unique(pairs[pairs[:, :newid] .== id, :EEZ]))
    end
    return Neez_ids.^(-α)
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
        prot_species_number[t_prot] += 1
    end
    return prot_species_number, prot_species_times
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



## 

α = 1
print("Higher payoff first (a=$α)")
rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
protected_times1, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="higher")
protected_species_number1, protected_species_times1 = protected_species(protected_number, protected_times1, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number1[1]) species are protected with $(protected_number[1]) individuals") 

##
p1 = plot(xlabel = "EEZs protecting", ylabel = "Fraction protected")
title!("Higher payoff first (a=$α)")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(protected_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
##
pcost = plot(xlabel = "EEZs protecting", ylabel = "Spents")
title!("Higher payoff first (a=$α)")
plot!(pcost, protected_cost, label="Per EEZ", color = "black")
plot!(pcost, cumsum(protected_cost), label="Cumulative", color = "blue")
plot!(pcost)
##
p3 = plot(xlabel = "# EEZs protecting", ylabel = "Cost/individuals")
p3twinx = twinx(p3)
plot!(p3twinx, ylabel = "Cost/species")
title!(p3, "Higher payoff first (a=$α)")
plot!(p3, cumsum(protected_cost[2:end])./(cumsum(protected_number[2:end])), label="individuals", color = "black")
plot!(p3twinx, cumsum(protected_cost[2:end])./(cumsum(protected_species_number[2:end])), label=false, color = "red")
plot!(p3, [-1], [0], label="species", color = "red")
plot!(p3)
##

# now do the same for α=-1

α = -1
print("Higher payoff first (a=$α)")
protected_times2, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="higher")
protected_species_number2, protected_species_times2 = protected_species(protected_number, protected_times2, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number2[1]) species are protected with $(protected_number[1]) individuals") 

##
p1 = plot(xlabel = "EEZs protecting", ylabel = "Fraction protected")
title!("Higher payoff first (a=$α)")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(protected_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
##
pcost = plot(xlabel = "EEZs protecting", ylabel = "Spents")
title!("Higher payoff first (a=$α)")
plot!(pcost, protected_cost, label="Per EEZ", color = "black")
plot!(pcost, cumsum(protected_cost), label="Cumulative", color = "blue")
plot!(pcost)
##
p3 = plot(xlabel = "# EEZs protecting", ylabel = "Cost/individuals")
p3twinx = twinx(p3)
plot!(p3twinx, ylabel = "Cost/species")
title!(p3, "Hiigher payoff first (a=$α)")
plot!(p3, cumsum(protected_cost[2:end])./(cumsum(protected_number[2:end])), label="individuals", color = "black")
plot!(p3twinx, cumsum(protected_cost[2:end])./(cumsum(protected_species_number[2:end])), label=false, color = "red")
plot!(p3, [-1], [0], label="species", color = "red")
plot!(p3)


# repeat for order = "lower"

α = 1
print("Lower payoff first (a=$α)")
protected_times3, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="lower")
protected_species_number3, protected_species_times3 = protected_species(protected_number, protected_times3, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number3[1]) species are protected with $(protected_number[1]) individuals") 

##
p1 = plot(xlabel = "EEZs protecting", ylabel = "Fraction protected")
title!("Lower payoff first (a=$α)")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(protected_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
##
pcost = plot(xlabel = "EEZs protecting", ylabel = "Spents")
title!("Lower payoff first (a=$α)")
plot!(pcost, protected_cost, label="Per EEZ", color = "black")
plot!(pcost, cumsum(protected_cost), label="Cumulative", color = "blue")
plot!(pcost)
##
p3 = plot(xlabel = "# EEZs protecting", ylabel = "Cost/individuals")
p3twinx = twinx(p3)
plot!(p3twinx, ylabel = "Cost/species")
title!(p3, "Lower payoff first (a=$α)")
plot!(p3, cumsum(protected_cost[2:end])./(cumsum(protected_number[2:end])), label="individuals", color = "black")
plot!(p3twinx, cumsum(protected_cost[2:end])./(cumsum(protected_species_number[2:end])), label=false, color = "red")
plot!(p3, [-1], [0], label="species", color = "red")
plot!(p3)


##
α = -1
print("Lower payoff first (a=$α)")
protected_times4, protected_number, protected_cost, protected_eezs = run_game(agg_data, start_protecting = rich, α=α, order="lower")
protected_species_number4, protected_species_times4 = protected_species(protected_number, protected_times4, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number4[1]) species are protected with $(protected_number[1]) individuals") 
##
p1 = plot(xlabel = "EEZs protecting", ylabel = "Fraction protected")
title!("Lower payoff first (a=$α)")
plot!(p1, cumsum(protected_number)./ N, label="individuals", color = "black")
plot!(p1, cumsum(protected_species_number)./N_species, label="species (50% of individuals)", color = "red")
plot!(p1)
##
pcost = plot(xlabel = "EEZs protecting", ylabel = "Spents")
title!("Lower payoff first (a=$α)")
plot!(pcost, protected_cost, label="Per EEZ", color = "black")
plot!(pcost, cumsum(protected_cost), label="Cumulative", color = "blue")
plot!(pcost)
# plot the cost of protecting one species at a time
# El coste por individuo es cuánto
##
p3 = plot(xlabel = "# EEZs protecting", ylabel = "Cost/individuals")
p3twinx = twinx(p3)
plot!(p3twinx, ylabel = "Cost/species")
title!(p3, "Lower payoff first (a=$α)")
plot!(p3, cumsum(protected_cost[2:end])./(cumsum(protected_number[2:end])), label="individuals", color = "black")
plot!(p3twinx, cumsum(protected_cost[2:end])./(cumsum(protected_species_number[2:end])), label=false, color = "red")
plot!(p3, [-1], [0], label="species", color = "red")
plot!(p3)


