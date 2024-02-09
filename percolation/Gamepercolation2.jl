using CSV, DataFrames, Random, Statistics, Plots, XLSX
# include("percolation_functions.jl")
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




##
# Write three functions to make the plots p1, pcost and p3

function plot_protected(prot_number, prot_species, n, n_species, title; xlabel = "EEZs cooperating", ylabel = "Fraction protected", savename = "none", title_location = :center)
    p1 = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(p1, cumsum(prot_number)./ n, label="individuals", color = "black")
    plot!(p1, cumsum(prot_species)./n_species, label="species (50%)", color = "red")
    title!(p1, title, title_location = title_location)
    if savename != "none"
        savefig(p1, "percolation/figures/$(savename).pdf")
        savefig(p1, "percolation/figures/$(savename).png")
    end
    plot!(p1)
    return p1
end

function plot_protection_cost(prot_cost, title; xlabel = "EEZs cooperating", ylabel = "Spents", savename = "none", title_location = :center)
    pcost = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(pcost, prot_cost, label="Per EEZ", color = "black")
    plot!(pcost, cumsum(prot_cost), label="Cumulative", color = "blue")
    title!(pcost, title, title_location = title_location)
    if savename != "none"
        savefig(pcost, "percolation/figures/$(savename).pdf")
        savefig(pcost, "percolation/figures/$(savename).png")
    end
    plot!(pcost)
    return pcost
end

function plot_protection_cost_per_individual(prot_cost, prot_number, prot_species, title; xlabel = "EEZs cooperating", ylabel = "Cost/individuals", savename = "none", title_location = :center)
    p3 = plot(xlabel = xlabel, ylabel = ylabel)
    p3twinx = twinx(p3)
    plot!(p3twinx, ylabel = "Cost/species")
    plot!(p3, cumsum(prot_cost[2:end])./(cumsum(prot_number[2:end])), label="individuals", color = "black")
    plot!(p3twinx, cumsum(prot_cost[2:end])./(cumsum(prot_species[2:end])), label=false, color = "red")
    plot!(p3, [-1], [0], label="species", color = "red")
    title!(p3, title, title_location = title_location)
    # plot!(p3, margin = 5) 
    if savename != "none"
        savefig(p3, "percolation/figures/$(savename).pdf")
        savefig(p3, "percolation/figures/$(savename).png")
    end
    plot!(p3)
    return p3
end


function plot_number_ids(data, prot_eez; xlabel = "EEZs cooperating", ylabel = "Number of individuals", label = "Individuals", color = "black", title="", title_location = ":center")
    ids_at_eez = [sum(data[:, :EEZ] .== eez) for eez in prot_eez]
    p4 = plot(ids_at_eez, xlabel = xlabel, ylabel = ylabel, label = label, color = color, title=title, title_location = title_location)
    return p4
end
###

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


#### Game functions

function compute_id_weight(pairs; α=1)
    ids = unique(pairs[:, :newid])
    Neez_ids = zeros(Float64, length(ids))
    for (i, id) in enumerate(ids)
        Neez_ids[i] = length(unique(pairs[pairs[:, :newid] .== id, :EEZ]))
    end
    return Neez_ids.^(-α)
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

function compute_final_payoff(pairs)
    eezs = unique(pairs[:, :EEZ])
    final_po = [sum(pairs[:, :EEZ] .== eez) for eez in eezs]
    return final_po
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



###
##

operation = -

##

rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
EEZ_neis_dict = compute_neighbors(agg_data)
α = 1
print("Incentives decrease with neighbors (a=$α)")
@time protected_times, protected_number, protected_cost, protected_eezs = run_game_incentives(agg_data, start_protecting = rich, α=α, operation = operation)
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species ($(round(100 * protected_species_number[1]/sum(protected_species_number), digits=1))%) and $(protected_number[1]) individuals ($(round(100 * protected_number[1]/sum(protected_number), digits=1))%) are protected") 

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "A", title_location = (-0.1,1.1))
pcost = plot_protection_cost(protected_cost, "B", title_location = (-0.1,1.1))
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, "C", title_location = (-0.1,1.1))

eez_ids = [sum(agg_data[:, :EEZ] .== eez) for eez in protected_eezs if !(eez ∈ rich)]
p4 = plot_number_ids(agg_data, eez_ids, title="D", title_location = (-0.1,1.1))
# eez_neighbors = [length(EEZ_neis_dict[eez]) for eez in protected_eezs if !(eez ∈ rich)]
# plot!(twinx(p4), eez_neighbors, label = "Neighbors", color = :blue, ylabel = "Number of neighbors")
# p4 = plot(eez_neighbors, color=:blue, label = "Neighbors", ylabel = "Number of neighbors")

p_title = plot(title = "Incentives decrease with neighbors (a=$α)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
ptot = plot(p_title, p1, pcost, p3, p4, layout = @layout([A{0.01h}; [B C]; [D E]]), size = (600, 600))
savefig(ptot, "percolation/figures/Game_Incentives_decrease_with_neighbors.pdf")
savefig(ptot, "percolation/figures/Game_Incentives_decrease_with_neighbors.png")
plot(ptot)

##

α = 0
print("Same incentive always (a=$α)")
@time protected_times, protected_number, protected_cost, protected_eezs = run_game_incentives(agg_data, start_protecting = rich, α=α, operation = operation)
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species ($(round(100 * protected_species_number[1]/sum(protected_species_number), digits=1))%) and $(protected_number[1]) individuals ($(round(100 * protected_number[1]/sum(protected_number), digits=1))%) are protected") 

p1 = plot_protected(protected_number, protected_species_number, N, N_species, "A", title_location = (-0.1,1.1))
pcost = plot_protection_cost(protected_cost, "B", title_location = (-0.1,1.1))
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, "C", title_location = (-0.1,1.1))

ids_at_eez = [sum(agg_data[:, :EEZ] .== eez) for eez in protected_eezs if !(eez ∈ rich)]
p4 = plot(ids_at_eez, xlabel = "EEZs cooperating", ylabel = "Size of the EEZ", label = "#Individuals", color = "black", title="D", title_location = (-0.1,1.1))


# Show the plot
p_title = plot(title = "Same incentive always (a=$α)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
ptot = plot(p_title, p1, pcost, p3, p4, layout = @layout([A{0.01h}; [B C]; [D E]]), size = (600, 600))
savefig(ptot, "percolation/figures/Game_Same_incentive_always.pdf")
savefig(ptot, "percolation/figures/Game_Same_incentive_always.png")
plot(ptot)
##

α = -1
print("Incentives increase with neighbors (a=$α)")
@time protected_times, protected_number, protected_cost, protected_eezs = run_game_incentives(agg_data, start_protecting = rich, α=α, operation = operation)
protected_species_number, protected_species_times = protected_species(protected_number, protected_times, id_to_species_int, newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) species ($(round(100 * protected_species_number[1]/sum(protected_species_number), digits=1))%) and $(protected_number[1]) individuals ($(round(100 * protected_number[1]/sum(protected_number), digits=1))%) are protected") 


p1 = plot_protected(protected_number, protected_species_number, N, N_species, "A", title_location = (-0.1,1.1))
pcost = plot_protection_cost(protected_cost, "B", title_location = (-0.1,1.1))
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, "C", title_location = (-0.1,1.1))

ids_at_eez = [sum(agg_data[:, :EEZ] .== eez) for eez in protected_eezs if !(eez ∈ rich)]
p4 = plot(ids_at_eez, xlabel = "EEZs cooperating", ylabel = "Number of individuals", label = "Individuals", color = "black", title="D", title_location = (-0.1,1.1))


p_title = plot(title = "Incentives increase with neighbors (a=$α)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
ptot = plot(p_title, p1, pcost, p3, p4, layout = @layout([A{0.01h}; [B C]; [D E]]), size = (600, 600))
savefig(ptot, "percolation/figures/Game_Incentives_increase_with_neighbors.pdf")
savefig(ptot, "percolation/figures/Game_Incentives_increase_with_neighbors.png")
plot(ptot)

##
"""
Incentives decrease with neighbors (a=1)  1.211168 seconds (1.62 M allocations: 2.421 GiB, 17.93% gc time)
Total cost: 3779.413636363636
Initially 93 species are protected with 9594 individuals

Same incentive always (a=0)  0.444331 seconds (697.43 k allocations: 680.187 MiB, 21.31% gc time)
Total cost: 5043.0
Initially 93 species are protected with 9594 individuals

Incentives increase with neighbors (a=-1)  3.629252 seconds (3.98 M allocations: 7.056 GiB, 18.90% gc time)
Total cost: 9533.0
Initially 93 species are protected with 9594 individuals
"""