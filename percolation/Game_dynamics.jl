using CSV, DataFrames, Random, Statistics, Plots, XLSX
include("percolation_functions.jl")
##

# En este caso miramos cuales son los individuos que son más fáciles de proteger, es decir,
# los que visitan menos EEZs, P. De las EEZs visitadas por los individuos P que sean más fáciles de proteger,
# se protegen primero las que más individuos P visitan. Cada vez que se protege una EEZ, se recalculan los individuos
# que son más fáciles de proteger, y se repite el proceso. El output es el numero de individuos protegidos en cada paso
# y el numero de EEZ protegidas antes de que cada individuo se proteja.

##
# open the full_data_inds.csv.gz file

# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))
eez_to_int = Dict(zip(eez_codes.EEZ, eez_codes.Int))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

# read the eez_to_iso3 file
eez_to_iso3_data = CSV.read("data/eez_to_iso3.csv", DataFrame)
eez_to_iso3 = Dict(zip(eez_to_iso3_data.Country, eez_to_iso3_data.ISO_3digit))
# read the economic data
economic_data = DataFrame(XLSX.readtable("data/CLASS.xlsx", "List of economies"))

mkpath("percolation/Game")

##
# Read the data, consisting on the newid, the species, the eez and the time spent in the eez    
agg_data = CSV.read("data/agg_data.csv", DataFrame)

id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs];

##

rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)


cum_p = CSV.read("percolation/comb_prob/cum_probs.csv", DataFrame)
cum_p = Matrix(cum_p)



function plot_protected(prot_number, prot_species, n, n_species, title; xlabel = "EEZs protecting", ylabel = "Fraction protected", savename = "none")
    p1 = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(p1, cumsum(prot_number)./ n, label="individuals", color = "black")
    plot!(p1, cumsum(prot_species)./n_species, label="species (50% of individuals)", color = "red")
    title!(p1, title)
    if savename != "none"
        savefig(p1, "percolation/figures/$(savename).pdf")
        savefig(p1, "percolation/figures/$(savename).png")
    end
    plot!(p1)
    return p1
end

function plot_protection_cost(prot_cost, title; xlabel = "EEZs protecting", ylabel = "Spents", savename = "none")
    pcost = plot(xlabel = xlabel, ylabel = ylabel)
    plot!(pcost, prot_cost, label="Per EEZ", color = "black")
    plot!(pcost, cumsum(prot_cost), label="Cumulative", color = "blue")
    title!(pcost, title)
    if savename != "none"
        savefig(pcost, "percolation/figures/$(savename).pdf")
        savefig(pcost, "percolation/figures/$(savename).png")
    end
    plot!(pcost)
    return pcost
end

function plot_protection_cost_per_individual(prot_cost, prot_number, prot_species, n, n_species, title; xlabel = "# EEZs protecting", ylabel = "Cost/individuals", savename = "none")
    p3 = plot(xlabel = xlabel, ylabel = ylabel)
    p3twinx = twinx(p3)
    plot!(p3twinx, ylabel = "Cost/species")
    plot!(p3, cumsum(prot_cost[2:end])./(cumsum(prot_number[2:end])), label="individuals", color = "black")
    plot!(p3twinx, cumsum(prot_cost[2:end])./(cumsum(prot_species[2:end])), label=false, color = "red")
    plot!(p3, [-1], [0], label="species", color = "red")
    title!(p3, title)
    # plot!(p3, margin = 5) 
    if savename != "none"
        savefig(p3, "percolation/figures/$(savename).pdf")
        savefig(p3, "percolation/figures/$(savename).png")
    end
    plot!(p3)
    return p3
end




"""
    compute_cooperating_probability(cumulative_prob::Matrix{Float64}, incentives::Vector{Float64}, max_reward::Int)

Computes the probability of cooperating for a given incentive for every zone.
The probability is computed using the probability of gettiing a reward higher
or equal to the incentive. 
"""
function compute_cooperating_probability(cumulative_prob::Matrix{Float64}, incentives::Vector{Float64}, max_reward::Int = size(cumulative_prob, 1))
    eezs_coop_odds = [incentive <= max_reward ?
        1. - 0.5*(cumulative_prob[floor(Int, incentive+1)] +
                  cumulative_prob[ceil(Int, incentive+1)]) :
        1.0
        for (incentive, cump) in zip(incentives, eachcol(cumulative_prob))
            ]
    return eezs_coop_odds
end

function incetives_game(data, prob; start_protecting = [0,8], α=0, order="higher", cooperating_threshold_p = 0.5)
    ids::Vector{Int64} = unique(data[:, :newid])
    eezs::Vector{Int64} = unique(data[:, :EEZ])
    iterated_eezs::Vector{Int64} = setdiff(eezs, start_protecting)
    Neez::Int = length(iterated_eezs)
    max_reward::Int = size(cum_p, 1)

    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    # this change respect to others. Instead of initialize the protecting
    # times to zeros, that means protected at the beginning, we initialize
    # them to the last iteration+1. Individuals that are not protected will
    # have this value.
    prot_times  = fill(Neez+1, length(ids)) 
    prot_number = zeros(Int64, Neez+1)
    prot_cost   = zeros(Float64, Neez+1)
    prot_eezs   = copy(start_protecting)

    # Protect the first EEZs
    unprotected_prob = prob[:, eezs .∈ (iterated_eezs, )]
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
    for ii in 2:Neez+1
        unique_pairs = unique(data[:, ["newid", "Species", "EEZ"]])
        unprotected_eezs = unique(unique_pairs[:, :EEZ])
        unprotected_prob = prob[:, eezs .∈ (unprotected_eezs, )]

        # 1- compute the weights of the individuals
        id_weights = compute_id_weight(unique_pairs; α=α)
        # 2- compute the incentives for each EEZ
        eezs_incentives = compute_eez_payoff(unique_pairs, id_weights)
        # write the incentives to a line in the file
        # 3- compute the probability of getting that incentive or more
        eezs_coop_odds = compute_cooperating_probability(
            unprotected_prob,
            eezs_incentives,
            max_reward)

        
        
        # 4- from those EEZs that have a probability of getting the incentive 
        # higher than the threshold, select the one with the highest incentive.
        # If all the probabilities are below the threshold, stop the game.
        willing_to_cooperate = eezs_coop_odds .> cooperating_threshold_p
        open("percolation/Game_incentives.csv", "a") do io
            write(io, join(string.(reverse(eezs_incentives[willing_to_cooperate])), ",")*"\n")
        end
        open("percolation/Game_coop_prob.csv", "a") do io
            write(io, join(string.(eezs_coop_odds[willing_to_cooperate]), ",")*"\n")
        end
        open("percolation/Game_willing_to_cooperate.csv", "a") do io
            sss = join(string.(eezs_incentives[willing_to_cooperate]), ",")*"\n"
            write(io,sss)
        end
        if all(!, willing_to_cooperate)
            println("All the probabilities are below the threshold, 
            stopping the game at iteration $ii.")
            println(eezs_incentives)
            println(eezs_coop_odds)
            println(unprotected_eezs)
            break
        end

        # 5- Protect that EEZ
        eezs_incentives = eezs_incentives[willing_to_cooperate]  
        if order == "higher"
            arg = argmax(eezs_incentives)
        elseif order == "lower"
            arg = argmin(eezs_incentives)
        end
        protect_eez = unprotected_eezs[willing_to_cooperate][arg]
        push!(prot_eezs, protect_eez)
        prot_cost[ii] = eezs_incentives[arg]

        open("percolation/Game_EEZSwilling_to_cooperate.csv", "a") do io
            sss = join(string.(unprotected_eezs[willing_to_cooperate]), ",")*":$protect_eez"*"\n"
            write(io,sss)
        end


        # Update the lists of protected and unprotected individuals
        data = data[data[:, :EEZ] .!= protect_eez, :]        
        new_unprotected_ids = unique(data[:, :newid])
        new_protected_ids = setdiff(unprotected_ids, new_unprotected_ids) # get the list of the individuals that are now protected

        if length(new_protected_ids) > 0
            prot_number[ii] = length(new_protected_ids)
            prot_times[ids .∈ (new_protected_ids,)] .= ii-1 # add the time at which they were protected
        end
        unprotected_ids = new_unprotected_ids # update the list of unprotected individuals

    end
    return prot_times, prot_number, prot_cost, prot_eezs
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
    prot_number, prot_times,
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


##
α = 0
protected_times, protected_number, protected_cost, protected_eezs =
 incetives_game(
    agg_data,
    cum_p;
    start_protecting = rich,
    α=α,
    order="higher",
    cooperating_threshold_p = 0.5)
    
protected_species_number, protected_species_times =
 protected_species(
    protected_number,
    protected_times,
    id_to_species_int,
    newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) of the $(N_species) species are protected with $(protected_number[1]) individuals")
println("At the end $(sum(protected_species_number[1:end])) of the $(N_species) species are protected with $(sum(protected_number[1:end])) individuals")
##
p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Game dynamics (a=$α)", savename = "game_dynamics_a_$(α)")
pcost = plot_protection_cost(protected_cost, "Game dynamics (a=$α)", savename = "game_dynamics_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Game dynamics (a=$α)", savename = "game_dynamics_cost_a_$(α)")
plot(p1,pcost,p3)
plot(p1)



##
α = 1
protected_times, protected_number, protected_cost, protected_eezs =
 incetives_game(
    agg_data,
    cum_p;
    start_protecting = rich,
    α=α,
    order="higher",
    cooperating_threshold_p = 0.5)
    
protected_species_number, protected_species_times =
 protected_species(
    protected_number,
    protected_times,
    id_to_species_int,
    newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) of the $(N_species) species are protected with $(protected_number[1]) individuals")
println("At the end $(sum(protected_species_number[1:end])) of the $(N_species) species are protected with $(sum(protected_number[1:end])) individuals")
##
p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Game dynamics (a=$α)", savename = "game_dynamics_a_$(α)")
pcost = plot_protection_cost(protected_cost, "Game dynamics (a=$α)", savename = "game_dynamics_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Game dynamics (a=$α)", savename = "game_dynamics_cost_a_$(α)")
plot(p1,pcost,p3)
plot(p1)

##
α = -1
protected_times, protected_number, protected_cost, protected_eezs =
 incetives_game(
    agg_data,
    cum_p;
    start_protecting = rich,
    α=α,
    order="higher",
    cooperating_threshold_p = 0.5)
    
protected_species_number, protected_species_times =
 protected_species(
    protected_number,
    protected_times,
    id_to_species_int,
    newids)
println("Total cost: ", sum(protected_cost))
println("Initially $(protected_species_number[1]) of the $(N_species) species are protected with $(protected_number[1]) individuals")
println("At the end $(sum(protected_species_number[1:end])) of the $(N_species) species are protected with $(sum(protected_number[1:end])) individuals")
##
p1 = plot_protected(protected_number, protected_species_number, N, N_species, "Game dynamics (a=$α)", savename = "game_dynamics_a_$(α)")
pcost = plot_protection_cost(protected_cost, "Game dynamics (a=$α)", savename = "game_dynamics_a_$(α)")
p3 = plot_protection_cost_per_individual(protected_cost, protected_number, protected_species_number, N, N_species, "Game dynamics (a=$α)", savename = "game_dynamics_cost_a_$(α)")
plot(p1,pcost,p3)
plot(p1)
