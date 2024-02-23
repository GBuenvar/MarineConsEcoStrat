using CSV, DataFrames, Random, Statistics, Plots, XLSX, ArgParse
include("percolation_functions.jl")

##

# Create a parser for the command line arguments

s = ArgParseSettings()
@add_arg_table! s begin
    "--nreps", "-n"
        help = "number of repetitions for the stochastic outcomes"
        default = 1
        arg_type = Int64
    "--saving", "-s"
        help = "save the results to a file"
        default = false
        arg_type = Bool
end

p = parse_args(ARGS, s)

nreps = p["nreps"]
saving = p["saving"]

##
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))
eez_to_int = Dict(zip(eez_codes.EEZ, eez_codes.Int))

species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))

# read the eez_to_iso3 file
eez_to_iso3_data = CSV.read("data/eez_to_iso3.csv", DataFrame)
eez_to_iso3 = Dict(zip(eez_to_iso3_data.Country, eez_to_iso3_data.ISO_3digit))
# add High Seas
eez_to_iso3["-1"] = "-1"

mkpath("percolation/comb_prob")

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


# Each individual i visits a number of EEZs n_i. The probability of haunting that 
# animal is uniformly distributed among the visited EEZs. The EEZ that succeeds in
# haunting gets a reward of 1. 

# For a EEZ that is visited by N_a animals, the probability of haunting all the animals
# gets a reward of N_a is

# P(N_a) = \prod_{i=1}^{N_a} p_i <= 1

# where p_i is the probability of haunting the i-th animal.

# Then, if it is proposed to conserve in exchange of a reward N_a, it should accept
# since it is the maximum reward that it can get.

# However, if the reward is r<N_a, then the EEZ will accept the proposal deppending
# on the probability of obtaining the reward r or more if it continues haunting.

# To compute this probability, we have to compute the number of situations in which
# the reward is r or more and divide it by the total number of situations.

# the probability of haunting none is:

# P(0) = \prod_{i=1}^{N_a} (1-p_i)

# the probability of haunting one is:

# P(1) = \sum_{i=1}^{N_a} p_i \prod_{j=1, j\neq i}^{N_a} (1-p_j)

# the probability of haunting two is:

# P(2) = \sum_{i=1}^{N_a} \sum_{j=i+1}^{N_a} p_i p_j \prod_{k=1, k\neq i, k\neq j}^{N_a} (1-p_k)

# and so on.

# The probability of haunting r is:

# P(r) = \sum_{i_1=1}^{N_a} \sum_{i_2=i_1+1}^{N_a} ... \sum_{i_r=i_{r-1}+1}^{N_a} p_{i_1} p_{i_2} ... p_{i_r} \prod_{j=1, j\neq i_1, j\neq i_2, ..., j\neq i_r}^{N_a} (1-p_j)


# Numerically, we can compute all situations and then compute the probability of getting
# a reward r or more.

# For this, each situation is encoded as a binary number of N_a bits, where 1 means that
# the animal is haunted and 0 means that it is not haunted. 

# The vector P is the vector of probabilities of haunting each animal, and Ϙ
# is the vector of the probability of not haunting each animal (1.-P).

# we compute all the possible combinations of 1s and 0s of length N_a, and then we
# compute the probability of each combination.



function bitwise_sum(bitvec1::BitVector, bitvec2::BitVector=BitVector(digits(1, base=2, pad=length(bitvec1))[end:-1:1]))::BitVector
    # Check if both bitvectors have the same length
    if length(bitvec1) != length(bitvec2)
        throw(ArgumentError("Bitvectors must have the same length"))
    end
    
    carry = false
    result = BitVector(zeros(length(bitvec1)))
    
    for i in reverse(1:length(bitvec1))
        # XOR operation for the sum
        sum_bit = bitvec1[i] ⊻ bitvec2[i] ⊻ carry
        # Update carry for the next iteration
        carry = (bitvec1[i] & bitvec2[i]) | (bitvec1[i] & carry) | (bitvec2[i] & carry)
        # Update result bit
        result[i] = sum_bit
    end
    
    return result
end

##

function possible_outcomes(data, eezs, neez)

    # for neez in 1:length(eezs)
    println("Processing $(eezs[neez])")
    # 1- find the newid that visits neez
    visiting_ids = unique(data[data.EEZ .== eezs[neez], :newid])
    N_a = length(visiting_ids)
    # 2- find the number of visited EEZ for each visiting_ids and compute their p_i form it
    nvisits = combine(groupby(data[data.newid .∈ [visiting_ids,], :], :newid), nrow)
    
    P = 1 ./ nvisits.nrow
    Ϙ = 1 .- P

    # 3- compute all the possible combinations of 1s and 0s of length N_a
    probs_N = zeros(Float64, N_a + 1) # probability of getting 0, 1, 2, ..., N_a
    situation = BitVector(digits(0, base=2, pad=N_a))
    nii = 0
    println("N_a = $N_a possible situations = $(2^N_a)")
    nn = 0
    while nii < N_a # all situations except all 1s

        nii = count(situation) # number of individuals haunted
        prob_situation = prod(P[situation]) * prod(Ϙ[.!situation])
        probs_N[nii + 1] += prob_situation
        situation = bitwise_sum(situation) # next situation
        nn += 1
        if nn == 10000000
            break
        end


    end

    # all 1s situation
    probs_N[N_a+1] += prod(P)
    
    println("done")
    # write the probabilities to a file
    # open("percolation/comb_prob/probs.csv", "a") do f
    #     write(f, "$(eezs[neez]), $(join(probs_N, ","))\n")
    # end
end


# @time p = possible_outcomes(agg_data, eezs, neez)
##
"""
    N_haunting_events(newid::Int64, data::DataFrame; rng=Xoshiro(newid), nreps::Int64=1000)::Vector{Int64}

Compute `nrep` haunting events for the individual `newid`. The haunting events considers
that all the EEZs visited by the individual have the same probability of haunting it.

# Arguments
- `newid::Int64`: the id of the individual
- `data::DataFrame`: any dataframe that contains the columns `newid` and `EEZ`, with
    the id of the individual and the id of the EEZ, respectively, with one unique entry
    for each visit of the individual to the EEZ. (normally only 1 per individual and EEZ)
- `rng::AbstractRNG=Xoshiro(newid)`: the random number generator to use. The default is
    to use the Xoshiro generator with the seed equal to the `newid`.
- `nreps::Int64=1000`: the number of repetitions to compute.
"""
function N_haunting_events(newid::Int64, data::DataFrame; rng=Xoshiro(newid), nreps::Int64=1000)::Vector{Int64}
    eez_i = data[data.newid .== newid, :EEZ]
    return rand(rng, eez_i, nreps)
end


"""
    succesfull_histogram(events, eezs, N_as)::Matrix{Float64}

Compute the histogram of the number of successful haunting events for each EEZ.

# Arguments
- `events::Matrix{Int64}`: a matrix of size `N x nreps` where `N` is the number of
    individuals and `nreps` is the number of repetitions. Each column `i` contains the
    number of haunting events for the individual `i`.
- `eezs::Vector{Int64}`: the vector of the EEZs.
- `N_as::Vector{Int64}`: the vector of the number of animals visiting each EEZ.

# Returns
- `probs_N::Matrix{Float64}`: a matrix of size `max(N_as) x length(eezs)` where each
    column `i` contains the probability of getting 'x' successful haunting events for
    the EEZ `i`.
"""
function succesfull_histogram(events, eezs, N_as)
    max_haunting = maximum(N_as)
    nreps = size(events, 2)
    probs_N = zeros(Float64, max_haunting + 1, length(eezs))
    @views for (i, eez) in enumerate(eezs)
        # count the number of success at each repetition
        success = [count(events[:, i] .== eez) for i in 1:nreps]
        prob_n = [count(success .== j) for j in 0:N_as[i]] 
        probs_N[1:length(prob_n), i] = prob_n
    end 
    return probs_N ./ nreps
end


"""
    stochastic_outcomes(agg_data::DataFrame; nreps::Int64=1000)

Compute the stochastic outcomes of the haunting events for each EEZ.

# Arguments
- `agg_data::DataFrame`: any dataframe that contains the columns `newid` and `EEZ`, with
    the id of the individual and the id of the EEZ, respectively, with one unique entry
    for each visit of the individual to the EEZ. (normally only 1 per individual and EEZ)
- `nreps::Int64=1000`: the number of repetitions to compute.

# Returns
- `probs_N::Matrix{Float64}`: a matrix of size `max(N_as) x length(eezs)` where each
    column `i` contains the probability of getting 'x' successful haunting events for
    the EEZ `i`.
- `cum_p::Matrix{Float64}`: a matrix of size `max(N_as) x length(eezs)` where each
    column `i` contains the cumulative probability of getting 'x' or more successful
    haunting events for the EEZ `i`.
"""
function stochastic_outcomes(agg_data::DataFrame; nreps::Int64=1000)
    eezs = unique(agg_data[:, :EEZ])
    newids = unique(agg_data[:, :newid])
    N_as = [count(agg_data.EEZ .== eez) for eez in eezs] # number of animals visiting each eez

    # 1- Compute nreps haunting events for each newid, 
    println("Computing haunting events")
    events = zeros(Int64, length(newids), nreps)
    @views for (i, id) in enumerate(newids)
        events[i, :] = N_haunting_events(id, agg_data, nreps=nreps)
    end

    # 2- Compute the number of successful haunting events for each eez
    println("Computing successful haunting events")
    probs_N = succesfull_histogram(events, eezs, N_as)
    cum_p =  1. .- cumsum(probs_N, dims=1)

    # 3- plot the results and save the figure
    println("Plotting the results")
    sort_idx = sortperm(N_as, rev=true) # sort the eezs by the number of animals visiting them
    p1 = plot(labels = false, title="prob of haunting r or more", xlabel="r/r_max", ylabel="P(X>r)", palette=:batlow, dpi=300)
    p2 = plot(labels = false, title="prob of haunting r", xlabel="r/r_max", ylabel="P(X=r)", palette=:batlow, dpi=300)
    for i in sort_idx
        p_i = probs_N[:, i]
        cum_p_i = cum_p[:, i]
        x = collect(0:(length(p_i)-1))./ length(p_i)
        if any((1. .- cum_p_i) .>= 1.)
            idx = findfirst((1. .- cum_p_i) .>= 1.)
            p_i = p_i[1:idx]
            cum_p_i = cum_p_i[1:idx]
            x = collect(0:(idx-1))./idx
        end
        plot!(p2, x, p_i, labels=false)
        plot!(p1, x, cum_p_i, labels=false)
    end
    p3 = plot(p2,p1, size=(1200, 400), leftmargin=25Plots.px, bottommargin=20Plots.px, legend=false, dpi=300)
    println("done!")
    savefig(p3, "percolation/comb_prob/haunting_probability.pdf")
    savefig(p3, "percolation/comb_prob/haunting_probability.png")
    return probs_N, cum_p, p3
end

@time p, cum_p, p3 = stochastic_outcomes(agg_data, nreps=nreps)

##
