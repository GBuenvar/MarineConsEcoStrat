using CSV, DataFrames, Random, Statistics, Plots, XLSX, ArgParse, LaTeXStrings
include("percolation_functions.jl")


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
# eez_to_iso3["-1"] = "-1"

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


function ind_protection_odds(N, m)
    # the number of ways to choose the m visited eezs
    # out of the n picked eezs divided by the
    # number of ways to choose the m visited eezs
    # out of the N total eezs
    P = [n < m ? 0.0 :
    binomial(big(n), big(m)) / binomial(big(N), big(m))
    for n in 0:N]
    P = convert(Vector{Float64}, P)
    return P
end

function proxy_odds(N, m)
    # the number of ways to choose the m visited eezs
    # out of the n picked eezs divided by the
    # number of ways to choose the m visited eezs
    # out of the N total eezs
    p=m/N
    P = [n < m ? 0.0 :
    (n/N)^m
    for n in 1:N]
    return P
end

"""
    protection_odds(agg_data)

Compute the expected number of protected eezs and its variance

Protection for an individual that visits m zones when n are protecting 
is a binomial process with two states, protected and not protected, 
with a probability p=P^m(n) of being protected given by the `proxy_odds` function.
the expected value and its variance are those of a binomial distribution.

μ_m = n_m p_m

σ_m^2 = n_m p_m(1-p_m)

To compute the total expected number of protected individuals we sum over all
populations of individuals that visited the same number of eezs.

μ = Σ μ_m = Σ n_m p_m

σ^2 = Σ σ_m^2 = Σ n_m p_m(1-p_m)
"""
function proxy_protection_odds(agg_data)
    # newids = unique(agg_data[:, :newid])
    ntot = length(eezs)
    μ = zeros(ntot+1) # μ=μ(n)=<s^m>(n) = s^m * P^m(n) 
    σ2 = zeros(ntot+1) # σ^2=σ^2(n)= s^m * P^m(n) (1-P^m(n)
    m_of_ids = combine(groupby(agg_data, :newid), nrow)[!, :nrow]
    for m in range(1, maximum(m_of_ids))
        Size_m = sum( m_of_ids .== m)
        # P^m(n)
        P_m_at_n = proxy_odds(ntot, m)
        # moments
        μ .+= Size_m .* P_m_at_n
        σ2 .+= Size_m .* P_m_at_n .* (1 .- P_m_at_n)
    end
    return μ, σ2
end



"""
    protection_odds(agg_data)

Compute the expected number of protected eezs and its variance

Protection for an individual that visits m zones when n are protecting 
is a binomial process with two states, protected and not protected, 
with a probability p=P^m(n) of being protected given by the `ind_protection_odds` function.
the expected value and its variance are those of a binomial distribution.

μ_m = n_m p_m

σ_m^2 = n_m p_m(1-p_m)

To compute the total expected number of protected individuals we sum over all
populations of individuals that visited the same number of eezs.

μ = Σ μ_m = Σ n_m p_m

σ^2 = Σ σ_m^2 = Σ n_m p_m(1-p_m)
"""
function protection_odds(agg_data)
    eezs = unique(agg_data[:, :EEZ])
    ntot = length(eezs)
    μ = zeros(ntot+1) # μ=μ(n)=<s^m>(n) = s^m * P^m(n) 
    σ2 = zeros(ntot+1) # σ^2=σ^2(n)= s^m * P^m(n) (1-P^m(n)
    m_of_ids = combine(groupby(agg_data, :newid), nrow)[!, :nrow]
    for m in range(1, maximum(m_of_ids))
        Size_m = sum(m_of_ids .== m)
        # P^m(n)
        P_m_at_n = ind_protection_odds(ntot, m)
        # moments
        μ .+= Size_m .* P_m_at_n
        σ2 .+= Size_m .* P_m_at_n .* (1 .- P_m_at_n)
    end
    return μ, σ2
end

eezsvisited = combine(groupby(agg_data, :newid), nrow)[!, :nrow]
##



# esto solo es para mostrar que el proxy funciona
p1 = plot(xlabel = "Protecting EEZs", ylabel="Expected number of protected Individuals")
p2 = plot()

protected_number_mean, protected_number_var = protection_odds(agg_data)
protected_number_std = sqrt.(protected_number_var) 
top_bound = protected_number_mean .+ protected_number_std
bottom_bound = protected_number_mean .- protected_number_std
plot!(p1, bottom_bound ./ N, fillrange = top_bound ./ N, label=L"\mu\pm \sigma^2", c=:blue, fillalpha=0.3, alpha=0)
plot!(p1, protected_number_mean ./ N, label=L"\mu", c=:black)
# plot!(p1, bottom_bound, c=:black, label=nothing)
savefig(p1, "percolation/figures/uncorrelated_random.png")
savefig(p1, "percolation/figures/uncorrelated_random.pdf")



# plot!(p1, sum([proxy_odds(180, n) .* sum(eezsvisited .== n) for n in 1:180]), label="<proxy>", c=:red)

plot(p1)

##

p2 = plot()
for m in 1:maximum(m_of_ids)
    plot!(p2, sum(m_of_ids .== m)^2 .* ind_protection_odds(180, m), label="<m=$m>, $(sum(m_of_ids .== m))", c=:black)
end
plot!(p2)
# plot!(p2, [proxy_odds(180, n) for n in 1:180], label=nothing, c=:red)

# plot!(p2, [[],[]], label=["P^m(n)"], c=[:blue])

# plot(p1, p2)? 


##

function species_protection_odds(agg_data, threshold=0.5)
    eezs = unique(agg_data[:, :EEZ])
    
    species = unique(agg_data[:, :Species])
    species_threshold = [ceil(
        sum(
            unique(
                agg_data[:, [:Species, :newid]]
                )[:, :Species] .== s
            ) * threshold
            ) for s in species]
    prot_species = zeros(length(eezs)+1)
    for (s, st) in zip(species, species_threshold)
        s_data = @view agg_data[agg_data[:, :Species] .== s, :]
        members_protected_mean, members_protected_var = protection_odds(s_data)
        prot_species .+= members_protected_mean .>= st
    end
    return prot_species    
end

prot_species = species_protection_odds(agg_data)

pspecies = plot(xlabel="EEZ", ylabel="Protected Species")
plot!(pspecies, prot_species, label="Protected Species", c=:black)
##


