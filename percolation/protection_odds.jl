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
# read the economic data
economic_data = DataFrame(XLSX.readtable("data/CLASS.xlsx", "List of economies"))


mkpath("percolation/uncorrelated_random/")

##
# Read the data, consisting on the newid, the species, the eez and the time spent in the eez    
agg_data = CSV.read("data/agg_data.csv", DataFrame)

id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs];


rich, poor = Rich_Poor_lists(eezs, iso3_eez, economic_data)
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
function protection_odds(agg_data; start_protecting=[], species=nothing)
    eezs = unique(agg_data[:, :EEZ])
    ntot = length(eezs) - length(start_protecting)
    μ = zeros(ntot+1) # μ=μ(n)=<s^m>(n) = s^m * P^m(n) 
    σ2 = zeros(ntot+1) # σ^2=σ^2(n)= s^m * P^m(n) (1-P^m(n)
    
    #filter the species
    if species != nothing
        agg_data = agg_data[agg_data[:, :Species] .== species, :]
    end
    #protect the first eezs
    tot_ids = length(unique(agg_data[:, :newid]))
    agg_data = agg_data[agg_data[:, :EEZ] .∉ (start_protecting,), :]
    intially_protected = tot_ids - length(unique(agg_data[:, :newid]))
    # println("Initially protected: ", intially_protected)
    μ .+= intially_protected
    if size(agg_data, 1) == 0
        return μ, σ2
    end
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

function species_protection_odds(agg_data; threshold=0.5, start_protecting=[])
    eezs = unique(agg_data[:, :EEZ])
    
    species = unique(agg_data[:, :Species])
    species_threshold = [ceil(
        sum(
            unique(
                agg_data[:, [:Species, :newid]]
                )[:, :Species] .== s
            ) * threshold
            ) for s in species]
    prot_species = zeros(length(eezs) +1 - length(start_protecting))
    for (s, st) in zip(species, species_threshold)
        members_protected_mean, members_protected_var = protection_odds(
            agg_data,
            species=s,
            start_protecting=start_protecting)
        # println(size(prot_species))
        # println(members_protected_mean)
        prot_species .+= members_protected_mean .>= st
    end
    return prot_species    
end


##
protected_number_mean, protected_number_var = protection_odds(
    agg_data,
    start_protecting=rich)
protected_number_std = sqrt.(protected_number_var) 
top_bound = protected_number_mean .+ protected_number_std
bottom_bound = protected_number_mean .- protected_number_std

prot_species = species_protection_odds(
    agg_data,
    threshold=0.5,
    start_protecting=rich
    )

# save the protected number mean, variance and the species protected
CSV.write("percolation/uncorrelated_random/protected_rich.csv",
    DataFrame(
        protected_number_mean=protected_number_mean,
        protected_number_var=protected_number_var,
        prot_species=prot_species
        )
        )


##
p1 = plot(
    xlabel = "Cooperating EEZs",
    ylabel="Expected number of protected Individuals"
    )
plot!(p1, collect(0:length(protected_number_mean)-1), bottom_bound ./ N, fillrange = top_bound ./ N, label=L"\mu\pm \sigma^2", c=:blue, fillalpha=0.3, alpha=0)
plot!(p1, collect(0:length(protected_number_mean)-1), protected_number_mean ./ N, label=L"\mu", c=:black)

pspecies = plot(xlabel="cooperating EEZ", ylabel="Protected Species")
plot!(pspecies, collect(0:length(prot_species)-1), prot_species, label=nothing, c=:black)

ptot = plot(p1, pspecies, layout=(2,1), size=(800, 800))
##
savefig(ptot, "percolation/figures/uncorrelated_random.png")
savefig(ptot, "percolation/figures/uncorrelated_random.pdf")
plot(ptot)

##



##
protected_number_mean, protected_number_var = protection_odds(
    agg_data,
    )
protected_number_std = sqrt.(protected_number_var) 
top_bound = protected_number_mean .+ protected_number_std
bottom_bound = protected_number_mean .- protected_number_std

prot_species = species_protection_odds(
    agg_data,
    threshold=0.5,
    )

# save the protected number mean, variance and the species protected
CSV.write("percolation/uncorrelated_random/protected.csv",
    DataFrame(
        protected_number_mean=protected_number_mean,
        protected_number_var=protected_number_var,
        protected_species=prot_species
        )
        )

##
p1 = plot(
    xlabel = "Cooperating EEZs",
    ylabel="Expected number of protected Individuals"
    )
plot!(p1, collect(0:length(protected_number_mean)-1), bottom_bound ./ N, fillrange = top_bound ./ N, label=L"\mu\pm \sigma^2", c=:blue, fillalpha=0.3, alpha=0)
plot!(p1, collect(0:length(protected_number_mean)-1), protected_number_mean ./ N, label=L"\mu", c=:black)

pspecies = plot(xlabel="cooperating EEZ", ylabel="Protected Species")
plot!(pspecies, collect(0:length(prot_species)-1), prot_species, label=nothing, c=:black)

ptot = plot(p1, pspecies, layout=(2,1), size=(800, 800))



