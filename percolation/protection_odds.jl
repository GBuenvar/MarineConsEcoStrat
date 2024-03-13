using CSV, DataFrames, Random, Statistics, Plots, XLSX, ArgParse
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


function ind_protection_odds(N, m)
    # the number of ways to choose the m visited eezs
    # out of the n picked eezs divided by the
    # number of ways to choose the m visited eezs
    # out of the N total eezs
    P = [n < m ? 0.0 :
    binomial(big(n), big(m)) / binomial(big(N), big(m))
    for n in 1:N]
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




function protection_odds(agg_data)
    # newids = unique(agg_data[:, :newid])
    ntot = length(eezs)
    μ = zeros(ntot)
    σ2 = zeros(ntot)
    # P_ids = zeros(length(newids), length(eezs))
    m_of_ids = combine(groupby(agg_data, :newid), nrow)[!, :nrow]
    for m in range(1, maximum(m_of_ids))
        ids_with_m = m_of_ids .== m
        Size_m = sum(ids_with_m)
        # P^m(n)
        P_m_at_n = ind_protection_odds(ntot, m) 
        # println(size(repeat(P_m, number_m)))
        # println(size(P_ids[idx, 2:end]))
        # store P_m in all the rows of P_ids[idx, 2:end]
        # P_ids[idx, 2:end] .= repeat(P_m, number_m)

        # first moment
        μ .+= Size_m .* P_m_at_n
        # second moment
        σ2 .+= (Size_m .^ 2.) .* P_m_at_n
    end
    σ2 .-= μ .^ 2.
    return μ, σ2
end

function proxy_protection_odds(agg_data)
    # newids = unique(agg_data[:, :newid])
    ntot = length(eezs)
    μ = zeros(ntot)
    σ2 = zeros(ntot)
    # P_ids = zeros(length(newids), length(eezs))
    m_of_ids = combine(groupby(agg_data, :newid), nrow)[!, :nrow]
    for m in range(1, maximum(m_of_ids))
        ids_with_m = m_of_ids .== m
        Size_m = sum(ids_with_m)
        # P^m(n)
        P_m_at_n = proxy_odds(ntot, m) 
        # println(size(repeat(P_m, number_m)))
        # println(size(P_ids[idx, 2:end]))
        # store P_m in all the rows of P_ids[idx, 2:end]
        # P_ids[idx, 2:end] .= repeat(P_m, number_m)

        # first moment
        μ .+= Size_m .* P_m_at_n
        # second moment
        σ2 .+= (Size_m .^ 2.) .* P_m_at_n
    end
    σ2 .-= μ .^ 2.
    return μ, σ2
end

eezsvisited = combine(groupby(agg_data, :newid), nrow)[!, :nrow]
##
p1 = plot()
p2 = plot()

plot!(p1, sum([ind_protection_odds(180, n) .* sum(eezsvisited .== n) for n in 1:180]), label="<inds>")
plot!(p1, sum([proxy_odds(180, n) .* sum(eezsvisited .== n) for n in 1:180]), label="<proxy>", c=:red)

plot(p1)

plot!(p2, [ind_protection_odds(180, n)  for n in 1:180], label=nothing, c=:blue)
plot!(p2, [proxy_odds(180, n) for n in 1:180], label=nothing, c=:red)

plot!(p2, [[],[]], label=["P^m(n)" "P^m(n) proxy"], c=[:blue :red])

plot(p1, p2)