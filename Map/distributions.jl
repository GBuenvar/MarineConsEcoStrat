# %%
using CSV, DataFrames, Plots, LaTeXStrings

# Read in the data

trajectories = CSV.read("data/full_data_inds.csv.gz", DataFrame)
# Read the codes dictionaries of the species and the eezs
eez_codes = CSV.read("data/eez_to_int.csv", DataFrame)
int_to_eez = Dict(zip(eez_codes.Int, eez_codes.EEZ))
species_codes = CSV.read("data/species_to_int.csv", DataFrame)
int_to_species = Dict(zip(species_codes.Int, species_codes.Species))



agg_data = unique(trajectories[:, [:newid, :Species, :EEZ]])


unique_species = unique(agg_data.Species)
unique_eez = unique(agg_data.EEZ)


counts_sp_eez = combine(groupby(agg_data, [:Species, :EEZ]), nrow)


# %%
function hill_number(p, q)
    @assert abs(sum(p) - 1 < 10^-6) "probabilities are not normalized"
    @assert isinteger(q) "q must be a float type integer"
    p = p[p .> 0.]
    q == 1. ? (h = exp(-sum(p .* log.(p)))) : (h = sum(p.^q)^(1/(1-q)))
    return h
end

# %%
bins = 0:2:80
p = plot(xlabel="Species", ylabel ="% of EEZs")
for (q, st) in zip([0, 1], ["Number", "Effective number (q=1)"])
    div_q = zeros(Float64, length(unique_eez))
    for (i, eez) in enumerate(unique_eez)
        sizes = Vector{Float64}(counts_sp_eez[counts_sp_eez.EEZ .== eez, :nrow])
        sizes ./= sum(sizes)
        div_q[i] = hill_number(sizes, q)
    end
    histogram!(p, div_q, alpha=0.5, label= "$st", bins=bins, normalize=:probability)
end
yaxis!(p, yticks =(0:0.20:0.60, [0,20,40]))

p

# %%

bins = 0:2:60
p = plot(xlabel="EEZs", ylabel ="% of Species")
for (q, st) in zip([0, 1], ["Visited", "Effective visited (q=1)"])
    div_q = zeros(Float64, length(unique_species))
    for (i, sp) in enumerate(unique_species)
        sizes = Vector{Float64}(counts_sp_eez[counts_sp_eez.Species .== sp, :nrow])
        sizes ./= sum(sizes)
        div_q[i] = hill_number(sizes, q)
    end
    histogram!(p, div_q, alpha=0.5, label= "$st", normalize=:probability, bins=bins)
    

end
yaxis!(p, yticks =(0:0.20:0.60, [0,20,40]))
xaxis!(p, xticks = (0:20:60))
p