# La idea es proyectar la red sobre animales (o especeies??). Tomar el degree de cada nodo
# y Rankearlos acorde a ello. Del que esta en el top, proteger todas las zonas que visita. 
# Recalcular el degree de cada nodo y volver a rankearlos. Repetir el proceso hasta que todos
# los nodos esten protegidos.

# Alternativamente se puede hacer lo mismo dividiendo el degree entre el numero de zonas que visita
# Para tener en cuenta que proteger muchas zonas es mas dificil que proteger pocas.

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

mkpath("percolation/individuals")

##
# Read the data, consisting on the newid, the species, the eez and the time spent in the eez    
agg_data = CSV.read("data/agg_data.csv", DataFrame)

id_to_species_int = Dict(zip(agg_data.newid, agg_data.Species))
newids = unique(agg_data[:, :newid])
N = length(newids)
N_species = length(unique(agg_data[:, :Species]))

eezs = unique(agg_data[:, :EEZ])
iso3_eez = [eez_to_iso3[int_to_eez[eez]] for eez in eezs];

#
# Consider the column newid as the nodes of the partition 1 and the column EEZ as the nodes of the partition 2
# The weight of the edges is the time spent in the EEZ as the colum "timestay (1/30days)"
# create a bipartite graph

# create a DataFrame to store projected graph
##

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

    for node in unique(bipartite_df[:, projected_partition])
        node_EEZs = bipartite_df[bipartite_df[:, projected_partition] .== node, complementary_partition]
        neighbouring_subset = bipartite_df[bipartite_df[:, complementary_partition] .∈ (node_EEZs,), :]
        # count the number of times each newid appears in the subset
        neighbors_times = combine(
            groupby(neighbouring_subset, :newid),
            nrow
            )
        # remove the rows of node in the neighbors_times
        neighbors_times = neighbors_times[neighbors_times[:, :newid] .!= node, :]
        # add the node to the projected graph
        node_A = fill(node, size(neighbors_times, 1))
        node_B = neighbors_times[:, :newid]
        weight = neighbors_times[:, :nrow]
        
        projected_df = vcat(
            projected_df,
            DataFrame(
                node_A = node_A,
                node_B = node_B,
                weight = weight
            )
        )

    end
    return projected_df
end


@time projected_graph(agg_data, :newid, :EEZ)


## 
using Graphs
function degree_ranking(
        projected_df,
        node_weight,
)

    # calculate the degree of each node
    degree_df = combine(
        groupby(projected_df, :node_A),
        :weight => sum
    )
    rename!(degree_df, :weight_sum => :degree)
    # sort the nodes by degree
    sort!(degree_df, :degree, rev=true
    
end

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
for node in nodes

    node_EEZ_list = bipartite_df[
        bipartite_df[:, projected_partition] .== node, complementary_partition
        ]
    neighbouring_subset = @view bipartite_df[
        bipartite_df[:, complementary_partition] .∈ (node_EEZ_list,), :]
    # remove the node from the subset
    neighbouring_subset = neighbouring_subset[
        neighbouring_subset[:, projected_partition] .!= node, :]
    # count the number of times each newid appears in the subset
    neighbors_times = combine(
        groupby(neighbouring_subset, :newid),
        nrow
        )
    # add the node to the projected graph
    node_A = fill(node, size(neighbors_times, 1))
    node_B = neighbors_times[:, :newid]
    weight = neighbors_times[:, :nrow]
    
    projected_df = vcat(
        projected_df,
        DataFrame(
            node_A = node_A,
            node_B = node_B,
            weight = weight
        )
    )

    projected_df = vcat(
        projected_df,
        DataFrame(
            node_A = node_B,
            node_B = node_A,
            weight = weight
        )
    )
    # remove the node from the bipartite graph
    bipartite_df = bipartite_df[bipartite_df[:, projected_partition] .!= node, :]
end
# sort the projected graph by node_A, node_B
#sort!(projected_df, [:node_A, :node_B])

return projected_df
end


