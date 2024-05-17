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
        for row in eachrow(neighbors_times)
            push!(projected_df, [node, row.newid, row.nrow])
            push!(projected_df, [row.newid, node, row.nrow])
        end
        # remove the node from the bipartite graph
        bipartite_df = bipartite_df[bipartite_df[:, projected_partition] .!= node, :]
    end
    return projected_df
    
end

@time projected_df = projected_graph(agg_data, :newid, :EEZ)

##

function degree_ranking(projected_df, node_N=nothing)
    # calculate the degree of each node
    # by summing all the weights of the same node
    degree_df = combine(
        groupby(projected_df, :node_A),
        :weight => sum
    )
    rename!(degree_df, :weight_sum => :degree)
    if node_N != nothing
        degree_df[!, :degree] = degree_df[!, :degree] ./ node_N
    end
    # sort the degree_df
    sort!(degree_df, :degree, rev=true)
    return degree_df
    
end

@time degree_df = degree_ranking(projected_df)

nodes_N = combine(groupby(agg_data, :newid), nrow)[:, :nrow]
@time degree_df_N = degree_ranking(projected_df, nodes_N)


function ranked_ids_remove_eezs(data, start_protecting = [0,8], include_eez_resistant=false)
    ids::Vector{Int64} = unique(data[:, :newid])
    eezs::Vector{Int64} = unique(data[:, :EEZ])
    iterated_eezs::Vector{Int64} = setdiff(eezs, start_protecting)
    Neez::Int = length(iterated_eezs)

    # Initialize the variables
    unprotected_ids  = unique(data[:, :newid])
    prot_times  = fill(Neez+1, length(ids)) 
    prot_number = zeros(Int64, Neez+1)
    prot_eezs   = copy(start_protecting)
    

    #######
    # Ver bien como implementar esto, porque en lugar de incorporar
    # las EEZ de una en una como en otros metodos,
    #  en este caso se incorportan todas las que visita el primero del 
    #  ranking a a la vez. 
    #  una opcion es guardar en prot_times el numero de EEZ que 
    #  se han protegido antes de proteger al individuo. En prot_number
    #  que haya saltos que indiquen que se han protegido varias a la vez, \
    #  y a la vez que prot_eezs se guarde en que paso se ha protegido cada eez
    #  ya que todas las que se protegen a la vez tienen el mismo tiempo.
    #######


    # compute the degree of each node and rank them, then remove its EEZs
    if include_eez_resistant
        higher_rank = degree_ranking(
            projected_graph(data, :newid, :EEZ),
            combine(groupby(data[:, :newid], :newid), nrow)
            )[1]
    else
        degree_df = degree_ranking(
            projected_graph(data, :newid, :EEZ)
            )[1]
    end

    