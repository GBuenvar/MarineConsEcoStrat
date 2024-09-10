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

mkpath("percolation/bipartite_rank")

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


#
# Consider the column newid as the nodes of the partition 1 and the column EEZ as the nodes of the partition 2
# The weight of the edges is the time spent in the EEZ as the colum "timestay (1/30days)"
# create a bipartite graph

"""
Rich protected first
"""
# create a DataFrame to store projected graph
##


"""
not include number of eez visited
not include time spent
"""

include_eezs_visited = false
weight_column = nothing
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    rich,
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_no_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_no_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_no_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_no_eezs_visited_unweighted_rich.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_no_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_no_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_species_times=prot_species_times))


p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)



##


"""
include number of eez visited
not include time spent
"""

include_eezs_visited = true
weight_column = nothing
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    rich,
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_eezs_visited_unweighted_rich.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_eezs_visited_unweighted_rich.csv.gz", DataFrame(protected_species_times=prot_species_times))



p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)

##

"""
not include number of eez visited
include time spent
"""


include_eezs_visited = false
weight_column = "timestay (1/30days)"
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    rich,
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_no_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_no_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_no_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_no_eezs_visited_weighted_rich.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_no_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_no_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_species_times=prot_species_times))


p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)

##

"""
include number of eez visited
include time spent
"""


include_eezs_visited = true
weight_column = "timestay (1/30days)"
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    rich,
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_eezs_visited_weighted_rich.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_eezs_visited_weighted_rich.csv.gz", DataFrame(protected_species_times=prot_species_times))


p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)




##
"""
Rich do not protect first
"""

##

"""
not include number of eez visited
include time spent
"""


include_eezs_visited = false
weight_column = nothing
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    [0, 8],
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_no_eezs_visited_unweighted.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_no_eezs_visited_unweighted.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_no_eezs_visited_unweighted.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_no_eezs_visited_unweighted.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_no_eezs_visited_unweighted.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_no_eezs_visited_unweighted.csv.gz", DataFrame(protected_species_times=prot_species_times))


p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)



##

include_eezs_visited = true
weight_column = nothing
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    [0, 8],
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_eezs_visited_unweighted.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_eezs_visited_unweighted.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_eezs_visited_unweighted.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_eezs_visited_unweighted.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_eezs_visited_unweighted.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_eezs_visited_unweighted.csv.gz", DataFrame(protected_species_times=prot_species_times))



p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)

##

include_eezs_visited = false
weight_column = "timestay (1/30days)"
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    [0, 8],
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_no_eezs_visited_weighted.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_no_eezs_visited_weighted.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_no_eezs_visited_weighted.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_no_eezs_visited_weighted.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_no_eezs_visited_weighted.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_no_eezs_visited_weighted.csv.gz", DataFrame(protected_species_times=prot_species_times))


p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)

##

include_eezs_visited = true
weight_column = "timestay (1/30days)"
@time prot_times, prot_number, prot_eezs, eez_prot_times = ranked_ids_remove_eezs(
    agg_data,
    [0, 8],
    include_eezs_visited,
    weight_column)
prot_species_number, prot_species_times = protected_species(
    prot_number,
    prot_times,
    id_to_species_int,
    newids
)


CSV.write("percolation/bipartite_rank/protected_times_eezs_visited_weighted.csv.gz", DataFrame(protected_times=prot_times))
CSV.write("percolation/bipartite_rank/protected_number_eezs_visited_weighted.csv.gz", DataFrame(protected_number=prot_number))
CSV.write("percolation/bipartite_rank/protected_eezs_eezs_visited_weighted.csv.gz", DataFrame(protected_eezs=prot_eezs))
CSV.write("percolation/bipartite_rank/eez_prot_times_eezs_visited_weighted.csv.gz", DataFrame(eez_prot_times=eez_prot_times))
CSV.write("percolation/bipartite_rank/protected_species_number_eezs_visited_weighted.csv.gz", DataFrame(protected_species_number=prot_species_number))
CSV.write("percolation/bipartite_rank/protected_species_times_eezs_visited_weighted.csv.gz", DataFrame(protected_species_times=prot_species_times))


p1 = plot(
    cumsum(prot_number)./N,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = plot(
    cumsum(prot_species_number)./N_species,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
# histogram(eez_prot_times, bins=30)
##
p1 = histogram(
    prot_times,
    bins=30,
    label="Number of individuals protected",
    xlabel="Step",
    ylabel="Number of individuals",
    legend=:topleft
    )
p2 = histogram(
    prot_species_times,
    bins=30,
    label="Number of species protected",
    xlabel="Step",
    ylabel="Number of species",
    legend=:topleft
    )
plot(p1, p2)
