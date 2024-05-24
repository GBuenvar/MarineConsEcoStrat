# Inspirado por las ideas de Pietronero Cristelli et al., M.A.Muñoz y Virginia Dominguez_García presentan el algoritom MusRank.

# Este algoritmo sirve para la ordenacion de redes mutualistas de la manera mas nested posiible.
# En Redes mutualistas, en las que se relacionan dos conjuntos de nodos, como pueden ser bienes exportados por paises, semillas
# dispersadas por pajaros, etc., se pueden encontrar estructuras nested en las que los nodos especialistas de un conjuntos
# se relacionan con los nodos generalistas del otro conjunto. En tal caso, definen uno de los conjuntos como el 
# "activo" (países, dospersadores, peces) y otro como el "pasivo" (bienes, semillas, anemonas). El rol de cada conjunto
# no es una mera cuestion de nomenclatura, sino que caracteriza la estructura de la red. 
# Por ejemplo, en el caso de la red de exportaciones de productos por paises, se observa una ordenacion en la que los productos
# que son exportados por pocos paises, son exportados por los paises que exportan muchos productos, mientras que los paises que 
# exportan pocos productos exportan productos que son exportados por muchos paises. No hay especializacion en la red

# - los paises generalistas exportan todo tipo de productos
# - los productos generalistas son exportados por todo tipo de paises
# - los productos que son exportados por pocos paises son exportados por los paises generalistas
# - los paises que exportan pocos productos exportan productos generalistas

# Este tipo de diseño confiere robusted a la red frente a perdidas de nodos, ya que los nodos generalistas pueden ser sustituidos
# por otros nodos generalistas, mientras que los nodos especialistas no pueden ser sustituidos por otros nodos especialistas.
using CSV, DataFrames, Random, Statistics, Plots, XLSX

##
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

function animal_eez_matrix(agg_data, weighted_col=nothing)
    n_animals = length(unique(agg_data[:, :newid]))
    n_eez = length(unique(agg_data[:, :EEZ]))
    matrix = zeros(n_animals, n_eez)
    eezs = unique(agg_data[:, :EEZ])
    newdis = unique(agg_data[:, :newid])
    for row in eachrow(agg_data)
        if isnothing(weighted_col)
            matrix[newdis .== row.newid, eezs .== row.EEZ] .= 1
        else
            matrix[newdis .== row.newid, eezs .== row.EEZ] .= row[weighted_col]
        end
    end
    return matrix
end

function species_eez_matrx(agg_data, species_codes, weighted_col=nothing)
    n_species = length(unique(agg_data[:, :Species]))
    n_eez = length(unique(agg_data[:, :EEZ]))
    matrix = zeros(n_species, n_eez)
    eezs = unique(agg_data[:, :EEZ])
    species = unique(agg_data[:, :Species])

    for row in eachrow(agg_data)
        if isnothing(weighted_col)
            matrix[species .== row.Species, eezs .== row.EEZ] .+= 1
        else
            matrix[species .== row.Species, eezs .== row.EEZ] .+= row[weighted_col]
        end
    end
    matrix = matrix ./ sum(matrix, dims=2)
end

##
function importance_n(matrix, prev_vulnerability)
    importance = sum(matrix .* prev_vulnerability, dims=1)
    importance = importance ./ sum(importance)
    return importance
end

function vulnerability_n(matrix, prev_importance)
    vulnerabilty = 1 ./ sum(matrix .* (1 ./ prev_importance)', dims=1)
    vulnerabilty = vulnerabilty ./ sum(vulnerabilty)
    return vulnerabilty
end

function fixed_point(matrix, max_iter=1000, tol=1e-6)
    importance = size(matrix, 2)
    vulnerability = size(matrix, 1)
    for i in 1:max_iter
        new_importance = importance_n(matrix, vulnerability)
        new_vulnerability = vulnerability_n(matrix, importance)
        if norm(new_importance - importance) < tol && norm(new_vulnerability - vulnerability) < tol
            break
        end
        importance = new_importance
        vulnerability = new_vulnerability
    end
    return importance, vulnerability
end


animal_eez = animal_eez_matrix(agg_data)
I,V = fixed_point(animal_eez)
