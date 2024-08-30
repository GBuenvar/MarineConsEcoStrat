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

mkpath("percolation/MusRank")


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

function interaction_matrix(agg_data, active_col, passive_col, weighted_col=nothing)
    n_pas = length(unique(agg_data[:, passive_col]))
    n_ac = length(unique(agg_data[:, active_col]))
    matrix = zeros(n_pas, n_ac)
    actives = unique(agg_data[:, active_col])
    pasives = unique(agg_data[:, passive_col])
    for row in eachrow(agg_data)
        if isnothing(weighted_col)
            matrix[pasives .== row[passive_col], actives .== row[active_col]] .= 1
        else
            matrix[pasives .== row[passive_col], actives .== row[active_col]] .= row[weighted_col]
        end
    end
    return matrix
end

"""
species_eez_matrix(agg_data, weighted_col=nothing, intensive=false)

Returns the species-EEZ matrix, where the rows are the species and the columns are the EEZs. The matrix is normalized by the total number of interactions of each species.

Parameters:
    - agg_data: DataFrame with the data
    - weighted_col: Column to weight the interactions. If nothing, each interaction is counted as 1
    - intensive: Whether to normalize the total number of interactions of 
"""
function species_eez_matrix(agg_data; weighted_col=nothing, intensive=false)
    # intensive := binary
    n_species = length(unique(agg_data[:, :Species]))
    n_eez = length(unique(agg_data[:, :EEZ]))
    matrix = zeros(n_species, n_eez)
    eezs = unique(agg_data[:, :EEZ])
    species = unique(agg_data[:, :Species])

    # for row in eachrow(agg_data)
    #     if isnothing(weighted_col)
    #         matrix[species .== row.Species, eezs .== row.EEZ] .+= 1
    #     else
    #         matrix[species .== row.Species, eezs .== row.EEZ] .+= row[weighted_col]
    #     end
    # end
    # matrix = matrix ./ sum(matrix, dims=2)

    # desplegar las cuatro opciones:

    # weighted col se usa para el valor del flujo. Si no se especifica, cada interaccion tiene un valor de 1,
    # y global_weight es el numero de interacciones. sp_total es el numero/pesp de interacciones de la especie
    # sp_eez es el numero/peso de interacciones de la especie en la eez.
    # M_cp es la fraccion de interacciones de la especie en la eez. Si intensive es true, se normaliza por el numero
    


    global_weight = nothing
    if intensive
        global_weight = isnothing(weighted_col) ? size(agg_data, 1) : sum(agg_data[!, weighted_col])
    end
    println("Global weight: $global_weight")
    @views for sp in species
        sp_total = isnothing(weighted_col) ? sum(agg_data.Species .== sp) : sum(agg_data[agg_data.Species .== sp, weighted_col])
        @views for eez in eezs
            cols = (agg_data.Species .== sp) .& (agg_data.EEZ .== eez)
            sp_eez = isnothing(weighted_col) ? sum(cols) : sum(agg_data[cols, weighted_col])
            M_cp = sp_eez / sp_total
            if intensive
                eez_total = isnothing(weighted_col) ? sum(agg_data.EEZ .== eez) : sum(agg_data[agg_data.EEZ .== eez, weighted_col])
                M_cp *= eez_total / global_weight # esto estaba al reves antes, no se cual es el que esta bien
            end
            matrix[species .== sp, eezs .== eez] .= M_cp
        end
    end
    return matrix
end

##

# Es posible que el error este en la normaliizacon
function importance_n(matrix, prev_vulnerability)
    importance = sum(matrix .* prev_vulnerability, dims=1)[1, :]
    importance = importance ./ mean(importance)
    return importance
end

function vulnerability_n(matrix, prev_importance)
    vulnerabilty = 1 ./ sum(matrix .* (1 ./ prev_importance)', dims=2)[:, 1]
    vulnerabilty = vulnerabilty ./ mean(vulnerabilty)
    return vulnerabilty
end

function fixed_point(matrix, max_iter=10000, tol=1e-6)
    vulnerability = ones(size(matrix, 1))
    importance = ones(size(matrix, 2))
    last_i = 1
    for i in 1:max_iter
        new_importance = importance_n(matrix, vulnerability)
        new_vulnerability = vulnerability_n(matrix, importance)
        if maximum(abs.(new_importance .- importance)) < tol && maximum(abs.(new_vulnerability .- vulnerability)) < tol
            last_i = i
            break
        end
        importance = new_importance
        vulnerability = new_vulnerability
    end
    last_i = max_iter
    println("Fixed point reached in $last_i iterations")
    return importance, vulnerability
end


function fixed_point_sort(matrix, I, V, act_list, pas_list)
    act_list = act_list[sortperm(I, rev=true)]
    pas_list = pas_list[sortperm(V, rev=true)]
    matrix = matrix[:, sortperm(I, rev=true)]
    matrix = matrix[sortperm(V, rev=true), :]
    I = I[sortperm(I, rev=true)]
    V = V[sortperm(V, rev=true)]
    return matrix, I, V, act_list, pas_list
    
end

animal_eez = interaction_matrix(agg_data, :EEZ, :newid)
I, V = fixed_point(animal_eez)
sorted_animal_eez, sorted_I, sorted_V, sorted_countries_animals, sorted_animals = fixed_point_sort(animal_eez, I, V, [int_to_eez[eez] for eez in eezs], newids)


##
heatmap(sorted_animal_eez',
        xflip=true,
        yflip=true, 
        yticks=(1:length(sorted_countries_animals), sorted_countries_animals),
        c=cgrad([:white, :orange]),
        xlabel="EEZ",
        ylabel="Animal",
        size=(800, 800)
        )

##
ps = []
importances = []
vulnerabilities = []
Countries_sorted = []
Species_sorted = []
for weighted_col in [nothing, "timestay (1/30days)"]
    for intensive in [false, true] 
    species_eez = species_eez_matrix(agg_data, weighted_col, intensive)
    Isp, Vsp = fixed_point(species_eez)
    sorted_species_eez,sorted_I, sorted_V, sorted_countries_species, sorted_species = fixed_point_sort(species_eez, Isp, Vsp, [int_to_eez[eez] for eez in eezs], [int_to_species[sp] for sp in unique(agg_data.Species)])
    push!(importances, sorted_I)
    push!(vulnerabilities, sorted_V)
    push!(Countries_sorted, sorted_countries_species)
    push!(Species_sorted, sorted_species)
    p = heatmap(sorted_species_eez',
            xflip=true,
            yflip=true, 
            xticks=(1:length(sorted_species), sorted_species),
            xtickfontsize=2,
            yticks=(1:length(sorted_countries_species), sorted_countries_species),
            xrotation=90,
            ytickfontsize=2,
            c=cgrad([:white, :orange]),
            ylabel="EEZ",
            xlabel="Species",
            # legendfontsize = 200,
            # xguidefontsize = 14, 
            # yguidefontsize = 14,
            # size=(800, 800),
            # fontsize = 5,
            bottom_margin=5Plots.mm,
            dpi=300)
    w = isnothing(weighted_col) ? "appearence" : "timestay"
    i = intensive ? "intensive" : "extensive"
    CSV.write("percolation/MusRank/species_eez_matrix_$(w)_$i.csv", DataFrame(sorted_species_eez, :auto))
    CSV.write("percolation/MusRank/sorted_species_$(w)_$i.csv", DataFrame(sorted_species=sorted_species, vulnerabilities=sorted_V))
    CSV.write("percolation/MusRank/sorted_countries_species_$(w)_$i.csv", DataFrame(sorted_countries_species=sorted_countries_species, importances=sorted_I))
    push!(ps, p)
    title!("Species-EEZ matrix$(w), $(i)")

    savefig("Nestedness/species_eez_matrix_$(w)_$(i).png")
    # remove title
    title!("$(w), $(i)")
    end
end
p_tot = plot(ps[1], ps[2], ps[3], ps[4], layout=(2, 2), size=(800, 800), dpi=300)
# write a 
savefig(p_tot, "Nestedness/species_eez_matrix_all.png")
p_tot

plot(vulnerabilities, labels = ["binary extensive" "binary intensive" "weighted extensive" "weighted intensive"])