# file contains some utility functions. 

using Dates
using Graphs

# module MyUtil


# finds the first free filename (that would not be overwriting) with
# the expression: filename * "_$i", where i is integer starting from 1.
# adds extension at the end.
# returns proper filename
function get_filename(filename, extension)
    i = 1
    #if extension does not have a "." at the beginning we should add one.
    if !startswith(extension, ".")
        extension = "." * extension
    end

    #filename should not end with the extension.
    if endswith(filename, extension)
        filename = filename[1:(end-length(extension))]
    end

    while isfile(filename * "_$i" * extension)
        i += 1
    end
    filename * "_$i" * extension
end
export get_filename

#get k max values from array a
function maxk(a, k)
    b = partialsortperm(a, 1:k, rev = true)
    return collect(zip(b, a[b]))
end
export maxk

macro Name(arg)
    string(arg)
end

#function for displaying progress
function display_res(names, vars)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    str = ""
    for i = 1:1:length(names)
        str = str * "$(names[i])" * " = " * "$(vars[i])" * "; "
    end
    println("Just finished: $str" * "Time: $cur")
    return nothing
end

function display_res(ntuple::NamedTuple)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    str = savename(ntuple, connector = "; ", sort = false)

    println("Just finished: $str" * ". Time: $cur")
    return nothing
end

function display_res(; kwargs...)
    ntuple = NamedTuple(kwargs)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    str = string(ntuple)

    println("Just finished: $str. Time: $cur")
    return nothing
end

function display_res(gamma::Float64, z::Int)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    println("Just finished: z = $z; gamma = $gamma; " * "Time: $cur")
    return nothing
end

# Function for in-place creating a symmetrical array of array x. 
# x should be an array of size `nxn`. 
function Sym!(x_sim, x, n)
    for i = 1:n
        x_sim[i, i] = x[i, i]
        for j = (i+1):n
            x_sim[i, j] = x[i, j]
            x_sim[j, i] = x[i, j]
        end
    end
end
export Sym!

#returns list of triads from adjacency matrix A
function get_triads(A)
    N = size(A)[1]

    triads = []

    for i in 1:N
        for j in (i+1):N
            if A[i,j] > 0
                for k in (j+1):N
                    if A[i,k] > 0 && A[j,k] > 0
                        push!(triads, (i,j,k))
                    end
                end
            end
        end
    end
    
    return triads
end
export get_triads

#returns all links that belong to triads
function get_links_in_triads(all_triads)
    links = []
    for triad in all_triads
        i, j, k = triad
        push!(links, (i,j))
        push!(links, (i,k))
        push!(links, (j,k))
    end
    links = unique(links)
    # sort according to y and then x. 
    sort!(links,
            lt = (x, y) -> (x[2] < y[2] || (x[2] == y[2] && x[1] < y[1])))
    return links
end
export get_links_in_triads

#creates an asjacency matrix based on given links and number of nodes. 
#It is used to create adjacency matrix with links forming triads. 
function get_adj_necessary_links(n, links; typ = Int)
    A = zeros(typ,n, n)

    for link in links
        ind1, ind2 = link
        A[ind1, ind2] = 1
    end

    return A
end
export get_adj_necessary_links

function get_edgelist_it(g::Graph)
    return edges(g)
end

function get_edgelist_it(A::Matrix)
    g = Graph(A)
    return edges(g)
end

# see, wikipedia line graph
# works for simle graphs (not directed)
function create_line_graph(g::Graph)
    # number of nodes in LG
    ln = g.ne 
    lg = Graph()
    add_vertices!(lg, ln)

    edgelist = collect(get_edgelist_it(g))
    dict_e = Dict(edgelist .=> 1:ln)

    for (ei, edge) in enumerate(edgelist)
        a1, a2 = edge.src, edge.dst

        # println(collect(edges(g)))
        for ai in [a1, a2]
            nais = deepcopy(neighbors(g, ai))
            # println(collect(edges(g)))
            filter!(x-> !(x in [a1, a2]), nais)
            # println(collect(edges(g)))
            #get edge numbers
            neis = [dict_e[Edge(min(ai, nai), max(ai, nai))] for nai in nais]
            # println(collect(edges(g)))
            #add edges
            map(x -> add_edge!(lg, x, ei), neis)
            # println(collect(edges(g)))
            # println(lg)
        end

    end
    return lg
end

# Returns a dictionary with links as keys (Tuple{Int64,Int64})
# and list of triads. Each triad contains a list of other links. 
# So the values of the dictionary are of type Vector{Vector{Tuple{Int64, Int64}}}
function get_triangles_around_links(g::Graph)
    n = nv(g)

    # edgelist = get_edgelist_it(g)
    # edgelist2 = [(edge.src, edge.dst) for edge in edgelist]
    dict_e = Dict{Tuple{Int64, Int64}, Vector{Vector{Tuple{Int64, Int64}}}}()

    for agent in 1:n
        neighs = neighbors(g, agent)

        for neigh in neighs
            if neigh < agent
                continue
            end
            neighs2 = neighbors(g, neigh)

            com_neighs = intersect(neighs, neighs2)

            for com_neigh in com_neighs
                if com_neigh < agent || com_neigh < neigh
                    continue
                end

                if !haskey(dict_e, (agent, neigh))
                    dict_e[(agent, neigh)] = [[(agent, com_neigh), (neigh, com_neigh)]]
                else
                    append!(dict_e[(agent, neigh)], [[(agent, com_neigh), (neigh, com_neigh)]])
                end
                if !haskey(dict_e, (agent, com_neigh))
                    dict_e[(agent, com_neigh)] = [[(agent, neigh), (neigh, com_neigh)]]
                else
                    append!(dict_e[(agent, com_neigh)], [[(agent, neigh), (neigh, com_neigh)]])
                end
                if !haskey(dict_e, (neigh, com_neigh))
                    dict_e[(neigh, com_neigh)] = [[(agent, neigh), (agent, com_neigh)]]
                else
                    append!(dict_e[(neigh, com_neigh)], [[(agent, neigh), (agent, com_neigh)]])
                end
            end
        end
    end
    return dict_e
end

function get_triangles_around_links(triads::Vector{Any})
    dict_e = Dict{Tuple{Int64, Int64}, Vector{Vector{Tuple{Int64, Int64}}}}()

    for (ti, triad) in enumerate(triads)
        agent, neigh, com_neigh = triad
        # tlinks = [(ai, aj), (ai, ak), (aj, ak)]

        if !haskey(dict_e, (agent, neigh))
            dict_e[(agent, neigh)] = [[(agent, com_neigh), (neigh, com_neigh)]]
        else
            append!(dict_e[(agent, neigh)], [[(agent, com_neigh), (neigh, com_neigh)]])
        end
        if !haskey(dict_e, (agent, com_neigh))
            dict_e[(agent, com_neigh)] = [[(agent, neigh), (neigh, com_neigh)]]
        else
            append!(dict_e[(agent, com_neigh)], [[(agent, neigh), (neigh, com_neigh)]])
        end
        if !haskey(dict_e, (neigh, com_neigh))
            dict_e[(neigh, com_neigh)] = [[(agent, neigh), (agent, com_neigh)]]
        else
            append!(dict_e[(neigh, com_neigh)], [[(agent, neigh), (agent, com_neigh)]])
        end
    end

    return dict_e
end
export get_triangles_around_links

# returns list whose elements correspond to links array
function get_triangles_around_links(dict_e::Dict, links::Vector)
    nl = length(links)

    dict_l = Dict(links .=> 1:nl)

    return [[(dict_l[pair[1]], dict_l[pair[2]]) for pair in dict_e[link]] 
        for link in links]
end

function link_triangles_count(dict_e::Dict; links = [])
    if isempty(links)
        return [length(value) for (key, value) in dict_e]
    else
        counts = zeros(size(links))
        for (i, link) in enumerate(links)
            if haskey(dict_e, link)
                counts[i] = length(dict_e[link])
            end
        end
        return counts
    end
end
export link_triangles_count

# creates a matrix with values indicating number of triads 
# a link (i,j) is involved in. 
# Output matrix has nonzero values only in upper matrix triangle. 
# To have the proper matrix run `Symmetric` afterwards.
function link_triangles_mat(n, links, counts)
    mat = zeros(n, n)

    map(i -> mat[links[i]...] = counts[i], 1:length(links))
    return mat
end

# Same as `link_triangles_mat` but it contains zeros and 
# inverse of triad counts. 
function link_triangles_mat_inv(n, links, counts)
    mat = zeros(n, n)

    map(i -> mat[links[i]...] = 1/counts[i], 1:length(links))
    return mat
end