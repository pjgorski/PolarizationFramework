# Windsurfer interactions dataset topology
# Preventing polarization from forming as a function of number of attributes. 
# Figure 1
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework
using Graphs, GraphIO
using LinearAlgebra

# ns = [5, 9]#[9,15,25]
file = datadir("windsurfers-interactions", "out.moreno_beach_beach")
gr = loadgraph(file, "graph_key", NETFormat())

ns = [nv(gr)]
A = adjacency_matrix(gr)
all_triads = get_triads(A)
all_links = get_links_in_triads(all_triads)
A = get_adj_necessary_links(size(A)[1], all_links; typ = Float64);

# Heider9!
link_indices = findall(triu(A, 1)[:] .> 0)
triads_around_links_dict = get_triangles_around_links(all_triads)
link_pairs = get_triangles_around_links(triads_around_links_dict, all_links)
link_pairs_triad_cnt = [length(link) for link in link_pairs];

reps = [100]
reps_dict = Dict(zip(ns, reps))
# two parts of simulations
gs = [1:2:21...]
# gs = [25, 29, 33, 37, 41, 45, 49, 55, 61, 67, 73, 81, 89, 97]
threshold = 0.5;
vs = [4, @onlyif("attr_types" == "OA", 1000)] #includes CA

attr_types = ["BA", "UA", "OA", "UPA"]
# attr_types = ["BA", "OA", "UPA"]
# attr_types = ["OA"]

gammas = [0.5, 1.5, 4]

all_params = @strdict(ns, gs, threshold, vs, attr_types)
dicts = dict_list(all_params)
[d["reps"] = reps_dict[d["ns"]] for d in dicts] #adding reps
# [d["gammas"] = gammas for d in dicts] #adding gammas

for params in dicts
    n, g, attr_type, v, rep = let
        @unpack ns, threshold, reps, gs, attr_types, vs = params
        ns, gs, attr_types, vs, reps
    end

    println("Started n=$n and g=$g and attr_type=", attr_type, " and v=$v.")

    if attr_type == "UA"
        attr = UnorderedAttributes(g, threshold, v)
    elseif attr_type == "BA"
        attr = BinaryAttributes(g)
    elseif attr_type == "OA"
        attr = OrderedAttributes(g, threshold, v)
    elseif attr_type == "UPA"
        attr = UnorderedPositiveAttributes(g, threshold, v)
    else
        throw(attr_type)
    end

    r = using_heider_attr(
        n,
        attr,
        gammas,
        rep,
        3000.0,
        "Heider92!";
        disp_each = 0,
        disp_more_every = 600,
        save_each = 600,
        files_folder = ["data", "windsurfers-interactions-sims"],
        filename_prefix = "NumKarG",
        all_links_mat = A,
        all_triads = all_triads,
        link_indices = link_indices,
        link_pairs = link_pairs,
        link_pairs_triad_cnt = link_pairs_triad_cnt,
    )
end
