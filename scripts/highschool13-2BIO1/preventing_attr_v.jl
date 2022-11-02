# Preventing polarization from forming as a function of number of attributes. 
# Figure 2
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework

using PolarizationFramework
using Graphs
using LinearAlgebra

# ns = [5, 9]#[9,15,25]
file = datadir("highschool13", "Highschool13-class2BIO1.lg")
gr = loadgraph(file)
ns = [nv(gr)]
A = adjacency_matrix(gr)
all_triads = get_triads(A)
all_links = get_links_in_triads(all_triads)
A = get_adj_necessary_links(size(A)[1], all_links; typ = Float64);

all_triads = get_triads(A);

#Heider72!
all_links = get_links_in_triads(all_triads)
triads_around_links_dict = get_triangles_around_links(all_triads)
counts = link_triangles_count(triads_around_links_dict; links = all_links)
triads_count_mat = PolarizationFramework.link_triangles_mat_inv(ns[1], all_links, counts)

reps = [100]
reps_dict = Dict(zip(ns, reps))
gs = [1, 3, 5, 7, 11]
# gs = [5]
threshold = 0.5;
vs = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1000];

attr_types = ["UA", "OA", "UPA"]

# sg = [-1.5, -1.1, -1.01, -0.99, -0.5, -0.1, -0.01]
# sg = [0, 0.01, 0.1, 0.5, 1.5]
# sg = [0.1, 0.3, 1.2]

gammas = [0.5, 1.5, 6]

all_params = @strdict(ns, gs, threshold, vs, attr_types)
dicts = dict_list(all_params)
[d["reps"] = reps_dict[d["ns"]] for d in dicts] #adding reps

for params in dicts
    n, g, attr_type, v, rep = let
        @unpack ns, threshold, reps, gs, attr_types, vs = params
        ns, gs, attr_types, vs, reps
    end

    println("Started Karate graph sims with g=$g and attr_type=", attr_type, " and v=$v.")
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
        "Heider722!";
        disp_each = 0,
        disp_more_every = 600,
        save_each = 600,
        files_folder = ["data", "highschool13-sims"],
        filename_prefix = "NumKarG",
        all_links_mat = A,
        all_triads = all_triads,
        triads_count_mat = triads_count_mat,
    )
end
