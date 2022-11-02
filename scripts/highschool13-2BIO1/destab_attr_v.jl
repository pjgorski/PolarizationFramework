# Highschool13 2BIO1 topology
# Destabilizing polarization as a function of number of attributes. 
# Figure 2
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework
using Graphs

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

gammas = [0.5, 1.5, 6]

larger_sizes = [@onlyif(i <= "ns" < 2 * i, i) for i = 1:maximum(ns)]

all_params = @strdict(ns, gs, threshold, vs, attr_types, larger_sizes)
dicts = dict_list(all_params)
[d["reps"] = reps_dict[d["ns"]] for d in dicts] #adding reps
# [d["gammas"] = gammas for d in dicts] #adding gammas

for params in dicts
    n, g, attr_type, v, rep, larger_size = let
        @unpack ns, threshold, reps, gs, attr_types, vs, larger_sizes = params
        ns, gs, attr_types, vs, reps, larger_sizes
    end
    println("Started H13", @ntuple(n, g, attr_type, v, larger_size))

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


    r = using_heider_attr_destab(
        n,
        attr,
        gammas,
        larger_size,
        rep,
        3000.0,
        "Heider722!",
        disp_each = 0,
        disp_more_every = 600,
        save_each = 600,
        files_folder = ["data", "highschool13-sims"],
        filename_prefix = "DesKarG",
        all_links_mat = A,
        all_triads = all_triads,
        triads_count_mat = triads_count_mat,
    )
end
