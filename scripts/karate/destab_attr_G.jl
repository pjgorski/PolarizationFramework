# Destabilizing polarization as a function of number of attributes. 
# Figure 1
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework
using Graphs

gr = KarateGraph()

specified_division = [
    [3, 14, 20, 2, 1, 4, 8, 18, 22, 12, 13, 5, 11, 6, 7, 17],
    [21, 30, 19, 24, 16, 27, 15, 23, 33, 34, 28, 29, 32, 26, 25, 31, 10, 9],
]

ns = [nv(gr)]
# ns = [5, 9]#[9,15,25]
reps = [900]
reps_dict = Dict(zip(ns, reps))
gs = [1:2:21...]
gs = [25, 29, 33, 37, 41, 45, 49, 55, 61, 67, 73, 81, 89, 97]

threshold = 0.5;
vs = [4, @onlyif("attr_types" == "OA", 1000)] #includes CA

attr_types = ["BA", "UA", "OA", "UPA"]


# sg = [-1.5, -1.1, -1.01, -0.99, -0.5, -0.1, -0.01]
# gammas=sort([sg..., 0, (-sg)...]);
# gammas = [2, 2.5, 3, 3.5]
gammas = [0.5, 1.5, 2, 2.5, 3, 3.5, 4, 6]
# gammas = [0.5]
# gammas = [6]

all_params = @strdict(ns, gs, threshold, vs, attr_types)
dicts = dict_list(all_params)
[d["reps"] = reps_dict[d["ns"]] for d in dicts] #adding reps
# [d["gammas"] = gammas for d in dicts] #adding gammas

for params in dicts
    n, g, attr_type, v, rep = let
        @unpack ns, threshold, reps, gs, attr_types, vs = params
        ns, gs, attr_types, vs, reps
    end
    println("Started ", @ntuple(g, attr_type, v))

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
        -1,
        rep,
        3000.0,
        "Heider9!";
        disp_each = 0,
        disp_more_every = 600,
        save_each = 600,
        files_folder = ["data", "karate-sims"],
        filename_prefix = "DesKarG",
        graph = gr,
        specified_division = specified_division,
    )

end
