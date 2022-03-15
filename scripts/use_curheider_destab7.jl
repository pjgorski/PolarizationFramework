# Destabilizing polarization as a function of number of attributes. 
# Figure 4
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework

ns = [3, 5, 7, 9, 11, 13, 15]#[9,15,25]
# ns = [3, 7]
reps = Int.([10000, 1000, 100, 100, 100] / 10)
# reps = [10000, 1000]
# reps = Int.(reps/10*9)
reps_dict = Dict(zip(ns, reps))
gs = [5, 11]
threshold = 0.5;
vs = [4, @onlyif("attr_types" == "OA", 1000)] #includes CA

attr_types = ["BA", "UA", "OA", "UPA"]

gammas = [0.5, 1.5, 4, 6]

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
    println("Started ", @ntuple(n, g, attr_type, v, larger_size))

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
    r = using_curheider_attr_destab(
        n,
        attr,
        gammas,
        larger_size,
        rep,
        3000.0,
        "Heider7!",
        disp_each = 0,
        disp_more_every = 600,
        save_each = 600,
        files_folder = ["data", "sims"],
        filename_prefix = "DestabFig4",
    )

end
