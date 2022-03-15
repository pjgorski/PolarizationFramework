# Preventing polarization from forming as a function of number of attributes. 
# Figure 4
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework

ns = [3, 5, 7, 9, 11, 13, 15]#[9,15,25]
# ns = [11,13,15]
#for test simulations number of repetition is decreased. 
reps = Int.([10000, 1000, 1000, 1000, 500, 500, 500] / 10)
# reps = [500, 500, 500]
reps_dict = Dict(zip(ns, reps))
gs = [5, 11]

threshold = 0.5;
vs = [4, @onlyif("attr_types" == "OA", 1000)] #includes CA

attr_types = ["BA", "UA", "OA", "UPA"]

gammas = [0.1, 0.5, 1.5, 2.5, 4]

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
        "Heider7!";
        disp_each = 0,
        disp_more_every = 600,
        save_each = 600,
        files_folder = ["data", "sims"],
        filename_prefix = "NumerFig4",
    )
end
