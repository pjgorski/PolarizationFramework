# Destabilizing polarization as a function of number of attributes. 
# Figure 3
using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework

ns = [5, 9]#[9,15,25]
# ns = [9]
reps = [10, 10]
reps_dict = Dict(zip(ns, reps))
# gs = [1, 3, 5, 7,11]
gs = [5, 11]
threshold = 0.5;
vs = [4, @onlyif("attr_types" == "OA",  1000)] #includes CA

attr_types = ["BA", "UA", "OA", "UPA"]

# sg = [-1.5, -1.1, -1.01, -0.99, -0.5, -0.1, -0.01]
# sg = [0, 0.01, 0.1, 0.5, 1.5]
# sg = [0.1, 0.3, 1.2]
sg = [0:0.1:1.5...]

larger_sizes = [@onlyif(i <="ns" < 2*i,  i) for i in 1:maximum(ns)]

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

    #gammas=sort([sg..., 0, (-sg)...])*g;
    gammas = sg*g
    # println("Started n=$n and g=$g and attr_type=", attr_type, " and v=$v.")
    if attr_type == "UA"
        attr = UnorderedAttributes(g, threshold, v)
        gammas = [sg..., gammas...]
    elseif attr_type == "BA"
        attr = BinaryAttributes(g)
        gammas = [sg..., gammas...]
    elseif attr_type == "OA"
        attr = OrderedAttributes(g, threshold, v)
    elseif attr_type == "UPA"
        attr = UnorderedPositiveAttributes(g, threshold, v)
    else
        throw(attr_type)
    end

    r = using_curheider_attr_destab(n, attr, gammas, larger_size,
        rep, 3000., "Heider7!",
        disp_each = 0, disp_more_every = 600, save_each = 600, files_folder = ["data", "sims"], 
        filename_prefix = "DestabFig3")
end
