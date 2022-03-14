# Preventing polarization from forming as a function of number of attributes. 
# Figure 3

using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework

ns = [5, 9]#[9,15,25]
reps = [100, 100]
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

all_params = @strdict(ns, gs, threshold, vs, attr_types)
dicts = dict_list(all_params)
[d["reps"] = reps_dict[d["ns"]] for d in dicts] #adding reps
# [d["gammas"] = gammas for d in dicts] #adding gammas

for params in dicts

    n, g, attr_type, v, rep = let
        @unpack ns, threshold, reps, gs, attr_types, vs = params
        ns, gs, attr_types, vs, reps
    end

    gammas = sg*g
    println("Started n=$n and g=$g and attr_type=", attr_type, " and v=$v.")
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
    r = using_curheider_attr(n, attr, gammas, rep, 3000., "Heider7!";
        disp_each = 0, disp_more_every = 600, save_each = 600, files_folder = ["data", "sims"], 
        filename_prefix = "NumerFig3")
end
