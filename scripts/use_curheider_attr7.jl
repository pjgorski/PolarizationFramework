# Figure 1
#Unordered attributes
using MeasureHeiderBalance
using Attributes

ns = [3,5,7,9,11,13,15]#[9,15,25]
ns = [11,13,15]
reps = [10000,1000,1000,1000,100,100,100]
reps = [500, 500, 500]
reps_dict = Dict(zip(ns, reps))
gs = [5, 11]
# gs = [25,29,33,37,41,45,49,55,61,67,73,81,89,97]
# gs = 3;
threshold = 0.5;
vs = [4];
vs = [1000]

attr_types = ["BA", "UA", "OA", "UPA"]
# attr_types = ["BA", "OA", "UPA"]
attr_types = ["OA"]

#its = its[4:end]

sg = [-1.5, -1.1, -1.01, -0.99, -0.5, -0.1, -0.01]
gammas=sort([sg..., 0, (-sg)...]);
gammas = [0.1, 0.5, 1.5, 2.5, 4]
# gammas = [1.5]

bal = [""]#, "init_random_balanced_enhanced"]

its = [Iterators.product(ns, gs, bal, vs, attr_types)...]

for (n, g, creator, v, attr_type) in its
    rep = reps_dict[n]
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
    r = using_curheider_attr(n, g, gammas, rep, 3000., "Heider7!", creator,
        attr, length(gammas)*rep/10, 600, 600,
        "attr_curheider_results_undef_balanced_paradise_deltas")
end
