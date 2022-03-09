# Figure 2
#Unordered attributes
using MeasureHeiderBalance
using Attributes

ns = [5, 9]#[9,15,25]
reps = [1000, 1000]
reps_dict = Dict(zip(ns, reps))
gs = [1, 3, 5, 7,11]
gs = [5]
threshold = 0.5;
vs = [2, 4,8, 16, 32, 64, 128, 256, 512, 1000];
# vs=2;
#vs = [100, 1000]

attr_types = ["UA", "OA", "UPA"]
# attr_types = ["UPA"]

#its = its[4:end]

sg = [-1.5, -1.1, -1.01, -0.99, -0.5, -0.1, -0.01]
sg = [0, 0.01, 0.1, 0.5, 1.5]
sg = [0.1, 0.3, 1.2]

bal = [""]#, "init_random_balanced_enhanced"]

its = [Iterators.product(ns, gs, bal, vs, attr_types)...]

for (n, g, creator, v, attr_type) in its
    rep = reps_dict[n]
    #gammas=sort([sg..., 0, (-sg)...])*g;
    gammas = sg*g
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
