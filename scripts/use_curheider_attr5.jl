#figure 3
using MeasureHeiderBalance
using Attributes

ns = [3]#[9,15,25]
reps = [10000]
reps_dict = Dict(zip(ns, reps))
gs = [5:2:11...]
threshold = 0.5;
vs = [32];
#vs = [100, 1000]

#its = its[4:end]

sg = [0.01:0.01:0.21...]

bal = [""]#, "init_random_balanced_enhanced"]

its = [Iterators.product(ns, gs, bal, vs)...]

for (n, g, creator, v) in its
    rep = reps_dict[n]
    gammas=sort(sg)*g;
    println("Started n=$n and g=$g and creator=", creator, " and v=$v.")
    attr = UnorderedAttributes(g, threshold, v)
    r = using_curheider_attr(n, g, gammas, rep, 3000., "Heider7!", creator,
        attr, length(gammas)*rep/10, 600, 600, "attr_curheider_results_undef_balanced")
end
