# Figure 3!!!
#destabilization
using MeasureHeiderBalance
using Attributes

ns = [5, 9]#[9,15,25]
# ns = [9]
reps = [900, 900]
reps_dict = Dict(zip(ns, reps))
# gs = [1, 3, 5, 7,11]
gs = [5, 11]
# gs = [11]
threshold = 0.5;
vs = [4]
vs = [1000]

attr_types = ["BA", "UA", "OA", "UPA"]
attr_types = ["UA", "OA", "UPA"]
attr_types = ["OA"]

# sg = [-1.5, -1.1, -1.01, -0.99, -0.5, -0.1, -0.01]
# sg = [0, 0.01, 0.1, 0.5, 1.5]
# sg = [0.1, 0.3, 1.2]
sg = [0:0.1:1.5...]

its = [Iterators.product(ns, gs, vs, attr_types)...]

for (n, g, v, attr_type) in its
    rep = reps_dict[n]
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

    larger_sizes = Int.(ceil(n/2):n)

    for larger_size in larger_sizes

        println("Started n=$n and g=$g and attr_type=", attr_type,
            " and v=$v and larger_size=$larger_size.")

        r = using_curheider_attr_destab(n, g, gammas, larger_size,
            rep, 3000., "Heider7!",
            attr, length(gammas)*rep/10, 600, 600,
            "")
    end
end
