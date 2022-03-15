using Base.Test
using Attributes

g = 5;
threshold = 0.55
v = 3
b = UnorderedAttributes(g, threshold, v)
nodes = 3;
attr = get_attributes(b, nodes)

@test min(attr...) == 1 #test rozmiaru
@test max(attr...) == v #test wartosci powinny byc +-1

attr = ones(Int, nodes, g)
attr[1, 1:4] = 3
attr[2, 1:3] = 3
attr[2, 4:5] = 2
#attr is:
# [3  3  3  3  1
#  3  3  3  2  2
#  1  1  1  1  1]
#This should give following weights:
# [0.0  0.1  -0.7
#  0.0  0.0  -1.0
#  0.0  0.0   0.0]
exp_weights = zeros(3, 3)
exp_weights[1, 2] = 0.1;
exp_weights[1, 3] = -0.7;
exp_weights[2, 3] = -1.1;

weights = get_attribute_layer_weights(b, attr)

@test size(weights) == (nodes, nodes)
@test weights ≈ exp_weights
