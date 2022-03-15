using Base.Test
using Attributes

g = 5;
b = BinaryAttributes(g)
nodes = 3;
attr = get_attributes(b, nodes)

@test size(attr) == (nodes, g) #test rozmiaru
@test all(abs.(attr) .== 1) #test wartosci powinny byc +-1

attr = ones(nodes, g)
attr[1, :] = -1
attr[2, 1:2] = -1
#attr is: [-1, -1, -1, -1, -1;
#          -1, -1, 1, 1, 1;
#           1, 1, 1, 1, 1]
#This should give following weights:
# [0.0  -0.2  -1.0
#  0.0   0.0   0.2
#  0.0   0.0   0.0]
exp_weights = zeros(3, 3)
exp_weights[1, 2] = -0.2;
exp_weights[1, 3] = -1.0;
exp_weights[2, 3] = 0.2;

weights = get_attribute_layer_weights(b, attr)

@test size(weights) == (nodes, nodes)
@test weights == exp_weights
# @test false
