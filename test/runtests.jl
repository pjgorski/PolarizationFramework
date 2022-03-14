# file for unit tests

using DrWatson
quickactivate(@__DIR__)

using PolarizationFramework
using Test

using DifferentialEquations

@testset "BinaryAttributes" begin
    g = 5;
    b = BinaryAttributes(g)
    nodes = 3;
    attr = get_attributes(b, nodes)

    @test size(attr) == (nodes, g) #test rozmiaru
    @test all(abs.(attr).==1) #test wartosci powinny byc +-1

    attr = ones(nodes, g)
    attr[1,:] .= -1
    attr[2,1:2] .= -1
    #attr is: [-1, -1, -1, -1, -1;
    #          -1, -1, 1, 1, 1;
    #           1, 1, 1, 1, 1]
    #This should give following weights:
    # [0.0  -0.2  -1.0
    #  0.0   0.0   0.2
    #  0.0   0.0   0.0]
    exp_weights = zeros(3,3)
    exp_weights[1,2] = -0.2;
    exp_weights[1,3] = -1.0;
    exp_weights[2,3] = 0.2;

    weights = get_attribute_layer_weights(b, attr)

    @test size(weights) == (nodes, nodes)
    @test weights == exp_weights
end

@testset "OrderedAttributes" begin
    g = 5;
    threshold = 0.55
    v = 3
    b = OrderedAttributes(g, threshold, v)
    nodes = 3;
    attr = get_attributes(b, nodes)

    @test min(attr...) == 1 #test rozmiaru
    @test max(attr...) == v #test wartosci powinny byc +-1

    attr = ones(Int, nodes, g)
    attr[1,:] .= 3
    attr[2,1:2] .= 3
    attr[2,3:5] .= 2
    #attr is:
    # [3  3  3  3  3
    #  3  3  2  2  2
    #  1  1  1  1  1]
    #This should give following weights:
    # [0.0  0.5  -0.9
    #  0.0  0.0  -0.3
    #  0.0  0.0   0.0]
    exp_weights = zeros(3,3)
    exp_weights[1,2] = 0.5;
    exp_weights[1,3] = -0.9;
    exp_weights[2,3] = -0.3;

    weights = get_attribute_layer_weights(b, attr)

    @test size(weights) == (nodes, nodes)
    @test weights ≈ exp_weights
end

@testset "UnorderedAttributes" begin
    g = 5;
    threshold = 0.55
    v = 3
    b = UnorderedAttributes(g, threshold, v)
    nodes = 3;
    attr = get_attributes(b, nodes)

    @test min(attr...) == 1 #test rozmiaru
    @test max(attr...) == v #test wartosci powinny byc +-1

    attr = ones(Int, nodes, g)
    attr[1,1:4] .= 3
    attr[2,1:3] .= 3
    attr[2,4:5] .= 2
    #attr is:
    # [3  3  3  3  1
    #  3  3  3  2  2
    #  1  1  1  1  1]
    #This should give following weights:
    # [0.0  0.1  -0.7
    #  0.0  0.0  -1.0
    #  0.0  0.0   0.0]
    exp_weights = zeros(3,3)
    exp_weights[1,2] = 0.1;
    exp_weights[1,3] = -0.7;
    exp_weights[2,3] = -1.1;

    weights = get_attribute_layer_weights(b, attr)

    @test size(weights) == (nodes, nodes)
    @test weights ≈ exp_weights
end

@testset "Simple simulation run" begin
    g = 5
    attr = BinaryAttributes(g)
    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider_attr(45, attr, 0., 1000., 
        Heider7!, AutoTsit5(Rodas5(autodiff = false)), true)
end

@testset "Simulation with initial conditions" begin
    n = 15
    g = 5
    gamma = 2.7
    attr = BinaryAttributes(g)

    rl_bal = init_random_balanced_relations(n)

    val0_attr = get_attributes(attr, n);
    al_w = get_attribute_layer_weights(attr, val0_attr);

    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider_attr(n, attr, 0., 1000., 
        Heider7!, AutoTsit5(Rodas5(autodiff = false)), false);
    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider_attr(n, attr, gamma, 1000., 
        Heider7!, AutoTsit5(Rodas5(autodiff = false)), false);
    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider_attr(n, attr, gamma, 1000., 
        Heider7!, AutoTsit5(Rodas5(autodiff = false)), false, rl_bal);
    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider_attr(n, attr, gamma, 1000., 
        Heider7!, AutoTsit5(Rodas5(autodiff = false)), false, (), al_w);
    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider_attr(n, attr, gamma, 1000., 
        Heider7!, AutoTsit5(Rodas5(autodiff = false)), false, rl_bal, al_w);
end

@testset "Series of simulations" begin
    n = 5
    g = 5
    attr = BinaryAttributes(g)
    gammas = [0.5, 1.5, 3.5]
    zmax = 10
    maxtime = 1000.
    ode_fun_name = "Heider7!"
    disp_each = 0.5
    disp_more_every = 60
    save_each = 60
    files_folder = ["test"]
    filename_prefix = "SoSNumerTest"

    using_curheider_attr(n, attr, gammas, zmax, maxtime, ode_fun_name; 
        disp_each = disp_each, disp_more_every = disp_more_every, 
        save_each = save_each, files_folder = files_folder, filename_prefix = filename_prefix)

    filename_prefix = "SoSDestabTest"
    larger_size = 3

    using_curheider_attr_destab(n, attr, gammas, larger_size, zmax, maxtime, ode_fun_name; 
        disp_each = disp_each, disp_more_every = disp_more_every, 
        save_each = save_each, files_folder = files_folder, filename_prefix = filename_prefix)
end

@testset "Longer series of simulations" begin
    n = 9
    g = 5
    attr = BinaryAttributes(g)
    gammas = [0.5, 1.5, 2., 3.5]
    zmax = 1000
    maxtime = 1000.
    ode_fun_name = "Heider7!"
    disp_each = 0.5
    disp_more_every = 60
    save_each = 60
    files_folder = ["test"]
    filename_prefix = "LSoSNumerTest"

    using_curheider_attr(n, attr, gammas, zmax, maxtime, ode_fun_name; 
        disp_each = disp_each, disp_more_every = disp_more_every, 
        save_each = save_each, files_folder = files_folder, filename_prefix = filename_prefix)

    filename_prefix = "LSoSDestabTest"
    larger_size = Int((n+1)/2 + 2)
    zmax = 300

    using_curheider_attr_destab(n, attr, gammas, larger_size, zmax, maxtime, ode_fun_name; 
        disp_each = disp_each, disp_more_every = disp_more_every, 
        save_each = save_each, files_folder = files_folder, filename_prefix = filename_prefix)
end