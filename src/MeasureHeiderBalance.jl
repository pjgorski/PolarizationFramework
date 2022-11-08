# This file contains main function to generate network ans solve differential equations 
# and to run the series of repetitions of given model. 
# Apart from that, functions also make some initial and final states analysis and save the results to Matlab file. 
# Code is optimized

# module MeasureHeiderBalance

using DrWatson
using DifferentialEquations
using Plots
using LinearAlgebra
using Dates
using Graphs

# function to measure HB in the case of a single-layered system with attributes
# n - number of agents
# attr - attribute type; contains all information about attributes
# gamma - coupling strength
# maxtime - maximal time of calculating differential equations
# ode_fun - function calculating derivatives
# solver - function solving different equations
# show_plot - parameter whether the solution will be plotted. If not transition values are not saved. 
# input is optional extra parameters with initial conditions of RL and AL (u0, xy_attr). 
# If input is empty, random initial conditions are generated. 
# If not, the first value (possibly empty) should contain initial conditions for RL. 
# The 2nd value of input may contain initial conditions for AL. 
# 
# Returns the tuple of: 
# ishb_sim_par - array of Boolean values [is it balanced state, is AL in balanced state,
# similarity between RL and AL, balanced ratio of RL, is it paradise state,
# is it hell state, counts od different triads (Deltas: Delta0, Delta1, Delta2, Delta3), 
# is it weakly balanced state, local polarization measure, is the system globally polarized],
# time of simulation to reach stable solution (or maxtime if not reached),
# values of RL weights at the end of simulation, 
# initial conditions of RL, 
# initial conditions of AL, 
# whole sol (if transition values are not saved, see show_plot, then they are not here).
function calc_heider_attr(
    n::Int,
    attr::AbstractAttributes,
    gamma::Float64,
    maxtime::Float64,
    ode_fun::Function,
    solver,
    show_plot::Bool,
    input...;
    all_links_mat = [], #if existing, this should be of size (n,n) and this would show which relations may change and which cannot because they do not exist
    all_triads = [], # list of all triads
    # all_links = [], # list of all links
    # triads_around_links_dict = [], #dict with which links form triads around other links
    triads_count_mat = [], # triad count matrix for each link. This is used in calculation of derivatives
    link_indices = [], # indices from adjacency matrix of links belonging to triads, i.e., links that are interesting in terms dynamics
    link_pairs = [], # vector (Vector{Vector{Tuple{Int64, Int64}}}) corresponding to link_indices. Each element is the vector with tuples of two link indices that close the triad with the given link
    link_pairs_triad_cnt = [], # vector corresponding to link_indices. Contains number of triads for each link. 
)

    if !isempty(all_links_mat)
        prob,
        cb,
        u0,
        xy_attr,
        mask,
        all_triads,
        triads_count_mat,
        link_indices,
        link_pairs,
        link_pairs_triad_cnt = initialize_calc_heider_attr_incomplete(
            n,
            attr,
            gamma,
            maxtime,
            ode_fun,
            all_links_mat,
            input...;
            all_triads = all_triads,
            triads_count_mat = triads_count_mat,
            link_indices = link_indices,
            link_pairs = link_pairs,
            link_pairs_triad_cnt = link_pairs_triad_cnt,
        )

        if ode_fun in [Heider9!, Heider92!]
            link_indices = mask
            mask = ones(size(link_indices))
            # mask .*= all_links_mat
        else
            mask .*= all_links_mat
        end
    else
        #initial conditions
        need_init_u0 = true
        need_init_xy = true

        if length(input) > 0
            if !isempty(input[1])
                need_init_u0 = false
                u0 = input[1]
            end
            if length(input) == 2
                if !isempty(input[2])
                    need_init_xy = false
                    xy_attr = input[2]
                end
            end
        end

        if need_init_u0
            u0 = triu((rand(n, n) * 2) .- 1, 1)
        end
        if need_init_xy
            val0_attr = get_attributes(attr, n)
            xy_attr = get_attribute_layer_weights(attr, val0_attr)
        end

        # help variable
        mask = triu(trues(size(u0)), 1)

        condition_here = (u, t, integrator) -> condition2(u, t, integrator, mask)
        affect!(integrator) = terminate!(integrator)
        cb = DiscreteCallback(condition_here, affect!)

        #ode and model parameters
        tspan = (0.0, maxtime)

        lay1mul = zeros(n, n)
        x_sim = zeros(n, n)

        p = (n, gamma .* xy_attr, lay1mul, x_sim, mask)
        if ode_fun in [Heider72!, Heider722!]
            p = (p..., triads_count_mat)
            if ode_fun == Heider722!
                p = (p..., zeros(Bool, n, n))
            end
        elseif ode_fun == Heider73!
            p[2] .*= triads_count_mat
        end

        #solving
        prob = ODEProblem(ode_fun, u0, tspan, p)
    end
    sol = solve(
        prob,
        solver,
        reltol = 1e-6,
        abstol = 1e-12,
        callback = cb,
        isoutofdomain = (u, p, t) -> any(x -> abs.(x) >= 1, u),
        save_everystep = show_plot,
    )

    #estimating output
    end_signs = sol.u[end] .* mask
    # println(end_signs)

    paradise = is_paradise(end_signs, n)
    hell = is_hell(end_signs, n)
    if isempty(all_links_mat)
        Deltas = get_triad_counts(end_signs, n; hlp = x_sim)
        # hb_in_rels = is_hb(sol.u[end], n)
        hb_in_attrs = is_hb(xy_attr, n)
        br = get_balanced_ratio_efficient(
            end_signs,
            Int(n * (n - 1) * (n - 2) / 6);
            hlp = x_sim,
        )
        br2 = get_balanced_ratio_efficient(
            sign.(end_signs),
            Int(n * (n - 1) * (n - 2) / 6);
            hlp = x_sim,
        )
        sim = get_similarity(end_signs, xy_attr, n)
    else
        x_sim = zeros(size(end_signs))
        if end_signs isa Vector
            Deltas =
                get_triad_counts(end_signs, link_pairs, link_pairs_triad_cnt; hlp = x_sim)
            br = get_balanced_ratio_not_complete(
                end_signs,
                link_pairs,
                link_pairs_triad_cnt;
                hlp = x_sim,
            )
            br2 = get_balanced_ratio_not_complete(
                sign.(end_signs),
                link_pairs,
                link_pairs_triad_cnt;
                hlp = x_sim,
            )

            sim = get_similarity2(end_signs, xy_attr[link_indices], length(end_signs))
            hb_in_attrs =
                get_balanced_ratio_not_complete(
                    sign.(xy_attr[link_indices]),
                    link_pairs,
                    link_pairs_triad_cnt;
                    hlp = x_sim,
                ) == 1
        else
            Deltas = get_triad_counts(end_signs, all_triads; hlp = x_sim)
            br = get_balanced_ratio_not_complete(end_signs, length(all_triads); hlp = x_sim)
            br2 = get_balanced_ratio_not_complete(
                sign.(end_signs),
                length(all_triads);
                hlp = x_sim,
            )
            sim = get_similarity2(end_signs, xy_attr, sum(sign.(end_signs) .> 0))

            hb_in_attrs =
                get_balanced_ratio_not_complete(
                    sign.(xy_attr .* mask),
                    length(all_triads);
                    hlp = x_sim,
                ) == 1
        end
    end
    weak_balance_in_complete_graph = Deltas[1+1] == 0
    local_polarization = get_local_polarization(Deltas)
    global_polarization = weak_balance_in_complete_graph && (!paradise)

    hb_in_rels = br2 == 1


    ishb_sim_par = [
        hb_in_rels,
        hb_in_attrs,
        sim,
        br,
        paradise,
        hell,
        Deltas,
        weak_balance_in_complete_graph,
        local_polarization,
        global_polarization,
    ]

    if show_plot

        if end_signs isa Vector
            h = plot(
                sol,
                linewidth = 5,
                title = "Solution to the linear ODE with a thick line",
                xaxis = "Time (t)",
                yaxis = "weights(t)",
                ylim = [-1, +1],
            )
        else
            h = plot(
                sol,
                linewidth = 5,
                title = "Solution to the linear ODE with a thick line",
                xaxis = "Time (t)",
                yaxis = "weights(t)",
                ylim = [-1, +1],
                vars = reshape(1:n^2, n, n)[mask],
            )
        end
        legend = false
        display(h)
        #gui()
    end

    return (ishb_sim_par, sol.t[end], end_signs, u0, xy_attr, sol)
end
export calc_heider_attr

function initialize_calc_heider_attr_incomplete(
    n::Int,
    attr::AbstractAttributes,
    gamma::Float64,
    maxtime::Float64,
    ode_fun::Function,
    all_links_mat::Matrix, #this should be of size (n,n) and this would show which relations may change and which cannot because they do not exist
    input...;
    all_triads = [], # list of all triads
    # all_links = [], # list of all links
    # triads_around_links_dict = [], #dict with which links form triads around other links
    triads_count_mat = [], # triad count matrix for each link. This is used in calculation of derivatives
    link_indices = [], # indices from adjacency matrix of links belonging to triads, i.e., links that are interesting in terms dynamics
    link_pairs = [], # vector (Vector{Vector{Tuple{Int64, Int64}}}) corresponding to link_indices. Each element is the vector with tuples of two link indices that close the triad with the given link
    link_pairs_triad_cnt = [], # vector corresponding to link_indices. Contains number of triads for each link. 
)

    # help variable
    mask = triu(trues(n, n), 1)
    mask .*= all_links_mat

    if isempty(all_triads)
        all_triads = get_triads(all_links_mat)
    end


    if ode_fun in [Heider72!, Heider722!, Heider73!]
        if isempty(triads_count_mat)
            all_links = get_links_in_triads(all_triads)

            triads_around_links_dict = get_triangles_around_links(all_triads)
            counts = link_triangles_count(triads_around_links_dict; links = all_links)

            if ode_fun in [Heider72!, Heider722!]
                triads_count_mat = link_triangles_mat_inv(n, all_links, counts)
            elseif ode_fun == Heider73!
                triads_count_mat = link_triangles_mat(n, all_links, counts)
            end
        end
    end
    if ode_fun in [Heider9!, Heider92!]
        if isempty(link_indices)
            link_indices = findall(triu(all_links_mat, 1)[:] .> 0)
        end
        nl = length(link_indices)

        if isempty(link_pairs)
            dict_e = get_triangles_around_links(all_triads)
            all_links = get_links_in_triads(all_triads)
            link_pairs = get_triangles_around_links(dict_e, all_links)
        end

        if isempty(link_pairs_triad_cnt)
            link_pairs_triad_cnt = [length(link) for link in link_pairs]
        end
    end


    need_init_u0 = true
    need_init_xy = true
    if length(input) > 0
        if !isempty(input[1])
            need_init_u0 = false
            u0 = input[1]
        end
        if length(input) == 2
            if !isempty(input[2])
                need_init_xy = false
                xy_attr = input[2]
            end
        end
    end

    if need_init_u0
        u0 = triu((rand(n, n) * 2) .- 1, 1)
    end
    if need_init_xy
        val0_attr = get_attributes(attr, n)
        xy_attr = get_attribute_layer_weights(attr, val0_attr)
    end

    u0 .*= mask
    xy_attr .*= mask

    if ode_fun in [Heider9!, Heider92!]
        u0_inc = u0[link_indices]
        xy_attr_inc = xy_attr[link_indices]
    end

    if ode_fun in [Heider9!, Heider92!]
        condition_here = (u, t, integrator) -> condition2(u, t, integrator, 1:nl)
    else
        condition_here = (u, t, integrator) -> condition2(u, t, integrator, mask)
    end

    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition_here, affect!)

    #ode and model parameters
    tspan = (0.0, maxtime)

    if ode_fun in [Heider9!, Heider92!]
        lay1mul = zeros(nl)

        p = (gamma .* xy_attr_inc, lay1mul, link_pairs, link_pairs_triad_cnt)

        if ode_fun in [Heider92!]
            p = (p..., zeros(Bool, nl))
        end
    else

        lay1mul = zeros(n, n)
        x_sim = zeros(n, n)

        p = (n, gamma .* xy_attr, lay1mul, x_sim, mask)
        if ode_fun in [Heider72!, Heider722!]
            p = (p..., triads_count_mat)
            if ode_fun == Heider722!
                p = (p..., zeros(Bool, n, n))
            end
        elseif ode_fun == Heider73!
            p[2] .*= triads_count_mat
        end
    end

    #solving
    if ode_fun in [Heider9!, Heider92!]
        prob = ODEProblem(ode_fun, u0_inc, tspan, p)
        return prob,
        cb,
        u0,
        xy_attr,
        link_indices,
        all_triads,
        triads_count_mat,
        link_indices,
        link_pairs,
        link_pairs_triad_cnt
    else
        prob = ODEProblem(ode_fun, u0, tspan, p)
        return prob,
        cb,
        u0,
        xy_attr,
        mask,
        all_triads,
        triads_count_mat,
        link_indices,
        link_pairs,
        link_pairs_triad_cnt
    end
end
export initialize_calc_heider_attr_incomplete

# Function creating results file. It parses parameters creating a unique file name and saves initially this file. 
# For parameters description see `using_heider_attr`. 
# Returns a Tuple with:
#   object of `Result` DataType,
#   `filename`
#   filename without extension
function initialize_file(
    n::Int,
    attr::AbstractAttributes,
    gammas::Vector{Float64},
    maxtime::Float64,
    ode_fun_name::String,
    files_folder::Vector{String},
    filename_prefix::String,
)

    r = Result(n, attr, gammas, maxtime, ode_fun_name)

    prefix = filename_prefix * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    file_params = savename(prefix, r, "mat", sort = false)
    file_params = replace(file_params, "attr_degeneracy" => "v")

    filename = projectdir(files_folder..., file_params)
    save_result(r, filename) #''allocating'' place
    return r, filename
end

# Function using `calc_heider_attr` function to simulate a number of repetitions
# of solving the system having given parameters and random initial conditions:
# Thus, this simulates the influence of attributes on preventing from forming 
# a polarized state. 
# 
# Parameters: 
# n - number of agents
# attr - attribute type; contains all information about attributes
# gammas - a vector of coupling strengths to simulate
# zmax - number of repetitions
# maxtime - maximal time of calculating differential equations
# ode_fun_name - function name calculating derivatives (`string`)
# 
# Additional named arguments are:
# disp_each - display status every how many ratio of realizations (if 0 then do not display) (default 0.5)
# disp_more_every - display status every certain number of seconds (if 0, then do not display) (default 600)
# save_each - how often (in seconds) should results be saved (if 0, then save at the end) (default 600)
# files_folder - folder name inside the project the results should be saved (default "data"). 
#       This can be also an array of folder names in the correct folder structure. 
# filename_prefix - string that should start each simulation results' file (default "")
# all_links_mat - adjacency matrix of links that should be considered. This should be used if not the complete network should be used. 
# kwargs - optional arguments used in the case of not complete graph topology (not empty `all_links_mat`). 
#       kwargs contain variables related to the topology. Goal is to speed up calculations. 
#       See `calc_heider_attr` description for details. 
# 
# `disp_more_every` and `save_each` work like that, that if the specified time has past
# then sth is displayed/saved. But it doesn't mean exact time of action.
# When the specified time has past means when the current function
# (heider function with parameters) finishes.
#  
# Returns `Result` object. 
function using_heider_attr(
    n::Int,
    attr::AbstractAttributes,
    gammas::Vector{Float64},
    zmax::Int,
    maxtime::Float64,
    ode_fun_name::String;
    disp_each = 0.5,
    disp_more_every = 600,
    save_each = 600,
    files_folder::Vector{String} = ["data"],
    filename_prefix::String = "",
    graph::Graph = Graph(), # `graph` or `all_links_mat` should be given if the considered network is not complete
    all_links_mat = [], # `graph` or `all_links_mat` should be given if the considered network is not complete
    kwargs...,
)

    ode_fun = getfield(PolarizationFramework, Symbol(ode_fun_name))
    solver = AutoTsit5(Rodas5(autodiff = false))

    if nv(graph) > 0 || ~isempty(all_links_mat)
        if nv(graph) > 0
            assumed_n = nv(graph)
            @assert nv(graph) == n "Wrong specified number of nodes `n` and given number of nodes in `graph."
        else
            assumed_n = size(all_links_mat)[1]
        end
        @assert assumed_n == n "Wrong specified number of nodes `n` and given number of nodes in `graph` or `all_links_mat`."

        kwargs_dict = Dict(kwargs)

        if isempty(all_links_mat)
            all_links_mat = adjacency_matrix(graph)
            all_triads = get_triads(all_links_mat)
            all_links = get_links_in_triads(all_triads)
            all_links_mat = get_adj_necessary_links(n, all_links; typ = Float64)

            kwargs = (kwargs..., all_triads = all_triads)
            kwargs_dict = Dict(pairs(kwargs))
            # kwargs_dict[:all_triads] = all_triads
        elseif !haskey(kwargs_dict, :all_triads)
            all_triads = get_triads(all_links_mat)

            kwargs = (kwargs..., all_triads = all_triads)
            # kwargs_dict[:all_triads] = all_triads
            kwargs_dict = Dict(pairs(kwargs))
        end

        if ode_fun_name in ["Heider72!", "Heider722!"]
            if !haskey(kwargs_dict, :triads_count_mat)
                triads_around_links_dict =
                    get_triangles_around_links(kwargs_dict[:all_triads])
                all_links = get_links_in_triads(kwargs_dict[:all_triads])
                counts = link_triangles_count(triads_around_links_dict; links = all_links)
                triads_count_mat = link_triangles_mat_inv(n, all_links, counts)

                kwargs = (kwargs..., triads_count_mat = triads_count_mat)
                kwargs_dict[:triads_count_mat] = triads_count_mat
            end
        elseif ode_fun_name in ["Heider9!", "Heider92!"]
            if !haskey(kwargs_dict, :link_indices)
                link_indices = findall(triu(all_links_mat, 1)[:] .> 0)

                kwargs = (kwargs..., link_indices = link_indices)
                kwargs_dict[:link_indices] = link_indices
            end

            if !haskey(kwargs_dict, :link_pairs)
                triads_around_links_dict =
                    get_triangles_around_links(kwargs_dict[:all_triads])
                all_links = get_links_in_triads(kwargs_dict[:all_triads])

                link_pairs = get_triangles_around_links(triads_around_links_dict, all_links)

                kwargs = (kwargs..., link_pairs = link_pairs)
                kwargs_dict[:link_pairs] = link_pairs
            end

            if !haskey(kwargs_dict, :link_pairs_triad_cnt)
                link_pairs_triad_cnt = [length(link) for link in kwargs_dict[:link_pairs]]

                kwargs = (kwargs..., link_pairs_triad_cnt = link_pairs_triad_cnt)
                kwargs_dict[:link_pairs_triad_cnt] = link_pairs_triad_cnt
            end
        end

    end

    r, filename = initialize_file(
        n,
        attr,
        gammas,
        maxtime,
        ode_fun_name,
        files_folder,
        filename_prefix,
    )

    # time preparation
    if disp_more_every != 0
        time_disp = time()
    end

    if save_each != 0
        time_save = time()
    end

    firstline = 1
    realization_counter = 0
    for i = 1:length(gammas)
        gamma1 = gammas[i]

        #zeroing data arrays
        HB = zeros(zmax)
        HB_x = zeros(zmax)
        BR = zeros(zmax) #ratio of balanced triads
        HB_attr = zeros(zmax)
        HB_only_weights = zeros(zmax)
        paradise = zeros(zmax)
        hell = zeros(zmax)
        Deltas = zeros(4, zmax)
        weak_balance_in_complete_graph = zeros(zmax)
        local_polarization = zeros(zmax)
        global_polarization = zeros(zmax)
        initial_neg_links_count = zeros(zmax)
        links_destab_changed = zeros(4, zmax)

        sim = zeros(zmax)
        x_attr_sim = zeros(zmax)
        stab = zeros(zmax)
        times = zeros(zmax)
        for rep = 1:zmax
            #simulation
            (ishb_sim_par, t, u, u0, xy_attr, sol) = calc_heider_attr(
                n,
                attr,
                gamma1,
                maxtime,
                ode_fun,
                solver,
                false;
                all_links_mat = all_links_mat,
                kwargs...,
            )
            realization_counter += 1

            #work on results
            HB_x[rep],
            HB_attr[rep],
            x_attr_sim[rep],
            BR[rep],
            paradise[rep],
            hell[rep],
            Deltas[:, rep],
            weak_balance_in_complete_graph[rep],
            local_polarization[rep],
            global_polarization[rep] = ishb_sim_par

            if u isa Vector
                link_indices = (; kwargs...).link_indices
                initial_neg_links_count[rep] = sum(u0[link_indices] .< 0)
                links_destab_changed[3, rep] = sum(u[u0[link_indices].>0] .< 0) #number of initial pos links that changed to negative
                links_destab_changed[4, rep] = sum(u[u0[link_indices].<0] .> 0) #number of initial neg links that changed to positive
            else
                initial_neg_links_count[rep] = sum(u0 .< 0)
                links_destab_changed[3, rep] = sum(u[u0.>0] .< 0) #number of initial pos links that changed to negative
                links_destab_changed[4, rep] = sum(u[u0.<0] .> 0) #number of initial neg links that changed to positive
            end

            if t < maxtime #we have stability
                HB[rep] = ishb_sim_par[1]
                stab[rep] = 1
            end

            times[rep] = t

            #checking time/repetitions and eventually informing about progress and saving
            if disp_more_every != 0
                if disp_more_every < time() - time_disp
                    time_disp = time()
                    # displaying
                    g = attr.g
                    v = r.attr_degeneracy
                    a = r.attr_name
                    # display_res(@ntuple(rep, a, n, g, v, gamma1))
                    display_res(; rep, a, n, g, v, gamma1)
                end
            end

            if save_each != 0
                if save_each < time() - time_save
                    time_save = time()
                    # partial saving
                    fields = (
                        HB,
                        HB_x,
                        HB_attr,
                        sim,
                        x_attr_sim,
                        BR,
                        paradise,
                        hell,
                        initial_neg_links_count,
                        links_destab_changed,
                        Deltas,
                        weak_balance_in_complete_graph,
                        local_polarization,
                        global_polarization,
                        stab,
                        times,
                        i,
                        rep,
                        firstline,
                    )
                    update_result!(r, fields)
                    save_result(r, filename)
                    save_result(r, filename, ext = "jld2")
                end
            end

            if disp_each != 0
                if disp_each <= realization_counter / zmax
                    realization_counter = 0
                    # displaying
                    g = attr.g
                    v = r.attr_degeneracy
                    a = r.attr_name
                    # display_res(@ntuple(rep, a, n, g, v, gamma1))
                    display_res(; rep, a, n, g, v, gamma1)
                end
            end
        end
        fields = (
            HB,
            HB_x,
            HB_attr,
            sim,
            x_attr_sim,
            BR,
            paradise,
            hell,
            initial_neg_links_count,
            links_destab_changed,
            Deltas,
            weak_balance_in_complete_graph,
            local_polarization,
            global_polarization,
            stab,
            times,
            i,
            zmax,
            firstline,
        )
        update_result!(r, fields)
        firstline += 2
    end

    save_result(r, filename)
    save_result(r, filename, ext = "jld2")

    return r
end
export using_heider_attr

# Function using `calc_heider_attr` function to simulate a number of repetitions
# of solving the system having given parameters and almost a balanced RL.
# Agents are divided into two groups and form positive links (+0.99) inside those groups
# and negative (-0.99) to agents outside their group. 
# Thus, this simulates destabilization of the polarized, initial state. 
# One can also specify a network that is simulated. It can be done in two ways:
#   by assigning an adjacency matrix to variable `all_links_mat`.
#       For speed purposes, it is recommended that only links that are in triads 
#       are given in the adjacency matrix. 
#   by assigning a keyword argument `graph` that is of the type `Graph` from `Graphs` package. 
# One may also specify the division of nodes into two groups 
# by specifying the keyword argument `specified_division`.
# 
# Parameters: 
# n - number of agents
# attr - attribute type; contains all information about attributes
# gammas - a vector of coupling strengths to simulate
# larger_size - size of the larger group in the initially balanced network
# zmax - number of repetitions
# maxtime - maximal time of calculating differential equations
# ode_fun_name - function name calculating derivatives (`string`)
# 
# Additional named arguments are:
# disp_each - display status every how many ratio of realizations (if 0 then do not display) (default 0.5)
# disp_more_every - display status every certain number of seconds (if 0, then do not display) (default 600)
# save_each - how often (in seconds) should results be saved (if 0, then save at the end) (default 600)
# files_folder - folder name inside the project the results should be saved (default "data")
#       This can be also an array of folder names in the correct folder structure. 
# filename_prefix - string that should start each simulation results' file (default "")
# graph - (type `Graphs.Graph`) specifies connections if sparse networks are considered. 
# all_links_mat - specifies adjacency matrix if sparse networks are considered.
# specified_division - (Vector of Vectors) specifies the initial division of nodes into groups. 
#     If this parameter is given, `larger_size` does not matter. 
# 
# `disp_more_every` and `save_each` work like that, that if the specified time has past
# then sth is displayed/saved. But it doesn't mean exact time of action.
# When the specified time has past means when the current function
# (heider function with parameters) finishes.
#  
# Returns `Result` object. 
function using_heider_attr_destab(
    n::Int,
    attr::AbstractAttributes,
    gammas::Vector{Float64},
    larger_size::Int,
    zmax::Int,
    maxtime::Float64,
    ode_fun_name::String;
    disp_each = 0.5,
    disp_more_every = 600,
    save_each = 600,
    files_folder::Vector{String} = ["data"],
    filename_prefix::String = "",
    graph::Graph = Graph(), # `graph` or `all_links_mat` should be given if the considered network is not complete
    all_links_mat = [], # `graph` or `all_links_mat` should be given if the considered network is not complete
    specified_division = [], # if incomplete network is considered, here a specific initial division of nodes into two groups can be given. 
    kwargs...,
)

    ode_fun = getfield(PolarizationFramework, Symbol(ode_fun_name))
    solver = AutoTsit5(Rodas5(autodiff = false))

    if nv(graph) > 0 || ~isempty(all_links_mat)
        if nv(graph) > 0
            assumed_n = nv(graph)
            @assert nv(graph) == n "Wrong specified number of nodes `n` and given number of nodes in `graph."
        else
            assumed_n = size(all_links_mat)[1]
        end
        @assert assumed_n == n "Wrong specified number of nodes `n` and given number of nodes in `graph` or `all_links_mat`."

        kwargs_dict = Dict(kwargs)

        if isempty(all_links_mat)
            all_links_mat = adjacency_matrix(graph)
            all_triads = get_triads(all_links_mat)
            all_links = get_links_in_triads(all_triads)
            all_links_mat = get_adj_necessary_links(n, all_links; typ = Float64)

            kwargs = (kwargs..., all_triads = all_triads)
            kwargs_dict = Dict(pairs(kwargs))
            # kwargs_dict[:all_triads] = all_triads
        elseif !haskey(kwargs_dict, :all_triads)
            all_triads = get_triads(all_links_mat)

            kwargs = (kwargs..., all_triads = all_triads)
            # kwargs_dict[:all_triads] = all_triads
            kwargs_dict = Dict(pairs(kwargs))
        end

        if ode_fun_name in ["Heider72!", "Heider722!"]
            if !haskey(kwargs_dict, :triads_count_mat)
                triads_around_links_dict =
                    get_triangles_around_links(kwargs_dict[:all_triads])
                all_links = get_links_in_triads(kwargs_dict[:all_triads])
                counts = link_triangles_count(triads_around_links_dict; links = all_links)
                triads_count_mat = link_triangles_mat_inv(n, all_links, counts)

                kwargs = (kwargs..., triads_count_mat = triads_count_mat)
                kwargs_dict[:triads_count_mat] = triads_count_mat
            end
        elseif ode_fun_name in ["Heider9!", "Heider92!"]
            if !haskey(kwargs_dict, :link_indices)
                link_indices = findall(triu(all_links_mat, 1)[:] .> 0)

                kwargs = (kwargs..., link_indices = link_indices)
                kwargs_dict[:link_indices] = link_indices
            end

            if !haskey(kwargs_dict, :link_pairs)
                triads_around_links_dict =
                    get_triangles_around_links(kwargs_dict[:all_triads])
                all_links = get_links_in_triads(kwargs_dict[:all_triads])

                link_pairs = get_triangles_around_links(triads_around_links_dict, all_links)

                kwargs = (kwargs..., link_pairs = link_pairs)
                kwargs_dict[:link_pairs] = link_pairs
            end

            if !haskey(kwargs_dict, :link_pairs_triad_cnt)
                link_pairs_triad_cnt = [length(link) for link in kwargs_dict[:link_pairs]]

                kwargs = (kwargs..., link_pairs_triad_cnt = link_pairs_triad_cnt)
                kwargs_dict[:link_pairs_triad_cnt] = link_pairs_triad_cnt
            end
        end

    end

    if ~isempty(all_links_mat) && ~isempty(specified_division)
        rl_weights = init_balanced_relations(n, specified_division)
    else
        rl_weights = init_random_balanced_relations(n, larger_size)
    end

    if ~isempty(all_links_mat)
        rl_weights .*= all_links_mat
    end

    r, filename = initialize_file(
        n,
        attr,
        gammas,
        maxtime,
        ode_fun_name,
        files_folder,
        filename_prefix,
    )

    # time prepariation
    if disp_more_every != 0
        time_disp = time()
    end

    if save_each != 0
        time_save = time()
    end

    firstline = 1
    realization_counter = 0

    art_attr = [ones(larger_size, 1); -ones(n - larger_size, 1)]
    for i = 1:length(gammas)
        gamma1 = gammas[i]

        #zeroing data arrays
        HB = zeros(zmax)
        HB_x = zeros(zmax)
        BR = zeros(zmax) #ratio of balanced triads
        HB_attr = zeros(zmax)
        # HB_only_weights = zeros(zmax);
        paradise = zeros(zmax)
        hell = zeros(zmax)
        Deltas = zeros(4, zmax)
        weak_balance_in_complete_graph = zeros(zmax)
        local_polarization = zeros(zmax)
        global_polarization = zeros(zmax)
        initial_neg_links_count = zeros(zmax)
        if isempty(all_links_mat)
            initial_neg_links_count .= ones(zmax) * larger_size * (n - larger_size)
        else
            initial_neg_links_count .= ones(zmax) .* sum(rl_weights .< 0)
            if issymmetric(all_links_mat)
                initial_neg_links_count ./= 2
            end
        end
        links_destab_changed = zeros(4, zmax)

        sim = zeros(zmax)
        x_attr_sim = zeros(zmax)
        stab = zeros(zmax)
        times = zeros(zmax)

        for rep = 1:zmax
            if !isempty(all_links_mat) && isempty(specified_division) # if this is true,then in each rep new division should be generated. 
                init_random_balanced_relations!(rl_weights, n, larger_size; art_attr = art_attr)
                rl_weights .*= all_links_mat

                initial_neg_links_count[zmax] = sum(rl_weights .< 0)
            end
            val0_attr = get_attributes(attr, n)
            al_weights = get_attribute_layer_weights(attr, val0_attr)

            if isempty(all_links_mat)
                (pos_destab, neg_destab) =
                    get_destabilized_links_count(rl_weights, al_weights, gamma1)
            elseif ode_fun in [Heider72!, Heider722!]
                (pos_destab, neg_destab) = get_destabilized_links_count(
                    rl_weights,
                    al_weights,
                    gamma1,
                    kwargs_dict[:triads_count_mat],
                )
            elseif ode_fun in [Heider9!, Heider92!]
                (pos_destab, neg_destab) = get_destabilized_links_count(
                    rl_weights[kwargs_dict[:link_indices]],
                    al_weights[kwargs_dict[:link_indices]],
                    gamma1,
                    kwargs_dict[:link_pairs],
                    kwargs_dict[:link_pairs_triad_cnt],
                )
            end
            links_destab_changed[1, :] .= pos_destab
            links_destab_changed[2, :] .= neg_destab
            if pos_destab + neg_destab == 0 #no destabilization => no sense of further calculations
                HB[rep] = 1
                HB_x[rep] = 1
                BR[rep] = 1

                paradise[rep] = is_paradise(rl_weights, n)
                hell[rep] = 0
                if isempty(all_links_mat)
                    HB_attr[rep] = is_hb(al_weights, n)
                    # HB_only_weights = zeros(zmax);
                    Deltas[:, rep] = get_triad_counts(rl_weights, n)
                    x_attr_sim[rep] = get_similarity(rl_weights, al_weights, n)
                else
                    HB_attr[rep] =
                        get_balanced_ratio_not_complete(
                            sign.(al_weights .* all_links_mat),
                            length(kwargs_dict[:all_triads]),
                        ) == 1

                    Deltas[:, rep] = get_triad_counts(rl_weights, kwargs_dict[:all_triads])
                    x_attr_sim[rep] =
                        get_similarity2(rl_weights, al_weights, sum(sign.(rl_weights) .> 0))
                end
                local_polarization[rep] = get_local_polarization(Deltas[:, rep])
                global_polarization[rep] = 1 - paradise[rep]
                weak_balance_in_complete_graph[rep] = 1

                stab[rep] = 1
                times[rep] = 0

                links_destab_changed[3, rep] = 0
                links_destab_changed[4, rep] = 0
            else
                #simulation
                (ishb_sim_par, t, u, u0, xy_attr, sol) = calc_heider_attr(
                    n,
                    attr,
                    gamma1,
                    maxtime,
                    ode_fun,
                    solver,
                    false,
                    rl_weights,
                    al_weights;
                    all_links_mat = all_links_mat,
                    kwargs...,
                )

                #work on results
                HB_x[rep],
                HB_attr[rep],
                x_attr_sim[rep],
                BR[rep],
                paradise[rep],
                hell[rep],
                Deltas[:, rep],
                weak_balance_in_complete_graph[rep],
                local_polarization[rep],
                global_polarization[rep] = ishb_sim_par

                if u isa Vector
                    links_destab_changed[3, rep] =
                        sum(u[u0[kwargs_dict[:link_indices]].>0] .< 0) #number of initial pos links that changed to negative
                    links_destab_changed[4, rep] =
                        sum(u[u0[kwargs_dict[:link_indices]].<0] .> 0) #number of initial neg links that changed to positive
                else
                    links_destab_changed[3, rep] = sum(u[u0.>0] .< 0) #number of initial pos links that changed to negative
                    links_destab_changed[4, rep] = sum(u[u0.<0] .> 0) #number of initial neg links that changed to positive
                end

                if t < maxtime #we have stability
                    HB[rep] = ishb_sim_par[1]
                    stab[rep] = 1
                end

                times[rep] = t
            end
            realization_counter += 1

            #checking time/repetitions and eventually informing about progress and saving
            if disp_more_every != 0
                if disp_more_every < time() - time_disp
                    time_disp = time()
                    # displaying
                    # display_res(gamma1, rep)
                    g = attr.g
                    v = r.attr_degeneracy
                    a = r.attr_name
                    # display_res(@ntuple(rep, a, n, g, v, gamma1, larger_size))
                    display_res(; rep, a, n, g, v, gamma1, larger_size)
                end
            end

            if save_each != 0
                if save_each < time() - time_save
                    time_save = time()
                    # partial saving
                    fields = (
                        HB,
                        HB_x,
                        HB_attr,
                        sim,
                        x_attr_sim,
                        BR,
                        paradise,
                        hell,
                        initial_neg_links_count,
                        links_destab_changed,
                        Deltas,
                        weak_balance_in_complete_graph,
                        local_polarization,
                        global_polarization,
                        stab,
                        times,
                        i,
                        rep,
                        firstline,
                    )
                    update_result!(r, fields)
                    save_result(r, filename)
                    save_result(r, filename, ext = "jld2")
                end
            end

            if disp_each != 0
                if disp_each <= realization_counter / zmax
                    realization_counter = 0
                    # displaying
                    g = attr.g
                    v = r.attr_degeneracy
                    a = r.attr_name
                    # display_res(@ntuple(rep, a, n, g, v, gamma1, larger_size))
                    display_res(; rep, a, n, g, v, gamma1, larger_size)
                end
            end
        end
        fields = (
            HB,
            HB_x,
            HB_attr,
            sim,
            x_attr_sim,
            BR,
            paradise,
            hell,
            initial_neg_links_count,
            links_destab_changed,
            Deltas,
            weak_balance_in_complete_graph,
            local_polarization,
            global_polarization,
            stab,
            times,
            i,
            zmax,
            firstline,
        )
        update_result!(r, fields)
        firstline += 2
    end

    save_result(r, filename)
    save_result(r, filename, ext = "jld2")

    return r
end
export using_heider_attr_destab

# end
