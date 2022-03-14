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
function calc_curheider_attr(n::Int, attr::AbstractAttributes, gamma::Float64, maxtime::Float64, ode_fun::Function, solver,
        show_plot::Bool, input...)
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
        u0 = triu((rand(n,n)*2).-1,1)
    end
    if need_init_xy
        val0_attr = get_attributes(attr, n);
        xy_attr = get_attribute_layer_weights(attr, val0_attr);
    end

    # help variable
    mask = triu(trues(size(u0)),1)

    condition_here(u, t, integrator) = condition2(u, t, integrator, mask)
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition_here,affect!)

    #ode and model parameters
    tspan = (0.0,maxtime)

    lay1mul = zeros(n,n);
    x_sim = zeros(n,n);

    p = (n, gamma.*xy_attr, lay1mul, x_sim, mask)

    #solving
    prob = ODEProblem(ode_fun,u0,tspan,p)
    sol = solve(prob,solver,reltol=1e-6,abstol=1e-12, callback=cb,
        isoutofdomain = (u,p,t)->any(x->abs.(x)>=1,u), save_everystep=show_plot)

    #estimating output
    paradise = is_paradise(sol.u[end], n)
    Deltas = get_triad_counts(sol.u[end], n)
    weak_balance_in_complete_graph = Deltas[1 + 1] == 0
    local_polarization = get_local_polarization(Deltas)
    global_polarization = weak_balance_in_complete_graph && (!paradise)
    ishb_sim_par = [is_hb(sol.u[end], n), is_hb(xy_attr, n),
        get_similarity(sol.u[end], xy_attr, n),
        get_balanced_ratio(sol.u[end], n),
        paradise,
        is_hell(sol.u[end], n),
        Deltas, weak_balance_in_complete_graph, local_polarization, global_polarization]

    if show_plot

        h = plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
             xaxis="Time (t)",yaxis="weights(t)", ylim=[-1,+1],
             vars = reshape(1:n^2,n,n)[mask]); legend=false
        display(h);
        #gui()
    end

    return (ishb_sim_par, sol.t[end], sol.u[end], u0, xy_attr, sol)
end
export calc_curheider_attr

# Function creating results file. It parses parameters creating a unique file name and saves initially this file. 
# For parameters description see `using_curheider_attr`. 
# Returns a Tuple with:
#   object of `Result` DataType,
#   `filename`
#   filename without extension
function initialize_file(n::Int, attr::AbstractAttributes, gammas::Vector{Float64}, maxtime::Float64, ode_fun_name::String,
    files_folder::String, filename_prefix::String)

    r = Result(n, attr, gammas, maxtime, ode_fun_name)

    prefix = filename_prefix*Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    file_params = savename(prefix, r, "mat", sort=false)
    file_params = replace(file_params, "attr_degeneracy"=>"v")

    filename = projectdir(files_folder, file_params)
    save_result(r, filename); #''allocating'' place
    return r, filename
end

# Function using `calc_curheider_attr` function to simulate a number of repetitions
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
# files_folder - folder name inside the project the results should be saved (default "data")
# filename_prefix - string that should start each simulation results' file (default "")
# 
# `disp_more_every` and `save_each` work like that, that if the specified time has past
# then sth is displayed/saved. But it doesn't mean exact time of action.
# When the specified time has past means when the current function
# (heider function with parameters) finishes.
#  
# Returns `Result` object. 
function using_curheider_attr(n::Int, attr::AbstractAttributes, gammas::Vector{Float64}, zmax::Int, maxtime::Float64, ode_fun_name::String;
    disp_each = 0.5, disp_more_every = 600, save_each = 600, files_folder::String = "data", filename_prefix::String = "")

    ode_fun = getfield(PolarizationFramework, Symbol(ode_fun_name))
    solver = AutoTsit5(Rodas5(autodiff = false))

    r, filename = initialize_file(n, attr, gammas, maxtime, ode_fun_name, files_folder, filename_prefix)

    # time preparation
    if disp_more_every != 0
        time_disp = time();
    end

    if save_each != 0
        time_save = time();
    end

    firstline = 1;
    realization_counter = 0;
    for i in 1:length(gammas)
        gamma1 = gammas[i];

        #zeroing data arrays
        HB = zeros(zmax);
        HB_x = zeros(zmax);
        BR = zeros(zmax); #ratio of balanced triads
        HB_attr = zeros(zmax);
        HB_only_weights = zeros(zmax);
        paradise = zeros(zmax)
        hell = zeros(zmax)
        Deltas = zeros(4, zmax)
        weak_balance_in_complete_graph = zeros(zmax)
        local_polarization = zeros(zmax)
        global_polarization = zeros(zmax)
        initial_neg_links_count = zeros(zmax)
        links_destab_changed = zeros(4, zmax)

        sim = zeros(zmax);
        x_attr_sim = zeros(zmax);
        stab = zeros(zmax);
        times = zeros(zmax);
        for rep in 1:zmax
            #simulation
            (ishb_sim_par, t, u, u0, xy_attr, sol) =
                calc_curheider_attr(n, attr, gamma1, maxtime, ode_fun, solver, false)
            realization_counter += 1

            #work on results
            HB_x[rep], HB_attr[rep], x_attr_sim[rep], BR[rep],
                paradise[rep], hell[rep], Deltas[:, rep],
                weak_balance_in_complete_graph[rep], local_polarization[rep], global_polarization[rep] = ishb_sim_par

            initial_neg_links_count[rep] = sum(u0 .< 0)
            links_destab_changed[3, rep] = sum(u[u0.>0].<0) #number of initial pos links that changed to negative
            links_destab_changed[4, rep] = sum(u[u0.<0].>0) #number of initial neg links that changed to positive

            if t < maxtime #we have stability
                HB[rep] = ishb_sim_par[1]
                stab[rep] = 1
            end

            times[rep] = t;

            #checking time/repetitions and eventually informing about progress and saving
            if disp_more_every != 0
                if disp_more_every < time() - time_disp
                    time_disp = time();
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
                    time_save = time();
                    # partial saving
                    fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
                        paradise, hell, initial_neg_links_count,
                        links_destab_changed, Deltas, weak_balance_in_complete_graph, local_polarization, global_polarization,
                        stab, times, i, rep, firstline);
                    update_result!(r, fields);
                    save_result(r, filename);
                    save_result(r, filename, ext = "jld2");
                end
            end

            if disp_each != 0
                if disp_each <= realization_counter / zmax
                    realization_counter = 0;
                    # displaying
                    g = attr.g
                    v = r.attr_degeneracy
                    a = r.attr_name
                    # display_res(@ntuple(rep, a, n, g, v, gamma1))
                    display_res(; rep, a, n, g, v, gamma1)
                end
            end
        end
        fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
            paradise, hell, initial_neg_links_count,
            links_destab_changed, Deltas, weak_balance_in_complete_graph, local_polarization, global_polarization,
            stab, times, i, zmax, firstline);
        update_result!(r, fields);
        firstline+=2;
    end

    save_result(r, filename);
    save_result(r, filename, ext = "jld2");

    return r;
end
export using_curheider_attr

# Function using `calc_curheider_attr` function to simulate a number of repetitions
# of solving the system having given parameters and almost a balanced RL.
# Agents are divided into two groups and form positive links (+0.99) inside those groups
# and negative (-0.99) to agents outside their group. 
# Thus, this simulates destabilization of the polarized, initial state. 
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
# files_folder - folder name inside the project the results should be saved (default "data")
# filename_prefix - string that should start each simulation results' file (default "")
# 
# `disp_more_every` and `save_each` work like that, that if the specified time has past
# then sth is displayed/saved. But it doesn't mean exact time of action.
# When the specified time has past means when the current function
# (heider function with parameters) finishes.
#  
# Returns `Result` object. 
function using_curheider_attr_destab(n::Int, attr::AbstractAttributes, gammas::Vector{Float64}, larger_size::Int, zmax::Int,
    maxtime::Float64, ode_fun_name::String;
    disp_each = 0.5, disp_more_every = 600, save_each = 600, files_folder::String = "data", filename_prefix::String = "")

    ode_fun = getfield(PolarizationFramework, Symbol(ode_fun_name))
    solver = AutoTsit5(Rodas5(autodiff = false))

    rl_weights = init_random_balanced_relations(n, larger_size)

    r, filename = initialize_file(n, attr, gammas, maxtime, ode_fun_name, files_folder, filename_prefix)

    # time prepariation
    if disp_more_every != 0
        time_disp = time();
    end

    if save_each != 0
        time_save = time();
    end

    firstline = 1;
    realization_counter = 0;
    for i in 1:length(gammas)
        gamma1 = gammas[i];

        #zeroing data arrays
        HB = zeros(zmax);
        HB_x = zeros(zmax);
        BR = zeros(zmax); #ratio of balanced triads
        HB_attr = zeros(zmax);
        # HB_only_weights = zeros(zmax);
        paradise = zeros(zmax)
        hell = zeros(zmax)
        Deltas = zeros(4, zmax)
        weak_balance_in_complete_graph = zeros(zmax)
        local_polarization = zeros(zmax)
        global_polarization = zeros(zmax)
        initial_neg_links_count = ones(zmax)*larger_size*(n-larger_size)
        links_destab_changed = zeros(4, zmax)

        sim = zeros(zmax);
        x_attr_sim = zeros(zmax);
        stab = zeros(zmax);
        times = zeros(zmax);

        for rep in 1:zmax
            val0_attr = get_attributes(attr, n);
            al_weights = get_attribute_layer_weights(attr, val0_attr);

            (pos_destab, neg_destab) = get_destabilized_links_count(rl_weights, al_weights, gamma1)
            links_destab_changed[1,:] .= pos_destab
            links_destab_changed[2,:] .= neg_destab
            if pos_destab + neg_destab == 0 #no destabilization => no sense of further calculations
                HB[rep] = 1
                HB_x[rep] = 1
                BR[rep] = 1
                HB_attr[rep] = is_hb(al_weights,n)
                # HB_only_weights = zeros(zmax);
                paradise[rep] = is_paradise(rl_weights, n)
                hell[rep] = 0
                Deltas[:, rep] = get_triad_counts(rl_weights, n)
                local_polarization[rep] = get_local_polarization(Deltas[:, rep])
                global_polarization[rep] = 1-paradise[rep]
                weak_balance_in_complete_graph[rep] = 1

                x_attr_sim[rep] = get_similarity(rl_weights, al_weights, n);
                stab[rep] = 1;
                times[rep] = 0

                links_destab_changed[3, rep] = 0
                links_destab_changed[4, rep] = 0
            else
                #simulation
                (ishb_sim_par, t, u, u0, xy_attr, sol) =
                    calc_curheider_attr(n, attr, gamma1, maxtime, ode_fun, solver, false,
                        rl_weights, al_weights)

                #work on results
                HB_x[rep], HB_attr[rep], x_attr_sim[rep], BR[rep],
                    paradise[rep], hell[rep], Deltas[:, rep],
                    weak_balance_in_complete_graph[rep], local_polarization[rep], global_polarization[rep] = ishb_sim_par

                links_destab_changed[3, rep] = sum(u[u0.>0].<0) #number of initial pos links that changed to negative
                links_destab_changed[4, rep] = sum(u[u0.<0].>0) #number of initial neg links that changed to positive

                if t < maxtime #we have stability
                    HB[rep] = ishb_sim_par[1]
                    stab[rep] = 1
                end

                times[rep] = t;
            end
            realization_counter += 1

            #checking time/repetitions and eventually informing about progress and saving
            if disp_more_every != 0
                if disp_more_every < time() - time_disp
                    time_disp = time();
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
                    time_save = time();
                    # partial saving
                    fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
                        paradise, hell, initial_neg_links_count,
                        links_destab_changed, Deltas, weak_balance_in_complete_graph, local_polarization, global_polarization,
                        stab, times, i, rep, firstline);
                    update_result!(r, fields);
                    save_result(r, filename);
                    save_result(r, filename, ext = "jld2");
                end
            end

            if disp_each != 0
                if disp_each <= realization_counter / zmax
                    realization_counter = 0;
                    # displaying
                    g = attr.g
                    v = r.attr_degeneracy
                    a = r.attr_name
                    # display_res(@ntuple(rep, a, n, g, v, gamma1, larger_size))
                    display_res(; rep, a, n, g, v, gamma1, larger_size)
                end
            end
        end
        fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
            paradise, hell, initial_neg_links_count,
            links_destab_changed, Deltas, weak_balance_in_complete_graph, local_polarization, global_polarization,
            stab, times, i, zmax, firstline);
        update_result!(r, fields);
        firstline+=2;
    end

    save_result(r, filename);
    save_result(r, filename, ext = "jld2");

    return r;
end
export using_curheider_attr_destab

# end
