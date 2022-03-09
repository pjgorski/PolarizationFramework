# This file contains main function to generate network ans solve differential equations 
# and to run the series of repetitions of given model. 
# Apart from that, functions also make some initial and final states analysis and save the results to Matlab file. 
# Code is optimized

# module MeasureHeiderBalance

using DifferentialEquations
using Plots
using LinearAlgebra
using Dates

# function to measure HB in the case of a single-layered system with attributes
# n - number of agents
# gamma - coupling strength
# maxtime - maximal time of calculating differential equations
# ode_fun - function calculating derivatives
# solver - function solving different equations
# show_plot - parameter whether the solution will be plotted. If not transition values are not saved. 
# attr - attribute type; contains all information about attributes
# input is optional extra parameters with initial conditions of RL and AL (u0, xy_attr). 
# If input is empty, random initial conditions are generated. 
# If not, the first value (possibly empty) should contain initial conditions for RL. 
# The 2nd value of input may contain initial conditions for AL. 
# 
# Returns the tuple of: 
# ishb_sim_par - array of Boolean values [is it balanced state, is AL in balanced state,
        # similarity between RL and AL, balanced ratio of RL, is it paradise state,
        # is it hell state, counts od different triads (Deltas: Delta0, Delta1, Delta2, Delta3), 
        # is it weakly balanced state],
# time of simulation to reach stable solution (or maxtime if not reached),
# values of RL weights at the end of simulation, 
# initial conditions of RL, 
# initial conditions of AL, 
# whole sol (if transition values are not saved, see show_plot, then they are not here).
function calc_curheider_attr(n, gamma, maxtime::Float64, ode_fun, solver,
        show_plot::Bool, attr::AbstractAttributes, input...)
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
    Deltas = get_triad_counts(sol.u[end], n)
    weak_balance_in_complete_graph = Deltas[1 + 1] == 0
    ishb_sim_par = [is_hb(sol.u[end], n), is_hb(xy_attr, n),
        get_similarity(sol.u[end], xy_attr, n),
        get_balanced_ratio(sol.u[end], n),
        is_paradise(sol.u[end], n),
        is_hell(sol.u[end], n),
        Deltas, weak_balance_in_complete_graph]

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

# Function using `calc_curheider_attr` function to simulate a number of repetitions
# of solving the system having given parameters:
# n, betas are number of nodes (scalar) and betas (array) used
#  zmax - number of repetitions
#  disp_each - display status every how many realizations (if 0 then do
#      not display)
#  disp_more_every - display more detailed status (n, betas) every
#      certain number of sec (if 0, then do not display)
#  save_each - how often (in sec) should results be saved (if 0, then save
#      at the end)
# %
#  Last two arguments work like that, that if the specified time has past
#  then sth is displayed/saved. But it doesn't mean exact time of action.
#  When the specified time has past means when the current function
#  (heider function with parameters) finishes.


# function using calc_curheider (SINGLE LAYER!) function to simulate a number of repetitions:
# ode_fun_name should be a string, that is later converted into the function with the same name.


# function using calc_curheider (SINGLE LAYER!) function to simulate a number of repetitions:
# ode_fun_name should be a string, that is later converted into the function with the same name.
# start_2nd_iter_from is whether iterating other gammas should start from "1" or "i".



# function using calc_curheider_attr (SINGLE LAYER!) function to simulate a number of repetitions:
# ode_fun_name should be a string, that is later converted into the function with the same name.
function using_curheider_attr(n, g, gammas, zmax, maxtime, ode_fun_name::String,
    input_fun_name::String, attr::AbstractAttributes, disp_each, disp_more_every, save_each, filename)

    ode_fun = getfield(HBUtil, Symbol(ode_fun_name))
    if ode_fun_name == "Heider6!"
        solver =AutoTsit5(Rodas5(autodiff = false))
    else
        solver =AutoTsit5(Rodas5(autodiff = false))
    end

    if any(input_fun_name .== ["", "random", "rand"]) #do nothing
        case = 1
        n_start = ()
        #input_fun(n,g) = ()
    elseif input_fun_name == "init_balanced"
        input_fun2 = getfield(HBUtil, Symbol(input_fun_name))
        case = 2
        n_start = ()
        #input_fun(n,g) = ((), input_fun2(n,g))
        #input = ((), init_balanced(n,g))
    elseif input_fun_name == "init_random_balanced_enhanced"
        input_fun2 = getfield(HBUtil, Symbol(input_fun_name))
        n_start = find_nstart(attr, n)
        println("Started from n_start=$n_start.")
        case = 3
    end

    function input_fun(case, n_start=())
        if case == 1
            return ()
        elseif case == 2
            return ((), input_fun2(attr, n))
        elseif case == 3
            return ((), input_fun2(attr, n,n_start))
        end
    end

    r = CurResult(n, g, gammas, maxtime, ode_fun_name, attr);
    if isempty(filename)
        filename = "data/attr_curheider_results";
    else
        filename = "data/"*filename;
    end
    filename = get_filename(filename, ".mat");
    save_result(r, filename); #''allocating'' place

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
        HB_only_weights = zeros(zmax);
        paradise = zeros(zmax)
        hell = zeros(zmax)
        Deltas = zeros(4, zmax)
        weak_balance_in_complete_graph = zeros(zmax)
        initial_neg_links_count = zeros(zmax)
        links_destab_changed = zeros(4, zmax)

        sim = zeros(zmax);
        x_attr_sim = zeros(zmax);
        stab = zeros(zmax);
        times = zeros(zmax);
        for rep in 1:zmax
            #simulation
            (ishb_sim_par, t, u, u0, xy_attr, sol) =
                calc_curheider_attr(n, g, gamma1, maxtime, ode_fun, solver, false,
                    attr, input_fun(case, n_start)...)
            realization_counter += 1

            #work on results
            HB_x[rep], HB_attr[rep], x_attr_sim[rep], BR[rep],
                paradise[rep], hell[rep], Deltas[:, rep],
                weak_balance_in_complete_graph[rep] = ishb_sim_par

            initial_neg_links_count[rep] = sum(u0 .< 0)
            links_destab_changed[3, rep] = sum(u[u0.>0].<0) #number of initial pos links that changed to negative
            links_destab_changed[4, rep] = sum(u[u0.<0].>0) #number of initial neg links that changed to positive
            # if !isempty(sol2)
            #     time_inf[rep] = 1
            #     t = maxtime # in the case of :MaxIters, t is smaller than maxtime,
            #                 # which can lead to wrong results
            # end

            if t < maxtime #we have stability
                HB[rep] = ishb_sim_par[1]
                stab[rep] = 1
            end

            #sim[rep] = get_layer_similarity(u, n);
            times[rep] = t;

            #checking time/repetitions and eventually informing about progress and saving
            if disp_more_every != 0
                if disp_more_every < time() - time_disp
                    time_disp = time();
                    # displaying
                    display_res(gamma1, rep)
                end
            end

            if save_each != 0
                if save_each < time() - time_save
                    time_save = time();
                    # partial saving
                    fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
                        paradise, hell, initial_neg_links_count,
                        links_destab_changed, Deltas, weak_balance_in_complete_graph,
                        stab, times, i, rep, firstline);
                    update_result!(r, fields);
                    save_result(r, filename);
                end
            end

            if disp_each != 0
                if disp_each <= realization_counter
                    realization_counter = 0;
                    # displaying
                    display_res(gamma1, rep)
                end
            end
        end
        fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
            paradise, hell, initial_neg_links_count,
            links_destab_changed, Deltas, weak_balance_in_complete_graph,
            stab, times, i, zmax, firstline);
        update_result!(r, fields);
        firstline+=2;
    end

    save_result(r, filename);

    return r;
end
export using_curheider_attr

# function using calc_curheider_attr (SINGLE LAYER!) function to simulate
# the outcome of a destabilization of a balanced system.
# ode_fun_name should be a string, that is later converted into the function with the same name.
function using_curheider_attr_destab(n, g, gammas, larger_size::Int, zmax,
    maxtime, ode_fun_name::String,
    attr::AbstractAttributes,
    disp_each, disp_more_every, save_each, filename)

    ode_fun = getfield(HBUtil, Symbol(ode_fun_name))
    if ode_fun_name == "Heider6!"
        solver =AutoTsit5(Rodas5(autodiff = false))
    else
        solver =AutoTsit5(Rodas5(autodiff = false))
    end

    # if any(input_fun_name .== ["", "random", "rand"]) #do nothing
    #     case = 1
    #     n_start = ()
    #     #input_fun(n,g) = ()
    # elseif input_fun_name == "init_balanced"
    #     input_fun2 = getfield(HBUtil, Symbol(input_fun_name))
    #     case = 2
    #     n_start = ()
    #     #input_fun(n,g) = ((), input_fun2(n,g))
    #     #input = ((), init_balanced(n,g))
    # elseif input_fun_name == "init_random_balanced_enhanced"
    #     input_fun2 = getfield(HBUtil, Symbol(input_fun_name))
    #     n_start = find_nstart(attr, n)
    #     println("Started from n_start=$n_start.")
    #     case = 3
    # elseif input_fun_name == "init_random_balanced_relations"
    #     input_fun2 = getfield(HBUtil, Symbol(input_fun_name))
    #     case = 4
    #     n_start = ()
    # end
    # function input_fun(case, n_start=())
    #     if case == 1
    #         return ()
    #     elseif case == 2
    #         return ((), input_fun2(attr, n))
    #     elseif case == 3
    #         return ((), input_fun2(attr, n,n_start))
    #     elseif case == 4
    #         return (input_fun2(n, larger_size))
    #     end
    # end

    rl_weights = init_random_balanced_relations(n, larger_size)

    # val0_attr = get_attributes(attr, n);
    # al_weights = get_attribute_layer_weights(attr, val0_attr);
    #
    # (pos_destab, neg_destab, gammas2) = get_destabilized_links_count(rl, al)

    r = CurResult(n, g, gammas, maxtime, ode_fun_name, attr);
    if isempty(filename)
        filename = "data/attr_curheider_destab_results";
    else
        filename = "data/"*filename;
    end
    filename = get_filename(filename, ".mat");
    save_result(r, filename); #''allocating'' place

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
            if pos_destab + neg_destab == 0 #no destabilization
                HB[rep] = 1
                HB_x[rep] = 1
                BR[rep] = 1
                HB_attr[rep] = is_hb(al_weights,n)
                # HB_only_weights = zeros(zmax);
                paradise[rep] = is_paradise(rl_weights, n)
                hell[rep] = 0
                Deltas[:, rep] = get_triad_counts(rl_weights, n)
                weak_balance_in_complete_graph[rep] = 1

                x_attr_sim[rep] = get_similarity(rl_weights, al_weights, n);
                stab[rep] = 1;
                times[rep] = 0

                links_destab_changed[3, rep] = 0
                links_destab_changed[4, rep] = 0
            else
                #simulation
                (ishb_sim_par, t, u, u0, xy_attr, sol) =
                    calc_curheider_attr(n, g, gamma1, maxtime, ode_fun, solver, false,
                        attr, rl_weights, al_weights)

                #work on results
                HB_x[rep], HB_attr[rep], x_attr_sim[rep], BR[rep],
                    paradise[rep], hell[rep], Deltas[:, rep],
                    weak_balance_in_complete_graph[rep] = ishb_sim_par

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
                    display_res(gamma1, rep)
                end
            end

            if save_each != 0
                if save_each < time() - time_save
                    time_save = time();
                    # partial saving
                    fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
                        paradise, hell, initial_neg_links_count,
                        links_destab_changed, Deltas, weak_balance_in_complete_graph,
                        stab, times, i, rep, firstline);
                    update_result!(r, fields);
                    save_result(r, filename);
                end
            end

            if disp_each != 0
                if disp_each <= realization_counter
                    realization_counter = 0;
                    # displaying
                    display_res(gamma1, rep)
                end
            end
        end
        fields = (HB, HB_x, HB_attr, sim, x_attr_sim, BR,
            paradise, hell, initial_neg_links_count,
            links_destab_changed, Deltas, weak_balance_in_complete_graph,
            stab, times, i, zmax, firstline);
        update_result!(r, fields);
        firstline+=2;
    end

    save_result(r, filename);

    return r;
end
export using_curheider_attr_destab

function test()
    ishb, t, u, u0, sol = calc_heider(45, 0., 1000., true);
end

function test2()
    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider(3, 5, -1.5, 2000., Heider6!, AutoTsit5(Rodas5(autodiff = false)), false)

    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider(3, 5, -1.5, 2000., Heider6!, AutoTsit5(Rodas5(autodiff = false)), false, ())

    (ishb_sim, t, u, u0, xy_attr, sol) = calc_curheider(3, 5, -1.5, 2000., Heider6!, AutoTsit5(Rodas5(autodiff = false)), false, (), init_balanced(3,5))
end
export test2

function test3()
    n = 3
    g = 3
    v = 4
    gamma = 0.5
    maxtime = 1000.
    ode_fun = Heider7!
    solver = AutoTsit5(Rodas5(autodiff = false))
    show_plot = false
    attr = OrderedAttributes(g, 0.5, v)
    input = ()
    return n, g, v, gamma, maxtime, ode_fun, solver, show_plot, attr, input
end
export test3
# n, g, v, gamma, maxtime, ode_fun, solver, show_plot, attr, input = test3()
# (ishb_sim_par, t, u, u0, xy_attr, sol) = calc_curheider_attr(n, g, gamma, maxtime, ode_fun, solver, show_plot, attr, input)

# end
