# This file defines struct type for keeping simulation data. 
# It also defines functions for saving data in MAT (Matlab) file. 

# module Results

using DrWatson
using MAT
using Statistics
using Suppressor

# using Attributes

#This type is to store simulation results
struct Result
    #parameters
    n::UInt32 #number of nodes
    g::UInt32 #number of attributes
    gammas::Array{Float64,1} #gammas
    maxtime::Float64 #simulation stopping time
    ode_fun_name::String #used ode_function
    larger_size::Int # This is -1 in the case this parameter is not important (Scenario B) or it shows the size of larger group in the case of destabilization (Scenario A)

    attr_name::String #atribute type
    attr_threshold::Float64 #threshold (if any)
    attr_degeneracy::Int #number of possible values

    #true results
    HB::Array{Float64,1} #number of times HB and condition were achieved divided by number of repetitions
    HB_x::Array{Float64,1} #number of times HB was achieved (no matter the condition) divided by number of repetitions
    HB_attr::Array{Float64,1} #number of times HB in attribute similarity matrix was achieved divided by number of repetitions
    HB_only_weights::Array{Float64,1} #number of times HB was achieved but attribute similarity matrix is not in HB (divided by number of times HB was not achieved in ASM)
    BR::Array{Float64,1} #mean value of ratio of balanced triads
    paradise::Array{Float64,1} #number of times paradise was achieved divided by number of repetitions
    hell::Array{Float64,1} #number of times hell was achieved divided by number of repetitions
    weak_balance_in_complete_graph::Array{Float64,1} #number of times weak SB was achieved divided by number of repetitions
    local_polarization::Array{Float64,1} #local polarization measure
    global_polarization::Array{Float64,1} #number of times system was globally polarized

    initial_neg_links_count::Array{Float64,1} #mean number of initial negative links
    initial_neg_links_count_std::Array{Float64,1} #std of initial negative links
    links_destab_changed::Array{Float64,2} #array of mean numbers of positive links destabilized, negative links destab., positive links changes and neg links changed
    links_destab_changed_std::Array{Float64,2} #array of std of positive links destabilized, negative links destab., positive links changes and neg links changed
    Deltas::Array{Float64,2} # array with mean numbers of triad counts
    Deltas_std::Array{Float64,2} # std of triad counts

    sim::Array{Float64,1} #mean similarity of layers at the end of simulations
    x_attr_sim::Array{Float64,1} #mean similarity between weights and attributes at the end of simulations
    stab::Array{Float64,1} #mean stability at the end of simulations
    times::Array{Float64,1} #mean times to reach the end of simulation

    BR_std::Array{Float64,1} #std deviation of ratio of balanced triads
    local_polarization_std::Array{Float64,1} #std deviation of local polarization measure
    sim_std::Array{Float64,1} #std deviation of similarity of layers at the end of simulations
    x_attr_sim_std::Array{Float64,1} #std deviation of similarity of weights and attributes at the end of simulations
    times_std::Array{Float64,1} #std deviation of times to reach the end of simulation

    zmax::Array{UInt32,1} #number of repetitions

    data::Array{Float64,2} #all data listed in an array
    interpretation::Array{String,1} #interpretation of columns of above field
    # interpretation = ["N", "G", "gamma", "HB_in_attr", "larger_size", 
    # "zmax", "HB", "HB_x", "paradise", "hell", "weak_balance_in_complete_graph",
    # "sim", "sim_std",
    # "x_attr_sim", "x_attr_sim_std",
    # "BR", "BR_std", "LP", "LP_std", "GP", "stab",
    # "times", "times_std",
    # "pos_links_destab", "pos_links_destab_std", "neg_links_destab", "neg_links_destab_std",
    # "pos_links_changed", "pos_links_changed_std", "neg_links_changed", "neg_links_changed_std",
    # "initial_neg_links_count", "initial_neg_links_count_std",
    # "Delta0", "Delta0_std", "Delta1", "Delta1_std", "Delta2",
    # "Delta2_std", "Delta3", "Delta3_std"]

    # Following is the different way of keeping and later saving data. 
    # Previously (above) data was separated into AL-balanced and AL-unbalanced cases. Below it is not. 
    data2::Array{Any,2} #all data listed in an array
    interpretation2::Array{String,1} #interpretation of columns of above field
    # interpretation = ["N", "G", "gamma", "maxtime", "ode_fun_name", "larger_size", 
    # "attr", "attr_name", "attr_threshold",  "attr_degeneracy", 
    # "zmax", "HB", "HB_x", "paradise", "hell", "weak_balance_in_complete_graph",
    # "sim", "sim_std",
    # "x_attr_sim", "x_attr_sim_std",
    # "BR", "BR_std", "LP", "LP_std", "GP", "stab",
    # "times", "times_std",
    # "pos_links_destab", "pos_links_destab_std", "neg_links_destab", "neg_links_destab_std",
    # "pos_links_changed", "pos_links_changed_std", "neg_links_changed", "neg_links_changed_std",
    # "initial_neg_links_count", "initial_neg_links_count_std",
    # "Delta0", "Delta0_std", "Delta1", "Delta1_std", "Delta2",
    # "Delta2_std", "Delta3", "Delta3_std"]
end
export Result

# Constructor
function Result(
    n::Int,
    attr::AbstractAttributes,
    gammas::Vector{Float64},
    maxtime::Float64,
    ode_fun_name::String,
    larger_size::Int = -1,
)
    nb = length(gammas)

    HB = zeros(Float64, nb)
    HB_x = copy(HB)
    HB_attr = copy(HB)
    HB_only_weights = copy(HB)
    BR = copy(HB)
    paradise = copy(HB)
    hell = copy(HB)
    weak_balance_in_complete_graph = copy(HB)
    local_polarization = copy(HB)
    global_polarization = copy(HB)

    initial_neg_links_count = copy(HB)
    initial_neg_links_count_std = copy(HB)
    Deltas = zeros(Float64, 4, nb)
    links_destab_changed = copy(Deltas)

    sim = zeros(Float64, nb)
    x_attr_sim = copy(sim)
    stab = copy(sim)
    times = copy(sim)

    BR_std = copy(sim)
    local_polarization_std = copy(sim)
    links_destab_changed_std = copy(Deltas)
    Deltas_std = copy(Deltas)
    sim_std = copy(sim)
    x_attr_sim_std = copy(sim)
    times_std = copy(sim)

    zmax = zeros(UInt32, nb)

    interpretation = [
        "N",
        "G",
        "gamma",
        "HB_in_attr",
        "larger_size",
        "zmax",
        "HB",
        "HB_x",
        "paradise",
        "hell",
        "weak_balance_in_complete_graph",
        "sim",
        "sim_std",
        "x_attr_sim",
        "x_attr_sim_std",
        "BR",
        "BR_std",
        "LP",
        "LP_std",
        "GP",
        "stab",
        "times",
        "times_std",
        "pos_links_destab",
        "pos_links_destab_std",
        "neg_links_destab",
        "neg_links_destab_std",
        "pos_links_changed",
        "pos_links_changed_std",
        "neg_links_changed",
        "neg_links_changed_std",
        "initial_neg_links_count",
        "initial_neg_links_count_std",
        "Delta0",
        "Delta0_std",
        "Delta1",
        "Delta1_std",
        "Delta2",
        "Delta2_std",
        "Delta3",
        "Delta3_std",
    ]
    data = zeros(Float64, nb * 2, length(interpretation))

    interpretation2 = [
        "N",
        "G",
        "gamma",
        "maxtime",
        "ode_fun_name",
        "larger_size",
        "attr_name",
        "attr_threshold",
        "attr_degeneracy",
        "zmax",
        "HB",
        "HB_x",
        "paradise",
        "hell",
        "weak_balance_in_complete_graph",
        "sim",
        "sim_std",
        "x_attr_sim",
        "x_attr_sim_std",
        "BR",
        "BR_std",
        "LP",
        "LP_std",
        "GP",
        "stab",
        "times",
        "times_std",
        "pos_links_destab",
        "pos_links_destab_std",
        "neg_links_destab",
        "neg_links_destab_std",
        "pos_links_changed",
        "pos_links_changed_std",
        "neg_links_changed",
        "neg_links_changed_std",
        "initial_neg_links_count",
        "initial_neg_links_count_std",
        "Delta0",
        "Delta0_std",
        "Delta1",
        "Delta1_std",
        "Delta2",
        "Delta2_std",
        "Delta3",
        "Delta3_std",
    ]
    data2 = Array{Any,2}(undef, nb, length(interpretation2))
    data2[:, 1:4] .= 0
    data2[:, [5, 7]] .= ""
    data2[:, 6] .= 0
    data2[:, 7:end] .= 0

    Result(
        n,
        attr.g,
        gammas,
        maxtime,
        ode_fun_name,
        larger_size,
        get_name(attr),
        get_threshold(attr),
        get_degeneracy(attr),
        HB,
        HB_x,
        HB_attr,
        HB_only_weights,
        BR,
        paradise,
        hell,
        weak_balance_in_complete_graph,
        local_polarization,
        global_polarization,
        initial_neg_links_count,
        initial_neg_links_count_std,
        links_destab_changed,
        links_destab_changed_std,
        Deltas,
        Deltas_std,
        sim,
        x_attr_sim,
        stab,
        times,
        BR_std,
        local_polarization_std,
        sim_std,
        x_attr_sim_std,
        times_std,
        zmax,
        data,
        interpretation,
        data2,
        interpretation2,
    )
end
export Result

function Result(
    n::Int,
    g::Int,
    gammas::Vector{Float64},
    maxtime::Float64,
    ode_fun_name::String,
    larger_size::Int = -1,
)
    Result(n, BinaryAttributes(g), gammas, maxtime, ode_fun_name, larger_size)
end
export Result

function Result(n::Int, g::Int, gammas::Vector{Float64}, maxtime::Float64)
    Result(n, g, gammas, maxtime, "Heider5!")
end
export Result

function Result(n::Int, g::Int, gammas::Vector{Float64}, ode_fun_name::String)
    Result(n, g, gammas, 1000.0, ode_fun_name)
end
export Result

function Result(n::Int, g::Int, gammas::Vector{Float64})
    Result(n, g, gammas, 1000.0)
end
export Result

#save results to MATLAB or jld2 file
function save_result(res, filename; ext = "mat")
    if ext == "mat"
        if !endswith(filename, "mat")
            filename = filename * ".mat"
        end
        file = matopen(filename, "w")
        write(file, "result", res)
        close(file)
    elseif ext == "jld2"
        if endswith(filename, r"\.[0-9a-z]+$")
            filename = replace(filename, (r"\.[0-9a-z]+$" => ""))
        end
        for (i, gamma) in enumerate(res.gammas)
            # This was not yet simulated. End of saving
            if res.zmax[i] == 0
                break
            end

            #creating filename
            local filename_g
            @suppress filename_g = savename(filename, (gamma = gamma,))
            filename_g2 = filename_g * "_1.jld2"

            # data_row = res.data[2*i-1, :]
            data_row = res.data2[i, :]
            @suppress @tagsave(
                filename_g2,
                # Dict(res.interpretation .=> data_row)
                Dict(res.interpretation2 .=> data_row)
            )
            # filename_g2 = filename_g*"_2.jld2"
            # data_row = res.data[2*i, :]
            # @suppress @tagsave(
            #     filename_g2,
            #     Dict(res.interpretation .=> data_row)
            # )
        end
    end

end
export save_result

function read_result(filename, varname)
    if !endswith(filename, ".mat")
        filename *= ".mat"
    end
    file = matopen(filename)
    out = read(file, varname) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    out
end
export read_result

#function to update result fields based on elements of fields tuple
function update_result!(res::Result, fields)
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
    firstline = fields
    res.HB[i] = sum(HB) / rep
    res.HB_x[i] = sum(HB_x) / rep
    res.HB_attr[i] = sum(HB_attr) / rep
    res.HB_only_weights[i] = sum(HB[HB_attr.==0]) / ((1 - res.HB_attr[i]) * rep)
    res.BR[i] = sum(BR) / rep
    res.paradise[i] = sum(paradise) / rep
    res.hell[i] = sum(hell) / rep
    res.weak_balance_in_complete_graph[i] = sum(weak_balance_in_complete_graph) / rep
    res.local_polarization[i] = sum(local_polarization) / rep
    res.global_polarization[i] = sum(global_polarization) / rep

    res.initial_neg_links_count[i] = sum(initial_neg_links_count) / rep
    res.links_destab_changed[:, i] = sum(links_destab_changed, dims = 2) / rep
    res.Deltas[:, i] = sum(Deltas, dims = 2) / rep

    res.sim[i] = sum(sim) / rep
    res.x_attr_sim[i] = sum(x_attr_sim) / rep
    res.stab[i] = sum(stab) / rep
    res.times[i] = sum(times) / rep

    res.BR_std[i] = std(BR[1:rep])
    res.local_polarization_std[i] = std(local_polarization[1:rep])
    res.initial_neg_links_count_std[i] = std(initial_neg_links_count[1:rep])
    res.Deltas_std[:, i] = std(Deltas[:, 1:rep], dims = 2)
    res.links_destab_changed_std[:, i] = std(links_destab_changed[:, 1:rep], dims = 2)

    res.sim_std[i] = std(sim[1:rep])
    res.x_attr_sim_std[i] = std(x_attr_sim[1:rep])
    res.times_std[i] = std(times[1:rep])

    res.zmax[i] = rep

    #line when HB is in attribute matrix
    mask_true = (HB_attr .== 1)
    smt = sum(mask_true)
    mask_false = .!mask_true
    smf = rep - smt
    # try
    res.data[firstline, :] = [
        res.n,
        res.g,
        res.gammas[i],
        true,
        res.larger_size,
        smt,
        sum(HB[mask_true]) / smt,
        sum(HB_x[mask_true]) / smt,
        sum(paradise[mask_true]) / smt,
        sum(hell[mask_true]) / smt,
        sum(weak_balance_in_complete_graph[mask_true]) / smt,
        sum(sim[mask_true]) / smt,
        std(sim[1:rep][mask_true[1:rep]]),
        sum(x_attr_sim[mask_true]) / smt,
        std(x_attr_sim[1:rep][mask_true[1:rep]]),
        sum(BR[mask_true]) / smt,
        std(BR[1:rep][mask_true[1:rep]]),
        sum(local_polarization[mask_true]) / smt,
        std(local_polarization[1:rep][mask_true[1:rep]]),
        sum(global_polarization[mask_true]) / smt,
        sum(stab[mask_true]) / smt,
        sum(times[mask_true]) / smt,
        std(times[1:rep][mask_true[1:rep]]),
        sum(links_destab_changed[1, mask_true]) / smt,
        std(links_destab_changed[1, 1:rep][mask_true[1:rep]]),
        sum(links_destab_changed[2, mask_true]) / smt,
        std(links_destab_changed[2, 1:rep][mask_true[1:rep]]),
        sum(links_destab_changed[3, mask_true]) / smt,
        std(links_destab_changed[3, 1:rep][mask_true[1:rep]]),
        sum(links_destab_changed[4, mask_true]) / smt,
        std(links_destab_changed[4, 1:rep][mask_true[1:rep]]),
        sum(initial_neg_links_count[mask_true]) / smt,
        std(initial_neg_links_count[1:rep][mask_true[1:rep]]),
        sum(Deltas[1, mask_true]) / smt,
        std(Deltas[1, 1:rep][mask_true[1:rep]]),
        sum(Deltas[2, mask_true]) / smt,
        std(Deltas[2, 1:rep][mask_true[1:rep]]),
        sum(Deltas[3, mask_true]) / smt,
        std(Deltas[3, 1:rep][mask_true[1:rep]]),
        sum(Deltas[4, mask_true]) / smt,
        std(Deltas[4, 1:rep][mask_true[1:rep]]),
    ]
    # catch y
    #     display([firstline, rep])
    #     display([length(mask_true), sum(mask_true), smt])
    #     display([length(HB_x), sum(HB_x[mask_true]), smt])
    #     display([length(sim), sum(sim[mask_true]), smt])
    #     display([length(HB_x), sum(HB_x[mask_true]), smt])
    #     display([length(sim), sum(sim[mask_true]), smt])
    #     throw(y)
    # end
    #line when HB is not in attribute matrix
    res.data[firstline+1, :] = [
        res.n,
        res.g,
        res.gammas[i],
        false,
        res.larger_size,
        smf,
        sum(HB[mask_false]) / smf,
        sum(HB_x[mask_false]) / smf,
        sum(paradise[mask_false]) / smf,
        sum(hell[mask_false]) / smf,
        sum(weak_balance_in_complete_graph[mask_false]) / smf,
        sum(sim[mask_false]) / smf,
        std(sim[1:rep][mask_false[1:rep]]),
        sum(x_attr_sim[mask_false]) / smf,
        std(x_attr_sim[1:rep][mask_false[1:rep]]),
        sum(BR[mask_false]) / smf,
        std(BR[1:rep][mask_false[1:rep]]),
        sum(local_polarization[mask_false]) / smf,
        std(local_polarization[1:rep][mask_false[1:rep]]),
        sum(global_polarization[mask_false]) / smf,
        sum(stab[mask_false]) / smf,
        sum(times[mask_false]) / smf,
        std(times[1:rep][mask_false[1:rep]]),
        sum(links_destab_changed[1, mask_false]) / smf,
        std(links_destab_changed[1, 1:rep][mask_false[1:rep]]),
        sum(links_destab_changed[2, mask_false]) / smf,
        std(links_destab_changed[2, 1:rep][mask_false[1:rep]]),
        sum(links_destab_changed[3, mask_false]) / smf,
        std(links_destab_changed[3, 1:rep][mask_false[1:rep]]),
        sum(links_destab_changed[4, mask_false]) / smf,
        std(links_destab_changed[4, 1:rep][mask_false[1:rep]]),
        sum(initial_neg_links_count[mask_false]) / smf,
        std(initial_neg_links_count[1:rep][mask_false[1:rep]]),
        sum(Deltas[1, mask_false]) / smf,
        std(Deltas[1, 1:rep][mask_false[1:rep]]),
        sum(Deltas[2, mask_false]) / smf,
        std(Deltas[2, 1:rep][mask_false[1:rep]]),
        sum(Deltas[3, mask_false]) / smf,
        std(Deltas[3, 1:rep][mask_false[1:rep]]),
        sum(Deltas[4, mask_false]) / smf,
        std(Deltas[4, 1:rep][mask_false[1:rep]]),
    ]

    # "maxtime", "ode_fun_name", 
    # "attr", "attr_name", "attr_threshold",  "attr_degeneracy", 
    res.data2[i, :] = [
        res.n,
        res.g,
        res.gammas[i],
        res.maxtime,
        res.ode_fun_name,
        res.larger_size,
        res.attr_name,
        res.attr_threshold,
        res.attr_degeneracy,
        rep,
        sum(HB) / rep,
        sum(HB_x) / rep,
        sum(paradise) / rep,
        sum(hell) / rep,
        sum(weak_balance_in_complete_graph) / rep,
        sum(sim) / rep,
        std(sim[1:rep]),
        sum(x_attr_sim) / rep,
        std(x_attr_sim[1:rep]),
        sum(BR) / rep,
        std(BR[1:rep]),
        sum(local_polarization) / rep,
        std(local_polarization[1:rep]),
        sum(global_polarization) / rep,
        sum(stab) / rep,
        sum(times) / rep,
        std(times[1:rep]),
        sum(links_destab_changed[1, mask_false]) / rep,
        std(links_destab_changed[1, 1:rep]),
        sum(links_destab_changed[2, mask_false]) / rep,
        std(links_destab_changed[2, 1:rep]),
        sum(links_destab_changed[3, mask_false]) / rep,
        std(links_destab_changed[3, 1:rep]),
        sum(links_destab_changed[4, mask_false]) / rep,
        std(links_destab_changed[4, 1:rep]),
        sum(initial_neg_links_count) / rep,
        std(initial_neg_links_count[1:rep]),
        sum(Deltas[1, mask_false]) / rep,
        std(Deltas[1, 1:rep]),
        sum(Deltas[2, mask_false]) / rep,
        std(Deltas[2, 1:rep]),
        sum(Deltas[3, mask_false]) / rep,
        std(Deltas[3, 1:rep]),
        sum(Deltas[4, mask_false]) / rep,
        std(Deltas[4, 1:rep]),
    ]
    # "zmax", "HB", "HB1", "HB2", "HB12", "sim", "sim_std", "stab",
    # "times", "times_std"]
end
export update_result!

DrWatson.allaccess(::Result) = (:attr_name, :n, :g, :attr_degeneracy)
# end
