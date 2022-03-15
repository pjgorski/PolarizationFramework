# file contains some utility functions. 

using Dates

# module MyUtil


# finds the first free filename (that would not be overwriting) with
# the expression: filename * "_$i", where i is integer starting from 1.
# adds extension at the end.
# returns proper filename
function get_filename(filename, extension)
    i = 1
    #if extension does not have a "." at the beginning we should add one.
    if !startswith(extension, ".")
        extension = "." * extension
    end

    #filename should not end with the extension.
    if endswith(filename, extension)
        filename = filename[1:(end-length(extension))]
    end

    while isfile(filename * "_$i" * extension)
        i += 1
    end
    filename * "_$i" * extension
end
export get_filename

#get k max values from array a
function maxk(a, k)
    b = partialsortperm(a, 1:k, rev = true)
    return collect(zip(b, a[b]))
end
export maxk

macro Name(arg)
    string(arg)
end

#function for displaying progress
function display_res(names, vars)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    str = ""
    for i = 1:1:length(names)
        str = str * "$(names[i])" * " = " * "$(vars[i])" * "; "
    end
    println("Just finished: $str" * "Time: $cur")
    return nothing
end

function display_res(ntuple::NamedTuple)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    str = savename(ntuple, connector = "; ", sort = false)

    println("Just finished: $str" * ". Time: $cur")
    return nothing
end

function display_res(; kwargs...)
    ntuple = NamedTuple(kwargs)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    str = string(ntuple)

    println("Just finished: $str. Time: $cur")
    return nothing
end

function display_res(gamma::Float64, z::Int)
    cur = Dates.format(now(), "HH:MM:SS dd.mm.yyyy")
    println("Just finished: z = $z; gamma = $gamma; " * "Time: $cur")
    return nothing
end

# Function for in-place creating a symmetrical array of array x. 
# x should be an array of size `nxn`. 
function Sym!(x_sim, x, n)
    for i = 1:n
        x_sim[i, i] = x[i, i]
        for j = (i+1):n
            x_sim[i, j] = x[i, j]
            x_sim[j, i] = x[i, j]
        end
    end
end
export Sym!
