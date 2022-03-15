# This file contains type Attributes and some functions related with it.
# The idea of Attributes class structure is that the instances do not contain data
# (i.e. attributes), but they define how the attributes are supposed to be treated
# (created, analyzed). The attributes itself should be kept in seperate Array objects.
#
# Attributes types inherit over AbstractAttributes.
# A few attribute types are defined (see PHD thesis, PJ Górski, p.94).
# The most basic is BinaryAttributes. Agents possess a vector of +-1, e.g. [+1, -1, -1].
# OrderedAttributes, UnorderedAttributes, UnorderedPositiveAttributes are for agents
# with vectors of integer attributes of 1, 2, 3, ..., e.g. [1, 4, 2, 3, 1].
# (The max value is defined.)
# OrderedAttributes assume attributes can be placed in certain order.
# UnorderedAttributes and UnorderedPositiveAttributes assume attributes either are the same
# or different (there is no many different levels of difference between different attributes).
# UnorderedAttributes assume different attributes have negative contribution to attribute's distance.
# UnorderedPositiveAttributes assume different attributes have neutral contribution to attribute's distance.
#
# To run requires to add to LOAD_PATH the filepath to this file, e.g.
# > push!(LOAD_PATH, pwd())

# module Attributes

using LinearAlgebra

abstract type AbstractAttributes end
export AbstractAttributes

struct BinaryAttributes <: AbstractAttributes
    g::Int
end
export BinaryAttributes

struct OrderedAttributes <: AbstractAttributes
    g::Int# number of attributes
    threshold::Float64#threshold when to accept the nodes are friends
    v::Int# number of different attribute values
end
export OrderedAttributes

struct UnorderedAttributes <: AbstractAttributes
    g::Int# number of attributes
    threshold::Float64#threshold when to accept the nodes are friends
    v::Int# number of different attribute values
end
export UnorderedAttributes

struct UnorderedPositiveAttributes <: AbstractAttributes
    g::Int# number of attributes
    threshold::Float64#threshold when to accept the nodes are friends
    v::Int# number of different attribute values
end
export UnorderedPositiveAttributes

#Returns the Array of attributes of type BinaryAttributes (attributes are +-1)
get_attributes(b::BinaryAttributes, n::Int) = sign.(rand(n, b.g) .- 0.5)
export get_attributes

#Returns the Array of attributes of type (Un)OrderedAttributes (attributes are 1, 2, ..., v)
get_attributes(b::AbstractAttributes, n::Int) = rand(1:b.v, n, b.g)
export get_attributes

#Returns weights of the AL in the case of BinaryAttributes
get_attribute_layer_weights(b::BinaryAttributes, attr::AbstractArray{<:Number}) =
    triu(attr * attr' / b.g, 1)
export get_attribute_layer_weights

function get_attribute_layer_weights(b::AbstractAttributes, attr::AbstractArray{<:Number})
    n = size(attr, 1)
    weights = zeros(n, n)
    for i = 1:n, j = (i+1):n
        if isa(b, OrderedAttributes)
            weights[i, j] =
                2(b.threshold - sum(abs.(attr[i, :] - attr[j, :])) / (b.v - 1) / b.g)
        elseif isa(b, UnorderedAttributes)
            weights[i, j] = 2(sum(attr[i, :] .== attr[j, :]) / b.g - b.threshold)
        elseif isa(b, UnorderedPositiveAttributes)
            weights[i, j] = sum(attr[i, :] .== attr[j, :]) / b.g
        end
    end
    return weights
end
export get_attribute_layer_weights

get_name(b::AbstractAttributes) = string(typeof(b))
export get_name

# Base.string(b::AbstractAttributes) = get_name(b)


get_threshold(b::AbstractAttributes) = b.threshold
get_threshold(b::BinaryAttributes) = 0.5
export get_threshold

get_degeneracy(b::AbstractAttributes) = b.v
get_degeneracy(b::BinaryAttributes) = 2
export get_degeneracy

# end
