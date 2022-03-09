# This file contains functions related to analysing the structural balance dynamcis and measuring the state of the system,
# among others whether the state is in paradise, structurally balanced etc. 

using DifferentialEquations
using LinearAlgebra

# using Attributes

f_inner(v1, v2, p) = triu((v1.^2-.1).*(v1*v1 .- v2.*p[1]).*(-1),1)
f_inner!(v, v1, v2, p) = triu!(v.=(v1.^2-.1).*(v1*v1 .- v2.*p[1]).*(-1),1)

# ode function for HB with attributes in single layer
# This is an ode function, so it follows the correct structure accepted by `DifferentialEquations` package. 
# dx - change of values
# x - current values
# p - help variable to optimize speed by limiting memory allocations. It should contain 4 elements:
#     number of nodes, gamma_attr (array with gamma already multiplied by attribute similarity matrix),
#     two help arrays for storing partial calculations. 
function Heider5!(dx, x, p, t)
  n, gamma_attr, lay1mul, x_sim = p
  # Symmetric3d!(x)
  x_sim .= Symmetric(x)
  # v1 = view(x,:,:,1)
  # v2 = view(x,:,:,2)
  A_mul_B!(lay1mul, x_sim, x_sim)
  # A_mul_B!(lay2mul, v2, v2)
  dx .= triu((x_sim.^2 .-1).*(lay1mul./(n-2) .+ gamma_attr).*(-1),1)
  # dx[:,:,1] .= triu((v1.^2 .-1).*(lay1mul./(n-2) .+ v2.*beta[1]).*(-1),1)
  # dx[:,:,2] .= triu((v2.^2 .-1).*(lay2mul./(n-2) .+ v1.*beta[2]).*(-1),1)
end
export Heider5!

# Following function is the same as `Heider5!`, but it was made even more efficient. 
# Help variable `p` should contain 5 elements. 
# The 5th element is a `mask`, which allows quick choosing of array elements above the diagonal. 
# `mask` should be prepared beforehand. 
function Heider7!(dx, x, p, t)
  n, gamma_attr, lay1mul, x_sim, mask = p
  Sym!(x_sim,x,n)
  mul!(lay1mul, x_sim, x_sim)
  dx .= (x_sim.^2 .-1).*(lay1mul./(n-2) .+ gamma_attr).*(-1).*mask
end
export Heider7!

# Following function is the same as `Heider5!`, used for speed comparison purposes. 
function Heider8!(dx, x, p, t)
  n, gamma_attr, lay1mul, x_sim, mask = p
  Sym!(x_sim,x,n)
  A_mul_B!(lay1mul, x_sim, x_sim)
  dx .= triu((x_sim.^2 .-1).*(lay1mul./(n-2) .+ gamma_attr).*(-1),1)
end

# # EVENT function used by `DifferentialEquations` package. 
# function condition(u,t,integrator) # Event when event_f(u,t) == 0
#     fc = abs.(u[mask]).>0.99
#     if !all(fc)
#         return false
#     end
#     fc .= view((get_du(integrator).*u) .> 0, mask);

#     if all(fc)
#         return true
#     else return false
#     end
# end
# export condition

# EVENT function used by `DifferentialEquations` package. 
# It contains mask variable for correct assessing end of simulations
function condition2(u,t,integrator, mask) # Event when event_f(u,t) == 0
    #println("$get_du(integrator)")

    fc = abs.(u[mask]).>0.99
    if !all(fc)
        return false
    end
    fc .= ((get_du(integrator).*u)[mask] .> 0);

    if all(fc)
        return true
    else return false
    end
end
export condition2

# checks whether multiplication of weights in all triads is positive.
# It does not have to be close to +-1. Then returns true.
# If true, the meaning is the system is strongly balanced. 
# This function considers a bilayer system. 
function is_hb(x::Array{Float64,3}, n)
    xs = sign.(Symmetric3d(x, n));
    for l=1:size(x,3)
        v = view(xs, :,:,l);
        if sum(v*v.*v) != (n-1)*n*(n-2)
            return false
        end
    end
    return true;
end
export is_hb

# Case with only one layer. 
# checks whether multiplication of weights in all triads is positive.
# It does not have to be close to +-1. Then returns true.
# If true, the meaning is the system is strongly balanced. 
function is_hb(x::Array{Float64,2}, n)
    xs = sign.(Symmetric(x))
    if sum(xs*xs.*xs) != (n-1)*n*(n-2)
      return false
    end
    return true;
end
export is_hb

# Case with only one layer
# checks whether paradise is achieved.
# Checks whether all weights in upper triangle are nonnegative.
# (They don't have to be close to 1).
function is_paradise(x::Array{Float64,2}, n)
    if all(x .>= 0.)
      return true
    else
      return false
    end
end
export is_paradise

# Case with only one layer
# checks whether hell is achieved.
# Checks whether all weights in upper triangle are nonpositive.
# (They don't have to be close to 1).
function is_hell(x::Array{Float64,2}, n)
    if all(x .<= 0.)
      return true
    else
      return false
    end
end
export is_hell

# checks whether multiplication of weights in all triads is positive.
# It does not have to be close to +-1. Then returns true.
# If true, the meaning is the system is strongly balanced. 
# This function considers a bilayer system. 
# It returns three bool values: HB in the whole network, HB in 1st layer and HB in the 2nd layer
function is_hb2(x::Array{Float64,3}, n)
    xs = sign.(Symmetric3d(x, n));
    hb = true;
    hb_l = [true, true];
    for l=1:size(x,3)
        v = view(xs, :,:,l);
        if sum(v*v.*v) != (n-1)*n*(n-2)
            hb = false;
            hb_l[l] = false;
        end
    end
    (hb, hb_l...)
end
export is_hb2

# Calculates ratio of balanced triads in a single-layered network. 
function get_balanced_ratio(x::Array{Float64,2}, n)
    xs = sign.(Symmetric(x))
    BN = (sum(xs*xs.*xs) + n*(n-1)*(n-2))/2/6 #number of balanced triads
    return BN / n/(n-1)/(n-2)*6
end
export get_balanced_ratio

# Case with only one layer.
# Returns the counts of triad types (Delta0, Delta1, Delta2, Delta3)
# Weights in upper triangle don't have to be close to 1.
function get_triad_counts(x::Array{Float64,2}, n)
    Deltas = zeros(4)

    xs = sign.(x) .< 0 # gets Boolean array with True, where a link is negative

    for i=1:1:n, j=(i+1):1:n, k=(j+1):1:n
      # triad = [xs[i,j],xs[i,k],xs[j,k]]

      Deltas[xs[i,j] + xs[i,k] + xs[j,k] + 1] += 1
    end

    return Deltas
end
export get_triad_counts

# calculates similarity parameter between the layers
# works only for 2 layers.
function get_layer_similarity(x::Array{Float64,3}, n)
    xs = sign.(triu3d(x, 1));
    s = sum(xs[:,:,1].*xs[:,:,2]) / (n*(n-1)/2)
end
export get_layer_similarity

# calculates similarity parameter between two layers (or layer and similarity matrix)
function get_similarity(x1::Array{Float64,2}, x2::Array{Float64,2}, n)
  xs1 = sign.(triu(x1, 1))
  xs2 = sign.(triu(x2, 1))

  s = sum(xs1.*xs2) / (n*(n-1)/2)
end
export get_similarity

# calculates correlation parameter between two layers, where
# xpm1 is the layer where values are +-1 and x2 is the other layer.
# The alg. is my alg.: it is the sum of absolute differences of sign(x2) and xpm1
# weighted by the abs(x2) and normalized by sum(abs(x2)).
function get_correlation(xpm1::Array{Float64,2}, x2::Array{Float64,2})
  u = (triu(xpm1, 1))
  L = (triu(x2, 1))

  s = 1-sum(abs.(sign.(L).-sign.(u)).*abs.(L)) / sum(abs.(L)) / 2 #2 is because sign(L)-u returns 0 or 2
end
export get_correlation

# generates a balanced RL network with predefined size of larger group. 
# n - total number of nodes
# larger_size - size of the larger group
# dist_to_1 - weights don't have to be exactly +-1, but are as close as this value
# Returns the weights.
function init_random_balanced_relations(n, larger_size::Int, dist_to_1 = 0.01)
  # larger_size = Int(ceil(rand(n/2:n)))

  art_attr = [ones(larger_size, 1); -ones(n-larger_size, 1)]
  xy_attr = triu(art_attr*art_attr',1)

  xy_attr = xy_attr.*(1-dist_to_1)
  return xy_attr
end
export init_random_balanced_relations

# generates a balanced RL network with random size of larger group
# n - total number of nodes
# larger_size - size of the larger group
# dist_to_1 - weights don't have to be exactly +-1, but are as close as this value
# Returns the weights.
function init_random_balanced_relations(n, dist_to_1 = 0.01)
  # larger_size = Int(ceil(rand(n/2:n)))

  art_attr = sign.(rand(n,1).-0.5)
  xy_attr = triu(art_attr*art_attr',1)

  xy_attr = xy_attr.*(1-dist_to_1)
  return xy_attr
end

# returns counts of pos and neg links that are destabilized based only al_weights and gamma. 
function get_destabilized_links_count(rl_weights, al_weights, gamma)
  n = size(rl_weights, 1)
  rl_sim = Symmetric(rl_weights)
  lay1mul = rl_sim * rl_sim ./ (n-2)
  sgn_rl = sign.(lay1mul)
  sgn_al = sign.(al_weights)
  difs = (sgn_rl .* sgn_al) .< 0
  if any(difs)
    links = lay1mul[difs] .+ al_weights[difs] .* gamma
    sgn_links = sign.(links)

    destab = (sgn_rl[difs] .* sgn_links) .< 0

    pos = sum(sgn_rl[difs][destab] .> 0)
    neg = sum(destab) - pos
  else
    pos = 0
    neg = 0
  end
  return (pos, neg)
end
export get_destabilized_links_count

# end
