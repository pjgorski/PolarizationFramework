using DrWatson

quickactivate(@__DIR__)

using PolarizationFramework
using DifferentialEquations

g = 5
attr = BinaryAttributes(g)

calc_heider_attr(45, 0.0, 1000.0, Heider7!, AutoTsit5(Rodas5(autodiff = false)), true, attr)
