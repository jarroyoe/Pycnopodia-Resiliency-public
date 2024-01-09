#using HomotopyContinuation
using JLD2
include("pycno_system.jl")
include("parameters.jl")

# Define Dynamical System
## kelp is at 5% of carrying capacity
## r reduced to 10% of original value: 584.26 -> 58.426
## unkown parameters:
### ϵ_U = 0.10
### γ_U = 1.0
### β = 1.0
### ϵ_S = 0.1
### δ_S = 0.0001
### δ_D = 0.3
pycno_sys1 = buildSystem([1.0,80.0,166.0], p)


# Find Attractors
## define range
xg = yg = range(0, 100; length = 20)
zg = range(0,4000; length = 20)

## Find basins of attraction
mapper = AttractorsViaRecurrences(pycno_sys1,(xg,yg,zg);diffeq = (alg=RadauIIA5(), reltol=1e-30, abstol=1e-30))
basins, attractors = basins_of_attraction(mapper)
array_attractors = [Matrix(attractors[i]) for i in 1:length(attractors)]

## create plot
####plot(basins) # basins is not a plottable object. need to transform.

## saving results
jldsave("attractors.jld2";array_attractors,basins)
