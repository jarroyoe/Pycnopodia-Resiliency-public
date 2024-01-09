using DelimitedFiles
using Distributions
using Random
using JLD2
#Random.seed!(11002458920032)

include("pycno_system.jl")
include("parameters.jl")
include("quasipotentials_in_julia.jl")

# number of samples to take
m = 1200

# Create matrix of samples for parameters of system
p_samples = hcat(
    rand(truncated(Normal(2.5,1.5); lower = 0.001),m), # r
    rand(truncated(Normal(10000,5000); lower = 0.001),m), # K
    rand(truncated(Normal(2,1); lower = 0.001),m) .* rand(Bernoulli(0.9),m), # κ_D
    rand(LogNormal(-3.7,1),m), # α_A
    rand(LogNormal(-0.2,1),m) .* rand(Bernoulli(0.9),m), # κ_S
    rand(LogUniform(exp(-4),5),m), # δ_A
    rand(Uniform(0,1),m), # ϵ_D
    rand(LogUniform(1e-5,1e-1),m), # ϵ_U
    rand(LogNormal(-2.7,1),m), # α_D
    rand(Uniform(0,10),m), # α_U
    rand(LogUniform(exp(-2),exp(3)),m), # γ_U
    rand(LogNormal(-8,1),m) .* rand(Bernoulli(0.9),m), # κ_A
    rand(LogNormal(-8,4),m), # δ_U
    rand(LogUniform(exp(-2),exp(2)),m) .* rand(Bernoulli(0.9),m), # β
    rand(LogUniform(1e-5,1e-1),m), # ϵ_S
    rand(LogUniform(exp(-9),exp(-2)),m), # δ_S
    rand(LogUniform(.001,2),m) # δ_D
    )


# array to store results for each simulation
action_data = Dict([("dim",Int64[]),("adj_n",Int64[]),("Bs",Array{Matrix{Float64}}[]),
    ("cs",Array{Matrix{Float64}}[]),("quadrature_w",Array{Float64}[]),("alpha0",Array{Float64}[]),("p",Array{Float64}[])])

#Threads.@threads 
for i in 1:m
    for j in 1:(length(p)-1)
        p[j] = p_samples[i,j]
    end
    #Find fixed points
    f(x) = pycno_roots_with_levers(x,p)[3]
    fps = try 
            find_zeros(f,0,1.1*p[2])
        catch e
            Inf
            continue
        end
    if length(fps)<2
        continue
    end
    edges = [minimum(fps),maximum(fps)]

    fp = [[pycno_roots(edges[i],p)[1],
        pycno_roots(edges[i],p)[2],
        edges[i]] for i in 1:2]
    #Compute the quasipotential
    action = py"pyritz.interpolation.Action"(lagrangian!, n, nq, fp[1], fp[2])
    alpha0 = py"pyritz.interpolation.utils.linear_path"(fp[1], fp[2], n)
    alpha0 += [rand() for i in 1:length(alpha0)]
    
    push!(action_data["dim"],action.dim)
    push!(action_data["adj_n"],action.adj_n)
    push!(action_data["cs"],action.cs)
    push!(action_data["Bs"],action.Bs)
    push!(action_data["quadrature_w"],action.quadrature_w)
    push!(action_data["alpha0"],alpha0)
    push!(action_data["p"],p[1:(length(p)-1)])
    println(i)
    if mod(i,100)==0
	JLD2.save("action_parameters2.jld2",action_data)
    end
end

JLD2.save("action_parameters2.jld2",action_data)

