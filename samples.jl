using DelimitedFiles
using Distributions
using Random
Random.seed!(11002093)

include("pycno_system.jl")
include("parameters.jl")

# number of samples to take
n = 2000

# Create matrix of samples for initial conditions
N_samples = hcat(
    rand(Uniform(0,20), n), # S
    rand(Uniform(0,20), n), # U
    rand(Uniform(0,4000), n) # A
)

# Create matrix of samples for parameters of system
p_samples = hcat(
    rand(truncated(Normal(2.5,1.5); lower = 0.001),n), # r
    rand(truncated(Normal(10000,5000); lower = 0.001),n), # K
    rand(truncated(Normal(2,1); lower = 0.001),n) .* rand(Bernoulli(0.9),n), # κ_D
    rand(LogNormal(-3.7,1),n), # α_A
    rand(LogNormal(-0.2,1),n) .* rand(Bernoulli(0.9),n), # κ_S
    rand(LogUniform(exp(-4),5),n), # δ_A
    rand(Uniform(0,1),n), # ϵ_D
    rand(LogUniform(1e-5,1e-1),n), # ϵ_U
    rand(LogNormal(-2.7,1),n), # α_D
    rand(Uniform(0,10),n), # α_U
    rand(LogUniform(exp(-2),exp(3)),n), # γ_U
    rand(LogNormal(-8,1),n) .* rand(Bernoulli(0.9),n), # κ_A
    rand(LogNormal(-8,4),n), # δ_U
    rand(LogUniform(exp(-2),exp(2)),n) .* rand(Bernoulli(0.9),n), # β
    rand(LogUniform(1e-5,1e-1),n), # ϵ_S
    rand(LogUniform(exp(-9),exp(-2)),n), # δ_S
    rand(LogUniform(.001,2),n) # δ_D
)


# array to store results for each simulation
traj_data = Array{Float64}(undef,n,3)
# number of weeks to simulate
t = 100

Threads.@threads for i in 1:n
    for j in 1:(length(p)-1)
        p[j] = p_samples[i,j]
    end
    N = N_samples[i,:]
    pycno_sys_i = ContinuousDynamicalSystem(pycno_with_levers,N,p)
    traj = trajectory(pycno_sys_i,t)[end,:]
    global traj_data[i,:] = [traj[1] traj[2] traj[3]]
end

sim_data = hcat(N_samples, p_samples, traj_data)
writedlm("sim_data.csv",sim_data)
