using DelimitedFiles, ForwardDiff
using Distributions
using LinearAlgebra
using Random
#Random.seed!(11014036)

include("pycno_system.jl")
include("parameters.jl")

# create grid for initial conditions
sg = range(0,2; length = 5)
ug = range(0,100; length = 5)
ag = range(0,4000; length = 5)
N_grid = reshape(collect(Iterators.product(sg,ug,ag)),125,1)
N_grid = [collect(N_grid[i]) for i in 1:125]

# number of parameter samples to take
n = 4000

#Create matrix of samples for parameters of system
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

# number of weeks to simulate
t = 1000

# array to store results
sim_data_stab = Array{Float64}(undef, n, (size(p_samples,2)+2))

# loop to run simulations
Threads.@threads for i in 1:n
    p_copy = copy(p)
    p_copy[1:(end-1)]=p_samples[i,:]
    pycno_sys_ij = ContinuousDynamicalSystem(pycno_with_levers,[1.,1.,1.],p_copy)
    f(x) = pycno_roots_with_levers(x,p_copy)[3]
#    fps = find_zeros(f,0,1.1*p[2])
    threshold = try 
        fps = find_zeros(f,0,1.1*p_copy[2])
        maximum(fps)
    catch e
        Inf
        continue
    end
    forest = 0
   if threshold != Inf
    	for j in 1:(size(N_grid,1))
        	N = N_grid[j,:][1]
        	traj, T = trajectory(pycno_sys_ij, t, N) # columns: [1] is S, [2] is U, [3] is A
        	forest += traj[:,3][end] ≥ threshold*0.95
    	end
    end
    kelpForest = [pycno_roots_with_levers(threshold,p_copy)[1],pycno_roots_with_levers(threshold,p_copy)[2],threshold]
    prop_forest = forest/size(N_grid,1)
    eig = maximum(real.(-eigvals(ForwardDiff.jacobian(x->pycno_with_levers(x,p_copy,0),kelpForest))))
    global sim_data_stab[i,:] = vcat(p_samples[i,:],prop_forest,eig)
    writedlm("sim_data_"*string(Threads.threadid())*".csv",sim_data_stab)
end

# save results as csv
writedlm("sim_data_stab.csv",sim_data_stab)
