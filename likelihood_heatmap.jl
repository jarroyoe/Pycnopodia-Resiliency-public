using DelimitedFiles, ForwardDiff
using Distributions
using FLoops
using LinearAlgebra
#Random.seed!(11014036)

include("pycno_system.jl")
include("parameters.jl")

# create grid for initial conditions
sg = range(0,2; length = 5)
ug = range(0,100; length = 5)
ag = range(0,4000; length = 5)
N_grid = reshape(collect(Iterators.product(sg,ug,ag)),125,1)
N_grid = [collect(N_grid[i]) for i in 1:125]


# number of weeks to simulate
t = 1000

#epsilon_U_range = 10 .^(-5:0.05:0)
#alpha_D_range = 0.01:0.1:10
dA_range = 0:0.02:5
r_range = 0.1:0.1:10
#K_range = 1000:100:10000
a = length(r_range)
b = length(dA_range)

# array to store results
sim_data_likelihood = zeros(a)

# loop to run simulations
@floop for i in 1:a
    p_copy = copy(p)
    p_copy[1] = r_range[i]
    #p_copy[6] = dA_range[i]
    pycno_sys_ij = ContinuousDynamicalSystem(pycno_with_levers,[1.,1.,1.],p_copy)
    f(x) = pycno_roots_with_levers(x,p_copy)[3]
    #fps = find_zeros(f,0,1.1*p[2])
    threshold = try 
        fps = find_zeros(f,0,1.1*p_copy[2])
        maximum(fps)
    catch e
        global sim_data_likelihood[i[1],i[2]] = 0
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
    #eig = maximum(real.(-eigvals(ForwardDiff.jacobian(x->pycno_with_levers(x,p_copy,0),kelpForest))))
    global sim_data_likelihood[i] = prop_forest
    writedlm("sim_data_"*string(Threads.threadid())*".csv",sim_data_likelihood)
    #writedlm("sim_data_"*string(Threads.threadid())*".csv",sim_data_rate)
end

# save results as csv
writedlm("r_likelihood_data.csv",sim_data_likelihood)
#writedlm("rate_data.csv",sim_data_rate)
