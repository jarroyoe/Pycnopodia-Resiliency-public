using DelimitedFiles
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

kappa_A_range = 0:1e-5:1e-3
kappa_D_range = 0:0.1:10
#base_kappa_A = p[12]
#base_kappa_D = p[3]

# array to store results
sim_data_likelihood = Array{Float64}(undef, length(kappa_A_range),2)
sim_data_rate = Array{Float64}(undef, length(kappa_A_range),2)

# loop to run simulations
for k in 1:2
    @Threads.threads for i in 1:length(kappa_A_range)
   	p_copy = copy(p)    
	if k == 1
		p_copy[12] = kappa_A_range[i]
	else
	    	p_copy[3] = kappa_D_range[i]
	end
    	pycno_sys_ij = ContinuousDynamicalSystem(pycno_with_levers,[1,1,1],p_copy)
    	f(x) = pycno_roots_with_levers(x,p_copy)[3]
    #fps = find_zeros(f,0,1.1*p[2])
    	threshold = try 
        	fps = find_zeros(f,0,1.1*p_copy[2])
        	maximum(fps)
    	catch e
        	Inf
        	continue
    	end
    	forest = 0
    	if threshold != Inf
    		for j in 1:length(N_grid)
        		N = N_grid[j,:][1]
        		traj = trajectory(pycno_sys_ij, t, N)[end,:] # columns: [1] is S, [2] is U, [3] is A
        		forest += traj[3] ≥ threshold*0.95
    		end
    	end
    	kelpForest = [pycno_roots_with_levers(threshold,p_copy)[1],pycno_roots_with_levers(threshold,p_copy)[2],threshold]
    	prop_forest = forest/size(N_grid,1)
    	eig = maximum(real.(-eigvals(jacobian(pycno_sys_ij,kelpForest,p_copy))))
    	global sim_data_likelihood[i,k] = prop_forest
    	global sim_data_rate[i,k] = eig
    	writedlm("sim_data_"*string(Threads.threadid())*".csv",sim_data_likelihood)
    	writedlm("sim_data_"*string(Threads.threadid())*".csv",sim_data_rate)
    end
end

# save results as csv
writedlm("likelihood_data_otherParams.csv",sim_data_likelihood)
writedlm("rate_data_otherParams.csv",sim_data_rate)
