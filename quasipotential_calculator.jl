using DelimitedFiles
using Distributions
using JLD2

include("pycno_system.jl")
include("parameters.jl")
include("quasipotentials_in_julia.jl")
action_data = JLD2.load("action_parameters.jld2")
m = min(length(action_data["dim"]),2000)

quasi_data = Array{Float64}(undef,m,length(p))

@Threads.threads for i in 1:m
    dim = action_data["dim"][i]
    adj_n = action_data["adj_n"][i]
    Bs = action_data["Bs"][i]
    cs = action_data["cs"][i]
    quadrature_w = action_data["quadrature_w"][i]
    alpha0 = action_data["alpha0"][i]
    ps = action_data["p"][i]
    function g(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21)
        return compute(dim,adj_n,Bs,cs,quadrature_w,[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21],p)
    end

    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_NEWUOA)
    set_optimizer_attribute(model, "xtol_rel", 1e-8)
    @variable(model, x[i = 1:length(alpha0)])
    register(model, :g, length(alpha0), g; autodiff = true)
    @NLobjective(model,Min,g(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21]))
    for i in 1:length(alpha0)
        set_start_value(x[i],alpha0[i])
    end
    JuMP.optimize!(model)
    qp = objective_value(model)
    println(qp)
    global quasi_data[i,:] = vcat(ps,qp)
    writedlm("quasi_data_"*string(Threads.threadid())*".csv",quasi_data)
end

writedlm("quasi_data.csv",quasi_data)