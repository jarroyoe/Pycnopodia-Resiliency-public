using Plots

include("pycno_system.jl")
include("parameters.jl")
include("quasipotentials_in_julia.jl")

paramToModify = [10,5]
parameterRange = collect(Iterators.product(0.1:0.1:10,0:0.01:1))
qps = Array{Float64}(undef,size(parameterRange))

for i in CartesianIndices(parameterRange)
    p[paramToModify] .= parameterRange[i]
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
    #alpha0 += [rand() for i in 1:length(alpha0)]
    function g(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21)
        return compute(action.dim,action.adj_n,action.Bs,action.cs,action.quadrature_w,[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21],p)
    end

    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_NEWUOA)
    set_optimizer_attribute(model, "xtol_rel", 1e-12)
    @variable(model, x[i = 1:length(alpha0)])
    register(model, :g, length(alpha0), g; autodiff = true)
    @NLobjective(model,Min,g(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21]))
    for i in 1:length(alpha0)
        set_start_value(x[i],alpha0[i])
    end
    JuMP.optimize!(model)
    qp = objective_value(model)
    println(qp)
    qps[i] = qp
end