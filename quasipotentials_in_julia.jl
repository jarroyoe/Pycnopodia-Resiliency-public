using PyCall
using LinearAlgebra
using OMEinsum
using JuMP, NLopt
include("parameters.jl")
include("pycno_system.jl")

#Load action parameters
n = 8
nq = 80

#Load Python PyRitz package
#=    py"""
    import os, sys

    #PyRitz directory goes here
    pyritz_dir = "~/PyRitz"
    nlopt_lib_dir = "%s/nlopt/lib" % pyritz_dir
    nlopt_py_dir = "%s/nlopt/nlopt_py" % pyritz_dir

    if "LD_LIBRARY_PATH" in os.environ:
        paths = os.environ["LD_LIBRARY_PATH"].split(":")
        paths.append(nlopt_lib_dir)
        os.environ["LD_LIBRARY_PATH"] = ":".join(paths)
    else:
        os.environ["LD_LIBRARY_PATH"] = ":%s" % nlopt_lib_dir
    
    sys.path.insert(0, nlopt_py_dir)
    sys.path.insert(0, pyritz_dir)

    # Import statements

    import pyritz

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import nlopt
    """
=#
#Define Lagrangian
function lagrangian(path,p)
    xs, vs = path
    tsas = [pycno_with_levers(xs[i,:],p,0) for i in 1:size(xs,1)]
    sas = [tsas[i][j] for i in 1:length(tsas), j in 1:3]
    v_norms = [norm(vs[i,:]) for i in 1:size(vs,1)]
    sa_norms = [norm(sas[i,:]) for i in 1:size(sas,1)]
    vs_dot_sas = ein"ij,ij->j"(vs',sas')
    return v_norms .* sa_norms - vs_dot_sas
end
function lagrangian!(ls, dxls, dvls, path, us, args)
    xs, vs = path
    tsas = [pycno_with_levers_JuMP(xs[:,i],p,0) for i in 1:size(xs,2)]
    sas = [tsas[i][j] for i in 1:length(tsas), j in 1:3]
    v_norms = [norm(vs[:,i]) for i in 1:size(vs,2)]
    sa_norms = [norm(sas[i,:]) for i in 1:size(sas,1)]
    vs_dot_sas = ein"ij,ij->j"(vs,sas')
    ls[:] = v_norms .* sa_norms - vs_dot_sas
end

function compute(dim,adj_n,Bs,cs,quadrature_w,alpha,p)
    alpha_reshaped = reshape(alpha,(dim,adj_n))'   

    newPath = [Bs[i]*alpha_reshaped+cs[i] for i in 1:2]

    L = lagrangian(newPath,p)
    #action.path = newPath
    return log(sum(quadrature_w .* L))
end

function quasipotential(dim,adj_n,Bs,cs,quadrature_w,alpha0,p)
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
    return model
end
