using DynamicalSystems
using DifferentialEquations
using Roots

# Model Definition
function pycno(N,p,t)
    
    ## Parameters
    r = p[1]    # growth rate of kelp at low density
    K = p[2]    # carrying capacity of kelp
    κ_D = p[3]  # effect of drift presence on urchin grazing rate
    α_A = p[4]  # urchin grazing rate at low kelp densities
    κ_S = p[5]  # effect of fear of predation on urchin grazing rate
    δ_A = p[6]  # conversion rate of live kelp biomass to drift kelp
    ϵ_D = p[7]  # proportion of drift kelp retained in the system
    ϵ_U = p[8]  # conversion efficiency from grams of grazed kelp to urchins
    α_D = p[9] # drift kelp consumption rate by urchins
    α_U = p[10] # seastar predation rate at low urchin counts
    γ_U = p[11] # seastar predation saturation constant
    κ_A = p[12] # effect of urchins starvation on seastar predation rate
    δ_U = p[13] # natural death rate of urchins
    β = p[14]   # impact of kelp density on nutritional value of urchins
    ϵ_S = p[15] # conversion proportion from urchins predated to seastars
    δ_S = p[16] # natural death rate of seastars
    δ_D = p[17] # natural (death rate) of drift kelp

    ## Variables
    S = max(N[1],0) # seastar population
    U = max(N[2],0) # urchin population
    A = max(N[3],0) # kelp population

    ## Drift Kelp Equation
    Dstar = A != 0 ? dStar(A,U,ϵ_D,δ_A,κ_D,α_D,δ_D) : 0

    ## Differential Equations
    dS = ((ϵ_S * (1 + β*(A+Dstar)) * α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_S*S
    dU = (1/(1 + κ_D*Dstar)) * (ϵ_U * α_A*A*U) * (1 / (1 + κ_S*S)) + (1 - (1/(1 + κ_D*Dstar))) * 
    (ϵ_U*α_D*Dstar*U) - ((α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_U*U
    dA = r*A*(1 - (A/K)) - (1/(1 + κ_D*Dstar)) * (α_A*A*U) * (1 / (1 + κ_S*S)) - δ_A*A

    return SVector{3}(dS,dU,dA)
end

# Model to find fixed points
function pycno_roots(A,p)
    
    ## Parameters
    r = p[1]    # growth rate of kelp at low density
    K = p[2]    # carrying capacity of kelp
    κ_D = p[3]  # effect of drift presence on urchin grazing rate
    α_A = p[4]  # urchin grazing rate at low kelp densities
    κ_S = p[5]  # effect of fear of predation on urchin grazing rate
    δ_A = p[6]  # conversion rate of live kelp biomass to drift kelp
    ϵ_D = p[7]  # proportion of drift kelp retained in the system
    ϵ_U = p[8]  # conversion efficiency from grams of grazed kelp to urchins
    α_D = p[9] # drift kelp consumption rate by urchins
    α_U = p[10] # seastar predation rate at low urchin counts
    γ_U = p[11] # seastar predation saturation constant
    κ_A = p[12] # effect of urchins starvation on seastar predation rate
    δ_U = p[13] # natural death rate of urchins
    β = p[14]   # impact of kelp density on nutritional value of urchins
    ϵ_S = p[15] # conversion proportion from urchins predated to seastars
    δ_S = p[16] # natural death rate of seastars
    δ_D = p[17] # natural (death rate) of drift kelp

    ## Equations for A
    Ustar(D) = D != 0 ? (ϵ_D*δ_A*A-δ_D*D)/(α_D*(1-1/(1+κ_D*D))*D) : 0
    dS(D) = ϵ_S*(1+β*(A+D))*α_U*Ustar(D)/(1+γ_U*Ustar(D))*1/(1+κ_A*(A+D))-δ_S
    Dstar = A != 0 ? max(find_zero(dS,A,Order2()),0) : 0
    U = Ustar(Dstar)
    dU(S) = (1/(1 + κ_D*Dstar)) * (ϵ_U * α_A*A*U) * (1 / (1 + κ_S*S)) + (1 - (1/(1 + κ_D*Dstar))) * 
    (ϵ_U*α_D*Dstar*U) - ((α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_U*U
    S = A != 0 ? max(find_zero(dU,1.,Order2()),0) : 0
    dA = r*A*(1 - (A/K)) - (1/(1 + κ_D*Dstar)) * (α_A*A*U) * (1 / (1 + κ_S*S)) - δ_A*A

    return SVector{3}(S,U,dA)
end

# Function for Building Dynamical System
function buildSystem(N,p)
    return ContinuousDynamicalSystem(pycno,N,p)
end

function dStar(A,U,ϵ_D,δ_A,κ_D,α_D,δ_D)
    a = δ_D*κ_D.+κ_D*α_D*U
    b = ϵ_D*δ_A*κ_D*A.-δ_D
    c = ϵ_D*δ_A*A

    return (b.+sqrt.(b .^2 .+4*a.*c))./(2*a)
end

function pycno_with_levers(N,p,t)
    #Levers for nonconsumptive feedbacks
    #1: κ_D
    #2: κ_S
    #3: κ_A
    #4: β
    levers = p[end]
    
    ## Parameters
    r = p[1]    # growth rate of kelp at low density
    K = p[2]    # carrying capacity of kelp
    κ_D = levers[1]*p[3]  # effect of drift presence on urchin grazing rate
    α_A = p[4]  # urchin grazing rate at low kelp densities
    κ_S = levers[2]*p[5]  # effect of fear of predation on urchin grazing rate
    δ_A = p[6]  # conversion rate of live kelp biomass to drift kelp
    ϵ_D = p[7]  # proportion of drift kelp retained in the system
    ϵ_U = p[8]  # conversion efficiency from grams of grazed kelp to urchins
    α_D = p[9] # drift kelp consumption rate by urchins
    α_U = p[10] # seastar predation rate at low urchin counts
    γ_U = p[11] # seastar predation saturation constant
    κ_A = levers[3]*p[12] # effect of urchins starvation on seastar predation rate
    δ_U = p[13] # natural death rate of urchins
    β = levers[4]*p[14]   # impact of kelp density on nutritional value of urchins
    ϵ_S = p[15] # conversion proportion from urchins predated to seastars
    δ_S = p[16] # natural death rate of seastars
    δ_D = p[17] # natural (death rate) of drift kelp

    ## Variables
    S = max(N[1],0) # seastar population
    U = max(N[2],0) # urchin population
    A = max(N[3],0) # kelp population

    ## Drift Kelp Equation
    Dstar = ((κ_D != 0) & (A != 0)) ? dStar(A,U,ϵ_D,δ_A,κ_D,α_D,δ_D) : 0

    ## Differential Equations
    dS = ((ϵ_S * (1 + β*(A+Dstar)) * α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_S*S
    dU = (1 ./(1 .+ κ_D*Dstar)) .* (ϵ_U * α_A*A*U) .* (1 ./ (1 .+ κ_S*S)) .+ (1 .- (1 ./(1 .+ κ_D*Dstar))) .* 
    (ϵ_U*α_D*Dstar*U) .- ((α_U*U.*S) ./ (1 .+ γ_U*U)) .* (1 ./ (1 .+ κ_A*(A+Dstar))) .- δ_U*U
    dA = r*A.*(1 .- (A/K)) .- (1 ./(1 .+ κ_D*Dstar)) .* (α_A*A.*U) .* (1 ./ (1 .+ κ_S*S)) .- δ_A*A

    return SVector{3}(dS,dU,dA)
end

function pycno_roots_with_levers(A,p)
    #Levers for nonconsumptive feedbacks
    #1: κ_D
    #2: κ_S
    #3: κ_A
    #4: β
    levers = p[end]
    
    ## Parameters
    r = p[1]    # growth rate of kelp at low density
    K = p[2]    # carrying capacity of kelp
    κ_D = levers[1]*p[3]  # effect of drift presence on urchin grazing rate
    α_A = p[4]  # urchin grazing rate at low kelp densities
    κ_S = levers[2]*p[5]  # effect of fear of predation on urchin grazing rate
    δ_A = p[6]  # conversion rate of live kelp biomass to drift kelp
    ϵ_D = p[7]  # proportion of drift kelp retained in the system
    ϵ_U = p[8]  # conversion efficiency from grams of grazed kelp to urchins
    α_D = p[9] # drift kelp consumption rate by urchins
    α_U = p[10] # seastar predation rate at low urchin counts
    γ_U = p[11] # seastar predation saturation constant
    κ_A = levers[3]*p[12] # effect of urchins starvation on seastar predation rate
    δ_U = p[13] # natural death rate of urchins
    β = levers[4]*p[14]   # impact of kelp density on nutritional value of urchins
    ϵ_S = p[15] # conversion proportion from urchins predated to seastars
    δ_S = p[16] # natural death rate of seastars
    δ_D = p[17] # natural (death rate) of drift kelp

    ## Equations for A
    Ustar(D) = D != 0 ? max((ϵ_D*δ_A*A-δ_D*D)/(α_D*(1-1/(1+κ_D*D))*D),0) : 0
    dS(D) = ϵ_S*(1+β*(A+D))*α_U*Ustar(D)/(1+γ_U*Ustar(D))*1/(1+κ_A*(A+D))-δ_S
    Dstar = (A != 0) ? try max(find_zero(dS,A,Order2()),0) catch e 0 end : 0
    U = Ustar(Dstar)
    dU(S) = (1/(1 + κ_D*Dstar)) * (ϵ_U * α_A*A*U) * (1 / (1 + κ_S*S)) + (1 - (1/(1 + κ_D*Dstar))) *
    (ϵ_U*α_D*Dstar*U) - ((α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_U*U
    S = A != 0 ? max(find_zero(dU,1.,Order2()),0) : 0
    dA = r*A*(1 - (A/K)) - (1/(1 + κ_D*Dstar)) * (α_A*A*U) * (1 / (1 + κ_S*S)) - δ_A*A

    return SVector{3}(S,U,dA)
end

function pycno_with_levers_JuMP(N,p,t)
    #Levers for nonconsumptive feedbacks
    #1: κ_D
    #2: κ_S
    #3: κ_A
    #4: β
    levers = p[end]
    
    ## Parameters
    r = p[1]    # growth rate of kelp at low density
    K = p[2]    # carrying capacity of kelp
    κ_D = levers[1]*p[3]  # effect of drift presence on urchin grazing rate
    α_A = p[4]  # urchin grazing rate at low kelp densities
    κ_S = levers[2]*p[5]  # effect of fear of predation on urchin grazing rate
    δ_A = p[6]  # conversion rate of live kelp biomass to drift kelp
    ϵ_D = p[7]  # proportion of drift kelp retained in the system
    ϵ_U = p[8]  # conversion efficiency from grams of grazed kelp to urchins
    α_D = p[9] # drift kelp consumption rate by urchins
    α_U = p[10] # seastar predation rate at low urchin counts
    γ_U = p[11] # seastar predation saturation constant
    κ_A = levers[3]*p[12] # effect of urchins starvation on seastar predation rate
    δ_U = p[13] # natural death rate of urchins
    β = levers[4]*p[14]   # impact of kelp density on nutritional value of urchins
    ϵ_S = p[15] # conversion proportion from urchins predated to seastars
    δ_S = p[16] # natural death rate of seastars
    δ_D = p[17] # natural (death rate) of drift kelp

    ## Variables
    S = N[1] # seastar population
    U = N[2] # urchin population
    A = N[3] # kelp population

    ## Drift Kelp Equation
    Dstar = ((κ_D != 0) & (A != 0)) ? dStar(A,U,ϵ_D,δ_A,κ_D,α_D,δ_D) : 0

    ## Differential Equations
    dS = ((ϵ_S * (1 + β*(A+Dstar)) * α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_S*S
    dU = (1/(1 + κ_D*Dstar)) * (ϵ_U * α_A*A*U) * (1 / (1 + κ_S*S)) + (1 - (1/(1 + κ_D*Dstar))) * 
    (ϵ_U*α_D*Dstar*U) - ((α_U*U*S) / (1 + γ_U*U)) * (1 / (1 + κ_A*(A+Dstar))) - δ_U*U
    dA = r*A*(1 - (A/K)) - (1/(1 + κ_D*Dstar)) * (α_A*A*U) * (1 / (1 + κ_S*S)) - δ_A*A

    return SVector{3}(dS,dU,dA)
end