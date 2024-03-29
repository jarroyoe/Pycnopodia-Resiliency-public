r = 2.5   # growth rate of kelp at low density
K = 10000    # carrying capacity of kelp
κ_D = 1.95  # effect of drift presence on urchin grazing rate
α_A = 0.025  # urchin grazing rate at low kelp densities
κ_S = 0.81  # effect of fear of predation on urchin grazing rate
δ_A = 1.8  # conversion rate of live kelp biomass to drift kelp
ϵ_D = 0.7  # proportion of drift kelp retained in the system
ϵ_U = 0.1  # conversion efficiency from grams of grazed kelp to urchins
α_D = 0.062 # drift kelp consumption rate by urchins
α_U = 4.77 # seastar predation rate at low urchin counts
γ_U = 3.42 # seastar predation saturation constant
κ_A = 0.00013 # effect of urchins starvation on seastar predation rate
δ_U = 0.0004 # natural death rate of urchins
β = 0.1   # impact of kelp density on nutritional value of urchins
ϵ_S = 0.1 # conversion proportion from urchins predated to seastars
δ_S = 0.0001 # natural death rate of seastars
δ_D = 0.3 # natural (death rate) of drift kelp

levers = [1,1,1,1]
p = [r,K,κ_D,α_A,κ_S,δ_A,ϵ_D,ϵ_U,α_D,α_U,γ_U,κ_A,δ_U,β,ϵ_S,δ_S,δ_D,levers]