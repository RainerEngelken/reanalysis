## first install julia from julialang.org/downloads/
## Then start julia and  install the modules Distributions and PyPlot by the following command:
## Pkg.add("Distributions"); Pkg.add("PyPlot")
## The following parameter generate a quick low resolution version of figure 1d.
## For higher precision, use a dV = 0.1; tWarmup = 5000; tSim = 5000

using Distributions, PyPlot
include("ratePerturbation.jl")

JRange = (1:10) # J: coupling strength in mV
dV = 1		# dV: integration step size for MFT
dt = 1		# dt:time step in rate simulation
tSim = 200	# tSim: simulation time in ms
tWarmup = 5000	# tWarmup: warmup time in ms
pertsize = 100	# pertsize: along least stable direction
seed = 18	# seed: topology and initial conditions

for J in JRange
    println("J = ", J)
    perturbeRatenetworkLIF(J/10,dV,dt,tSim,tWarmup,pertsize,seed)
end
plot()
