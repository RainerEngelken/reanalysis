## first install julia from julialang.org/downloads/
## Then start julia and  install the modules Distributions and PyPlot by the following command:
## Pkg.add("Distributions"); Pkg.add("PyPlot")
## The following parameter generate a quick low resolution version of figure 1d.
## For higher precision, use a dV = 0.1; tWarmup = 5000; tSim = 5000

using Distributions#, PyPlot
include("spikingNetworkPerturbation.jl")

JRange = (4)/10 # J: coupling strength in mV
dV = 1		# dV: integration step size for MFT
dt = 1		# dt:time step in rate simulation
tSim = 100	# tSim: simulation time in ms
tWarmup = 100	# tWarmup: warmup time in ms
pertsize = 100	# pertsize: along least stable direction
seed = 18	# seed: topology and initial conditions

N = 10000
C = 1000
nPert = 1000
J = 0.6
tauRef = 0.5/1000
tauDelay = 0.55/1000
seedIC = 1
seedTopo = 18
pertsize = 100
for J in JRange
    println("J = ", J)
    spikingNetworkPerturbation(N,C,nPert,J,tauRef,tauDelay,seedIC,seedTopo,pertsize)
end

