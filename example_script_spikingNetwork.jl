## first install julia from julialang.org/downloads/
## Then start julia and  install the modules Distributions and PyPlot by the following command:
## Pkg.add("Distributions"); Pkg.add("PyPlot")
## The following parameter generate a quick low resolution version of figure 1a
## For higher precision, use tSim = 200
using Distributions, PyPlot
include("spikingNetwork.jl")

JRange = 1:10 # J: coupling strength in mV/10
τ_delayRange = [0, 0.00055]
N = 10000
C = 1000
tSim = 1	# tSim: simulation time in s
J = 0.6
τ_ref = 0.0005
τ_delay = 0.00055
seed = 1
rateSim = Float64[]
for τ_delay in τ_delayRange
    for J in JRange
        println("J = ", J/10)
        push!(rateSim,spikingNetwork(N,C,tSim,J/10,τ_ref,τ_delay,seed))
    end
end
rateSim = reshape(rateSim,length(JRange),length(τ_delayRange))
plot(JRange/10,rateSim)
xlabel("coupling J (mV)")
ylabel("rate (1/s)")
legend(["delay = 0 ms","delay = 0.55 ms"])
