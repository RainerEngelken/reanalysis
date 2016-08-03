##################################################################
###     	rate network of Ostojic Nature Neuroscience 2014	##
### 	Copyright Rainer Engelken, 2016 rainer@nld.ds.mpg.de	##
##################################################################
function perturbeRatenetworkLIF(J,dV,dt,tSim,tWarmup,pertsize,seed)
############ Input parameters:###################
# J: coupling strength (mV)			#
# dV: integration step size for MFT		#
# dt:time step in rate simulation (ms)		#
# tSim: simulation time	(ms)			#
# tWarmup: warmup time (ms)			#
# pertsize: along least stable direction	#
# seed: topology and initial conditions		#
#################################################
# example: using Distributions;include("rc06ratePerturbation.jl");perturbeRatenetworkLIF(0.3,0.1,1,300,300,1,2)

##################
### PARAMETERS ###
##################
nStepSim = round(Int,tSim/dt) # number of simulation steps after perturbation
nStepWarmup = round(Int,tWarmup/dt) # number of warmup steps before perturbation
sPrint = 10 	# frequency of output
f = 0.8 	# fraction of excitation in the network
g = 5 		# relative inhibitory strength
τ_m = 20 	# membrane time constant (ms)
Vth = 20 	# threshold voltage (mV)
Vre = 10 	# reset voltage (mV)
μ0 = 24 	# constant external input (mV)
C = 1000 	# number of synapses per neuron
N = 2000 	# number of neurons in the network
τRef = 0.5 	# refractory time (ms)
Vl = -50	# the lower bound for the voltage range
tPerturbe = 2	# perturbation time (ms)
ε = 1e-5 	# slope
#################
### FUNCTIONS ###
#################
include("makeConnectivity.jl")
include("getMaxRealEigenvector.jl")
include("lifMFT_Richardson.jl")

function runSimulation(nStep,sPrint,dt,ν,μ0Perturbed)
    t0 = time_ns()
    νAll = zeros(N,nStep)
    νAll[:,1] = νInit
    F = ones(N)
    μ0Unperturbed = μ0*ones(N)
    for s = 2 : nStep 
        (s-1)*dt <= tPerturbe ? μ0Vec = μ0Perturbed : μ0Vec = μ0Unperturbed # perturbe μ0 during tPerturbe
        μ = μ0Vec + τ_m * JMatrix * ν 	# equation 18 in Ostojic 2014
        σ2 = τ_m * JMatrix2 * ν 	# equation 19 in Ostojic 2014
        σ = sqrt(σ2/2)
        for n = 1 : N
            F[n] = LIF0dVref(τ_m,σ[n],μ[n],Vth,Vre,τRef,dV,Vl) # equ. 1 via Richardson method
        end
        ν = ν + dt/τ_m * (-ν + F) 	# equation 17 in Ostojic 2014
        νAll[:,s] = ν
        if mod(s,sPrint)==0
            t=(time_ns()-t0)/1e9;
            println(round(100s/nStep),"% @ ",round(Int,t),"s. Left: ",round(Int,(nStep-s)*t/(s-1)),"s, ν = ",1000mean(ν),"Hz")
        end
    end
    return νAll
end


##########################
### NETWORK SIMULATION ###
##########################
ParaString = ("N=$N\_C=$C\_tSim=$tSim\_tWarmup=$tWarmup\_J=$J\_tauRef=$τRef\_dV=$dV\_dt=$dt\_pertsize=$pertsize\_seed=$seed") ;println(ParaString)
ParaStringVec = ("N=$N\_C=$C\_J=$J\_tauRef=$τRef\_dV=$dV\_seed=$seed")
dir = "dataRatePert/";~iswritable(dir)?mkdir(dir):true
FileOut = dir*ParaString*"rateDevPerturbationDirectionHz.dat"

if ~isfile(FileOut)

    println("generating topology")
    JMatrix = makeConnectivity(N,C,f,g,J,seed)
    JMatrix2 = JMatrix.^2

    ### WARMUP ###
    ν0 = getSelfconsistentRateLIF(τ_m,Vth,Vre,g,f,J,C,μ0,τRef,dV,Vl)
    println("warmup with ", nStepWarmup, " steps at initial rate: ",1000ν0, " Hz")
    νInit = ν0 * ones(N)# + 0.001 * randn(N)) # start with random rates near mean field prediction
    μ0Perturbed =  μ0*ones(N)
    νWarmup = runSimulation(nStepWarmup,sPrint,dt,νInit,μ0Perturbed) 
    ν0WarmupTimeaverage = mean(νWarmup[:,round(2end/3):end],2)
    ν0WarmupPopulation = mean(νWarmup,1)

    ### SIMULATE WITH PERTURBATION ###
    vReal = getMaxRealEigenvector(JMatrix,J,μ0,τ_m,ν0,C,f,g,Vth,Vre,ε,τRef,dV,Vl,dir*ParaStringVec)
    println("simulate network with perturbation in the direction of the largest eigenvector")
    νInit = νWarmup[:,end]

    μ0Perturbed = μ0 + pertsize * vReal / norm(vReal)
    println("μ0Perturbed",μ0Perturbed[1:100])
    println("vReal",vReal[1:100])
    νPerturbed = runSimulation(nStepSim,sPrint,dt,νInit,μ0Perturbed)
    νDev = νPerturbed - mean(νPerturbed[:,end],2) * ones(size(νPerturbed,2))'
    νDevPerturbationDirectionHz = νDev' * vReal * 1000

    ### SAVE RESULTS ###
    println("simulation finished, saving results")
    writedlm(dir*ParaString*"nu0WarmupTimeaverage.dat",ν0WarmupTimeaverage)
    writedlm(dir*ParaString*"nu0Warmup.dat",mean(ν0WarmupTimeaverage))
    writedlm(dir*ParaString*"nu0WarmupPopulationrate.dat",ν0WarmupPopulation)
    writedlm(FileOut,νDevPerturbationDirectionHz)
    writedlm(dir*ParaString*"rateWarmup.dat",mean(νWarmup,1))
    writedlm(dir*ParaString*"rateDev.dat",mean(νDev,1))
    writedlm(dir*ParaString*"ratePerturbed.dat",mean(νPerturbed,1))
    println("saving finished, everything finished")
else
    println("Results already exists: ",FileOut)
end
end
