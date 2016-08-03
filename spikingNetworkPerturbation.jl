##################################################################
###	spiking network of Ostojic Nature Neuroscience 2014	##
###	Copyright Rainer Engelken, 2015 rainer@nld.ds.mpg.de	##
##################################################################

# 1. warmup of network
# 2. save warmup rates (trow away beginning of warmup)
# 3. calculate stability matrix G
# 4. perturbe every 100*ϕ ms for 2ms along least stable direction of G
# 5. readout the perturbed: save spiketimes with corresponding weighted
# 6. average over perturbations

using Distributions
include("makeConnectivity.jl")
include("getMaxRealEigenvector.jl")
include("lifMFT_Richardson.jl")
include("weightedhist.jl")

hostn = readall(`hostname`);mypid=getpid();print("running pid ",mypid," on ",hostn)
function spikingNetworkPerturbation(N,C,nPert,J,τ_ref ,τ_delay,seedIC,seedTopo,pertsize)
## This julia function is an event-based implementation of Nature Neuro Ostojic 2014
############ Input parameters:###################
# N: number of neurons in network		#
# C: number of synapses				#
# TC: simulation time in s			#
# J: synaptic strength in mV			#
# τ_delay: syn. delay in s			#
# τ_ref : refractory period in s		#
# seedTopo/seedIC for topology and IC		#
# pertsize: norm of perturbation vector		#
#################################################
# example:  Data = spikingNetworkPerturbation(10000,1000,10,0.8,0.0005,0.00055,1)
#############################
# neuron network parameters #
#############################
f = 0.8 # fraction of excitatory cells
g = 5.0 # relative strength of inhibition
Vre = 10  # threshold voltage
Vth = 20 # reset of voltage
τ_m = 20/1000# membrane time constant ( seconds )
μ0 = 24 # external input
Iext = μ0-Vth # external input beyond rheobase current (mV)
println("N=",N," C=",C," nPert=",nPert," J=",J," τ_ref =",τ_ref ," τ_delay=",τ_delay," seedIC=",seedIC," seedTopo=",seedTopo)
# simulation parameters
TW = 2 # warmup time in s.
TWthrowAway = TW/4 # don't include beginning in mean rate calculation
dV = 0.1 # step size in voltage integration
Vl = -200 # lower bound for integration
vReal = zeros(N)
ε = 1e-5 
ParaStringVec = ("N=$N\_C=$C\_J=$J\_tauRef=$τ_ref \_seedTopo=$seedTopo") #saving
ParaString =    ("N=$N\_C=$C\_nPert=$nPert\_tauRef=$τ_ref \_J=$J\_pertsize=$pertsize\_tauDelay=$τ_delay\_seedIC=$seedIC\_seedTopo=$seedTopo") #saving Hash of File
dir = "dataSpikingPert/";~iswritable(dir)?mkdir(dir):true
FileOut = dir*ParaString*"FileOut.dat"
if  ~isfile(FileOut)
    ###################################################################
    println(" initialization, network setup, memory preallocation ") ##
    ###################################################################
    Iext = Iext*ones(N) # vector because Iext is heterogenous during perturbation
    I = Vre - Iext - Vth # change of variables for convenience
    ω = - 1./log( - Iext./I )  # phase velocity of LIF
    c = - J ./ I # effective coupling in phase response courve
    ϕTh = exp( - 1./ω) # treshold in ϕ to check wether neuron spikes
    ## set up mixed network topology with fixed indegree
    JMatrix = makeConnectivity(N,C,f,g,J,seedTopo)
    JMatrix2 = JMatrix.^2
    srand(seedIC) # setting seed of random number generator different for IC
    if isfile("PofV/MFT_C$C\_J$J\Vsample.dat")
        Vsample = readdlm("PofV/MFT_C$C\_J$J\Vsample.dat")
        ϕSample = -ω[1]*log((Vsample-Vth-Iext[1])/(Vre-Vth-Iext[1]))
        ϕ = ϕSample[rand(1:length(ϕSample),N)]
        rateMFT = readdlm("PofV/MFT_C$C\_J$J\Rate.dat")
    else
        rateMFT = 78; ϕ = randn(N )-1; ϕ[ϕ.>1] = -1; # neurons  initial phases
    end
    
    srand(seedIC)
    nQstart = round(Int,rateMFT[1]*τ_delay*N); 
    spikeIdxQ = rand(1:N,nQstart); spikeTimeQ = τ_delay*sort(rand(nQstart))
    
    #preallocate memory and initialize some variables
    LastSpikeOfThisNeuron = - 2τ_ref *ones(N)
    nDisplay = 10; nDisplayW = 10000
    NSyncSpikesTimes = Float64[]; NSyncSpikes = Int[]
    trainTimes = Float64[]; trainNeuron = Int[]; nSpikes = 0; nUniqueSpiketime = 0
    SynchronyWarning = false; measureIdx = 0; s = 0; t  = 0
    
    # preallocate some stuff for perturbation:
    PerturbationPeriod = 0.1 * (1 + sqrt(5)) / 2 # in seconds. Using golden ratio to avoid resonance ;-)
    PerturbationTime  = 0.002 # in s
    ωUnpert = ω # perturbation in Iext also affect phase velocity ω and ϕTh
    IUnpert = I
    ϕThUnpert = ϕTh
    pertDefined = false
    Pert = false

    pertNo = 0 # perturbation number
    println("PerturbationPeriod = ",PerturbationPeriod)
    println("PerturbationTime = ",PerturbationTime)
    TWrand = 0.5rand() # add random time to avoid any artefacts from start
    writedlm(dir*ParaString*"TWrand.dat",TWrand) # save this time
    TW = TW + TWrand 
    println("TWrand = ",TWrand)

    #################################
    println("network simulation ") ##
    #################################
    tStart = time_ns()
    while pertNo < nPert
        s += 1 # s count number of events
        # check whether it is time for a perturbation:
        if t  > TW
            if !pertDefined # if perturbation fector is not yet defined, define perturbation vector
                ν0 = 1000*getSelfconsistentRateLIF( τ_m *1000,Vth,Vre,g,f,J,C,μ0,τ_ref*10000,dV,Vl) #1000 for Hz
                println("ν0 = ",ν0)
                vReal = getMaxRealEigenvector(JMatrix,J,μ0,τ_m,ν0,C,f,g,Vth,Vre,ε,τ_ref,dV,Vl,dir*ParaStringVec)
                println("simulate network with perturbation in the direction of the largest eigenvector")
                IextP = Iext + pertsize*vReal/norm(vReal) 
                println("minimum(IextP)",minimum(IextP))
                if minimum(IextP) <= 0 # this should never happen
                    IextP[IextP.<=0]=0.001
		    println("perturbation too strong")
                end
                IextP = IextP[:] # turn into vector otherwise Nx1 array
                IPert = Vre - IextP - Vth
                ωPert = - 1./log( - IextP./IPert )
                ϕThPert = exp( - 1./ωPert)
                pertDefined = true
                tStartSim = time_ns()
            end


            tMod = mod(t -TW,PerturbationPeriod)
            if tMod < PerturbationTime 
                if ~Pert
                    pertNo = pertNo + 1
                    if mod(Pert,nDisplay)==0
                        println("pertNo = ",pertNo)
                        tNow = (time_ns()-tStartSim)/1e9;# RateNow = nSpikes /t  /N;#print("\r")
                        RateNow = sum(trainTimes.>TW)/N / (t  - TW)
                        println(round(pertNo*100/nPert),"% done after ",round(Int,tNow)," s, Tsim = ",round(100*t )/100,"s Remaining: ",round(Int,(nPert - pertNo)*tNow/pertNo),"s ",s," events"," @ " ,round(RateNow,3)," Hz")
                    end
                    ω = ωPert
                    I = IPert
                    ϕTh = ϕThPert
                    Pert = true
                end
            else
                if Pert
                    ω = ωUnpert
                    I = IUnpert
                    ϕTh = ϕThUnpert
                    Pert = false
                end
            end
        else
        
            if mod(s,nDisplayW) == 0
                tNow = (time_ns()-tStart)/1e9
                RateNow = sum(trainTimes.>TWthrowAway)/N / (t  - TWthrowAway)
                println(round(t *100/TW),"% of TW done after ",round(Int,tNow)," s, Tsim = ",round(100*t )/100,"s Remaining: ",round((TW - t )*tNow/t ),"s ",s," events"," @ " ,round(RateNow,3)," Hz")
            end
        end
        ## Find next event in the network
        dtAll = ( 1 - ϕ )./ ω
        dtMin = minimum(dtAll)

        j = find(dtAll .== dtMin)         #detecting all sync spikes
        if (length(j) > N/10 && !SynchronyWarning ); SynchronyWarning = true; println("network synchrony");end # throw warning if more than 10% of network spikes sync.
            dt = dtMin #( 1 - ϕMax )./ w[j[1]]  # calculate next spike time xy CHECK!! FIXME

      
            if ~isempty(spikeTimeQ) && (spikeTimeQ[1] == dt); println("next phase update and next spike coincide, this should almost never happen!");end # this should almost never happen
            if (~isempty(spikeTimeQ) && spikeTimeQ[1]<= dt)
                dt = spikeTimeQ[1]  # in case next measurement is before next spike
      
                NextEventIsPhaseUpdate = true
            else
                NextEventIsPhaseUpdate = false
            end
            ## evolve phases till next event (spiketime, phase update or phase measurement)
            NotInTauRef = τ_ref  .<= (t  - LastSpikeOfThisNeuron)
            ϕ[NotInTauRef] +=  ω[NotInTauRef] .* dt
            t  += dt * τ_m  # update time
            NotInTauRefPost =  τ_ref .<=collect(t  - LastSpikeOfThisNeuron)
            ϕRefOver =  NotInTauRef $ NotInTauRefPost
            if sum(ϕRefOver)>0
                dt2 = (t  - LastSpikeOfThisNeuron[ϕRefOver] - τ_ref )/τ_m
                ϕ[ϕRefOver] += ω[ϕRefOver] .* dt2
            end
            spikeTimeQ -=  dt  # update time of Q which is always relative to now in units of τ_m like dt
                
                ## if next event is spike, save spiketimes and put spike into queue.
            if ~(NextEventIsPhaseUpdate) 
                append!(trainTimes,t *ones(length(j)))
                append!(trainNeuron,j)
                push!(NSyncSpikes,length(j))
                push!(NSyncSpikesTimes,t )
                LastSpikeOfThisNeuron[j] = t 
                append!(spikeIdxQ,j)
                append!(spikeTimeQ,τ_delay/τ_m*ones(length(j)))
                nSpikes += length(j)
                nUniqueSpiketime += 1
                ϕ[j] = 0 #ϕReset
                        
                        ## if next event is phase update , get spiking neuron index, find postsynaptic neurons, sum PSP and propagate phase
            elseif NextEventIsPhaseUpdate
                j = Int64[]
                while ~isempty(spikeTimeQ) && spikeTimeQ[1]==0
                    push!(j,spikeIdxQ[1])
                    shift!(spikeIdxQ)
                    shift!(spikeTimeQ)
                end
                Jtemp = sum(JMatrix[:,j],2) # sum PSP of all spiking neurons.
                Jtemp[~NotInTauRefPost,1] = 0 # neurons in τ_ref  receive no input
                Jtemp=Jtemp./I # every I entry can be different
                postUnique = find(Jtemp.!=0) # find unique postsynaptic neurons
                InLog = exp( - ϕ[ postUnique ]./ω[ postUnique ]) + Jtemp[postUnique]
                ϕLogPos = InLog .> ϕTh[postUnique] # find neurons whose phase get beyond 1
                ϕ[postUnique[ϕLogPos]] = - ω[ postUnique[ϕLogPos] ].*log(InLog[ϕLogPos]) #update postsyn. phases
                ϕ[postUnique[~ϕLogPos]] = 1 # set ϕ to 1 in case excitatory input kicks it beyond threshold

            else
                error("Neither spike nor phase update nor phase measure")
            end
        end
        ##########
        # output #
        ##########
        "postprocessing"
        meanRate = sum(trainTimes.>TW)/N / (t  - TW);tNow = (time_ns()-tStart)/1e9;
        #print("\r");
        println("\n","#"^20, "Simulation finished \n Mean rate : ",  round(meanRate,2) , "Hz from  ", nSpikes, " spikes after ",round(Int,tNow),"s", ", writing rate to ",FileOut)
        writedlm(FileOut,meanRate)
        UniqueSpikeFraction = nUniqueSpiketime/nSpikes
        writedlm(dir*ParaString*"UniqueSpikeFraction.dat",UniqueSpikeFraction)
        trainTimes2 = trainTimes[trainTimes.>TW]
        trainNeuron2 = trainNeuron[trainTimes.>TW]
        tBins1ms = 0:0.001:PerturbationPeriod
        
        w = vReal[trainNeuron2];
        histUni1ms = weightedhist(mod(trainTimes2-TW,PerturbationPeriod),tBins1ms)
        histWeighted1ms = weightedhist(mod(trainTimes2-TW,PerturbationPeriod),tBins1ms;weights=w)
        writedlm(dir*ParaString*"histTimes1ms.dat",histUni1ms[1])
        writedlm(dir*ParaString*"histUni1ms.dat",histUni1ms[2])
        writedlm(dir*ParaString*"histWeighted1ms.dat",histWeighted1ms[2])
        MeanNeuronNumberPerSpiketime = mean(NSyncSpikes[NSyncSpikesTimes.>TW])
        writedlm(dir*ParaString*"MeanNeuronNumberPerSpiketime.dat",MeanNeuronNumberPerSpiketime)
        writedlm(dir*ParaString*"NSyncSpikes.dat",NSyncSpikes)
        writedlm(dir*ParaString*"NSyncSpikesTimes.dat",NSyncSpikesTimes)
    else
        println("Simulation was already done, loading results from ",FileOut)
        Rate = readdlm(FileOut);meanRate = Rate[1] # load array to float
    end
    println("#### End of simulation #####")
    return meanRate
end
