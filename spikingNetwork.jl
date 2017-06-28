##################################################################
###    	spiking network of Ostojic Nature Neuroscience 2014	##
### 	Copyright Rainer Engelken, 2016 rainer@nld.ds.mpg.de	##
##################################################################
using Distributions
hostn=readall(`hostname`);mypid=getpid();print("running pid ",mypid," on ",hostn)
function spikingNetwork(N,C,TC,J,τ_ref ,τ_delay,seed)
## This julia function is in event-based implementation of Nature Neuro Ostojic 2014
# input parameters:
# N: number of neurons in network
# C: number of synapses
# TC: simulation time in s
# J: synaptic strength in mV
# τ_delay: syn. delay in s
# τ_ref : refractory period in s
# seed for topology and initial conditions
# example:  spikingNetwork(10000,1000,10,0.8,0.0005,0.00055,1)
#############################
# neuron network parameters #
#############################
f = 0.8 # fraction of excitatory cells
g = 5 # relative strength of inhibition
Vr = 10  # threshold voltage
Ith = 20 # reset of voltage
τ_m = 20/1000 # membrane time constant ( seconds )
Iext = 4 # external input beyond rheobase current (mV)
println("N=",N," C=",C," TC=",TC," J=",J," τ_ref =",τ_ref ," τ_delay=",τ_delay," seed=",seed)
plotFlat = false
# simulation parameters
ϕN = N;ϕN>N ?ϕN=N : ϕN=ϕN  # number of phases to be measured (smaller than N)
nMeasure = 1000 # number of phase measurements
TW = 2 # warmup time in s. Phases will be measured between TW and TC
ParaString = ("N=$N\_C=$C\_TC=$TC\_J=$J\_tauRef=$τ_ref\_tauDelay=$τ_delay\_seed=$seed") #saving Hash of File
dir = "dataSpikingNet/";~iswritable(dir)?mkdir(dir):true
FileOut = dir*ParaString*"FileOut.dat"
if  ~isfile(FileOut)
    ###################################################################
    println(" initialization, network setup, memory preallocation ") ##
    ###################################################################
    Ieff = Vr - Iext - Ith
    ω = - 1/log( - Iext/Ieff )  # phase velocity LIF
    c = - J/Ieff # effective coupling in PRC
    ϕTh = exp( - 1/ω)

    ## set up mixed network topology with fixed indegree
    TimeStmpBeforeTopo = tic()
    srand(seed)
    Ni = round(Int,N - N*f)	
    Ne = N - Ni
    Ce = round(Int,f*C)
    Ci = round(Int,C - C*f)
    JMatrix = zeros(N,N)
    for n=1:N
      for m = sample([1:min(n-1,Ne);(n+1):Ne],Ce,replace=false) # using Distributions
	  JMatrix[n,m] -= c
      end
      for m = sample([Ne+1:n-1;max(n+1,Ne+1):N],Ci,replace=false)
	  JMatrix[n,m] += g*c
      end
    end   #JMatrix = sparse(JMatrix) # this takes time and is only recommended for large N


    TimeForTopo = toq(); println("Building topology took " ,TimeForTopo," s")
    srand(seed)

    if isfile("PofV/MFT_C$C\_J$J\Vsample.dat")#optionally use p(V) from MFT for IC
        Vsample = readdlm("PofV/MFT_C$C\_J$J\Vsample.dat")
        ϕSample = -ω[1]*log((Vsample-Ith-Iext)/(Vr-Ith-Iext))
        ϕ = ϕSample[rand(1:length(ϕSample),N)]::Array{Float64,1}
        rateMFT = readdlm("PofV/MFT_C$C\_J$J\Rate.dat")
    else
        rateMFT = 15.; ϕ = randn(N )-1; ϕ[ϕ.>1] = -1; # neurons  initial phases
    end
    
    srand(seed)
    nQstart = round(Int,rateMFT*τ_delay*N); println("nQstart=",nQstart)
    spikeIdxQ = rand(1:N,nQstart); spikeTimeQ = τ_delay*sort(rand(nQstart));    #spikeIdxQ = Int[]; spikeTimeQ = Float64[];
    measureTimes = TW + TC*sort(rand(nMeasure)) # measure at random times
    #measureTimes = TW + (TC - TW)*[1:nMeasure]/nMeasure # measure at equally space times

    #preallocate memory and initialize some variables
    LastSpikeOfThisNeuron = - 2*τ_ref *ones(N)
    nDisplay = 10000;RateNow = 0
    NSyncSpikesTimes = Float64[]; NSyncSpikes = Int[]
    trainTimes = Float64[]; trainNeuron = Int[]; nSpikes = 0; nUniqueSpiketime = 0; ϕAll = zeros(ϕN,nMeasure); tAll = zeros(nMeasure)  
    SynchronyWarning = false; measureIdx = 0; s = 0; t  = 0.
    #################################
    println("network simulation ") ##
    #################################
    tStart = time_ns()
    while t  <(TC+TW) && RateNow<10000
	s += 1
	if mod(s,nDisplay) == 0
	tNow = (time_ns()-tStart)/1e9; RateNow = nSpikes /t  /N;#print("\r")
	println(round(t *100/(TC+TW)),"% done after ",round(Int,tNow)," s, Tsim = ",round(100*t )/100,"s Remaining: ",round(Int,((TC+TW) - t )*tNow/t ),"s ",s," events"," @ " ,round(RateNow,3)," Hz")
	end

	## Find next event in the network
	ϕMax = maximum( ϕ ) # find next spiking neuron j
	j = find(ϕ .== ϕMax)         #detecting all sync spikes
	if (length(j) > N/10 && !SynchronyWarning ); SynchronyWarning = true; println("network synchrony");end # throw warning if more than 10% of network spikes sync.

	dt = ( 1 - ϕMax )/ ω  # calculate next spike time
	if ~isempty(spikeTimeQ) && (spikeTimeQ[1] == dt); println("next phase update and next spike coincide");end # this should almost never happen
	if (~isempty(spikeTimeQ) && spikeTimeQ[1]<= dt)
	    dt = spikeTimeQ[1]  # in case next measurement is before next spike
	    NextEventIsPhaseUpdate = true
	else
	    NextEventIsPhaseUpdate = false
	end
	if ~isempty(measureTimes) && measureTimes[1] <= t  + dt*τ_m
	    NextEventIsPhaseUpdate = false
	    NextEventIsPhaseMeasure = true
	    dt = (measureTimes[1] - t )/τ_m
	    shift!(measureTimes)
	else
	    NextEventIsPhaseMeasure = false
	end
	## evolve phases till next event (spiketime, phase update or phase measurement)
	NotInTauRef = τ_ref  .<= (t  - LastSpikeOfThisNeuron)
	ϕ[NotInTauRef] +=  ω * dt 
	t  += dt * τ_m  # update time
	NotInTauRefPost =  τ_ref .<= t  - LastSpikeOfThisNeuron
	ϕRefOver =  NotInTauRef $ NotInTauRefPost
	if sum(ϕRefOver)>0
	    dt2 = (t  - LastSpikeOfThisNeuron[ϕRefOver] - τ_ref )/τ_m
	    ϕ[ϕRefOver] += ω * dt2 
	end
	spikeTimeQ -=  dt  # update time of Q which is always relative to now in units of τ_m like dt
	
	## if next event is spike, save spiketimes and put spike into queue.
	if ~(NextEventIsPhaseUpdate) && ~(NextEventIsPhaseMeasure)
	    append!(trainTimes,t *ones(length(j)))
	    append!(trainNeuron,j)
	    push!(NSyncSpikes,length(j))
	    push!(NSyncSpikesTimes,t )
	    LastSpikeOfThisNeuron[j] = t 
	    append!(spikeIdxQ,j)
	    append!(spikeTimeQ,τ_delay/τ_m*ones(length(j)))
	    nSpikes += length(j)
	    nUniqueSpiketime += 1
	    ϕ[j ] = 0 #ϕReset
	
	## if next event is phase update , get spiking neuron index, find postsynaptic neurons, sum PSP and propagate phase
	elseif NextEventIsPhaseUpdate
	    j=Int64[]
	    while ~isempty(spikeTimeQ) && spikeTimeQ[1]==0
		push!(j,spikeIdxQ[1])
		shift!(spikeIdxQ)
		shift!(spikeTimeQ)
	    end
	    Jtemp = sum(JMatrix[:,j],2) # sum PSP of all spiking neurons.
	    Jtemp[~NotInTauRefPost,1] = 0 # neurons in τ_ref  receive no input
	    postUnique = find(Jtemp.!=0) # find unique postsynaptic neurons
	    InLog = exp( - ϕ[ postUnique ]/ω) + Jtemp[postUnique]
	    ϕLogPos = InLog .> ϕTh # find neurons whose phase get beyond 1
	    ϕ[postUnique[ϕLogPos]] =  - ω*log(InLog[ϕLogPos]) #update postsyn. phases
	    ϕ[postUnique[~ϕLogPos]] = 1 # set ϕ to 1 in case excitatory input kicks it beyond threshold
	## in case next event is phase update measure phases
	elseif NextEventIsPhaseMeasure
	    measureIdx += 1
	    ϕAll[:,measureIdx] = ϕ[1 : ϕN]
	    tAll[measureIdx] = t 
        else
           throw("Neither spike nor phase update nor phase measure")
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
UniqueSpikeFraction=nUniqueSpiketime/nSpikes
writedlm(dir*ParaString*"UniqueSpikeFraction.dat",UniqueSpikeFraction)
χ = std(mean(ϕAll,1),2)/mean(std(ϕAll,2),1)
println( "χ = ",  χ)
writedlm(dir*ParaString*"chi.dat",χ)
cv = NaN*zeros(N); rates = NaN*zeros(N)
for n = 1:N
    trainTemp = trainTimes[trainNeuron.==n]
    if length(trainTemp)>6 # to few spikes contaminate the statistics
       ISItemp = diff(trainTemp)
       cv[n] = std(ISItemp)/mean(ISItemp)
       rates[n] = length(trainTemp)/(trainTemp[end]-trainTemp[1])
    end
end
cvMean = mean(cv[!isnan(cv)])
writedlm(dir*ParaString*"cvMean.dat",cvMean)
println( "Mean cv = ",  cvMean)


trainTimesShort = trainTimes[trainTimes.<1+TW]
trainNeuronShort = trainNeuron[trainTimes.<1+TW]
writedlm(dir*ParaString*"trainTimes.dat",trainTimesShort)
writedlm(dir*ParaString*"trainNeuron.dat",trainNeuronShort)
MeanNeuronNumberPerSpiketime = mean(NSyncSpikes[NSyncSpikesTimes.>TW])
VoltageAll = (Vr-Ith-Iext)*(exp(-ϕAll/ω))+Ith+Iext;
MeanVoltage = mean(VoltageAll);
MedianVoltage = median(VoltageAll);
writedlm(dir*ParaString*"MeanNeuronNumberPerSpiketime.dat",MeanNeuronNumberPerSpiketime)
writedlm(dir*ParaString*"MeanVoltage.dat",MeanVoltage)
writedlm(dir*ParaString*"MedianVoltage.dat",MedianVoltage)
minV = minimum(VoltageAll[:]); maxV = maximum(VoltageAll[:]); 
histVx = collect(round(Int,minV):round(Int,maxV))
histV = hist(VoltageAll[:], histVx)
writedlm(dir*ParaString*"histVx.dat",collect(histV[1]))
writedlm(dir*ParaString*"histVy.dat",histV[2])
writedlm(dir*ParaString*"maxV.dat",maxV)
writedlm(dir*ParaString*"minV.dat",minV)
writedlm(dir*ParaString*"cvAll.dat",cv)
writedlm(dir*ParaString*"rates.dat",rates)
#if seed == 1; writedlm(dir*ParaString*"phiAll.dat",ϕAll) ; end
else
println("Simulation was already done, loading results from ",FileOut) 
Rate = readdlm(FileOut);meanRate = Rate[1] # load array to float
end
println("%%%% End of simulation #####")
return meanRate
end
