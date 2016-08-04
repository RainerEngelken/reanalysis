# reanalysis
Code for reanalysis of “Two types of asynchronous activity in networks of  excitatory and inhibitory spiking neurons”

## Installation
first install julia 0.4 from [julialang.org/downloads/](http://www.julialang.org/downloads/)
Then start julia and  install the modules Distributions and PyPlot by the following command:
~~~julia
julia> Pkg.add("Distributions"); Pkg.add("PyPlot")
~~~

##Usage
`example_script_spikingNetwork.jl` generates a quick low resolution version of figure 1a for N=10000, K=1000.
It runs `spikingNetwork.jl`, which is an numerically exact event-based implementation of the spiking network in
Ostojic, S. Two types of asynchronous activity in networks of excitatory and inhibitory spiking neurons. Nat Neurosci 17, 594–600 (2014).

For higher precision, use longer simulation time and more coupling strengths. 

`spikingNetwork.jl` returns the population firing rate and saves it to `LOCALDIRECTORY/dataSpikingNet/FileOut.dat`
Similarly, spike times and indices of spiking neurons are stored at `trainTimes.dat`, `trainNeuron.dat`.
Further, it stores the mean and median voltage, a histogram of all voltages including its minimum and maximum, the coefficient of variation for all neurons and all individual time-averaged firing rates. 
chi.dat saves for each simulation the common synchrony measure Χ (Hansel, D. & Mato, 2003), which requires no binning.
To measures how many exactly synchronous spikes occur in the network, `UniqueSpikeFraction.dat` gives the number of unique spike times divided by the number of total network spikes, so a unique spike fraction of 0.1 means that on average there are 10 neurons spiking at each spike time.
One reason for the overabundance of exactly synchronous spikes is the hard threshold of the LIF model: 
When an excitatory neuron spikes, it can push all postsynaptic neurons that are less then J below Vth into threshold.

Figures 1b, 1c and S1 were generated using the scripts `Ostojic_RateNetwork_simulation.py` and `Ostojic_SpikingNetwork_simulation_NEST.py`.

Figure 1d were done using `ratePerturbation.jl` and `spikingNetworkPerturbation.jl`. 
After a warmup, the external current μ₀ is perturbed into the least stable direction of the linearised dynamics of the rate network.
This direction is calculated using `getMaxRealEigenvector.jl`. 
The input-output function of the single units for the rate network was solved using the efficient Richardson method (Richardson PRE 2007) using
the code in `lifMFT_Richardson.jl`.
The resulting rate perturbation was projected on the least stable direction of the linearized rate dynamics.

In the case of the spiking network, exactly the same network parameters and topology were used. The exactly same perturbation was applied to
$\mu_0$. To obtain good statistics, this was repeated fo over 1.42 Mio trials. In each trial, the resulting rate vector (binned into 1ms bins) was projected onto the least stable direction.
The average rate deviation after the end of the external perturbation is depicted in figure 1d.
Different perturbation durations and perturbation strength had similar effect.

All results from figure 2 are generated with `spikingNetwork.jl`
Please report any bugs or problems running the code to rainer <at> nld <dot> ds <dot> mpg <dot> de (add @ before nld and dots otherwise).
