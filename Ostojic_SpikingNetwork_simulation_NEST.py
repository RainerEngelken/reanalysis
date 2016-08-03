# -*- coding: utf-8 -*-
#
# Ostojic_SpikingNetwork_simulation_NEST.py
#
# Copyright (C) 2015 Farzad Farkhooi 
# 
# This code is provided based on GNU GPLv3 license.
# Please concult with farzad@bccn-berlin.de in case of problem running it.


def local(src):
	""" returns GIDs in src which are local to this process """
	l = nest.GetStatus(src,'local')
	# convert to indexable 
	idx = nonzero(l)
	return array(src)[idx].tolist()



from scipy.optimize import fsolve

import nest
import nest.raster_plot
import nest.voltage_trace

from numpy import sqrt

import numpy
from numpy import exp

import time

nest.ResetKernel()

startbuild= time.time()

dt      = 0.1    # the resolution in ms
simtime = 2000.0 # Simulation time in ms
delay   = 0.55    # synaptic delay in ms

# Parameters for asynchronous irregular firing
g       = 5.0
epsilon = 0.1    # connection probability

factor = 1

order     = 2000*factor
NE        = 4*order
NI        = 1*order
N_neurons = NE+NI
N_rec     = int(N_neurons) # record from 50 neurons

mu_0 = 24.0

CE    = 800*factor  # number of excitatory synapses per neuron
CI    = 200*factor  # number of inhibitory synapses per neuron  
C_tot = int(CI+CE)  # total number of synapses per neuron

# Initialize the parameters of the integrate and fire neuron 
CMem = 20.0
tauMem = 20.0
theta  = 20.0
V_reset = 10.0
J      = 0.3*sqrt(1000./(C_tot))# postsynaptic amplitude in mV

J_ex  = J 
J_in  = -g*J_ex


nest.SetKernelStatus({"resolution": dt, 
	"print_time": True, 
	'off_grid_spiking': True,
	'total_num_virtual_procs': 4})

print "Building network"

neuron_params= {"C_m"       : CMem,
                "tau_m"     : tauMem,
                "t_ref"     : 0.5,
                "E_L"       : 0.0,
                "V_m"       : 0.0,
		"I_e"	    : mu_0,
                "V_th"      : theta,
		"V_reset"   : V_reset}

nest.SetDefaults("iaf_psc_delta", neuron_params)

nodes_ex=nest.Create("iaf_psc_delta",NE)
nodes_in=nest.Create("iaf_psc_delta",NI)

rate_node = nest.Create("iaf_psc_delta",1)
nest.SetStatus(rate_node, 'V_th', +Inf)
nest.SetStatus(rate_node, 'tau_m', .00001)

nest.SetStatus(nodes_ex, 'V_m', 10.*np.random.random(NE))
nest.SetStatus(nodes_in, 'V_m', 10.*np.random.random(NI))



voltmeter = nest.Create("voltmeter")


espikes=nest.Create("spike_detector")
ispikes=nest.Create("spike_detector")

nest.SetStatus(espikes,[{"start":0.0,"label": "Ostoijc-py-ex",
                   "withtime": True,
                   "withgid": True}])

nest.SetStatus(ispikes,[{"start":0.0, "label": "Ostoijc-py-in",
                   "withtime": True,
                   "withgid": True}])

print "Connecting devices."

nest.CopyModel("static_synapse","excitatory",{"weight":J_ex, "delay":delay})
nest.CopyModel("static_synapse","inhibitory",{"weight":J_in, "delay":delay})

nest.ConvergentConnect(nodes_ex,espikes,model="excitatory")
nest.ConvergentConnect(nodes_in,ispikes,model="excitatory")

nest.ConvergentConnect(nodes_ex,rate_node, weight=1.0, delay=dt)
nest.ConvergentConnect(nodes_in,rate_node, weight=1.0, delay=dt)
nest.Connect(voltmeter, rate_node)


print "Connecting network."
numpy.random.seed(1234)

sources_ex = numpy.random.random_integers(1,NE,(N_neurons,CE))
sources_in = numpy.random.random_integers(NE+1,N_neurons,(N_neurons,CI))

# We now iterate over all neuron IDs, and connect the neuron to
# the sources from our array. The first loop connects the excitatory neurons
# and the second loop the inhibitory neurons.


for n in xrange(N_neurons):
    nest.ConvergentConnect(list(sources_ex[n]),[n+1],model="excitatory")


for n in xrange(N_neurons):
    nest.ConvergentConnect(list(sources_in[n]),[n+1],model="inhibitory")
endbuild=time.time()

print "Simulating."

nest.Simulate(simtime)

endsimulate= time.time()


events_ex = nest.GetStatus(espikes,"n_events")[0]
rate_ex   = events_ex/simtime*1000.0/NE
events_in = nest.GetStatus(ispikes,"n_events")[0]
rate_in   = events_in/simtime*1000.0/NI

num_synapses = nest.GetDefaults("excitatory")["num_connections"]+\
nest.GetDefaults("inhibitory")["num_connections"]

build_time = endbuild-startbuild
sim_time   = endsimulate-endbuild

print "Ostojic network simulation (Python)"
print "Number of neurons :", N_neurons
print "Number of synapses:", num_synapses
print "       Exitatory  :", int(CE*N_neurons)+N_neurons
print "       Inhibitory :", int(CI*N_neurons)
print "Excitatory rate   : %.2f Hz" % rate_ex
print "Inhibitory rate   : %.2f Hz" % rate_in
print "Building time     : %.2f s" % build_time
print "Simulation time   : %.2f s" % sim_time

nest.raster_plot.from_device(espikes, hist=True)
nest.raster_plot.show()
