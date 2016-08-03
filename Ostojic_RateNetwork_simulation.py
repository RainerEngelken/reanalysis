# -*- coding: utf-8 -*-
#
# Ostojic_RateNetwork_simulation.py
#
# Copyright (C) 2015 Farzad Farkhooi
#
# This code is provided based on GNU GPLv3 license.
# Please concult with farzad@bccn-berlin.de in case of problem running it.

from ricchiardi_lb import ricciardi_lb as nu
from scipy.optimize import fsolve
from scipy.sparse import csr_matrix
import numpy as np

def mu(tau_m, nu0, C, J, g ,f, mu0):
	return tau_m*nu0*C*J*(f - (1-f)*g)+mu0

def sigma2(tau_m, nu0, C, J, g ,f):
	return tau_m*nu0*C*J**2*(f+(1-f)*g**2)

# Ostojic parameters
tau_m =	 0.020
tau_r =	0.0005
mu_0 =	 24.
f =	 0.8
g =	 5.0
v_th =	 20.
v_re =	 10.
J = 	0.5375
C = 	1000.
N = 10000


SEED=1234

N_rec = 10000 

# In this version: The connectivity is generated in BRIAN

from brian import *
def gen_conn_brian(Seed):
	eqs = "dv/dt = ((-v+ mu_0*mV)) /(tau_m*second)  : volt"
	P = NeuronGroup(N, eqs, threshold=v_th*mV, reset=v_re*mV , refractory=tau_r*ms)
	Pe = P.subgroup(int(round(N*f)))
	Pi = P.subgroup(int(round(N*(1-f))))
	
	K = Connection(P, P, state='v')
	K.connect_random(Pe,P,weight=1000*J*mV,sparseness=C/N,seed=Seed,fixed=True)
	K.connect_random(Pi,P,weight=-g*J*1000*mV,sparseness=C/N,seed=Seed,fixed=True)
	
	K2 = K.W.multiply(K.W) 
	K = csr_matrix(K.W)
	return K.T,  K2.T


Conn, Conn2 = gen_conn_brian(SEED)


# Here the equilibrium rate is calc to start with!
# This step is not needed!
mu_x = lambda x :  mu(tau_m, x, C, J, g, f, mu_0)
sigma_x = lambda s : sqrt( sigma2(tau_m, s, C, J, g ,f) ) 
nu_eq = fsolve(lambda x: x-nu(mu_x(x), sigma_x(x), tau_m, tau_r, v_th, v_re, -Inf), 13.82)
m = np.ones(N)*mu_x(nu_eq)
s = np.ones(N)*sigma_x(nu_eq)
r = np.random.lognormal(log(nu_eq), .1, N)
vnu=np.vectorize(nu)
x = np.array([])



t = 0.0
dt = 0.0001
delta = 55


while t < 20.:
	r += (dt/tau_m)*(-r + vnu(m, s, tau_m, tau_r, v_th, v_re, -Inf))
	m = mu_0 + tau_m* Conn.dot(r) 
	s = sqrt(tau_m* Conn2.dot(r))
	t += dt
	x= np.append(x, r[0:N_rec])
	if int(t/dt) % 100. == 0:
		print t	, mean(mean(r))
	
rates = x.reshape(len(x)/N_rec, N_rec)

rates.dump('rates_C'+str(C)+'_N'+str(N)+'_J'+str(J)+'_seed'+str(SEED)+'_T'+str(t))

