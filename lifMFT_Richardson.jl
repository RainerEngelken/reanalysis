function LIF0dVref(tau,sig,μ,Vth,Vre,tauRef,dV,Vl)
#solves the LIF assuming Gaussian white noise. Method from Richardson, PRE(2007)
    V = Vl:dV:Vth		 # set up the lattice for the discretized integration
    n = length(V)
    kre = round(Int,(Vre-Vl)/dV) 	# Vre must fall on a lattice point!
    G = (V-μ)/sig^2
    A = exp(dV*G) 		# modified Euler method (Phys. Rev. E 2007, Eqs A1-A6)
    B = (A-1)./(sig^2*dV*G)
    B[G.==0] = 1/sig^2
    j0 = ones(n) 		# allocate mem for scaled current and probability density
    p0 = zeros(n) 		# initial conditions at j0(Vth)=1, p0(n)=0 already.
    for k = n:-1:2 		# integrate backwards from threshold
        #j0[k-1] = j0[k] - round(Int,k.==kre+1)
        k==kre+1? j0[k-1] = j0[k] - 1 : j0[k-1] = j0[k]
        p0[k-1] = p0[k]*A[k] + dV*B[k]*tau*j0[k]
    end
    r0::Float64 = 1/(tauRef+dV*sum(p0)) # steady-state firing rate (in kHz)
end


function get_γ_μ(τ,σ,μ,Vth,Vre,ε,τ_ref,dV,Vl) # derivative of dF/dμ
# calculates derivative of transfer function w.r.t mean
    ν0Plus = LIF0dVref(τ,σ,μ+ε,Vth,Vre,τ_ref,dV,Vl)
    ν0Minus = LIF0dVref(τ,σ,μ-ε,Vth,Vre,τ_ref,dV,Vl)
    dν0dMu = 1000(ν0Plus-ν0Minus)/ε/2 	# equation 13 in Ostojic 2014
end

function get_γ_σ2(τ,σ,μ,Vth,Vre,ε,τ_ref,dV,Vl) # derivative of dF/dσ2
# calculates derivative of transfer function w.r.t variance
    σ2 = σ^2*sqrt(2)
    ν0Plus = LIF0dVref(τ,sqrt(σ2+ε)/sqrt(2),μ,Vth,Vre,τ_ref,dV,Vl)
    ν0Minus = LIF0dVref(τ,sqrt(σ2-ε)/sqrt(2),μ,Vth,Vre,τ_ref,dV,Vl)
    dνdMu = 1000(ν0Plus-ν0Minus)/ε/2 	# equation 14 in Ostojic 2014
end

function getSelfconsistentRateLIF(τ_m,Vth,Vre,g,f,J,C,μ0,τ_ref,dV,Vl)
# calculates self-consistent firing rate assuming Gaussian white noise currents
    ν0 = J/(-μ0/τ_m/(f-(1-f)*g))
    ν0Up = 2 # 2000 Hz upper bound
    ν0Down = 0
    ν0Out = 0
    pR = 0.00001
    i = 0
    μ = μ0 + τ_m*ν0*C*J*(f-(1-f)*g)
    σ = sqrt(τ_m*ν0*C*J^2*(f+(1-f)*g^2)/2)
    while (abs(1-ν0/ν0Out)>pR && i<45)
        (ν0Up>ν0 && i<40) || error("ν0>ν0Up or i>40 in getSelfconsistentRateLIF")
        ν0Out<ν0 ? ν0Up = ν0 : ν0Down = ν0
        ν0 = (ν0Down+ν0Up)/2
        μ = μ0 + τ_m*ν0*C*J*(f-(1-f)*g)
        σ = sqrt(τ_m*ν0*C*J^2*(f+(1-f)*g^2))/sqrt(2)
        ν0Out = LIF0dVref(τ_m,σ,μ,Vth,Vre,τ_ref,dV,Vl)
        i+=1
    end
return ν0Out
end
