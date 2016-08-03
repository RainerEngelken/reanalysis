function getMaxRealEigenvector(JMatrix,J,μ0,τ_m,ν0,C,f,g,Vth,Vre,ε,τ_ref,dV,Vl,path) 
    if !isfile(path*"vReal.dat")
        JMatrix2 = JMatrix.^2
        println("calculating largest eigenvector of linearized dynamics")
        μ = μ0 + τ_m * ν0 * C * J * (f - (1 - f) * g) # equation 5 in Ostojic 2014
        σ2 = τ_m * ν0 * C * J^2 * (f + (1 - f) * g^2) # equation 6 in Ostojic 2014
        σ = sqrt(σ2)/sqrt(2)
        #ε = 0.001 # stepsize to calculate slope of transfer function
        γ_μ = get_γ_μ(τ_m,σ,μ,Vth,Vre,ε,τ_ref,dV,Vl) # equation 13 in Ostojic 2014
        γ_σ2 = get_γ_σ2(τ_m,σ,μ,Vth,Vre,ε,τ_ref,dV,Vl) # equation 14 in Ostojic 2014
        println("γ_σ2 = ", round(γ_σ2,5), ", γ_μ = ", round(γ_μ,5) )
        Stabilitymatrix = τ_m / 1000 * (JMatrix * γ_μ + JMatrix2 * γ_σ2) # equation 12 in Ostojic 2014
        TimeStampBeforeTopo = time_ns() 
        d,v,nconv,niter,nmult,resid = eigs(Stabilitymatrix; nev = 1,which=:LR,ritzvec=true,maxiter=1300)
        println("calculating eigenvalue finished, it took ",round((time_ns()-TimeStampBeforeTopo)/1e9),"s")
        println(" nconv = ",(nconv)," niter = ",(niter)," nmult = ",(nmult)," norm of resid = ",norm(resid))
        # equation 16 in Ostojic 2014:
        λMax = τ_m / 1000 * J * sqrt(C) * sqrt( f * (γ_μ + J * γ_σ2)^2 + (1 - f) * g^2 * ( - γ_μ + g * J * γ_σ2)^2 )
        println(" λMaxMFT = ",λMax," abs(d) = ",abs(d)," d = ",(d))
        writedlm(path*"vReal.dat",real(v))
        writedlm(path*"d.dat",d)
        writedlm(path*"lambdaMax.dat",λMax)
        writedlm(path*"Arnoldi.dat",[nconv,niter,nmult,norm(resid)])
        vReal = real(v)
    else
        vReal = readdlm(path*"vReal.dat")
    end
    return vReal
end 
