function makeConnectivity(N,C,f,g,J,seed)
# generates a network of N neurons with fixed indegree C
    Ni = round(Int,N - N*f)
    Ne = N - Ni
    Ce = round(Int,f*C)
    Ci = round(Int,C - C*f)
    srand(seed)
    JMatrix = zeros(N,N)
    for n = 1:N # function sample is from module Distributions
        for m = sample([1:min(n-1,Ne);(n+1):Ne],Ce,replace = false) 
            JMatrix[n,m] = J
        end
        for m = sample([Ne+1:n-1;max(n+1,Ne+1):N],Ci,replace = false)
            JMatrix[n,m] = -g*J
        end
    end
    return JMatrix
end
