function sturges(n)  # Sturges' formula
    n==0 && return one(n)
    iceil(log2(n))+1
end

function weightedhist!{HT}(h::AbstractArray{HT}, v::AbstractVector, edg::AbstractVector; init::Bool=true, weights::AbstractVector = ones(HT,length(v)))
# computes a weighted histogram, from stackoverflow unser Rory
    n = length(edg) - 1
    length(weights) == length(v) || error("length(weights) must equal length(v)")
    length(h) == n || error("length(h) must equal length(edg) - 1.")
    if init
        fill!(h, zero(HT))
    end
    for j=1:length(v)
        i = searchsortedfirst(edg, v[j])-1
        if 1 <= i <= n
            h[i] += weights[j]
        end
    end
    edg, h
end

weightedhist(v::AbstractVector, edg::AbstractVector; weights::AbstractVector = ones(Int,length(v))) = weightedhist!(Array(eltype(weights), length(edg)-1), v, edg; weights=weights)
weightedhist(v::AbstractVector, n::Integer; weights::AbstractVector = ones(Int,length(v))) = weightedhist(v, histrange(v,n); weights=weights)
weightedhist(v::AbstractVector; weights::AbstractVector = ones(Int,length(v))) = weightedhist(v, sturges(length(v)); weights=weights)
