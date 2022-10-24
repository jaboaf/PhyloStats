# should I make these ⨂ ? it wont conflict with actual product operation
⊗(A::Array{T},B::Array{U}) where {T<: Union{String,Char},U<: Union{String,Char}} = map(x-> *(x...),Base.product(A,B))
⊗(A::Array{T},B::Array{T}) where T<: Number = prod.(Base.product(A,B))
⊗(V::Vararg{Array{T}}) where T<: Number =reduce(⊗,V) 


⊗(a::NTuple{T},b::NTuple{T}) where T  = (a...,b...)
⊗(a::NTuple{N,T},b::NTuple{M,T}) where {T,N,M}  = (a...,b...)

# Why duplicate?
⊕(V::Vararg{Array{T,N} where N}) where T = [ sum(filter(x->size(x)==d,V)) for d in unique(size.(V)) ] 
⊕(V::Vararg{Array{T,N}}) where {T,N} = [ sum(filter(x->size(x)==d,V)) for d in unique(size.(V)) ] 


# maybe
⨁(A::Array{Array{T,N} where N,1}) where T <: Number = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨁(A::Array{Array{T},1}) where T<:Number = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨁(A::Array{Array{T,N},1}) where {T<:Number,N} = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨂(A::Array{Array{T,N} where N,1}) where T <: Number = reduce(⊗,A)

