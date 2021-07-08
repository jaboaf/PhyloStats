⊗(A::Array{T},B::Array{U}) where {T<: Union{String,Char},U<: Union{String,Char}} = map(x-> *(x...),Base.product(A,B))
