# Exploratory Analysis of the GENETIC peices of all the genomic observations data GSAID had on 07/12/21


# The keys of this dict are "genes" and the values will be frequencies over sequences
Genes = ["NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP11","NSP12","NSP13","NSP14","NSP15","NSP16","Spike","NS3","E","M","NS6","NS7a","NS7b","NS8","N","NS9b","NS9c"]
GeneSeqs = Dict([ G => Dict{String,Int32}() for G in Genes])
GeneSeqs[""] = Dict{String,Int32}()

M = FiniteSuppFn(Dict())
for line in eachline("Data/GeneticObs.txt")
	M[split(line,"|",limit=11)[end]...] += 1
end


# A = Dict([G=>collect(values(GeneSeqs[G])) for G in keys(GeneSeqs)])

#=
L = ['A','C','G','T']
function M_L(s::String)
	N = countmap(s)
	return [ l in keys(N) ? N[l] : 0 for l in letters ]
end

uniquemonomials = map(uniqueparams) do 

struct FinitelySupportedFunction{It,Vt}
	D::Dict{It,Vt}
end

#
FiniteSuppFn(F::Dict{It,Vt}) = FiniteSuppFn(filter(p->p.last != zero(Vt),F))
# if i is in index set return D[i] else return 0
Base.getindex(F::FiniteSuppFn{It,Vt},i::It) where {It,Vt} = get(D, i, zero(Vt))
# if y is not 0 D[x] = y
Base.setindex!(F::FiniteSuppFn{It,Vt},y::Vt,x::It) where {It,Vt} = if !(iszero(y)) D[x] = y end
# If 
Base.getindex(F::FiniteSuppFn{It,Vt},I::Vector{It}) where {It,Vt} = map(i->get(D, i, zero(Vt)),I)

+(A::FiniteSuppFn{It,Vt},B::FiniteSuppFn{It,Vt}) where {It,Vt} = 
=#


# uniqueparams is ultimately H_n
# M_geneseqs is counts of uniqueparams
# countvals is counts of integer vals of M_geneseqs