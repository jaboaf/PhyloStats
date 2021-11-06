import Base.+
import Base.*

struct FiniteSuppFn{It,Vt}
	D::Dict{It,Vt}
end

#FiniteSuppFn(D::Dict{It,Vt}) where {It,Vt}= FiniteSuppFn(filter(p->p.last != zero(Vt),))
# if i is in index set return D[i] else return 0
Base.getindex(F::FiniteSuppFn{It,Vt},i::It) where {It,Vt} = if i in keys(F.D) F.D[i] else zero(Vt) end
# if y is not 0 then let D[x] = y
Base.setindex!(F::FiniteSuppFn{It,Vt},y::Vt,x::It) where {It,Vt} = if !(iszero(y)) F.D[x] = y end
# If 

Base.getindex(F::FiniteSuppFn{It,Vt},I::Vector{It}) where {It,Vt} = map(F[i],I)
Base.getindex(F::FiniteSuppFn{It,Vt},I::Set{It}) where {It,Vt} = map(i->get(F.D, i, zero(Vt)),I)

function +(A::FiniteSuppFn{It,Vt},B::FiniteSuppFn{It,Vt}) where {It,Vt}
	if A.D.count >= B.D.count
		+(B,A)
	else
		for k in union(A.D.keys,B.D.keys) A.D[k] += B[k] end
	end
end
function *(A::FiniteSuppFn{It,Vt},B::FiniteSuppFn{It,Vt}) where {It,Vt}
	if A.D.count >= B.D.count
		*(B,A)
	else
		# I should remove A.D.keys\B.D.keys from A.D
		for k in union(A.D.keys,B.D.keys) A.D.keys[k] *= B[k] end
	end
end

# inv(A::FiniteSuppFn{It,Vt}) where {It,Vt} = Dict([])


# lineseqs = mapreduce(l->LToPowerB[l],(x,y)->prod.(Base.product(x,y)), line)
# global seqs = prod.(Base.product(seqs,lineseqs))

#=
k = length(seqs)
for s in seqs 
	try
		GeneSeqs[s] += 1/k
	catch
		GeneSeqs[s] = 1/k
	end
end
global seqs = [""]
=#
