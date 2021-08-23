Genes = ["NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP11","NSP12","NSP13","NSP14","NSP15","NSP16","Spike","NS3","E","M","NS6","NS7a","NS7b","NS8","N","NS9b","NS9c"]
GeneSeqs = Dict([ G => Dict{String,Int32}() for G in Genes])
GeneSeqs[""] = Dict{String,Int32}()

gene = ""
seq = ""
for line in eachline("Data/allnuc0801.fasta")
	if line[1] == '>'
		if seq in keys(GeneSeqs[gene]) GeneSeqs[gene][seq] += 1
		else GeneSeqs[gene][seq] = 1
		end
		global gene = split(line[2:end],'|')[1]
		global seq = ""
	else
		global seq = seq * line
	end
end

function Newt(N::Int32,K::Int32,tol::Float64)
	F(x) = K/x - sum([1/(x+j-1) for j in 1:N]) 
	f(x) = -K/x^2 + sum([1/(x+j-1)^2 for j in 1:N])
	x = 1
	while abs(F(x))>=tol
		x = x - F(x)/f(x)
	end
	return x
end

println("Gene & N & K & Newt(N,K,0.0001)")
for g in Genes 
	N = sum(values(GeneSeqs[g]))
	K = length(GeneSeqs[g])
	println(g," & ",N," & ",K," & ",Newt(N,K,0.0001))
end

# for g in Genes println(Newt(sum(values(GeneSeqs[g])),length(GeneSeqs[g]),0.0001))

# using GRUtils

Nmax = maximum([sum(values(GeneSeqs[g])) for g in Genes])
Kmax = maximum([length(GeneSeqs[g]) for g in Genes])
Vmax = maximum([ maximum(values(GeneSeqs[g])) for g in Genes])

#=
GRUtils.plot(sort(collect(values(GeneSeqs[Genes[1]]))),label=Genes[1],xlim=(0,Kmax),ylim=(0,Vmax))
for g in Genes[2:end] GRUtils.oplot(sort(collect(values(GeneSeqs[g]))),label=g,xlim=(0,Kmax),ylim=(0,Vmax)) end

GRUtils.plot(sort(collect(values(GeneSeqs[Genes[1]]))),label=Genes[1],xlim=(0,Kmax),ylim=(0,Vmax))
for g in Genes[2:end] GRUtils.oplot(sort(collect(values(GeneSeqs[g]))),label=g,xlim=(0,Kmax),ylim=(0,Vmax)) end

GRUtils.plot(sort(collect(values(GeneSeqs[Genes[1]]))),label=Genes[1],xlim=(0,Kmax),ylim=(0,Vmax))
for g in Genes[2:end] GRUtils.oplot(sort(collect(values(GeneSeqs[g]))),label=g,xlim=(0,Kmax),ylim=(0,Vmax)) end
xlabel("Index of Sequence for Gene, sorted by frequency")
ylabel("Frequency")


GRUtils.plot(LinRange(0,1,length(GeneSeqs[Genes[1]])),cumsum(sort(collect(values(GeneSeqs[Genes[1]])))),label=Genes[1],ylim=(0,Vmax))
for g in Genes[2:end] GRUtils.oplot(LinRange(0,1,length(GeneSeqs[g])),cumsum(sort(collect(values(GeneSeqs[g])))),ylim=(0,Vmax)) end
xlabel("Index of Order Statistics Embedded in [0,1]")
ylabel("Cumulative Frequency")
title("Empirical Cumulative Densities of Order Statistics By Gene")


GRUtils.plot(LinRange(0,1,length(GeneSeqs[Genes[1]])),cumsum(sort(collect(values(GeneSeqs[Genes[1]])))) / sum(values(GeneSeqs[Genes[1]])),label=Genes[1],ylim=(0,1))
for g in Genes[2:end] GRUtils.oplot(LinRange(0,1,length(GeneSeqs[g])),cumsum(sort(collect(values(GeneSeqs[g])))) / sum(values(GeneSeqs[g])),ylim=(0,1)) end
xlabel("Index of Order Statistics Embedded in [0,1]")
ylabel("Probability")
title("Empirical Cumulative DF of Order Statistics By Gene")
=#