GENES = Set()
YYYY = Set()
MM = Set()
DD = Set()
TYPE = Set()
SUBLOC = Set()
LOC = Set()
for line in eachline("Data/allnuc0801.fasta")
	if line[1] == '>'
		gene,yyyy_mm_dd,passage,type_subloc,host,loc = split(line[2:end],'|')[[1,3,5,6,7,end]]
		if passage == "Original"
			push!(GENES,gene)
			yyyy, mm, dd = split(yyyy_mm_dd,'-')
			push!(YYYY,yyyy)
			push!(MM,mm)
			push!(DD,dd)
			type, subloc = split(type_subloc,"^^")
			push!(TYPE,type)
			push!(SUBLOC,subloc)
			push!(LOC,loc)
		end
	end
end
println("Length of GENES: $(length(GENES)) ")
println("Length of YYYY: $(length(YYYY)) ")
println("Length of MM: $(length(MM)) ")
println("Length of DD: $(length(DD)) ")
println("Length of TYPE: $(length(TYPE)) ")
println("Length of SUBLOC: $(length(SUBLOC)) ")
println("Length of LOC: $(length(LOC)) ")

			


