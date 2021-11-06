# C:\Users\joabo\Documents\GitHub\PhyloStats\Data\allnuc0712.tar\allnuc0712
# See https://docs.julialang.org/en/v1/base/io-network/#Base.eachline

# n is the number of genetic observations
n = 0
seq = ""
date = ""
loc = ""
newform = open("Data/GeneticObs.dat","w")
GenomeSeqs = Dict{String,Int32}()

fileNames = readdir(pwd()*"\\RawData")
for name in fileNames
	for line in eachline("RawData\\"*name)
	    if line[1] != '>'
	    	global seq = seq * line
	    else
			if seq in keys(GenomeSeqs) GenomeSeqs[seq] += 1
			else GenomeSeqs[seq] = 1
			end
			println(newform, date * "|" * loc * "|" * seq)
			global n = n+1
			global date = join(split(split(line,'|')[end],'-'),'|')
			global loc = split(line,'/')[2]
	    	global seq = ""
		end
	end
	println("done with $name")
end
close(newform)