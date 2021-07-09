using StatsBase
@enum NucleotideBase A=1 C G T

sf = read("Data/sequences.fasta",String);
println("read in data")
samp = split(sf,'>')[2:end];
println("separated obs")
samp = map(d->split(d,'\n',limit=2),samp);

subsamp = samp[1:10000];
println("sub sampled first 10000 obs")
INFOs = first.(subsamp);
SEQs = map(s->replace(s,'\n'=>""),last.(subsamp));

C_L(s::String) = [count(==(l),s) for l in ['A','C','G','T']];
function C_N(itr)
    M = maximum(itr)
    cnts = zeros(Int,M)
   	for v in itr cnts[v] += 1 end
    return cnts
end

heC_N(values(countmap(C_L.(SEQs))))
C_N(values(countmap(sort.(C_L.(SEQs)))))


