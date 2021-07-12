using StatsBase
using GRUtils
#=
C_L(s::String) = [count(==(l),s) for l in ['A','C','G','T']];
function C_N(itr)
    M = maximum(itr)
    cnts = zeros(Int,M)
   	for v in itr cnts[v] += 1 end
    return cnts
end

heC_N(values(countmap(C_L.(SEQs))))
C_N(values(countmap(sort.(C_L.(SEQs)))))
=#
# A large amount of Genomic and Genetic data is widely available. An observation typically consists of a sequence of nucleotide bases, and any information associated with it; for example, specie, time, or place. We call this a genomic observation if the sequence is a full length DNA sequence, and a genetic one if the sequence corresponds to a gene identified by a sub-sequence of DNA. The structure of sequence data is quite simple: pick any number of 'A's,'C's,'G's, and 'T's, and arrange them in any way. Frequently, sequences contain the letter 'N', denoting that the sequencing machines was unable to determine the nucleotide at a position. These machines are also subject to incorrectly reading a nucleotide.

# Here are some references:
# Pfeiffer, F., Gröber, C., Blank, M. et al. Systematic evaluation of error rates and causes in short samples in next-generation sequencing. Sci Rep 8, 10950 (2018). https://doi.org/10.1038/s41598-018-29325-6
# Ma, X., Shao, Y., Tian, L. et al. Analysis of error profiles in deep next-generation sequencing data. Genome Biol 20, 50 (2019). https://doi.org/10.1186/s13059-019-1659-6

# Statistical methods should be able to accomodate sequence data alone, and more importantly, genomic or genetic observations. These must be computationally feasible in context of the quantity of data available. Otherwise researchers are forced to use less data than what is available, which leaves the world less informed. Perhaps we ought to be more informed.

#( Its oddly hard to find the number of observations in a bunch of the covid-19 research, and many of the papers sample from their data).

# \section{Data}
# We've been data-less, so for concreteness I've been playing with covid data.  This was retrieved from, and publicly available at the National Center for Biotechnology Innovation (NCBI).

# I have a hard copy of all of the "complete" sequences as of ~ 2 weeks ago on my pc (there are ~ M = 356,471 of these), and i put a "toy" sample of start of covid to 5/31/20 in the github to try different stuff out with. Generally the following variables are available.
# \begin{itemize}\item Accesion: This seems to uniquely identify a sequence
# \item Submitters: This is a list of names of the people attached to the submission of the sequence.
# \item Release Date: This is the date the submission was released I think.
# \item Pangolin: unsure.
# \item Species: This takes the value "Severe acute respiratory syndrome-related coronavirus"
# \item Molecule Type: This takes the values "ssRNA(+)"
# \item Length: Some Integer. I am wondering why these are not all multiples of 3.
# \item Geo Location: Of the form "Country:City"
# \item USA: If given, this is the state of the sample.
# \item Host: This takes the values "Homo sapiens", "Felis Catus", and "Canis Lupus Familiaris"
# \item Isolation Source: Where/how the sample was collected from the source, e.g. swab.
# \item Collection Date: Date of collection, in any of the following formats, YYYY, YYYY-MM, or YYYY-MM-DD
#\item Sequence: This is some word on the alphabet: 'A','C','G','T','N'. 'N' denotes "unknown nucleic acid residue". There seem to be quite a bit of 'N's. Approximately, 10\% of the first sequence is 'N'.
#\end{itemize}

# In the toy example, I have chosen the following variables based on my interests: collection date, location, sequence. I also used length, but this isn't a variable, it is statistic (right?).

# \section{Information Space}
# I'd like to break this section up into a couple parts.
# 1) Preliminary Definitions and Notations:
# 2) Sequence Space: This accomodates any observed sequence with no 'N's.
# 3) Ambigious Sequence Space: This accomodates any observed sequence.
# 4) Observation Space: This combines ambiious sequence space with other information.
# 5) Data Space: I don't know yet, this may be better as its own section.

#\section{Preliminary Definitions and Notations}
B = Set(['A','C','G','T'])
B̄ = ['A';'C';'G';'T'] # mathematically this is equal to a 4-tuple
toydata = read("some_sequences.fasta",String);
n = count(==('>'),toydata)
yyyys = String[]; mms=String[]; dds=String[];
geos=String[]; lens= Int[]; seqs = String[];
SAMPLE = Tuple{String,String,String,String,Int64,String}[]
for x in split(toydata,'>')[2:end]
	x = replace(x,'\n'=>"")
	x = split(x,'|')
	# date,location,length,seq

	geo = String(x[3]); push!(geos,geo)
	len = parse(Int,x[2]); push!(lens,len)
	seq = String(x[5]); push!(seqs,seq)

	k = length(split(x[1],"-"))
	if k ==0 yyyy = ""; mm = ""; dd = ""
	elseif k==1 yyyy=split(x[1],"-")[1]; mm= ""; dd = ""
	elseif k==2 yyyy,mm = String.(split(x[1],"-")); dd = ""
	else yyyy,mm,dd = String.(split(x[1],"-"))
	end
	push!(yyyys, String(yyyy)); push!(mms, String(mm)); push!(dds, String(dd));

	push!(SAMPLE, (yyyy,mm,dd,geo,len,seq))
end; # note that the first element of split(toydata,'>') is ""
sort!(SAMPLE);
YYYY = sort(unique(yyyys)); MM = sort(unique(mms)); DD = sort(unique(dds));
GEO=sort(unique(geos)); LEN =sort(unique(lens)); SEQ = sort(unique(seqs));
println("These are $(length(GEO)) unique geolocations ")
println("These are $(length(LEN)) unique lengths of sequences ")
println("These are $(length(SEQ)) unique sequences")

Nice = filter(s->!('N' in s),SEQ);

# \section{Sequence Space}

# Lets start by identifying the measurable space that each (completely specified) sequence is contained in because then our sample will be in the n-fold product of this space (where n in the number of observations in our sample).
# So, we are in need of space that accomodates all of our observations. There are a couple of weird things that come up in biological data that we should (be able to) deal with such as 'N's in sequences denoting that the nucleotide in a position was not able to be observed. Other things include alignment (though I am not sure if this is a scientific procedure, under the biologists juristiction, or something of a statistical nature). Regardless, these are not a result of the data-generating process; they are a result of a (possibly informed) observation process. Maybe lets focus on the "point data", and then return to some of these intricacies? (where "point data" is data of points, and a point is a complete sequence with no 'N's)

# It appears that the length of a sequence and the counts of nucleotides in a sequence vary across the sample. So our sequence space should accomodate sequences with different lengths and different counts of nucleotides. 

# Without further chattiness, denote the set of nucleotide bases by B. i.e. B = \{'A','C','G','T'\}.
# That means the set of sequences of length l is Bˡ = \{ (w¹,…,wˡ) ∣ 1≤i≤m, wⁱ ∈ B \}.
# So, the set of all sequences is the union over l ∈ ℕ of sequences of length l i.e. ⋃_{l∈ℕ} Bˡ. We will denote this set B^⋆.

# Bada bing! Any sequence of 'A','C','G','T' is an element of B^⋆, so it seems like B^⋆ is a good candidate for the set of possible outcomes for an observation. We can take 𝒫(B^⋆) to be the sigma-algebra. We could endow this space with a measure, and if we were being honest bayesians we'd do this before looking at the data.
# So, what probability measure would you like?
# It is easy to pick ~a~ probability measure on this discrete space, however I don't enjoy doing this because I know very very little about covid so I'm just making arbitrary choices. Also, I find it way more fun to explore data, than to arbitrarily specify some probability measure that ought to reflect how shit actually works in the world.

# So maybe we could take a more adventurous approach that allows us to perform inference (whatever that means) and estimate probabilities in a consistent manner (in the precise logical definition). Any probability measure that reflects the data generating process seems to satisfy these desires.
# Recall that we are still in the world of "point data", 𝐱 = (x_1,…,x_n) where x_i ∈ B^⋆
# One such measure, is the empirical measure: ℙₙ = 1/n ∑_i δ_{x_i}
# ℙₙ:𝒫(B^⋆)→[0,1]
# In the "non-point data" world we will have to be more careful because ∫_{B^⋆}δ_{X}(w) dw ≥ 1 for X ∈ 𝒫(B^⋆)\∅. We will return to this.

# We have not explitly defined a parameter space or a parametric family of functions to define the empirical measure ℙₙ, however we can do this using "parameters". Below are some ways of doing this:
# 0.0) Let Ω = B^⋆. A density (or function) on B^⋆ is f:B^⋆→ℝ, so f = ∑_{w ∈ B^⋆} f(w)1_w = ∑_{w ∈ B^⋆} f_w 1_w. Hence, the collection of densities is a "parametric family" with the parameter (f_w)_{w∈B^⋆}. NOTE: Technially, for (f_w)_{w∈B^⋆} to be well defined we must (be able to) order B^⋆; one can do this by partially ordering B^⋆ by length and then lexicographically w.r.t. B̄. The order of application of these partial orders gives a total order, <, and (B^⋆,<) looks like: ()<(A)<(C)<(G)<(T)<(A,A)<(A,C)<…<(T,G)<(T,T)<(A,A,A)<(A,A,C)<…<(T,T,G)<(T,T,T)…. One could also partially order B^⋆ lexicographically w.r.t. B̄ and then by length, so (B^⋆,<') would look like: ()<(A)<(A,A)<(A,A,A)<…<(A,T)<…<(B)<(B,B)
# 0.1) A probability density on B^⋆ is a density f:B^⋆→ℝ such that ∑_w |f_w| = 1. Hence, we may have "parametric family" of probability densities with the parameter (f_w)_{w∈B^⋆}.
# 0.2) A probability distribution on B^⋆ is a density f:B^⋆→ℝ taking non-negative values, i.e. f_w ≥ 0, such that ∑_w f_w = 1. Hence, we may have "parametric family" of probability distributions with the parameter (f_w)_{w∈B^⋆}.
# 1) Let Θ = ℕ⁴ be the parameter space. Each θ ∈ Θ determines an equivalence class of B^⋆ where there are θ^1=θ_A 'A's, θ^2=θ_C 'C's,θ^3=θ_G 'G's, and θ^4=θ_T 'T's. The equivalence class is the set of all rearrangements of (A)^{⊗θ_A}⊗(C)^{⊗θ_C}⊗(G)^{⊗θ_G}⊗(T)^{⊗θ_T} = B̄^{θ_{B̄}} is 𝔖^{|θ|} B̄^{θ_{B̄}}. And Θ partitions B^⋆ because for any w ∈ B^⋆, N_{B̄}w ∈ ℕ^4 means that N_{B̄}w = ϑ for some ϑ∈Θ.
# Note: You can replace (b)^{⊗θ_b} with {b}^{×θ_b} to get an e}; thaquivalent expression. The latter is a more set theorhetic approcach. I will avoid this for a very particular reason that I hope to get to later.
# 2) Let Φ = \{ ϕ∈ℕ⁴ ∣ 1≤i≤j≤4 ⟹ϕᵢ≤ϕⱼ\}. Each ϕ ∈ Φ determines a equivalence class which is a union of equivalence classes from 1. Namely, the set of w in B^⋆ such that w has ϕ_1 b₁s, ϕ_2 b₂s, ϕ_3 b₃s, and ϕ_4 b₄s, and {b₁,b₂,b₃,b₄}=B. Namely, ϕ∈Φ determines the equivalence class, B_ ⋃_{g∈𝔖^4} 𝔖^{|ϕ|} B̄^{g ϕ_{B̄}}. This partitions B^⋆ because it is a coarsening of the partiion given in 1. Also, ⋃_{ϕ∈Φ} 𝔖ϕ = ℕ⁴.
# 3) Let ψ = (ψ^{(l)})_{l∈ℕ} where ψ^{(l)} ∈B^l. Then for any w ∈ B^l we may find an element of 𝔖_B^l, such that gψ^{(l)} = w. There may be multiple such g, so we can actually find a subgroup H of 𝔖_B^l, such that Hψ^{(l)} = \{w\}.
# 4) This one is a riff on 3) and the finite length sequence approach to everything. First, let Ψ ∈B^∞ be some fixed element. Now, for the rest of this bullet, let w ∈ B^l. Instead of thinking of w as some finite sequence, we could think of it as an equivlence class of sequences w̃ ⊂ B^∞, where every v ∈ w̃ has its first \# w elements equal to w, so w̃  is an element of B^∞ / ⟨e^{(1,…,\# v)} v' = v ∣ v∈B^⋆⟩=\tilde{B^⋆}. Formally, ̃ is a function, ̃:B^⋆→𝒫(B^∞), which is defined by w↦\{v∈B^∞∣ (v^1,…,v^{\# w}) = w \}). Since, B^∞ = 𝔖_B^∞ Ψ, there is a subgroup H = (H_1,H_2,…) ∈ 𝔖_B^∞ such that HΨ = w̃. The structure of H is fairly simple: h ∈ H_1 maps Ψ^1 to w^1,…,h ∈ H_{\# w} maps Ψ^{\# w} to w^{\# w},h∈H_{\# w + 1} maps Ψ^{\# w +1} to any element of B,…. So for 1≤i≤\# w, Hᵢ = (Ψ^i \: w^i)𝔖_{ B \\ \{Ψ^i,w^i\}}I_B, and for i>\#w, Hᵢ = 𝔖_B.

# Some topic-specific thoughts on each parameterization before we get concrete.
# 0) This basic probability on the combinitorial structure of sequences. I would avoid using this as THE approach. This considers the dual space of B^⋆, aka, the space of linear functionals from B^⋆ into a field of choice. After looking at it this way, the dual space gives a better perspective on measures, and is the same as (V_B^⋆)^* when \{e_A,e_C,e_G,e_T\} is an orthonormal basis for V_B. If \{E_A,E_C,E_G,E_T\} is orthogonal (so ||E_b|| need not be 1), then we can consider densities over the basis E_w = ||E_w|| e_w; these elements can get very big or very small, very slowly or very quickly.
# 1) This one is probably the most intuitive and very useful. For any sequence its not so hard to count the number of 'A's,'C's,'G's, and 'T's.
function M_B̄(s::String)
	cnts = zeros(Int,4)
	for c in s
		cnts += B̄ .== c
	end
	return cnts
end

𝔸 = M_B̄.(seqs);
sum.(𝔸)
# note !! n on x axis,
plot(sort(sum.(𝔸)))
#e.g. plot(sort(sum.(𝔸))[1:10000])
maximum(𝔸)

#=
Polya enumeration.
Y^X is set of functions from X to Y
|𝔖||Y^X/𝔖| = ∑_{g∈𝔖} |Y|^{c(g)}

Burnsides lemma
X^g = {x∈X|gx=x}
|𝔖||X/𝔖| = ∑_{g∈𝔖} |X^g|

𝔖_x = {g∈𝔖| gx=x}
|𝔖x| = [𝔖:𝔖_x] = |𝔖|/|𝔖_x| 
=#
# % Statiticians and probabalists sometimes go about this by trying out different models, which are specified by probability measures, and determining one that best fits the data.


#lets just consider our data to be stochastic process,
# So the first observation x₁
#=
π:Θ→[0,1] is a prior distribution

π(θ|x) = P(x|θ)π(θ) / m(x)

P(x and θ) = P(x*1_{𝔖θ}) = (x*1_{𝔖θ})^*|_1
P(x|θ) =  = x*1_{𝔖θ}/π(θ)

P(x|θ) = ∂_θ(x^*)|_0
P(θ|x) = ∂_x(θ^*)|_0

m(x) = ∫_Θ P(x|θ) π(θ) dθ = ∫_Θ P(x|θ) dπ = ∫_Θ P(x|θ) dπ
P(θ) = ∫_Ω P(θ|ω) P(ω) dω = ∫_Ω P(θ|ω) dω =

NOTE: P(x|θ) = P_θ(x)

P(𝐗|θ) = ∂_θ  𝐗^* |_0


P(θ|𝐗) = P(θ|⨁_{l∈ℕ} 𝐗_l) = ⨁_l P(θ|𝐗_l)
P(θ|𝐗) = P(θ|⨁_{θ∈Θ} 𝐗*1_{𝔖θ}) = ⨁_{θ∈Θ} P(θ|𝐗*1_{𝔖θ})


P(θ|𝐗) = P(𝐗|θ)P(θ) / P(𝐗)
P(𝐗) = ∫_Θ P(𝐗|θ)P(θ) dθ
=#
# I am going to be a bad bayesian. Let P = ℙ and π = ℼ. ℼ = 1/n∑_i δ_{N_B̄(X_i)}
#(θ|x) = P(x|θ)π(θ) / m(x)


# eⁱ(n_1e_1+…+n_ke_l) = n_i
# for some orthogonal (?) e_b \in V_B and E = (e_A,e_C,e_G,e_T), the dual of E is \hat{E} = (e^1,e^2,e^3,e^4).

# EᵀÊ = I_l = J_l*1_l
# If we let E_w = (e_{w¹},…,e_{wˡ}); 
# maybe let ⊗E_w = e_{(w¹,…,wˡ)} and ⊗Ê^ℕ = e^{(1,…,l)}
# ⊗Eᵀ⊗Ê = λ_{≤l} = ∑_i e_{w^i}
# e.g. E_w = e_

# A bit of a radical take on this is to condider each observation to be an equivalence class of sequences
# 𝔖 B^* / 𝔖_{1^l,∞} = 𝔖^l𝔖^B B̄

# Orr just say we have the following : 𝔾^l(𝔾_B θ_1,…,𝔾_B θ_k) = 𝔾^l 𝔾_B^{×k} θ

# orr maybe this is just the entire space: 𝔾^l (𝔾_{B,1}θ_1,…,𝔾_{B,k}θ_k)

# wait no, it is just (𝔾_{B,1}θ_1,…,𝔾_{B,l}θ_l) ∀ l∈ℕ and for some defined θ ∈ B^{∞}. (it should be θ = B̄). maybe we ought to say 𝔾_{B}^* θ_*
# I really like this one

# Another way i like to paramaterize this is ⨁_{α∈ℕ^k} 𝔾^{|α|} θ^{⊗α} = ⨁_l ⨁_{π ∈ Π_l^k} 𝔾^{|π|} θ^{⊗π}
# this has been my favorite for a while, its transparent


# which we could parameterize as ⨁_{α∈ℕ^k s.t. i<j⟹αᵢ<αⱼ} 𝔾^{|α|} (𝔾_B)^{⊗|α|} θ^{⊗α}, where (𝔾_B)^{⊗|α|} is meant as the \{(g,\dots,g)|g ∈ 𝔾_B \}
# the above is equal to  ⨁_{α∈ℕ^k s.t. i<j⟹αᵢ<αⱼ} 𝔾^{|α|} θ^{⊗𝔾^kα} = ⨁_l ⨁_{y∈Υ_l^k} 𝔾^{l} θ^{⊗𝔾^ky}
# 




# \section{Kingman's Approach}
# Kingman has an interesting article, "Random Partitions in Population Genetics", in which he defines a partition structure, describes its connection with Ewens' sampling formula, and provides some probablistic tools. In this paper he views a partition as an integer partition of the frequency distribution of the sample:
# Let M_{B^⋆}: CtsTime → ℕ^{B^⋆} be a counting process (using the ordering described in Sequence Space.0), and let 𝕄_{B^⋆} be the empirical analouge. Now define another counting process M_N:ℕ^{B^⋆} → ϖ_N, where ϖ_N = \{a∈ℕ^N ∣∑_{j=1}^N ja_j = N\} is the set of ordered integer partitions of N. Kingman's random partition is π = M_{|𝕄_{B^⋆}(\text{right now})|}∘𝕄_{B^⋆}(\text{right now}). If π=(a_1,a_2,…,a_n), then a_1 sequences appeared once, a_2 sequences appeared twice, ..., a_n sequences appeared n times.
# Kingman then defines the set of all probability distributions on ϖ_n as Π_n, so any P ∈ Π_n satisfies P(π)≥0 and ∑_{π∈ϖ_n} P(π) = 1. So you might wonder what the probability of π=(a_1,a_2,…,a_n) is. It could be anything. If you also specified that for every m<n (these are variables for sample size), P_m∈Π_m must be the distribution arising from sub sampling without replacement from P_n ∈ Π_n, you can consider the probability of π :
# "Consider an infinite population divided into a countable number of sub-populations labelled (say) by colours 1,2,3,..., and suppose that the proportion of the population which is of colour i is x_i, where x_i ≥ 0, ∑_{i=1}^∞ x_i = 1" (p.4). Then P_n(π) = ϕ_π(y) = n!∏_{j=1}^n (j!)^{-a_j} ∑_{v∈[0:n]^∞ ∣ M_n(v) = π \} x_1^{v_1}x_2^{v_2}… . For the particular choice of y=x, this is a number, otherwise this is a function of the variable y.

# I don't want to get too deep into this because I don't think this is the right perspective. It seems that the motivation for this view of a partition (and the resulting importance of a partition structure) comes from considering "a geneticist that examines a sample of n gametes from a population of size N, and whose experimental techniques do not permit a labelling of the alleles he can distinguish"(p.3). I challenge you to distinguish things without being able to label them.

hashSEQ = Dict([s=>i for (i,s) in enumerate(SEQ)]);
countingSEQprocess = zeros(Int,(n,length(SEQ)));
for (i,s) in enumerate(seqs)
	countingSEQprocess[i,hashSEQ[s]] += 1
end

# Letters or fancy letters used thus far:
# 'A', B, 'C', 'G', 'T', 'N',
# N for random n,
# M for counting
# ^⋆,𝒫,ℕ,∫Φ,ϕ