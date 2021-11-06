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

# Statistical methods should be able to accomodate sequence data alone, and more importantly, genomic or genetic observations. These must be computationally feasible in context of the quantity of data available. Otherwise researchers are forced to use less data than what is available, which leaves the world less informed. This is not good.

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
data = read("some_sequences.fasta",String);
n = count(==('>'),data)
yyyys = String[]; mms=String[]; dds=String[];
ctrys=String[]; lens= Int[]; seqs = String[];
SAMPLE = Tuple{String,String,String,String,Int64,String}[]
for x in split(data,'>')[2:end]
	x = split(strip(x),'|') # gives date,location,length,seq
	
	ctry = String(split(x[3],':')); push!(ctrys,ctry)
	len = parse(Int,x[2]); push!(lens,len)
	seq = String(x[5]); push!(seqs,seq)

	k = length(split(x[1],"-"))
	if k ==0 yyyy = ""; mm = ""; dd = ""
	elseif k==1 yyyy=split(x[1],"-")[1]; mm= ""; dd = ""
	elseif k==2 yyyy,mm = String.(split(x[1],"-")); dd = ""
	else yyyy,mm,dd = String.(split(x[1],"-"))
	end
	push!(yyyys, String(yyyy)); push!(mms, String(mm)); push!(dds, String(dd));

	push!(SAMPLE, (yyyy,mm,dd,ctry,len,seq))
end; # note that the first element of split(data,'>') is ""
sort!(SAMPLE);
YYYY = sort(unique(yyyys)); MM = sort(unique(mms)); DD = sort(unique(dds));
CTRY=sort(unique(ctrys)); LEN =sort(unique(lens)); SEQ = sort(unique(seqs));
println("These are $(length(CTRY)) unique countries ")
println("These are $(length(LEN)) unique lengths of sequences ")
println("These are $(length(SEQ)) unique sequences")

# \section{Sequence Space}

# Lets start by identifying the measurable space that each (completely specified) sequence is contained in because then our sample will be in the n-fold product of this space (where n in the number of observations in our sample).
# So, we are in need of space that accomodates all of our observations. There are a couple of weird things that come up in biological data that we should (be able to) deal with such as 'N's in sequences denoting that the nucleotide in a position was not able to be observed. Other things include alignment (though I am not sure if this is a scientific procedure, under the biologists juristiction, or something of a statistical nature). These are not reflective of the biological information, rather, the (erronous or possibly informed)data-generating process. I am not a biologist or geneticist, so I will ignore the question of alignment until later. On the other hand, I do feel well equipped to handle 'N's because they are simple: each instance of 'N' is a placeholder for a nucleotide base. Also, it seems important to include them in our understanding of the data because they show up quite a bit.

numberOfSequencesWithNs = count(s ->'N' in s, seqs)
println(" $numberOfSequencesWithNs of the $n sequences contain at least one N, which is $(numberOfSequencesWithNs/n) percent of the sequences")

# Upon looking into the Ns I found out that there are other letters that occur in the sequences: 
print(union(unique.(seqs)...))

M_N(s::String) = count(==('N'),s)

numNsInseqs = M_N.(seqs);

#= https://www.qmul.ac.uk/sbcs/iubmb/misc/naseq.html
3.1. Guanine. adenine, thymine, cytosine: G,A,T,C
3.2. Purine (adenine or guanine): R
3.3. Pyrimidine (thymine or cytosine): Y
3.4. Adenine or thymine: W
3.5. Guanine or cytosine: S
3.6. Adenine or cytosine: M
3.7. Guanine or thymine: K
3.8. Adenine or thymine or cytosine: H
3.9. Guanine or cytosine or thymine: B
3.10. Guanine or adenine or cytosine: V
3.11. Guanine or adenine or thymine: D
3.12. Guanine or adenine or thymine or cytosine: N
so we need a map
'A'=>'A'
'B'=>['C','G','T']
'C'=>'C'
'D'=>['A','G','T']
'E'=>???
'F'=>???
'G'=>'G'
'H'=>['A','C','T']
'I'=>???
'J'=>???
'K'=>['G','T']
'L'=>???
'M'=>['A','C']
'N'=>['A','C','G','T']
'O'=>???
'P'=>???
'R'=>['A','G']
'S'=>['G','S']
'T'=>'T'
'U'=>???
'V'=>['A','C','G']
'W'=>['A','T']
'X'=>???
'Y'=>['C','T']
'Z'=>???
=#


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
# ℙₙ:B^⋆→[0,1]
# In the "non-point data" world we will have to be more careful because ∫_{B^⋆}δ_{X}(w) dw ≥ 1 for X ∈ 𝒫(B^⋆)\∅. We will return to this.

# We have not explitly defined a parameter space or a parametric family of functions to define the empirical measure ℙₙ, however we can do this using "parameters". Below are some ways of doing this:
# 0.0) Let $ Ω = B^⋆ $. A density (or function) on $ B^⋆ $ is $ f:B^⋆→ℝ $, so $ f = ∑_{w ∈ B^⋆} f(w)1_w = ∑_{w ∈ B^⋆} f_w 1_w $. Hence, the collection of densities is a "parametric family" with the parameter (f_w)_{w ∈ B^⋆}. NOTE: Technially, for $ (f_w)_{w∈ B^⋆} $ to be well defined we must (be able to) order $ B^⋆ $; one can do this by partially ordering $ B^⋆ $ by length and then lexicographically w.r.t. $ \bar{B}$. The order of application of these partial orders gives a total order, $<$, and $(B^⋆,<)$ looks like: $ ()<(A)<(C)<(G)<(T)<(A,A)<(A,C)<…<(T,G)<(T,T)<(A,A,A)<(A,A,C)<…<(T,T,G)<(T,T,T)… $. One could also partially order $ B^⋆ $ lexicographically w.r.t. \bar{B} and then by length, so $(B^⋆,<') $ would look like: $ ()<'(A)<'(A,A)<'(A,A,A)<'… <'(A,T)<'… <'(B)<(B,B) $
# 0.1) A probability density on $ B^⋆ $ is a density $ f:B^⋆→ℝ $ such that $ ∑_w |f_w| = 1 $. Hence, we may have "parametric family" of probability densities with the parameter $(f_w)_{w∈B^⋆}$.
# 0.2) A probability distribution on $ B^⋆ $ is a density $ f:B^⋆→ℝ $ taking non-negative values, i.e. $ f_w ≥ 0 $, such that $ ∑_w f_w = 1 $. Hence, we may have "parametric family" of probability distributions with the parameter $(f_w)_{w∈B^⋆}$.
# 1) Let $Θ = ℕ⁴$ be the parameter space. Each $ θ ∈ Θ $ determines an equivalence class of $ B^⋆ $ where there are $ θ^1=θ_A $ 'A's, $ θ^2=θ_C $ 'C's,$ θ^3=θ_G $ 'G's, and $ θ^4=θ_T $ 'T's. The equivalence class is the set of all rearrangements of $ (A)^{⊗θ_A}⊗(C)^{⊗θ_C}⊗(G)^{⊗θ_G}⊗(T)^{⊗θ_T} = \bar{B}^{θ_{\bar{B}}} $ is $ 𝔖^{|θ|} \bar{B}^{θ_{\bar{B}}} $. And $Θ $ partitions $ B^⋆ $ because for any $ w ∈ B^⋆ , M_{\bar{B}}w ∈ ℕ^4 $ means that $ M_{\bar{B}}w = ϑ $ for some $ ϑ∈Θ $.
# Note: You can replace $(b)^{⊗θ_b} $ with $ {b}^{×θ_b}$  to get an equivalent expression. The latter is a more set theorhetic approcach. I will avoid this for a very particular reason that I hope to get to later.
# 2) Let $ Φ = \{ ϕ∈ℕ⁴ ∣ 1 ≤ i ≤ j ≤ 4 ⟹ϕᵢ≤ϕⱼ\} $. Each $ ϕ ∈ Φ $ determines a equivalence class which is a union of equivalence classes from 1. Namely, the set of w in $ B^⋆ $ such that w has $ϕ_1 b₁s, ϕ_2 b₂s, ϕ_3 b₃s$, and $ϕ_4 b₄s$, and ${b₁,b₂,b₃,b₄}=B$. Namely, $ ϕ∈Φ $ determines the equivalence class, $ ⋃_{g∈𝔖^4} 𝔖^{|ϕ|} \bar{B}^{g ϕ} $. This partitions $ B^⋆ $ because it is a coarsening of the partiion given in 1. Also, $ ⋃_{ϕ∈Φ} 𝔖ϕ = ℕ⁴ $.
# 2.1) This equivalence class may be given in another way. Let $\widehat{𝔖_B}:B^⋆→ B^⋆$ be given by $⋃_l ⋃_{g∈𝔖_B}g^{×l}$. This is the action of the union over l of the diagonal subgroups of  $𝔖_B^{× l}$ acting on the union of their domains $B^l$. Then the equivalence class given by $ϕ$ from (2),  $⋃_{g∈𝔖_B^4} 𝔖_B^{|\phi|} \bar{B}^{gϕ}$, is the same thing as $𝔖^{|ϕ|}\widehat{𝔖_B}\bar{B}^ϕ$. I love this one for a lot of reasons: $\mathfrak{S}^{|\phi|}$ and $\widehat{𝔖_B}$ commute, the klien 4 group is a subgroup of $ \widehat{𝔖_B}$ if you're willing to associate $\mathbb{Z}_2^2$ with nucleotide bases (i have a very very fun, perhaps computationally useful, idea here with $ℂ$), every subgroup of $\widehat{𝔖_B}$ is  isomorphic to a torus of some dimension. There are more reasons.... 
# 3) Let $ψ = (ψ^{(l)})_{l∈ℕ}$ where $ψ^{(l)} ∈B^l$. Then for any $w ∈ B^l$ we may find an element of $ 𝔖_B^l $, such that $ gψ^{(l)} = w $. There may be multiple such g, so we can actually find a subgroup H of $ 𝔖_B^l $, such that $ Hψ^{(l)} = \{w\} $.
# 4) This one is a riff on 3) and the finite length sequence approach to everything. First, let $ Ψ ∈B^∞ $ be some fixed element. Now, for the rest of this bullet, let w ∈ B^l. Instead of thinking of w as some finite sequence, we could think of it as an equivlence class of sequences $ \tilde{w} ⊂ B^∞ $, where every $ v ∈ \tilde{w} $ has its first $ \# w $  elements equal to w, so $\tilde{w}$ is an element of $ B^∞ / ⟨e^{(1,…,\# v)} v' = v ∣ v∈B^⋆⟩=\tilde{B^⋆} $. Formally, $\tilde{}$ is a function, $\tilde:B^⋆→𝒫(B^∞)$, which is defined by $ w↦\{v∈B^∞∣ (v^1,…,v^{\# w}) = w \}$. Since, $B^∞ = 𝔖_B^∞ Ψ $, there is a subgroup $ H = (H_1,H_2,…) ∈ 𝔖_B^∞ $ such that $ HΨ = \tilde{w}$. The structure of H is fairly simple: $ h ∈ H_1 $ maps $ Ψ^1 $ to $ w^1,…,h ∈ H_{\# w}$ maps $ Ψ^{\# w} $ to $ w^{\# w},h∈H_{\# w + 1} $ maps $ Ψ^{\# w +1}$ to any element of B,…. So for $ 1 ≤ i ≤ \# w, Hᵢ = (Ψ^i \: w^i)𝔖_{ B / \{Ψ^i,w^i\}}I_B $, and for $ i>\#w, Hᵢ = 𝔖_B $.

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

# \section{On Models}
# Acccording to Se Yoon Lee, Bayesian models are given by \{P(y|θ),π(θ)\}. This is supposed to be read as \{P(⋅|θ):Ω→ℝ∣ θ ∈ Θ\}∪\{π:Θ→ℝ\} and interpreted as the collection of conditional probability densities and a prior.
# Se Yoon Lee, Gibbs sampler and coordinate ascent variational inference: A set-theoretical review, https://arxiv.org/pdf/2008.01006.pdf

# According to Sullivant, a parametric statistical model given by a parameter spcae with a family of conditional probability densities, (Θ,\{P_θ:Ω→ℝ∣θ∈Θ\}).

# In either case, Ω is called the outcome space.

# Lets get one thing out of the way:
# \[ P_θ(⋅) = P(⋅|θ) \]

# Can we all just get along and dance?:
# P_θ(⋅) = "Traditional Statistics" = "Left θ-ists"
# P(⋅|θ) = "Bayesian Statistics" = "Right θ-ists" 
# Using, P_θ(⋅) = P(⋅|θ), we arrive at the following:
# \[ "Left θ statistics" = "Traditional Statistics" = P_θ(⋅) = P(⋅|θ) = "Bayesian" = "Right θ statistics" \]

# Why can't we define a model using ℳ:Ω×Θ→ℝ_+? I am not sure, lets see what happens.
# Suppose we define a model by a function ℳ:Ω×Θ→ℝ_+,
# A conditional density given by θ is P_θ(⋅)=P(⋅|θ)=∫_Ω ℳ(⋅,θ)dω = \frac{∂ℳ(⋅,ϑ)}{∂ϑ}|_{ϑ=θ}
# The family of conditional probability densities is \{∫_Ω ℳ(⋅,θ)dω ∣ θ∈Θ\}
# The prior density is π(⋅) = ∫_Ω ℳ(ω,⋅) dω

# The "evidence" of the model is m(⋅) = ∫_Θ P(⋅|θ)π(θ)dθ, huh, ∫_Θ P(⋅|θ)π(θ)dθ = ∫_Θ ∫_Ω ℳ(⋅,θ)dω ∫_Ω ℳ(ω,θ) dω dθ = ∫_Θ ∫_Ω ∫_Ω ℳ(⋅,θ)ℳ(ω,θ)dω dω dθ.

# Maybe consider, m(⋅) = "∫_Θ P(⋅|θ)π(θ)dθ" = "∫_Θ P(⋅,θ) dθ" = ∫_Θℳ(⋅,θ)dθ. This may be easily intereted as a probability density because it is a function from Ω to ℝ_+. Lets take a moment to reflect, we simply specified ℳ as a non-negative function over Ω×Θ and it appears as though the common usage of "P" is precisely this arbitrary function. Huh, we never specified anything more about ℳ, so the number, ∫_Ω∫_Θ ℳ(ω,θ)dωdθ, could really be any number in ℝ_+, which includes 0. (As a side note, lets stop using notation that makes the most basic letter for the most basic concept an arbitrary non-negative function.)

# We observed that ∫_Ω∫_Θ ℳ(ω,θ)dωdθ is a number, and gave names to the following integrals: ∫_Ω ℳ(⋅,θ)dω, ∫_Ω ℳ(ω,⋅)dω, and ∫_Θ ℳ(⋅,θ) dθ.

# There is one more integral of ℳ to consider: ∫_Θ ℳ(ω,⋅)dθ.
# This is called the "posterior distribution", typically written as "π(θ∣ω)", and holds two honors: the most desired thing to compute and the biggest pain in statistics' butt. Lets first observe that we only supposed a parametric statistical model (bayesian or not) was provided by a non-negative function on the product of an (arbitrary) outcome space (Ω) and a (arbitary) parameter space (Θ). The statisticain often knows the parameter space. So to determine the posterior, given some observation, the statistician should evaluate the integral.
# Often, statisticians have multiple observations x_1,…,x_n. The data consisting of these observations is typically denoted, 𝐗. So if the statistician wishes to compute the "posterior given the data", they may compute it ∫_Θ ℳ(𝐗,⋅)dθ; that is, if and only if 𝐗 ∈ Ω.
# In many cases, the data consistsing of observations is not in the form of an observation. The form that 𝐗 takes is up to the modeller's discretion.
# The data, 𝐗, could be:
# a set S = \{x_1,…,x_n\}
# a set function f:S→ℕ counting the number of times an element of S occured in the sample, i.e. f(x) = ∑_{i=1}^n 1_{x=x_i}
# or a n-tuple, e.g. (x_1,…,x_n), e.g.(x_n,…,x_1)
# these are some possibilites when Ω is a set.
# In the case that Ω is a set, we could consider the data to be given by any sum or product or really any binary operation, e.g. +,*,⨁,⨂,∨,∧,∘,[⋅,⋅].

# It seems like most of the difficulties encountered in the bayesian notation come from blatant misuse of bayes' rule. If P:Ω→[0,1], then the probability of A given B may be computed using the formula for conditional probability:
# \[ P(A∩B) = P(A|B)P(B) \]
# If you apply this formula to $P(A∩B)$ and $P(B∩A)$, you end up with Bayes' rule:
# \[ A∩B=B∩A ⟹ P(A∩B)=P(B∩A) ⟹ P(A|B)P(B)=P(B|A)P(A) ⟹ P(A|B)=\frac{P(B|A)P(A)}{P(B)} \]
# Bayesians are understandably excited to use Bayes' rule. So when tasked with computing the posterior distribution, typically written as "P(θ|𝐗)", they jump to compute \frac{ P(𝐗|θ)P(θ)}{P(𝐗)}. If you replace A with θ, replace B with 𝐗, and reverse the chain of implications (using the modus pomens rule in logic) you find:
#\[ P(θ|𝐗)= \frac{P(𝐗|θ)P(θ)}{P(𝐗)} ⟹ P(θ|𝐗)P(𝐗)= P(𝐗|θ)P(θ) ⟹ P(θ∩𝐗)=P(𝐗∩θ) ⟹ θ∩𝐗=𝐗∩θ\]
# In bayesian models, θ and 𝐗 tend to live in different spaces, so it tends tobe the case that $ θ∩𝐗=𝐗∩θ=∅$.

# I am going to be a bad bayesian. Let P = ℙ and π = ℼ. ℼ = 1/n∑_i δ_{N_\bar{B}(X_i)}
#(θ|x) = P(x|θ)π(θ) / m(x)

# eⁱ(n_1e_1+…+n_ke_l) = n_i
# for some orthogonal (?) e_b \in V_B and E = (e_A,e_C,e_G,e_T), the dual of E is \hat{E} = (e^1,e^2,e^3,e^4).

# EᵀÊ = I_l = J_l*1_l
# If we let E_w = (e_{w¹},…,e_{wˡ}); 
# maybe let ⊗E_w = e_{(w¹,…,wˡ)} and ⊗Ê^ℕ = e^{(1,…,l)}
# ⊗Eᵀ⊗Ê = λ_{≤l} = ∑_i e_{w^i}
# e.g. E_w = e_


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