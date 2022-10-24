Statistics with Genetic and Genomic Data
========================================
Data: "all" GSAID data accessible on 07/12/2021, NCBI covid data
Stats: see WhatIveDone.txt
Visuals: viz folder
Languages: Julia


# Must Do List
- Make FieldVals.jl produce an output file with everything it is currently printing (lol this is what happens when you write files by copying the code you write in the REPL )
- (GeneticEDA.jl) clean up, add comments, produce an output with computations, and make the commented code produce visuals in the viz folder (FYI)
- Merge AllDataEDA.jl into GeneticEDA.jl (also "AllDataEDA.jl" is a bad name)
- Combine tex documents into one.
- Write GenomicEDA.jl
- Do you need OperatorDefns.jl ?
- Utility file?
	- Do you need NewtonsMethod.jl ? It is also in GeneticEDA.jl.
	- Maybe merge NewtonsMethod.jl with FiniteSupportedFunction.jl, which tbh idk if you need...
	- Every F::FiniteSuppFn{It,Vt} satisfies all( !.(iszero.(keys(F))) ) and supports the same indexing patterns as Dict{It,Vt} except...
		- if x isa It and iszero(y), F will be unchanged by F[x]=y
		- if z isa IT and not a key, F[z] = zero(Vt)
		- in REPL
			> x isa & iszero(y)
			true
			> G = F
			> F[x] = y;
			> F[x]
			


# Data
This project uses all of the data available from GSAID on July 12th, 2021.
- You cannot download all of their genomic data at once.
- You can download at most 1,000 genomic observation at a time.

- I downloaded what they consider "all of it".
	- This is a zipped .fasta file is the output of a process that the original observations.
	- The process chops up a genomic observation into many peices.
	- Each peice has the demographic data of the genomic observation and a subsequence of the genome associated with the genomic observation.
	- Some bases in the observed genome do not make it into any of the resulting peices.
- I downloaded over 230,000 genomic observations.

I asked GSAID to let me download all of the (unprocessed) genomic observations at once. GSAID did not respond to me. GSAID has bad policies.

The data is not included in this repository because one must agree not to share data in order to access it. If you GSAID credentials, I am happy to share the 230,000 genomic observations with you.

My code will reflect the fact I have rewritten .fasta files into a different format (GenomeFastaToObs.jl). The resulting file requires half as much storage.

This project also uses from NCBI (National Center for Biotechnology Innovation)



