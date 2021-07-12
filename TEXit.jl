import REPL

typeface = ["\\bf","\\bb","\\scr","\\frak"]

symbols_latex = Pair{String,String}[]
for p in REPL.REPLCompletions.latex_symbols
	if p[1][2] in ["_","^"]
		push!(symbols_latex, p[2]=>p[1][2:end])
	elseif any(occursin.(typeface,p[1]))
		f = findfirst(x->occursin(x,p[1]),typeface)
		i = length(typeface[f])-1
		texcmd = "\\math"*p[1][2:(2+i-1)]*"{"*p[1][(2+i):end]*"}"
		push!(symbols_latex, p[2]=>texcmd)
	else
		push!(symbols_latex, p[2]=>p[1])
	end
end

function old_new(x)
	if x in first.(symbols_latex)
		return last.(symbols_latex)[findfirst(==(x),first.(symbols_latex))]
	else
		return x
	end
end

file = read(ARGS[1],String);
texedFile = mapreduce(old_new,*, "" .* collect(file));
write(ARGS[1][1:(end-2)] * "tex", texedFile)



