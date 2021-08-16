const K = 2,805,095
const N = 61,065,458
const tol = 0.0001
F(x) = K/x - sum([1/(x+j-1) for j in 1:N]) 
f(x) = -K/x^2 + sum([1/(x+j-1)^2 for j in 1:N])
x = 1
while abs(f(x))>=tol
	x = x - F(x)/f(x)
end

