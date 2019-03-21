f = @(x) cos(3*x) - cos(7*x);
a = 0; b = 1; tol = 0.5d-6; 
[x,n] = bisection(f,a,b,tol)
