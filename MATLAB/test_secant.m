f = @(x) cos(3*x) - cos(7*x);
a = 1/2; b = 1; tol = 0.5d-6; 
[x,n] = secant(f,a,b,tol)
