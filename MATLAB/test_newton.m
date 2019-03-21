f = @(x) cos(3*x) - cos(7*x);
fp = @(x) 7*sin(7*x) - 3*sin(3*x);
x0 = 1; tol = 0.5d-6; 
[x,n] = newton(f,fp,x0,tol)
