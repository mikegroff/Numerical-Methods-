function [x,n] = secant(f,x0,x1,tol)
error = 1; n = 0; x = x1;
while error > tol
    n = n+1;
    x = x - f(x)*(x-x0)/(f(x) - f(x0));
    error = abs(x-x1);
    x0 = x1; x1 = x; 
end
