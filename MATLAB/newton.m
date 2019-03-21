function [x,n] = newton(f,fp,x0,tol)
error = 1; n = 0; x = x0;
 
while error > tol
    n = n+1;
    x = x - f(x)/fp(x);
    error = abs(x-x0);  
    x0 = x;
end

