function [x,n] = bisection(f,a,b,tol)

if f(a)*f(b) > 0
    x = -1000; disp('Wrong end points'); return; 
end
error = 1; fa = f(a); fb = f(b); n = 0;

while error > tol
    n = n+1; c = (a+b)/2;
    if fa*f(c) > 0
        a = c; fa = f(a);
    else 
        b = c; fb = f(b);
    end
    error = (b-a)/2;  
end
x = (a+b)/2;