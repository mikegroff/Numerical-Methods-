% testing RI
f = @(x) sin(x); fp = @(x) cos(x); fpp = @(x) -sin(x);
x = pi/3; h = 0.1; n = 4;
phi = @(h) ( f(x+h) -2*f(x) + f(x-h) )/(h^2);
D = RI(phi,h,n);
disp([abs(fpp(x) - D(1,:))])