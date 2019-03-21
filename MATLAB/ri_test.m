f = @(x) sin(x); fp = @(x) cos(x); fpp = @(x) -sin(x);
x = pi/3; h = 0.1; n =4;
phi = @(x) (f(x+h)-f(x-h))/(2*h);
D = RI(phi,h,n);
disp([abs(fp(x)-D(1,:))]);