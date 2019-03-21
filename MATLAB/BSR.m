%Simpson's rule
f = @(x) sin(x);
A =2; h = pi/2;
SR = h/3*(f(0) +4*f(h) + f(pi));
disp(abs(A-SR));
N = 2; h = pi/N; error =1; tol = .5d-6;
while error > tol
    N= 2*N; h = pi/N;
    sr = h/3*f(0) + f(pi) +4*sum(sin( (1:2:N-1)*h))+...
        2*sum(sin( (1:2:N-2)*h));
    error = abs(2-sr);
end

