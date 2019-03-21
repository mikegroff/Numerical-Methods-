f = @(x) sech(x); fp = @(x) sec(x)*tan(x); g = @(x) (abs(x).^2).*x; L = 20; N = 2^7; 
h = L/N; k = 1d-4; err = 0; tol = 1d-5; nt = 1/k; U = zeros(1,N+1); levelcap = 15;
d = -2*ones(1,N+1); a = ones(1,N);
A = diag(d) + diag(a,-1) + diag(a,1);
A(N+1,1) = 1; A(1,N+1) = 1;
M = [];

for i = 0:N
    p = -L/2 + i*h;
    U(1,i+1) = f(p);
end
lvl=1;
while err < tol
    prev = U(lvl,:)';
    if lvl >= levelcap; break; end
        for j = 1:nt
            prev = prev + 1j*((k/(h^2))*A*prev + (2*k)*g(prev));
        end
        
        lvl = lvl+1;
        U(lvl,:) = prev';
end

mesh(abs(U));

