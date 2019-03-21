import containers.*

f = @(x) sech(x); fp = @(x) sec(x)*tan(x); g = @(x) (abs(x).^2).*x; s = @(x) (abs(x).^2);
L = 20; N = 2^7; emp = [];
h = L/N; k = 1d-4; err = 0; tol = 1d-1; nt = 1/k; U = zeros(1,N+1); levelcap = 15;
d = -2*ones(1,N+1); a = ones(1,N);
A = diag(d) + diag(a,-1) + diag(a,1);
A(N+1,1) = 1; A(1,N+1) = 1;
Ax = diag(d/2) + diag(a,1);
Ax(N+1,1) = -1; Ax(1,N+1) = -1;
Mm = zeros(1,levelcap); Nn = zeros(1,levelcap);

for i = 0:N
    p = -L/2 + i*h;
    U(1,i+1) = f(p);
end
lvl=1;
Ux = Ax*(U(lvl,:)')/h;
Nf = s(s(U(lvl,:))) - s(Ux');
G = containers.Map((-L/2:h:L/2),Nf);
H = containers.Map((-L/2:h:L/2),s(U(lvl,:))); 
Mo = ASimpson(H,-L/2,L/2,tol, 0, 6, emp,emp);
No = ASimpson(G,-L/2,L/2,tol, 0, 6, emp,emp);
Mm(1) = Mo; Nn(1) = No;
while err < tol
    prev = U(lvl,:)';
    if lvl >= levelcap; break; end
        for j = 1:nt
            prev = prev + 1j*((k/(h^2))*A*prev + (2*k)*g(prev));
        end
        
        lvl = lvl+1;
        U(lvl,:) = prev';
        H = containers.Map((-L/2:h:L/2),s(U(lvl,:))); 
        Ux = Ax*(U(lvl,:)');
        Nf = s(s(U(lvl,:))) - s(Ux');
        G = containers.Map((-L/2:h:L/2),Nf);
        Mm(lvl) = ASimpson(H,-L/2,L/2,tol, 0, 6, emp,emp);
        Nn(lvl) = ASimpson(G,-L/2,L/2,tol, 0, 6, emp,emp);
        %err = max(abs(No-Nn(lvl)),abs(Mo-Mm(lvl)));
end
subplot(2,1,1);
mesh(abs(U));
subplot(2,1,2);
c = 1:lvl;
hold on;
semilogy(c,abs(Mo-Mm(c)),'r');
semilogy(c,abs(No-Nn(c)),'b');
Hold off;