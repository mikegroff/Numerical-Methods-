f = @(x) sech(x); 
%f = @(x) 0.5*(ones(size(x))+0.1*cos(2*pi*x/L) + 0.0375*cos(4*pi*x/L)); 
g = @(x) (abs(x).^2).*x; s = @(x) (abs(x).^2);
L = 20; N = 2^7; emp = [];
h = L/N; k = 1d-4; err = 0; tol = 1d-8; nt = 1/k; U = zeros(1,N+1); levelcap = 30;
d = -2*ones(1,N+1); a = ones(1,N); xn = zeros(1,N+1);
A = diag(d) + diag(a,-1) + diag(a,1);
A(N+1,1) = 1; A(1,N+1) = 1;
Ax = diag(d/2) + diag(a,1);
Ax(N+1,1) = 1; Ax(1,N+1) = -1;
Mm = zeros(1,ceil(levelcap*nt)); Nn = zeros(1,ceil(levelcap*nt));

for i = 0:N
    p = -L/2 + i*h;
    xn(i+1) = p;
    U(1,i+1) = f(p);
end
lvl=1; count = 1;
Ux = Ax*(U(lvl,:)')/(2*h);
Nf = s(s(U(lvl,:))) - s(Ux');
Mo = simpsons(s(U(lvl,:)),-L/2,L/2,N);
No = simpsons(Nf/2,-L/2,L/2,N);
Mm(1) = Mo; Nn(1) = No;
while err < tol
    prev = U(lvl,:)';
    if lvl >= levelcap; break; end
        for j = 1:nt
            count = count + 1;
            prev = prev + 1j*((k/(h^2))*A*prev + (2*k)*g(prev));
            Vx = Ax*(prev)/(2*h);
            Nf = s(s(prev')) - s(Vx');
            Mm(count) = simpsons(s(prev'),-L/2,L/2,N);
            Nn(count) = simpsons(Nf/2,-L/2,L/2,N);
        end
        
        lvl = lvl+1;
        U(lvl,:) = prev';
        disp(lvl);
        %err = max(abs(No-Nn(lvl)),abs(Mo-Mm(lvl)));
end
figure;
set(gcf,'Color', 'w');
mesh( xn, 1:levelcap, abs(U));
title('Numerical approxiamtion of NLS');
xlabel('x'); ylabel('t'); zlabel('|U(x,t)|');

c = 1:count;
disp(count);
figure;
set(gcf,'Color', 'w');
yyaxis left
semilogy(c/nt,abs(Mo-Mm(c)),'r');
yyaxis right
semilogy(c/nt,abs(No-Nn(c)),'b');
title('Error of Constants of Motion');
xlabel('t'); ylabel('error');
legend('M constant', 'N constant');

figure;
set(gcf,'Color', 'w');
yyaxis left
plot(c/nt,Mm(c))
yyaxis right
plot(c/nt,Nn(c));

Maxvals = zeros(1,lvl);
for r = 1:lvl
    Maxvals(r) = max(abs(U(r,:)));
end

a = 1:lvl;

figure;
set(gcf,'Color', 'w');
plot( a,Maxvals(a));

function [retval] = simpsons(f,a,b,n)
    h = (b-a)/n; tot = 0;
    for k= 1:n-1
        tot = tot+ h*(f(k)+4*f(k+1)+ f(k+2))/6;
    end
    retval = tot;
end