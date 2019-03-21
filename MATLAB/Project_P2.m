%f = @(x) 0.5*(ones(size(x))+0.1*cos(2*pi*x/L)); 
f = @(x) sech(x); 
g = @(x) x.*(abs(x).^2); s = @(x) (abs(x).^2);
pp = 7; L = 20; N = 2^pp; emp = [];
h = L/N; k = (1/4)*h^2; err = 0; tol = 1d-8; nt = 1/k; 
%The V matrix will hold each nt step of U
V = zeros(1,N+1); levelcap = 30; xn = zeros(1,N+1);
d = -2*ones(1,N+1); a = ones(1,N); i = ones(1,N+1);
%need identity to calculate the B
I = diag(i);
A = diag(d) + diag(a,-1) + diag(a,1);
A(N+1,1) = 1; A(1,N+1) = 1;
B = I - 1i*(k/(h^2))*A;
B = inv(B);

Ax = diag(-a,1) + diag(a,1);
Ax(N+1,1) = 1; Ax(1,N+1) = -1;
%arrays to hold constants of motion
Mm2 = zeros(1,ceil(levelcap*nt)); Nn2 = zeros(1,ceil(levelcap*nt));

for i = 0:N
    p = -L/2 + i*h;
    xn(i+1) = p;
    V(1,i+1) = f(p);
end
lvl=1; 
Vx = Ax*(V(lvl,:)')/(2*h);
Nf = s(s(V(lvl,:))) - s(Vx');
Mo = simpsons(s(V(lvl,:)),-L/2,L/2,N);
No = simpsons(Nf/2,-L/2,L/2,N);
Mm2(1) = Mo; Nn2(1) = No;
count2 = 1; 

while err < 2
    prev = V(lvl,:)'; 
    if lvl >= levelcap; break; end    
    for l = 1:nt    
        hold = prev; errr = 1; lvlv = 1; count2 = count2+1;
        while errr > tol
            comp = prev;
            prev = B*(hold + 1i*(2*k)*g(prev));
            lvlv = lvlv+1;
            errr = norm(prev - comp);
            if lvlv >= levelcap; break; end
        end
        %evaluating the constants at every k step
        Vx = Ax*(prev)/(2*h);
        Nf = s(s(prev')) - s(Vx');
        Mm2(count2) = simpsons(s(prev'),-L/2,L/2,N);
        Nn2(count2) = simpsons(Nf/2,-L/2,L/2,N);
        
    end
        
        lvl = lvl+1;
        %saving scheme after every nt steps
        V(lvl,:) = prev';
        %err = max(abs(No-Nn(lvl)),abs(Mo-Mm(lvl)));
        err = max(V(lvl,:));
        disp(lvl);
end
figure;
set(gcf,'Color', 'w');
mesh( xn, 1:levelcap, abs(V));
title('Crank-Nicholson approxiamtion of NLS');
xlabel('x'); ylabel('t'); zlabel('|U(x,t)|');

c = 1:count2;

figure;
set(gcf,'Color', 'w');
yyaxis left
semilogy(c/nt,abs(Mo-Mm2(c)),'r');
yyaxis right
semilogy(c/nt,abs(No-Nn2(c)),'b');
title('Error of Constants of Motion');
xlabel('t'); ylabel('error');
legend('M constant', 'N constant');

figure;
set(gcf,'Color', 'w');
yyaxis left
plot(c/nt,Mm2(c))
yyaxis right
plot(c/nt,Nn2(c));

Maxvals = zeros(1,lvl);
for r = 1:lvl
    Maxvals(r) = max(abs(V(r,:)));
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




