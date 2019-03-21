%f = @(x) 0.5*(ones(size(x))+0.1*cos(2*pi*x/L)); 
f = @(x) sech(x); 
g = @(x) (abs(x).^2).*x; s = @(x) (abs(x).^2);
pp = 7; L = 20; N = 2^pp; emp = [];
h = L/N; k = (1/4)*h^2; err = 0; tol = 1d-15; nt = 1/k; 
%The V matrix will hold each nt step of U
w = zeros(1,N+1); levelcap = 30; xn = zeros(1,N+1);
d = -2*ones(1,N+1); a = ones(1,N); i = ones(1,N+1);
%need identity to calculate the B
I = diag(i);
A = diag(d) + diag(a,-1) + diag(a,1);
A(N+1,1) = 1; A(1,N+1) = 1;
D1 = I - (1j*k/(2*h^2))*A;
D2 = I + (1j*k/(2*h^2))*A;
D1 = D1^(-1);
Ax = diag(-a,1) + diag(a,-1);
Ax(N+1,1) = 1; Ax(1,N+1) = -1;
%arrays to hold constants of motion
Mm3 = zeros(1,ceil(levelcap*nt)); Nn3 = zeros(1,ceil(levelcap*nt));

for i = 0:N
    p = -L/2 + i*h;
    xn(i+1) = p;
    w(1,i+1) = f(p);
end
lvl=1; 
Vx = Ax*(w(lvl,:)')/(2*h);
Nf = s(s(w(lvl,:))) - s(Vx');
Mo = simpsons(s(w(lvl,:)),-L/2,L/2,N);
No = simpsons(Nf,-L/2,L/2,N);
Mm3(1) = Mo; Nn3(1) = No;
count3 = 1; 
while err < 2
    if lvl >= levelcap; break; end
    prev = w(lvl,:)'; 
    for l = 1:nt
        hold = prev; count3 = count3+1;
        errr = 1; lvlv = 1; 
        while errr > tol
            comp = prev;
            prev = D1*(D2*hold + (1j*k)*g(prev+hold)/4);
            lvlv = lvlv+1;
            %disp(lvlv);
            errr = norm(prev - comp);
            if lvlv >= levelcap; break; end
        end
        %evaluating the constants at every k step
        Vx = Ax*(prev)/(2*h);
        Nf = s(s(prev')) - s(Vx');
        Mm3(count3) = simpsons(s(prev'),-L/2,L/2,N);
        Nn3(count3) = simpsons(Nf,-L/2,L/2,N);
        
    end
        
        lvl = lvl+1;
        %saying scheme after every nt steps
        w(lvl,:) = prev';
        %err = max(abs(No-Nn(lvl)),abs(Mo-Mm(lvl)));
        err = max(w(lvl,:));
        disp(lvl);
end
figure;
set(gcf,'Color', 'w');
mesh( xn, 1:levelcap, abs(w));
title('Crank-Nicholson approxiamtion of NLS');
xlabel('x'); ylabel('t'); zlabel('|U(x,t)|');

c = 1:count3;

figure;
set(gcf,'Color', 'w');
yyaxis left
semilogy(c/nt,abs(Mo-Mm3(c)),'r');
yyaxis right
semilogy(c/nt,abs(No-Nn3(c)),'b');
title('Error of Constants of Motion');
xlabel('t'); ylabel('error');
legend('M constant', 'N constant');

figure;
set(gcf,'Color', 'w');
yyaxis left
plot(c/nt,Mm3(c))
yyaxis right
plot(c/nt,Nn3(c));

Maxvals = zeros(1,lvl);
for r = 1:lvl
    Maxvals(r) = max(abs(w(r,:)));
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