clear;
close all;

L = 4*sqrt(2)*pi;
f = @(x) 0.5*(ones(size(x))+0.1*cos(2*pi*x/L)+0.0375*cos(4*pi*x/L)); 
%f = @(x) sech(x); 
g = @(x) (abs(x).^2).*x; s = @(x) (abs(x).^2);
pp = 7; N = 2^pp; emp = [];
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
        holds = prev; count3 = count3+1;
        errr = 1; lvlv = 1; 
        while errr > tol
            comp = prev;
            prev = D1*(D2*holds + (1j*k)*g(prev+holds)/4);
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
hold on;




h = L/N; k = (1/4)*h^2; err = 0; tol = 1d-15; nt = 1/k; 
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
        holds = prev; errr = 1; lvlv = 1; count2 = count2+1;
        while errr > tol
            comp = prev;
            prev = B*(holds + 1i*(2*k)*g(prev));
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

ntt = nt;

h = L/N; k = 1d-4; err = 0; tol = 1d-8; nt = 1/k; U = zeros(1,N+1); levelcap = 10;
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

c = 1:count; cc = 1:count2; ccc = 1:count3;

figure;
set(gcf,'Color', 'w');
semilogy(c/nt,abs(Mm(1)-Mm(c)),'r',cc/ntt,abs(Mm2(1)-Mm2(cc)),'g',ccc/ntt,abs(Mm3(1)-Mm3(ccc)),'b');

title('Log Error of M Constant for Two Bump');
xlabel('t'); ylabel('error');
legend('NLS1', 'NLS2','C-N');

figure;
set(gcf,'Color', 'w');
semilogy(c/nt,abs(Nn(1)-Nn(c)),'r',cc/ntt,abs(Nn2(1)-Nn2(cc)),'g',ccc/ntt,abs(Nn3(1)-Nn3(ccc)),'b');
title('Log Error of N Constant for Two Bump');
xlabel('t'); ylabel('error');
legend('NLS1', 'NLS2','C-N');

hold off;

function [retval] = simpsons(f,a,b,n)
    h = (b-a)/n; tot = 0;
    for k= 1:n-1
        tot = tot+ h*(f(k)+4*f(k+1)+ f(k+2))/6;
    end
    retval = tot;
end