clear; close all
for q = 1:3
    for gq=1:2
        t0 = 5*q-15;
        f = @(x) per1(x,t0);
        %f = @(x) 0.5*(ones(size(x))+0.1*cos(2*pi*x/L) + 0.0375*cos(4*pi*x/L));
        g = @(x) (abs(x).^2).*x; s = @(x) (abs(x).^2);
        L = 80/gq; N = 2^7; emp = []; 
        h = L/N; k = 1d-4; err = 0; tol = 1d-8; nt = 1/k; U = zeros(1,N+1); levelcap = 15;
        d = -2*ones(1,N+1); a = ones(1,N); xn = zeros(1,N+1);
        A = diag(d) + diag(a,-1) + diag(a,1);
        A(N+1,1) = 1; A(1,N+1) = 1;
        Ax = diag(d/2) + diag(a,1);
        Ax(N+1,1) = 1; Ax(1,N+1) = -1;
        Mm = zeros(1,ceil(levelcap*nt)); Nn = zeros(1,ceil(levelcap*nt));
        
        p = L*cos(pi*(0:N)/N);
        xn = p;
        U(1,:) = f(p);
        
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
        title(['Numerical approxiamtion of NLS1 (L = ' num2str(L) ', t0 = ' num2str(t0) ]);
        xlabel('x'); ylabel('t'); zlabel('|U(x,t)|');
        
        c = 1:count;
        disp(count);
        figure;
        set(gcf,'Color', 'w');
        yyaxis left
        semilogy(c/nt,abs(Mo-Mm(c)),'r');
        yyaxis right
        semilogy(c/nt,abs(No-Nn(c)),'b');
        title(['Log Error of Constants of Motion t0 = ' num2str(t0)]);
        xlabel('t'); ylabel('error');
        legend('M constant', 'N constant');
        
        figure;
        set(gcf,'Color', 'w');
        yyaxis left
        plot(c/nt,Mm(c))
        yyaxis right
        plot(c/nt,Nn(c));
        title(['Constants of Motion t0 = ' num2str(t0)]);
        xlabel('t'); 
        legend('M constant', 'N constant');
        Maxvals = zeros(1,lvl);
        for r = 1:lvl
            Maxvals(r) = max(abs(U(r,:)));
        end
        
        a = 1:lvl;
        
        figure;
        set(gcf,'Color', 'w');
        plot( a,Maxvals(a));
    end
end

function [retval] = simpsons(f,a,b,n)
h = (b-a)/n; tot = 0;
for k= 1:n-1
    tot = tot+ h*(f(k)+4*f(k+1)+ f(k+2))/6;
end
retval = tot;
end

function u = per1(x,t)
a = 1/2;  t = t;
%u = a*exp(2*a^2*i*t).*(1 - 4*( 1 + 4*i*t )./( 1 + 16*t.^2 + 4*x.^2 ));
% if q(x,t) is a solution then u(x,t) = a q(ax,a^2*t) is a solution too
% From Dysthe & Trulsen:
% Physica Scripta (T82) 48-52 (1999)
%
% q(x,t) = e^{2it}[ 1 - 4(1+4it)/(1 + 16t^2 + 4x^2) ]
%
% so
u = a*exp( 2*i*(a^2*t) ).*(1 - 4*( 1 + 4*i*(a^2*t) )./( 1 + 16*(a^2*t).^2 + 4*(a*x).^2 ));
end
function u = per2(x,t)
  a = 1/2; t = a^2*t; x = a*x;
g = -12*( x.^4 + 6*(t.^2 + 1).*x.^2 + 5*t.^4 + 18*t.^2 - 3 );
h = -12*( x.^4 + 2*(t.^2 - 3).*x.^2 + (t.^2 + 5).*(t.^2 - 3) );
f = x.^6 + 3*(t.^2 + 1).*x.^4 + 3*(t.^2 - 3).^2.*x.^2 + t.^6 + 27*t.^4 ...
    + 99*t.^2 + 9;
u = a*exp(2*i*t).*( 1 + (g + i*t.*h)./f );
end