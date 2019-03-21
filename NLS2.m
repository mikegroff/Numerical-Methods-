clear; close all
for q = 1:3
    for gq=1:2
        t0 = 5*q-15;
        f = @(x) per1(x,t0);
        g = @(x) x.*(abs(x).^2); s = @(x) (abs(x).^2);
        pp = 8; L = 80/gq; N = 2^pp; emp = [];
        h = L/N; k = h^2/4; err = 0; tol = 1d-8; nt = 1/k;
        %The V matrix will hold each nt step of U
        V = zeros(1,N+1); levelcap = 15; xn = zeros(1,N+1);
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
        
        p = L*cos(pi*(0:N)/N);
        xn = p;
        V(1,:) = f(p);
        
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
        title(['Numerical approxiamtion of NLS2 (L = ' num2str(L) ', t0 = ' num2str(t0) ]);
        xlabel('x'); ylabel('t'); zlabel('|U(x,t)|');
        
        c = 1:count2;
        
        figure;
        set(gcf,'Color', 'w');
        yyaxis left
        semilogy(c/nt,abs(Mo-Mm2(c)),'r');
        yyaxis right
        semilogy(c/nt,abs(No-Nn2(c)),'b');
        title(['Log Error of Constants of Motion t0 = ' num2str(t0)]);
        xlabel('t'); ylabel('error');
        legend('M constant', 'N constant');
        
        figure;
        set(gcf,'Color', 'w');
        yyaxis left
        plot(c/nt,Mm2(c))
        yyaxis right
        plot(c/nt,Nn2(c));
        title(['Constants of Motion t0 = ' num2str(t0)]);
        xlabel('t'); 
        legend('M constant', 'N constant');
        
        Maxvals = zeros(1,lvl);
        for r = 1:lvl
            Maxvals(r) = max(abs(V(r,:)));
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