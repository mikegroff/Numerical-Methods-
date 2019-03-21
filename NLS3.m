clear; close all
for q = 1:3
    for gq=1:2
        t0 = 5*q-15;
        f = @(x) per1(x,t0);
        g = @(x) (abs(x).^2).*x; s = @(x) (abs(x).^2);
        pp = 8; L = 80/gq; N = 2^pp; emp = []; Y = 2;
        h = L/N; k = h^2/(2^Y); err = 0; tol = 1d-15; nt = 1/k;
        %The V matrix will hold each nt step of U
        w = zeros(1,N+1); levelcap = 15; xn = zeros(1,N+1);
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
        
        p = L*cos(pi*(0:N)/N);
        xn = p;
        w(1,:) = f(p);
        
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
         title(['Numerical approxiamtion of NLS3 (L = ' num2str(L) ', t0 = ' num2str(t0) ]);
        xlabel('x'); ylabel('t'); zlabel('|U(x,t)|');
        
        c = 1:count3;
        
        figure;
        set(gcf,'Color', 'w');
        yyaxis left
        semilogy(c/nt,abs(Mo-Mm3(c)),'r');
        yyaxis right
        semilogy(c/nt,abs(No-Nn3(c)),'b');
        title(['Log Error of Constants of Motion t0 = ' num2str(t0)]);
        xlabel('t'); ylabel('error');
        legend('M constant', 'N constant');
       
        
        figure;
        set(gcf,'Color', 'w');
        yyaxis left
        plot(c/nt,Mm3(c))
        yyaxis right
        plot(c/nt,Nn3(c));
        title(['Constants of Motion t0 = ' num2str(t0)]);
        xlabel('t'); 
        legend('M constant', 'N constant');
        
        
        Maxvals = zeros(1,lvl);
        for r = 1:lvl
            Maxvals(r) = max(abs(w(r,:)));
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