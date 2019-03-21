f = @(t,x) x+t; 
f1 = @(t,x) f(t,x) + 1; %f2 = @(t,x) f1(t,x); f3= @(t,x) f2(t,x);
%T1 = @(t,x,h) x + h*f(t,x); 
T2 = @(t,x,h) x + h*f(t,x) + h^2/2*f1(t,x);
%T3 = @(t,x,h) x + h*f(t,x) + h^2/2*f1(t,x)+ h^3/6*f2(t,x);
%T4 = @(t,x,h) x + h*f(t,x) + h^2/2*f1(t,x)+ h^3/6*f2(t,x)+h^4/24*f3(t,x);
x0 = 0; k = 0; E = zeros(6,4); h = 0.1; t = 0; rk = 0; s2 = x0;%exact soln
    for n = 1:1/h
        s2(n+1) = T2(t(n),s2(n),h);
        t(n+1) = t(n)+h;
        k1 = f(t(n),rk(n));
        k2 = f(t(n) + h, rk(n) + h*k1);
        rk(n+1) = rk(n) + h/2*(k1+k2);
    end
s = -t -1 + exp(t);
subplot(2,1,1), plot(t,s,'r--',t,s2,t,rk);
title("The Solutions");
xlabel('t'); ylabel('x(t)');
subplot(2,1,2), semilogy(t,abs(s-s2),'o',t,abs(s-rk));
title("The Errors");
xlabel('t'); ylabel('x(t)');