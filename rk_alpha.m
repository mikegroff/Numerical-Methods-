aa = 100;
f = @(t,x) aa*(sin(t) -x); 
f1 = @(t,x) aa*(cos(t)-f(t,x)); 
T2 = @(t,x,h) x + h*f(t,x) + h^2/2*f1(t,x);
x0 = 0; k = 0; h = 0.01; t = 0; rk = 0; s2 = x0; 
for alpha = 0.1:0.1:1
beta = alpha; b = 1/2/alpha; a = 1-b;
    for n = 1:1/h
        s2(n+1) = T2(t(n),s2(n),h);
        t(n+1) = t(n)+h;
        k1 = h*f(t(n),rk(n));
        k2 = h*f(t(n) + alpha*h, rk(n) + beta*k1);
        rk(n+1) = rk(n) +a*k1+b*k2;
    end
s = aa/(1+aa^2)*(aa*sin(t) - cos(t) + exp(-aa*t));
subplot(2,1,1), plot(t,s,'r--',t,s2,t,rk);
title("The Solutions");
xlabel('t'); ylabel('x(t)');
subplot(2,1,2), semilogy(t,abs(s-s2),'r--',t,abs(s-rk));
title("The Errors");
xlabel('t'); ylabel('x(t)');
pause(1)
end