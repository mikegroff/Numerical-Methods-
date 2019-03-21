% testing TR
f = @(x) 2/sqrt(pi)*exp(-x.^2); a = 0; b = 1; N = sqrt((1d6)/3);
h = (b-a)/N; xn = a:h:b; fn = f(xn);
s = TR(fn,h); [s erf(b)],
% ii) check how the error decreases when h = 1/2,1/4,...
h = 1; s0 = erf(b);
for n=1:10
    h = h/2; xn = a:h:b; fn = f(xn); s = TR(fn,h);
    H(n) = h; E(n) = abs(s - s0);
end
% if E ~ CH^2 then to match the graphs we need C = E(1)/H(1)^2;
C = E(1)/H(1)^2;
loglog(H,E,H,C*H.^2,'ro')