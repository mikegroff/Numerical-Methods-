n = 1000;
h = 2*pi/n;
x = 0:h:2*pi; f = sin(x) ;
d = -2*ones(n+1,1)/(h^2); a = ones(n,1)/(h^2); c = a;
b = f;
F = tri(a,b,c,d);
plot(x, -sin(x), x , F);