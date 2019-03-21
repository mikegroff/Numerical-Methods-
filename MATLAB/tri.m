function x = tri(a,b,c,d)

n = length(d);

for k=2:n
    m = a(k-1)/d(k-1);
    d(k) = d(k) - m*c(k-1);
    b(k) = b(k)- m*b(k-1);
end

x(n) = b(n)/d(n);

for k=n-1:-1:1
    x(k) = (b(k)-c(k)*x(k+1))/d(k);
end