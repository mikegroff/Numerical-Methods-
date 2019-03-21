function x = NG(A,b)
% Gauss elimination: forward elimination and backward substitution
n = length(b);
% forward elimination
for k = 1:n-1
    for i = k+1:n
        m = A(i,k)/A(k,k);
        for j = k:n
            A(i,j) = A(i,j) - m*A(k,j);
        end
        b(i) = b(i) - m*b(k);
    end
end
% back substitution
x(n) = b(n)/A(n,n);
for i = n-1:-1:1
    x(i) = (b(i) - sum(A(i,i+1:n).*x(i+1:n)))/A(i,i);
end
x = x(:);