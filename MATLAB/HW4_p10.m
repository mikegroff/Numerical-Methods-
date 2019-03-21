f = @(x) 1/x; N = 16; a = 1; b = 2; A = zeros(N,N);
h = (b-a)/2;
A(1,1) = h*(f(b)+f(a));

for i=2:N
    h = h/2;
    S = 0;
    for j = 1:2:2^i-1
        S =  S + f(a + j*h);
    end
    A(i,1) = A(i-1,1)/2 + h*S;
    
    for k = 2:i
        A(i,k) = A(i,k-1) + (A(i,k-1)-A(i-1,k-1))/(4^k -1);
    end
    
end

disp(A(N,N));