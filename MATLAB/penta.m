m = 20;% must be an even number
n = m/2;
d = 1:m; a = 1:m-1; c = a; e = 1:m-2; f = e;x = ones(m,1);
D = zeros(2,2,n); A = zeros(2,2,n-1); C = zeros(2,2,n-1); X = zeros(2,n); B =  zeros(2,n);
M = diag(d) + diag(a,-1) + diag(c,1) + diag(e,-2) + diag(f,2);
b = M*x;



for i=1:n-1
    D(:,:,i) = [d(2*i-1),c(2*i-1); a(2*i-1),d(2*i)];
    A(:,:,i) = [e(2*i-1),a(2*i);0,e(2*i)];
    C(:,:,i) = [f(2*i-1),0;c(2*i),f(2*i)]; 
    B(:,i) = [b(2*i-1);b(2*i)];
    X(:,i) = [x(2*i-1);x(2*i)];
end
D(:,:,n) = [d(2*n-1),c(2*n-1); a(2*n-1),d(2*n)];
B(:,n) = [b(2*n-1);b(2*n)];
X(:,n) = [x(2*n-1);x(2*n)];
%since the A and C arrays are size n-1 need to initialize the nth term
%disp(M);
%disp(A);
%disp(C);

for k=2:n 
    aa = D(1,1,k-1); bb = D(1,2,k-1); cc = D(2,1,k-1); dd = D(2,2,k-1);
    P = A(:,:,k-1)*[dd, -bb; -cc, aa]/(aa*dd-cc*bb);
    D(:,:,k) = D(:,:,k) - P*C(:,:,k-1);
    B(:,k) = B(:,k)- P*B(:,k-1);
end

aa = D(1,1,n); bb = D(1,2,n); cc = D(2,1,n); dd = D(2,2,n);
X(:,n) = [dd, -bb; -cc, aa]*B(:,n)/(aa*dd-cc*bb);
%disp(X(:,n));

for k=n-1:-1:1
    aa = D(1,1,k); bb = D(1,2,k); cc = D(2,1,k); dd = D(2,2,k);
    X(:,k) = [dd, -bb; -cc, aa]*(B(:,k)- C(:,:,k)*X(:,k+1))/(aa*dd-cc*bb);
end

disp(X);


    
