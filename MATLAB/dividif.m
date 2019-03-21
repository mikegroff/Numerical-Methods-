function a = dividif(xn,yn)
n = length(xn);
A = zeros(n);
A(:,1) = yn(:); 
for j=2:n
    for i=1:n-j+1
        A(i,j) = (A(i+1,j-1) -A(i,j-1))/(xn(j+i-1)-xn(i));               
    end
end
a = A(1,:);