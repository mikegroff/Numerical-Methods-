function [Ans] = NaiveGE(M,B)
A = zeros(length(B),1);

n = length(B);

for i=2:n
    m = M(i-1,i-1);

    for o=i:n
        r= M(o,i-1)/m;
        B(o) = B(o) - r*B(i-1);
        
        for j=1:n
            M(o,j) = M(o,j) - r*M(i-1,j);
        end
    end
    %disp(M);
  
end

%disp(M);
%disp(B);


k = n;



while (k>0)

    for l=k:n
        A(k) = A(k)+A(l)*M(k,l);
    end    
    A(k)= (B(k) - A(k))/M(k,k);
        k = k-1;
    
end

Ans = A;
disp(A);
end

