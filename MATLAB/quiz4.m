n=15;

a = zeros(n,n);
b = zeros(1,n);


for i=1:n
    for j=1:n
        if i>j
            a(i,j)= -1 +2*i;
        else
            a(i,j)= -1 +2*j;
        end      
    end
end

for k=1:n
    for l=1:n
    b(k) = b(k) + a(l,k);
    end
end

NaiveGE(a,b);
    
    

