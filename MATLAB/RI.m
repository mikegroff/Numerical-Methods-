function D = RI(phi,h,n)

D(1,1) = phi(h);
if(n==1); return; end;

for i = 2:n
    D(i,1) = phi(h/(2^(i-1)));
    for j = 2:i
        D(i+1-j,j) = ( 4^(j-1)*D(i+2-j,j-1)- D(i+1-j,j-1))/(4^(j-1)-1);
    end
end

        
    