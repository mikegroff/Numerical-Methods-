A = zeros(4,4); y = linspace(0,10); v = 1/20:1/20:10;
z = zeros(1,length(v)); N = length(v);
f = @(x) sin(x);
for k=1:97
    xn = [k/10, (k+1)/10, (k+2)/10, (k+3)/10];
    yn = [f(k/10), f((k+1)/10), f((k+2)/10), f((k+3)/10)];
    n = length(xn)*2;
    A(:,1) = yn;
    A(1,1) = yn(1);
    p = A(1,1).*ones(N,1);
    
    for i = 2:4
        for j = 1:i-1
            A(i-j,j+1) = ( A(i-j+1,j) - A(i-j,j) )/( xn(i) - xn(i-j));
        end
    end
    

    
    K = newti(v,xn,i);
    p = p + A(1,i).*K(1,:);
    
    z(2*k) = p(1,2*k);
    z(2*k+1) = p(1,2*k+1);
    
    if(k==97)
        for l=195:200
            z(l) = p(1,l);
        end
    end
    
end

    
%subplot(2,1,1);
%plot(v,f(v),v,z);
%subplot(2,1,2);
plot(v,abs(f(v)-z)); 



