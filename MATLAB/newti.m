function ppi = newti(x,xn,i)
% computes the i-th newton polynomial
ppi = ones(1,length(x)); if i == 1; return; end
for k=1:i-1
    ppi = ppi.*(x - xn(k));
end
