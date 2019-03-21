e=exp(1.0);
for k=1:10
    ee = (1+1/(8^k))^(8^k);
    error = abs(e-ee);
    disp(e);
    disp(error);
end
