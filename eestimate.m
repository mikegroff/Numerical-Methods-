e = exp(1);
for k=1:8
    ee = (1+ 1/(8^k))^(8^k)
    error = abs(e-ee)
end
