S1 = 0; S2 = 0;
for k=1:8
    S1 = S1 + (-1)^(k+1)/k;
end

for k=1:2:8
    S2 = S2 +3^(-k)/(k);
end


retval1 = abs(log(2)- S1)
retval2 = abs(log(2)-2*S2)
