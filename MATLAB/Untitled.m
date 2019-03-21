p = .55;
T = 0;
for n=525:1000
    T = T + factorial(1000)*(p^i)*(1-p)^(1000-i)/(factorial(i)*factorial(1000-i));
end
disp(T);