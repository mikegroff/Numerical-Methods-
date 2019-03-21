function [] = sequencelarge(a)
for n=2:1000
    b = n*a;
    a = b/n;
end
disp(a);
disp(b);

end