clear; close all;
f1 = @(x) log(x);f2 = @(x) log((x+1)/(x-1));
x  =  1.1:1:100.1;
y1 = f1(x); y2 = f2(x);

%plot(x,y1);

plot(x,y2); 