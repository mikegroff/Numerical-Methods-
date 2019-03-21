x = linspace(0,2*pi); g = @(x) cos(2*x)./(exp(x)); 
plot(x,g(x),'o');