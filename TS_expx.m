clear; close all;
f1 = @(x) 1; f2 = @(x) 1+x ; f3 = @(x) f2(x) + x.^2/2;
x= -1:1/128:1; y = exp(x); y1 = f1(x);  y2 = f2(x); y3 = f3(x);

subplot(2,1,1), plot(x,y,x,y1,x,y2, x,y3 )

legend('e^x', 'S_1', 'S_2','S_3')
title('e^x) and its first three Taylor approximations')
xlabel('x'); ylabel('y');
subplot(2,1,2),semilogy(x,abs(y-y1),x,abs(y-y2),x,abs(y-y3))
legend('|y-y_1|','|y-y_2|','|y-y_3|')
title('Errors')
xlabel('x'); ylabel('ln(y)');