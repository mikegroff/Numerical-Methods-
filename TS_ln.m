clear; close all;
f1 = @(x) x; f2 = @(x) x-x.^2/2 ; f3 = @(x) f2(x) + x.^3/3;  f4 = @(x) f3(x) - x.^4/4; f5 = @(x) f4(x) + x.^5/5;
x= 0:1/128:3; y = log(1+x); y1 = f1(x);  y2 = f2(x); y3 = f3(x); y4= f4(x); y5 = f5(x);

subplot(2,1,1), plot(x,y,x,y1,x,y2, x,y3 ,x,y4,x,y5 )

legend('ln(1+x)', 'S_1', 'S_2','S_3','S_4','S_5')
title('Ln(1+x) and its first five Taylor approximations')
xlabel('x'); ylabel('y');
subplot(2,1,2),semilogy(x,abs(y-y1),x,abs(y-y2),x,abs(y-y3),x,abs(y-y4),x,abs(y-y5))
legend('|y-y_1|','|y-y_2|','|y-y_3|','|y-y_4|','|y-y_5|')
title('Errors')
xlabel('x'); ylabel('ln(y)');