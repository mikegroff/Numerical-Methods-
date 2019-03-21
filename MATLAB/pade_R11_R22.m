x = linspace(-1,1); f = exp(x);
R11 = (2+x)./(2-x);
R22 = (2+x+(1/6)*x.^2)./(2-x+(1/6)*x.^2);
subplot(2,1,1);
plot(x,f,x, R11, '--',x , R22,':');
xlabel('X'); ylabel('Y')
legend('exp(x)','R11' , 'R22');
subplot(2,1,2),plot(x, log(abs(f-R11)), x , log(abs(f-R22)),'--');
xlabel('X'); ylabel('Log(Y)')
legend('R11' , 'R22');      