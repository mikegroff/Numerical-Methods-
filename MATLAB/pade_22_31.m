x = linspace(-1+eps,1); f = log(1+x);
R22 = (x+(1/2)*x.^2)./(1-(1/6)*x - (1/36)*x.^2);
R31 = (x+x.^2/4 - x.^3/24)./(1+3*x./4);
subplot(2,1,1), plot(x,f, x,R22,'--',x,R31,':') 
xlabel('X'); ylabel('Y')
legend('ln(1+x)','R22' , 'R31');
subplot(2,1,2), plot(x , log(abs(f-R22)),x , log(abs(f-R31)),'--')
xlabel('X'); ylabel('Log(Y)')
legend('R22' , 'R31');