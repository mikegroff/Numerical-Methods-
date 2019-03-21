f = @(x) 2*x^3 + 13*x^2 - 67*x + 30;
fp = @(x) 6*x^2 + 26*x - 67;
fpp = @(x) 12*x + 26;
fppp = @(x) 12;
tol = 0.5d-6; b = -30; c = 30;
%using the rational roots theorem we see that the smallest and largest possible zeros are -30 and 30
%bisection 
[x1,n1] = bisection(fpp,b,c,tol);
a = x1;
[x1,n2] = bisection(fp,b,a,tol);
[x2,n3] = bisection(fp,a,c,tol);
a = x1; b = x2; c = -30; d = 30;

[x1,n1] = bisection(f,c,a,tol);
[x2,n2] = bisection(f,a,b,tol);
[x3,n3] = bisection(f,b,d,tol);
Xb = [x1,n1;x2,n2;x3,n3];

%newtons
[x1,n] = newton(fpp,fppp,0,tol);
a = x1;
[x1,n] = newton(fp,fpp,a-50,tol);
[x2,n] = newton(fp,fpp,a+50,tol);
c = -30; d = 30;
[x1,n1] = newton(f,fp,c,tol);
[x2,n2] = newton(f,fp,a,tol);
[x3,n3] = newton(f,fp,d,tol);
Xn = [x1,n1;x2,n2;x3,n3];
%secant
c = -30; d = 30;
[x1,n] = secant(fpp,c,d,tol);
a = x1; c = -30; d = 30;
[x1,n] = secant(fp,c,a,tol);
[x2,n] = secant(fp,a,d,tol);
a = x1; b = x2; c = -30; d = 30;
[x1,n1] = secant(f,c,a,tol);
[x2,n2] = secant(f,a,b,tol);
[x3,n3] = secant(f,b,d,tol);
Xs = [x1,n1;x2,n2;x3,n3];

disp(Xb);
disp(Xn);
disp(Xs);

