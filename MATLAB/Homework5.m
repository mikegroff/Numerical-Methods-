f = @(x) cos(2.*x)./(exp(x));
x = linspace(0,5);
a = 0; b = 5*pi/4; epsi = .0005;
levl = 0; mxlevl = 4;

[n,C,F] = Simpson(f,a,b,epsi,levl,mxlevl);
disp(n)
disp(C)
disp(F)

Z = zeros(1,length(F));
P1 = [C;C]; P2 = [Z;F];

plot(x,f(x),'r',C,F,'ob',P1(:,:),P2(:,:),'k',x,zeros(100,1),'k')
title('Simpson''s Rule of a Function');
legend('f(x)=cos(2x)/e^x','Partition Points');
xlabel('x'); ylabel('y');

function [n,C,F] = Simpson(f,a,b,epsi,levl,mxlevl)
  levl = 1 + levl;
  h = b-a;
  c = (a+b)/2;
  d = (a+c)/2;
  e = (c+b)/2;
  C=c;
  F=f(c);
  simpsonone = h*(f(a) + 4*f(c) + f(b))/6;
  simpsontwo = h*(f(a) + 4*f(d) + 2*f(c) + 4*f(e) + f(b))/12;
  
  if levl>=mxlevl
      n = simpsontwo;
  elseif abs(simpsontwo - simpsonone) < (15*epsi)
      n = simpsontwo + (simpsontwo - simpsonone)/15;
  else
      [lsimpson,C1,F1] = Simpson(f,a,c,epsi/2,levl,mxlevl);
      [rsimpson,C2,F2] = Simpson(f,c,b,epsi/2,levl,mxlevl);
      C = [C,C1,C2];
      F = [F,F1,F2];
      n = lsimpson + rsimpson;
  end
 end
      
    
  