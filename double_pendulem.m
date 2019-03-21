%double pendulum problem
clear; close all;
m1 = 1.1; m2 = 1; l1 = 1; l2 = 1; g = 9.81;
f1 = @(x) m2*l2*x(4)^2*sin(x(1)-x(2)) - (m1+m2)*g*sin(x(1));
f2 = @(x) -l1*x(3)^2*sin(x(1)-x(2)) - g*sin(x(2));
F = @(t,x) [ x(3);x(4);
    (12*f1(x)-m2*l2*cos(x(1)-x(2))*f2(x))/( (m1+m2)*l1*l2 - m2*l1*l2*cos(x(1)-x(2))^2 );
    ( -l1*cos(x(1)-x(2))*f1(x) + (m1+m2)*l1*f2(x))/( (m1+m2)*l1*l2 - m2*l1*l2*cos(x(1)-x(2))^2 )];

x(1) = 0 ;x(2) = pi; x(3) = 0.001; x(4) = 0; h = 0.01; t = 0;
X = zeros(4,1/h+1);
X(:,1) = [x(1);x(2);x(3);x(4)];
figure(1);
 plot([0 l1*sin(X(1,1))],[0 -l1*cos(X(1,1))],'o',...
     [l1*sin(X(1,1)), l1*sin(X(1,1))+l2*sin(X(2,1))],[-l1*cos(X(1,1)),-l1*cos(X(1,1))-cos(X(2,1))],'o',...
     'MarkerSize',5,'LineWidth',3);
 hold; %circles(l1*sin(X(1,1)),-l1*cos(X(1,1)),0.075,'color','b');
 %circles(l1*sin(X(1,1)) + l2*sin(X(2,1)), -l1*cos(X(1,1))-l2*cos(X(2,1)),0.075,'color','r');hold
 title('Initial Configuration')
 axis([-2 2 -2 2]); grid; axis equal; pause(1) 
 
for k = 1:10/h
    X(:,k+1) = rk4_step(F,t(k),X(:,k),h);
    t(k+1) = t(k) + h;
    plot([0 l1*sin(X(1,k+1))],[0 -l1*cos(X(1,k+1))],'o',[l1*sin(X(1,k+1)), l1*sin(X(1,k+1)) + l2*sin(X(2,k+1))],...
       [-l1*cos(X(1,k+1)), -l1*cos(X(1,k+1))-l2*cos(X(2,k+1))],'o','MarkerSize',5,'LineWidth',3); 
   hold; %circles(l1*sin(X(1,k+1)),-l1*cos(X(1,k+1)),0.075,'color','b');
% circles(l1*sin(X(1,k+1)) + l2*sin(X(2,k+1)), -l1*cos(X(1,k+1))-l2*cos(X(2,k+1)),0.075,'color','r');hold
    axis([-(l1+l2), l1+l2, -(l1+l2), l1+l2]); grid; axis equal;
    title(sprintf('Evolution at t = %g', t(k+1)));
    pause(0.005)
end

figure(2); plot(t,X(1,:)); hold; plot(t,X(2,:)); hold
