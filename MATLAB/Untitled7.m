clf;
xn = -1:0.5:1; yn = sin(xn);
x = linspace(-1,1); x= x(:); l = zeros(length(x),length(xn)-1);
hold;
l(:,1) = (x -xn(1))/(xn(2) - xn(1));
v = {'k','b','r','g','c'};
p = zeros(length(x),1);
ppi(:,1) = ones(length(x),1);
plot(x,ppi(:,1),v(1));
for k = 2:5
    %l(:,k) = ones(length(x),1);
    ppi(:,k) = ones(length(x),1);
    for j=1:k-1
        l(:,k) = ppi(:,k).*(x-xn(j));
    end
    plot(x,ppi(:,k),v(k));
    p = p +yn(k)*l(:,k);
end
hold
%plot(xn,sin(xn),'k*',x(1:5:end),sin(x(1:5:end)),'ro',x,p)

