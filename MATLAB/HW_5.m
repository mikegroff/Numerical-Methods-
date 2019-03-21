import containers.*
s = 0; f = 2*pi; x = linspace(s,f); g = @(x) cos(2.*x)./(exp(x)); emp = [];
hold on;
[A, P, fP] = ASimpson(g,s,f,5.0d-5,0,5,emp,emp);
disp(A);
plot(x,g(x),'g');
plot(x,zeros(length(x)));

for i =1:length(P)
    j = P(i);
    if fP(j) > 0
        plot([P(i).',P(i).'],[0.',fP(j).'],'b');  
    else
        plot([P(i).',P(i).'],[fP(j).',0.'],'b');
    end
    plot(P(i),fP(j),'or');
end

hold off;




%H = containers.Map((-L/2:h:L/2),t(U(lvl,:))); 
        %M(lvl) = ASimpson(H,-L/2,L/2,tol, 0, 6, emp,emp);