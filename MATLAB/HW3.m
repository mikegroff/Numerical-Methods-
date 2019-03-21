J = @(x,y) [3*(x^2)-3*(y^2), -6*x*y; 6*x*y, 3*(x^2)-2*(y^2)];
F = @(x,y) [x^3-3*x*y^2-1; 3*x^2*y-y^3];

hold on;
for i =-1:2/100:1
    disp(i);
    for j=-1:2/100:1
        
        if(i == j)
            continue;
        end
        
        X = [i;j]; error = 1; tol = 0.5d-6; n = 0; X0 = X;
        while error > tol
            n = n+1;
            A = J(X(1),X(2)); b = F(X(1),X(2));
            H = NG(A,b);
            X = X - H;
            error = norm(X-X0);
            X0 = X;
        end
        
        
        if abs(X(1)-1) < tol % clearly z = 1 is a root
            scatter(i,j,'b');
            continue;
        end
        if X(2) < 0 % two complex root one a+bi and a-bi
            scatter(i,j,'g');
            continue;
        end
        if X(2) > 0 % the other one will have a b <0
            scatter(i,j,'r');
            continue;
        end
    end
end
hold off;
