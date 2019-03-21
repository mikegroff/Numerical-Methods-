x = linspace(-1,1); n=1000;
N = length(x); f = @(x) exp(x).*sin(10*x);
a = zeros(N,N); xn = zeros(N,1); yn = zeros(N,1); k=0;
while k<1
    w = input('Fist node ','s');  
    z = str2num(w); %#ok<*ST2NM>
    %to break out of loop when a good enough approximation is found
    if strcmp(w,'exit') == 1; return; end 
    % allows loop to reject values without continuing iteration
    if (abs(z) > 1); disp('Outside of domain'); continue; end  
    if strcmp(w,'') == 1; disp('Outside of domain'); continue; end  
    k = k+1;
end
xn(1) = z; yn(1) = f(xn(1));
a(1,1) = yn(1); p = a(1,1).*ones(N,1);
%disp(xn(1));
subplot(2,1,1);
plot(xn(1),yn(1),'o',x,p(1,:),'g',x,f(x));
legend('Data Points','Approximation','Function');
xlabel('x'); ylabel('y'); title('Newtons Interpolation');
subplot(2,1,2);
plot(x,abs(f(x) - p(1,:)),'r');
legend('Error');
xlabel('x'); ylabel('y'); title('Approximation Error');
i = 2;
while i <= N
    w = input('Next node ','s'); 
    z = str2num(w); 
    %to break out of loop when a good enough approximation is found
    if strcmp(w,'exit') == 1; return; end 
    % allows loop to reject values without continuing iteration
    if (abs(z) > 1); disp('Outside of domain'); continue; end  
    if strcmp(w,'') == 1; disp('Outside of domain'); continue; end  
    % disallows repeat values which cause errors
    if ismember(z,xn(1:i-1)); disp('Number used previously'); continue; end
    xn(i) = z; yn(i) = f(xn(i));
    a(i,1) = yn(i);
    for j = 1:i-1
        a(i-j,j+1) = ( a(i-j+1,j) - a(i-j,j) )/( xn(i) - xn(i-j));
    end

    K = newti(x,xn,i);
    p = p + a(1,i).*K(1,:);
    subplot(2,1,1);
    plot(xn(1:i),yn(1:i),'o',x,p(1,:),'g',x,f(x));
    legend('Data Points','Approximation','Function');
    xlabel('x'); ylabel('y'); title('Newtons Interpolation');
    subplot(2,1,2);
    plot(x,abs(f(x) - p(1,:)),'r');
    legend('Error');
    xlabel('x'); ylabel('y'); title('Approximation Error');
    i = i+1;
end



