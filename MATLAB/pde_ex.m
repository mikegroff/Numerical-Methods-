L = pi; N = 100; dx = L/N; dt = 1d-1; x=dx:dx:pi-dx; x= x(:); u = sin(x);
r = dt/dx^2; B = diag((1-2*r)*ones(length(x),1))+ ...
    diag(ones(length(x)-1,1)*r,1)+ diag(ones(length(x)-1,1)*r,-1);
tfinal = 1; dtout = 0.1; nout = dtout/dt;
nsteps = tfinal/dtout; usave=[0;u;0]; uex = usave; t =0;
X=[0;x;pi]; %plot(X, usave, 'g'); 
hold
for k =1:nsteps
    for n = 1:nout
        t = t + dt; u = B*u;
    end
    usave = [usave [0;u;0]];
    uex = [uex [0;exp(-t)*sin(x);0]];
    plot(X,usave,X,uex,'ro');pause;
    %semilogy(X,abs(usave-uex));
    
end
hold;
