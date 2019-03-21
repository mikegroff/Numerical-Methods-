%Differentiation of two discrete functions at the same time
N = 24; h = 2*pi/N; x = h*(1:N)';
vh = max(0,1-abs(x-pi)/2); vp= exp(sin(x));
v = [vh'; vp'];
v_hat = fft(v');
w_hat = 1i*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
w = real(ifft(w_hat)); w = w'; clf
subplot(3,2,1), plot(x,v(1,:),'.-','markersize',13)
axis([0 2*pi -.5 1.5]), grid on, title('function')
subplot(3,2,2), plot(x,w(1,:),'.-','markersize',13)
axis([0 2*pi -1 1]), grid on, title('spectral derivative')

subplot(3,2,3), plot(x,v(2,:),'.-','markersize',13)
axis([0 2*pi 0 3]), grid on
subplot(3,2,4), plot(x,w(2,:),'.-','markersize',13)
axis([0 2*pi -2 2]), grid on