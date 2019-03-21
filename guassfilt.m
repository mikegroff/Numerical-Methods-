N =13;
%%horizontal and vertical coordinates
coordinates = -floor(N/2) : floor(N/2);
[X Y] = meshgrid(coordinates, coordinates);
sigma=2;

gfilter = exp(-(X.^2 + Y.^2)/(2*sigma.^2));
gfilter = gfilter/sum(gfilter(:))

gridInv = (-1).^(X+Y)
hfilter = gfilter.*gridInv;

subplot(2,2,1);
mesh(gfilter);
subplot(2,2,2);
mesh(hfilter);

subplot(2,2,3);
imagesc(abs(fftshift(fft2(gfilter))));
subplot(2,2,4);
imagesc(abs(fftshift(fft2(hfilter))));