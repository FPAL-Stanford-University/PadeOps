function dudy = ddy_hit(u)
[nx,ny,nz] = size(u);
kx = fftshift(-nx/2:1:nx/2-1); ky = fftshift(-ny/2:1:ny/2-1); kz = fftshift(-nz/2:1:nz/2-1);
[~,k2,~] = ndgrid(kx,ky,kz);
k2(:,end/2+1,:) = 0;
dudy = ifft(1i*k2.*fft(u,[],2),[],2,'symmetric');