function dudz = ddz_hit(u)
[nx,ny,nz] = size(u);
kx = fftshift(-nx/2:1:nx/2-1); ky = fftshift(-ny/2:1:ny/2-1); kz = fftshift(-nz/2:1:nz/2-1);
[~,~,k3] = ndgrid(kx,ky,kz);
k3(:,:,end/2+1) = 0;
dudz = ifft(1i*k3.*fft(u,[],3),[],3,'symmetric');