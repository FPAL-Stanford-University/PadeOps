function dudx = ddx_hit(u)
[nx,ny,nz] = size(u);
kx = fftshift(-nx/2:1:nx/2-1); ky = fftshift(-ny/2:1:ny/2-1); kz = fftshift(-nz/2:1:nz/2-1);
[k1,~,~] = ndgrid(kx,ky,kz);
k1(end/2+1,:,:) = 0;
dudx = ifft(1i*k1.*fft(u,[],1),[],1,'symmetric');