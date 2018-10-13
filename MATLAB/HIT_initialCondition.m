clear, clc, close all

seed = 12345;
nx = 256; ny = 256; nz = 256;
[u,v,w] = generate_isotropic_turbulence(nx,ny,nz,seed);

dat = zeros(nx*ny*nz,3);
dat(:,1) = real(u(:));
dat(:,2) = real(v(:));
dat(:,3) = real(w(:));

figure, surface(real(u(:,:,45))','edgecolor','none'), daspect([1 1 1])

save('PadeOps_HIT.dat','dat','-ascii')