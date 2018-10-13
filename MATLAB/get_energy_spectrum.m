function [E,kline] = get_energy_spectrum(uMesh,Nspect)
[nx,ny,nz] = size(uMesh);
uHat = fftn(uMesh);
k_ind = fftshift(-nx/2 : nx/2 - 1);
[k1,k2,k3] = ndgrid(k_ind);
uHat = uHat(:)/(nx*ny*nz);
k_abs_use = sqrt(k1(:).^2 + k2(:).^2 + k3(:).^2);
[k_abs_use,idx] = sort(k_abs_use);
uHat = uHat(idx);
kline = linspace(min(k_abs_use),max(k_abs_use),Nspect);
[N,BIN] = histc(k_abs_use,kline);
E = zeros(length(N),1);
for i = 1:nx*ny*nz
    E(BIN(i)) = E(BIN(i)) + (abs(uHat(i)).^2);
end
E = E/(kline(2)-kline(1));