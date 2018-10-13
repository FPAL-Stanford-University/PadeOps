function ffil = sharp_filter(f,kco)
n = size(f,1);
k_ind = fftshift(-n/2 : n/2 - 1);
[k1,k2,k3] = ndgrid(k_ind);
kabs_sq = k1.*k1 + k2.*k2 + k3.*k3;
fhat = fftn(f);
kco_sq = kco^2;
fhat(kabs_sq >= kco_sq) = 0;
ffil = real(ifftn(fhat));