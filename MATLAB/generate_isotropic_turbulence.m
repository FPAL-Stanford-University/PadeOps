function [u,v,w] = generate_isotropic_turbulence(nx,ny,nz,seed)

% inputs: 
% 1. N - number of grid points (NxNxN), 
% 2. Integer Seed (for deterministic generation)

C_vk = 1.4528;
i0 = nx/2+1 ;
j0 = ny/2+1 ;
k0 = nz/2+1 ;

Mx = nx/2-1 ;
My = ny/2-1 ;
Mz = nz/2-1 ;

u = zeros(nx,ny,nz) ; v = u ; w = u ;

rng(seed)

for k1= 0:Mx
    for k2=-My:My
        for k3=-Mz:Mz
            k = sqrt(k1^2+k2^2+k3^2) ;
            if k > 0
                E = C_vk*(k.^4)./((1 + k.^2).^(17/6)).*exp(-5.2.*k.*0.0245);
                phi = rand(1,3)*2*pi ;
                a = sqrt( E/(2*pi*k^2) ) * exp(1i*phi(1)) * cos(phi(3)) ;
                b = sqrt( E/(2*pi*k^2) ) * exp(1i*phi(2)) * sin(phi(3)) ;
                k12 = sqrt(k1^2+k2^2) ;
                if k12 == 0
                    uhs = a ;
                    vhs = b ;
                    whs = 0 ;
                else
                    uhs = ( a*k*k2 + b*k1*k3 ) / ( k*k12 ) ;
                    vhs = ( b*k2*k3 - a*k*k1 ) / ( k*k12 ) ;
                    whs = ( - b*k12 ) / k ;
                end
                u(i0+k1,j0+k2,k0+k3) = uhs ;
                v(i0+k1,j0+k2,k0+k3) = vhs ;
                w(i0+k1,j0+k2,k0+k3) = whs ;
                u(i0-k1,j0-k2,k0-k3) = conj( uhs ) ;
                v(i0-k1,j0-k2,k0-k3) = conj( vhs ) ;
                w(i0-k1,j0-k2,k0-k3) = conj( whs ) ;
            end
        end
    end
end
u = fftn(fftshift(u)) ; v = fftn(fftshift(v)) ; w = fftn(fftshift(w))  ;
%u = real(u) ; v = real(v) ; w = real(w) ;

