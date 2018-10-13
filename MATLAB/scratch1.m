clear; clc;
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');

for n=1:7
u = read_fortran_box(['Run01_uVel_t00' num2str(n) '000.out'], 256, 256, 256, 'double');
v = read_fortran_box(['Run01_vVel_t00' num2str(n) '000.out'], 256, 256, 256, 'double');
w = read_fortran_box(['Run01_wVel_t00' num2str(n) '000.out'], 256, 256, 256, 'double');
fid =fopen(['Run01_info_t00' num2str(n) '000.out']);
t=fscanf(fid,'%f'); t=t(1);

% Plot u
figure, surface(real(u(:,:,45))','edgecolor','none'), daspect([1 1 1]);
title(['t=' num2str(t)]);
saveas(gcf,[num2str(n) '000u.pdf']);
close;

% Plot energy spectrum
title(['t=' num2str(t)]);
[Eu,~] = get_energy_spectrum(u,500);
[Ev,~] = get_energy_spectrum(v,500);
[Ew,kline] = get_energy_spectrum(w,500);
figure;
loglog(kline, 3/2*Eu);
hold on;
loglog(kline, 1.4528*(kline.^4)./((1 + kline.^2).^(17/6)).*exp(-5.2.*kline.*0.0245));
hold off;
saveas(gcf,[num2str(n) '000E.pdf']);
close;

% Energy dissipation rate
dudx=ddx_hit(u); dudy=ddy_hit(u); dudz=ddz_hit(u);
dvdx=ddx_hit(v); dvdy=ddy_hit(v); dvdz=ddz_hit(v);
dwdx=ddx_hit(w); dwdy=ddy_hit(w); dwdz=ddz_hit(w);

epsilon=dudx.*dudx+dudy.*dudy+dudz.*dudz...
    +dvdx.*dvdx+dvdy.*dvdy+dvdz.*dvdz...
    +dwdx.*dwdx+dwdy.*dwdy+dwdz.*dwdz;
    
epsilon=1/34*mean(epsilon(:));

% Save
fid1=fopen('result', 'w');
fprintf(fid1, [num2str(t) ' ' num2str(epsilon)]);

end