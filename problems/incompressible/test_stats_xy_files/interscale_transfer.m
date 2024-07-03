clc, clearvars
addpath('/home1/06632/ryanhass/MATLAB/MATLABpostProcessing-PhD/matlabFunctions')
addpath('/home1/06632/ryanhass/MATLAB/MATLABpostProcessing-PhD/problems/HITdecay/shearlessMixing_Bodart/SMcommon')
echo on

tid = 10000;
runID = 78;

kc = 1.5;

nx = 128;
ny = 128;
nz = 128;

Lx = 6;
Ly = 6;
zmin = -4;
zmax = 2;

dx = Lx/nx;
dy = Ly/ny;
dz = (zmax-zmin)/nz;

inputdir = '/home1/06632/ryanhass/codes/PadeOps/problems/incompressible/test_stats_xy_files/data/'; 

[u,v,w] = readInVelocityData(inputdir,nx,ny,nz,runID,tid);
T    = read_fortran_box([inputdir,'Run',sprintf('%02d',runID),'_potT_t',sprintf('%06d',tid),'.out'],...
  nx, ny, nz, 'double');

U = squeeze(mean(mean(u,1),2));
V = squeeze(mean(mean(v,1),2));
W = squeeze(mean(mean(w,1),2));

u = u - mean(mean(u,1),2);
v = v - mean(mean(v,1),2);
w = w - mean(mean(w,1),2);
T = T - mean(mean(T,1),2);

duidxj = zeros(nx,ny,nz,3,3);
dTdxj  = zeros(nx,ny,nz,3);
[duidxj(:,:,:,1,1),duidxj(:,:,:,1,2),duidxj(:,:,:,1,3)] = getGradient(u,Lx,Ly,dz);
[duidxj(:,:,:,2,1),duidxj(:,:,:,2,2),duidxj(:,:,:,2,3)] = getGradient(v,Lx,Ly,dz);
[duidxj(:,:,:,3,1),duidxj(:,:,:,3,2),duidxj(:,:,:,3,3)] = getGradient(w,Lx,Ly,dz);
[dTdxj(:,:,:,1),dTdxj(:,:,:,2),dTdxj(:,:,:,3)] = getGradient(T,Lx,Ly,dz);

tau11 = zeros(nx,ny,nz); 
tau12 = zeros(nx,ny,nz); 
tau13 = zeros(nx,ny,nz); 
tau22 = zeros(nx,ny,nz); 
tau23 = zeros(nx,ny,nz); 
tau33 = zeros(nx,ny,nz); 
fx    = zeros(nx,ny,nz);
fy    = zeros(nx,ny,nz);
fz    = zeros(nx,ny,nz);

% Scale splitting
[q,~,~,duidxj] = scale_splitting(u,v,w,T,duidxj,...
  tau11,tau12,tau13,tau22,tau23,tau33,fx,fy,fz,kc,Lx,Ly);
dTdxj_ss = zeros(nx,ny,nz,3,2);
dTdxj_ss(:,:,:,1,1) = lowpass_filter(dTdxj(:,:,:,1),kc,Lx,Ly);
dTdxj_ss(:,:,:,2,1) = lowpass_filter(dTdxj(:,:,:,2),kc,Lx,Ly);
dTdxj_ss(:,:,:,3,1) = lowpass_filter(dTdxj(:,:,:,3),kc,Lx,Ly);
dTdxj_ss(:,:,:,1,2) = dTdxj(:,:,:,1) - dTdxj_ss(:,:,:,1,1);
dTdxj_ss(:,:,:,2,2) = dTdxj(:,:,:,2) - dTdxj_ss(:,:,:,2,1);
dTdxj_ss(:,:,:,3,2) = dTdxj(:,:,:,3) - dTdxj_ss(:,:,:,3,1);

% Do scale 1 ("r") first
for scale = 1:2

    % Turbulent transport
    statsX = -ddzCD06(...
      squeeze(mean(mean(w.*q(:,:,:,1,scale).^2,1),2)),...
      dz,-1,-1);
    statsY = -ddzCD06(...
      squeeze(mean(mean(w.*q(:,:,:,2,scale).^2,1),2)),...
      dz,-1,-1);
    statsZ = -ddzCD06(...
      squeeze(mean(mean(w.*q(:,:,:,3,scale).^2,1),2)),...
      dz,-1,-1);
    statsT = -ddzCD06(...
      squeeze(mean(mean(w.*q(:,:,:,4,scale).^2,1),2)),...
      dz,-1,-1);
    statswT = -ddzCD06(...
      squeeze(mean(mean(w.*q(:,:,:,4,scale).*q(:,:,:,3,scale),1),2)),...
      dz,-1,-1);

    save_1Dascii(statsX,['turb_trans_X_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statsY,['turb_trans_Y_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statsZ,['turb_trans_Z_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statsT,['turb_trans_T_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statswT,['turb_trans_wT_scl',sprintf('%02d',scale),'.dat'])

    % Turbulent transport mixed scale
    statsX = -2*ddzCD06(...
      squeeze(mean(mean(q(:,:,:,1,1).*q(:,:,:,1,2).*q(:,:,:,3,1+mod(scale,2)),1),2)),...
      dz,-1,-1);
    statsY = -2*ddzCD06(...
      squeeze(mean(mean(q(:,:,:,2,1).*q(:,:,:,2,2).*q(:,:,:,3,1+mod(scale,2)),1),2)),...
      dz,-1,-1);
    statsZ = -2*ddzCD06(...
      squeeze(mean(mean(q(:,:,:,3,1).*q(:,:,:,3,2).*q(:,:,:,3,1+mod(scale,2)),1),2)),...
      dz,-1,-1);
    statsT = -2*ddzCD06(...
      squeeze(mean(mean(q(:,:,:,4,1).*q(:,:,:,4,2).*q(:,:,:,3,1+mod(scale,2)),1),2)),...
      dz,-1,-1);
    statswT = -ddzCD06(...
      squeeze(mean(mean(q(:,:,:,3,1).*q(:,:,:,4,2).*q(:,:,:,3,1+mod(scale,2)),1),2)) + ...
      squeeze(mean(mean(q(:,:,:,3,2).*q(:,:,:,4,1).*q(:,:,:,3,1+mod(scale,2)),1),2)),...
      dz,-1,-1);

    save_1Dascii(statsX,['turb_trans_mixed_X_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statsY,['turb_trans_mixed_Y_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statsZ,['turb_trans_mixed_Z_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statsT,['turb_trans_mixed_T_scl',sprintf('%02d',scale),'.dat'])
    save_1Dascii(statswT,['turb_trans_mixed_wT_scl',sprintf('%02d',scale),'.dat'])

end
% Interscale transfer
statsX = squeeze(mean(mean(...
  q(:,:,:,1,1).*(q(:,:,:,1,1).*duidxj(:,:,:,1,1,2) + q(:,:,:,2,1).*duidxj(:,:,:,1,2,2) + q(:,:,:,3,1).*duidxj(:,:,:,1,3,2)),1),2)) - ...
  squeeze(mean(mean(...
  q(:,:,:,1,2).*(q(:,:,:,1,2).*duidxj(:,:,:,1,1,1) + q(:,:,:,2,2).*duidxj(:,:,:,1,2,1) + q(:,:,:,3,2).*duidxj(:,:,:,1,3,1)),1),2));
  
statsY = squeeze(mean(mean(...
  q(:,:,:,2,1).*(q(:,:,:,1,1).*duidxj(:,:,:,2,1,2) + q(:,:,:,2,1).*duidxj(:,:,:,2,2,2) + q(:,:,:,3,1).*duidxj(:,:,:,2,3,2)),1),2)) - ...
  squeeze(mean(mean(...
  q(:,:,:,2,2).*(q(:,:,:,1,2).*duidxj(:,:,:,2,1,1) + q(:,:,:,2,2).*duidxj(:,:,:,2,2,1) + q(:,:,:,3,2).*duidxj(:,:,:,2,3,1)),1),2));
  
statsZ = squeeze(mean(mean(...
  q(:,:,:,3,1).*(q(:,:,:,1,1).*duidxj(:,:,:,3,1,2) + q(:,:,:,2,1).*duidxj(:,:,:,3,2,2) + q(:,:,:,3,1).*duidxj(:,:,:,3,3,2)),1),2)) - ...
  squeeze(mean(mean(...
  q(:,:,:,3,2).*(q(:,:,:,1,2).*duidxj(:,:,:,3,1,1) + q(:,:,:,2,2).*duidxj(:,:,:,3,2,1) + q(:,:,:,3,2).*duidxj(:,:,:,3,3,1)),1),2));
  
statsT = squeeze(mean(mean(...
  q(:,:,:,4,1).*(q(:,:,:,1,1).*dTdxj_ss(:,:,:,1,2) + q(:,:,:,2,1).*dTdxj_ss(:,:,:,2,2) + q(:,:,:,3,1).*dTdxj_ss(:,:,:,3,2)),1),2)) - ...
  squeeze(mean(mean(...
  q(:,:,:,4,2).*(q(:,:,:,1,2).*dTdxj_ss(:,:,:,1,1) + q(:,:,:,2,2).*dTdxj_ss(:,:,:,2,1) + q(:,:,:,3,2).*dTdxj_ss(:,:,:,3,1)),1),2));
  
statswT = squeeze(mean(mean(...
  q(:,:,:,4,1).*(q(:,:,:,1,1).*duidxj(:,:,:,3,1,2) + q(:,:,:,2,1).*duidxj(:,:,:,3,2,2) + q(:,:,:,3,1).*duidxj(:,:,:,3,3,2)),1),2)) - ...
  squeeze(mean(mean(...
  q(:,:,:,4,2).*(q(:,:,:,1,2).*duidxj(:,:,:,3,3,1) + q(:,:,:,2,2).*duidxj(:,:,:,3,2,1) + q(:,:,:,3,2).*duidxj(:,:,:,3,3,1)),1),2)) + ...
  squeeze(mean(mean(...
  q(:,:,:,3,1).*(q(:,:,:,1,1).*dTdxj_ss(:,:,:,1,2) + q(:,:,:,2,1).*dTdxj_ss(:,:,:,2,2) + q(:,:,:,3,1).*dTdxj_ss(:,:,:,3,2)),1),2)) - ...
  squeeze(mean(mean(...
  q(:,:,:,3,2).*(q(:,:,:,1,2).*dTdxj_ss(:,:,:,3,1) + q(:,:,:,2,2).*dTdxj_ss(:,:,:,2,1) + q(:,:,:,3,2).*dTdxj_ss(:,:,:,3,1)),1),2)); 
  
save_1Dascii(statsX,['interscale_X.dat'])
save_1Dascii(statsY,['interscale_Y.dat'])
save_1Dascii(statsZ,['interscale_Z.dat'])
save_1Dascii(statsT,['interscale_T.dat'])
save_1Dascii(statswT,['interscale_wT.dat'])
