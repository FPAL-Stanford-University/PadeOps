program CoriolisBudget
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank
   use decomp_2d_io 
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4
   real(rkind), dimension(:,:,:), allocatable :: dudx, dudy, dudz
   real(rkind), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz
   real(rkind), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T,  Tfluct, P, kSGS, nSGS, Pfluct, wfluct, wE
   real(rkind), dimension(:,:,:), allocatable :: tau_11, tau_12, tau_13, tau_23, tau_22, tau_33
   real(rkind) :: dx, dy, dz, Rib = 0.05d0, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile, tempname, fname
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0, Ro, Fr, omega_2, omega_3, G_1, G_2
   real(rkind) :: latitude, frameAngle
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false., isTurbines = .true., noisy = .FALSE. 
   real(rkind), dimension(:,:), allocatable :: R11, R12, R13, R22, R23, R33, TT, T13, T23
   real(rkind), dimension(:,:,:), allocatable :: S11, S12, S13, S22, S23, S33
   real(rkind), dimension(:,:,:), allocatable :: S11m, S12m, S13m, S22m, S23m, S33m
   real(rkind), dimension(:,:), allocatable :: umn_set, Tmn_set, wT, uT, vT, vmn_set, Rij_mean, Tij_mean
   real(rkind), dimension(:,:), allocatable :: P_mke, D_mke, C_mke, transport_mke, transportSGS_mke 
   real(rkind), dimension(:,:), allocatable :: budget_MKE, eps, bou_tke, prod_tke, transport_tke, transportSGS_tke, transport_p_tke
   integer :: VizDump_Schedule = 0
   integer, dimension(:), allocatable :: timesteps
   real(rkind), dimension(:), allocatable :: times, buff1d_1, buff1d_2, buff1d_3, dudzM, dvdzM, buff1d_4
   real(rkind), dimension(:), allocatable :: turbine_u, turbine_v, turbine_w
   real(rkind), dimension(:,:), allocatable :: timewrite
   real(rkind), dimension(:,:,:), allocatable :: budget_MKE_time, R_time
   integer :: nt 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, NumericalSchemeVert, VizDump_Schedule, Ro, Fr, G_1, G_2, latitude, frameAngle, isTurbines, noisy
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx =     Lx/real(nx,rkind) 
   dy =     Ly/real(ny,rkind) 
   dz =     Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(v)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(wfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)
   call ops%allocate3DField(Tfluct)
   call ops%allocate3DField(P)
   call ops%allocate3DField(Pfluct)
   call ops%allocate3DField(kSGS)       
   call ops%allocate3DField(nSGS)
   call ops%allocate3Dfield(dudx)
   call ops%allocate3Dfield(dudy)
   call ops%allocate3Dfield(dudz)
   call ops%allocate3Dfield(dvdx)
   call ops%allocate3Dfield(dvdy)
   call ops%allocate3Dfield(dvdz)
   call ops%allocate3Dfield(dwdx)
   call ops%allocate3Dfield(dwdy)
   call ops%allocate3Dfield(dwdz)
   call ops%allocate3Dfield(tau_13)
   call ops%allocate3Dfield(tau_23)
   call ops%allocate3Dfield(S11)
   call ops%allocate3Dfield(S12)
   call ops%allocate3Dfield(S13)
   call ops%allocate3Dfield(S22)
   call ops%allocate3Dfield(S23)
   call ops%allocate3Dfield(S33)
   call ops%allocate3Dfield(S11m)
   call ops%allocate3Dfield(S12m)
   call ops%allocate3Dfield(S13m)
   call ops%allocate3Dfield(S22m)
   call ops%allocate3Dfield(S23m)
   call ops%allocate3Dfield(S33m)   
 
   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
  
   allocate(dudzM(nz)) 
   allocate(dvdzM(nz))
   allocate(buff1d_1(nz))
   allocate(buff1d_2(nz))
   allocate(buff1d_3(nz))
   allocate(buff1d_4(nz))
   allocate(turbine_u(nz))
   allocate(turbine_v(nz))
   allocate(turbine_w(nz))
   !allocate(wE(ops%gpE%xsz(1),ops%gpE%xsz(2),ops%gpE%xsz(3)))
   turbine_u = 0.
   turbine_v = 0.
   turbine_w = 0.

   if (VizDump_Schedule == 1) then
      call ops%Read_VizSummary(times, timesteps)
      nt = size(timesteps)
   else
      nt = (tstop - tstart)/tstep
   end if

   

   call message(0,"Number of snapshots to read:", nt)

   !if (nrank == 0) then
      ! Means
      allocate(umn_set(nz,nt))  
      allocate(Tmn_set(nz,nt))  
      allocate(vmn_set(nz,nt))
      ! Reynolds stresses
      allocate(R11(nz,nt))
      allocate(R12(nz,nt))
      allocate(R13(nz,nt))
      allocate(R22(nz,nt))
      allocate(R23(nz,nt))
      allocate(R33(nz,nt))
      allocate(T13(nz,nt))
      allocate(T23(nz,nt))
      allocate(TT(nz,nt))
      ! Velocity temperature corr
      allocate(uT(nz,nt))
      allocate(vT(nz,nt))
      allocate(wT(nz,nt))
      allocate(timewrite(nt,1))

      ! MKE
      allocate(P_mke(nz,nt))     
      allocate(D_mke(nz,nt))
      allocate(C_mke(nz,nt)) 
      allocate(transport_mke(nz,nt))
      allocate(transportSGS_mke(nz,nt))

      ! TKE
      allocate(transport_tke(nz,nt))
      allocate(transport_p_tke(nz,nt))
      allocate(transportSGS_tke(nz,nt))
      allocate(eps(nz,nt))
    
    
    !end if 
    allocate(budget_MKE(nz,16))
    allocate(budget_MKE_time(nz,nt,6))
    allocate(R_time(nz,nt,2))
    allocate(Rij_mean(nz,6))
    allocate(Tij_mean(nz,2))

   ! Initialize the turbine array
   if (isTurbines) then    
       call ops%create_turbine_array(inputfile) 
   end if

   idx = 1
   
   do while(idx <= nt)
      
      if (VizDump_Schedule == 1) then
         tidx = timesteps(idx)
      else
         tidx = tstart + tstep * (idx - 1)
      end if
      
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)
      call ops%ReadField3D(T,"potT",TIDX)
      !call ops%ReadField3D(kSGS,'kSGS',TIDX)
      call ops%ReadField3D(nSGS,'nSGS',TIDX)
      call ops%ReadField3D(P,'prss', TIDX)
      !call ops%ReadField3DEdge(wE, 'wVel', TIDX)
      !call message(0, "Read simulation data at time:", times(idx))
      ! Read edge velocity
      !write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",RunID, "_w.",TIDX
      !fname = ops%InputDir(:len_trim(ops%InputDir))//"/"//trim(tempname)
      !call decomp_2d_read_one(1, wE, fname, ops%gpE)


      ! STEP 2: Compute velocity gradients
      call ops%GetGradient(u, dudx, dudy, dudz, 0, 1)
      call ops%GetGradient(v, dvdx, dvdy, dvdz, 0, 1)
      call ops%GetGradient(w, dwdx, dwdy, dwdz, -1, -1)
      tau_11 = - 2. * nSGS * dudx
      tau_12 = - nSGS * (dudy + dvdx)
      tau_13 = - nSGS * (dudz + dwdx)
      tau_22 = - 2. * nSGS * dvdy
      tau_23 = - nSGS * (dvdz + dwdy)
      tau_33 = - 2. * nSGS * dwdz
      call ops%TakeMean_xy(tau_13, T13(:,idx)) !T13 
      call ops%TakeMean_xy(tau_23, T23(:,idx)) !T23
      S11m = 0.5 * (dudx + dudx)
      S12m = 0.5 * (dudy + dvdx)
      S13m = 0.5 * (dudz + dwdx)
      S22m = 0.5 * (dvdy + dvdy)
      S23m = 0.5 * (dvdz + dwdy)
      S33m = 0.5 * (dwdz + dwdz)    

 
      ! STEP 0: Compute means, mean gradients and fluctuations
      call ops%TakeMean_xy(u,umn_set(:,idx)) ! umn_set
      call ops%TakeMean_xy(T,Tmn_set(:,idx)) ! Tmn_set
      call ops%TakeMean_xy(v,vmn_set(:,idx)) ! vmn_set
      
      call ops%getFluct_from_MeanZ(v,vfluct)
      call ops%getFluct_from_MeanZ(u,ufluct)
      call ops%getFluct_from_MeanZ(w,wfluct)
      call ops%getFluct_from_MeanZ(T,Tfluct)
      call ops%getFluct_from_MeanZ(P,Pfluct)

      ! STEP 1: Compute correlations
      buff2 = ufluct*ufluct
      call ops%TakeMean_xy(buff2,R11(:,idx)) ! R11
      buff2 = ufluct*vfluct
      call ops%TakeMean_xy(buff2,R12(:,idx)) ! R12
      buff2 = ufluct*wfluct
      call ops%TakeMean_xy(buff2,R13(:,idx)) ! R13
      buff2 = vfluct*vfluct
      call ops%TakeMean_xy(buff2,R22(:,idx)) ! R22
      buff2 = vfluct*wfluct
      call ops%TakeMean_xy(buff2,R23(:,idx)) ! R23
      buff2 = wfluct*wfluct
      call ops%TakeMean_xy(buff2,R33(:,idx)) ! R33
      buff2 = Tfluct*Tfluct
      call ops%TakeMean_xy(buff2,TT(:,idx)) ! TT
      buff2 = ufluct*Tfluct
      call ops%TakeMean_xy(buff2,uT(:,idx)) ! uT
      buff2 = vfluct*Tfluct
      call ops%TakeMean_xy(buff2,vT(:,idx)) ! vT
      buff2 = wfluct*Tfluct
      call ops%TakeMean_xy(buff2,wT(:,idx)) ! wT

      ! STEP 2: Compute velocity gradients
      ! call ops%GetGradient(ufluct, dudx, dudy, dudz, -1, 0)
      ! call ops%GetGradient(vfluct, dvdx, dvdy, dvdz, -1, 0)
      ! call ops%GetGradient(wfluct, dwdx, dwdy, dwdz, -1, 0)
      call ops%ddz_1d(umn_set(:,idx), dudzM, 0, 1)
      call ops%ddz_1d(vmn_set(:,idx), dvdzM, 0, 1)

      ! STEP 3: Production/Destruction of MKE terms
      P_mke(:,idx) = R13(:,idx) * dudzM + R23(:,idx) * dvdzM ! P_mke
      D_mke(:,idx) = T13(:,idx) * dudzM + T23(:,idx) * dvdzM ! D_mke
 
      ! STEP 4: Forcing 
      omega_3 = sin(latitude*pi/180.d0)
      omega_2 = cos(latitude*pi/180.d0)*cos(frameAngle*pi/180.d0)
      C_mke(:,idx) = -(2./Ro) * omega_3 * umn_set(:,idx) * G_2 + (2./Ro) * omega_3 * vmn_set(:,idx) * G_1 ! C_mke

      ! STEP 5: MKE LHS
      call ops%ddz_1d(umn_set(:,idx) * R13(:,idx) + vmn_set(:,idx) * R23(:,idx), transport_mke(:,idx), 0, -1) ! transport_mke
      call ops%ddz_1d(umn_set(:,idx) * T13(:,idx) + vmn_set(:,idx) * T23(:,idx), transportSGS_mke(:,idx), 0, -1) ! transportSGS_mke
     
      budget_MKE_time(:,idx,1) = transport_mke(:,idx)
      budget_MKE_time(:,idx,2) = transportSGS_mke(:,idx)
      budget_MKE_time(:,idx,3) = P_mke(:,idx)
      budget_MKE_time(:,idx,4) = D_mke(:,idx)
      budget_MKE_time(:,idx,5) = C_mke(:,idx)    
 
      ! Get turbine forcing
      if (isTurbines == .true.) then
          call ops%get_turbine_RHS(u, v, w, buff1, buff2, buff3)
          call ops%TakeMean_xy(buff1, buff1d_1)
          turbine_u = turbine_u + buff1d_1
          budget_MKE_time(:,idx,6) = umn_set(:,idx) * buff1d_1
          call ops%TakeMean_xy(buff2, buff1d_1)
          turbine_v = turbine_v + buff1d_1
          budget_MKE_time(:,idx,6) = budget_MKE_time(:,idx,6) + vmn_set(:,idx) * buff1d_1
      end if

      ! Store diagonal Reynolds stress terms for all the time steps
      R_time(:,idx,1) = R12(:,idx)
      R_time(:,idx,2) = R13(:,idx)

      if (idx == 1) then 
          call ops%WriteField3D(buff1, "turb", idx)
          if (nrank == 0) then
          end if
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! TKE
      ! buff1 = 0.5 * (ufluct*ufluct + vfluct*vfluct + wfluct*wfluct)
      ! call ops%GetGradient(ufluct, dudx, dudy, dudz, 0, 1)
      ! call ops%GetGradient(vfluct, dvdx, dvdy, dvdz, 0, 1)
      ! call ops%GetGradient(wfluct, dwdx, dwdy, dwdz, -1, -1)
      ! tau_13 = - nSGS * (dudz + dwdx)
      ! tau_23 = - nSGS * (dvdz + dwdy)
      ! ! TKE LHS
      ! call ops%TakeMean_xy(wfluct*buff1, buff1d_1)
      ! call ops%ddz_1d(buff1d_1, transport_tke(:,idx), 0, 1) ! transport_tke
      ! call ops%TakeMean_xy(wfluct*Pfluct, buff1d_1)
      ! call ops%ddz_1d(buff1d_1, transport_p_tke(:,idx), 0, 0) ! transport_p_tke
      ! call ops%TakeMean_xy(ufluct*tau_13 + vfluct*tau_23, buff1d_1)
      ! call ops%ddz_1d(buff1d_1, transportSGS_tke(:,idx), 0, 0) ! transportSGS_tke

      ! ! TKE RHS
      ! ! Production
      ! call ops%TakeMean_xy(-ufluct*wfluct, buff1d_1)
      ! call ops%TakeMean_xy(-vfluct*wfluct, buff1d_2)
      ! prod_tke(:,idx) = buff1d_1 * dudzM + buff1d_2 * dvdzM ! prod_tke
      ! ! Bouancy
      ! call ops%TakeMean_xy(wfluct*Tfluct / Fr**2., bou_tke(:,idx)) ! bou_tke
      ! ! Dissipation
      ! S11 = 0.5 * (dudx + dudx) 
      ! S12 = 0.5 * (dudy + dvdx) 
      ! S13 = 0.5 * (dudz + dwdx) 
      ! S22 = 0.5 * (dvdy + dvdy) 
      ! S23 = 0.5 * (dvdz + dwdy) 
      ! S33 = 0.5 * (dwdz + dwdz)
      ! buff1 = -2. * nSGS * (S11*S11 + 2.*S12*S12 + 2.*S13*S13 + S22*S22 + 2.*S23*S23 + S33*S33) 
      ! buff2 = -2. * nSGS * (S11*S11m + 2.*S12*S12m + 2.*S13*S13m + S22*S22m + 2.*S23*S23m + S33*S33m)
      ! call ops%TakeMean_xy(buff1 + buff2, eps(:,idx)) ! eps
      
      idx = idx + 1
      call toc()
   end do 
   ! Destrio
   if (isTurbines == .true.) then
       call ops%destroy_turbine_array()
   end if

   ! STEP 1: Compute MKE > TKE, SGS
   call ops%ddz_1d(sum(umn_set,2)/real(idx,rkind),buff1d_1,0,1)
   budget_MKE(:,1) = (sum(R13,2)/real(idx,rkind))*buff1d_1
   budget_MKE(:,3) = (sum(T13,2)/real(idx,rkind))*buff1d_1
   
   call ops%ddz_1d(sum(vmn_set,2)/real(idx,rkind),buff1d_1,0,1)
   budget_MKE(:,2) = (sum(R23,2)/real(idx,rkind))*buff1d_1
   budget_MKE(:,4) = (sum(T23,2)/real(idx,rkind))*buff1d_1
    
   ! STEP 2: Compute Transfer 
   if (noisy == .TRUE.) then
       buff1d_1 = (sum(R13,2)/real(idx,rkind))*(sum(umn_set,2)/real(idx,rkind))
       call ops%ddz_1d(buff1d_1,budget_MKE(:,5),0,-1)
       buff1d_1 = (sum(R23,2)/real(idx,rkind))*(sum(vmn_set,2)/real(idx,rkind))
       call ops%ddz_1d(buff1d_1,budget_MKE(:,6),0,-1)
       buff1d_1 = (sum(T13,2)/real(idx,rkind))*(sum(umn_set,2)/real(idx,rkind))
       call ops%ddz_1d(buff1d_1,budget_MKE(:,7),0,-1)
       buff1d_1 = (sum(T23,2)/real(idx,rkind))*(sum(vmn_set,2)/real(idx,rkind))
       call ops%ddz_1d(buff1d_1,budget_MKE(:,8),0,-1)
   else 
       call ops%ddz_1d(sum(R13,2)/real(idx,rkind), buff1d_1, 0, -1) 
       budget_MKE(:,5) = budget_MKE(:, 1) + buff1d_1 * sum(umn_set,2)/real(idx,rkind)
       call ops%ddz_1d(sum(R23,2)/real(idx,rkind), buff1d_1, 0, -1) 
       budget_MKE(:,6) = budget_MKE(:, 2) + buff1d_1 * sum(vmn_set,2)/real(idx,rkind)
       call ops%ddz_1d(sum(T13,2)/real(idx,rkind), buff1d_1, 0, -1) 
       budget_MKE(:,7) = budget_MKE(:, 3) + buff1d_1 * sum(umn_set,2)/real(idx,rkind)
       call ops%ddz_1d(sum(T23,2)/real(idx,rkind), buff1d_1, 0, -1) 
       budget_MKE(:,8) = budget_MKE(:, 4) + buff1d_1 * sum(vmn_set,2)/real(idx,rkind)
   end if


   budget_MKE(:,9) = sum(C_mke, 2) / real(idx, rkind)
   budget_MKE(:,10) = sum(umn_set,2) / real(idx,rkind) * turbine_u / real(idx,rkind)
   budget_MKE(:,10) = budget_MKE(:,10) + sum(vmn_set,2) / real(idx,rkind) * turbine_v / real(idx,rkind)

   budget_MKE(:,11) = (sum(R13,2)/real(idx,rkind))*(sum(umn_set,2)/real(idx,rkind))
   budget_MKE(:,12) = (sum(R23,2)/real(idx,rkind))*(sum(vmn_set,2)/real(idx,rkind))
   budget_MKE(:,13) = sum(umn_set,2)/real(idx,rkind)
   budget_MKE(:,14) = sum(vmn_set,2)/real(idx,rkind)
   budget_MKE(:,15) = sum(Tmn_set,2)/real(idx,rkind)
   budget_MKE(:,16) = sum(wT,2)/real(idx,rkind)
       

   ! Reynolds stresses
   Rij_mean(:,1) = sum(R11,2) / real(idx, rkind) 
   Rij_mean(:,2) = sum(R12,2) / real(idx, rkind) 
   Rij_mean(:,3) = sum(R13,2) / real(idx, rkind) 
   Rij_mean(:,4) = sum(R22,2) / real(idx, rkind) 
   Rij_mean(:,5) = sum(R23,2) / real(idx, rkind) 
   Rij_mean(:,6) = sum(R33,2) / real(idx, rkind) 
   ! SGS stresses
   Tij_mean(:,1) = sum(T13,2) / real(idx, rkind) 
   Tij_mean(:,2) = sum(T23,2) / real(idx, rkind) 

   if (nrank == 0) then
      call ops%WriteASCII_2D(budget_MKE, "mkeB")
      call ops%WriteASCII_2D(budget_MKE_time(:,:,1), 'mke1')
      call ops%WriteASCII_2D(budget_MKE_time(:,:,2), 'mke2')
      call ops%WriteASCII_2D(budget_MKE_time(:,:,3), 'mke3')
      call ops%WriteASCII_2D(budget_MKE_time(:,:,4), 'mke4')
      call ops%WriteASCII_2D(budget_MKE_time(:,:,5), 'mke5')
      call ops%WriteASCII_2D(budget_MKE_time(:,:,6), 'mke6')
      call ops%WriteASCII_2D(Rij_mean, 'rijM')
      call ops%WriteASCII_2D(Tij_mean, 'tijM')
      call ops%WriteASCII_2D(umn_set, 'umnT')
      call ops%WriteASCII_2D(vmn_set, 'vmnT')
   end if 

   !if (nrank == 0) then
   !   call ops%WriteASCII_2D(transport_tke, "transport_tke")
   !   call ops%WriteASCII_2D(transport_p_tke, "transport_p_tke")  
   !   call ops%WriteASCII_2D(transportSGS_tke, "transportSGS_tke")   
   !   call ops%WriteASCII_2D(prod_tke, "prod_tke") 
   !   call ops%WriteASCII_2D(bou_tke, "bou_tke")
   !   call ops%WriteASCII_2D(eps, "eps")

   !   call ops%WriteASCII_2D(umn_set, "umn")
   !   call ops%WriteASCII_2D(Tmn_set, "Tmn")
   !   call ops%WriteASCII_2D(R11, "R11")
   !   call ops%WriteASCII_2D(R12, "R12")
   !   call ops%WriteASCII_2D(R13, "R13")
   !   call ops%WriteASCII_2D(R22, "R22")
   !   call ops%WriteASCII_2D(R23, "R23")
   !   call ops%WriteASCII_2D(R33, "R33")
   !   
   !   call ops%WriteASCII_2D(TT, "TTm")
   !   call ops%WriteASCII_2D(uT, "uTm")
   !   call ops%WriteASCII_2D(vT, "vTm")
   !   call ops%WriteASCII_2D(wT, "wTm")

   !   call ops%WriteASCII_2D(P_mke, "P_mke")
   !   call ops%WriteASCII_2D(D_mke, "D_mke")
   !   call ops%WriteASCII_2D(C_mke, "C_mke")
   !   call ops%WriteASCII_2D(transport_mke, "transport_mke")
   !   call ops%WriteASCII_2D(transportSGS_mke, "transportSGS_mke")

   !   timewrite(:,1) = times
   !   call ops%WriteASCII_2D(timewrite, "time")
   !end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

