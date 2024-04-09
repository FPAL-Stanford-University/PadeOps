module spectralForcingLayerMod
  ! This module implements a forcing layer similar to Briggs et al. (1996), but with some modifications
    use kind_parameters, only: rkind, clen, mpirkind
    use reductions,      only: p_maxval, p_sum
    use gridtools,       only: onThisRank, getStEndIndices
    use constants,       only: im0, pi
    use spectralMod,     only: spectral
    use PadePoissonMod,  only: Padepoisson 
    use random,          only: uniform_random
    use exits,           only: message
    use decomp_2d
    use decomp_2d_io
    use PadeDerOps,      only: Pade6Stagg
    use fortran_assert,  only: assert
    use fringeMethod,    only: S_fringe 
    use mpi
    use seedGen,         only: initializeLogMap, incrementLogisticMap
    use sgsmod_igrid,    only: sgs_igrid  
    !use scalar_igridMod, only: scalar_igrid
    implicit none

    integer, parameter :: zbc_bot = 0, zbc_top = 0
    integer, parameter :: gaussian = 1, fringeFunction = 2

    type :: spectForcingLayer

      real(rkind), dimension(:,:,:), allocatable :: rmaskC, cmaskC, cmaskE, jC, jE
      real(rkind), dimension(:,:,:), allocatable :: integralMask
      real(rkind), dimension(:),     allocatable :: integralMask_1D
      real(rkind) :: seedFact, dV, tgtKE, tgtKEon3, ampFact_x, ampFact_y, &
        ampFact_z, ampFact_T, tgtDissipation, integralTime, tgtMPE, ampFact_time 
      integer :: seed
      real(rkind) :: gain, buoyancyFact
      type(spectral), pointer :: spectC, spectE
      type(decomp_info), pointer :: gpC => null(), gpE => null(), &
        sp_gpC => null(), sp_gpE => null()
      integer :: nz
      complex(rkind), dimension(:,:,:), allocatable :: fxhat, fyhat, fzhat, fThat
      real(rkind), dimension(:,:,:), allocatable :: fx, fy, fz, fT
      logical :: dumpForce, projectDivergenceFree, isStratified
      type(Pade6Stagg), pointer :: Pade6opZ
      real(rkind) :: maxDiv, maxDivAllTime, avgFact, Re, Pr, oneOnRe, oneOnRePr
      real(rkind) :: dz, Lx, Ly
      character(len=clen) :: outputdir
      real(rkind) :: lambda ! Sets the temperature forcing time scale
      real(rkind) :: Ttgt ! Target temp for forcing layer
      real(rkind), dimension(:,:), allocatable :: zbuff, covBuff
      real(rkind), dimension(:),   allocatable :: meanT, dTdz
      real(rkind), dimension(:,:,:), pointer :: rbuffxC, rbuffyC, rbuffzC, &
        rbuffyE, rbuffzE
      complex(rkind), dimension(:,:,:), pointer :: cbuffyC, cbuffyE, cbuffzC, cbuffzE

      ! Intrinsic scales of the forcing
      real(rkind) :: len_f, tau_f, u_f
      
      ! Variables needed for "duty-cycle" forcing
      real(rkind) :: torigin, ton, toff, rampfact
      logical :: fon, dutyCycleForcing

      contains
        procedure          :: init
        procedure          :: destroy
        procedure          :: updateRHS
        procedure          :: updateScalarRHS
        procedure, private :: interpE2C_real
        procedure, private :: interpE2C_cmplx
        procedure          :: dumpField
        procedure          :: getFringeMask
        generic            :: interpE2C => interpE2C_real, interpE2C_cmplx
        procedure, private :: mean
        procedure, private :: ddz
        procedure, private :: covariance
        procedure, private :: get_time_ampfact
    end type

    contains
      
      subroutine init(this,inputfile,tsim,spectC,spectE,mesh,zEin,gpC,gpE,Pade6opZ,&
          outputdir,sgsModel,buoyancyFact,buoyancyDirection,isStratified,Re,Pr,&
          rbuffxC, rbuffyC, rbuffzC, rbuffyE, rbuffzE, cbuffyC, cbuffyE, &
          cbuffzC, cbuffzE)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), intent(in) :: tsim
          character(len=*), intent(in) :: inputfile
          type(spectral), intent(inout), target :: spectC, spectE
          real(rkind), dimension(:,:,:,:), intent(in), target :: mesh
          real(rkind), dimension(:,:,:), intent(in), target :: zEin
          type(decomp_info), intent(inout), target :: gpC, gpE
          type(Pade6Stagg), intent(inout), target :: Pade6opZ
          characteR(len=*), intent(in) :: outputdir
          type(sgs_igrid), intent(in) :: sgsModel
          real(rkind), intent(in) :: buoyancyFact, Re, Pr
          integer, intent(in) :: buoyancyDirection
          logical, intent(in) :: isStratified
          real(rkind), dimension(:,:,:), intent(in), target :: rbuffxC, rbuffyC, &
            rbuffzC, rbuffyE, rbuffzE
          complex(rkind), dimension(:,:,:), intent(in), target :: cbuffyC, cbuffyE, &
            cbuffzC, cbuffzE
          integer :: ioUnit, ierr, nx, ny
          real(rkind) :: dx, dy
          real(rkind), dimension(:,:,:), pointer :: zC => null(), zE => null()
          real(rkind), dimension(:,:,:), allocatable, target :: zinY
          real(rkind) :: tgtKE = 0.6d0, tgtDissipation = 0.1d0, tgtMPE = 1.d0
          real(rkind) :: zmid = 0.d0, lf = 0.1d0, kmin = 3.d0, kmax = 10.d0
          real(rkind) :: gain = 60.d0
          logical :: dumpForce = .false., projectDivergenceFree = .true.
          integer :: maskType = gaussian
          real(rkind) :: fringe_delta = 0.5d0, lambdaFact_temp = 0.d0
          real(rkind) :: Temp_in_forcing_layer = 1.d0, onoffratio = 1.d0
          integer :: nForceTimeScalesOn = 0 ! The number of forcing time scales to have forcing turned "on". Set to 0 for forcing always on
          real(rkind) :: rampfact = 0.25d0 ! Sets the rate at which the forcing turns on/off

          namelist /spectForceLayer/ zmid, lf, kmin, kmax, tgtKE, tgtMPE,tgtDissipation, &
            gain, dumpForce, maskType, fringe_delta, projectDivergenceFree, lambdaFact_temp, &
            Temp_in_forcing_layer, nForceTimeScalesOn, onoffratio, rampfact
          
          ioUnit = 123
          open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
          read(unit=ioUnit, NML=spectForceLayer)
          close(ioUnit)

          call message(0,'Initializing the spectral forcing layer')
  
          ! Nullify all pointers to be safe
          this%spectC   => null()
          this%spectE   => null()
          this%gpC      => null()
          this%gpE      => null()
          this%Pade6opZ => null()

          this%tgtKE                 = tgtKE
          this%tgtMPE                = tgtMPE
          this%tgtKEon3              = tgtKE/3.d0
          this%tgtDissipation        = tgtDissipation
          this%integralTime          = tgtKE/tgtDissipation
          this%gain                  = gain
          this%dumpForce             = dumpForce
          this%projectDivergenceFree = projectDivergenceFree
          this%spectC => spectC
          this%spectE => spectE
          this%sp_gpC => spectC%spectdecomp
          this%sp_gpE => spectE%spectdecomp
          this%gpC => gpC
          this%gpE => gpE
          this%Pade6opZ => Pade6opZ
          this%outputdir = trim(outputdir)
          this%buoyancyFact = buoyancyFact
          this%isStratified = isStratified
          this%lambda = lambdaFact_temp
          this%Ttgt = Temp_in_forcing_layer

          ! Forcing layer intrinsic scales
          this%len_f = 2.d0*pi/kmin
          this%u_f   = sqrt(tgtKE)
          this%tau_f = this%len_f/this%u_f

          ! If using on/off "duty-cycle" forcing these need to be defined
          if (nForceTimeScalesOn > 0) then
              this%ton = real(nForceTimeScalesOn)*this%tau_f
              this%toff = this%ton/onoffratio
              this%dutyCycleForcing = .TRUE.
          else
              this%ton  = huge(1.d0)
              this%toff = 0.d0
              this%dutyCycleForcing = .FALSE.
              this%ampFact_time = 1.d0
          end if
          this%torigin  = tsim
          this%fon      = .TRUE.
          this%rampfact = rampfact

          nx      = gpC%xsz(1)
          ny      = gpC%ysz(2)
          this%nz = gpC%zsz(3)
          this%avgFact = 1.d0/(real(nx,rkind)*real(ny,rkind))

          dx      = mesh(2,1,1,1) - mesh(1,1,1,1)
          dy      = mesh(1,2,1,2) - mesh(1,1,1,2)
          this%dz = mesh(1,1,2,3) - mesh(1,1,1,3)

          this%Lx = dx*real(nx,rkind)
          this%Ly = dy*real(ny,rkind)

          this%Re        = Re
          this%Pr        = Pr
          this%oneOnRe   = 1.d0/Re
          this%oneOnRePr = 1.d0/(Re*Pr)

          this%rbuffxC => rbuffxC
          this%rbuffyC => rbuffyC
          this%rbuffzC => rbuffzC
          this%rbuffyE => rbuffyE
          this%rbuffzE => rbuffzE
          this%cbuffyC => cbuffyC
          this%cbuffzC => cbuffzC
          this%cbuffyE => cbuffyE
          this%cbuffzE => cbuffzE

          ! User safe-guards
          ! First, make sure the decompositions are the same size in z (assumed in this module)
          call assert(this%sp_gpC%ysz(3) == gpC%ysz(3),&
            'this%sp_gpC%ysz(3) == gpC%ysz(3) -- forcingLayer_spectral.F90')
          call assert(this%sp_gpE%ysz(3) == gpE%ysz(3),&
            'this%sp_gpE%ysz(3) == gpE%ysz(3) -- forcingLayer_spectral.F90')

          ! The amplitude conroller for the forcing assumes the SGS model is an
          ! eddy viscosity model
          call assert(sgsModel%IamEddyViscosityModel(),'The force controller'//&
            ' assumes the SGS term is closed via an EVM')

          ! This class assumes the buoyancy direction is in z
          if (isStratified) call assert(BuoyancyDirection == 3,'BuoyancyDirection == 3 -- forcingLayer_spectral.F90')

          ! Allocate memory
          allocate(this%cmaskC(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3) ))
          allocate(this%jC(   this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3) ))
          allocate(this%cmaskE(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3) ))
          allocate(this%jE(   this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3) ))
          allocate(this%rmaskC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
         
          allocate(this%integralMask(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
          allocate(this%integralMask_1D(this%nz))
          
          allocate(this%fxhat(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3) ))
          allocate(this%fyhat(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3) ))
          allocate(this%fzhat(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3) ))
          allocate(this%fThat(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3) ))

          allocate(this%fx(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
          allocate(this%fy(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
          allocate(this%fz(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
          allocate(this%fT(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))

          allocate(this%zbuff(this%nz,2))
          allocate(this%covBuff(this%nz,2))
          allocate(this%meanT(this%nz),this%dTdz(this%nz))

          if (maskType == gaussian) then
              ! Define Gaussian envelopes
              zC => mesh(:,:,:,3)
              this%rmaskC = gaussianMask(zC,zmid,lf)

              allocate(zinY(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
              call transpose_x_to_y(zC,zinY,gpC)
              nullify(zC)
              zC => zinY(1:this%sp_gpC%ysz(1),1:this%sp_gpC%ysz(2),:)
              !this%cmaskC = exp(-((zC - zmid)/(0.5d0*lf))**2.d0)
              this%cmaskC = gaussianMask(zC,zmid,lf)
              deallocate(zinY)

              zE => zEin
              allocate(zinY(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))
              call transpose_x_to_y(zE,zinY,gpE)
              nullify(zE)
              zE => zinY(1:this%sp_gpE%ysz(1),1:this%sp_gpE%ysz(2),:)
              !this%cmaskE = exp(-((zE - zmid)/(0.5d0*lf))**2.d0)
              this%cmaskE = gaussianMask(zE,zmid,lf)
              
              nullify(zC,zE)
          elseif (maskType == fringeFunction) then
              call this%getFringeMask(mesh(:,:,:,3),gpC,zmid,lf,Fringe_delta,this%rmaskC,'x')
              call this%getFringeMask(mesh(:,:,:,3),gpC,zmid,lf,Fringe_delta,this%cmaskC,'y') 
              call this%getFringeMask(zEin         ,gpE,zmid,lf,Fringe_delta,this%cmaskE,'y') 
          end if
          
          ! For verification purposes
          call this%dumpField(this%cmaskC,this%sp_gpC,'cmaskC',2)
          call this%dumpField(this%cmaskE,this%sp_gpE,'cmaskE',2)

          ! Define spectral envelopes
          this%jC = 1.d0
          this%jE = 1.d0
          where(sqrt(spectC%kabs_sq) < kmin) this%jC = 0.d0
          where(sqrt(spectE%kabs_sq) < kmin) this%jE = 0.d0
          where(sqrt(spectC%kabs_sq) > kmax) this%jC = 0.d0
          where(sqrt(spectE%kabs_sq) > kmax) this%jE = 0.d0
          
          ! Define cutoff mask for integral calculations used for forcing control
          this%integralMask = 1.d0
          where(mesh(:,:,:,3) < -0.5d0*Lf) this%integralMask = 0.d0
          where(mesh(:,:,:,3) >  0.5d0*Lf) this%integralMask = 0.d0

          call transpose_x_to_y(this%integralMask,this%rbuffyC,gpC)
          call transpose_y_to_z(this%rbuffyC,this%rbuffzC,gpC)
          this%integralMask_1D = rbuffzC(1,1,:)

          ! Compute differential volume element
          this%dV = (mesh(2,1,1,1) - mesh(1,1,1,1)) * &
                    (mesh(1,2,1,2) - mesh(1,1,1,2)) * &
                    (mesh(1,1,2,3) - mesh(1,1,1,3))

          if (associated(zC))     nullify(zC)
          if (associated(zE))     nullify(zE)

          call message(1,'Forcing layer initialized successfully!')
          call message(2,'Target KE',this%tgtKE)
          call message(2,'Gain divided by time scale',this%gain/this%integralTime)
      end subroutine

      subroutine updateRHS(this,uhat,vhat,what,u,v,wC,ThatE,T,duidxjC,nSGS,tsim,dt,&
          poiss,urhs,vrhs,wrhs,Trhs)!,scalars)
          class(spectForcingLayer), intent(inout) :: this
          complex(rkind), dimension(:,:,:), intent(in) :: uhat, vhat, what
          real(rkind), dimension(:,:,:), intent(in) :: u, v, wC
          complex(rkind), dimension(:,:,:), intent(in) :: ThatE
          type(PadePoisson), intent(inout) :: poiss
          complex(rkind), dimension(:,:,:), intent(inout) :: urhs, vrhs, wrhs, Trhs
          real(rkind), dimension(:,:,:,:), intent(in), target :: duidxjC
          real(rkind), dimension(:,:,:), intent(in) :: nSGS, T
          real(rkind), intent(in) :: dt, tsim
          !type(scalar_igrid), dimension(:), allocatable :: scalars
          real(rkind) :: maxDiv
          real(rkind) :: KE_x, forceWork_x, eps_x
          real(rkind) :: KE_y, forceWork_y, eps_y
          real(rkind) :: KE_z, forceWork_z, eps_z, BV
          real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, &
            dvdz, dwdx, dwdy, dwdz
          integer :: sca_id
          
          ! Step 1: Construct the force
          this%fxhat = uhat*this%cmaskC*this%jC
          this%fyhat = vhat*this%cmaskC*this%jC
          this%fzhat = what*this%cmaskE*this%jE
          
          ! Step 2: Project out divergence
          call poiss%divergenceCheck(this%fxhat, this%fyhat, this%fzhat, &
            this%rbuffxC, fixDiv = this%projectDivergenceFree, printMessage = .false.)
          this%maxDiv = p_maxval(maxval(abs(this%rbuffxC)))
          this%maxDivAllTime = max(this%maxDivAllTime,this%maxDiv)

          ! Step 3: Remove mean component of the force
          if (this%spectC%carryingZeroK) then
              this%fxhat(this%spectC%zeroK_i,this%spectC%zeroK_j,:) = im0
              this%fyhat(this%spectC%zeroK_i,this%spectC%zeroK_j,:) = im0
              this%fzhat(this%spectC%zeroK_i,this%spectC%zeroK_j,:) = im0
          end if
          
          ! Step 4: Get physical space forcing
          call this%spectC%ifft(this%fxhat,this%fx)
          call this%spectC%ifft(this%fyhat,this%fy)
          call this%spectE%ifft(this%fzhat,this%fz)
         
          ! Step 5: Compute the amplitude
          associate( &
              dudx => duidxjC(:,:,:,1), dudy => duidxjC(:,:,:,2), dudz => duidxjC(:,:,:,3), &
              dvdx => duidxjC(:,:,:,4), dvdy => duidxjC(:,:,:,5), dvdz => duidxjC(:,:,:,6), &
              dwdx => duidxjC(:,:,:,7), dwdy => duidxjC(:,:,:,8), dwdz => duidxjC(:,:,:,9))
         
              ! Integrated kinetic energy 
              KE_x = this%dV*0.5d0*(p_sum(sum(this%integralMask*u *u )))
              KE_y = this%dV*0.5d0*(p_sum(sum(this%integralMask*v *v )))
              KE_z = this%dV*0.5d0*(p_sum(sum(this%integralMask*wC*wC)))

              ! Integrated force work 
              forceWork_x = this%dV*p_sum(sum(this%integralMask*(u *this%fx)))
              forceWork_y = this%dV*p_sum(sum(this%integralMask*(v *this%fy)))
              forceWork_z = this%dV*p_sum(sum(this%integralMask*(wC*this%fz)))

              ! NOTE: This assumes an eddy-viscosity SGS closure
              ! Integrated disspation
              eps_x = this%dV*p_sum(sum(this%integralMask*(this%oneOnRe + nSGS)*(dudx*dudx + dudy*dudy + dudz*dudz)))
              eps_y = this%dV*p_sum(sum(this%integralMask*(this%oneOnRe + nSGS)*(dvdx*dvdx + dvdy*dvdy + dvdz*dvdz)))
              eps_z = this%dV*p_sum(sum(this%integralMask*(this%oneOnRe + nSGS)*(dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)))

              ! Buoyancy transfer
              if (this%isStratified) then
                  this%cbuffyE = this%buoyancyFact*ThatE

                  ! Consistent with the implemented governing equations, subtract the mean
                  ! temperature
                  if (this%spectE%carryingZeroK) then
                      this%cbuffyE = im0
                  end if 
                  call this%interpE2C(this%cbuffyE,this%cbuffyC,this%cbuffzE,this%cbuffzC)
                  call this%spectC%ifft(this%cbuffyC,this%rbuffxC)

                  BV = this%dV*p_sum(sum(this%integralMask*wC*this%rbuffxC))
              else
                  BV = 0.d0
              end if

              ! Compute the amplifiation factor
              this%ampFact_x = (eps_x - this%gain*(KE_x - this%tgtKEon3)/this%integralTime     )/(forceWork_x + 1.d-14)
              this%ampFact_y = (eps_y - this%gain*(KE_y - this%tgtKEon3)/this%integralTime     )/(forceWork_y + 1.d-14)
              this%ampFact_z = (eps_z - this%gain*(KE_z - this%tgtKEon3)/this%integralTime - BV)/(forceWork_z + 1.d-14)
          end associate

          call this%get_time_ampfact(tsim,this%ampFact_time)
          this%ampFact_x = this%ampFact_x*this%ampFact_time
          this%ampFact_y = this%ampFact_y*this%ampFact_time
          this%ampFact_z = this%ampFact_z*this%ampFact_time
          
          ! Step 6: update RHS
          urhs = urhs + this%ampFact_x*this%fxhat
          vrhs = vrhs + this%ampFact_y*this%fyhat
          wrhs = wrhs + this%ampFact_z*this%fzhat

          if (this%isStratified) call this%updateScalarRHS(T,dt,Trhs)
          !if (allocated(scalars)) then
          !    do sca_id = 1,size(scalars)
          !        call this%updateScalarRHS(scalars(sca_id)%F,dt,scalars(sca_id)%rhs)
          !    end do
          !end if
      end subroutine

      subroutine updateScalarRHS(this,T,dt,Trhs)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: T
          real(rkind), intent(in) :: dt
          complex(rkind), dimension(:,:,:), intent(inout) :: Trhs

          if (this%isStratified .and. dt > 0.d0) then
              ! Penalty term
              this%ampFact_T = 1.d0
              this%fT = this%rmaskC*this%lambda/dt*(this%Ttgt - T)
              call this%spectC%fft(this%fT,this%fThat)
              Trhs = Trhs + this%ampFact_T*this%fThat
          end if
      end subroutine

      subroutine get_time_ampfact(this,tsim,ampFact)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), intent(in) :: tsim
          real(rkind), intent(out) :: ampFact
          real(rkind), parameter :: tol = 1.d-10

          if (this%dutyCycleForcing) then
              if (this%fon) then
                  ampFact = 1.d0 - 0.5d0*(1.d0 + tanh((tsim - (this%torigin + this%ton))/this%rampFact))
                  if ((tsim > this%torigin + this%ton) .and. (ampFact < tol)) then
                      this%fon = .FALSE.
                      this%torigin = this%torigin + this%ton
                  end if
              else
                  ampFact = 0.5d0*(1.d0 + tanh((tsim - (this%torigin + this%toff))/this%rampFact))
                  if ((tsim > this%torigin + this%toff) .and. ((1.d0 - ampFact) < tol)) then
                      this%fon = .TRUE.
                      this%torigin = this%torigin + this%toff
                  end if
              end if
          else
              ampFact = 1.d0
          end if
      end subroutine
      
      subroutine covariance(this,f,g,fg,rbuffxC,cbuffyC,cbuffzC)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: f, g
          real(rkind), dimension(:), intent(out) :: fg
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffxC
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffyC, cbuffzC

          rbuffxC = f*g
          call this%mean(rbuffxC,fg,cbuffyC,cbuffzC)
          call this%mean(f,this%covBuff(:,1),cbuffyC,cbuffzC)
          call this%mean(g,this%covBuff(:,2),cbuffyC,cbuffzC)

          fg = fg - this%covBuff(:,1)*this%covBuff(:,2)
          
      end subroutine

      subroutine interpE2C_cmplx(this,fE,fC,cbuffzE,cbuffzC)
          class(spectForcingLayer), intent(inout) :: this
          complex(rkind), dimension(:,:,:), intent(in) :: fE
          complex(rkind), dimension(:,:,:), intent(out) :: fC
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffzE, cbuffzC
          call transpose_y_to_z(fE,cbuffzE,this%sp_gpE)
          call this%pade6opz%interpz_E2C(cbuffzE,cbuffzC,zbc_bot,zbc_top)
          call transpose_z_to_y(cbuffzC,fC,this%sp_gpC)
      end subroutine
      
      subroutine interpE2C_real(this,fE,fC,rbuffyC,rbuffzC,rbuffyE,rbuffzE)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: fE
          real(rkind), dimension(:,:,:), intent(out) :: fC
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffyC, rbuffyE
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffzC, rbuffzE
          
          call transpose_x_to_y(fE,rbuffyE,this%gpE)
          call transpose_y_to_z(rbuffyE,rbuffzE,this%gpE)
          call this%pade6opZ%interpz_E2C(rbuffzE, rbuffzC, zbc_bot, zbc_top)
          call transpose_z_to_y(rbuffzC, rbuffyC, this%gpC)
          call transpose_y_to_x(rbuffyC, fC, this%gpC)
      end subroutine
      
      subroutine destroy(this)
          class(spectForcingLayer), intent(inout) :: this
          if (allocated(this%cmaskC)) deallocate(this%cmaskC)
          if (allocated(this%cmaskE)) deallocate(this%cmaskE)
          if (allocated(this%jC)) deallocate(this%jC)
          if (allocated(this%jE)) deallocate(this%jE)
          if (allocated(this%fxhat)) deallocate(this%fxhat)
          if (allocated(this%fyhat)) deallocate(this%fyhat)
          if (allocated(this%fzhat)) deallocate(this%fzhat)
          if (allocated(this%fx)) deallocate(this%fx)
          if (allocated(this%fy)) deallocate(this%fy)
          if (allocated(this%fz)) deallocate(this%fz)
          if (allocated(this%integralMask)) deallocate(this%integralMask)

          if (associated(this%spectC)) nullify(this%spectC)
          if (associated(this%spectE)) nullify(this%spectE)
          if (associated(this%gpC)) nullify(this%gpC)
          if (associated(this%gpE)) nullify(this%gpE)
          if (associated(this%Pade6opZ)) nullify(this%Pade6opZ)

      end subroutine

      subroutine dumpField(this,field,gp,fieldName,coord)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: field
          type(decomp_info), intent(inout) :: gp
          character(len=*), intent(in) :: fieldName
          integer, intent(in) :: coord
          character(len=clen) :: fname

          write(fname,'(A,I2.2,A)') trim(this%outputdir)//'/forcingLayer_'//trim(fieldName)//'.out'
          call decomp_2d_write_one(coord,field,trim(fname),gp)
      end subroutine

      function gaussianMask(z,zmid,lf) result (mask)
          real(rkind), dimension(:,:,:), intent(in) :: z
          real(rkind), intent(in) :: zmid, lf
          real(rkind), dimension(size(z,1),size(z,2),size(z,3)) :: mask

          mask = exp(-((z - zmid)/(0.5d0*lf))**2.d0)

      end function
      
      subroutine getFringeMask(this,zinX,gp,zmid,lf,Fringe_delta,mask,pencil)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: zinX
          type(decomp_info), intent(inout) :: gp
          real(rkind), intent(in) :: zmid, lf, Fringe_delta
          real(rkind), dimension(:,:,:), intent(out) :: mask
          character(len=1), intent(in) :: pencil
          real(rkind), dimension(:), allocatable :: z1D, z1, z2, &
            S1, S2, fringeFunc
          real(rkind), dimension(:,:,:), allocatable :: maskInX, maskInY
          real(rkind) :: Fringe_zst, Fringe_zen
          integer :: i, j, k

          allocate(z1D(gp%xsz(3)))
          allocate(z1(gp%xsz(3)), S1(gp%xsz(3)))
          allocate(z2(gp%xsz(3)), S2(gp%xsz(3)))
          allocate(fringeFunc(gp%xsz(3)))
          allocate(maskInX(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
          allocate(maskInY(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
         
          z1D = zinX(1,1,:)
          
          Fringe_zst = zmid - lf/2.d0
          Fringe_zen = zmid + lf/2.d0
          z1 = ((z1D -  Fringe_zst)/Fringe_delta)
          z2 = ((z1D -  Fringe_zen)/Fringe_delta) + 1.d0
          call S_fringe(z1, S1)
          call S_fringe(z2, S2)
          fringeFunc = S1 - S2
          
          maskInX = 0.d0

          !do k = 1,gp%xsz(3)
              do j = 1,gp%xsz(2)
                  do i = 1,gp%xsz(1)
                      maskInX(i,j,:) = fringeFunc(:)
                      !maskInX(i,j,k) = maskInX(i,j,k) + fringeFunc(k)
                  end do 
              end do
          !end do

          where(maskInX > 1.d0) maskInX = 1.d0

          if (pencil == 'x') then
              mask = maskInX
          elseif (pencil == 'y') then
              call transpose_x_to_y(maskInX,maskInY,gp)
              call assert(size(mask,3) == size(maskInY,3),'size(mask,3) == size(maskInY,3) -- forcingLayer_spectral.F90')
              do j = 1,size(mask,2)
                  do i = 1,size(mask,1)
                      mask(i,j,:) = maskInY(1,1,:)
                  end do
              end do
          else
              call assert(.false.,'Only x or y pencils supported for this routine -- forcingLayer_spectral.F90')
          end if

          deallocate(z1D, z1, z2, S1, S2, fringeFunc, maskInX, maskInY)
      end subroutine

      subroutine mean(this, f, fmean, cbuffyC, cbuffzC)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:), intent(out) :: fmean
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffyC, cbuffzC
          integer :: ierr

          call this%spectC%fft(f, cbuffyC)
          call transpose_y_to_z(cbuffyC, cbuffzC, this%sp_gpC)
          if (nrank == 0) then
              fmean = real(cbuffzC(1,1,:),rkind)*this%avgFact
          else
              fmean = 0.d0 ! Only 0 processor has the actual mean  
          end if 
          call MPI_BCAST(fmean,this%nz,MPIRKIND,0,MPI_COMM_WORLD, ierr)

      end subroutine 

      subroutine ddz(this, f, ddz_f, bc_bot, bc_top)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:), intent(in)  :: f
          real(rkind), dimension(:), intent(out) :: ddz_f
          integer, intent(in) :: bc_bot, bc_top 

          this%zbuff(:,1) = f
          call this%Pade6opZ%ddz_1d_C2C(this%zbuff(:,1),this%zbuff(:,2), bc_bot, bc_top)
          ddz_f = this%zbuff(:,2)

      end subroutine 
         
end module spectralForcingLayerMod
