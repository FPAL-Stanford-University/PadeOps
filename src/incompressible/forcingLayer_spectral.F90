module spectralForcingLayerMod
  ! This module implements a forcing layer similar to Briggs et al. (1996), but with some modifications
    use kind_parameters, only: rkind, clen, mpirkind
    use reductions,      only: p_maxval, p_sum
    use gridtools,       only: onThisRank, getStEndIndices
    use constants,       only: rmaxInt, im0
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
    implicit none

    integer, parameter :: zbc_bot = 0, zbc_top = 0
    integer, parameter :: gaussian = 1, fringeFunction = 2
    type :: spectForcingLayer

      real(rkind), dimension(:,:,:), allocatable :: maskC, maskE, jC, jE
      real(rkind), dimension(:,:,:), allocatable :: integralMask
      real(rkind) :: seedFact, dV, tgtKE, ampFact, tgtDissipation, integralTime 
      integer :: seed
      real(rkind) :: gain
      type(spectral), pointer :: spectC, spectE
      type(decomp_info), pointer :: gpC, gpE
      complex(rkind), dimension(:,:,:), allocatable :: fxhat, fyhat, fzhat
      real(rkind), dimension(:,:,:), allocatable :: fx, fy, fz
      logical :: dumpForce, projectDivergenceFree
      type(Pade6Stagg), pointer :: Pade6opZ
      real(rkind) :: maxDiv, maxDivAllTime
      character(len=clen) :: outputdir

      contains
        procedure :: init
        procedure :: destroy
        procedure :: updateRHS
        procedure :: interpE2C
        procedure :: dumpField
        procedure :: getFringeMask
    end type

    contains
      
      subroutine init(this,inputfile,spectC,spectE,mesh,zEin,gpC,gpE,Pade6opZ,outputdir)
          class(spectForcingLayer), intent(inout) :: this
          character(len=*), intent(in) :: inputfile
          type(spectral), intent(inout), target :: spectC, spectE
          real(rkind), dimension(:,:,:,:), intent(in), target :: mesh
          real(rkind), dimension(:,:,:), intent(in), target :: zEin
          type(decomp_info), intent(inout), target :: gpC, gpE
          type(Pade6Stagg), intent(inout), target :: Pade6opZ
          characteR(len=*), intent(in) :: outputdir
          integer :: ioUnit, ierr
          real(rkind), dimension(:,:,:), pointer :: zC => null(), zE => null()
          real(rkind), dimension(:,:,:), allocatable, target :: zinY
          type(decomp_info), pointer :: sp_gpC => null(), sp_gpE => null()
          real(rkind) :: tgtKE = 0.6d0, tgtDissipation = 0.1d0
          real(rkind) :: zmid = 0.d0, lf = 0.1d0, kmin = 3.d0, kmax = 10.d0
          real(rkind) :: gain = 60.d0
          logical :: dumpForce = .false., projectDivergenceFree = .true.
          integer :: maskType = gaussian
          real(rkind) :: fringe_delta = 0.5d0

          namelist /spectForceLayer/ zmid, lf, kmin, kmax, tgtKE, tgtDissipation, &
            gain, dumpForce, maskType, fringe_delta, projectDivergenceFree
          
          ioUnit = 123
          open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
          read(unit=ioUnit, NML=spectForceLayer)
          close(ioUnit)
  
          ! Nullify all pointers to be safe
          this%spectC   => null()
          this%spectE   => null()
          this%gpC      => null()
          this%gpE      => null()
          this%Pade6opZ => null()

          this%tgtKE                 = tgtKE
          this%tgtDissipation        = tgtDissipation
          this%integralTime          = tgtKE/tgtDissipation
          this%gain                  = gain
          this%dumpForce             = dumpForce
          this%projectDivergenceFree = projectDivergenceFree
          this%spectC => spectC
          this%spectE => spectE
          sp_gpC => spectC%spectdecomp
          sp_gpE => spectE%spectdecomp
          this%gpC => gpC
          this%gpE => gpE
          this%Pade6opZ => Pade6opZ
          this%outputdir = trim(outputdir)

          ! First, make sure the decompositions are the same size in z (assumed in this module)
          call assert(sp_gpC%ysz(3) == gpC%ysz(3),&
            'sp_gpC%ysz(3) == gpC%ysz(3) -- forcingLayer_spectral.F90')
          call assert(sp_gpE%ysz(3) == gpE%ysz(3),&
            'sp_gpE%ysz(3) == gpE%ysz(3) -- forcingLayer_spectral.F90')

          ! Allocate memory
          allocate(this%maskC(sp_gpC%ysz(1), sp_gpC%ysz(2), sp_gpC%ysz(3) ))
          allocate(this%jC(   sp_gpC%ysz(1), sp_gpC%ysz(2), sp_gpC%ysz(3) ))
          allocate(this%maskE(sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3) ))
          allocate(this%jE(   sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3) ))
         
          allocate(this%integralMask(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
          
          allocate(this%fxhat(sp_gpC%ysz(1), sp_gpC%ysz(2), sp_gpC%ysz(3) ))
          allocate(this%fyhat(sp_gpC%ysz(1), sp_gpC%ysz(2), sp_gpC%ysz(3) ))
          allocate(this%fzhat(sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3) ))

          allocate(this%fx(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
          allocate(this%fy(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
          allocate(this%fz(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))

          if (maskType == gaussian) then
              ! Define Gaussian envelopes
              zC => mesh(:,:,:,3)
              allocate(zinY(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
              call transpose_x_to_y(zC,zinY,gpC)
              nullify(zC)
              zC => zinY(1:sp_gpC%ysz(1),1:sp_gpC%ysz(2),:)
              this%maskC = exp(-((zC - zmid)/(0.5d0*lf))**2.d0)
              deallocate(zinY)

              zE => zEin
              allocate(zinY(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))
              call transpose_x_to_y(zE,zinY,gpE)
              nullify(zE)
              zE => zinY(1:sp_gpE%ysz(1),1:sp_gpE%ysz(2),:)
              this%maskE = exp(-((zE - zmid)/(0.5d0*lf))**2.d0)
              
              nullify(zC,zE)
          elseif (maskType == fringeFunction) then
              call this%getFringeMask(mesh(:,:,:,3),gpC,sp_gpC,zmid,lf,Fringe_delta,this%maskC) 
              call this%getFringeMask(zEin         ,gpE,sp_gpE,zmid,lf,Fringe_delta,this%maskE) 
          end if
          
          ! For verification purposes
          call this%dumpField(this%maskC,sp_gpC,'maskC',2)
          call this%dumpField(this%maskE,sp_gpE,'maskE',2)

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

          ! Compute differential volume element
          this%dV = (mesh(2,1,1,1) - mesh(1,1,1,1)) * &
                    (mesh(1,2,1,2) - mesh(1,1,1,2)) * &
                    (mesh(1,1,2,3) - mesh(1,1,1,3))

          if (associated(zC))     nullify(zC)
          if (associated(zE))     nullify(zE)
          if (associated(sp_gpC)) nullify(sp_gpC)
          if (associated(sp_gpE)) nullify(sp_gpE)
      end subroutine

      subroutine updateRHS(this,uhat,vhat,what,u,v,wC,duidxjC,nSGS,rbuffxC,poiss,Re,urhs,vrhs,wrhs)
          class(spectForcingLayer), intent(inout) :: this
          complex(rkind), dimension(:,:,:), intent(in) :: uhat, vhat, what
          real(rkind), dimension(:,:,:), intent(in) :: u, v, wC
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffxC
          type(PadePoisson), intent(inout) :: poiss
          complex(rkind), dimension(:,:,:), intent(inout) :: urhs, vrhs, wrhs
          real(rkind), dimension(:,:,:,:), intent(in), target :: duidxjC
          real(rkind), dimension(:,:,:), intent(in) :: nSGS
          real(rkind), intent(in) :: Re
          real(rkind) :: maxDiv, KE, forceWork, eps
          real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, &
            dvdz, dwdx, dwdy, dwdz
          
          ! Step 1: Construct the force
          this%fxhat = uhat*this%maskC*this%jC
          this%fyhat = vhat*this%maskC*this%jC
          this%fzhat = what*this%maskE*this%jE
          
          ! Step 2: Project out divergence
          call poiss%divergenceCheck(this%fxhat, this%fyhat, this%fzhat, &
            rbuffxC, fixDiv = this%projectDivergenceFree, printMessage = .false.)
          this%maxDiv = p_maxval(maxval(abs(rbuffxC)))
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
              KE = this%dV*0.5d0*(p_sum(sum(this%integralMask*(u*u + v*v + wC*wC))))

              ! Integrated force work
              forceWork = this%dV*p_sum(sum(this%integralMask*(u*this%fx + v*this%fy + wC*this%fz)))

              ! Integrated disspation
              eps = this%dV*p_sum(sum(this%integralMask*(1.d0/Re + nSGS)*(&
                dudx*dudx + dudy*dudy + dudz*dudz + &
                dvdx*dvdx + dvdy*dvdy + dvdz*dvdz + &
                dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)))

              ! Compute the amplifiation factor
              this%ampFact = (eps - this%gain*(KE - this%tgtKE)/this%integralTime)/(forceWork + 1.d-14)
          end associate
          
          ! Step 6: update RHS
          urhs = urhs + this%ampFact*this%fxhat
          vrhs = vrhs + this%ampFact*this%fyhat
          wrhs = wrhs + this%ampFact*this%fzhat

      end subroutine
      
      subroutine interpE2C(this,fE,fC,rbuffyC,rbuffzC,rbuffyE,rbuffzE)
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
          if (allocated(this%maskC)) deallocate(this%maskC)
          if (allocated(this%maskE)) deallocate(this%maskE)
          if (allocated(this%jC)) deallocate(this%jC)
          if (allocated(this%jE)) deallocate(this%jE)
          if (allocated(this%fxhat)) deallocate(this%fxhat)
          if (allocated(this%fyhat)) deallocate(this%fyhat)
          if (allocated(this%fzhat)) deallocate(this%fzhat)
          if (allocated(this%fx)) deallocate(this%fx)
          if (allocated(this%fy)) deallocate(this%fy)
          if (allocated(this%fz)) deallocate(this%fz)

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
      
      subroutine getFringeMask(this,zinX,gp,sp_gp,zmid,lf,Fringe_delta,mask)
          class(spectForcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: zinX
          type(decomp_info), intent(inout) :: gp, sp_gp
          real(rkind), intent(in) :: zmid, lf, Fringe_delta
          real(rkind), dimension(:,:,:), intent(out) :: mask
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

          do k = 1,gp%xsz(3)
              do j = 1,gp%xsz(2)
                  do i = 1,gp%xsz(1)
                      maskInX(i,j,k) = maskInX(i,j,k) + fringeFunc(k)
                  end do 
              end do
          end do

          where(maskInX > 1.d0) maskInX = 1.d0

          call transpose_x_to_y(maskInX,maskInY,gp)
          call assert(size(mask,3) == size(maskInY,3),'size(mask,3) == size(maskInY,3) -- forcingLayer_spectral.F90')
          do j = 1,size(mask,2)
              do i = 1,size(mask,1)
                  mask(i,j,:) = maskInY(1,1,:)
              end do
          end do

          deallocate(z1D, z1, z2, S1, S2, fringeFunc, maskInX, maskInY)
      end subroutine
         
end module spectralForcingLayerMod
