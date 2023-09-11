module jet_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    use FiltersMod,       only: filters
    implicit none
    
    real(rkind) :: gam = 1.4_rkind
    real(rkind) :: Rgas = one
    real(rkind), allocatable, dimension(:,:) :: u_left_bc

    !real(rkind) :: pL = one
    !real(rkind) :: pR = 0.1_rkind

    ! Gaussian filter for sponge
    type(filters) :: mygfil
    logical :: periodicx = .false., periodicy = .true., periodicz = .true.
end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind, clen
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info
    use jet_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k,ioUnit
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    character(clen) :: inputfile='input.dat'

    real(rkind) :: xcore, AA, D, U0, Uco, ycen, zcen
    real(rkind) :: narg, uc, raddist, radialfac, sigbyDsq
    real(rkind) :: Lx, Ly, Lz, Re, Pr, Sc, T_ref

    namelist /PROBINPUT/ Lx, Ly, Lz, U0, Uco, Re, Pr, Sc, &
                        T_ref, D, ycen, zcen

    ioUnit = 11
    open(unit=ioUnit, file=inputfile, form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)
    
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! Base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nx-1,rkind)
        dy = Ly/real(ny-0,rkind)   ! periodic in y
        dz = Lz/real(nz-0,rkind)   ! periodic in z
        
        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use IdealGasEOS,                 only: idealgas
    use SutherLandViscosityMod,      only: sutherlandViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use SutherLandConductivityMod,   only:SutherLandConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
 
    use jet_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim,tstop,dt,tviz

    type(sutherlandViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(sutherlandConductivity) :: thermcond
    !real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: Lx, Ly, Lz, Re, Pr, Sc, T_ref, S, Sk, T0
    real(rkind) :: xcore, AA, D, U0, Uco, ycen, zcen, mu_ref
    real(rkind) :: xl, yl, zl, sig
    real(rkind) :: narg, uc, raddist, radialfac, sigbyDsq
    integer :: icore, i, j, k, ioUnit

        namelist /PROBINPUT/ Lx, Ly, Lz, U0, Uco, Re, Pr, Sc, &
                            T_ref, D, ycen, zcen

    associate( rho => fields(:,:,:,rho_index), u => fields(:,:,:,u_index), &
                 v => fields(:,:,:,  v_index), w => fields(:,:,:,w_index), &
                 p => fields(:,:,:,  p_index), T => fields(:,:,:,T_index), &
                 e => fields(:,:,:,  e_index),                             &
                 Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),           &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        

        !!!!--------------------read input values--------------------------!!!!
        ! default values!, zcen
        xcore = 0.5d0; AA = 4.0d0; D = 1.0d0; U0 = 1.0d0;
        ycen = Ly/2.0d0; zcen = Lz/2.0d0; Uco = 0.05d0
        if(decomp%ysz(3)==1) then
            !! there is only 1 point along z
            zcen = zero
        endif

        ioUnit = 11
        open(unit=ioUnit, file=inputfile, form='FORMATTED')
        read(unit=ioUnit, NML=PROBINPUT)
        close(ioUnit)
    
        !tmp = half * ( one+tanh(x/(dx)) )

        !!!!--------------------set velocity values------------------------!!!!
        icore = minloc(abs(x(:,1,1)-xcore),1)
        do k=1,size(u,3)
         do j=1,size(u,2)
          do i=1,size(u,1)
              xl = x(i,j,k); yl = y(i,j,k); zl = z(i,j,k)
              narg = 2.0d0 + AA*erfc(xl-xcore)
              if(xl < xcore) then
                  uc = U0
                  sig = 0.32d0*D
              else
                  uc = U0 * xcore / xl
                  sig = 0.32d0*D + 0.15d0* (xl - xcore)
              endif
              raddist   = sqrt((yl-ycen)**2 + (zl-zcen)**2)/D
              sigbyDsq  = (sig/D)**2
              radialfac = exp(-0.5d0*raddist**narg/sigbyDsq)
              u(i,j,k) = uc * radialfac
              if(i==1) then
                write(*,'(7(e19.12,1x))') y(1,j,1), narg, raddist, sigbyDsq, radialfac, xl, uc
              endif
              if (xl > xcore) then
                 u(i,j,k) = 0.0
              endif
          enddo
         enddo
        enddo
        u = u + Uco
        v   = zero
        w   = zero
        print *, '--Done initialization-'
        !!!!--------------------done setting velocity values---------------!!!!

        !!!!--------------------set thermodynamic properties---------------!!!!
        T0 =  101325.0_rkind/(287.0_rkind*1.2_rkind)
        S = 110.4_rkind/T0
        Sk = 194.0_rkind/T0
        mu_ref = one/Re

        ! Set material's transport coefficient object
        shearvisc = sutherlandViscosity(mu_ref, T_ref, 1.5_rkind, S)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = sutherlandConductivity(Pr, T_ref, 1.5_rkind, Sk)
        call mix%set_material( 1, idealgas(gam, Rgas), &
             shearvisc = shearvisc, & 
             bulkvisc  = bulkvisc, &
             thermcond = thermcond  )

        rho = one;    p = one;    T = p/(rho*Rgas) 
        !!!!-----------done setting thermodynamic properties----------------!!!!

        !!!!------------------set mixture properties------------------------!!!!
        Ys(:,:,:,1)  = one;  call mix%update(Ys)
        !!!!------------done setting mixture properties---------------------!!!!

        !!!!------------store left bc for using later-----------------------!!!!
        if(decomp%yst(1) == 1) then 
            allocate(u_left_bc(size(u,2), size(u,3)))
            !! left-boundary processor
            u_left_bc = u(1,:,:)
            write(*,*) '----In initfields----'
            do j=1,size(u,2)
                write(*,'(3(e19.12,1x))') y(1,j,1), u_left_bc(j,1), u(1,j,1)
            enddo
        endif

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, &
                        "gaussian", "gaussian", "gaussian" )
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use jet_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile
    integer :: i

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        !! perhaps write jet-centerline properties later
        !!write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/jet_", vizcount, ".dat"

        !!open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        !!do i=1,decomp%ysz(1)
        !!    write(outputunit,'(8ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
        !!                                   mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        !!
        !!end do
        !!close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use jet_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    real(rkind) :: dx, filpt, thickT
    integer :: i, j

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! set Dirichlet BC for velocity, density, pressure
        if(decomp%yst(1) == 1) then 
            !! left-boundary processor
            u(1,:,:) = u_left_bc
            v(1,:,:) = zero
            w(1,:,:) = zero
            rho(1,:,:) = one
            p(1,:,:) = one
         endif

        ! Sponge+bulk for exit bc
        ! Gradually apply the exit boundary conditions
        dx = x(2,1,1) - x(1,1,1)
        filpt = 3.0_rkind/dx 
        thickT = real(5.0D0, rkind)
        
        ! Gussian Filter for 
        do i=1,decomp%ysz(1)
            dumT(i,:,:)=half*(one-tanh( (real(decomp%xsz(1) - (decomp%yst(1) - 1 + i - 1), rkind)-filpt) / thickT ))
        end do
            
        dumF = u
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        u = u + dumT*(dumF-u) 

        dumF = v
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        v = v + dumT*(dumF-v)

        dumF = w
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        w = w + dumT*(dumF-w)

        dumF = p
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        p = p + dumT*(dumF-p)

        dumF = rho
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        rho = rho + dumT*(dumF-rho)

        !write(*,*) '----In hook_bc----'
        !do j=1,size(u,2)
        !    write(*,'(3(e19.12,1x))') y(1,j,1), u_left_bc(j,1), u(1,j,1)
        !enddo

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim,sgsmodel)
!subroutine hook_timestep(decomp,mesh,fields,mix,tsim)
    use kind_parameters,  only: rkind
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use exits,            only: message
    use sgsmod_cgrid,     only: sgs_cgrid
    use reductions,       only: P_MAXVAL,P_MINVAL

    use jet_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: tsim
    integer,                         intent(in) :: step
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid), optional,       intent(in) :: sgsmodel

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Minimum bulk viscosity",P_MINVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,   only: mixture

    use jet_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
