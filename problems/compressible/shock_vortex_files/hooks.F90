module Shock_vortex_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info, nrank
    use basic_io,         only: read_2d_ascii 
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message, message_min_max
    use mpi
    implicit none
    integer, parameter :: ns = 1

    real(rkind), allocatable, dimension(:,:) :: u_left_bc, p_left_bc, rho_left_bc, T_left_bc
    ! Problem parameters
    real(rkind) :: Sc = 1._rkind              ! Schmidt number
    real(rkind) :: u_max   = 0.9301_rkind     ! max vel
    real(rkind) :: p_max   = 0.3041_rkind     ! max press
    real(rkind) :: T_max   = 1.1426_rkind     ! max temp
    real(rkind) :: rho_max = 0.2661_rkind     ! max density
    real(rkind) :: rho_ref = 0.125_rkind      ! ref density
    real(rkind) :: p_ref   = 0.1_rkind        ! ref pressure
    
    ! Parameters for the 2 materials
    real(rkind):: gam=1.4_rkind
    real(rkind):: Pr=0.7_rkind 
    real(rkind):: nu_ref=0.02_rkind 
    real(rkind):: Rgas=one

    ! Domain size data
    real(rkind) :: Ly = 7._rkind, Lx=7._rkind, Lz=7._rkind
    real(rkind) :: x1, y1, z1
    logical :: periodicx = .false., periodicy = .false., periodicz = .false.

    ! Gaussian filter for sponge
    type(filters) :: mygfil

end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use Shock_vortex_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    character(clen) :: inputfile='input.dat'

    namelist /PROBINPUT/ Lx, Ly, Lz, Pr, Sc, nu_ref, p_ref, rho_ref , &
                        T_max, p_max, rho_max, u_max
    ioUnit = 11
    open(unit=ioUnit, file=inputfile, form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)
    
    ! Global domain size 
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)

    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        if (nrank == 0) then
            print *, "Domain size: ",Lx,Ly,Lz,ix1,iy1,iz1,ixn,iyn,izn
        end if
        dx = Lx/real(nx-1,rkind)
        dy = Ly/real(ny-1,rkind)
        dz = Lz/real(nz-1,rkind)   ! Change here for 2d
  
        x1 = 0._rkind! -Lx/2._rkind
        y1 = -Ly/2._rkind
        z1 = -Lz/2._rkind
        
        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 + k - 1, rkind ) * dz
                end do
            end do
        end do
    end associate
end subroutine


subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,four,five,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,&
                                           p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info,nrank
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use SutherLandViscosityMod,      only: sutherlandViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use SutherLandConductivityMod,   only: SutherLandConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  

    use Shock_vortex_data
    use mpi

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    type(sutherlandViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(sutherlandConductivity) :: thermcond
    real(rkind) :: S, Sk, T0, mu_ref , T_ref! U0, rad
    integer :: i,j, k, iounit, nx, ny, nz, nxl, nyl, nzl
    character(len=clen) :: outputfile
    real(rkind), dimension(decomp%ysz(1)) :: x_new
    real(rkind), dimension(decomp%ysz(2)) :: y_new
    real(rkind), dimension(decomp%ysz(3)) :: z_new
    
    namelist /PROBINPUT/ Lx, Ly, Lz, Pr, Sc, nu_ref, p_ref, rho_ref , &
                        T_max, p_max, rho_max, u_max

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Global domain sizes
    nx = decomp%xsz(1);     ny = decomp%ysz(2);    nz = decomp%zsz(3)

    ! Local domain sizes
    nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)
    
    !print *, "checking nx,ny,nz,n ",nx,ny,nz,nxl,nyl,nzl
    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), &
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)
        
        T0     =  101325.0_rkind/(287.0_rkind*1.2_rkind)        
        S      = 110.4_rkind/T0
        Sk     = 194.0_rkind/T0

        mu_ref = nu_ref * rho_ref

        T_ref = p_ref / rho_ref
        ! Set each material's transport coefficient object
        shearvisc = sutherlandViscosity(mu_ref, T_ref, 1.5_rkind, S)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = sutherlandConductivity(Pr, T_ref, 1.5_rkind, Sk)
        call mix%set_material( 1, idealgas( gam, Rgas ),&
             shearvisc = shearvisc, & 
             bulkvisc  = bulkvisc, &
             thermcond = thermcond  )

        Ys(:,:,:,1)  = one 
        call mix%update(Ys)     

        
        ! Add base flow profiles
        u = u + zero
        v = v + zero
        w = w + zero
        p = p + p_ref
        rho = rho + rho_ref
        T = p/(rho*Rgas) 

        !!!!------------store left bc for using later-----------------------!!!!
        !if(decomp%yst(1) == 1) then 
        !    allocate(u_left_bc(size(u,2), size(u,3)))
        !    allocate(p_left_bc(size(p,2), size(p,3)))
        !    allocate(rho_left_bc(size(rho,2), size(rho,3)))
        !    allocate(T_left_bc(size(T,2), size(T,3)))
        !    !! left-boundary processor
        !    u_left_bc   = u(1,:,:)
        !    p_left_bc   = p(1,:,:)
        !    rho_left_bc = rho(1,:,:)
        !    T_left_bc   = T(1,:,:)
        !endif

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, &
                        "gaussian", "gaussian", "gaussian" )
    end associate
end subroutine


subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture
    use reductions,       only: P_MEAN
    use Shock_vortex_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    character(len=clen) :: outputfile,str
    integer :: i,outputunit=229

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Shock_vortex_t.dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(8ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        
        end do
        close(outputunit)
    end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info, nrank
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use Shock_vortex_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i, j, k
    real(rkind) :: dx, rad, filpt, thickT, U0, P0, rho0, T0, ta = 0.006_rkind, tb = 0.38_rkind
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    character(len=clen) :: outputfile

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        !U_vel
        !if (tsim .le. (ta + (tb/2))) then
        !   U0 = u_max * (0.499_rkind * (tanh((34.83_rkind/ta) * (tsim-ta)) + one))
        !   P0 = p_max * (0.334_rkind * (tanh((43.37_rkind/ta) * (tsim-ta)) + 1.985_rkind))
        !   rho0 = rho_max * (0.264_rkind * (tanh((41.346_rkind/ta) * (tsim-ta)) + 2.778_rkind))
        !   T0 = T_max * (0.149_rkind * (tanh((42.194_rkind/ta) * (tsim-ta)) + 5.678_rkind))
        !else
        !   U0 =  u_max * (0.499_rkind * (tanh((34.83_rkind/ta) * (ta+tb-tsim)) + one))
        !   P0 =  p_max * (0.334_rkind * (tanh((43.37_rkind/ta) * (ta+tb-tsim)) + 1.985_rkind))
        !   rho0 =  rho_max * (0.264_rkind * (tanh((41.346_rkind/ta) * (ta+tb-tsim)) + 2.778_rkind))
        !   T0 =  T_max * (0.149_rkind * (tanh((42.194_rkind/ta) * (ta+tb-tsim)) + 5.678_rkind))
        !endif


        U0   = 0.487_rkind + (0.97538_rkind - 0.487_rkind) * tanh(246.78_rkind * (tsim - 0.0095_rkind));
        P0   = 0.209_rkind + (0.31814_rkind - 0.209_rkind) * tanh(277.79_rkind * (tsim - 0.0113_rkind));
        rho0 = 0.201_rkind + (0.27652_rkind - 0.201_rkind) * tanh(263.71_rkind * (tsim - 0.0110_rkind));
        T0   = 0.979_rkind + (1.15650_rkind - 0.979_rkind) * tanh(283.38_rkind * (tsim - 0.0145_rkind));
       
        !if (tsim .le. (ta + (tb/2))) then
        !    U0   = 0.487_rkind + (0.97538_rkind - 0.487_rkind) * tanh(246.78_rkind * (tsim - ta - 0.0095_rkind));
        !    P0   = 0.209_rkind + (0.31814_rkind - 0.209_rkind) * tanh(277.79_rkind * (tsim - ta - 0.0113_rkind));
        !    rho0 = 0.201_rkind + (0.27652_rkind - 0.201_rkind) * tanh(263.71_rkind * (tsim - ta - 0.0110_rkind));
        !    T0   = 0.979_rkind + (1.15650_rkind - 0.979_rkind) * tanh(283.38_rkind * (tsim - ta - 0.0145_rkind));
        !else
        !    U0   = 0.487_rkind + (0.97538_rkind - 0.487_rkind) * tanh(246.78_rkind * (ta + tb - tsim - 0.0095_rkind));
        !    P0   = 0.209_rkind + (0.31814_rkind - 0.209_rkind) * tanh(277.79_rkind * (ta + tb - tsim - 0.0113_rkind));
        !    rho0 = 0.201_rkind + (0.27652_rkind - 0.201_rkind) * tanh(263.71_rkind * (ta + tb - tsim - 0.0110_rkind));
        !    T0   = 0.979_rkind + (1.15650_rkind - 0.979_rkind) * tanh(283.38_rkind * (ta + tb - tsim - 0.0145_rkind));
        !endif



        if(nrank == zero) then
        print*, "--values", tsim, U0, P0, rho0, T0
        endif
        ! set Dirichlet BC for velocity, density, pressure
        if(decomp%yst(1) == 1) then 
            !! left-boundary processor
            !u(1,:,:)   = u_left_bc
            do k = 1, decomp%ysz(3) 
               do j = 1, decomp%ysz(2)
                   rad = sqrt(y(1,j,1)**2 + z(1,1,k)**2)
                   u(1,j,k)   = U0 * (half - half * tanh(1.0_rkind  * (rad-1.0_rkind))) 
                   p(1,j,k)   = P0 * (half - half * tanh(1.0_rkind  * (rad-1.0_rkind))) 
                   rho(1,j,k) = rho0 * (half - half * tanh(1.0_rkind  * (rad-1.0_rkind))) 
                   T(1,j,k)   = T0 * (half - half * tanh(1.0_rkind  * (rad-1.0_rkind))) 
               enddo
            enddo
            v(1,:,:)   = zero
            w(1,:,:)   = zero
         endif


        ! Sponge+bulk for exit bc
        ! Gradually apply the exit boundary conditions
        !dy = Ly/real(decomp%ysz(2)-1,rkind)
        dx = Lx/real(decomp%xsz(1)-1,rkind)
        filpt = 2.0_rkind/dx 
        thickT = real(2.0D0, rkind)
        
        ! Gussian Filter for right side of domain 
        do i=1,decomp%ysz(1)
            dumT(i,:,:)=half*(one-tanh( (real(decomp%xsz(1)- (decomp%yst(1) - 1 + i - 1), rkind)-filpt) / thickT ))
        end do

        !write(outputfile, '(a,i3.3,a)') 'dumT_', nrank, '.dat'
        !open(10,file=outputfile,status='unknown')
        !do i=1,decomp%ysz(1)
        !   write(10,'(2(e19.12),1x)') x(i,1,1), dumT(i,1,1)
        !end do
        !close(10)

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

    end associate
end subroutine


subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim,sgsmodel)
!subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info, nrank
    use MixtureEOSMod,    only: mixture
    use sgsmod_cgrid,     only: sgs_cgrid
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use Shock_vortex_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid), optional,       intent(in) :: sgsmodel

    real(rkind) :: dx, Ythick, oob
    integer :: ny  , j
    integer :: iounit = 229
    character(len=clen) :: outputfile
    !real(rkind), dimension(decomp%ysz(2)) :: cmodel_loc, cmodel_loc_Qjsgs, cmodel_loc_tke

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))


    end associate

end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,   only: mixture

    use shock_vortex_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
