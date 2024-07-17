module shocktube_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none
    
    integer     :: ns     = 1
    real(rkind) :: Lx     = 3.0_rkind
    real(rkind) :: Ly     = 4.0_rkind
    real(rkind) :: Lz     = 0.01_rkind
    real(rkind) :: Pr     = 0.7_rkind 
    real(rkind) :: Sc     = 1.0_rkind  
    real(rkind) :: gam    = 1.4_rkind
    real(rkind) :: Rgas   = 287.053_rkind
    real(rkind) :: mu_ref = 0.00001849_rkind 
    real(rkind) :: pLef   = 1013250.0_rkind
    real(rkind) :: pRig   = 101325.0_rkind
    real(rkind) :: rhoLef = 11.851_rkind
    real(rkind) :: rhoRig = 1.1851_rkind
    real(rkind) :: x1, y1, z1

end module

!subroutine meshgen(decomp, dx, dy, dz, mesh)
subroutine meshgen(decomp, dx, dy, dz, mesh, inputfile, xmetric, ymetric, zmetric, xi, eta, zeta, dxs, dys, dzs, xbuf, zbuf)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y

    use shocktube_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn, iounit
    character(len=*),                intent(in)    :: inputfile
    logical,                         intent(in   ) :: xmetric, ymetric, zmetric
    real(rkind), dimension(:,:,:  ), intent(inout) :: xi, eta, zeta
    real(rkind), dimension(:      ), intent(inout) :: dxs, dys, dzs
    real(rkind), dimension(:,:,:,:), target,intent(in):: xbuf, zbuf
    real(rkind), dimension(:,:,:), pointer :: xtmp1, xtmp2, ztmp1, ztmp2

    namelist /PROBINPUT/ ns, Lx, Ly, Lz, Pr, Sc, gam, Rgas, mu_ref, pLef, pRig, rhoLef, rhoRig

    ioUnit = 15
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    xtmp1 => xbuf(:,:,:,1); xtmp2 => xbuf(:,:,:,2) 
    ztmp1 => zbuf(:,:,:,1); ztmp2 => zbuf(:,:,:,2) 

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nx-1,rkind)
        dy = Ly/real(ny-1,rkind)
        dz = dx
        
        x1 = -0.285_rkind   ! location of diaphragm
        !print*, nx, ny,nz,ix1,ixn,iy1,iyn,iz1,izn

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx + x1
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do
    !For now let's consider it as a uniform mesh..
    dxs = dx
    dys = dy
    dzs = dz
   
    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,&
                                           p_index,T_index,e_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use IdealGasEOS,                 only: idealgas
    use SutherLandViscosityMod,      only: sutherlandViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use SutherLandConductivityMod,   only: SutherLandConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use exits,                       only: GracefulExit, message, nancheck
 
    use shocktube_data

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
    real(rkind) :: S, Sk, T0
    integer :: i,j, k, iounit, nx, ny, nz, nxl, nyl, nzl
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp

    namelist /PROBINPUT/ ns, Lx, Ly, Lz, Pr, Sc, gam, Rgas, mu_ref, pLef, pRig, rhoLef, rhoRig

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), &
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)
        T0     =  pRig/(Rgas*rhoRig)        
        S      = 110.4_rkind
        Sk     = 194.0_rkind
        
        !print*, T0, mu_ref, Pr, S, Sk, Rgas
        !!!! Set each material's transport coefficient object
        shearvisc = sutherlandViscosity(mu_ref, T0, 1.5_rkind, S)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = sutherlandConductivity(Pr, T0, 1.5_rkind, Sk)
        call mix%set_material( 1, idealgas( gam, Rgas ), shearvisc = shearvisc, bulkvisc  = bulkvisc, thermcond = thermcond  )

        Ys(:,:,:,1)  = one 
        call mix%update(Ys)     
     
        ! Initialise the flow field
        tmp = half * ( one+tanh(x/(dx)) )

        rho = (one-tmp)*rhoLef + tmp*rhoRig
        u   = zero
        v   = zero
        w   = zero
        p   = (one-tmp)*pLef + tmp*pRig
        T   = p/(rho*Rgas) 
           
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use shocktube_data

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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/shocktube_", vizcount, ".dat"

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
    use constants,        only: zero
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture

    use shocktube_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
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

    use shocktube_data

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

    use shocktube_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
