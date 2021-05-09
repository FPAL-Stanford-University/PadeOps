module NCtemplate_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, one
    use FiltersMod,       only: filters

    implicit none

    integer, parameter :: ns = 2

    ! Problem parameters

    ! Domain size data
    real(rkind) :: L_z = 1.0_rkind
    real(rkind) :: L_x = 1.0_rkind, L_y = 1.0_rkind
    real(rkind) :: x1, y1, z1, tmp_input

    ! Filtering mask parameters (for nonperiodic boundaries)
    logical :: useSponge = .false.            ! Use sponge region at boundaries?
    real(rkind) :: spongesize = 0.1_rkind     ! Fraction of domain used for filtering mask
    real(rkind) :: spongewidth = 0.1_rkind    ! Transition width for the filtering mask

    logical :: periodicx = .true., periodicy = .true., periodicz = .false.

    ! Gaussian filter for sponge
    type(filters) :: mygfil

contains


    subroutine get_volumefractions(mix,Ys,Xs)
        use MixtureEOSMod,               only: mixture
        type(mixture),                   intent(in)  :: mix
        real(rkind), dimension(:,:,:,:), intent(in)  :: Ys    ! Species mass fractions
        real(rkind), dimension(size(Ys,1), size(Ys,2), size(Ys,3), size(Ys,4)), intent(out) :: Xs    ! Species volume fractions
        
        integer :: n

        ! Get volume fractons
        do n = 1,mix%ns
            Xs(:,:,:,n) = mix%material(n)%mat%Rgas * Ys(:,:,:,n) / mix%Rgas
        end do
   
    end subroutine

end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two
    use decomp_2d,        only: decomp_info

    use NCtemplate_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k,ioUnit
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        dx = L_x/real(nx  ,rkind)
        dy = L_y/real(ny  ,rkind)
        dz = L_z/real(nz-1,rkind)


        x1 = zero
        y1 = zero
        z1 = zero

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do


    end associate

end subroutine


subroutine initfields(decomp,der,dx,dy,dz,inputfile,mesh,fields,Wbackground,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,pi,eight
    use CompressibleGridNC,          only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index,mom_index,TE_index
    use decomp_2d,                   only: decomp_info,transpose_x_to_y,transpose_y_to_x,transpose_y_to_z,transpose_z_to_y
    use DerivativesMod,              only: derivatives
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use PowerLawViscosityMod,        only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use io_hdf5_stuff,               only: io_hdf5
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use mpi
    
    use NCtemplate_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: Wbackground
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    type(powerLawViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond
    integer :: i, j, k, l, mode, iounit
    integer :: nx, ny, nz

    namelist /PROBINPUT/ tmp_input
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,u_index),  &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,w_index),  &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,T_index),  &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:,mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap =>fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),       &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)


    end associate

end subroutine


subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use NCtemplate_data

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
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/rayleightaylor_", vizcount, ".dat"

    end associate
end subroutine


subroutine hook_filmask(decomp,mesh,Wcnsrv,mask,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, pi
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info

    use NCtemplate_data

    implicit none
    type(decomp_info),               intent(in)      :: decomp
    real(rkind),                     intent(in)      :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)      :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)      :: Wcnsrv
    real(rkind), dimension(:,:,:),   intent(inout)   :: mask
    integer, dimension(2),           intent(in)      :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: tmp
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        ! makes the mask function
        tmp = one-half*erf((z-L_z*spongesize)/spongewidth)+half*erf((z-L_z+L_z*spongesize)/spongewidth)
        where (tmp < 0.00001_rkind)
            mask = zero
        elsewhere
            mask = tmp
        end where
   end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, pi
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use NCtemplate_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind) :: dx
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


    end associate
end subroutine


subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two,eps
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use NCtemplate_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    integer :: ny
    integer :: iounit = 229
    character(len=clen) :: outputfile

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))

    end associate
end subroutine


subroutine hook_source(decomp,der,mesh,fields,Wcnsrv,Wbackground,mix,tsim,dt,rhs,rhsg)
    use kind_parameters,    only: rkind
    use constants,          only: zero,half,one,two,four
    use decomp_2d,          only: decomp_info,transpose_x_to_y,transpose_y_to_x,transpose_y_to_z,transpose_z_to_y
    use DerivativesMod,     only: derivatives
    use MixtureEOSMod,      only: mixture
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index,mom_index,TE_index

    use NCtemplate_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind),                     intent(in)    :: dt
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(in)    :: Wcnsrv
    real(rkind), dimension(:,:,:,:), intent(in)    :: Wbackground
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), optional, intent(inout) ::rhsg
    
    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

    end associate
end subroutine
