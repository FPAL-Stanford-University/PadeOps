module VDHIT_data
    use kind_parameters,  only: rkind
    use constants,        only: one,half,four
    implicit none
    
    real(rkind) :: Mt = 0.1_rkind, Re_lambda = 100._rkind, visc_exp = 0.75_rkind, gam = 1.4_rkind, k0 = four
    real(rkind) :: tke0, enstrophy0, dilatation_var0, u_rms0, lambda0, tau

    real(rkind) :: mu_ref
    real(rkind) :: T_ref = 1._rkind / 1.4_rkind
    real(rkind) :: Pr = 0.70_rkind
    logical     :: useRestart = .false.

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: two,pi
    use decomp_2d,        only: decomp_info

    use VDHIT_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)

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
    use mpi
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,three,pi
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info, nrank, nproc
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use PowerLawViscosityMod,        only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use exits,                       only: GracefulExit
    
    use VDHIT_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim,tstop,dt,tviz
    
    integer :: ioUnit
    integer :: counter, rank, ierr
    integer :: nx, ny, nz
    integer :: ix_dat, iy_dat, iz_dat, ix, iy, iz
    real(rkind) :: u_dat, v_dat, w_dat
    character(len=clen) :: datafile

    namelist /PROBINPUT/  Mt, Re_lambda, visc_exp, Pr, gam, k0, datafile, useRestart
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),                    &
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),                    &
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),                    &
                 e => fields(:,:,:,  e_index), Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        rho = one
        p   = one / gam
        u   = zero
        v   = zero
        w   = zero

        Ys(:,:,:,1) = half  ! equilibrium mass fraction of species 1
        Ys(:,:,:,2) = half  ! equilibrium mass fraction of species 2

        mu_ref = Mt * half / (Re_lambda * sqrt(three))
        T_ref = 1._rkind / gam

        ! Set the material property
        call mix%set_material(1, idealgas(gam,one), shearvisc = powerLawViscosity(mu_ref, T_ref, visc_exp), &
                                                    bulkvisc = constRatioBulkViscosity(zero), &
                                                    thermcond = constPrandtlConductivity(Pr) )

        if (.not. useRestart) then
            ioUnit = 27
            ! Make all processes read the file one by one to avoid file access conflicts
            do rank = 0,nproc-1
                call MPI_Barrier(MPI_COMM_WORLD, ierr)
                if (nrank == rank) then
                    open(unit=ioUnit, file=trim(datafile), form='FORMATTED')
                    read(ioUnit,*) ix_dat, iy_dat, iz_dat

                    if ((ix_dat /= nx) .or. (iy_dat /= ny) .or. (iz_dat /= nz)) then
                        call GracefulExit("Grid size in datafile does not match simulation grid size. Exiting...", 4652)
                    end if

                    counter = 0

                    ! Read in velocity data from file (u_rms = 1)
                    do while (counter < decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3))
                      read(ioUnit,*) ix_dat, iy_dat, iz_dat, u_dat, v_dat, w_dat
                      ix_dat = ix_dat + 1; iy_dat = iy_dat + 1; iz_dat = iz_dat + 1 ! Convert to 1 based indexing

                      if ( (ix_dat >= decomp%yst(1)) .and. (ix_dat <= decomp%yen(1)) ) then
                          if ( (iy_dat >= decomp%yst(2)) .and. (iy_dat <= decomp%yen(2)) ) then
                              if ( (iz_dat >= decomp%yst(3)) .and. (iz_dat <= decomp%yen(3)) ) then
                                  ix = ix_dat - decomp%yst(1) + 1
                                  iy = iy_dat - decomp%yst(2) + 1
                                  iz = iz_dat - decomp%yst(3) + 1
                                  u(ix,iy,iz) = u_dat
                                  v(ix,iy,iz) = v_dat
                                  w(ix,iy,iz) = w_dat
                                  counter = counter + 1
                              end if
                          end if
                      end if
                    end do
                    close(ioUnit)
                end if
            end do
        end if

        u_rms0 = Mt / sqrt(three)
        u = u_rms0 * u
        v = u_rms0 * v
        w = u_rms0 * w

        lambda0 = two / k0
        tau = lambda0 / u_rms0
        tstop = tstop * tau
        dt    = dt * tau
        tviz  = tviz * tau

    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info, nrank
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture
    use operators,        only: curl, divergence
    use reductions,       only: P_MEAN

    use VDHIT_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    
    integer                                     :: outputunit=229
    character(len=clen) :: outputfile, str
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tke, enstrophy, dilatation
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),3) :: vorticity
    real(rkind) :: tke_mean, enstrophy_mean, dilatation_var

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! Get TKE and Enstrophy to output to file
        tke = u*u + v*v + w*w
        call curl(decomp, der, u, v, w, vorticity)
        enstrophy = vorticity(:,:,:,1)**2 + vorticity(:,:,:,2)**2 + vorticity(:,:,:,3)**2

        call divergence(decomp, der, u, v, w, dilatation)

        write(str,'(I4.4)') decomp%ysz(2)
        write(outputfile,'(2A)') trim(outputdir),"/VDHIT_"//trim(str)//".dat"

        if (vizcount == 0) then
            tke0 = P_MEAN( tke )
            enstrophy0 = P_MEAN( vorticity(:,:,:,1)**2 + vorticity(:,:,:,2)**2 + vorticity(:,:,:,3)**2 )
            dilatation_var0 = P_MEAN( dilatation**2 )
            if (nrank == 0) then
                open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
                write(outputunit,'(4A26)') "Time", "Velocity Var", "Enstrophy", "Dilatation Var"
            end if
        else
            if (nrank == 0) then
                open(unit=outputunit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
            end if
        end if

        tke_mean = P_MEAN( tke )
        enstrophy_mean = P_MEAN( vorticity(:,:,:,1)**2 + vorticity(:,:,:,2)**2 + vorticity(:,:,:,3)**2 )
        dilatation_var = P_MEAN( dilatation**2 )

        if (nrank == 0) then
            ! write(outputunit,'(4ES26.16)') tsim, tke_mean/tke0, enstrophy_mean/(u_rms0**2 / lambda0**2), dilatation_var/(u_rms0**2 / lambda0**2)
            write(outputunit,'(4ES26.16)') tsim, tke_mean, enstrophy_mean, dilatation_var
            close(outputunit)
        end if
        
        ! write(str,'(I4.4,A,I4.4,A,I6.6)') decomp%ysz(2), "_", vizcount, "_", nrank
        ! write(outputfile,'(2A)') trim(outputdir),"/VDHIT_"//trim(str)//".dat"
        ! open(unit=outputunit, file=trim(outputfile), form='UNFORMATTED', status='REPLACE')
        ! write(outputunit) tsim
        ! write(outputunit) decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)
        ! write(outputunit) decomp%yst(1), decomp%yst(2), decomp%yst(3)
        ! write(outputunit) decomp%yen(1), decomp%yen(2), decomp%yen(3)
        ! write(outputunit) rho
        ! write(outputunit) u
        ! write(outputunit) v
        ! write(outputunit) w
        ! write(outputunit) vorticity(:,:,:,1)
        ! write(outputunit) vorticity(:,:,:,2)
        ! write(outputunit) vorticity(:,:,:,3)
        ! write(outputunit) dilatation
        ! write(outputunit) p
        ! close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture

    use VDHIT_data

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

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: half
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use exits,            only: message
    use reductions,       only: P_MAXVAL, P_MEAN

    use VDHIT_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    real(rkind) :: tke

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
       
        tke = P_MEAN( half * rho * (u*u + v*v + w*w) )

        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"TKE",tke/tke0)

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,rhsg)
    use kind_parameters,    only: rkind
    use decomp_2d,          only: decomp_info
    use MixtureEOSMod,      only: mixture
    use CompressibleGrid,   only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index

    use VDHIT_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), optional, intent(inout) ::rhsg

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),                    &
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),                    &
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),                    &
                 e => fields(:,:,:,  e_index), Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


             ! TODO: add focing schemes
    end associate
end subroutine
