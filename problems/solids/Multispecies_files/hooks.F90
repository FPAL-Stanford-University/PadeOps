module Multispecies_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: minVF = 0.2_rkind, thick = one
    real(rkind) :: rhoRatio = one
    logical     :: sharp = .FALSE.
    integer     :: tviz_counter = -1

end module

subroutine meshgen(decomp, dx, dy, dz, mesh, xcentered)
    use kind_parameters,  only: rkind
    use constants,        only: one, half, zero
    use decomp_2d,        only: decomp_info

    use Multispecies_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    logical,                         intent(in)    :: xcentered

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xfst

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,1)x[0,1)x[0,1) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = one/real(nx,rkind)
        dy = dx
        dz = dx

        if(xcentered) then
            xfst = half*dx
        else
            xfst = zero
        endif

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx + xfst
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two
    use SolidGrid,        only: u_index,v_index,w_index
    use decomp_2d,        only: decomp_info
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    
    use Multispecies_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, sharp, thick, minVF, rhoRatio
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if

        ! Set materials (both are the same material here)
        ! call mix%set_material(1,stiffgas(gamma,Rgas,p_infty),sep1solid(rho_0,zero,1.0D30,1.0D-10))
        call mix%set_material(1,stiffgas(gamma,Rgas,p_infty),sep1solid(rho_0,mu,1.0D30,1.0D-10))
        call mix%set_material(2,stiffgas(gamma,Rgas/rhoRatio,p_infty),sep1solid(rhoRatio*rho_0,mu,1.0D30,1.0D-10))

        u   = half
        v   = zero
        w   = zero

        if ( sharp ) then
            tmp = half * ( erf( (x-half+0.1_rkind)/(thick*dx) ) - erf( (x-half-0.1_rkind)/(thick*dx) ) )
        else
            tmp = exp(-((x-half)/0.1_rkind)**2)
        end if

        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        mix%material(1)%p  = p_amb
        mix%material(2)%p  = p_amb

        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        mix%material(2)%VF = one - mix%material(1)%VF

        tmp = rho_0*(mix%material(1)%VF+rhoRatio*(one-mix%material(1)%VF)) ! Mixture density
        mix%material(1)%Ys = mix%material(1)%VF * rho_0 / tmp
        mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity

    end associate

    tviz_counter = -1

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use DerivativesMod,   only: derivatives

    use Multispecies_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str,fname
    integer :: i
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, xnew
    real(rkind), dimension(decomp%ysz(1),6) :: exsoln
    real(rkind) :: err2, err3, err4, err5, err6, mat1vf, mat2vf, mat1ys, mat2ys, mxtrho

    tviz_counter = tviz_counter + 1

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (sharp) then
            write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2,A)') decomp%ysz(1), "_", minVF, "_", rhoRatio, "_sharp"
        else
            write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2,A)') decomp%ysz(1), "_", minVF, "_", rhoRatio, "_smooth"
        end if

        if (mix%use_gTg) then
            str = trim(str)//'_gTg'
        else
            str = trim(str)//'_g'
        end if

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Multispecies_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        write(outputunit,'(4ES27.16E3)') tsim, minVF, thick, rhoRatio
        do i=1,decomp%ysz(1)
            write(outputunit,'(23ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           mix%material(1)%p (i,1,1), mix%material(2)%p (i,1,1), &
                                           mix%material(1)%Ys(i,1,1), mix%material(2)%Ys(i,1,1), &
                                           mix%material(1)%VF(i,1,1), mix%material(2)%VF(i,1,1), &
                                           mix%material(1)%eh(i,1,1), mix%material(2)%eh(i,1,1), &
                                           mix%material(1)%T (i,1,1), mix%material(2)%T (i,1,1), &
                                           mix%material(1)%g11(i,1,1), mix%material(2)%g11(i,1,1), &
                                           mu(i,1,1), bulk(i,1,1), mix%material(1)%kap(i,1,1), mix%material(2)%kap(i,1,1), &
                                           mix%material(1)%diff(i,1,1), mix%material(2)%diff(i,1,1)
        end do
        close(outputunit)

        ! compute exact solution
        if ( sharp ) then
            tmp = half * ( erf( (x-half+0.1_rkind)/(thick*dx) ) - erf( (x-half-0.1_rkind)/(thick*dx) ) )
        else
            xnew = x-half-half*tsim
            where(xnew > 0.5d0)
              xnew = xnew - one
            endwhere
            where(xnew < -0.5d0)
              xnew = xnew + one
            endwhere
            tmp = exp(-(xnew/0.1_rkind)**2)
        end if

        err2 = 0.0d0; err3 = 0.0d0; err4 = 0.0d0; err5 = 0.0d0; err6 = 0.0d0
        do i=1,decomp%ysz(1)
!print *, i, xnew(i,1,1), tmp(i,1,1)
            mat1vf = minVF + (one-two*minVF)*tmp(i,1,1)
            mat2vf = one - mat1vf

            mxtrho = rho_0*(mat1vf+rhoRatio*(one-mat1vf)) ! Mixture density
            mat1ys = mat1vf * rho_0 / mxtrho
            mat2ys = one - mat1ys ! Enforce sum to unity

            exsoln(i,1) = x(i,1,1); exsoln(i,2) = mxtrho; 
            exsoln(i,3) = mat1vf;   exsoln(i,4) = mat2vf;  
            exsoln(i,5) = mat1ys;   exsoln(i,6) = mat2ys;  

            err2 = err2 + (exsoln(i,2) - rho(i,1,1))**2
            err3 = err3 + (exsoln(i,3) - mix%material(1)%VF(i,1,1))**2
            err4 = err4 + (exsoln(i,4) - mix%material(2)%VF(i,1,1))**2
            err5 = err5 + (exsoln(i,5) - mix%material(1)%Ys(i,1,1))**2
            err6 = err6 + (exsoln(i,6) - mix%material(2)%Ys(i,1,1))**2
        enddo
!stop
        err2 = sqrt(err2/decomp%ysz(1))
        err3 = sqrt(err3/decomp%ysz(1))
        err4 = sqrt(err4/decomp%ysz(1))
        err5 = sqrt(err5/decomp%ysz(1))
        err6 = sqrt(err6/decomp%ysz(1))

        if(tviz_counter==0) then
          open(10,file='exact_solution.dat',status='replace',action='write')
        else
          open(10,file='exact_solution.dat',status='unknown',action='write',position='append')
        endif
        write(10,'(a,e19.12,a,i5)') 'ZONE T="', tsim, '", F=POINT, I=', decomp%ysz(1)
        do i=1,decomp%ysz(1)
          write(10,'(7(e19.12,1x))') exsoln(i,:), tmp(i,1,1)
        enddo
        close(10)

        write(fname,'(a17,i4.4,a4)') 'errors_compiled_t', tviz_counter, '.dat'
        open(10,file=trim(fname),status='unknown',action='write',position='append')
        write(10,'(i5,1x,5(e19.12,1x))') decomp%ysz(1), err2, err3, err4, err5, err6
        close(10)
    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix

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
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL, P_MINVAL
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        call message("Min. VF1",P_MINVAL(mix%material(1)%VF))
        call message("Min. VF2",P_MINVAL(mix%material(2)%VF))
        call message("Max. VF1",P_MAXVAL(mix%material(1)%VF))
        call message("Max. VF2",P_MAXVAL(mix%material(2)%VF))
    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    type(solid_mixture),             intent(in)    :: mix

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_material_g_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine

subroutine hook_material_mass_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine

subroutine hook_material_energy_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine

subroutine hook_material_VF_source(decomp,hydro,elastic,x,y,z,tsim,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use Multispecies_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
