module obliqueimpact_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none
    
    real(rkind) :: omega = 1._rkind
    real(rkind) :: dxdyfixed = 1._rkind
    real(rkind) :: uvel = 300._rkind, vvel = 100.0_rkind
    real(rkind) :: pinit   = real(1.0D5,rkind)
    real(rkind) :: rho_0

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use obliqueimpact_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xa, xb, yc, yd

    xa = 0.0D0; xb = 2.0D0
    yc = 0.0D0; yd = 1.0D0

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = (xb-xa)/real(nx-1,rkind)
        dy = (yd-yc)/real(ny-1,rkind)
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = xa + real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = yc + real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,rho0,mu,yield,gam,PInf,tau0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    
    use obliqueimpact_data

    implicit none
    character(len=*),                                               intent(in)    :: inputfile
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind),                                          optional, intent(inout) :: rho0, mu, gam, PInf, tstop, dt, tviz, yield, tau0
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit, i, j, k
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: stp_x,bkstp_x,bkstp_y,phiang,u1,u2,v1,v2,rad,stp_r1,bkstp_r2,regfrac,dxdy

    namelist /PROBINPUT/  uvel, vvel, dxdyfixed
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    rho_0 = rho0

    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), g11 => fields(:,:,:,g11_index), &
               g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), & 
               g23 => fields(:,:,:,g23_index), g31 => fields(:,:,:,g31_index), & 
               g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        u1 = -uvel; u2 = uvel
        v1 = -vvel; v2 = vvel

        dxdy = sqrt(dx**2+dy**2)
        if(dxdyfixed>zero) dxdy = dxdyfixed
             
        do k=1,decomp%ysz(3)
         do j=1,decomp%ysz(2)
          do i=1,decomp%ysz(1)

             stp_r1 = half*(tanh(( x(i,j,k) - (1.0d0 + 0.1d0*cos(6.0d0*pi*y(i,j,k))))/dxdy) + one)
             regfrac = stp_r1

             u(i,j,k)   = (u1+regfrac*(u2-u1))
             v(i,j,k)   = (v1+regfrac*(v2-v1))

          end do 
         end do 
        end do 
        w   = zero
        p   = pinit

        g11 = one;  g12 = zero; g13 = zero
        g21 = zero; g22 = one;  g23 = zero
        g31 = zero; g32 = zero; g33 = one

        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp

    end associate

end subroutine

subroutine hook_output(decomp,der,fil,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index, &
                                sxx_index,sxy_index,sxz_index,syy_index,syz_index,szz_index
    use decomp_2d,        only: decomp_info, transpose_x_to_y, transpose_y_to_x
    use DerivativesMod,   only: derivatives
    use operators,        only: curl
    use FiltersMod,       only: filters

    use obliqueimpact_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(filters),                   intent(in) :: fil
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    real(rkind), dimension(2),       intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    real(rkind), allocatable, dimension(:,:,:,:) :: curlg

    character(len=clen) :: outputfile, velstr
    integer :: i,j,k

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3),       &
               sxx => fields(:,:,:, sxx_index), sxy => fields(:,:,:, sxy_index), sxz => fields(:,:,:, sxz_index), &
               syy => fields(:,:,:, syy_index), syz => fields(:,:,:, syz_index), szz => fields(:,:,:, szz_index)  )


        ! do post-processing stuff
        allocate(curlg(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),9))

        call curl(decomp, der, g11, g12, g13, curlg(:,:,:,1:3))
        call curl(decomp, der, g21, g22, g23, curlg(:,:,:,4:6))
        call curl(decomp, der, g31, g32, g33, curlg(:,:,:,7:9))


        write(velstr,'(I3.3)') int(1.0_rkind)
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/vort_"//trim(velstr)//"_", vizcount, ".dat"

        if(vizcount==0) then
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
          write(outputunit,'(290a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p","g11","g12","g13","g21","g22","g23","g31","g32","g33","sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar","curlg11","curlg12","curlg13","curlg21","curlg22","curlg23","curlg31","curlg32","curlg33"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(36ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9)
          
            end do
           end do
          end do
          close(outputunit)
        else
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(33ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9)
          
            end do
           end do
          end do
          close(outputunit)
        endif

        deallocate(curlg)
    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use obliqueimpact_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(2),       intent(in)    :: x_bc, y_bc, z_bc

    integer :: nx, ny

    nx = decomp%ysz(1); ny = decomp%ysz(2)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! left x
        rho( 1,:,:) = rho_0
        u  ( 1,:,:) = zero
        v  ( 1,:,:) = zero
        w  ( 1,:,:) = zero
        p  ( 1,:,:) = pinit
        
        g11( 1,:,:) = one;  g12( 1,:,:) = zero; g13( 1,:,:) = zero
        g21( 1,:,:) = zero; g22( 1,:,:) = one;  g23( 1,:,:) = zero
        g31( 1,:,:) = zero; g32( 1,:,:) = zero; g33( 1,:,:) = one

        ! right x
        rho(nx,:,:) = rho_0
        u  (nx,:,:) = zero
        v  (nx,:,:) = zero
        w  (nx,:,:) = zero
        p  (nx,:,:) = pinit
        
        g11(nx,:,:) = one;  g12(nx,:,:) = zero; g13(nx,:,:) = zero
        g21(nx,:,:) = zero; g22(nx,:,:) = one;  g23(nx,:,:) = zero
        g31(nx,:,:) = zero; g32(nx,:,:) = zero; g33(nx,:,:) = one

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL

    use obliqueimpact_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,tsim,rhs,rhsg)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhsg

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_finalize

end subroutine 
