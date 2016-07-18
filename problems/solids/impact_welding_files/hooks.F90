module impact_welding_data
    use kind_parameters,  only: rkind
    use constants,        only: zero,one,two,eight,three,six
    use FiltersMod,       only: filters
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 10._rkind, rho_0_2 = one
    real(rkind) :: minVF = real(1.0D-8,rkind), thick = one
    real(rkind) :: rhoRatio = one, pRatio = two
    real(rkind) :: uimpact = one, theta = 30._rkind
    logical     :: sharp = .FALSE.
    real(rkind) :: rhoL, rhoR, YsL, YsR, VFL, VFR
    real(rkind) :: yield = one, yield2 = one, eta0k = 0.4_rkind
    logical     :: explPlast = .FALSE., explPlast2 = .FALSE.
    logical     :: plastic = .FALSE., plastic2 = .FALSE.
    real(rkind) :: Ly = one, Lx = six, interface_init = zero, kwave = 4.0_rkind

    type(filters) :: mygfil

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: three
    use decomp_2d,        only: decomp_info
    use exits,            only: warning

    use impact_welding_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,1)x[0,1)x[0,1) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny,rkind)
        dz = dx

        if(abs(dx-dy)>1.0d-13) then
          call warning("dx not equal to dy")
        endif

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1     + i - 1, rkind ) * dx - three  ! x \in (-3,3]
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi
    use SolidGrid,        only: u_index,v_index,w_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    
    use impact_welding_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    real(rkind), dimension(8) :: fparams
    integer, dimension(2) :: iparams

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rhoRatio, pRatio, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, plastic, explPlast, yield,   &
                          plastic2, explPlast2, yield2, interface_init, kwave, eta0k, uimpact, theta
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if

        if(rhoRatio > 0) then
          ! if rhoRatio is positive, only rho_0 is different. Rgas is set such
          ! that Temperature equilibrium condition is satisfied
          gamma_2 = gamma; Rgas_2 = Rgas/rhoRatio; p_infty_2 = p_infty; 
          rho_0_2 = rho_0*rhoRatio; mu_2 = mu
        else
          ! if rhoRatio is negative, all quantities except Rgas need to be
          ! specified in input file. Rgas is then set such
          ! that Temperature equilibrium condition is satisfied
          Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2
        endif

        ! write material properties
        if (nrank == 0) then
            print *, '---Material 1---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0, ', gam  = ', gamma, ', p_infty = ', p_infty
            write(*,'(3(a,e12.5))') 'mu    = ', mu,    ', Rgas = ', Rgas
            print *, '---Material 2---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0_2, ', gam  = ', gamma_2, ', p_infty = ', p_infty_2
            write(*,'(3(a,e12.5))') 'mu    = ', mu_2,    ', Rgas = ', Rgas_2
        end if

        ! Set materials
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield, 1.0D-10))
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10))

        ! set logicals for plasticity
        mix%material(1)%plast = plastic;  mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        tmp = half*(one + erf( (x-(interface_init+eta0k/(2.0_rkind*pi*kwave)*sin(2.0_rkind*kwave*pi*y)))/(thick*dx) )) - rho_0_2/(rho_0 + rho_0_2)

        theta = theta * pi / 180._rkind ! Convert angle from degrees to radians

        u   =-uimpact*cos(theta)*tmp
        v   =-uimpact*sin(theta)*tmp
        w   = zero

        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one
        
        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        mix%material(1)%p  = zero
        mix%material(2)%p  = mix%material(1)%p

        tmp = half*(one - erf( (x-(interface_init+eta0k/(2.0_rkind*pi*kwave)*sin(2.0_rkind*kwave*pi*y)))/(thick*dx) ))
        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        mix%material(2)%VF = one - mix%material(1)%VF

        tmp = rho_0*mix%material(1)%VF*mix%material(1)%g11 + rho_0_2*mix%material(2)%g11*(one-mix%material(1)%VF) ! Mixture density
        mix%material(1)%Ys = mix%material(1)%VF * rho_0 / tmp
        mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity

        rhoL = tmp(1,1,1)
        rhoR = tmp(decomp%ysz(1),1,1)
        YsL  = mix%material(1)%Ys(1,1,1)
        YsR  = mix%material(1)%Ys(decomp%ysz(1),1,1)
        VFL  = mix%material(1)%VF(1,1,1)
        VFR  = mix%material(1)%VF(decomp%ysz(1),1,1)

    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,four,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                sxx_index,syy_index,szz_index,sxy_index,sxz_index,syz_index,sos_index
    use decomp_2d,        only: decomp_info, nrank
    use DerivativesMod,   only: derivatives
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: curl
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL

    use impact_welding_data

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

    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),3) :: vort
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)  ) :: tmp
    real(rkind), dimension(decomp%ysz(1)) :: Ys1_mean, Ys2_mean
    real(rkind) :: vort_pos, vort_neg, mixwidth, Al_mass, xspike, xbubbl, xspike_proc, xbubbl_proc
    character(len=clen) :: outputfile, str
    integer :: i, j, k

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 sxx  => fields(:,:,:, sxx_index), syy => fields(:,:,:,syy_index), &
                 szz  => fields(:,:,:, szz_index), sxy => fields(:,:,:,sxy_index), &
                 sxz  => fields(:,:,:, sxz_index), syz => fields(:,:,:,syz_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3),       &
                 sos  => fields(:,:,:,sos_index) )

       if (rhoRatio > 0) then
           write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2)') nrank, "_", minVF, "_", rhoRatio
       else
           write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2)') nrank, "_", minVF, "_", rho_0_2/rho_0
       end if

       if (decomp%ysz(2) == 1) then
           write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/impact_welding_"//trim(str)//"_", vizcount, ".dat"

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
       end if

       call curl(decomp, der, u, v, w, vort, x_bc, y_bc, z_bc)
       
       tmp = zero
       where (vort(:,:,:,3) .GE. zero)
           tmp = vort(:,:,:,3)
       end where
       vort_pos = P_MEAN(tmp)*six*one
       
       tmp = zero
       where (vort(:,:,:,3) .LE. zero)
           tmp = vort(:,:,:,3)
       end where
       vort_neg = P_MEAN(tmp)*six*one

       Ys1_mean = SUM(mix%material(1)%Ys(:,:,1),2) / real(decomp%ysz(2),rkind)
       Ys2_mean = SUM(mix%material(2)%Ys(:,:,1),2) / real(decomp%ysz(2),rkind)

       Ys1_mean = four * Ys1_mean * Ys2_mean
       mixwidth = P_SUM(Ys1_mean) * dx

       Al_mass = P_MEAN(rho*mix%material(2)%Ys)*six*one

       ! Get bubble and spike extents
       xspike_proc = -two
       xbubbl_proc = four
       do j = 1,decomp%ysz(2)
           do i = 1,decomp%ysz(1)
               if (mix%material(1)%Ys(i,j,1) .GE. half) xspike_proc = max(xspike_proc,x(i,j,1))
               if (mix%material(1)%Ys(i,j,1) .LE. half) xbubbl_proc = min(xbubbl_proc,x(i,j,1))
           end do
       end do
       xspike = P_MAXVAL(xspike_proc)
       xbubbl = P_MINVAL(xbubbl_proc)
           
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/impact_welding_statistics.dat"

       if (vizcount == 0) then
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
           write(outputunit,'(7A27)') 'tsim', 'mixwidth', 'vort_pos', 'vort_neg', 'Al_mass', 'xspike', 'xbubbl'
       else
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', action='WRITE', status='OLD', position='APPEND')
       end if
       write(outputunit,'(7ES27.16E3)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl
       close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use impact_welding_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: nx, ax, i, j
    real(rkind) :: dx, xspng, tspng
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    
    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if(decomp%yst(1)==1) then
          if(x_bc(1)==0) then
              rho( 1,:,:) = rho( 2,:,:)
              u  ( 1,:,:) = u  ( 2,:,:)
              v  ( 1,:,:) = v  ( 2,:,:)
              w  ( 1,:,:) = zero
              mix%material(1)%p( 1,:,:) = mix%material(1)%p(2,:,:)
              mix%material(2)%p( 1,:,:) = mix%material(2)%p(2,:,:)
             
              mix%material(1)%g( 1,:,:,:) = mix%material(1)%g(2,:,:,:)
              mix%material(2)%g( 1,:,:,:) = mix%material(2)%g(2,:,:,:)
              
              ! mix%material(1)%Ys ( 1,:,:) = YsL
              ! mix%material(2)%Ys ( 1,:,:) = one - YsL
  
              mix%material(1)%VF ( 1,:,:) = mix%material(1)%VF ( 2,:,:)
              mix%material(2)%VF ( 1,:,:) = mix%material(2)%VF ( 2,:,:)
          end if
        endif

        xspng = -2.75_rkind
        tspng = 0.4_rkind
        dx = x(2,1,1) - x(1,1,1)
        dum = half*(one - tanh( (x-xspng)/(tspng) ))

        do i=1,4
            tmp = u
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            u = u + dum*(tmp - u)

            tmp = v
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            v = v + dum*(tmp - v)

            tmp = w
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            w = w + dum*(tmp - w)

            tmp = e
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            e = e + dum*(tmp - e)

            tmp = rho
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            rho = rho + dum*(tmp - rho)

            tmp = mix%material(1)%p
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            mix%material(1)%p = mix%material(1)%p + dum*(tmp - mix%material(1)%p)

            tmp = mix%material(2)%p
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            mix%material(2)%p = mix%material(2)%p + dum*(tmp - mix%material(2)%p)

            do j = 1,9
                tmp = mix%material(1)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

                tmp = mix%material(2)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(2)%g(:,:,:,j) = mix%material(2)%g(:,:,:,j) + dum*(tmp - mix%material(2)%g(:,:,:,j))
            end do
        end do

        nx = decomp%xsz(1)
        ax = decomp%ysz(1)
        if(decomp%yen(1)==nx) then
          if(x_bc(1)==0) then
              rho(ax,:,:) = rho(ax-1,:,:)
              u  (ax,:,:) = u  (ax-1,:,:)
              v  (ax,:,:) = v  (ax-1,:,:)
              w  (ax,:,:) = zero
              mix%material(1)%p(ax,:,:) = mix%material(1)%p(ax-1,:,:)
              mix%material(2)%p(ax,:,:) = mix%material(2)%p(ax-1,:,:)
       
              mix%material(1)%g(ax,:,:,:) = mix%material(1)%g(ax-1,:,:,:)
              mix%material(2)%g(ax,:,:,:) = mix%material(2)%g(ax-1,:,:,:)
  
              mix%material(1)%VF (ax,:,:) = mix%material(1)%VF (ax-1,:,:)
              mix%material(2)%VF (ax,:,:) = mix%material(2)%VF (ax-1,:,:)
          end if
        endif

        xspng = 2.75_rkind
        tspng = 0.4_rkind
        dx = x(2,1,1) - x(1,1,1)
        dum = half*(one + tanh( (x-xspng)/(tspng) ))

        do i=1,4
            tmp = u
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            u = u + dum*(tmp - u)

            tmp = v
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            v = v + dum*(tmp - v)

            tmp = w
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            w = w + dum*(tmp - w)

            tmp = e
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            e = e + dum*(tmp - e)

            tmp = rho
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            rho = rho + dum*(tmp - rho)

            tmp = mix%material(1)%p
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            mix%material(1)%p = mix%material(1)%p + dum*(tmp - mix%material(1)%p)

            tmp = mix%material(2)%p
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            mix%material(2)%p = mix%material(2)%p + dum*(tmp - mix%material(2)%p)

            do j = 1,9
                tmp = mix%material(1)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

                tmp = mix%material(2)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(2)%g(:,:,:,j) = mix%material(2)%g(:,:,:,j) + dum*(tmp - mix%material(2)%g(:,:,:,j))
            end do
        end do

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use impact_welding_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer                                     :: imin, ind(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use impact_welding_data

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

    use impact_welding_data

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

    use impact_welding_data

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

    use impact_welding_data

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

    use impact_welding_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
