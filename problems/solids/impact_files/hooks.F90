module impact_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    use FiltersMod,       only: filters
    implicit none
    
    real(rkind) :: uimpact = 100._rkind, vimpact = 0.0_rkind, wimpact = 0.0_rkind
    real(rkind) :: pinit   = real(1.0D5,rkind)
    real(rkind) :: rho_0
    real(rkind) :: thick = 1.1_rkind
    type(filters) :: mygfil

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use impact_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xa, xb

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)
    xa = -0.5_rkind;    xb = 1.5_rkind

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = (xb-xa)/real(nx-1,rkind)
        dy = dx
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx + xa
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,eostype,eosparams,rho0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    
    use impact_data

    implicit none
    character(len=*),                                               intent(in)    :: inputfile
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    integer,                                                        intent(in)    :: eostype
    real(rkind), dimension(:),                                      intent(inout) :: eosparams
    real(rkind),                                          optional, intent(inout) :: rho0, tstop, dt, tviz
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: mu, gam, PInf, yield, tau0

    namelist /PROBINPUT/  uimpact, vimpact, wimpact, pinit, thick
    
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
       
        if(thick<zero) then 
            tmp = tanh( (x-half)/(abs(thick)) )
        else
            tmp = tanh( (x-half)/(thick*dx) )
        endif

        u   = -uimpact*tmp
        v   = -vimpact*tmp
        w   = -wimpact*tmp
        p   = pinit

        g11 = one;  g12 = zero; g13 = zero
        g21 = zero; g22 = one;  g23 = zero
        g31 = zero; g32 = zero; g33 = one

        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp

        if(eostype==3) then
          T = pinit/(eosparams(7)*rho0*eosparams(5)) + eosparams(6) ! = pinit/(gamma*rho0*CV) + T0
        endif

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    end associate

end subroutine

subroutine hook_output(decomp,der,fil,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index, &
                                sxx_index,sxy_index,sxz_index,syy_index,syz_index,szz_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters

    use impact_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(derivatives),               intent(in) :: der
    type(filters),                   intent(in) :: fil
    real(rkind), dimension(2),       intent(in) :: x_bc,y_bc,z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, velstr
    integer :: i,j,k
    integer :: indx(1), nx, indhalf, numshocks
    real(rkind) :: xshock(2), betmax
    real(rkind), allocatable, dimension(:) :: bettmp(:)

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


        write(velstr,'(I3.3)') int(uimpact)
        write(*,*) uimpact
        write(*,*) trim(outputdir)
        write(*,*) "/tec_impact_"//trim(velstr)//".dat"
        write(outputfile,'(2A)') trim(outputdir),"/tec_impact_"//trim(velstr)//".dat"

        if(vizcount==0) then
          open(unit=outputunit, file=trim(outputfile), form='formatted', status='unknown')
          write(outputunit,'(200a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p","g11","g12","g13","g21","g22","g23","g31","g32","g33","sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar","T"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(28ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), T(i,j,k)
          
            end do
           end do
          end do
          close(outputunit)
        else
          open(unit=outputunit, file=trim(outputfile), form='formatted', status='old', action='write',position='append')
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(25ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), T(i,j,k)
          
            end do
           end do
          end do
          close(outputunit)
        endif

        !! determine shock speeds
        !nx = decomp%ysz(1); indhalf = nx/2
        !allocate(bettmp(nx))
        !bettmp = bulk(:,1,1)
        !betmax = maxval(bettmp(indhalf:nx));  
        !do i = nx, nx/2, -1
        !  if(bettmp(i)<0.01d0*betmax) bettmp(i) = zero
        !enddo
        !numshocks = 0; xshock = 0.5d0
        !do i = nx/2+1, nx-1
        !  if((bettmp(i-1) < bettmp(i)) .and. (bettmp(i)>bettmp(i+1))) then
        !    numshocks = numshocks+1
        !    xshock(numshocks) = x(i,1,1)
        !    if(numshocks>2) write(*,*) 'More than 2 shocks detected. Check details.'
        !  endif
        !enddo
        !write(111,*) tsim, xshock(1), xshock(2)
        !deallocate(bettmp)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one, half
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use impact_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc

    integer :: nx, i
    real(rkind) :: xspng, tspng, dx
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        !rho( 1,:,:) = rho_0
        !u  ( 1,:,:) = uimpact
        !v  ( 1,:,:) = vimpact
        !w  ( 1,:,:) = wimpact
        !p  ( 1,:,:) = pinit
        !
        !g11( 1,:,:) = one;  g12( 1,:,:) = zero; g13( 1,:,:) = zero
        !g21( 1,:,:) = zero; g22( 1,:,:) = one;  g23( 1,:,:) = zero
        !g31( 1,:,:) = zero; g32( 1,:,:) = zero; g33( 1,:,:) = one

        !rho(nx,:,:) = rho_0
        !u  (nx,:,:) = -uimpact
        !v  (nx,:,:) = -vimpact
        !w  (nx,:,:) = -wimpact
        !p  (nx,:,:) = pinit
        !
        !g11(nx,:,:) = one;  g12(nx,:,:) = zero; g13(nx,:,:) = zero
        !g21(nx,:,:) = zero; g22(nx,:,:) = one;  g23(nx,:,:) = zero
        !g31(nx,:,:) = zero; g32(nx,:,:) = zero; g33(nx,:,:) = one

        rho( 1,:,:) = rho(2,:,:)
        u  ( 1,:,:) = u(2,:,:)
        v  ( 1,:,:) = v(2,:,:)
        w  ( 1,:,:) = w(2,:,:)
        p  ( 1,:,:) = p(2,:,:)
        
        g11( 1,:,:) = g11(2,:,:);  g12( 1,:,:) = g12(2,:,:);  g13( 1,:,:) = g13(2,:,:)
        g21( 1,:,:) = g21(2,:,:);  g22( 1,:,:) = g22(2,:,:);  g23( 1,:,:) = g23(2,:,:)
        g31( 1,:,:) = g31(2,:,:);  g32( 1,:,:) = g32(2,:,:);  g33( 1,:,:) = g33(2,:,:)

        rho(nx,:,:) = rho(nx-1,:,:)
        u  (nx,:,:) = u(nx-1,:,:)
        v  (nx,:,:) = v(nx-1,:,:)
        w  (nx,:,:) = w(nx-1,:,:)
        p  (nx,:,:) = p(nx-1,:,:)
        
        g11(nx,:,:) = g11(nx-1,:,:);  g12(nx,:,:) = g12(nx-1,:,:);  g13(nx,:,:) = g13(nx-1,:,:)
        g21(nx,:,:) = g21(nx-1,:,:);  g22(nx,:,:) = g22(nx-1,:,:);  g23(nx,:,:) = g23(nx-1,:,:)
        g31(nx,:,:) = g31(nx-1,:,:);  g32(nx,:,:) = g32(nx-1,:,:);  g33(nx,:,:) = g33(nx-1,:,:)

        ! sponge
        dx = x(2,1,1)-x(1,1,1)
        xspng = x(nx,1,1)-0.1_rkind;  tspng = 0.15_rkind
        dum = half*(one - tanh((x-xspng)/tspng))          ! backstep at x2
        xspng = x(1,1,1)+0.1_rkind;  tspng = 0.15_rkind
        dum = dum*half*(one + tanh((x-xspng)/tspng))      ! step at x1
        dum = one - dum
        !do i = 1, nx
        !  write(*,*) x(i,1,1), dum(i,1,1)
        !enddo
        !stop

        do i=1,4
            tmp = u
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            u = u + dum*(tmp - u)

            tmp = v
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            v = v + dum*(tmp - v)

            tmp = w
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            w = w + dum*(tmp - w)

            tmp = e
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            e = e + dum*(tmp - e)

            tmp = rho
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            rho = rho + dum*(tmp - rho)

            tmp = p
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            p = p + dum*(tmp - p)

            tmp = g11(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g11(:,:,:) = g11(:,:,:) + dum*(tmp - g11(:,:,:))

            tmp = g12(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g12(:,:,:) = g12(:,:,:) + dum*(tmp - g12(:,:,:))

            tmp = g13(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g13(:,:,:) = g13(:,:,:) + dum*(tmp - g13(:,:,:))

            tmp = g21(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g21(:,:,:) = g21(:,:,:) + dum*(tmp - g21(:,:,:))

            tmp = g22(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g22(:,:,:) = g22(:,:,:) + dum*(tmp - g22(:,:,:))

            tmp = g23(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g23(:,:,:) = g23(:,:,:) + dum*(tmp - g23(:,:,:))

            tmp = g31(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g31(:,:,:) = g31(:,:,:) + dum*(tmp - g31(:,:,:))

            tmp = g32(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g32(:,:,:) = g32(:,:,:) + dum*(tmp - g32(:,:,:))

            tmp = g33(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g33(:,:,:) = g33(:,:,:) + dum*(tmp - g33(:,:,:))
        end do

    end associate
end subroutine

subroutine filter(decomp,myfil,arr,numtimes,x_bc_,y_bc_,z_bc_)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use FiltersMod,       only: filters

    implicit none
    type(decomp_info),                                                    intent(in) :: decomp
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(inout) :: arr
    type(filters),                                                        intent(in) :: myfil
    integer,                                                              intent(in) :: numtimes
    integer, dimension(2),                                                intent(in) :: x_bc_, y_bc_, z_bc_

    integer :: times2fil
    integer, dimension(2) :: x_bc, y_bc, z_bc
    real(rkind), dimension(:,:,:), pointer :: tmp_in_y, tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
    real(rkind), dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),2), target :: xbuf
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),1), target :: ybuf
    real(rkind), dimension(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),2), target :: zbuf
    integer :: lastx, lasty, lastz, idx


    !if (present(myfil)) then
    !    fil2use => myfil
    !else
    !    !fil2use => this%fil
    !end if 

    !if (present(numtimes)) then
        times2fil = numtimes
    !else
    !    times2fil = 1
    !end if

    ! Allocate pointers for the needed buffers 
    ! Atleast 2 buffers in x and z are assumed
    ! Last two buffers are occupied

    !lastx = size(this%xbuf,4)
    !lasty = size(this%ybuf,4)
    !lastz = size(this%zbuf,4)

    tmp1_in_x => xbuf(:,:,:,1)
    tmp2_in_x => xbuf(:,:,:,2)
    tmp_in_y => ybuf(:,:,:,1)
    tmp1_in_z => zbuf(:,:,:,1)
    tmp2_in_z => zbuf(:,:,:,2)
   
    x_bc = x_bc_
    y_bc = y_bc_
    z_bc = z_bc_
   

    !! Order y -> x 
    !  ! First filter in y
    !  call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
    !  ! Subsequent refilters 
    !  do idx = 1,times2fil-1
    !      arr = tmp_in_y
    !      call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
    !  end do
    !  
    !  ! Then transpose to x
    !  call transpose_y_to_x(tmp_in_y,tmp1_in_x,decomp)

    !  ! First filter in x
    !  call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
    !  ! Subsequent refilters
    !  do idx = 1,times2fil-1
    !      tmp1_in_x = tmp2_in_x
    !      call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
    !  end do 

    !  ! Now transpose back to y
    !  call transpose_x_to_y(tmp2_in_x,tmp_in_y,decomp)
    !! ---------Order y -> x 

    ! Order x -> y 
      ! Transpose to x
      call transpose_y_to_x(arr,tmp1_in_x,decomp)

      ! First filter in x
      call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
      ! Subsequent refilters
      do idx = 1,times2fil-1
          tmp1_in_x = tmp2_in_x
          call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
      end do 

      ! Now transpose back to y
      call transpose_x_to_y(tmp2_in_x,arr,decomp)

      ! First filter in y
      call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
      ! Subsequent refilters 
      do idx = 1,times2fil-1
          arr = tmp_in_y
          call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
      end do
    ! ---------Order x -> y 
      

    ! Now transpose to z
    call transpose_y_to_z(tmp_in_y,tmp1_in_z,decomp)

    !First filter in z
    call myfil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
    ! Subsequent refilters
    do idx = 1,times2fil-1
        tmp1_in_z = tmp2_in_z
        call myfil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
    end do 

    ! Now transpose back to y
    call transpose_z_to_y(tmp2_in_z,arr,decomp)

    ! Finished
end subroutine

subroutine hook_timestep(decomp,der,mesh,fields,step,tsim,dt,x_bc,y_bc,z_bc,hookcond)
    use kind_parameters,  only: rkind
    use SolidGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use exits,            only: message
    use reductions,       only: P_MAXVAL

    use impact_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim, dt
    real(rkind), dimension(2),       intent(in) :: x_bc,y_bc,z_bc
    logical,            optional, intent(inout) :: hookcond

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

subroutine hook_finalize()

end subroutine

