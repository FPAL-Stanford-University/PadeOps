module ibmgpmod
    use kind_parameters, only: rkind, clen
    use constants      , only: zero, one, third, half, two, kappa
    use exits          , only: GracefulExit
    use reductions     , only: p_maxval, p_minval
    use kdtree_wrapper , only: initialize_kdtree_ib, create_kdtree_ib, probe_nearest_points_kdtree_ib, finalize_kdtree_ib
    use decomp_2d

    implicit none

    private
    public :: ibmgp

    type :: ibmgp
        private 
        class(decomp_info), pointer :: gpC, gpE
        !class(spectral), pointer :: spectC, spectE

        integer     :: num_surfelem, num_gptsC, num_gptsE, ibwm
        real(rkind) :: ibwm_ustar
        real(rkind), allocatable, dimension(:,:,:) :: surfelem
        real(rkind), allocatable, dimension(:,:)   :: surfcent, surfnormal
        real(rkind), allocatable, dimension(:,:)   :: gptsC_xyz, gptsE_xyz, gptsC_bpt, gptsE_bpt
        real(rkind), allocatable, dimension(:,:)   :: gptsC_img, gptsE_img, gptsC_bnp, gptsE_bnp
        integer,     allocatable, dimension(:,:)   :: gptsC_ind, gptsE_ind
        integer,     allocatable, dimension(:)     :: gptsC_bpind, gptsE_bpind!, gptsC_ileft, gptsC_jleft, gptsC_kleft
        !integer,     allocatable, dimension(:)     :: gptsE_ileft, gptsE_jleft, gptsE_kleft
        !real(rkind), allocatable, dimension(:)     :: gptsC_ifacx, gptsC_jfacy, gptsC_kfacz
        !real(rkind), allocatable, dimension(:)     :: gptsE_ifacx, gptsE_jfacy, gptsE_kfacz
        integer,     allocatable, dimension(:)     ::  imptsC_numonproc, imptsE_numonproc
        integer,     allocatable, dimension(:,:,:) ::  imptsC_indices,   imptsE_indices
        real(rkind), allocatable, dimension(:,:)   ::  imptsC_multfac,   imptsE_multfac
        real(rkind), allocatable, dimension(:)     ::  imptsC_u, imptsC_v, imptsC_w, imptsE_u, imptsE_v, imptsE_w
        real(rkind), allocatable, dimension(:)     ::   gptsC_u,  gptsC_v,  gptsC_w,  gptsE_u,  gptsE_v,  gptsE_w
        real(rkind), allocatable, dimension(:)     ::  gptsC_dst, gptsE_dst

        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure          :: destroy
            procedure, private :: compute_levelset
            procedure, private :: mark_ghost_points
            procedure, private :: save_ghost_points
            procedure, private :: setup_interpolation
            procedure, private :: mark_boundary_points
            procedure, private :: compute_image_points
            procedure, private :: interp_imptsC
            procedure, private :: interp_imptsE
            procedure, private :: update_ghostptsCE
            procedure, private :: smooth_solidptsCE
            procedure, private :: set_interpfac_imptsC
            procedure, private :: set_interpfac_imptsE
    end type 

contains

subroutine init(this, inputDir, inputFile, gpC, gpE, mesh, Lx, Ly, zBot, zTop, dz)
  class(ibmgp),     intent(inout) :: this
  character(len=*), intent(in)    :: inputFile, inputDir
  class(decomp_info), intent(in), target :: gpC, gpE
  real(rkind), dimension(:,:,:,:), intent(in) :: mesh
  real(rkind), intent(in) :: Lx, Ly, zBot, zTop, dz

  character(len=clen)    :: surfaceMeshFile, fname, dumstr
  integer :: io, ii, jj, k, ioUnit, nlines, nlayers = 2, ibwm = 1
  real(rkind) :: rnum1, rnum2, rnum3, Lscale = one, invLscale, dotprod, dx, dy
  real(rkind) :: translate_x = zero, translate_y = zero, translate_z = zero
  real(rkind) :: xmax_ib, xmin_ib, ymax_ib, ymin_ib, zmax_ib, zmin_ib, xkmag
  real(rkind) :: solidpt_x = zero, solidpt_y = zero, solidpt_z = zero
  real(rkind), allocatable, dimension(:) :: xlinepart, ylinepart, zlinepart, zlinepartE
  integer,     allocatable, dimension(:,:,:) :: mapC, mapE, flagC, flagE
  real(rkind), allocatable, dimension(:,:,:) :: levelsetC, levelsetE
  real(rkind), dimension(3) :: vec1, vec2, xk
  real(rkind) :: ibwm_ustar = one

  namelist /IBMGP/ surfaceMeshFile, Lscale, translate_x, translate_y, translate_z, &
                   solidpt_x, solidpt_y, solidpt_z, nlayers, ibwm, ibwm_ustar

  ioUnit = 11
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
  read(unit=ioUnit, NML=IBMGP)
  close(ioUnit)

  this%gpC => gpC
  this%gpE => gpE
  this%ibwm = ibwm
  this%ibwm_ustar = ibwm_ustar

  ! -------Read surfaceMeshFile (x, y, z) for triangles
  ! -- First count the number of elements ---
  fname = inputDir(:len_trim(inputDir))//"/"//trim(surfaceMeshFile)
  ioUnit = 11
  open(unit=ioUnit, file=fname, form='FORMATTED')
  nlines = 0
  do 
      read(ioUnit, *, iostat=io)
      if(io .ne. 0) exit
      nlines = nlines+1
  enddo
  nlines = nlines-2
  this%num_surfelem = nlines/7
  close(ioUnit)

  allocate(this%surfelem(3, this%num_surfelem, 3))
  allocate(this%surfcent(   this%num_surfelem, 3))
  allocate(this%surfnormal( this%num_surfelem, 3))

  ! -- Now read in the vertex data ---
  ioUnit = 11
  open(unit=ioUnit, file=fname, form='FORMATTED')
  read(ioUnit, *)
  do jj = 1, this%num_surfelem 
      ! ignore two lines
      read(ioUnit, *);       read(ioUnit, *)

      ! read first vertex of element ii
      ii = 1
      read(ioUnit, *) dumstr, rnum1, rnum2, rnum3
      this%surfelem(ii,jj,1) = rnum1; this%surfelem(ii,jj,2) = rnum2; this%surfelem(ii,jj,3) = rnum3

      ! read second vertex of element ii
      ii = 2
      read(ioUnit, *) dumstr, rnum1, rnum2, rnum3
      this%surfelem(ii,jj,1) = rnum1; this%surfelem(ii,jj,2) = rnum2; this%surfelem(ii,jj,3) = rnum3

      ! read third vertex of element ii
      ii = 3
      read(ioUnit, *) dumstr, rnum1, rnum2, rnum3
      this%surfelem(ii,jj,1) = rnum1; this%surfelem(ii,jj,2) = rnum2; this%surfelem(ii,jj,3) = rnum3

      ! ignore two lines
      read(ioUnit, *);       read(ioUnit, *)
  enddo
  close(ioUnit)

  if(nrank==0) then
    print *, 'nlines = ', nlines
    print *, 'nelems = ', this%num_surfelem
    do jj = 1, this%num_surfelem 
        print *, '----Element ', jj, '-------------------------------------------'
        print '(a,3(1x,e19.12))', 'point 1:=', this%surfelem(1,jj,1:3)
        print '(a,3(1x,e19.12))', 'point 2:=', this%surfelem(2,jj,1:3)
        print '(a,3(1x,e19.12))', 'point 3:=', this%surfelem(3,jj,1:3)
        print *, '------------------------------------------------------'
    enddo
  endif

  ! IB may be generated at a scale and position that is inconsistent with the
  ! problem setup. Modify the IB surface elements to be consistent with problem
  ! --- First  :: scale the IB surface mesh
  invLscale = one/Lscale
  this%surfelem = invLscale*this%surfelem

  ! --- Second :: translate the IB surface mesh
  this%surfelem(:,:,1) = this%surfelem(:,:,1) - translate_x
  this%surfelem(:,:,2) = this%surfelem(:,:,2) - translate_y
  this%surfelem(:,:,3) = this%surfelem(:,:,3) - translate_z

  ! --- Now check if they are consistent
  xmax_ib = p_maxval(this%surfelem(:,:,1)); xmin_ib = p_minval(this%surfelem(:,:,1))
  ymax_ib = p_maxval(this%surfelem(:,:,2)); ymin_ib = p_minval(this%surfelem(:,:,2))
  zmax_ib = p_maxval(this%surfelem(:,:,3)); zmin_ib = p_minval(this%surfelem(:,:,3))
  if((xmax_ib > Lx) .or. (xmin_ib < zero)) then
      call GracefulExit("x dimension of IB is larger than domain", 111)
  endif
  if((ymax_ib > Ly) .or. (ymin_ib < zero)) then
      call GracefulExit("y dimension of IB is larger than domain", 111)
  endif
  if((zmax_ib > zTop) .or. (zmin_ib < zBot)) then
      call GracefulExit("z dimension of IB is larger than domain", 111)
  endif

  ! compute surface mesh element centroids
  do jj = 1, this%num_surfelem
    this%surfcent(jj,1) = third*sum(this%surfelem(:,jj,1))
    this%surfcent(jj,2) = third*sum(this%surfelem(:,jj,2))
    this%surfcent(jj,3) = third*sum(this%surfelem(:,jj,3))
  enddo

  ! compute surface mesh element normals
  do jj = 1, this%num_surfelem
      vec1(1) = this%surfelem(2,jj,1) - this%surfelem(1,jj,1)
      vec1(2) = this%surfelem(2,jj,2) - this%surfelem(1,jj,2)
      vec1(3) = this%surfelem(2,jj,3) - this%surfelem(1,jj,3)

      vec2(1) = this%surfelem(3,jj,1) - this%surfelem(1,jj,1)
      vec2(2) = this%surfelem(3,jj,2) - this%surfelem(1,jj,2)
      vec2(3) = this%surfelem(3,jj,3) - this%surfelem(1,jj,3)

      ! cross product
      xk(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
      xk(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
      xk(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

      xkmag = sqrt(sum(xk*xk))
      this%surfnormal(jj,1) = xk(1)/(xkmag + 1.0d-18)
      this%surfnormal(jj,2) = xk(2)/(xkmag + 1.0d-18)
      this%surfnormal(jj,3) = xk(3)/(xkmag + 1.0d-18)

      ! ensure it points into the fluid
      vec1(1) = solidpt_x-this%surfcent(jj,1)
      vec1(2) = solidpt_y-this%surfcent(jj,2)
      vec1(3) = solidpt_z-this%surfcent(jj,3)
      
      vec2(:) = this%surfnormal(jj,:)

      dotprod = sum(vec1*vec2)
      if(dotprod > zero) then
          ! normal is probably pointing into the solid so flip the sign
          this%surfnormal(jj,:) = -this%surfnormal(jj,:)
      endif

  enddo

  ! create local grid lines
  allocate(xlinepart(this%gpC%xsz(1)), ylinepart(this%gpC%xsz(2)) )
  allocate(zlinepart(this%gpC%xsz(3)), zlinepartE(this%gpE%xsz(3)))
  xlinepart = mesh(:,1,1,1)
  ylinepart = mesh(1,:,1,2)
  zlinepart = mesh(1,1,:,3)
  do k = 1,this%gpC%xsz(3)
    zlinepartE(k) = zlinepart(k) - half*dz
  enddo
  if(this%gpE%xsz(3) == this%gpC%xsz(3)+1) then
      k = this%gpE%xsz(3)
      zlinepartE(k) = zlinepart(k-1) + half*dz
  endif

  allocate(levelsetC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(levelsetE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
  allocate(mapC     (this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(mapE     (this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
  allocate(flagC    (this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(flagE    (this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))

  call this%compute_levelset(xlinepart, ylinepart, zlinepart, zlinepartE, mapC, mapE, levelsetC, levelsetE)


  call this%mark_ghost_points(gpC, nlayers, mapC, flagC);  this%num_gptsC = sum(flagC)
  call this%mark_ghost_points(gpE, nlayers, mapE, flagE);  this%num_gptsE = sum(flagE)

  allocate(this%gptsC_xyz(this%num_gptsC,3), this%gptsE_xyz(this%num_gptsE,3))
  allocate(this%gptsC_ind(this%num_gptsC,3), this%gptsE_ind(this%num_gptsE,3))
  allocate(this%gptsC_bpt(this%num_gptsC,3), this%gptsE_bpt(this%num_gptsE,3))
  allocate(this%gptsC_bpind(this%num_gptsC), this%gptsE_bpind(this%num_gptsE))

  allocate(this%imptsC_numonproc(this%num_gptsC), this%imptsE_numonproc(this%num_gptsE))
  allocate(this%imptsC_indices(3,8,this%num_gptsC), this%imptsE_indices(3,8,this%num_gptsE))
  allocate(this%imptsC_multfac(  8,this%num_gptsC), this%imptsE_multfac(  8,this%num_gptsE))
  allocate(this%imptsC_u(this%num_gptsC), this%imptsE_u(this%num_gptsE))
  allocate(this%imptsC_v(this%num_gptsC), this%imptsE_v(this%num_gptsE))
  allocate(this%imptsC_w(this%num_gptsC), this%imptsE_w(this%num_gptsE))
  allocate(this%gptsC_u(this%num_gptsC),  this%gptsE_u(this%num_gptsE))
  allocate(this%gptsC_v(this%num_gptsC),  this%gptsE_v(this%num_gptsE))
  allocate(this%gptsC_w(this%num_gptsC),  this%gptsE_w(this%num_gptsE))
  !allocate(this%gptsC_ileft(this%num_gptsC), this%gptsE_ileft(this%num_gptsE))
  !allocate(this%gptsC_jleft(this%num_gptsC), this%gptsE_jleft(this%num_gptsE))
  !allocate(this%gptsC_kleft(this%num_gptsC), this%gptsE_kleft(this%num_gptsE))
  !allocate(this%gptsC_ifacx(this%num_gptsC), this%gptsE_ifacx(this%num_gptsE))
  !allocate(this%gptsC_jfacy(this%num_gptsC), this%gptsE_jfacy(this%num_gptsE))
  !allocate(this%gptsC_kfacz(this%num_gptsC), this%gptsE_kfacz(this%num_gptsE))

  call this%save_ghost_points(flagC, flagE, xlinepart, ylinepart, zlinepart, zlinepartE)

  call this%mark_boundary_points()

  call this%compute_image_points()

  dx = xlinepart(2)-xlinepart(1); dy = ylinepart(2)-ylinepart(1)
  call this%setup_interpolation(dx,dy,dz, zBot, xlinepart, ylinepart, zlinepart)

  deallocate(flagE, flagC, mapC, mapE, levelsetC, levelsetE, xlinepart, ylinepart, zlinepart, zlinepartE)

end subroutine

subroutine update_ibmgp(this, u, v, w, uE, vE, wC)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: u, v, w, uE, vE, wC

  integer :: ii, itmp

  call this%interp_imptsC(u, v,  wC)
  call this%interp_imptsE(uE, vE, w)

  call this%update_ghostptsCE()

  call this%smooth_solidptsCE()

end subroutine

subroutine smooth_solidptsCE(this)
  class(ibmgp),     intent(inout) :: this

  integer :: it, i, j, k, num_smooth, nx, ny, nz
  real(rkind) :: diffcoedd

  nx = this%gpC%xsz(1);  ny = this%gpC%ysz(2);  nz = this%gpC%zsz(3)
  num_smooth = 10; diffcoeff = 1.0d0;

  ! In x decomp; Smooth along x
  diffcoeff_x = diffcoeff ! * (dt/dx^2)
  do it = 1, num_smooth
      ! for u
      do k = 1, this%gpC%xsz(3)
        do j = 1, this%gpC%xsz(2)
          rhs(1,j,k) = zero
          rhs(2:nx-1,j,k) = (u(3:nx,j,k) + u(1:nx-2) - two*u(2:nx-1,j,k))
          rhs(nx,j,k) = zero
        enddo
      enddo
      u = u + diffcoeff_x*rhs*this%solid_maskC_x
      
      ! for v
      do k = 1, this%gpC%xsz(3)
        do j = 1, this%gpC%xsz(2)
          rhs(1,j,k) = zero
          rhs(2:nx-1,j,k) = (v(3:nx,j,k) + v(1:nx-2) - two*v(2:nx-1,j,k))
          rhs(nx,j,k) = zero
        enddo
      enddo
      v = v + diffcoeff_x*rhs*this%solid_maskC_x
      
      ! for v
      do k = 1, this%gpE%xsz(3)
        do j = 1, this%gpE%xsz(2)
          rhsE(1,j,k) = zero
          rhsE(2:nx-1,j,k) = (w(3:nx,j,k) + w(1:nx-2) - two*w(2:nx-1,j,k))
          rhsE(nx,j,k) = zero
        enddo
      enddo
      w = w + diffcoeff_x*rhsE*this%solid_maskE_x
  enddo

end subroutine

subroutine update_ghostptsCE(this)
  class(ibmgp),     intent(inout) :: this

  integer :: ii, jj
  real(rkind) :: vec1(3), vec2(3), unrm(3), utan(3), dotpr, twodist
  real(rkind) :: ibwallstress, utan_gp(3), unrm_gp(3)


  if(this%ibwm==1) then
      ! no slip immersed boundary
      this%gptsC_u = -this%imptsC_u;   this%gptsC_v = -this%imptsC_v;   this%gptsC_w = -this%imptsC_w
      this%gptsE_u = -this%imptsE_u;   this%gptsE_v = -this%imptsE_v;   this%gptsE_w = -this%imptsE_w
  elseif(this%ibwm==2) then
      ! simple log-law immersed boundary
      do ii = 1, this%num_gptsC
          vec1(1) = this%imptsC_u(ii);  vec1(2) = this%imptsC_v(ii);  vec1(3) = this%imptsC_w(ii)
          jj = this%gptsC_bpind(ii)
          vec2(:) = this%surfnormal(jj,:)
          dotpr = sum(vec1*vec2)
          unrm = vec2*dotpr
          utan = vec1-unrm

          twodist = this%gptsC_dst(ii) ! distance from ghost to image point
          ibwallstress = this%ibwm_ustar/kappa/twodist*two
          utan_gp = utan-twodist*ibwallstress
          unrm_gp = -unrm

          this%gptsC_u(ii) = utan_gp(1) + unrm_gp(1)
          this%gptsC_v(ii) = utan_gp(2) + unrm_gp(2)
          this%gptsC_w(ii) = utan_gp(3) + unrm_gp(3)
      enddo

      do ii = 1, this%num_gptsE
          vec1(1) = this%imptsE_u(ii);  vec1(2) = this%imptsE_v(ii);  vec1(3) = this%imptsE_w(ii)
          jj = this%gptsE_bpind(ii)
          vec2(:) = this%surfnormal(jj,:)
          dotpr = sum(vec1*vec2)
          unrm = vec2*dotpr
          utan = vec1-unrm

          twodist = this%gptsE_dst(ii) ! distance from ghost to and image point
          ibwallstress = this%ibwm_ustar/kappa/twodist*two
          utan_gp = utan-twodist*ibwallstress
          unrm_gp = -unrm

          this%gptsE_u(ii) = utan_gp(1) + unrm_gp(1)
          this%gptsE_v(ii) = utan_gp(2) + unrm_gp(2)
          this%gptsE_w(ii) = utan_gp(3) + unrm_gp(3)
      enddo

  endif

end subroutine

subroutine interp_imptsE(this, uE, vE, w)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: uE, vE, w

  integer :: ii, jj, i, j, k
  real(rkind) :: fac, tmpval1, tmpval2, tmpval3

  this%imptsE_u = zero;  this%imptsE_v = zero;  this%imptsE_w = zero
  do ii = 1, this%num_gptsE
      ! how many of the 8 interpolation points are on this proc?
      tmpval1 = zero; tmpval2 = zero; tmpval3 = zero
      do jj = 1, this%imptsE_numonproc(ii)
          i   = this%imptsE_indices(1,jj,ii)
          j   = this%imptsE_indices(2,jj,ii)
          k   = this%imptsE_indices(3,jj,ii)
          fac = this%imptsE_multfac(jj, ii)
          tmpval1 = tmpval1 + fac*uE(i,j,k)
          tmpval2 = tmpval2 + fac*vE(i,j,k)
          tmpval3 = tmpval3 + fac* w(i,j,k)
      enddo
      this%imptsE_u(ii) = tmpval1
      this%imptsE_v(ii) = tmpval2
      this%imptsE_w(ii) = tmpval3
  enddo

end subroutine

subroutine interp_imptsC(this, u, v, wC)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: u, v, wC

  integer :: ii, jj, i, j, k
  real(rkind) :: fac, tmpval1, tmpval2, tmpval3

  this%imptsC_u = zero;  this%imptsC_v = zero;  this%imptsC_w = zero
  do ii = 1, this%num_gptsC
      ! how many of the 8 interpolation points are on this proc?
      tmpval1 = zero; tmpval2 = zero; tmpval3 = zero
      do jj = 1, this%imptsC_numonproc(ii)
          i   = this%imptsC_indices(1,jj,ii)
          j   = this%imptsC_indices(2,jj,ii)
          k   = this%imptsC_indices(3,jj,ii)
          fac = this%imptsC_multfac(jj, ii)
          tmpval1 = tmpval1 + fac* u(i,j,k)
          tmpval2 = tmpval2 + fac* v(i,j,k)
          tmpval3 = tmpval2 + fac*wC(i,j,k)
      enddo
      this%imptsC_u(ii) = tmpval1
      this%imptsC_v(ii) = tmpval2
      this%imptsC_w(ii) = tmpval3
  enddo

end subroutine

subroutine setup_interpolation(this, dx, dy, dz, zBot, xlinepart, ylinepart, zlinepart)
  class(ibmgp),              intent(inout) :: this
  real(rkind),               intent(in)    :: dx, dy, dz, zBot
  real(rkind), dimension(:), intent(in)    :: xlinepart, ylinepart, zlinepart

  integer :: ii, itmp, jtmp, ktmp, nxloc, nyloc, nzloc, point_number
  real(rkind) :: xdom_left, xdom_right, ydom_left, ydom_right, zdom_left, zdom_right
  real(rkind) :: xloc, yloc, zloc, facx, facy, facz, onemfacx, onemfacy, onemfacz

   ! set imptsC_numonproc (any number between 1 and 8)
   ! fill imptsC_indices(1:3,1:this%imptsc_numonproc)
   ! fill imptsC_multfac(1:this%imptsc_numonproc)    :: depending on whether point (i,j,k) is point 1, 2, 3, .. or 8, some combination of facx, facy, facz, (1-facx), (1-facy) and (1-facz)

  nxloc = this%gpC%xsz(1);  nyloc = this%gpC%xsz(2);  nzloc = this%gpC%xsz(3)

  ! note :: padded domain extents
  xdom_left = xlinepart(1)-dx;   xdom_right = xlinepart(nxloc)+dx
  ydom_left = ylinepart(1)-dy;   ydom_right = ylinepart(nyloc)+dy
  zdom_left = zlinepart(1)-dz;   zdom_right = zlinepart(nzloc)+dz

  do ii = 1, this%num_gptsC
      ! first check if this point is on this processor
      xloc = this%gptsC_img(ii,1);  yloc = this%gptsC_img(ii,2);  zloc = this%gptsC_img(ii,3);

      if( (xloc>=xdom_left) .and. (xloc<=xdom_right) .and. &
          (yloc>=ydom_left) .and. (yloc<=ydom_right) .and. &
          (zloc>=zdom_left) .and. (zloc<=zdom_right) ) then
      endif
      itmp = floor((xloc-xdom_left)/dx);        
      jtmp = floor((yloc-ydom_left)/dy);        
      ktmp = floor((zloc-zdom_left)/dz);        

      ! determine facx, facy, facz
      if(itmp==0) then
          facx = one - (xlinepart(1)-xloc)/dx
      else
          facx = one - (xloc-xlinepart(itmp))/dx
      endif

      if(jtmp==0) then
          facy = one - (ylinepart(1)-yloc)/dy
      else
          facy = one - (yloc-ylinepart(jtmp))/dy
      endif

      if(ktmp==0) then
          facz = one - (zlinepart(1)-zloc)/dz
      else
          facz = one - (zloc-zlinepart(ktmp))/dz
      endif

      onemfacx = one-facx;   onemfacy = one-facy;   onemfacz = one-facz

      ! now consider 8 points separately
      point_number = 0
      call this%set_interpfac_imptsC(itmp  , jtmp  , ktmp  , ii, point_number,     facx*    facy*    facz)
      call this%set_interpfac_imptsC(itmp+1, jtmp  , ktmp  , ii, point_number, onemfacx*    facy*    facz)
      call this%set_interpfac_imptsC(itmp  , jtmp+1, ktmp  , ii, point_number,     facx*onemfacy*    facz)
      call this%set_interpfac_imptsC(itmp+1, jtmp+1, ktmp  , ii, point_number, onemfacx*onemfacy*    facz)
      call this%set_interpfac_imptsC(itmp  , jtmp  , ktmp+1, ii, point_number,     facx*    facy*onemfacz)
      call this%set_interpfac_imptsC(itmp+1, jtmp  , ktmp+1, ii, point_number, onemfacx*    facy*onemfacz)
      call this%set_interpfac_imptsC(itmp  , jtmp+1, ktmp+1, ii, point_number,     facx*onemfacy*onemfacz)
      call this%set_interpfac_imptsC(itmp+1, jtmp+1, ktmp+1, ii, point_number, onemfacx*onemfacy*onemfacz)

      !this%gptsC_ileft(ii) = max(1, min(itmp+1, this%gpC%xsz(1)-1))
      !this%gptsC_ifacx(ii) = (this%gptsC_img(ii,1)-real(itmp, rkind)*dx)/dx

      !itmp = floor(this%gptsC_img(ii,2)/dy);
      !this%gptsC_jleft(ii) = max(1, min(itmp+1, this%gpC%ysz(2)-1))
      !this%gptsC_jfacy(ii) = (this%gptsC_img(ii,2)-real(itmp, rkind)*dy)/dy

      !itmp = floor((this%gptsC_img(ii,3)-zBot)/dz);
      !this%gptsC_kleft(ii) = max(1, min(itmp+1, this%gpC%zsz(3)-1))
      !this%gptsC_kfacz(ii) = (this%gptsC_img(ii,3)-real(itmp, rkind)*dz - zBot)/dz
  enddo

  do ii = 1, this%num_gptsE
      ! first check if this point is on this processor
      xloc = this%gptsE_img(ii,1);  yloc = this%gptsE_img(ii,2);  zloc = this%gptsE_img(ii,3);

      if( (xloc>=xdom_left) .and. (xloc<=xdom_right) .and. &
          (yloc>=ydom_left) .and. (yloc<=ydom_right) .and. &
          (zloc>=zdom_left) .and. (zloc<=zdom_right) ) then
      endif
      itmp = floor((xloc-xdom_left)/dx);        
      jtmp = floor((yloc-ydom_left)/dy);        
      ktmp = floor((zloc-zdom_left)/dz);        

      ! determine facx, facy, facz
      if(itmp==0) then
          facx = one - (xlinepart(1)-xloc)/dx
      else
          facx = one - (xloc-xlinepart(itmp))/dx
      endif

      if(jtmp==0) then
          facy = one - (ylinepart(1)-yloc)/dy
      else
          facy = one - (yloc-ylinepart(jtmp))/dy
      endif

      if(ktmp==0) then
          facz = one - (zlinepart(1)-zloc)/dz
      else
          facz = one - (zloc-zlinepart(ktmp))/dz
      endif

      onemfacx = one-facx;   onemfacy = one-facy;   onemfacz = one-facz

      ! now consider 8 points separately
      point_number = 0
      call this%set_interpfac_imptsE(itmp  , jtmp  , ktmp  , ii, point_number,     facx*    facy*    facz)
      call this%set_interpfac_imptsE(itmp+1, jtmp  , ktmp  , ii, point_number, onemfacx*    facy*    facz)
      call this%set_interpfac_imptsE(itmp  , jtmp+1, ktmp  , ii, point_number,     facx*onemfacy*    facz)
      call this%set_interpfac_imptsE(itmp+1, jtmp+1, ktmp  , ii, point_number, onemfacx*onemfacy*    facz)
      call this%set_interpfac_imptsE(itmp  , jtmp  , ktmp+1, ii, point_number,     facx*    facy*onemfacz)
      call this%set_interpfac_imptsE(itmp+1, jtmp  , ktmp+1, ii, point_number, onemfacx*    facy*onemfacz)
      call this%set_interpfac_imptsE(itmp  , jtmp+1, ktmp+1, ii, point_number,     facx*onemfacy*onemfacz)
      call this%set_interpfac_imptsE(itmp+1, jtmp+1, ktmp+1, ii, point_number, onemfacx*onemfacy*onemfacz)

      !itmp = floor(this%gptsE_img(ii,1)/dx);
      !this%gptsE_ileft(ii) = max(1, min(itmp+1, this%gpE%xsz(1)-1))
      !this%gptsE_ifacx(ii) = (this%gptsE_img(ii,1)-real(itmp, rkind)*dx)/dx

      !itmp = floor(this%gptsE_img(ii,2)/dy);
      !this%gptsE_jleft(ii) = max(1, min(itmp+1, this%gpE%ysz(2)-1))
      !this%gptsE_jfacy(ii) = (this%gptsE_img(ii,2)-real(itmp, rkind)*dy)/dy

      !itmp = floor((this%gptsE_img(ii,3)-zBot)/dz);
      !this%gptsE_kleft(ii) = max(1, min(itmp+1, this%gpE%zsz(3)-1))
      !this%gptsE_kfacz(ii) = (this%gptsE_img(ii,3)-real(itmp, rkind)*dz - zBot)/dz
  enddo

end subroutine

subroutine set_interpfac_imptsC(this, ip, jp, kp, ii, point_number, multfac)
  class(ibmgp),     intent(inout) :: this
  integer,     intent(in)     :: ip, jp, kp, ii
  integer,     intent(inout)  :: point_number
  real(rkind), intent(in)     :: multfac

  integer :: nxloc, nyloc, nzloc

  nxloc = this%gpC%xsz(1);  nyloc = this%gpC%xsz(2);  nzloc = this%gpC%xsz(3)

  if( (ip .ge. 1) .and. (ip .le. nxloc) .and. &
      (jp .ge. 1) .and. (jp .le. nyloc) .and. &
      (kp .ge. 1) .and. (kp .le. nzloc) ) then
          point_number = point_number + 1
          this%imptsC_numonproc(ii) = point_number
          this%imptsC_indices(1,point_number,ii) = ip
          this%imptsC_indices(2,point_number,ii) = jp
          this%imptsC_indices(3,point_number,ii) = kp
          this%imptsC_multfac(  point_number,ii) = multfac
  endif

end subroutine

subroutine set_interpfac_imptsE(this, ip, jp, kp, ii, point_number, multfac)
  class(ibmgp),     intent(inout) :: this
  integer,     intent(in)    :: ip, jp, kp, ii
  integer,     intent(inout) :: point_number
  real(rkind), intent(in)    :: multfac

  integer :: nxloc, nyloc, nzloc

  nxloc = this%gpC%xsz(1);  nyloc = this%gpC%xsz(2);  nzloc = this%gpC%xsz(3)

  if( (ip .ge. 1) .and. (ip .le. nxloc) .and. &
      (jp .ge. 1) .and. (jp .le. nyloc) .and. &
      (kp .ge. 1) .and. (kp .le. nzloc) ) then
          point_number = point_number + 1
          this%imptsE_numonproc(ii) = point_number
          this%imptsE_indices(1,point_number,ii) = ip
          this%imptsE_indices(2,point_number,ii) = jp
          this%imptsE_indices(3,point_number,ii) = kp
          this%imptsE_multfac(  point_number,ii) = multfac
  endif

end subroutine

subroutine compute_image_points(this)
  class(ibmgp),     intent(inout) :: this

  integer :: ii, jj
  real(rkind) :: vec1(3), vec2(3), dotpr

  do ii = 1, this%num_gptsC
      ! compute image of gptsC_xyz(ii,:) wrt gptsc_bpt(ii,:) and
      vec1(:) = this%gptsC_bpt(ii,:) - this%gptsC_xyz(ii,:)
      jj = this%gptsC_bpind(ii)
      vec2(:) = this%surfnormal(jj,:)
      dotpr = sum(vec1*vec2)
      this%gptsC_img(ii,:) = this%gptsC_xyz(ii,:) + two*dotpr*vec2
      this%gptsC_bnp(ii,:) = this%gptsC_xyz(ii,:) + dotpr*vec2
      this%gptsC_dst(ii)   = two*dotpr
  enddo

  do ii = 1, this%num_gptsE
      ! compute image of gptsC_xyz(ii,:) wrt gptsc_bpt(ii,:) and
      vec1(:) = this%gptsE_bpt(ii,:) - this%gptsE_xyz(ii,:)
      jj = this%gptsE_bpind(ii)
      vec2(:) = this%surfnormal(jj,:)
      dotpr = sum(vec1*vec2)
      this%gptsE_img(ii,:) = this%gptsE_xyz(ii,:) + two*dotpr*vec2
      this%gptsE_bnp(ii,:) = this%gptsE_xyz(ii,:) + dotpr*vec2
      this%gptsE_dst(ii)   = two*dotpr
  enddo

end subroutine

subroutine mark_boundary_points(this)
  class(ibmgp),     intent(inout) :: this

  integer :: ii, imin
  real(rkind), allocatable, dimension(:) :: distfn

  allocate(distfn(this%num_surfelem))

  do ii = 1, this%num_gptsC
      ! find the distance between this ghost point and all surface centroids
      distfn(:) = (this%gptsC_xyz(ii,1) - this%surfcent(:,1))**2 + &
                  (this%gptsC_xyz(ii,2) - this%surfcent(:,2))**2 + &
                  (this%gptsC_xyz(ii,3) - this%surfcent(:,3))**2 
      imin = minloc(distfn,1)
      this%gptsC_bpt(ii,1) = this%surfcent(imin,1)
      this%gptsC_bpt(ii,2) = this%surfcent(imin,2)
      this%gptsC_bpt(ii,3) = this%surfcent(imin,3)
      this%gptsC_bpind(ii) = imin
  enddo

  do ii = 1, this%num_gptsE
      ! find the distance between this ghost point and all surface centroids
      distfn(:) = (this%gptsE_xyz(ii,1) - this%surfcent(:,1))**2 + &
                  (this%gptsE_xyz(ii,2) - this%surfcent(:,2))**2 + &
                  (this%gptsE_xyz(ii,3) - this%surfcent(:,3))**2 
      imin = minloc(distfn,1)
      this%gptsE_bpt(ii,1) = this%surfcent(imin,1)
      this%gptsE_bpt(ii,2) = this%surfcent(imin,2)
      this%gptsE_bpt(ii,3) = this%surfcent(imin,3)
      this%gptsE_bpind(ii) = imin
  enddo

  deallocate(distfn)

end subroutine

subroutine save_ghost_points(this, flagC, flagE, xlinepart, ylinepart, zlinepart, zlinepartE)
  class(ibmgp),     intent(inout) :: this
  integer, dimension(:,:,:), intent(in)   :: flagC, flagE
  real(rkind), dimension(:), intent(in) :: xlinepart, ylinepart, zlinepart, zlinepartE

  integer :: i, j, k, ii

  ii = 0
  do k = 1, this%gpC%xsz(3)
    do j = 1, this%gpC%xsz(2)
      do i = 1, this%gpC%xsz(1)
          if(flagC(i,j,k)==1) then
             ii = ii+1
             
             this%gptsC_xyz(ii,1) = xlinepart(i)
             this%gptsC_xyz(ii,2) = ylinepart(j)
             this%gptsC_xyz(ii,3) = zlinepart(k)

             this%gptsC_ind(ii,1) = i
             this%gptsC_ind(ii,2) = j
             this%gptsC_ind(ii,3) = k
          endif
      enddo
    enddo
  enddo

  ii = 0
  do k = 1, this%gpE%xsz(3)
    do j = 1, this%gpE%xsz(2)
      do i = 1, this%gpE%xsz(1)
          if(flagE(i,j,k)==1) then
             ii = ii+1
             this%gptsE_xyz(ii,1) = xlinepart(i)
             this%gptsE_xyz(ii,2) = ylinepart(j)
             this%gptsE_xyz(ii,3) = zlinepartE(k)

             this%gptsE_ind(ii,1) = i
             this%gptsE_ind(ii,2) = j
             this%gptsE_ind(ii,3) = k
             
          endif
      enddo
    enddo
  enddo



end subroutine

subroutine mark_ghost_points(this, gp, nlayers, map, flag)
  class(ibmgp),     intent(in) :: this
  class(decomp_info), pointer, intent(in) :: gp
  integer, intent(in) :: nlayers
  integer, dimension(:,:,:), intent(in) :: map
  integer, dimension(:,:,:), intent(out) :: flag

  integer :: i, j, k, ii, jj
  integer, allocatable, dimension(:) :: locflags

  flag = 0
  allocate(locflags(3*nlayers))
  ! mark ghost points
  do k = 1, gp%xsz(3)
    do j = 1, gp%xsz(2)
      do i = 1, gp%xsz(1)
          if(map(i,j,k)==0) then
             ! point (i,j,k) is solid; check its neighbours
             ii = 0; locflags = 0
             do jj = 1, nlayers
                 ii = ii+1;   locflags(ii) = map(i+jj, j,    k   )
                 ii = ii+1;   locflags(ii) = map(i,    j+jj, k   )
                 ii = ii+1;   locflags(ii) = map(i,    j,    k+jj)
             enddo 
             if(sum(locflags) > 0) flag(i,j,k) = 1
          endif
      enddo
    enddo
  enddo
  deallocate(locflags)

end subroutine

subroutine compute_levelset(this, xlinepart, ylinepart, zlinepart, zlinepartE, mapC, mapE, levelsetC, levelsetE)
  class(ibmgp),     intent(in) :: this
  real(rkind), dimension(:), intent(in) :: xlinepart, ylinepart, zlinepart, zlinepartE
  real(rkind), dimension(:,:,:), intent(out) :: levelsetC, levelsetE
  integer, dimension(:,:,:), intent(out) :: mapC, mapE

  integer :: i, j, k, jj
  real(rkind) :: xk(3), xkmag
  integer, allocatable, dimension(:,:,:) :: lnearlistC, lnearlistE

  call initialize_kdtree_ib(this%num_surfelem, this%surfcent(:,1), this%surfcent(:,2), this%surfcent(:,3))
  call create_kdtree_ib()

  allocate(lnearlistC(1:this%gpC%xsz(1), 1:this%gpC%xsz(2), 1:this%gpC%xsz(3)))
  allocate(lnearlistE(1:this%gpE%xsz(1), 1:this%gpE%xsz(2), 1:this%gpE%xsz(3)))

  ! compute levelsetC using lnearlistC
  call probe_nearest_points_kdtree_ib(1, this%gpC%xsz(1), 1, this%gpC%xsz(2), 1, this%gpC%xsz(3), &
                                      xlinepart, ylinepart, zlinepart, lnearlistC)
  mapC = 0
  do k = 1, this%gpC%xsz(3)
    do j = 1, this%gpC%xsz(2)
      do i = 1, this%gpC%xsz(1)
          jj = lnearlistC(i,j,k)
          xk(1) = xlinepart(i)-this%surfcent(jj,1)
          xk(2) = ylinepart(i)-this%surfcent(jj,2)
          xk(3) = zlinepart(i)-this%surfcent(jj,3)
          xkmag = sqrt(sum(xk*xk))
          levelsetC(i,j,k) = sign(1.0_rkind, dot_product(xk, this%surfnormal(jj,:)))*xkmag
      enddo
    enddo
  enddo
  where (levelsetC > zero)
    mapC = 1
  endwhere

  ! compute levelsetE using lnearlistE
  call probe_nearest_points_kdtree_ib(1, this%gpE%xsz(1), 1, this%gpE%xsz(2), 1, this%gpE%xsz(3), &
                                      xlinepart, ylinepart, zlinepartE, lnearlistE)
  mapE = 0
  do k = 1, this%gpE%xsz(3)
    do j = 1, this%gpE%xsz(2)
      do i = 1, this%gpE%xsz(1)
          jj = lnearlistE(i,j,k)
          xk(1) = xlinepart(i)-this%surfcent(jj,1)
          xk(2) = ylinepart(i)-this%surfcent(jj,2)
          xk(3) = zlinepartE(i)-this%surfcent(jj,3)
          xkmag = sqrt(sum(xk*xk))
          levelsetE(i,j,k) = sign(1.0_rkind, dot_product(xk, this%surfnormal(jj,:)))*xkmag
      enddo
    enddo
  enddo
  where (levelsetE > zero)
    mapE = 1
  endwhere
  
  call finalize_kdtree_ib()
  deallocate (lnearlistC, lnearlistE)

end subroutine


subroutine destroy(this)
  class(ibmgp), intent(inout) :: this

  !deallocate(this%gptsC_ifacx, this%gptsE_ifacx)
  !deallocate(this%gptsC_jfacy, this%gptsE_jfacy)
  !deallocate(this%gptsC_kfacz, this%gptsE_kfacz)
  !deallocate(this%gptsC_ileft, this%gptsE_ileft)
  !deallocate(this%gptsC_jleft, this%gptsE_jleft)
  !deallocate(this%gptsC_kleft, this%gptsE_kleft)
  deallocate(this%gptsC_u,  this%gptsE_u)
  deallocate(this%gptsC_v,  this%gptsE_v)
  deallocate(this%gptsC_w,  this%gptsE_w)
  deallocate(this%imptsC_u, this%imptsE_u)
  deallocate(this%imptsC_v, this%imptsE_v)
  deallocate(this%imptsC_w, this%imptsE_w)
  deallocate(this%imptsC_numonproc, this%imptsE_numonproc)
  deallocate(this%imptsC_indices, this%imptsE_indices)
  deallocate(this%imptsC_multfac, this%imptsE_multfac)
  deallocate(this%gptsC_dst, this%gptsE_dst)
  deallocate(this%gptsC_bnp, this%gptsE_bnp)
  deallocate(this%gptsC_img, this%gptsE_img)
  deallocate(this%gptsC_xyz, this%gptsE_xyz)
  deallocate(this%gptsC_ind, this%gptsE_ind)
  deallocate(this%surfnormal)
  deallocate(this%surfcent)
  deallocate(this%surfelem)

end subroutine

end module
