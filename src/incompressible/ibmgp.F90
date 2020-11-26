module ibmgpmod
    use kind_parameters, only: rkind, clen
    use constants      , only: zero, one, third, half
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

        integer :: num_surfelem
        real(rkind), allocatable, dimension(:,:,:) :: surfelem
        real(rkind), allocatable, dimension(:,:)   :: surfcent, surfnormal

        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure          :: destroy
            procedure, private :: compute_levelset
            procedure, private :: mark_ghost_points
    end type 

contains

subroutine init(this, inputDir, inputFile, gpC, gpE, mesh, Lx, Ly, zBot, zTop, dz)
  class(ibmgp),     intent(inout) :: this
  character(len=*), intent(in)    :: inputFile, inputDir
  class(decomp_info), intent(in), target :: gpC, gpE
  real(rkind), dimension(:,:,:,:), intent(in) :: mesh
  real(rkind), intent(in) :: Lx, Ly, zBot, zTop, dz

  character(len=clen)    :: surfaceMeshFile, fname, dumstr
  integer :: io, ii, jj, k, ioUnit, nlines, nlayers = 2
  real(rkind) :: rnum1, rnum2, rnum3, Lscale = one, invLscale, dotprod
  real(rkind) :: translate_x = zero, translate_y = zero, translate_z = zero
  real(rkind) :: xmax_ib, xmin_ib, ymax_ib, ymin_ib, zmax_ib, zmin_ib, xkmag
  real(rkind) :: solidpt_x = zero, solidpt_y = zero, solidpt_z = zero
  real(rkind), allocatable, dimension(:) :: xlinepart, ylinepart, zlinepart, zlinepartE
  integer,     allocatable, dimension(:,:,:) :: mapC, mapE, flagC, flagE
  real(rkind), allocatable, dimension(:,:,:) :: levelsetC, levelsetE
  real(rkind), dimension(3) :: vec1, vec2, xk

  namelist /IBMGP/ surfaceMeshFile, Lscale, translate_x, translate_y, translate_z, &
                   solidpt_x, solidpt_y, solidpt_z, nlayers

  ioUnit = 11
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
  read(unit=ioUnit, NML=IBMGP)
  close(ioUnit)

  this%gpC => gpC
  this%gpE => gpE

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


  call this%mark_ghost_points(gpC, nlayers, mapC, flagC);  this%num_ghostptsC = sum(flagC)
  call this%mark_ghost_points(gpE, nlayers, mapE, flagE);  this%num_ghostptsE = sum(flagE)



  deallocate(flagE, flagC, mapC, mapE, levelsetC, levelsetE, xlinepart, ylinepart, zlinepart, zlinepartE)

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

  deallocate(this%surfnormal)
  deallocate(this%surfcent)
  deallocate(this%surfelem)

end subroutine

end module
