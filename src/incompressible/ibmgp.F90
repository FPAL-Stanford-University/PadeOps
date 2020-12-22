module ibmgpmod
    use kind_parameters, only: rkind, clen, mpiinteg, mpirkind
    use constants      , only: zero, one, third, half, two, kappa
    use exits          , only: GracefulExit
    use reductions     , only: p_maxval, p_minval
    use kdtree_wrapper , only: initialize_kdtree_ib, create_kdtree_ib, probe_nearest_points_kdtree_ib, finalize_kdtree_ib
    use spectralMod    , only: spectral  
    use decomp_2d
    use decomp_2d_io
    use mpi

    implicit none

    private
    public :: ibmgp

    type :: ibmgp
        private 
        type(decomp_info), pointer :: gpC, gpE
        class(spectral), pointer :: spectC, spectE

        integer     :: num_surfelem, num_gptsC, num_gptsE, ibwm, runID
        integer     :: num_imptsC_glob, num_imptsE_glob
        real(rkind) :: ibwm_ustar, ibwm_z0
        real(rkind), allocatable, dimension(:,:,:) :: surfelem
        real(rkind), allocatable, dimension(:,:)   :: surfcent, surfnormal
        real(rkind), allocatable, dimension(:,:)   :: gptsC_xyz, gptsE_xyz, gptsC_bpt, gptsE_bpt
        real(rkind), allocatable, dimension(:,:)   :: gptsC_img, gptsE_img, gptsC_bnp, gptsE_bnp
        integer,     allocatable, dimension(:,:)   :: gptsC_ind, gptsE_ind
        integer,     allocatable, dimension(:)     :: gptsC_bpind, gptsE_bpind!, gptsC_ileft, gptsC_jleft, gptsC_kleft
        integer,     allocatable, dimension(:)     :: num_imptsC_arr, imptsC_index_st
        integer,     allocatable, dimension(:)     :: num_imptsE_arr, imptsE_index_st
        !integer,     allocatable, dimension(:)     :: gptsE_ileft, gptsE_jleft, gptsE_kleft
        !real(rkind), allocatable, dimension(:)     :: gptsC_ifacx, gptsC_jfacy, gptsC_kfacz
        !real(rkind), allocatable, dimension(:)     :: gptsE_ifacx, gptsE_jfacy, gptsE_kfacz
        !real(rkind), allocatable, dimension(:)     ::   gptsC_u,  gptsC_v,  gptsC_w,  gptsE_u,  gptsE_v,  gptsE_w
        integer,     allocatable, dimension(:)     ::  imptsC_numonproc, imptsE_numonproc
        integer,     allocatable, dimension(:,:,:) ::  imptsC_indices,   imptsE_indices
        real(rkind), allocatable, dimension(:,:)   ::  imptsC_multfac,   imptsE_multfac, imptsC_xyz, imptsE_xyz
        real(rkind), allocatable, dimension(:,:)   ::  imptsC_xyz_loc, imptsE_xyz_loc
        real(rkind), allocatable, dimension(:)     ::  imptsC_u, imptsC_v, imptsC_w, imptsE_u, imptsE_v, imptsE_w
        real(rkind), allocatable, dimension(:)     ::  imptsC_u_tmp, imptsC_v_tmp, imptsC_w_tmp
        real(rkind), allocatable, dimension(:)     ::  imptsE_u_tmp, imptsE_v_tmp, imptsE_w_tmp
        real(rkind), allocatable, dimension(:)     ::  gptsC_dst, gptsE_dst

        real(rkind), dimension(:,:,:), pointer     :: rbuffxC1, rbuffxC2, rbuffyC1, rbuffyC2, rbuffzC1, rbuffzC2 
        real(rkind), dimension(:,:,:), pointer     :: rbuffxE1, rbuffxE2, rbuffyE1, rbuffyE2, rbuffzE1, rbuffzE2 
        real(rkind), allocatable, dimension(:,:,:) ::  mask_solid_xC, mask_solid_yC, mask_solid_zC
        real(rkind), allocatable, dimension(:,:,:) ::  mask_solid_xE, mask_solid_yE, mask_solid_zE

        character(len=clen) :: outputDir

        integer :: debug_flag = 1, debug_elemid = 10

        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure          :: destroy
            procedure          :: update_ibmgp
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
            procedure, private :: smooth_along_x
            procedure, private :: smooth_along_y
            procedure, private :: smooth_along_z
    end type 

contains

subroutine init(this, inputDir, inputFile, outputDir, runID, gpC, gpE, spectC, spectE, mesh, Lx, Ly, zBot, zTop, dz, rbuffxC, rbuffyC, rbuffzC, rbuffxE, rbuffyE, rbuffzE)
  class(ibmgp),     intent(inout) :: this
  character(len=*), intent(in)    :: inputFile, inputDir, outputDir
  integer, intent(in) :: runID
  type(decomp_info), intent(in), target :: gpC, gpE
  type(spectral),    intent(in), target :: spectC, spectE
  real(rkind), dimension(:,:,:,:), intent(in) :: mesh
  real(rkind), dimension(:,:,:,:), intent(in), target :: rbuffxC, rbuffyC, rbuffzC, rbuffxE, rbuffyE, rbuffzE
  real(rkind), intent(in) :: Lx, Ly, zBot, zTop, dz

  character(len=clen)    :: surfaceMeshFile, fname, dumstr
  integer :: io, ii, jj, k, ioUnit, nlines, nlayers = 2, ibwm = 1, ierr
  real(rkind) :: rnum1, rnum2, rnum3, Lscale = one, invLscale, dotprod, dx, dy
  real(rkind) :: translate_x = zero, translate_y = zero, translate_z = zero
  real(rkind) :: xmax_ib, xmin_ib, ymax_ib, ymin_ib, zmax_ib, zmin_ib, xkmag
  real(rkind) :: solidpt_x = zero, solidpt_y = zero, solidpt_z = zero
  real(rkind), allocatable, dimension(:) :: xlinepart, ylinepart, zlinepart, zlinepartE
  integer,     allocatable, dimension(:,:,:) :: mapC, mapE, flagC, flagE
  real(rkind), allocatable, dimension(:,:,:) :: levelsetC, levelsetE
  real(rkind), dimension(3) :: vec1, vec2, xk
  real(rkind) :: ibwm_ustar = one, ibwm_z0 = 1.0d-4

  namelist /IBMGP/ surfaceMeshFile, Lscale, translate_x, translate_y, translate_z, &
                   solidpt_x, solidpt_y, solidpt_z, nlayers, ibwm, ibwm_ustar, ibwm_z0

  this%runID = runID
  this%outputDir = outputDir

  ioUnit = 11
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
  read(unit=ioUnit, NML=IBMGP)
  close(ioUnit)

  this%gpC => gpC
  this%gpE => gpE
  this%spectC => spectC
  this%spectE => spectE
  this%ibwm = ibwm
  this%ibwm_ustar = ibwm_ustar
  this%ibwm_z0    = ibwm_z0

  if( (this%ibwm .ne. 1) .and. (this%ibwm .ne. 2)) then
      call GracefulExit("Wrong choice for IB Wall Model. Check input file.", 111)
  endif

  if( (this%ibwm==2) .and. (this%ibwm_ustar < zero)) then
      call GracefulExit("IB Wall Model ustar cannot be negative. Check input file.", 111)
  endif

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

  !if(nrank==0) then
  !  print *, 'nlines = ', nlines
  !  print *, 'nelems = ', this%num_surfelem
  !  do jj = 1, this%num_surfelem 
  !      print *, '----Element ', jj, '-------------------------------------------'
  !      print '(a,3(1x,e19.12))', 'point 1:=', this%surfelem(1,jj,1:3)
  !      print '(a,3(1x,e19.12))', 'point 2:=', this%surfelem(2,jj,1:3)
  !      print '(a,3(1x,e19.12))', 'point 3:=', this%surfelem(3,jj,1:3)
  !      print *, '------------------------------------------------------'
  !  enddo
  !endif

  ! IB may be generated at a scale and position that is inconsistent with the
  ! problem setup. Modify the IB surface elements to be consistent with problem
  ! --- First  :: scale the IB surface mesh
  invLscale = one/Lscale
  this%surfelem = invLscale*this%surfelem

  ! --- Second :: translate the IB surface mesh
  this%surfelem(:,:,1) = this%surfelem(:,:,1) + translate_x
  this%surfelem(:,:,2) = this%surfelem(:,:,2) + translate_y
  this%surfelem(:,:,3) = this%surfelem(:,:,3) + translate_z

  ! --- Now check if they are consistent
  xmax_ib = p_maxval(this%surfelem(:,:,1)); xmin_ib = p_minval(this%surfelem(:,:,1))
  ymax_ib = p_maxval(this%surfelem(:,:,2)); ymin_ib = p_minval(this%surfelem(:,:,2))
  zmax_ib = p_maxval(this%surfelem(:,:,3)); zmin_ib = p_minval(this%surfelem(:,:,3))
  if((this%debug_flag==1) .and. (nrank==0)) then
        print *, '------------------------------------------------------'
        print '(a,3(1x,e19.12))', '(xmin, xmax):=', xmin_ib, xmax_ib
        print '(a,3(1x,e19.12))', '(ymin, ymax):=', ymin_ib, ymax_ib
        print '(a,3(1x,e19.12))', '(zmin, zmax):=', zmin_ib, zmax_ib
        print *, '------------------------------------------------------'
  endif

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

  if(this%debug_flag==1) then
    if(nrank==0) then
      jj = this%debug_elemid
      print '(a,3(e19.12,1x))', 'pt. 1 : ', this%surfelem(1,jj,:)
      print '(a,3(e19.12,1x))', 'pt. 2 : ', this%surfelem(2,jj,:)
      print '(a,3(e19.12,1x))', 'pt. 3 : ', this%surfelem(3,jj,:)
      print '(a,3(e19.12,1x))', 'centr : ', this%surfcent(jj,:)
      print '(a,3(e19.12,1x))', 'normal: ', this%surfnormal(jj,:)
    endif
  endif

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

  this%rbuffxC1 => rbuffxC(:,:,:,1);  this%rbuffxC2 => rbuffxC(:,:,:,2)
  this%rbuffyC1 => rbuffyC(:,:,:,1);  this%rbuffyC2 => rbuffyC(:,:,:,2)
  this%rbuffzC1 => rbuffzC(:,:,:,1);  this%rbuffzC2 => rbuffzC(:,:,:,2)
  this%rbuffxE1 => rbuffxE(:,:,:,1);  this%rbuffxE2 => rbuffxE(:,:,:,2)
  this%rbuffyE1 => rbuffyE(:,:,:,1);  this%rbuffyE2 => rbuffyE(:,:,:,2)
  this%rbuffzE1 => rbuffzE(:,:,:,1);  this%rbuffzE2 => rbuffzE(:,:,:,2)
  allocate(this%mask_solid_xC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(this%mask_solid_yC(this%gpC%ysz(1), this%gpC%ysz(2), this%gpC%ysz(3)))
  allocate(this%mask_solid_zC(this%gpC%zsz(1), this%gpC%zsz(2), this%gpC%zsz(3)))
  allocate(this%mask_solid_xE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
  allocate(this%mask_solid_yE(this%gpE%ysz(1), this%gpE%ysz(2), this%gpE%ysz(3)))
  allocate(this%mask_solid_zE(this%gpE%zsz(1), this%gpE%zsz(2), this%gpE%zsz(3)))

  allocate(levelsetC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(levelsetE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
  allocate(mapC     (this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(mapE     (this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
  allocate(flagC    (this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
  allocate(flagE    (this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))

  call this%compute_levelset(xlinepart, ylinepart, zlinepart, zlinepartE, mapC, mapE, levelsetC, levelsetE)

  write(dumstr,"(A3,I2.2,A14)") "Run",this%runID,"_levelsetC.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  call decomp_2d_write_one(1, levelsetC, fname, this%gpC)

  write(dumstr,"(A3,I2.2,A14)") "Run",this%runID,"_levelsetE.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  call decomp_2d_write_one(1, levelsetE, fname, this%gpE)


  !!!print *, 'Out mapC  min = ', -p_maxval(maxval(-mapC))
  !!!print *, 'Out mapC  max = ', p_maxval(maxval( mapC))
  !!!!print '(a,5(i4,1x))', 'Decomp info C = ', nrank, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)
  !!!!print '(a,5(i4,1x))', 'Decomp info E = ', nrank, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)

  call this%mark_ghost_points(gpC, nlayers, mapC, flagC, this%rbuffxC1, this%rbuffyC1, this%rbuffzC1);  this%num_gptsC = sum(flagC)
  write(dumstr,"(A3,I2.2,A9)") "Run",this%runID,"_mapC.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  this%rbuffxC1 = mapC
  call decomp_2d_write_one(1, this%rbuffxC1, fname, this%gpC)

  call this%mark_ghost_points(gpE, nlayers, mapE, flagE, this%rbuffxE1, this%rbuffyE1, this%rbuffzE1);  this%num_gptsE = sum(flagE)
  write(dumstr,"(A3,I2.2,A9)") "Run",this%runID,"_mapE.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  this%rbuffxE1 = mapE
  call decomp_2d_write_one(1, this%rbuffxE1, fname, this%gpE)
  !!!print *, '   flagC  min = ', -p_maxval(maxval(-flagC))
  !!!print *, '   flagC  max = ', p_maxval(maxval( flagC))

! comehere

  allocate(this%num_imptsC_arr(nproc), this%num_imptsE_arr(nproc))
  allocate(this%imptsC_index_st(nproc), this%imptsE_index_st(nproc))

  this%num_imptsC_arr = 0
  call MPI_AllGather(this%num_gptsC, 1, mpiinteg, this%num_imptsC_arr, 1, mpiinteg, MPI_COMM_WORLD, ierr)
  this%num_imptsC_glob = sum(this%num_imptsC_arr)
  this%imptsC_index_st(1) = 1
  do ii = 2, nproc
    this%imptsC_index_st(ii) = this%imptsC_index_st(ii-1)
  enddo

  this%num_imptsE_arr = 0
  call MPI_AllGather(this%num_gptsE, 1, mpiinteg, this%num_imptsE_arr, 1, mpiinteg, MPI_COMM_WORLD, ierr)
  this%num_imptsE_glob = sum(this%num_imptsE_arr)
  this%imptsE_index_st(1) = 1
  do ii = 2, nproc
    this%imptsE_index_st(ii) = this%imptsE_index_st(ii-1)
  enddo

  allocate(this%imptsC_xyz(this%num_imptsC_glob,3), this%imptsE_xyz(this%num_imptsE_glob,3))
  allocate(this%imptsC_xyz_loc(this%num_imptsC_glob,3), this%imptsE_xyz_loc(this%num_imptsE_glob,3))
  allocate(this%imptsC_u_tmp(this%num_imptsC_glob), this%imptsE_u_tmp(this%num_imptsE_glob))
  allocate(this%imptsC_v_tmp(this%num_imptsC_glob), this%imptsE_v_tmp(this%num_imptsE_glob))
  allocate(this%imptsC_w_tmp(this%num_imptsC_glob), this%imptsE_w_tmp(this%num_imptsE_glob))
     ! imptsC_numonproc; imptsC_indices, imptsC_multfac
  allocate(this%imptsC_numonproc(this%num_imptsC_glob), this%imptsE_numonproc(this%num_imptsE_glob))
  allocate(this%imptsC_indices(3,8,this%num_imptsC_glob), this%imptsE_indices(3,8,this%num_imptsE_glob))
  allocate(this%imptsC_multfac(  8,this%num_imptsC_glob), this%imptsE_multfac(  8,this%num_imptsE_glob))

  allocate(this%gptsC_xyz(this%num_gptsC,3), this%gptsE_xyz(this%num_gptsE,3))
  allocate(this%gptsC_ind(this%num_gptsC,3), this%gptsE_ind(this%num_gptsE,3))
  allocate(this%gptsC_bpt(this%num_gptsC,3), this%gptsE_bpt(this%num_gptsE,3))
  allocate(this%gptsC_bpind(this%num_gptsC), this%gptsE_bpind(this%num_gptsE))
  allocate(this%gptsC_img(this%num_gptsC,3), this%gptsE_img(this%num_gptsE,3))
  allocate(this%gptsC_bnp(this%num_gptsC,3), this%gptsE_bnp(this%num_gptsE,3))
  allocate(this%gptsC_dst(this%num_gptsC),   this%gptsE_dst(this%num_gptsE))

  !allocate(this%imptsC_numonproc(this%num_gptsC), this%imptsE_numonproc(this%num_gptsE))
  !allocate(this%imptsC_indices(3,8,this%num_gptsC), this%imptsE_indices(3,8,this%num_gptsE))
  !allocate(this%imptsC_multfac(  8,this%num_gptsC), this%imptsE_multfac(  8,this%num_gptsE))
  allocate(this%imptsC_u(this%num_imptsC_glob), this%imptsE_u(this%num_imptsE_glob))
  allocate(this%imptsC_v(this%num_imptsC_glob), this%imptsE_v(this%num_imptsE_glob))
  allocate(this%imptsC_w(this%num_imptsC_glob), this%imptsE_w(this%num_imptsE_glob))

  if(nrank==7) then
      write(*,'(a,i4,1x,3(e19.12,1x))') 'testing: ', nrank, xlinepart(32), ylinepart(4), zlinepart(18)
      write(*,'(a,i4,1x,100(e19.12,1x))') 'ylinepart: ', nrank, ylinepart(:)
      write(*,'(a,i4,1x,100(i4,1x))') 'mapC: ', nrank, mapC(32,:,18)
      write(*,'(a,i4,1x,100(i4,1x))') 'numptsC, numptsE: ', this%num_gptsC, this%num_gptsE
  endif

  call this%save_ghost_points(flagC, flagE, mapC, mapE, xlinepart, ylinepart, zlinepart, zlinepartE)

  call this%mark_boundary_points()

  call this%compute_image_points(Lx, Ly, zBot, zTop)

  dx = xlinepart(2)-xlinepart(1); dy = ylinepart(2)-ylinepart(1)
  call this%setup_interpolation(dx,dy,dz, zBot, xlinepart, ylinepart, zlinepart)

  deallocate(flagE, flagC, mapC, mapE, levelsetC, levelsetE, xlinepart, ylinepart, zlinepart, zlinepartE)

end subroutine

subroutine update_ibmgp(this, u, v, w, uE, vE, wC, uhat, vhat, what)
  class(ibmgp),     intent(inout) :: this
  real(rkind),    dimension(:,:,:), intent(inout) :: u, v, w, uE, vE, wC
  complex(rkind), dimension(:,:,:), intent(out)   :: uhat, vhat, what

  integer :: ii, itmp
  real(rkind) :: umax, umin, vmax, vmin, wmax, wmin, uEmax, uEmin, vEmax, vEmin, wCmax, wCmin

  !umax = p_maxval(u); umin = p_minval(u);  vmax = p_maxval(v); vmin = p_minval(v);  wmax = p_maxval(w); wmin = p_minval(w);
  !uEmax = p_maxval(uE); uEmin = p_minval(uE);  vEmax = p_maxval(vE); vEmin = p_minval(vE);  wCmax = p_maxval(wC); wCmin = p_minval(wC);
  !if(nrank==0) then
  !    print *, '-----Before interp_impts------'
  !    print '(a,6(e19.12,1x))', 'Cell uvw:', umax, umin, vmax, vmin, wCmax, wCmin
  !    print '(a,6(e19.12,1x))', 'Edge uvw:', uEmax, uEmin, vEmax, vEmin, wmax, wmin
  !    print *, '------------------------------'
  !endif
  call this%interp_imptsC(u, v,  wC)
  call this%interp_imptsE(uE, vE, w)

  !umax = p_maxval(u); umin = p_minval(u);  vmax = p_maxval(v); vmin = p_minval(v);  wmax = p_maxval(w); wmin = p_minval(w);
  !uEmax = p_maxval(uE); uEmin = p_minval(uE);  vEmax = p_maxval(vE); vEmin = p_minval(vE);  wCmax = p_maxval(wC); wCmin = p_minval(wC);
  !if(nrank==0) then
  !    print *, '-----After  interp_impts------'
  !    print '(a,6(e19.12,1x))', 'Cell uvw:', umax, umin, vmax, vmin, wCmax, wCmin
  !    print '(a,6(e19.12,1x))', 'Edge uvw:', uEmax, uEmin, vEmax, vEmin, wmax, wmin
  !    print *, '------------------------------'
  !endif
  call this%update_ghostptsCE(u, v, w)

  !umax = p_maxval(u); umin = p_minval(u);  vmax = p_maxval(v); vmin = p_minval(v);  wmax = p_maxval(w); wmin = p_minval(w);
  !uEmax = p_maxval(uE); uEmin = p_minval(uE);  vEmax = p_maxval(vE); vEmin = p_minval(vE);  wCmax = p_maxval(wC); wCmin = p_minval(wC);
  !if(nrank==0) then
  !    print *, '-----After  update_gpts ------'
  !    print '(a,6(e19.12,1x))', 'Cell uvw:', umax, umin, vmax, vmin, wCmax, wCmin
  !    print '(a,6(e19.12,1x))', 'Edge uvw:', uEmax, uEmin, vEmax, vEmin, wmax, wmin
  !    print *, '------------------------------'
  !endif
  call this%smooth_solidptsCE(u, v, w)

  !umax = p_maxval(u); umin = p_minval(u);  vmax = p_maxval(v); vmin = p_minval(v);  wmax = p_maxval(w); wmin = p_minval(w);
  !uEmax = p_maxval(uE); uEmin = p_minval(uE);  vEmax = p_maxval(vE); vEmin = p_minval(vE);  wCmax = p_maxval(wC); wCmin = p_minval(wC);
  !if(nrank==0) then
  !    print *, '-----After  smooth_solidpts---'
  !    print '(a,6(e19.12,1x))', 'Cell uvw:', umax, umin, vmax, vmin, wCmax, wCmin
  !    print '(a,6(e19.12,1x))', 'Edge uvw:', uEmax, uEmin, vEmax, vEmin, wmax, wmin
  !    print *, '------------------------------'
  !endif

  ! Step 3: Take it to spectral fields
  call this%spectC%fft(u, uhat)
  call this%spectC%fft(v, vhat)
  call this%spectE%fft(w, what)

end subroutine

subroutine smooth_solidptsCE(this, u, v, w)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(inout) :: u, v, w

  integer :: num_smooth, ii, jj, kk
  real(rkind) :: diffcoeff
  real(rkind) :: umax, umin, vmax, vmin, wmax, wmin
  real(rkind) :: uEmax, uEmin, vEmax, vEmin, wCmax, wCmin

  num_smooth = 10; diffcoeff = 1.0d-1;

  !if(nrank==8) then
  !  jj = 1; kk = 16;
  !  open(10,file='smooth_along_x_before.dat',status='replace')
  !  do ii = 1, size(u,1)
  !      write(10,*) u(ii,jj,kk), this%mask_solid_xC(ii,jj,kk)
  !  enddo
  !  close(10)
  !endif
  !umax = p_maxval(u); umin = p_minval(u);  vmax = p_maxval(v); vmin = p_minval(v);  wmax = p_maxval(w); wmin = p_minval(w);
  !!uEmax = p_maxval(uE);uEmin = p_minval(uE); vEmax = p_maxval(vE);vEmin = p_minval(vE); wCmax = p_maxval(wC);wCmin = p_minval(wC);
  !if(nrank==0) then
  !    print *, '-----Before Smoothing ---'
  !    print '(a,6(e19.12,1x))', 'CeEd uvw:', umax, umin, vmax, vmin, wmax, wmin
  !    !print '(a,6(e19.12,1x))', 'Edge uvw:', uEmax, uEmin, vEmax, vEmin, wmax, wmin
  !    print *, '------------------------------'
  !endif

  call this%smooth_along_x(u,             this%rbuffxC1, this%mask_solid_xC, num_smooth, diffcoeff)
  !umax = p_maxval(u); umin = p_minval(u)
  !if(nrank==0) then
  !    print *, '-----Smoothing u velocity ---'
  !    print *, '-----After  smooth_along_x---'
  !    print '(a,6(e19.12,1x))', 'Cell u:', umax, umin
  !    print *, '------------------------------'
  !endif
  !if(nrank==8) then
  !  jj = 1; kk = 16;
  !  open(10,file='smooth_along_x_after.dat',status='replace')
  !  do ii = 1, size(u,1)
  !      write(10,*) u(ii,jj,kk), this%mask_solid_xC(ii,jj,kk)
  !  enddo
  !endif

  call transpose_x_to_y   (u,             this%rbuffyC1, this%gpC)
  !if(nrank==8) then
  !  ii = 1; kk = 16;
  !  open(10,file='smooth_along_y_before.dat',status='replace')
  !  do jj = 1, size(this%rbuffyC1,2)
  !      write(10,*) this%rbuffyC1(ii,jj,kk), this%mask_solid_yC(ii,jj,kk)
  !  enddo
  !  close(10)
  !endif
  call this%smooth_along_y(this%rbuffyC1, this%rbuffyC2, this%mask_solid_yC, num_smooth, diffcoeff)
  !if(nrank==8) then
  !  ii = 1; kk = 16;
  !  open(10,file='smooth_along_y_after.dat',status='replace')
  !  do jj = 1, size(this%rbuffyC1,2)
  !      write(10,*) this%rbuffyC1(ii,jj,kk), this%mask_solid_yC(ii,jj,kk)
  !  enddo
  !  close(10)
  !endif
  !umax = p_maxval(this%rbuffyC1); umin = p_minval(this%rbuffyC1)
  !if(nrank==0) then
  !    print *, '-----After  smooth_along_y---'
  !    print '(a,6(e19.12,1x))', 'Cell u:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_y_to_z   (this%rbuffyC1, this%rbuffzC1, this%gpC)
  !if(nrank==8) then
  !  print *, 'nrank 8 x = ', this%gpC%zst(1), this%gpC%zen(1)
  !  ii = 1; jj = 33;
  !  open(10,file='smooth_along_z_before.dat',status='replace')
  !  do kk = 1, size(this%rbuffzC1,3)
  !      write(10,*) this%rbuffzC1(ii,jj,kk), this%mask_solid_zC(ii,jj,kk)
  !  enddo
  !  close(10)
  !endif
  call this%smooth_along_z(this%rbuffzC1, this%rbuffzC2, this%mask_solid_zC, num_smooth, diffcoeff)
  !if(nrank==8) then
  !  ii = 1; jj = 33;
  !  open(10,file='smooth_along_z_after.dat',status='replace')
  !  do kk = 1, size(this%rbuffzC1,3)
  !      write(10,*) this%rbuffzC1(ii,jj,kk), this%mask_solid_zC(ii,jj,kk)
  !  enddo
  !  close(10)
  !endif
  !umax = p_maxval(this%rbuffzC1); umin = p_minval(this%rbuffzC1)
  !if(nrank==0) then
  !    print *, '-----After  smooth_along_z---'
  !    print '(a,6(e19.12,1x))', 'Cell u:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_z_to_y   (this%rbuffzC1, this%rbuffyC1, this%gpC)
  call transpose_y_to_x   (this%rbuffyC1, u,             this%gpC)

  call this%smooth_along_x(v,             this%rbuffxC1, this%mask_solid_xC, num_smooth, diffcoeff)
  !umax = p_maxval(v); umin = p_minval(v)
  !if(nrank==0) then
  !    print *, '-----Smoothing v velocity ---'
  !    print *, '-----After  smooth_along_x---'
  !    print '(a,6(e19.12,1x))', 'Cell v:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_x_to_y   (v,             this%rbuffyC1, this%gpC)
  call this%smooth_along_y(this%rbuffyC1, this%rbuffyC2, this%mask_solid_yC, num_smooth, diffcoeff)
  !umax = p_maxval(this%rbuffyC1); umin = p_minval(this%rbuffyC1)
  !if(nrank==0) then
  !    print *, '-----After  smooth_along_y---'
  !    print '(a,6(e19.12,1x))', 'Cell v:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_y_to_z   (this%rbuffyC1, this%rbuffzC1, this%gpC)
  call this%smooth_along_z(this%rbuffzC1, this%rbuffzC2, this%mask_solid_zC, num_smooth, diffcoeff)
  !umax = p_maxval(this%rbuffzC1); umin = p_minval(this%rbuffzC1)
  !if(nrank==0) then
  !    print *, '-----After  smooth_along_z---'
  !    print '(a,6(e19.12,1x))', 'Cell v:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_z_to_y   (this%rbuffzC1, this%rbuffyC1, this%gpC)
  call transpose_y_to_x   (this%rbuffyC1, v,             this%gpC)

  call this%smooth_along_x(w,             this%rbuffxE1, this%mask_solid_xE, num_smooth, diffcoeff)
  !umax = p_maxval(w); umin = p_minval(w)
  !if(nrank==0) then
  !    print *, '-----Smoothing w velocity ---'
  !    print *, '-----After  smooth_along_x---'
  !    print '(a,6(e19.12,1x))', 'Edge w:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_x_to_y   (w,             this%rbuffyE1, this%gpE)
  call this%smooth_along_y(this%rbuffyE1, this%rbuffyE2, this%mask_solid_yE, num_smooth, diffcoeff)
  !umax = p_maxval(this%rbuffyE1); umin = p_minval(this%rbuffyE1)
  !if(nrank==0) then
  !    print *, '-----After  smooth_along_z---'
  !    print '(a,6(e19.12,1x))', 'Edge w:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_y_to_z   (this%rbuffyE1, this%rbuffzE1, this%gpE)
  call this%smooth_along_z(this%rbuffzE1, this%rbuffzE2, this%mask_solid_zE, num_smooth, diffcoeff)
  !umax = p_maxval(this%rbuffzE1); umin = p_minval(this%rbuffzE1)
  !if(nrank==0) then
  !    print *, '-----After  smooth_along_z---'
  !    print '(a,6(e19.12,1x))', 'Edge w:', umax, umin
  !    print *, '------------------------------'
  !endif
  call transpose_z_to_y   (this%rbuffzE1, this%rbuffyE1, this%gpE)
  call transpose_y_to_x   (this%rbuffyE1, w,             this%gpE)

  !call GracefulExit("Stopping here for now", 111)

end subroutine


subroutine smooth_along_x(this, uarr, rhs, mask_solid_x, num_smooth, diffcoeff)
  class(ibmgp),                  intent(in)    :: this
  real(rkind), dimension(:,:,:), intent(inout) :: uarr, rhs
  real(rkind), dimension(:,:,:), intent(in   ) :: mask_solid_x
  integer,                       intent(in)    :: num_smooth
  real(rkind),                   intent(in)    :: diffcoeff

  integer :: it, i, j, k, nx, nyloc, nzloc
  real(rkind) :: diffcoeff_x

  nx = size(uarr, 1); nyloc = size(uarr, 2); nzloc = size(uarr, 3)

  ! In x decomp; Smooth along x
  diffcoeff_x = diffcoeff ! * (dt/dx^2)

  do it = 1, num_smooth
      do k = 1, nzloc
        do j = 1, nyloc
          rhs(1,j,k) = (uarr(2,j,k) + uarr(nx,j,k) - two*uarr(1,j,k))
          rhs(2:nx-1,j,k) = (uarr(3:nx,j,k) + uarr(1:nx-2,j,k) - two*uarr(2:nx-1,j,k))
          rhs(nx,j,k) = (uarr(1,j,k) + uarr(nx-1,j,k) - two*uarr(nx,j,k))
        enddo
      enddo
      uarr = uarr + diffcoeff_x*rhs*mask_solid_x
  enddo

end subroutine

subroutine smooth_along_y(this, uarr, rhs, mask_solid_y, num_smooth, diffcoeff)
  class(ibmgp),                  intent(in)    :: this
  real(rkind), dimension(:,:,:), intent(inout) :: uarr, rhs
  real(rkind), dimension(:,:,:), intent(in   ) :: mask_solid_y
  integer,                       intent(in)    :: num_smooth
  real(rkind),                   intent(in)    :: diffcoeff

  integer :: it, i, j, k, nxloc, ny, nzloc
  real(rkind) :: diffcoeff_y

  nxloc = size(uarr, 1); ny = size(uarr, 2); nzloc = size(uarr, 3)

  ! In y decomp; Smooth along y
  diffcoeff_y = diffcoeff ! * (dt/dy^2)

  do it = 1, num_smooth
      do k = 1, nzloc
        rhs(:,1,k) = (uarr(:,2,k) + uarr(:,ny,k) - two*uarr(:,1,k))
        do j = 2, ny-1
          rhs(:,j,k) = (uarr(:,j+1,k) + uarr(:,j-1,k) - two*uarr(:,j,k))
        enddo
        rhs(:,ny,k) = (uarr(:,1,k) + uarr(:,ny-1,k) - two*uarr(:,ny,k))
      enddo
      uarr = uarr + diffcoeff_y*rhs*mask_solid_y
  enddo

end subroutine

subroutine smooth_along_z(this, uarr, rhs, mask_solid_z, num_smooth, diffcoeff)
  class(ibmgp),                  intent(in)    :: this
  real(rkind), dimension(:,:,:), intent(inout) :: uarr, rhs
  real(rkind), dimension(:,:,:), intent(in   ) :: mask_solid_z
  integer,                       intent(in)    :: num_smooth
  real(rkind),                   intent(in)    :: diffcoeff

  integer :: it, i, j, k, nxloc, nyloc, nz
  real(rkind) :: diffcoeff_z, dcmax, dcmin, rhsmax, rhsmin
  integer :: maxmsk, minmsk

  nxloc = size(uarr, 1); nyloc = size(uarr, 2); nz = size(uarr, 3)

  ! In z decomp; Smooth along z
  diffcoeff_z = diffcoeff ! * (dt/dz^2)

  do it = 1, num_smooth
      do k = 2, nz-1
          rhs(:,:,k) = (uarr(:,:,k-1) + uarr(:,:,k+1) - two*uarr(:,:,k))
      enddo
      rhs(:,:,1) = rhs(:,:,2)
      rhs(:,:,nz) = rhs(:,:,nz-1)
      uarr = uarr + diffcoeff_z*rhs*mask_solid_z
  enddo
  dcmax = p_maxval(diffcoeff_z); dcmin = p_minval(diffcoeff_z)
  rhsmax = p_maxval(rhs);    rhsmin = p_minval(rhs)
  maxmsk = p_maxval(maxval(mask_solid_z)); minmsk = p_minval(minval(mask_solid_z))
  !if(nrank==0) then
  !    print *, 'diffcoeff: ', dcmax, dcmin
  !    print *, 'rhs      : ', rhsmax, rhsmin
  !    print *, 'mask     : ', maxmsk, minmsk
  !endif

end subroutine

subroutine update_ghostptsCE(this, u, v, w)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(inout) :: u, v, w

  integer :: ii, jj, kk, i, j, k
  real(rkind) :: vec1(3), vec2(3), unrm(3), utan(3), dotpr, dist
  real(rkind) :: ibwallvel, utan_gp(3), unrm_gp(3), vmax, vmin


  if(this%ibwm==1) then
      ! no slip immersed boundary
      do ii = 1, this%num_gptsC
          i = this%gptsC_ind(ii,1);  j = this%gptsC_ind(ii,2);  k = this%gptsC_ind(ii,3)
          jj = this%imptsC_index_st(nrank) + ii - 1
          u(i,j,k) = -this%imptsc_u(jj)
          v(i,j,k) = -this%imptsc_v(jj)
      enddo
      do ii = 1, this%num_gptsE
          i = this%gptsE_ind(ii,1);   j = this%gptsE_ind(ii,2);  k = this%gptsE_ind(ii,3)
          jj = this%imptsE_index_st(nrank) + ii - 1
          w(i,j,k) = -this%imptsE_w(jj)
      enddo
  elseif(this%ibwm==2) then
      ! simple log-law immersed boundary
      do ii = 1, this%num_gptsC
          kk = this%imptsC_index_st(nrank) + ii - 1
          vec1(1) = this%imptsC_u(kk);  vec1(2) = this%imptsC_v(kk);  vec1(3) = this%imptsC_w(kk)

          jj = this%gptsC_bpind(ii)
          vec2(:) = this%surfnormal(jj,:)
          dotpr = sum(vec1*vec2)
          unrm = vec2*dotpr
          utan = vec1-unrm

          dist = this%gptsC_dst(ii) ! distance from ghost pt to image point
          !ibwallvel = this%ibwm_ustar/kappa*log(dist/this%ibwm_z0)
          !utan_gp = -utan+two*ibwallvel   ! (WM01)
           utan_gp = utan*(one - two*kappa/log(dist/this%ibwm_z0)) ! (WM04) 
          !utan_gp = utan  ! (WM03)
          !utan_gp = utan - two*this%ibwm_ustar/kappa ! (WM02) 
          unrm_gp = -unrm

          i = this%gptsC_ind(ii,1);  j = this%gptsC_ind(ii,2);  k = this%gptsC_ind(ii,3)
          u(i,j,k) = utan_gp(1) + unrm_gp(1)
          v(i,j,k) = utan_gp(2) + unrm_gp(2)
          !this%gptsC_w(ii) = utan_gp(3) + unrm_gp(3)
          !if(v(i,j,k) > 17.2d0 ) then
          !    print *, '++++++ ii = ', ii, '++++++'
          !    print '(a,3(e19.12,1x))', 'vec1   : ', vec1
          !    print '(a,3(e19.12,1x))', 'vec2   : ', vec2
          !    print '(a,3(e19.12,1x))', 'dotpr  : ', dotpr
          !    print '(a,3(e19.12,1x))', 'unrm   : ', unrm
          !    print '(a,3(e19.12,1x))', 'utan   : ', utan
          !    print '(a,3(e19.12,1x))', 'unrm_gp: ', unrm_gp
          !    print '(a,3(e19.12,1x))', 'utan_gp: ', utan_gp
          !    print '(a,3(e19.12,1x))', 'u, v   : ', utan_gp(1:2) + unrm_gp(1:2)
          !    print *, '++++++---------------------------'
          !endif
      enddo

      !vmax = p_maxval(v); vmin = p_minval(v)
      !print '(a,i4.4,1x,2(e19.12,1x))', '(vmax vmin): ', nrank, maxval(v), minval(v)

      do ii = 1, this%num_gptsE
          kk = this%imptsC_index_st(nrank) + ii - 1
          vec1(1) = this%imptsE_u(kk);  vec1(2) = this%imptsE_v(kk);  vec1(3) = this%imptsE_w(kk)

          jj = this%gptsE_bpind(ii)
          vec2(:) = this%surfnormal(jj,:)
          dotpr = sum(vec1*vec2)
          unrm = vec2*dotpr
          utan = vec1-unrm

          dist = this%gptsE_dst(ii) ! distance from ghost pt to image point
          !ibwallvel = this%ibwm_ustar/kappa*log(dist/this%ibwm_z0)
          !utan_gp = -utan+two*ibwallvel ! (WM01)
           utan_gp = utan*(one - two*kappa/log(dist/this%ibwm_z0)) ! (WM04) 
          !utan_gp = utan  ! (WM03)
          !utan_gp = utan - two*this%ibwm_ustar/kappa ! (WM02)
          unrm_gp = -unrm

          i = this%gptsE_ind(ii,1);   j = this%gptsE_ind(ii,2);  k = this%gptsE_ind(ii,3)
          !this%gptsE_u(ii) = utan_gp(1) + unrm_gp(1)
          !this%gptsE_v(ii) = utan_gp(2) + unrm_gp(2)
          w(i,j,k) = utan_gp(3) + unrm_gp(3)
      enddo

  endif

end subroutine

subroutine interp_imptsE(this, uE, vE, w)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: uE, vE, w

  integer :: ii, jj, i, j, k, ierr
  real(rkind) :: fac, tmpval1, tmpval2, tmpval3

  this%imptsE_u_tmp = zero;  this%imptsE_v_tmp = zero;  this%imptsE_w_tmp = zero
  do ii = 1, this%num_imptsE_glob ! this%num_gptsE
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
      this%imptsE_u_tmp(ii) = tmpval1
      this%imptsE_v_tmp(ii) = tmpval2
      this%imptsE_w_tmp(ii) = tmpval3
  enddo

  call MPI_AllReduce(this%imptsE_u_tmp, this%imptsE_u, this%num_imptsE_glob, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(this%imptsE_v_tmp, this%imptsE_v, this%num_imptsE_glob, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(this%imptsE_w_tmp, this%imptsE_w, this%num_imptsE_glob, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine

subroutine interp_imptsC(this, u, v, wC)
  class(ibmgp),     intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: u, v, wC

  integer :: ii, jj, i, j, k, ierr
  real(rkind) :: fac, tmpval1, tmpval2, tmpval3

  this%imptsC_u_tmp = zero;  this%imptsC_v_tmp = zero;  this%imptsC_w_tmp = zero
  do ii = 1, this%num_imptsC_glob ! this%num_gptsC
      ! how many of the 8 interpolation points are on this proc?
      tmpval1 = zero; tmpval2 = zero; tmpval3 = zero
      do jj = 1, this%imptsC_numonproc(ii)
          i   = this%imptsC_indices(1,jj,ii)
          j   = this%imptsC_indices(2,jj,ii)
          k   = this%imptsC_indices(3,jj,ii)
          fac = this%imptsC_multfac(jj, ii)
          tmpval1 = tmpval1 + fac* u(i,j,k)
          tmpval2 = tmpval2 + fac* v(i,j,k)
          tmpval3 = tmpval3 + fac*wC(i,j,k)
      enddo
      this%imptsC_u_tmp(ii) = tmpval1
      this%imptsC_v_tmp(ii) = tmpval2
      this%imptsC_w_tmp(ii) = tmpval3
  enddo

  call MPI_AllReduce(this%imptsC_u_tmp, this%imptsC_u, this%num_imptsC_glob, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(this%imptsC_v_tmp, this%imptsC_v, this%num_imptsC_glob, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(this%imptsC_w_tmp, this%imptsC_w, this%num_imptsC_glob, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

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

  do ii = 1, this%num_imptsC_glob  !this%num_gptsC
      ! first check if this point is on this processor
      !xloc = this%gptsC_img(ii,1);  yloc = this%gptsC_img(ii,2);  zloc = this%gptsC_img(ii,3);
      xloc = this%imptsC_xyz(ii,1);  yloc = this%imptsC_xyz(ii,2);  zloc = this%imptsC_xyz(ii,3);

      if( (xloc>=xdom_left) .and. (xloc<=xdom_right) .and. &
          (yloc>=ydom_left) .and. (yloc<=ydom_right) .and. &
          (zloc>=zdom_left) .and. (zloc<=zdom_right) ) then

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
      endif

      !this%gptsC_ileft(ii) = max(1, min(itmp+1, this%gpC%xsz(1)-1))
      !this%gptsC_ifacx(ii) = (this%gptsC_img(ii,1)-real(itmp, rkind)*dx)/dx

      !itmp = floor(this%gptsC_img(ii,2)/dy);
      !this%gptsC_jleft(ii) = max(1, min(itmp+1, this%gpC%ysz(2)-1))
      !this%gptsC_jfacy(ii) = (this%gptsC_img(ii,2)-real(itmp, rkind)*dy)/dy

      !itmp = floor((this%gptsC_img(ii,3)-zBot)/dz);
      !this%gptsC_kleft(ii) = max(1, min(itmp+1, this%gpC%zsz(3)-1))
      !this%gptsC_kfacz(ii) = (this%gptsC_img(ii,3)-real(itmp, rkind)*dz - zBot)/dz
  enddo

  do ii = 1, this%num_imptsE_glob ! this%num_gptsE
      ! first check if this point is on this processor
      !xloc = this%gptsE_img(ii,1);  yloc = this%gptsE_img(ii,2);  zloc = this%gptsE_img(ii,3);
      xloc = this%imptsC_xyz(ii,1);  yloc = this%imptsC_xyz(ii,2);  zloc = this%imptsC_xyz(ii,3);

      if( (xloc>=xdom_left) .and. (xloc<=xdom_right) .and. &
          (yloc>=ydom_left) .and. (yloc<=ydom_right) .and. &
          (zloc>=zdom_left) .and. (zloc<=zdom_right) ) then

        itmp = floor((xloc-xdom_left)/dx);        
        jtmp = floor((yloc-ydom_left)/dy);        
        ktmp = floor((zloc-zdom_left)/dz);        

        ! determine facx, facy, facz
        if(itmp==0) then
            facx = one - (xlinepart(1)-xloc)/dx
        else
            facx = one - (xloc-xlinepart(itmp))/dx
        endif

        if(jtmp==5) then
            print *, '-----nrank = ', nrank 
            print '(a,2(i5,1x),3(e19.12,1x))', '--nrank: ', nrank, ii, xloc, yloc, zloc
            print '(a,1(i5,1x),3(e19.12,1x))', '--xdom : ', itmp, xdom_left, xdom_right, dx
            print '(a,1(i5,1x),3(e19.12,1x))', '--ydom : ', jtmp, ydom_left, ydom_right, dy
            print '(a,1(i5,1x),3(e19.12,1x))', '--zdom : ', ktmp, zdom_left, zdom_right, dz
            print '(a,       100(e19.12,1x))', '--xline: ', xlinepart
            print '(a,       100(e19.12,1x))', '--yline: ', ylinepart
            print '(a,       100(e19.12,1x))', '--zline: ', zlinepart
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

      endif
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

subroutine compute_image_points(this, Lx, Ly, zBot, zTop)
  class(ibmgp),     intent(inout) :: this
  real(rkind), intent(in) :: Lx, Ly, zBot, zTop

  integer :: ii, jj, kk, ierr
  real(rkind) :: vec1(3), vec2(3), dotpr, xmax_img, xmin_img, ymax_img, ymin_img, zmax_img, zmin_img
  character(len=clen) :: fname, dumstr

  kk = this%imptsC_index_st(nrank)-1

  do ii = 1, this%num_gptsC
      ! compute image of gptsC_xyz(ii,:) wrt gptsc_bpt(ii,:) and
      vec1(:) = this%gptsC_bpt(ii,:) - this%gptsC_xyz(ii,:)
      jj = this%gptsC_bpind(ii)
      vec2(:) = this%surfnormal(jj,:)
      dotpr = sum(vec1*vec2)

      this%gptsC_img(kk,:) = this%gptsC_xyz(ii,:) + two*dotpr*vec2
      this%gptsC_bnp(ii,:) = this%gptsC_xyz(ii,:) + dotpr*vec2
      this%gptsC_dst(ii)   = dotpr

      ! store in the global array
      kk = kk+1
      this%imptsC_xyz_loc(kk,:) = this%gptsC_img(kk,:)
  enddo

  do ii = 1, this%num_gptsE
      ! compute image of gptsC_xyz(ii,:) wrt gptsc_bpt(ii,:) and
      vec1(:) = this%gptsE_bpt(ii,:) - this%gptsE_xyz(ii,:)
      jj = this%gptsE_bpind(ii)
      vec2(:) = this%surfnormal(jj,:)
      dotpr = sum(vec1*vec2)

      this%gptsE_img(ii,:) = this%gptsE_xyz(ii,:) + two*dotpr*vec2
      this%gptsE_bnp(ii,:) = this%gptsE_xyz(ii,:) + dotpr*vec2
      this%gptsE_dst(ii)   = dotpr

      ! store in the global array
      kk = kk+1
      this%imptsE_xyz_loc(kk,:) = this%gptsE_img(kk,:)
  enddo

  !create global image point arrays; this%imptsC_glob(index,3);
  !this%imptsE_glob(index,3)
  call MPI_AllReduce(this%imptsC_xyz_loc, this%imptsC_xyz, this%num_imptsC_glob*3, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(this%imptsE_xyz_loc, this%imptsE_xyz, this%num_imptsE_glob*3, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

  write(dumstr,"(A3,I2.2,A15,I4.4,A4)") "Run",this%runID,"_ibm_gptsC_img_",nrank,".dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  open(11,file=fname,status='unknown',action='write')
  do ii = 1, this%num_gptsC
      write(11, '(7(e19.12,1x))') this%gptsC_img(ii,:), this%gptsC_bnp(ii,:), this%gptsC_dst(ii)
  enddo
  close(11)

  write(dumstr,"(A3,I2.2,A15,I4.4,A4)") "Run",this%runID,"_ibm_gptsE_img_",nrank,".dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  open(11,file=fname,status='unknown',action='write')
  do ii = 1, this%num_gptsE
      write(11, '(7(e19.12,1x))') this%gptsE_img(ii,:), this%gptsE_bnp(ii,:), this%gptsE_dst(ii)
  enddo
  close(11)

  xmax_img = p_maxval(maxval(this%gptsC_img(:,1)));   xmin_img = p_minval(minval(this%gptsC_img(:,1)))
  ymax_img = p_maxval(maxval(this%gptsC_img(:,2)));   ymin_img = p_minval(minval(this%gptsC_img(:,2)))
  zmax_img = p_maxval(maxval(this%gptsC_img(:,3)));   zmin_img = p_minval(minval(this%gptsC_img(:,3)))
  if((xmax_img > Lx) .or. (xmin_img < zero)) then
      call GracefulExit("Cell image points x coordinates outside domain", 111)
  endif
  if((ymax_img > Ly) .or. (ymin_img < zero)) then
      call GracefulExit("Cell image points y coordinates outside domain", 111)
  endif
  if((zmax_img > zTop) .or. (zmin_img < zBot)) then
      call GracefulExit("Cell image points z coordinates outside domain", 111)
  endif

  xmax_img = p_maxval(maxval(this%gptsE_img(:,1)));   xmin_img = p_minval(minval(this%gptsE_img(:,1)))
  ymax_img = p_maxval(maxval(this%gptsE_img(:,2)));   ymin_img = p_minval(minval(this%gptsE_img(:,2)))
  zmax_img = p_maxval(maxval(this%gptsE_img(:,3)));   zmin_img = p_minval(minval(this%gptsE_img(:,3)))
  if((xmax_img > Lx) .or. (xmin_img < zero)) then
      call GracefulExit("Edge image points x coordinates outside domain", 111)
  endif
  if((ymax_img > Ly) .or. (ymin_img < zero)) then
      call GracefulExit("Edge image points y coordinates outside domain", 111)
  endif
  if((zmax_img > zTop) .or. (zmin_img < zBot)) then
      call GracefulExit("Edge image points z coordinates outside domain", 111)
  endif


end subroutine

subroutine mark_boundary_points(this)
  class(ibmgp),     intent(inout) :: this

  integer :: ii, imin
  real(rkind), allocatable, dimension(:) :: distfn
  character(len=clen) :: fname, dumstr

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

  write(dumstr,"(A3,I2.2,A15,I4.4,A4)") "Run",this%runID,"_ibm_gptsC_bpt_",nrank,".dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  open(11,file=fname,status='unknown',action='write')
  do ii = 1, this%num_gptsC
      write(11, '(3(e19.12,1x), 3(i5,1x))') this%gptsC_bpt(ii,:), this%gptsC_bpind(ii)
  enddo
  close(11)

  write(dumstr,"(A3,I2.2,A15,I4.4,A4)") "Run",this%runID,"_ibm_gptsE_bpt_",nrank,".dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  open(11,file=fname,status='unknown',action='write')
  do ii = 1, this%num_gptsE
      write(11, '(3(e19.12,1x), 3(i5,1x))') this%gptsE_bpt(ii,:), this%gptsE_bpind(ii)
  enddo
  close(11)

end subroutine

subroutine save_ghost_points(this, flagC, flagE, mapC, mapE,  xlinepart, ylinepart, zlinepart, zlinepartE)
  class(ibmgp),     intent(inout) :: this
  integer, dimension(:,:,:), intent(in)   :: flagC, flagE, mapC, mapE
  real(rkind), dimension(:), intent(in) :: xlinepart, ylinepart, zlinepart, zlinepartE

  integer :: i, j, k, ii
  character(len=clen)    :: fname, dumstr

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

  ! fill mask_solid
  this%mask_solid_xC = zero
  do k = 1, this%gpC%xsz(3)
    do j = 1, this%gpC%xsz(2)
      do i = 1, this%gpC%xsz(1)
          if( (flagC(i,j,k)==0) .and. (mapC(i,j,k)==0) ) then
              this%mask_solid_xC(i,j,k) = one
          endif
      enddo
    enddo
  enddo
  call transpose_x_to_y(this%mask_solid_xC, this%mask_solid_yC, this%gpC)
  call transpose_y_to_z(this%mask_solid_yC, this%mask_solid_zC, this%gpC)

  this%mask_solid_xE = zero
  do k = 1, this%gpE%xsz(3)
    do j = 1, this%gpE%xsz(2)
      do i = 1, this%gpE%xsz(1)
          if( (flagE(i,j,k)==0) .and. (mapE(i,j,k)==0) ) then
              this%mask_solid_xE(i,j,k) = one
          endif
      enddo
    enddo
  enddo
  call transpose_x_to_y(this%mask_solid_xE, this%mask_solid_yE, this%gpE)
  call transpose_y_to_z(this%mask_solid_yE, this%mask_solid_zE, this%gpE)

  write(dumstr,"(A3,I2.2,A11,I4.4,A4)") "Run",this%runID,"_ibm_gptsC_",nrank,".dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  open(11,file=fname,status='unknown',action='write')
  do ii = 1, this%num_gptsC
      write(11, '(3(e19.12,1x), 3(i5,1x))') this%gptsC_xyz(ii,:), this%gptsC_ind(ii,:)
  enddo
  close(11)
  write(dumstr,"(A3,I2.2,A11,I4.4,A4)") "Run",this%runID,"_ibm_gptsE_",nrank,".dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  open(11,file=fname,status='unknown',action='write')
  do ii = 1, this%num_gptsE
      write(11, '(3(e19.12,1x), 3(i5,1x))') this%gptsE_xyz(ii,:), this%gptsE_ind(ii,:)
  enddo
  close(11)
  !if(nrank==0) then
      print ('(a,2(i5,1x))'), '--Done writing ghost points information to file--', nrank, this%num_gptsC
  !endif

  ! Write all mask fields to file
  write(dumstr,"(A3,I2.2,A12)") "Run",this%runID,"_mask_xC.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  call decomp_2d_write_one(1, this%mask_solid_xC, fname, this%gpC)

  write(dumstr,"(A3,I2.2,A12)") "Run",this%runID,"_mask_yC.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  call decomp_2d_write_one(2, this%mask_solid_yC, fname, this%gpC)

  write(dumstr,"(A3,I2.2,A12)") "Run",this%runID,"_mask_zC.dat"
  fname = this%outputDir(:len_trim(this%outputDir))//"/"//trim(dumstr)
  call decomp_2d_write_one(3, this%mask_solid_zC, fname, this%gpC)


end subroutine

subroutine mark_ghost_points(this, gp, nlayers, map, flag, rbuffx1, rbuffy1, rbuffz1)
  class(ibmgp),     intent(in) :: this
  type(decomp_info), pointer, intent(in) :: gp
  integer, intent(in) :: nlayers
  integer, dimension(:,:,:), intent(in) :: map
  integer, dimension(:,:,:), intent(out) :: flag
  real(rkind), dimension(:,:,:), pointer, intent(inout) :: rbuffx1, rbuffy1, rbuffz1

  integer :: i, j, k, ii, jj
  integer, allocatable, dimension(:) :: locflags
  integer, allocatable, dimension(:,:,:) :: map_yd, map_zd, flag_yd, flag_zd

  flag = 0

  allocate(locflags(2*nlayers))
  allocate(map_yd(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
  allocate(map_zd(gp%zsz(1), gp%zsz(2), gp%zsz(3)))
  allocate(flag_yd(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
  allocate(flag_zd(gp%zsz(1), gp%zsz(2), gp%zsz(3)))

  print '(a,5(i4,1x))', '-----mark_ghost_points x decomp---', nrank, gp%xsz(1), gp%xsz(2), gp%xsz(3)
  print '(a,5(i4,1x))', '-----mark_ghost_points y decomp---', nrank, gp%ysz(1), gp%ysz(2), gp%ysz(3)
  print '(a,5(i4,1x))', '-----mark_ghost_points z decomp---', nrank, gp%zsz(1), gp%zsz(2), gp%zsz(3)
  print '(a,5(i4,1x))', '-----size array map_yd         ---', nrank, size(map_yd,1), size(map_yd,2), size(map_yd,3)
  print '(a,5(i4,1x))', '-----size array map_zd         ---', nrank, size(map_zd,1), size(map_zd,2), size(map_zd,3)

  ! mark ghost points 
  ! first check in x direction
  do k = 1, gp%xsz(3)
    do j = 1, gp%xsz(2)
      do i = 1+nlayers, gp%xsz(1)-nlayers
          !if(nrank==7 .and. i==32 .and. k==18) then
          !    print '(a,3(i4,1x))', 'In mark_ghost_points', j, map(i,j,k), nlayers
          !endif
          if(map(i,j,k)==0) then
             ! point (i,j,k) is solid; check its x neighbours
             ii = 0; locflags = 0
             do jj = 1, nlayers
                 ii = ii+1;   locflags(ii) = map(i+jj, j, k)
                 ii = ii+1;   locflags(ii) = map(i-jj, j, k)
             enddo 
             if(sum(locflags) > 0) flag(i,j,k) = 1
          endif
      enddo
    enddo
  enddo

  ! now check in y direction
  rbuffx1 = map;   call transpose_x_to_y(rbuffx1,  rbuffy1,  gp);  map_yd = rbuffy1
  rbuffx1 = flag;  call transpose_x_to_y(rbuffx1,  rbuffy1,  gp); flag_yd = rbuffy1
  do k = 1, gp%ysz(3)
    do j = 1+nlayers, gp%ysz(2)-nlayers
      do i = 1, gp%ysz(1)
          !if(nrank==7 .and. i==32 .and. k==18) then
          !    print '(a,3(i4,1x))', 'In mark_ghost_points', j, map(i,j,k), nlayers
          !endif
          if(map_yd(i,j,k)==0) then
             ! point (i,j,k) is solid; check its x neighbours
             ii = 0; locflags = 0
             do jj = 1, nlayers
                 ii = ii+1;   locflags(ii) = map_yd(i, j+jj, k)
                 ii = ii+1;   locflags(ii) = map_yd(i, j-jj, k)
             enddo 
             if(sum(locflags) > 0) flag_yd(i,j,k) = flag_yd(i,j,k) + 1
          endif
      enddo
    enddo
  enddo

  ! now check in z direction
  rbuffy1 = map_yd;   call transpose_y_to_z(rbuffy1,  rbuffz1,  gp);  map_zd = rbuffz1
  rbuffy1 = flag_yd;  call transpose_y_to_z(rbuffy1,  rbuffz1,  gp); flag_zd = rbuffz1
  do k = 1+nlayers, gp%zsz(3)-nlayers
    do j = 1, gp%zsz(2)
      do i = 1, gp%zsz(1)
          !if(nrank==7 .and. i==32 .and. k==18) then
          !    print '(a,3(i4,1x))', 'In mark_ghost_points', j, map(i,j,k), nlayers
          !endif
          if(map_zd(i,j,k)==0) then
             ! point (i,j,k) is solid; check its x neighbours
             ii = 0; locflags = 0
             do jj = 1, nlayers
                 ii = ii+1;   locflags(ii) = map_zd(i, j, k+jj)
                 ii = ii+1;   locflags(ii) = map_zd(i, j, k-jj)
             enddo 
             if(sum(locflags) > 0) flag_zd(i,j,k) = flag_zd(i,j,k) + 1
          endif
      enddo
    enddo
  enddo

  rbuffz1 = flag_zd;   call transpose_z_to_y(rbuffz1,  rbuffy1,  gp)
                       call transpose_y_to_x(rbuffy1,  rbuffx1,  gp);   flag = rbuffx1

  flag = min(flag, 1)

  deallocate(locflags, flag_yd, flag_zd, map_yd, map_zd)

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
          xk(2) = ylinepart(j)-this%surfcent(jj,2)
          xk(3) = zlinepart(k)-this%surfcent(jj,3)
          xkmag = sqrt(sum(xk*xk))
          levelsetC(i,j,k) = sign(1.0_rkind, dot_product(xk, this%surfnormal(jj,:)))*xkmag
      enddo
    enddo
  enddo
  where (levelsetC > zero)
    mapC = 1
  endwhere
  print *, 'levelsetC max = ', p_maxval(levelsetC)
  print *, 'levelsetC min = ', p_minval(levelsetC)
  print *, '    mapC  min = ', -p_maxval(maxval(-mapC))
  print *, '    mapC  max = ', p_maxval(maxval( mapC))

  ! compute levelsetE using lnearlistE
  call probe_nearest_points_kdtree_ib(1, this%gpE%xsz(1), 1, this%gpE%xsz(2), 1, this%gpE%xsz(3), &
                                      xlinepart, ylinepart, zlinepartE, lnearlistE)
  mapE = 0
  do k = 1, this%gpE%xsz(3)
    do j = 1, this%gpE%xsz(2)
      do i = 1, this%gpE%xsz(1)
          jj = lnearlistE(i,j,k)
          xk(1) = xlinepart(i)-this%surfcent(jj,1)
          xk(2) = ylinepart(j)-this%surfcent(jj,2)
          xk(3) = zlinepartE(k)-this%surfcent(jj,3)
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

  deallocate(this%mask_solid_xC, this%mask_solid_yC, this%mask_solid_zC)
  deallocate(this%mask_solid_xE, this%mask_solid_yE, this%mask_solid_zE)
  nullify(this%rbuffxC1, this%rbuffxC2)
  nullify(this%rbuffyC1, this%rbuffyC2)
  nullify(this%rbuffzC1, this%rbuffzC2)
  nullify(this%rbuffxE1, this%rbuffxE2)
  nullify(this%rbuffyE1, this%rbuffyE2)
  nullify(this%rbuffzE1, this%rbuffzE2)
  !deallocate(this%gptsC_ifacx, this%gptsE_ifacx)
  !deallocate(this%gptsC_jfacy, this%gptsE_jfacy)
  !deallocate(this%gptsC_kfacz, this%gptsE_kfacz)
  !deallocate(this%gptsC_ileft, this%gptsE_ileft)
  !deallocate(this%gptsC_jleft, this%gptsE_jleft)
  !deallocate(this%gptsC_kleft, this%gptsE_kleft)
  !deallocate(this%gptsC_u,  this%gptsE_u)
  !deallocate(this%gptsC_v,  this%gptsE_v)
  !deallocate(this%gptsC_w,  this%gptsE_w)
  deallocate(this%num_imptsC_arr, this%num_imptsE_arr)
  deallocate(this%imptsC_index_st, this%imptsE_index_st)
  deallocate(this%imptsC_xyz, this%imptsE_xyz)
  deallocate(this%imptsC_xyz_loc, this%imptsE_xyz_loc)
  deallocate(this%imptsC_u_tmp, this%imptsE_u_tmp)
  deallocate(this%imptsC_v_tmp, this%imptsE_v_tmp)
  deallocate(this%imptsC_w_tmp, this%imptsE_w_tmp)
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
