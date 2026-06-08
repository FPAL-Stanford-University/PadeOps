module MultiBlockTopologyMod

    use kind_parameters, only: rkind,clen
    use constants,       only: half,one,zero
    use exits,           only: GracefulExit
    use reductions,      only: P_MAXVAL, P_MINVAL
    use decomp_2d
    use decomp_2d_io

    implicit none

    !private

    !! maximum of 10 domain/processor blocks allowed for now. Increase if needed.
    integer, parameter :: max_numbl = 10
    public :: multiblocktopol

    type :: multiblocktopol

        !private 
        integer :: x_num_blocks, y_num_blocks, z_num_blocks
        integer, allocatable, dimension(:,:) :: xst, yst, zst, xen, yen, zen
        integer, allocatable, dimension(:,:) :: ybclo_st, ybclo_en, ybchi_st, ybchi_en
        integer, allocatable, dimension(:,:) :: y_intbd_left_st, y_intbd_left_en
        integer, allocatable, dimension(:,:) :: y_intbd_rght_st, y_intbd_rght_en
        integer, allocatable, dimension(:)   :: y_num_intbd_left, y_num_intbd_rght
        real(rkind), allocatable, dimension(:,:,:) :: mask
        real(rkind), dimension(3,max_numbl) :: x_dombl_pt1, x_dombl_pt2   !! left and right ends of the domain blocks
        real(rkind), dimension(3,max_numbl) :: y_dombl_pt1, y_dombl_pt2
        real(rkind), dimension(3,max_numbl) :: z_dombl_pt1, z_dombl_pt2
        
        contains

            procedure :: init
            procedure :: destroy

    end type

contains

  subroutine init(this, decomp, mesh, inputfile, xbuf, zbuf)
      class(multiblocktopol), intent(inout) :: this
      class(decomp_info), intent(in), target :: decomp
      real(rkind), dimension(:,:,:,:), intent(in) :: mesh
      character(len=clen), intent(in) :: inputfile
      real(rkind), target, intent(in), dimension(:,:,:,:) :: xbuf, zbuf

      integer :: imb, ierr, nxl, nyl, nzl, indices_pt1(3), indices_pt2(3)
      character(len=clen) :: err_message, outdebugfile

      integer :: x_dom_numbl=1, y_dom_numbl=1, z_dom_numbl=1            !! number of domain blocks 
      real(rkind), dimension(3,max_numbl) :: x_dombl_pt1, x_dombl_pt2   !! left and right ends of the domain blocks
      real(rkind), dimension(3,max_numbl) :: y_dombl_pt1, y_dombl_pt2
      real(rkind), dimension(3,max_numbl) :: z_dombl_pt1, z_dombl_pt2
      integer,     dimension(3,max_numbl) :: xst, yst, zst, xen, yen, zen   !! for temporary storage
      integer,     dimension(3,max_numbl) :: ybclo_st, ybclo_en, ybchi_st, ybchi_en !! temporary storage for bc
      integer,     dimension(3,max_numbl) :: y_intbd_left_st, y_intbd_left_en !! temporary storage for internal (left/right) boundaries
      integer,     dimension(max_numbl)   :: y_num_intbd_left, y_num_intbd_rght !! number of internal (left/right) boundaries
      integer,     dimension(3,max_numbl) :: y_intbd_rght_st, y_intbd_rght_en !! temporary storage for internal (left/right) boundaries
      integer :: x_num_blocks=0, y_num_blocks=0, z_num_blocks=0             !! for temporary storage
      real(rkind), allocatable, dimension(:) :: xline_x, yline_x, zline_x
      real(rkind), allocatable, dimension(:) :: xline_y, yline_y, zline_y
      real(rkind), allocatable, dimension(:) :: xline_z, yline_z, zline_z
      integer :: iounit=123, i1, i2, j1, j2, k1, k2, jj, kk, i_intbd, imb1, imb2
      logical :: intersection_exists, detailed_debug = .false., file_exists
      real(rkind), dimension(3) :: intrbl_pt1, intrbl_pt2, x_procbl_pt1, x_procbl_pt2
      real(rkind), dimension(3) :: y_procbl_pt1, y_procbl_pt2, z_procbl_pt1, z_procbl_pt2
      real(rkind), dimension(:,:,:), pointer :: xtmp, ztmp

      !! domain block :: `global'
      namelist /MULTIBLOCK/ x_dom_numbl, y_dom_numbl, z_dom_numbl, &
                            x_dombl_pt1, x_dombl_pt2, y_dombl_pt1, &
                            y_dombl_pt2, z_dombl_pt1, z_dombl_pt2, &
                            detailed_debug
      !!---default---
      !!! only 1 block by default
      this%x_dombl_pt1 = x_dombl_pt1
      this%x_dombl_pt2 = x_dombl_pt2
      this%y_dombl_pt1 = y_dombl_pt1
      this%y_dombl_pt2 = y_dombl_pt2
      this%z_dombl_pt1 = z_dombl_pt1
      this%z_dombl_pt2 = z_dombl_pt2
      !!this%y_num_blocks = 1
      !!this%z_num_blocks = 1

      !!allocate(this%xst(3,this%x_num_blocks))
      !!allocate(this%xen(3,this%x_num_blocks))
      !!allocate(this%yst(3,this%y_num_blocks))
      !!allocate(this%yen(3,this%y_num_blocks))
      !!allocate(this%zst(3,this%z_num_blocks))
      !!allocate(this%zen(3,this%z_num_blocks))

      !!!! -- needed?? --this%decomp => decomp

      !!! the block spans the entire domain by default
      !!! x-decomp
      !!do ib = 1, this%x_num_blocks
      !!  this%xst(1, ib) = 1;          this%xen(1, ib) = decomp%xsz(1)
      !!  this%xst(2, ib) = 1;          this%xen(2, ib) = decomp%xsz(2)
      !!  this%xst(3, ib) = 1;          this%xen(3, ib) = decomp%xsz(3)
      !!end do

      !!! y-decomp
      !!do ib = 1, this%y_num_blocks
      !!  this%yst(1, ib) = 1;          this%yen(1, ib) = decomp%ysz(1)
      !!  this%yst(2, ib) = 1;          this%yen(2, ib) = decomp%ysz(2)
      !!  this%yst(3, ib) = 1;          this%yen(3, ib) = decomp%ysz(3)
      !!end do

      !!! z-decomp
      !!do ib = 1, this%z_num_blocks
      !!  this%zst(1, ib) = 1;          this%zen(1, ib) = decomp%zsz(1)
      !!  this%zst(2, ib) = 1;          this%zen(2, ib) = decomp%zsz(2)
      !!  this%zst(3, ib) = 1;          this%zen(3, ib) = decomp%zsz(3)
      !!end do
      !!---end default---

      ! STEP 1 :: Read Inputs
      open(unit=iounit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=iounit, NML=MULTIBLOCK)
      close(iounit)

      ! STEP 2 :: Exit if domain blocks (read from inputs) are not correct
      !! check if the x_dom_numbl <= 10. If not, exit. Similarly for y, z
      if( (x_dom_numbl > max_numbl) .or. (y_dom_numbl > max_numbl) .or. (z_dom_numbl > max_numbl)) then
          call GracefulExit("dom_numbl (no. of domain blocks) in each direction must be < max_numbl (usually 10).", 11)
      endif

      !! for each domain block, x_dombl_pt1(1:3) should be less than
      !or equal to x_dombl_pt2(1:3). Otherwise, exit. Similarly for y, z
      do imb = 1, x_dom_numbl
          if( (x_dombl_pt1(1, imb) >= x_dombl_pt2(1, imb)) .or. &
              (x_dombl_pt1(2, imb) >= x_dombl_pt2(2, imb)) .or. &
              (x_dombl_pt1(3, imb) >= x_dombl_pt2(3, imb)) ) then
            write(err_message, '(a,1x,i10)') "Error in x_dombl_pt block num", imb
            call GracefulExit(err_message, 11)
          endif
      enddo

      do imb = 1, y_dom_numbl
          if( (y_dombl_pt1(1, imb) >= y_dombl_pt2(1, imb)) .or. &
              (y_dombl_pt1(2, imb) >= y_dombl_pt2(2, imb)) .or. &
              (y_dombl_pt1(3, imb) >= y_dombl_pt2(3, imb)) ) then
            write(err_message, '(a,1x,i10)') "Error in y_dombl_pt block num", imb
            call GracefulExit(err_message, 11)
          endif
      enddo

      do imb = 1, z_dom_numbl
          if( (z_dombl_pt1(1, imb) >= z_dombl_pt2(1, imb)) .or. &
              (z_dombl_pt1(2, imb) >= z_dombl_pt2(2, imb)) .or. &
              (z_dombl_pt1(3, imb) >= z_dombl_pt2(3, imb)) ) then
            write(err_message, '(a,1x,i10)') "Error in z_dombl_pt block num", imb
            call GracefulExit(err_message, 11)
          endif
      enddo

      !! STEP 3  :: Define processor block
      xtmp => xbuf(:,:,:,1);    ztmp => zbuf(:,:,:,1)

      !! STEP 3a :: for x-decomposition
      allocate(xline_x(decomp%xsz(1)),  yline_x(decomp%xsz(2)),   zline_x(decomp%xsz(3)))
      call transpose_y_to_x(mesh(:,:,:,1), xtmp, decomp);  xline_x = xtmp(:,1,1)
      call transpose_y_to_x(mesh(:,:,:,2), xtmp, decomp);  yline_x = xtmp(1,:,1)
      call transpose_y_to_x(mesh(:,:,:,3), xtmp, decomp);  zline_x = xtmp(1,1,:)

      !print *, '++==', nrank, decomp%xsz(1), size(xline_x)
      x_procbl_pt1(1) = xline_x(1);      x_procbl_pt2(1) = xline_x(decomp%xsz(1))
      x_procbl_pt1(2) = yline_x(1);      x_procbl_pt2(2) = yline_x(decomp%xsz(2))
      x_procbl_pt1(3) = zline_x(1);      x_procbl_pt2(3) = zline_x(decomp%xsz(3))
      !! STEP 3a :: Done for x-decomposition

      !! STEP 3b :: for y-decomposition
      allocate(xline_y(decomp%ysz(1)),  yline_y(decomp%ysz(2)),   zline_y(decomp%ysz(3)))
      xline_y = mesh(:,1,1,1);          yline_y = mesh(1,:,1,2);  zline_y = mesh(1,1,:,3)

      y_procbl_pt1(1) = xline_y(1);      y_procbl_pt2(1) = xline_y(decomp%ysz(1))
      y_procbl_pt1(2) = yline_y(1);      y_procbl_pt2(2) = yline_y(decomp%ysz(2))
      y_procbl_pt1(3) = zline_y(1);      y_procbl_pt2(3) = zline_y(decomp%ysz(3))
      !! STEP 3b :: Done for y-decomposition

      !! STEP 3c :: for z-decomposition
      allocate(xline_z(decomp%zsz(1)),  yline_z(decomp%zsz(2)),   zline_z(decomp%zsz(3)))
      call transpose_y_to_z(mesh(:,:,:,1), ztmp, decomp); xline_z = ztmp(:,1,1)
      call transpose_y_to_z(mesh(:,:,:,2), ztmp, decomp); yline_z = ztmp(1,:,1)
      call transpose_y_to_z(mesh(:,:,:,3), ztmp, decomp); zline_z = ztmp(1,1,:)

      z_procbl_pt1(1) = xline_z(1);      z_procbl_pt2(1) = xline_z(decomp%zsz(1))
      z_procbl_pt1(2) = yline_z(1);      z_procbl_pt2(2) = yline_z(decomp%zsz(2))
      z_procbl_pt1(3) = zline_z(1);      z_procbl_pt2(3) = zline_z(decomp%zsz(3))
      !! STEP 3c :: Done for z-decomposition

      !! STEP 4  :: Set up multiblock topology
      !! STEP 4a :: for x-decomposition
      if(detailed_debug) then
          iounit = 100+nrank
          write(outdebugfile, '(a,i6.6,a)') 'debuginfo_', nrank,'.dat'
          inquire(file=outdebugfile, exist=file_exists)
          if(file_exists) then
            open(unit=iounit, file=outdebugfile, form='formatted', status='old', action='write', position='append')
          else
            open(unit=iounit, file=outdebugfile, form='formatted', status='new', action='write')
          endif
          write(iounit,*) xline_x
          write(iounit,*) yline_x
          write(iounit,*) zline_x
          write(iounit,*) x_procbl_pt1, x_procbl_pt2
          write(iounit,*) '---------------------'
          close(iounit)
      endif
      do imb = 1, x_dom_numbl
        ! set outputs to meaningless values
        intrbl_pt1(:) = -1.0d0; intrbl_pt2(:) = -1.0d0
        indices_pt1(:) = -1   ; indices_pt2 = -1; intersection_exists = .false.
        ! calculate intersection between domain block and processor block
        call get_intersection_box(x_dombl_pt1(:,imb), x_dombl_pt2(:,imb), x_procbl_pt1, x_procbl_pt2, &
               intrbl_pt1, intrbl_pt2, xline_x, yline_x, zline_x, indices_pt1, indices_pt2, intersection_exists)
        if(detailed_debug) then
            iounit = 100+nrank
            write(outdebugfile, '(a,i6.6,a)') 'debuginfo_', nrank,'.dat'
            inquire(file=outdebugfile, exist=file_exists)
            if(file_exists) then
              open(unit=iounit, file=outdebugfile, form='formatted', status='old', action='write', position='append')
            else
              open(unit=iounit, file=outdebugfile, form='formatted', status='new', action='write')
            endif
            write(iounit,*) intersection_exists, x_num_blocks, imb
            write(iounit,*) x_dombl_pt1(:,imb), x_dombl_pt2(:,imb)
            write(iounit,*) intrbl_pt1, intrbl_pt2
            write(iounit,*) indices_pt1, indices_pt2
            write(iounit,*) '    -----------------'
            close(iounit)
        endif
        if(intersection_exists) then
           ! store important information in temporary vars
           x_num_blocks = x_num_blocks + 1
           xst(:,x_num_blocks) = indices_pt1
           xen(:,x_num_blocks) = indices_pt2
        endif
      enddo
      ! transfer information from temporary vars to topology object
      this%x_num_blocks = x_num_blocks
      allocate(this%xst(3,this%x_num_blocks))
      allocate(this%xen(3,this%x_num_blocks))
      do imb = 1, this%x_num_blocks
          this%xst(:,imb) = xst(:,imb)
          this%xen(:,imb) = xen(:,imb)
      enddo
      !! STEP 4a :: Done

      !! STEP 4b :: for y-decomposition
      do imb = 1, y_dom_numbl
        ! calculate intersection between domain block and processor block
        call get_intersection_box(y_dombl_pt1(:,imb), y_dombl_pt2(:,imb), y_procbl_pt1, y_procbl_pt2, &
               intrbl_pt1, intrbl_pt2, xline_y, yline_y, zline_y, indices_pt1, indices_pt2, intersection_exists)
        if(intersection_exists) then
           ! store important information in temporary vars
           y_num_blocks = y_num_blocks + 1
           yst(:,y_num_blocks) = indices_pt1
           yen(:,y_num_blocks) = indices_pt2

           ! store lower boundary indices
           ybclo_st(:, y_num_blocks) = (/ indices_pt1(1), indices_pt1(2), indices_pt1(3)/)
           ybclo_en(:, y_num_blocks) = (/ indices_pt2(1), indices_pt1(2), indices_pt2(3)/)  !! note 2nd entry here

           ! store upper boundary indices
           ybchi_st(:, y_num_blocks) = (/ indices_pt1(1), indices_pt2(2), indices_pt1(3)/)  !! note 2nd entry here
           ybchi_en(:, y_num_blocks) = (/ indices_pt2(1), indices_pt2(2), indices_pt2(3)/)
        endif
      enddo


      !if(nrank==1) then
      !    print *, 'Topology: bclo_st:', ybclo_st(:,1:y_num_blocks), 'bclo_en:', ybclo_en(:,1:y_num_blocks)
      !endif

      ! figure out vertical internal boundaries --  only for 2D blocks (in x-y) as of now
      y_num_intbd_left(:) = 0;      y_num_intbd_rght(:) = 0;
      y_intbd_left_st(:,:) = -1;    y_intbd_rght_st(:,:) = -1
      y_intbd_left_en(:,:) = -1;    y_intbd_rght_en(:,:) = -1
      do imb1 = 1, y_num_blocks
        ! when imb1 is to the left of imb2
        i1 = ybclo_en(1, imb1);  j1 = ybclo_en(2, imb1);  k1 = ybclo_en(3, imb1)
        do imb2 = imb1+1, y_num_blocks
            if(imb1==imb2) then
                ! internal boundary cannot exist if both blocks are the same
                cycle
            endif
            i2 = ybclo_st(1, imb2);  j2 = ybclo_st(2, imb2);  k2 = ybclo_st(3, imb2)
            if((j1/=j2) .and. (i1+1==i2)) then      !! internal boundary exists only if j1/=j2
                !! internal boundary exists
                !! count this for imb1 or imb2; is it left or right boundary
                if(j1 > j2) then
                    !! block imb2 has a left internal boundary
                    y_num_intbd_left(imb2) = y_num_intbd_left(imb2) + 1
                    i_intbd = y_num_intbd_left(imb2)
                    y_intbd_left_st(:, i_intbd) = (/i1, j2, k2/)
                    y_intbd_left_en(:, i_intbd) = (/i1, j1-1, k1/)
                else
                    !! block imb1 has a right internal boundary
                    y_num_intbd_rght(imb1) = y_num_intbd_rght(imb1) + 1
                    i_intbd = y_num_intbd_rght(imb1)
                    y_intbd_rght_st(:, i_intbd) = (/i2, j1, k2/)
                    y_intbd_rght_en(:, i_intbd) = (/i2, j2-1, k1/)
                endif
            endif
        enddo

        ! when imb1 is to the right of imb2
        i1 = ybclo_st(1, imb1);  j1 = ybclo_st(2, imb1);  k1 = ybclo_st(3, imb1)
        do imb2 = imb1+1, y_num_blocks
            if(imb1==imb2) then
                ! internal boundary cannot exist if both blocks are the same
                cycle
            endif
            i2 = ybclo_en(1, imb2);  j2 = ybclo_en(2, imb2);  k2 = ybclo_en(3, imb2)
            if((j1/=j2) .and. (i2+1==i1)) then
                !! internal boundary exists
                !! count this for imb1 or imb2; is it left or right boundary
                if(j1 > j2) then
                    !! block imb2 has a right internal boundary
                    y_num_intbd_rght(imb2) = y_num_intbd_rght(imb2) + 1
                    i_intbd = y_num_intbd_rght(imb2)
                    y_intbd_rght_st(:, i_intbd) = (/i1, j2, k1/)
                    y_intbd_rght_en(:, i_intbd) = (/i1, j1-1, k2/)
                else
                    !! block imb1 has a left internal boundary
                    y_num_intbd_left(imb1) = y_num_intbd_left(imb1) + 1
                    i_intbd = y_num_intbd_left(imb1)
                    y_intbd_left_st(:, i_intbd) = (/i2, j1, k1/)
                    y_intbd_left_en(:, i_intbd) = (/i2, j2-1, k2/)
                endif
            endif
        enddo
      enddo

      ! transfer information from temporary vars to topology object
      ! first for the bulk
      this%y_num_blocks = y_num_blocks
      allocate(this%yst(3,this%y_num_blocks))
      allocate(this%yen(3,this%y_num_blocks))
      do imb = 1, this%y_num_blocks
          this%yst(:,imb) = yst(:,imb)
          this%yen(:,imb) = yen(:,imb)
      enddo
      ! next for the boundaries (needed only for ydecomposition)
      ! bottom and top boundaries
      allocate(this%ybclo_st(3,this%y_num_blocks), this%ybclo_en(3,this%y_num_blocks))
      allocate(this%ybchi_st(3,this%y_num_blocks), this%ybchi_en(3,this%y_num_blocks))
      do imb = 1, this%y_num_blocks
          this%ybclo_st(:,imb) = ybclo_st(:,imb)
          this%ybclo_en(:,imb) = ybclo_en(:,imb)
          this%ybchi_st(:,imb) = ybchi_st(:,imb)
          this%ybchi_en(:,imb) = ybchi_en(:,imb)
      enddo
      ! internal (left/right) boundaries
      allocate(this%y_num_intbd_left(this%y_num_blocks))
      allocate(this%y_num_intbd_rght(this%y_num_blocks))
      this%y_num_intbd_left(1:this%y_num_blocks) = y_num_intbd_left(1:this%y_num_blocks)
      this%y_num_intbd_rght(1:this%y_num_blocks) = y_num_intbd_rght(1:this%y_num_blocks)
      allocate(this%y_intbd_left_st(3,this%y_num_blocks))
      allocate(this%y_intbd_left_en(3,this%y_num_blocks))
      allocate(this%y_intbd_rght_st(3,this%y_num_blocks))
      allocate(this%y_intbd_rght_en(3,this%y_num_blocks))
      do imb = 1, this%y_num_blocks
        do i_intbd = 1, this%y_num_intbd_left(imb)
          this%y_intbd_left_st(:,i_intbd) = y_intbd_left_st(:,i_intbd)
          this%y_intbd_left_en(:,i_intbd) = y_intbd_left_en(:,i_intbd)
        enddo
        do i_intbd = 1, this%y_num_intbd_rght(imb)
          this%y_intbd_rght_st(:,i_intbd) = y_intbd_rght_st(:,i_intbd)
          this%y_intbd_rght_en(:,i_intbd) = y_intbd_rght_en(:,i_intbd)
        enddo
      enddo
      !! STEP 4b :: Done

      !! STEP 4c :: for z-decomposition
      do imb = 1, z_dom_numbl
        ! calculate intersection between domain block and processor block
        call get_intersection_box(z_dombl_pt1(:,imb), z_dombl_pt2(:,imb), z_procbl_pt1, z_procbl_pt2, &
               intrbl_pt1, intrbl_pt2, xline_z, yline_z, zline_z, indices_pt1, indices_pt2, intersection_exists)
        if(intersection_exists) then
           ! store important information in temporary vars
           z_num_blocks = z_num_blocks + 1
           zst(:,z_num_blocks) = indices_pt1
           zen(:,z_num_blocks) = indices_pt2
        endif
      enddo
      ! transfer information from temporary vars to topology object
      this%z_num_blocks = z_num_blocks
      allocate(this%zst(3,this%z_num_blocks))
      allocate(this%zen(3,this%z_num_blocks))
      do imb = 1, this%z_num_blocks
          this%zst(:,imb) = zst(:,imb)
          this%zen(:,imb) = zen(:,imb)
      enddo
      !! STEP 4c :: Done
      !! STEP 4  :: Done

      !! STEP 5 :: Create Mask
      allocate(this%mask(decomp%ysz(1),  decomp%ysz(2),   decomp%ysz(3)))
      this%mask = zero
      !! construct mask in y from x-decomposition
      xtmp = zero
      do imb = 1, this%x_num_blocks
          i1 = this%xst(1,imb);    j1 = this%xst(2,imb);    k1 = this%xst(3,imb)
          i2 = this%xen(1,imb);    j2 = this%xen(2,imb);    k2 = this%xen(3,imb)
          do kk = k1, k2
            do jj = j1, j2
               xtmp(i1:i2, jj, kk) = one
            enddo
          enddo
      enddo
      call transpose_x_to_y(xtmp, this%mask, decomp)
      call decomp_2d_write_one(2, this%mask, 'multiblock_mask_y.out', decomp)
      call decomp_2d_write_one(1, xtmp,      'multiblock_mask_x.out', decomp)
      !print '(a,e19.12,1x,e19.12,1x)', 'Mask-minmax: ', p_maxval(maxval(this%mask)), p_minval(minval(this%mask))
      !! STEP 5  :: Done Create Mask

      !! STEP 6 :: Debug
      write(outdebugfile, '(a,i6.6,a)') 'multiblock_info_', nrank,'.dat'
      open(unit=iounit, file=outdebugfile, form='formatted', status='replace', action='write')
      write(iounit,'(a)') '----X-Decomposition Full Domain----'
      write(iounit, '(6(i8,1x))') decomp%xst(1:3), decomp%xen(1:3)
      write(iounit, '(6(e19.12,1x))') xline_x(1), xline_x(decomp%xsz(1)), yline_x(1), yline_x(decomp%xsz(2)), zline_x(1), zline_x(decomp%xsz(3))
      write(iounit,*)
      write(iounit,'(a)') '----X-Decomposition Block Information----'
      write(iounit, '(i6)') this%x_num_blocks
      do imb = 1, this%x_num_blocks
        write(iounit, '(6(i8,1x))') this%xst(1:3,imb), this%xen(1:3,imb)
        write(iounit, '(6(e19.12,1x))') xline_x(this%xst(1,imb)), xline_x(this%xen(1,imb)), yline_x(this%xst(2,imb)), &
                                        yline_x(this%xen(2,imb)), zline_x(this%xst(3,imb)), zline_x(this%xen(3,imb))
      enddo
      write(iounit,*)
      write(iounit,'(a)') '----Y-Decomposition Full Domain----'
      write(iounit, '(6(i8,1x))') decomp%yst(1:3), decomp%yen(1:3)
      write(iounit, '(6(e19.12,1x))') xline_y(1), xline_y(decomp%ysz(1)), yline_y(1), yline_y(decomp%ysz(2)), zline_y(1), zline_y(decomp%ysz(3))
      write(iounit,*)
      write(iounit,'(a)') '----Y-Decomposition Block Information----'
      write(iounit, '(i6)') this%y_num_blocks
      do imb = 1, this%y_num_blocks
        write(iounit, '(6(i8,1x))') this%yst(1:3,imb), this%yen(1:3,imb)
        write(iounit, '(6(e19.12,1x))') xline_y(this%yst(1,imb)), xline_y(this%yen(1,imb)), yline_y(this%yst(2,imb)), &
                                        yline_y(this%yen(2,imb)), zline_y(this%yst(3,imb)), zline_y(this%yen(3,imb))
      enddo
      write(iounit,*)
      write(iounit,'(a)') '----Z-Decomposition Full Domain----'
      write(iounit, '(6(i8,1x))') decomp%zst(1:3), decomp%zen(1:3)
      write(iounit, '(6(e19.12,1x))') xline_z(1), xline_z(decomp%zsz(1)), yline_z(1), yline_z(decomp%zsz(2)), zline_z(1), zline_z(decomp%zsz(3))
      write(iounit,*)
      write(iounit,'(a)') '----Z-Decomposition Block Information----'
      write(iounit, '(i6)') this%z_num_blocks
      do imb = 1, this%z_num_blocks
        write(iounit, '(6(i8,1x))') this%zst(1:3,imb), this%zen(1:3,imb)
        write(iounit, '(6(e19.12,1x))') xline_z(this%zst(1,imb)), xline_z(this%zen(1,imb)), yline_z(this%zst(2,imb)), &
                                        yline_z(this%zen(2,imb)), zline_z(this%zst(3,imb)), zline_z(this%zen(3,imb))
      enddo

      ! write information about ybclo and ybchi
      write(iounit,*)
      write(iounit,'(a)') '----ybclo and ybchi in each blocks----'
      do imb = 1, this%y_num_blocks
        write(iounit,'(a,i2)') 'Block number = ',imb
        write(iounit,'(a)') '----ybclo start  and ybclo end in each blocks----'
        write(iounit, '(6(i8,1x))') this%ybclo_st(1:3,imb), this%ybclo_en(1:3,imb)
        write(iounit,'(a)') '----left internal boundary start and end in each blocks if exist----'
        if (this%y_num_intbd_left(imb) /= 0 )then
          write(iounit,'(a,i2)') 'number of left int boundary = ',this%y_num_intbd_left(imb)
          do i_intbd = 1, this%y_num_intbd_left(imb)
            write(iounit,'(a,i2)') 'left int boundary = ',i_intbd
            write(iounit,'(a)') '----left internal boundary start and end in each blocks----'
            write(iounit, '(6(i8,1x))') this%y_intbd_left_st(1:3,i_intbd), this%y_intbd_left_en(1:3,i_intbd)
          enddo
        endif
        write(iounit,'(a)') '----right internal boundary start and end in each blocks if exist----'
        if (this%y_num_intbd_rght(imb) /= 0 )then
          write(iounit,'(a,i2)') 'number of right int boundary = ',this%y_num_intbd_rght(imb)
          do i_intbd = 1, this%y_num_intbd_rght(imb)
            write(iounit,'(a,i2)') 'right int boundary = ',i_intbd
            write(iounit,'(a)') '----right internal boundary start and end in each blocks----'
            write(iounit, '(6(i8,1x))') this%y_intbd_rght_st(1:3,i_intbd), this%y_intbd_rght_en(1:3,i_intbd)
          enddo
        endif
      enddo
      ! write information about internal boundaries

      close(unit=iounit)
      !! STEP 6 :: Done

      !! STEP 7 :: Nullify and destroy
      nullify(ztmp, xtmp)
      deallocate(xline_z, yline_z, zline_z)
      deallocate(xline_y, yline_y, zline_y)
      deallocate(xline_x, yline_x, zline_x)
      !! STEP 7 :: Done

  end subroutine

  subroutine destroy(this)
      class(multiblocktopol), intent(inout) :: this

      deallocate(this%mask)
      deallocate(this%zen)
      deallocate(this%zst)
      deallocate(this%yen)
      deallocate(this%yst)
      deallocate(this%y_intbd_rght_en)
      deallocate(this%y_intbd_rght_st)
      deallocate(this%y_intbd_left_en)
      deallocate(this%y_intbd_left_st)
      deallocate(this%y_num_intbd_rght)
      deallocate(this%y_num_intbd_left)
      deallocate(this%ybclo_en, this%ybclo_st)
      deallocate(this%ybchi_en, this%ybchi_st)
      deallocate(this%xen)
      deallocate(this%xst)

  end subroutine

  !subroutine get_st_en_indices()
  subroutine  get_intersection_box(domn_block_pt1, domn_block_pt2, proc_block_pt1, proc_block_pt2, intr_block_pt1, intr_block_pt2, xline, yline, zline, indices_pt1, indices_pt2, intersection_exists)
      real(rkind), dimension(3), intent(in)  ::  domn_block_pt1, domn_block_pt2  !! domain block
      real(rkind), dimension(3), intent(in)  ::  proc_block_pt1, proc_block_pt2  !! processor block
      real(rkind), dimension(3), intent(out) ::  intr_block_pt1, intr_block_pt2  !! intersecting block
      real(rkind), dimension(:), intent(in)  ::  xline, yline, zline
      integer,     dimension(3), intent(out) ::  indices_pt1, indices_pt2
      logical,                   intent(out) ::  intersection_exists

      real(rkind) :: proc_x1, proc_y1, proc_z1, proc_x2, proc_y2, proc_z2
      real(rkind) :: domn_x1, domn_y1, domn_z1, domn_x2, domn_y2, domn_z2
      real(rkind) :: intr_x1, intr_y1, intr_z1, intr_x2, intr_y2, intr_z2
      integer :: imb, num_intr_elems
  
      !!%% unpack the inputs
      proc_x1 = proc_block_pt1(1); proc_y1 = proc_block_pt1(2); proc_z1 = proc_block_pt1(3)
      proc_x2 = proc_block_pt2(1); proc_y2 = proc_block_pt2(2); proc_z2 = proc_block_pt2(3)
      
      domn_x1 = domn_block_pt1(1); domn_y1 = domn_block_pt1(2); domn_z1 = domn_block_pt1(3)
      domn_x2 = domn_block_pt2(1); domn_y2 = domn_block_pt2(2); domn_z2 = domn_block_pt2(3)
      
      intr_x1 = max(proc_x1, domn_x1); intr_x2 = min(proc_x2, domn_x2);
      intr_y1 = max(proc_y1, domn_y1); intr_y2 = min(proc_y2, domn_y2);
      intr_z1 = max(proc_z1, domn_z1); intr_z2 = min(proc_z2, domn_z2);
      
      intr_block_pt1(1) = intr_x1;  intr_block_pt1(2) = intr_y1; intr_block_pt1(3) = intr_z1
      intr_block_pt2(1) = intr_x2;  intr_block_pt2(2) = intr_y2; intr_block_pt2(3) = intr_z2
      !!intr_block = [intr_x1  intr_y1  intr_z1  intr_x2  intr_y2  intr_z2];

      call get_closest_indices(intr_block_pt1, xline, yline, zline, indices_pt1, .true. ) !! left  point
      call get_closest_indices(intr_block_pt2, xline, yline, zline, indices_pt2, .false.) !! right point
  
      num_intr_elems = (indices_pt2(1)-indices_pt1(1)+1) * (indices_pt2(2)-indices_pt1(2)+1) * (indices_pt2(3)-indices_pt1(3)+1);
      if (num_intr_elems < 1) then
          intersection_exists = .false.
      else
          intersection_exists = .true.
      endif

  end subroutine

  subroutine get_closest_indices(refpt, xline, yline, zline, indices_pt, pt_on_left)
      real(rkind), dimension(3), intent(in)  :: refpt
      real(rkind), dimension(:), intent(in)  :: xline, yline, zline
      integer,     dimension(3), intent(out) :: indices_pt
      logical,                   intent(in)  :: pt_on_left

      integer :: ii, jj, kk

      ii = minloc(abs(xline-refpt(1)), 1)
      jj = minloc(abs(yline-refpt(2)), 1)
      kk = minloc(abs(zline-refpt(3)), 1)
      !print*,ii,'ii',jj,'jj',kk,'kk',xline(ii),'xline(ii)',yline(jj),'yline(jj)',zline(kk),'zline(kk)',refpt(1),'refpt(1)',refpt(2),'refpt(2)',refpt(3),'refpt(3)'
      if(pt_on_left) then
          !! (ii,jj,kk) must be greater than refpt
          if(xline(ii) < refpt(1)) ii = ii+1
          if(yline(jj) < refpt(2)) jj = jj+1
          if(zline(kk) < refpt(3)) kk = kk+1
      else
          !! (ii,jj,kk) must be less than refpt
          if(xline(ii) > refpt(1)) ii = ii-1
          if(yline(jj) > refpt(2)) jj = jj-1
          if(zline(kk) > refpt(3)) kk = kk-1
      endif
      indices_pt(1) = ii; indices_pt(2) = jj; indices_pt(3) = kk

  end subroutine

end module
