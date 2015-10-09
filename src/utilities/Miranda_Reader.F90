module miranda_tools

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize,       &
                               transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, &
                               nrank, nproc
    use exits,           only: GracefulExit, message
    use gridtools,       only: alloc_buffs
    implicit none

    private

    integer, parameter :: u_index=1                               ! x-velocity Index
    integer, parameter :: v_index=2                               ! y-velocity Index
    integer, parameter :: w_index=3                               ! z-velocity Index
    integer, parameter :: rho_index=4                             ! Density Index
    integer, parameter :: e_index=5                               ! Energy Index
    integer, parameter :: p_index=6                               ! Pressure Index
    integer, parameter :: T_index=7                               ! Temp Index
    integer, parameter :: c_index=8                               ! Speed of sound Index
    integer, parameter :: mu_index=9                              ! Shear visc Index
    integer, parameter :: bulk_index=10                           ! Bulk visc Index
    integer, parameter :: ktc_index=11                            ! Thermal cond Index
    integer, parameter :: Diff_index=12                           ! Species diffusion Index
    integer, parameter :: Ys_index=13                             ! Species mass-fraction Index

    public :: miranda_reader

    type miranda_reader

        type(decomp_info) :: gp

        integer :: nx, ny, nz
        integer :: ax, ay, az
        integer :: px, py, pz, mirprocs
        integer :: nvars
        integer :: ns
        integer :: nsteps

        integer :: prow, pcol
        integer :: ppx, ppz

        real(rkind) :: dx, dy, dz

        integer, dimension(:,:),   allocatable :: procmap
        integer, dimension(:,:,:), allocatable :: invprocmap

        real(rkind), dimension(:,:,:,:), allocatable :: mesh
        real(rkind), dimension(:,:,:,:), allocatable :: fields

        character(len=clen) :: jobdir

        real(rkind), dimension(:,:,:),   pointer :: x 
        real(rkind), dimension(:,:,:),   pointer :: y 
        real(rkind), dimension(:,:,:),   pointer :: z 
        
        real(rkind), dimension(:,:,:),   pointer :: u                                 ! x-velocity Index
        real(rkind), dimension(:,:,:),   pointer :: v                                 ! y-velocity Index
        real(rkind), dimension(:,:,:),   pointer :: w                                 ! z-velocity Index
        real(rkind), dimension(:,:,:),   pointer :: rho                               ! Density Index
        real(rkind), dimension(:,:,:),   pointer :: e                                 ! Energy Index
        real(rkind), dimension(:,:,:),   pointer :: p                                 ! Pressure Index
        real(rkind), dimension(:,:,:),   pointer :: T                                 ! Temp Index
        real(rkind), dimension(:,:,:),   pointer :: c                                 ! Speed of sound Index
        real(rkind), dimension(:,:,:),   pointer :: mu                                ! Shear visc Index
        real(rkind), dimension(:,:,:),   pointer :: bulk                              ! Bulk visc Index
        real(rkind), dimension(:,:,:),   pointer :: ktc                               ! Thermal cond Index
        real(rkind), dimension(:,:,:),   pointer :: Diff                              ! Species diffusion Index
        real(rkind), dimension(:,:,:,:), pointer :: Ys                                ! Species mass-fraction Index
        
        contains

        procedure :: init
        procedure :: destroy
        
        procedure :: read_grid
        procedure :: read_data

        procedure, private :: get_procmap
        procedure, private :: read_metadata

    end type

contains

    subroutine init(this, jobdir_, prow_, pcol_)
        class(miranda_reader), intent(inout) :: this
        character(len=*), intent(in) :: jobdir_
        integer, intent(in) :: prow_, pcol_

        this%jobdir = trim(jobdir_)
        this%prow = prow_
        this%pcol = pcol_

        ! Get processor to grid mapping
        call this%get_procmap()

        ! Read in Miranda meta data
        call this%read_metadata()

    end subroutine

    subroutine destroy(this)
        class(miranda_reader), intent(inout) :: this

        ! Deallocate the procmap arrays
        if ( allocated( this%procmap    ) ) deallocate( this%procmap    )
        if ( allocated( this%invprocmap ) ) deallocate( this%invprocmap )
        
        ! Deallocate mesh
        if ( allocated( this%mesh ) ) deallocate( this%mesh )

        ! Deallocate fields
        if ( allocated( this%fields ) ) deallocate( this%fields )

    end subroutine

    subroutine get_procmap(this)
        class(miranda_reader), intent(inout) :: this
        
        integer :: p, xp, yp, zp, i
        character(len=clen) :: procmapfile
        integer :: ioUnit = 17

        ! Read the processor map
        WRITE(procmapfile,'(2A)') TRIM(this%jobdir),'/procmap'
        OPEN(UNIT=ioUnit, FILE=TRIM(procmapfile), FORM='FORMATTED')
        
        ! Read in number of processors per direction in Miranda run
        READ(ioUnit,*) this%pz, this%py, this%px

        ! Skip header line
        READ(ioUnit,*)

        this%mirprocs = this%px * this%py * this%pz

        ! Allocate processor rank to cartesian processor coordinate mapping
        if ( allocated(this%procmap) ) deallocate(this%procmap)
        allocate( this%procmap(this%mirprocs,4) )
        
        ! Allocate inverse cartesian processor coordinate to processor rank mapping
        if ( allocated(this%invprocmap) ) deallocate(this%invprocmap)
        allocate( this%invprocmap(this%px,this%py,this%pz) )
        
        DO i=1,this%mirprocs
            READ(ioUnit,*) p,zp,yp,xp
            this%procmap(p+1,1) = p
            this%procmap(p+1,2) = xp
            this%procmap(p+1,3) = yp
            this%procmap(p+1,4) = zp
            this%invprocmap(xp+1,yp+1,zp+1) = p
        END DO
        CLOSE(ioUnit)

    end subroutine

    subroutine read_metadata(this)
        class(miranda_reader), intent(inout) :: this

        real(kind=4) :: dx_, dy_, dz_
        integer :: i, ioUnit = 18
        character(len=clen) :: plotmir, dumchar

        ! Read plot.mir file to get # of variables and materials
        WRITE(plotmir,'(2A)') TRIM(this%jobdir),'/plot.mir'
        OPEN(UNIT=ioUnit, FILE=TRIM(plotmir), FORM='FORMATTED')

        ! Read metadata
        DO i=1,5; READ(ioUnit,*); END DO                        ! Skip first 5 lines
        READ(ioUnit,*) dumchar,this%nx,this%ny,this%nz          ! Domain size
        DO i=1,2; READ(ioUnit,*); END DO                        ! Skip next 2 lines
        READ(ioUnit,*) dumchar,dx_,dy_,dz_                      ! Grid spacing
        READ(ioUnit,*) dumchar, this%nvars                      ! # of variables
        DO i=1,this%nvars; READ(ioUnit,*); END DO               ! Skip variable lines
        READ(ioUnit,*) dumchar, this%ns                         ! # of species
        DO i=1,this%ns; READ(ioUnit,*); END DO                  ! Skip material lines
        READ(ioUnit,*) dumchar, this%nsteps                     ! # of timesteps

        this%nvars = this%nvars + 2 ! Since one of the vars is velocity

        ! Convert from kind 4 to rkind
        this%dx = real(dx_,rkind); this%dy = real(dy_,rkind); this%dz = real(dz_,rkind)

        this%ppx = this%px / this%prow
        this%ppz = this%pz / this%pcol
        if (this%ppx*this%prow .ne. this%px) then
            call GracefulExit("Incompatible number of X processors: prow does not divide Miranda's px",991)
        end if
        if (this%ppz*this%pcol .ne. this%pz) then
            call GracefulExit("Incompatible number of Z processors: pcol does not divide Miranda's pz",992)
        end if

        ! Get Miranda grid per processor
        this%ax = this%nx / this%px
        this%ay = this%ny / this%py
        this%az = this%nz / this%pz

        if (this%ax*this%px .ne. this%nx) then
            call GracefulExit("Total X grid points not divisible by no. of X processors in Miranda metadata",1001)
        end if
        if (this%ay*this%py .ne. this%ny) then
            call GracefulExit("Total Y grid points not divisible by no. of Y processors in Miranda metadata",1002)
        end if
        if (this%az*this%pz .ne. this%nz) then
            call GracefulExit("Total Z grid points not divisible by no. of Z processors in Miranda metadata",1003)
        end if

        CLOSE(ioUnit) 
        
        ! Initialize decomp
        call decomp_2d_init(this%nx, this%ny, this%nz, this%prow, this%pcol)
        call get_decomp_info(this%gp)

    end subroutine

    subroutine read_grid(this)
        class(miranda_reader), target, intent(inout) :: this
        
        integer :: xp,yp,zp,proc
        integer :: xp1,xpn,yp1,ypn,zp1,zpn
        real(kind=4), dimension(:,:,:), allocatable :: procgrid
        character(len=clen) :: procfile
        integer :: pUnit = 27

        allocate( procgrid(this%ax,this%ay,this%az) )
        
        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        call alloc_buffs(this%mesh,3,'y',this%gp)

        xp1 = (this%gp%yst(1)-1) / this%ax + 1
        xpn = this%gp%yen(1) / this%ax
        
        yp1 = 1
        ypn = this%py
        
        zp1 = (this%gp%yst(3)-1) / this%az + 1
        zpn = this%gp%yen(3) / this%az


        do zp = zp1,zpn
            do yp = yp1,ypn
                do xp = xp1,xpn
                    
                    ! Get processor ID
                    proc = this%invprocmap(xp,yp,zp)

                    ! Read in the grid
                    WRITE(procfile,'(2A,I6.6)') TRIM(this%jobdir),'/grid/p',proc
                    OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')

                    ! Read x coordinate
                    READ(pUnit) procgrid
                    this%mesh( (xp-xp1)*this%ax+1:(xp-xp1)*this%ax+this%ax, (yp-yp1)*this%ay+1:(yp-yp1)*this%ay+this%ay, (zp-zp1)*this%az+1:(zp-zp1)*this%az+this%az, 1) = real(procgrid,rkind)

                    ! Read y coordinate
                    READ(pUnit) procgrid
                    this%mesh( (xp-xp1)*this%ax+1:(xp-xp1)*this%ax+this%ax, (yp-yp1)*this%ay+1:(yp-yp1)*this%ay+this%ay, (zp-zp1)*this%az+1:(zp-zp1)*this%az+this%az, 2) = real(procgrid,rkind)
                    
                    ! Read z coordinate
                    READ(pUnit) procgrid
                    this%mesh( (xp-xp1)*this%ax+1:(xp-xp1)*this%ax+this%ax, (yp-yp1)*this%ay+1:(yp-yp1)*this%ay+this%ay, (zp-zp1)*this%az+1:(zp-zp1)*this%az+this%az, 3) = real(procgrid,rkind)

                    CLOSE(pUnit)
                
                end do
            end do
        end do
        
        deallocate( procgrid )

        ! Associate pointers for ease of use
        this%x    => this%mesh(:,:,:,1)
        this%y    => this%mesh(:,:,:,2)
        this%z    => this%mesh(:,:,:,3)

    end subroutine

    subroutine read_data(this, step)
        class(miranda_reader), target, intent(inout) :: this
        integer, intent(in) :: step

        integer :: xp,yp,zp,proc
        integer :: xp1,xpn,yp1,ypn,zp1,zpn
        real(kind=4), dimension(:,:,:), allocatable :: procdata
        character(len=clen) :: vizdir
        character(len=clen) :: procfile
        integer :: idx, pUnit = 27

        allocate( procdata(this%ax,this%ay,this%az) )
        

        if (.not. allocated(this%fields)) then
            ! Allocate fields if not already allocated
            call alloc_buffs(this%fields,this%nvars+this%ns,'y',this%gp)
        else
            ! Check size for compatibility
            if (size(this%fields,1) .ne. this%gp%ysz(1)) then
                call GracefulExit("Incorrect size of fields in read_data.",993)
            end if
            if (size(this%fields,2) .ne. this%gp%ysz(2)) then
                call GracefulExit("Incorrect size of fields in read_data.",993)
            end if
            if (size(this%fields,3) .ne. this%gp%ysz(3)) then
                call GracefulExit("Incorrect size of fields in read_data.",993)
            end if
            if (size(this%fields,4) .ne. this%nvars+this%ns) then
                call GracefulExit("Incorrect size of fields in read_data.",993)
            end if
        end if

        WRITE(vizdir,'(2A,I4.4)') TRIM(this%jobdir),'/vis',step

        xp1 = (this%gp%yst(1)-1) / this%ax + 1
        xpn = this%gp%yen(1) / this%ax
        
        yp1 = 1
        ypn = this%py
        
        zp1 = (this%gp%yst(3)-1) / this%az + 1
        zpn = this%gp%yen(3) / this%az


        do zp = zp1,zpn
            do yp = yp1,ypn
                do xp = xp1,xpn
                    
                    ! Get processor ID
                    proc = this%invprocmap(xp,yp,zp)

                    ! Read in the grid
                    WRITE(procfile,'(2A,I6.6)') TRIM(vizdir),'/p',proc
                    OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')

                    do idx = 1,this%nvars+this%ns
                        READ(pUnit) procdata
                        this%fields( (xp-xp1)*this%ax+1:(xp-xp1)*this%ax+this%ax, &
                                     (yp-yp1)*this%ay+1:(yp-yp1)*this%ay+this%ay, &
                                     (zp-zp1)*this%az+1:(zp-zp1)*this%az+this%az, idx ) = real(procdata,rkind)
                    end do

                    CLOSE(pUnit)
                
                end do
            end do
        end do

        deallocate( procdata )
        
        ! Associate pointers for ease of use
        this%u     => this%fields(:,:,:,   u_index)
        this%v     => this%fields(:,:,:,   v_index)
        this%w     => this%fields(:,:,:,   w_index)
        this%rho   => this%fields(:,:,:, rho_index)
        this%e     => this%fields(:,:,:,   e_index)
        this%p     => this%fields(:,:,:,   p_index)
        this%T     => this%fields(:,:,:,   T_index)
        this%c     => this%fields(:,:,:,   c_index)
        this%mu    => this%fields(:,:,:,  mu_index)
        this%bulk  => this%fields(:,:,:,bulk_index)
        this%ktc   => this%fields(:,:,:, ktc_index)
        this%Diff  => this%fields(:,:,:,Diff_index)
        this%Ys    => this%fields(:,:,:,  Ys_index:Ys_index+this%ns-1)

    end subroutine

end module
