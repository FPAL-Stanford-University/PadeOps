module miranda_restart_mod

    use iso_fortran_env, only: iostat_end
    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
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

    public :: miranda_restart

    type miranda_restart

        type(decomp_info), pointer :: gp

        integer :: nx, ny, nz
        integer :: ax, ay, az
        integer :: px, py, pz, mirprocs
        integer :: nvars
        integer :: nres
        integer :: ns
        integer :: nsteps

        real(rkind) :: dx, dy, dz

        integer, dimension(:,:),   allocatable :: procmap
        integer, dimension(:,:,:), allocatable :: invprocmap

        character(len=clen) :: jobdir           ! Directory that contains all
                                                ! miranda output and restart files
        character(len=clen) :: resfile          ! Name of the miranda res file

        integer :: u_index, v_index, w_index, rho_index, e_index, Ys_index, p_index, T_index

        contains

        procedure :: init
        procedure :: destroy
        
        procedure :: read_grid
        procedure :: read_data

        procedure, private :: get_procmap
        procedure, private :: read_metadata
        procedure, private :: read_restart_times

    end type

contains

    subroutine init(this, gp_, jobdir_, resfile_)
        class(miranda_restart), intent(inout) :: this
        type(decomp_info), target, intent(in) :: gp_
        character(len=*), intent(in) :: jobdir_, resfile_

        this%jobdir  = trim(jobdir_)
        this%resfile = trim(resfile_)
        this%gp => gp_

        ! Get processor to grid mapping
        call this%get_procmap()

        ! Read in Miranda meta data
        call this%read_metadata()

        ! Check if grid size in gp matches the data set
        if (this%gp%xsz(1) /= this%nx) then
            call GracefulExit("Grid size in X in gp does not match that in the data", 347)
        end if
        if (this%gp%ysz(2) /= this%ny) then
            call GracefulExit("Grid size in Y in gp does not match that in the data", 347)
        end if
        if (this%gp%zsz(3) /= this%nz) then
            call GracefulExit("Grid size in Z in gp does not match that in the data", 347)
        end if

        u_index = 1; v_index = 2; w_index = 3; rho_index = 4; e_index = 5
        Ys_index = 6
        p_index  = Ys_index + this%ns
        T_index  = p_index + 1

    end subroutine

    subroutine destroy(this)
        class(miranda_restart), intent(inout) :: this

        ! Deallocate the procmap arrays
        if ( allocated( this%procmap    ) ) deallocate( this%procmap    )
        if ( allocated( this%invprocmap ) ) deallocate( this%invprocmap )
        
    end subroutine

    subroutine get_procmap(this)
        class(miranda_restart), intent(inout) :: this
        
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
        class(miranda_restart), intent(inout) :: this

        integer :: i, ioUnit = 18
        character(len=clen) :: plotmir, dumchar

        ! Read plot.mir file to get # of variables and materials
        WRITE(plotmir,'(2A)') TRIM(this%jobdir),'/plot.mir'
        OPEN(UNIT=ioUnit, FILE=TRIM(plotmir), FORM='FORMATTED')

        ! Read metadata
        DO i=1,5; READ(ioUnit,*); END DO                        ! Skip first 5 lines
        READ(ioUnit,*) dumchar,this%nx,this%ny,this%nz          ! Domain size
        DO i=1,2; READ(ioUnit,*); END DO                        ! Skip next 2 lines
        READ(ioUnit,*) dumchar,this%dx,this%dy,this%dz          ! Grid spacing
        READ(ioUnit,*) dumchar, this%nvars                      ! # of variables
        DO i=1,this%nvars; READ(ioUnit,*); END DO               ! Skip variable lines
        READ(ioUnit,*) dumchar, this%ns                         ! # of species
        DO i=1,this%ns; READ(ioUnit,*); END DO                  ! Skip material lines
        READ(ioUnit,*) dumchar, this%nsteps                     ! # of timesteps

        this%nvars = this%nvars + 2  ! Since one of the vars is velocity
        this%nres  = 5 + this%ns + 2 ! u, v, w, rho, e, Y's, p, T

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

    end subroutine

    subroutine read_restart_times(this, step, tsim, dt)
        class(miranda_restart), intent(inout) :: this
        integer,                intent(in)    :: step
        real(rkind),            intent(out)   :: tsim, dt

        real(rkind) :: mytsim, mytviz, mydt, mydtold, mydtolder
        character(len=clen) :: myioformat, myfileorder
        integer :: mystep, iounit = 18, reason = 0
        character(len=clen) :: charout

        ! Open the res file
        open(unit=iounit, file=trim(this%jobdir)//'/'//trim(this%resfile), form='FORMATTED', status='UNKNOWN')

        ! Find the line for current step
        do
            read(iounit, '(X,I4.4,5(3X,D22.15),2(x,a3))',iostat=reason) mystep, mytsim, mytviz, mydt, &
                                                                        mydtold, mydtolder, myioformat, myfileorder
            ! write(*, '(X,I,5D25.15,2(x,a3))') mystep, mytsim, mytviz, mydt, mydtold, mydtolder, trim(myioformat), trim(myfileorder)
            ! write(*,*) "reason = ", reason
            if  (reason == iostat_end) then ! EOF reached
                write(charout,'(A,I0.0,A,A)') "Step ", step, " not in res file ", trim(this%resfile)
                call GracefulExit(trim(charout), 349)
            else if (reason /= 0) then ! Error reading in line
                write(charout,'(A,A)') "Error reading res file ", trim(this%resfile)
                call GracefulExit(trim(charout), 349)
            else
                if (step /= mystep) then
                    continue
                else
                    tsim = mytsim
                    dt = mydt
                    exit
                end if
            end if
        end do

        close(iounit)

    end subroutine

    subroutine read_grid(this, mesh)
        class(miranda_restart), intent(inout) :: this
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),3), intent(out) :: mesh
        
        integer :: xp,yp,zp,proc
        integer, dimension(3) :: lo, hi, glo, ghi, plo, phi
        real(kind=4), dimension(:,:,:), allocatable :: procgrid
        character(len=clen) :: procfile
        integer :: pUnit = 27

        allocate( procgrid(this%ax,this%ay,this%az) )
        
        do zp = 1,this%pz
            do yp = 1,this%py
                do xp = 1,this%px
                    
                    ! Get processor ID
                    proc = this%invprocmap(xp,yp,zp)

                    ! Get lo and hi of this proc
                    lo = [ (xp-1)*this%ax+1, (yp-1)*this%ay+1, (zp-1)*this%az+1 ]
                    hi = [ xp*this%ax, yp*this%ay, zp*this%az ]

                    if ( (lo(1) > this%gp%yen(1)) .or. (lo(2) > this%gp%yen(2)) .or. (lo(3) > this%gp%yen(3)) ) then
                        cycle
                    end if
                    if ( (hi(1) < this%gp%yst(1)) .or. (hi(2) < this%gp%yst(2)) .or. (hi(3) < this%gp%yst(3)) ) then
                        cycle
                    end if

                    glo = [ 1+max(0,lo(1)-this%gp%yst(1)), 1+max(0,lo(2)-this%gp%yst(2)), 1+max(0,lo(3)-this%gp%yst(3)) ]
                    ghi = [ min(this%gp%ysz(1), hi(1)-this%gp%yst(1)+1), &
                            min(this%gp%ysz(2), hi(2)-this%gp%yst(2)+1), &
                            min(this%gp%ysz(3), hi(3)-this%gp%yst(3)+1) ]

                    plo = [ 1+max(0, this%gp%yst(1)-lo(1)), 1+max(0, this%gp%yst(2)-lo(2)), 1+max(0, this%gp%yst(3)-lo(3)) ]
                    phi = [ min( hi(1)-lo(1)+1, this%gp%yen(1)-lo(1)+1 ), &
                            min( hi(2)-lo(2)+1, this%gp%yen(2)-lo(2)+1 ), &
                            min( hi(3)-lo(3)+1, this%gp%yen(3)-lo(3)+1 ) ]

                    ! Read in the grid
                    WRITE(procfile,'(2A,I6.6)') TRIM(this%jobdir),'/grid/p',proc
                    OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')

                    ! Read x coordinate
                    READ(pUnit) procgrid
                    mesh( glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 1) = &
                        real( procgrid(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3)), rkind )

                    ! Read y coordinate
                    READ(pUnit) procgrid
                    mesh( glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 2) = &
                        real( procgrid(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3)), rkind )
                    
                    ! Read z coordinate
                    READ(pUnit) procgrid
                    mesh( glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3) = &
                        real( procgrid(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3)), rkind )

                    CLOSE(pUnit)
                
                end do
            end do
        end do
        
        deallocate( procgrid )

    end subroutine

    subroutine read_data(this, step, resdata, tsim, dt)
        class(miranda_restart), intent(inout) :: this
        integer, intent(in) :: step
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),this%nres), intent(out) :: resdata
        real(rkind), intent(out) :: tsim, dt

        integer :: xp,yp,zp,proc
        integer, dimension(3) :: lo, hi, glo, ghi, plo, phi
        double precision, dimension(:,:,:,:), allocatable :: procdata
        character(len=clen) :: resdir
        character(len=clen) :: procfile
        integer :: pUnit = 27

        call this%read_restart_times(step, tsim, dt)

        allocate( procdata(this%ax,this%ay,this%az,this%nres) )
        

        ! Check size for compatibility
        if (size(resdata,1) .ne. this%gp%ysz(1)) then
            call GracefulExit("Incorrect size of array resdata in read_data.",993)
        end if
        if (size(resdata,2) .ne. this%gp%ysz(2)) then
            call GracefulExit("Incorrect size of array resdata in read_data.",993)
        end if
        if (size(resdata,3) .ne. this%gp%ysz(3)) then
            call GracefulExit("Incorrect size of array resdata in read_data.",993)
        end if
        if (size(resdata,4) .ne. this%nres) then
            call GracefulExit("Incorrect size of array resdata in read_data.",993)
        end if

        WRITE(resdir,'(2A,I4.4)') TRIM(this%jobdir),'/res',step

        call message("Reading Miranda restart files from "//trim(resdir))

        do zp = 1,this%pz
            do yp = 1,this%py
                do xp = 1,this%px
                    
                    ! Get processor ID
                    proc = this%invprocmap(xp,yp,zp)

                    ! Get lo and hi of this proc
                    lo = [ (xp-1)*this%ax+1, (yp-1)*this%ay+1, (zp-1)*this%az+1 ]
                    hi = [ xp*this%ax, yp*this%ay, zp*this%az ]

                    if ( (lo(1) > this%gp%yen(1)) .or. (lo(2) > this%gp%yen(2)) .or. (lo(3) > this%gp%yen(3)) ) then
                        cycle
                    end if
                    if ( (hi(1) < this%gp%yst(1)) .or. (hi(2) < this%gp%yst(2)) .or. (hi(3) < this%gp%yst(3)) ) then
                        cycle
                    end if

                    glo = [ 1+max(0,lo(1)-this%gp%yst(1)), 1+max(0,lo(2)-this%gp%yst(2)), 1+max(0,lo(3)-this%gp%yst(3)) ]
                    ghi = [ min(this%gp%ysz(1), hi(1)-this%gp%yst(1)+1), &
                            min(this%gp%ysz(2), hi(2)-this%gp%yst(2)+1), &
                            min(this%gp%ysz(3), hi(3)-this%gp%yst(3)+1) ]

                    plo = [ 1+max(0, this%gp%yst(1)-lo(1)), 1+max(0, this%gp%yst(2)-lo(2)), 1+max(0, this%gp%yst(3)-lo(3)) ]
                    phi = [ min( hi(1)-lo(1)+1, this%gp%yen(1)-lo(1)+1 ), &
                            min( hi(2)-lo(2)+1, this%gp%yen(2)-lo(2)+1 ), &
                            min( hi(3)-lo(3)+1, this%gp%yen(3)-lo(3)+1 ) ]

                    ! Read in the grid
                    WRITE(procfile,'(2A,I6.6)') TRIM(resdir),'/p',proc
                    OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')

                    READ(pUnit) procdata
                    resdata( glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), :) = &
                        real( procdata(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), :), rkind )

                    CLOSE(pUnit)
                
                end do
            end do
        end do

        deallocate( procdata )
        
    end subroutine

end module
