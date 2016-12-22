module kspreprocessing
    use kind_parameters, only: rkind, clen
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral
    use gaussianstuff, only: gaussian
    use decomp_2d
    use decomp_2d_io
    use staggOpsMod, only: staggops
    use constants, only: two, pi, eight 
    use lstsqstuff, only: lstsq

    implicit none

    private
    public :: ksprep


    type :: ksprep
        private
        real(rkind), dimension(:,:,:), allocatable :: ffil,fdump, fxDnX, fxDnY, fxDnyDnY, fxDnyDnZ, fallDnZ, ztmp
        type(spectral), pointer :: spectE, spectC
        type(decomp_info), pointer :: gpE, gpC, sp_gpE
        type(decomp_info) :: gp_xDn, gp_xDnyDn, gp_allDn, gp_Dump0, gp_xDn8, gp_xDn8yDn8
        integer, dimension(:), allocatable :: planes2dumpC, planes2dumpF
        integer :: nxF, nyF, nzF
        character(len=clen) :: outputDir
        integer :: RunID, step
        complex(rkind), dimension(:,:,:), allocatable :: cbuffY, fCtmp1, fCtmp2
        real(rkind), dimension(:,:,:), allocatable :: fdumpFinal, fdump8, fdump8x1, fdump8x2, fdump8y
        real(rkind), dimension(:,:,:), allocatable :: ufil, vfil, wfil, ztmp8, fxDn8X, fxDn8Y, fxDn8yDn8Y, fxDn8yDn8Z
        type(lstsq ) :: zfil1, zfil2
        type(spectral) :: spectSmall
        type(staggops) :: OpsSmall
        real(rkind), dimension(:,:,:), allocatable :: Gspectral
        integer, dimension(:,:), pointer :: probes
        logical :: isAllocated = .false. 
        logical :: doZFilter = .true.  

        contains 
            procedure :: initFull
            procedure :: initLES2KSfilter 
            procedure :: destroy
            procedure :: LES_for_KS
            procedure :: LES_to_KS
            procedure :: applyFilterForKS
            procedure :: link_pointers
            generic :: init => initFull, initLES2KSfilter 
    end type


contains 

    subroutine link_pointers(this, uout, vout, wout) 
        class(ksprep), intent(inout), target :: this
        real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: uout, vout, wout 
        
        if (.not.this%isAllocated) then
            call GracefulExit("You cannot call LINK_POINTERS to KSPREP before initializing it",12)
        end if

        uout => this%ufil
        vout => this%vfil
        wout => this%wfil

    end subroutine

    pure subroutine DownsampleBy4_x(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: i, j,k, ifine, nxL, nxS
        nxL = size(arrIn,1)
        nxS = (nxL)/4

        do k = 1,size(arrIn,3)
            do j = 1,size(arrIn,2)
                ifine = 1
                do i = 1,nxS
                    arrOut(i,j,k) = arrIn(ifine,j,k)
                    ifine = ifine + 4
                end do 
            end do 
        end do 

    end subroutine

    pure subroutine DownsampleBy4_y(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, jfine, nyL, nyS
        nyL = size(arrIn,2)
        nyS = (nyL)/4

        do k = 1,size(arrIn,3)
            jfine = 1
            do j = 1,nyS
                arrOut(:,j,k) = arrIn(:,jfine,k)
                jfine = jfine + 4
            end do 
        end do 

    end subroutine

    pure subroutine DownsampleBy8_x(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: i, j,k, ifine, nxL, nxS
        nxL = size(arrIn,1)
        nxS = (nxL)/8

        do k = 1,size(arrIn,3)
            do j = 1,size(arrIn,2)
                ifine = 1
                do i = 1,nxS
                    arrOut(i,j,k) = arrIn(ifine,j,k)
                    ifine = ifine + 8
                end do 
            end do 
        end do 

    end subroutine

    pure subroutine DownsampleBy8_y(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, jfine, nyL, nyS
        nyL = size(arrIn,2)
        nyS = (nyL)/8

        do k = 1,size(arrIn,3)
            jfine = 1
            do j = 1,nyS
                arrOut(:,j,k) = arrIn(:,jfine,k)
                jfine = jfine + 8
            end do 
        end do 

    end subroutine

    pure subroutine DownsampleBy2_z_E2C(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: k, nzL, nzS
        nzL = size(arrIn,3)
        nzS = (nzL-1)/2

        do k = 1,nzS
            arrOut(:,:,k) = arrIn(:,:,2*k)
        end do 

    end subroutine


    subroutine initLES2KSfilter(this, spectC, gpC, dx, dy, outputdir, RunID, probes, FilFact, doZfilter)
        class(ksprep), intent(inout) :: this
        integer, intent(in) :: RunID
        type(decomp_info), target, intent(in) :: gpC 
        type(spectral), target, intent(in) :: spectC
        character(len=*), intent(in) :: outputdir
        real(rkind), intent(in) :: FilFact
        real(rkind), intent(in) :: dx, dy
        integer, dimension(:,:), intent(in), target, optional :: probes
        integer :: ierr, i, j, k 
        real(rkind) :: kdealiasx, kdealiasy
        logical, intent(in) :: doZfilter

        if (this%isAllocated) then
            call GracefulExit("You cannot allocate ksprep derived type if it has already been allocated",12)
        end if
        this%RunID = runID
        this%spectC => spectC
        this%gpC => gpC
        this%outputdir = outputdir
        this%doZfilter = doZfilter
        ierr = this%zfil1%init(gpC%zsz(3),.false.)
        if (allocated(this%cbuffY)) deallocate(this%cbuffY)
        call this%spectC%alloc_r2c_out(this%cbuffY)
        if (allocated(this%ffil)) deallocate(this%ffil)

        allocate(this%ufil(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
        allocate(this%vfil(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
        allocate(this%wfil(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))

        allocate(this%fCtmp1(spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2),spectC%spectdecomp%zsz(3)))
        allocate(this%fCtmp2(spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2),spectC%spectdecomp%zsz(3)))
       
        if (present(probes)) then
            this%probes => probes
        end if 

        allocate(this%Gspectral(size(this%cbuffY,1),size(this%cbuffY,2),size(this%cbuffY,3)))

        kdealiasx = (1.d0/FilFact)*pi/dx
        kdealiasy = (1.d0/FilFact)*pi/dy
        do k = 1,size(this%Gspectral,3)
            do j = 1,size(this%Gspectral,2)
                do i = 1,size(this%Gspectral,1)
                    if ((abs(this%spectC%k1(i,j,k)) < kdealiasx) .and. (abs(this%spectC%k2(i,j,k))< kdealiasy)) then
                        this%Gspectral(i,j,k) = 1.d0
                    else
                        this%Gspectral(i,j,k) = 0.d0
                    end if
                end do 
            end do 
        end do 


        this%isAllocated = .true. 

    end subroutine

    subroutine applyFilterForKS(this, u, v, w)
        class(ksprep), intent(inout) :: this
        real(rkind), intent(in), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: u, v, w

        ! Filter u into ufil 
        call this%spectC%fft(u,this%cbuffY)
        this%cbuffY = this%cbuffY*this%Gspectral
        if (this%doZfilter) then
            call transpose_y_to_z(this%cbuffY, this%fCtmp1,this%spectC%spectdecomp)
            call this%zfil1%filter3(this%fCtmp1,this%fCtmp2,size(this%fCtmp1,1),size(this%fCtmp1,2))
            call transpose_z_to_y(this%fCtmp2,this%cbuffY,this%spectC%spectdecomp)
        end if
        call this%spectC%ifft(this%cbuffY,this%ufil)

        
        ! Filter v into vfil 
        call this%spectC%fft(v,this%cbuffY)
        this%cbuffY = this%cbuffY*this%Gspectral
        if (this%doZfilter) then
            call transpose_y_to_z(this%cbuffY, this%fCtmp1,this%spectC%spectdecomp)
            call this%zfil1%filter3(this%fCtmp1,this%fCtmp2,size(this%fCtmp1,1),size(this%fCtmp1,2))
            call transpose_z_to_y(this%fCtmp2,this%cbuffY,this%spectC%spectdecomp)
        end if
        call this%spectC%ifft(this%cbuffY,this%vfil)


        ! Filter w into wfil 
        call this%spectC%fft(w,this%cbuffY)
        this%cbuffY = this%cbuffY*this%Gspectral
        if (this%doZfilter) then
            call transpose_y_to_z(this%cbuffY, this%fCtmp1,this%spectC%spectdecomp)
            call this%zfil1%filter3(this%fCtmp1,this%fCtmp2,size(this%fCtmp1,1),size(this%fCtmp1,2))
            call transpose_z_to_y(this%fCtmp2,this%cbuffY,this%spectC%spectdecomp)
        end if
        call this%spectC%ifft(this%cbuffY,this%wfil)

    end subroutine

    subroutine dumpKSfilteredFields(this)
        use decomp_2d_io
        class(ksprep), intent(in) :: this
        character(len=clen) :: tempname, fname
        character(len=4) :: label
        
        label = "u_ks"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%ufil,fname, this%gpC)

        label = "v_ks"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%vfil,fname, this%gpC)

        label = "w_ks"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%wfil,fname, this%gpC)


    end subroutine



    subroutine initFull(this, nx, ny, nz, spectE, gpE, outputdir, RunID, dx, dy, dz, planes2DumpC, planes2DumpF)
        class(ksprep), intent(inout) :: this
        real(rkind), intent(in) :: dx, dy, dz
        integer, intent(in) :: runID, nx, ny, nz
        type(decomp_info), target, intent(in) :: gpE 
        type(spectral), target, intent(in) :: spectE
        integer, dimension(:), intent(in) :: planes2dumpC, planes2dumpF
        character(len=*), intent(in) :: outputdir
        integer :: ierr

        if (this%isAllocated) then
            call GracefulExit("You cannot allocate ksprep derived type if it has already been allocated",12)
        end if
        this%outputdir = outputdir
        this%gpE => gpE
        this%sp_gpE => spectE%spectdecomp
        this%spectE => spectE
        this%RunID = runID

        call decomp_info_init(nx/4,ny  ,nz+1,this%gp_xDn   )
        call decomp_info_init(nx/4,ny/4,nz+1,this%gp_xDnyDn)
        call decomp_info_init(nx/4,ny/4,nz/2,this%gp_allDn )

        allocate(this%ffil(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
        call this%spectE%alloc_r2c_out(this%cbuffY)
        allocate(this%fxDnX(this%gp_xDn%xsz(1),this%gp_xDn%xsz(2),this%gp_xDn%xsz(3)))
        allocate(this%fxDnY(this%gp_xDn%ysz(1),this%gp_xDn%ysz(2),this%gp_xDn%ysz(3)))
        allocate(this%fxDnyDnY(this%gp_xDnyDn%ysz(1),this%gp_xDnyDn%ysz(2),this%gp_xDnyDn%ysz(3)))
        allocate(this%fxDnyDnZ(this%gp_xDnyDn%zsz(1),this%gp_xDnyDn%zsz(2),this%gp_xDnyDn%zsz(3)))
        allocate(this%ztmp(this%gp_xDnyDn%zsz(1),this%gp_xDnyDn%zsz(2),this%gp_xDnyDn%zsz(3)))
        allocate(this%fallDnZ(this%gp_allDn%zsz(1),this%gp_allDn%zsz(2),this%gp_allDn%zsz(3)))
        allocate(this%fdump(this%gp_allDn%zsz(1),this%gp_allDn%zsz(2),this%gp_allDn%zsz(3)))

        ierr = this%zfil1%init(nz+1,.false.)
        ierr = this%zfil2%init(nz/2,.false.)

        allocate(this%planes2DumpC(size(planes2dumpC)))
        allocate(this%planes2DumpF(size(planes2dumpF)))
        this%planes2dumpC = planes2dumpC
        this%planes2dumpF = planes2dumpF


        call decomp_info_init(nx/8, ny/8, nz/2, this%gp_Dump0)
        call decomp_info_init(nx/8, ny, nz+1, this%gp_xDn8)
        call decomp_info_init(nx/8, ny/8, nz+1, this%gp_xDn8yDn8)


        call this%spectSmall%init("x", nx/8, ny/8, nz/2, eight*dx, eight*dy, two*dz, "four", "2/3rd", 2, .false.)
        call this%OpsSmall%init( this%gp_Dump0, this%gp_Dump0, 2 , dx*eight, dy*eight, dz*two, &
                           this%spectSmall%spectdecomp, this%spectSmall%spectdecomp, .true., .true.)
    
        allocate(this%fdump8(this%gp_Dump0%zsz(1),this%gp_Dump0%zsz(2),this%gp_Dump0%zsz(3)))               
        allocate(this%fdumpFinal(this%gp_Dump0%zsz(1),this%gp_Dump0%zsz(2),this%gp_Dump0%zsz(3)))               
        allocate(this%fdump8x1(this%gp_Dump0%xsz(1),this%gp_Dump0%xsz(2),this%gp_Dump0%xsz(3)))               
        allocate(this%fdump8x2(this%gp_Dump0%xsz(1),this%gp_Dump0%xsz(2),this%gp_Dump0%xsz(3)))               
        allocate(this%fdump8y(this%gp_Dump0%ysz(1),this%gp_Dump0%ysz(2),this%gp_Dump0%ysz(3)))               
        call this%spectSmall%alloc_r2c_out(this%fCtmp1)
        call this%spectSmall%alloc_r2c_out(this%fCtmp2)

        allocate(this%ztmp8(this%gp_Dump0%zsz(1),this%gp_Dump0%zsz(2),this%gp_Dump0%zsz(3)))
        allocate(this%fxDn8X(this%gp_xDn8%xsz(1),this%gp_xDn8%xsz(2),this%gp_xDn8%xsz(3)))
        allocate(this%fxDn8Y(this%gp_xDn8%ysz(1),this%gp_xDn8%ysz(2),this%gp_xDn8%ysz(3)))
        allocate(this%fxDn8yDn8Y(this%gp_xDn8yDn8%ysz(1),this%gp_xDn8yDn8%ysz(2),this%gp_xDn8yDn8%ysz(3)))
        allocate(this%fxDn8yDn8Z(this%gp_xDn8yDn8%zsz(1),this%gp_xDn8yDn8%zsz(2),this%gp_xDn8yDn8%zsz(3)))
        this%isAllocated = .true.  

    end subroutine

    subroutine LES_to_KS(this, uE, vE, wE, tid)
        class(ksprep), intent(inout) :: this
        real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
        integer :: idx, pid, dirid = 3
        integer, intent(in) :: tid
        character(len=clen) :: fname, tempname
        character(len=4) :: flabel


        ! u field
        call this%spectE%fft(uE,this%cbuffY)
        call this%spectE%KSprepFilter2(this%cbuffY)
        call this%spectE%ifft(this%cbuffY,this%ffil)
        call DownsampleBy8_x(this%ffil,this%fxDn8X)
        call transpose_x_to_y(this%fxDn8X,this%fxDn8Y,this%gp_xDn8)
        call DownsampleBy8_y(this%fxDn8Y,this%fxDn8yDn8Y)
        call transpose_y_to_z(this%fxDn8yDn8Y,this%fxDn8yDn8Z,this%gp_xDn8yDn8)
        call this%zfil1%filter3(this%fxDn8yDn8Z, this%zTmp8, size(this%zTmp8,1),size(this%zTmp8,2))
        call DownsampleBy2_z_E2C(this%zTmp8, this%fdump8)
        call this%zfil2%filter3(this%fdump8,this%fdumpFinal,size(this%fdump8,1),size(this%fdump8,2))
        flabel = ".usw"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump8,dirid, pid,fname,this%gp_Dump0)
        end do 
        call transpose_z_to_y(this%fdumpFinal,this%fdump8y,this%gp_Dump0)
        call transpose_y_to_x(this%fdump8y,this%fdump8x1,this%gp_Dump0)
        call this%spectSmall%fft(this%fdump8x1,this%fCtmp1)
        call this%spectSmall%mtimes_ik1_oop(this%fCtmp1,this%fCtmp2)
        call this%spectSmall%ifft(this%fCtmp2,this%fdump8x1)
        flabel = ".udx"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,this%fdump8x1,dirid, pid,fname,this%gp_Dump0)
        end do 
        call this%spectSmall%mtimes_ik2_ip(this%fCtmp1)
        call this%spectSmall%ifft(this%fCtmp1,this%fdump8x1)
        flabel = ".udy"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,this%fdump8x1,dirid, pid,fname,this%gp_Dump0)
        end do 
        call this%OpsSmall%ddz_C2C(this%fdumpFinal,this%fdump8,.true.,.false.)
        flabel = ".udz"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump8,dirid, pid,fname,this%gp_Dump0)
        end do 


        ! v field
        call this%spectE%fft(vE,this%cbuffY)
        call this%spectE%KSprepFilter2(this%cbuffY)
        call this%spectE%ifft(this%cbuffY,this%ffil)
        call DownsampleBy8_x(this%ffil,this%fxDn8X)
        call transpose_x_to_y(this%fxDn8X,this%fxDn8Y,this%gp_xDn8)
        call DownsampleBy8_y(this%fxDn8Y,this%fxDn8yDn8Y)
        call transpose_y_to_z(this%fxDn8yDn8Y,this%fxDn8yDn8Z,this%gp_xDn8yDn8)
        call this%zfil1%filter3(this%fxDn8yDn8Z, this%zTmp8, size(this%zTmp8,1),size(this%zTmp8,2))
        call DownsampleBy2_z_E2C(this%zTmp8, this%fdump8)
        call this%zfil2%filter3(this%fdump8,this%fdumpFinal,size(this%fdump8,1),size(this%fdump8,2))
        flabel = ".vsw"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump8,dirid, pid,fname,this%gp_Dump0)
        end do 
        call transpose_z_to_y(this%fdumpFinal,this%fdump8y,this%gp_Dump0)
        call transpose_y_to_x(this%fdump8y,this%fdump8x1,this%gp_Dump0)
        call this%spectSmall%fft(this%fdump8x1,this%fCtmp1)
        call this%spectSmall%mtimes_ik1_oop(this%fCtmp1,this%fCtmp2)
        call this%spectSmall%ifft(this%fCtmp2,this%fdump8x1)
        flabel = ".vdx"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,this%fdump8x1,dirid, pid,fname,this%gp_Dump0)
        end do 
        call this%spectSmall%mtimes_ik2_ip(this%fCtmp1)
        call this%spectSmall%ifft(this%fCtmp1,this%fdump8x1)
        flabel = ".vdy"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,this%fdump8x1,dirid, pid,fname,this%gp_Dump0)
        end do 
        call this%OpsSmall%ddz_C2C(this%fdumpFinal,this%fdump8,.true.,.false.)
        flabel = ".vdz"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump8,dirid, pid,fname,this%gp_Dump0)
        end do 

        
        ! w field
        call this%spectE%fft(wE,this%cbuffY)
        call this%spectE%KSprepFilter2(this%cbuffY)
        call this%spectE%ifft(this%cbuffY,this%ffil)
        call DownsampleBy8_x(this%ffil,this%fxDn8X)
        call transpose_x_to_y(this%fxDn8X,this%fxDn8Y,this%gp_xDn8)
        call DownsampleBy8_y(this%fxDn8Y,this%fxDn8yDn8Y)
        call transpose_y_to_z(this%fxDn8yDn8Y,this%fxDn8yDn8Z,this%gp_xDn8yDn8)
        call this%zfil1%filter3(this%fxDn8yDn8Z, this%zTmp8, size(this%zTmp8,1),size(this%zTmp8,2))
        call DownsampleBy2_z_E2C(this%zTmp8, this%fdump8)
        call this%zfil2%filter3(this%fdump8,this%fdumpFinal,size(this%fdump8,1),size(this%fdump8,2))
        flabel = ".wsw"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump8,dirid, pid,fname,this%gp_Dump0)
        end do 
        call transpose_z_to_y(this%fdumpFinal,this%fdump8y,this%gp_Dump0)
        call transpose_y_to_x(this%fdump8y,this%fdump8x1,this%gp_Dump0)
        call this%spectSmall%fft(this%fdump8x1,this%fCtmp1)
        call this%spectSmall%mtimes_ik1_oop(this%fCtmp1,this%fCtmp2)
        call this%spectSmall%ifft(this%fCtmp2,this%fdump8x1)
        flabel = ".wdx"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,this%fdump8x1,dirid, pid,fname,this%gp_Dump0)
        end do 
        call this%spectSmall%mtimes_ik2_ip(this%fCtmp1)
        call this%spectSmall%ifft(this%fCtmp1,this%fdump8x1)
        flabel = ".wdy"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,this%fdump8x1,dirid, pid,fname,this%gp_Dump0)
        end do 
        call this%OpsSmall%ddz_C2C(this%fdumpFinal,this%fdump8,.true.,.false.)
        flabel = ".wdz"
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump8,dirid, pid,fname,this%gp_Dump0)
        end do 




    end subroutine

    subroutine LES_for_KS(this,uE,vE,wE, tid)
        class(ksprep), intent(inout) :: this
        real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
        integer :: idx, pid, dirid = 3
        integer, intent(in) :: tid
        character(len=clen) :: fname, tempname
        character(len=4) :: flabel

        ! u field
        flabel = ".plu"
        call this%spectE%fft(uE,this%cbuffY)
        call this%spectE%KSprepFilter1(this%cbuffY)
        call this%spectE%ifft(this%cbuffY,this%ffil)
        call DownsampleBy4_x(this%ffil,this%fxDnX)
        call transpose_x_to_y(this%fxDnX,this%fxDnY,this%gp_xDn)
        call DownsampleBy4_y(this%fxDnY,this%fxDnyDnY)
        call transpose_y_to_z(this%fxDnyDnY,this%fxDnyDnZ,this%gp_xDnyDn)
        call this%zfil1%filter3(this%fxDnyDnZ, this%zTmp, size(this%zTmp,1),size(this%zTmp,2))
        call DownsampleBy2_z_E2C(this%zTmp, this%fallDnZ)
        call this%zfil2%filter3(this%fallDnZ,this%fdump,size(this%fdump,1),size(this%fdump,2))
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump,dirid, pid,fname,this%gp_allDn)
        end do  
        do idx = 1,size(this%planes2dumpF)
            pid = this%planes2dumpF(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_F",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,uE,dirid, pid, fname, this%gpE)
        end do  
       
        ! v field
        flabel = ".plv"
        call this%spectE%fft(vE,this%cbuffY)
        call this%spectE%KSprepFilter1(this%cbuffY)
        call this%spectE%ifft(this%cbuffY,this%ffil)
        call DownsampleBy4_x(this%ffil,this%fxDnX)
        call transpose_x_to_y(this%fxDnX,this%fxDnY,this%gp_xDn)
        call DownsampleBy4_y(this%fxDnY,this%fxDnyDnY)
        call transpose_y_to_z(this%fxDnyDnY,this%fxDnyDnZ,this%gp_xDnyDn)
        call this%zfil1%filter3(this%fxDnyDnZ, this%zTmp, size(this%zTmp,1),size(this%zTmp,2))
        call DownsampleBy2_z_E2C(this%zTmp, this%fallDnZ)
        call this%zfil2%filter3(this%fallDnZ,this%fdump,size(this%fdump,1),size(this%fdump,2))
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump,dirid, pid,fname,this%gp_allDn)
        end do  
        do idx = 1,size(this%planes2dumpF)
            pid = this%planes2dumpF(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_F",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,vE,dirid, pid, fname, this%gpE)
        end do  


        ! w field
        flabel = ".plw"
        call this%spectE%fft(wE,this%cbuffY)
        call this%spectE%KSprepFilter1(this%cbuffY)
        call this%spectE%ifft(this%cbuffY,this%ffil)
        call DownsampleBy4_x(this%ffil,this%fxDnX)
        call transpose_x_to_y(this%fxDnX,this%fxDnY,this%gp_xDn)
        call DownsampleBy4_y(this%fxDnY,this%fxDnyDnY)
        call transpose_y_to_z(this%fxDnyDnY,this%fxDnyDnZ,this%gp_xDnyDn)
        call this%zfil1%filter3(this%fxDnyDnZ, this%zTmp, size(this%zTmp,1),size(this%zTmp,2))
        call DownsampleBy2_z_E2C(this%zTmp, this%fallDnZ)
        call this%zfil2%filter3(this%fallDnZ,this%fdump,size(this%fdump,1),size(this%fdump,2))
        do idx = 1,size(this%planes2dumpC)
            pid = this%planes2dumpC(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_C",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(3,this%fdump,dirid, pid,fname,this%gp_allDn)
        end do  
        do idx = 1,size(this%planes2dumpF)
            pid = this%planes2dumpF(idx)
            write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_F",pid,flabel
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,wE,dirid, pid, fname, this%gpE)
        end do  
     

    end subroutine 



    subroutine destroy(this)
        class(ksprep), intent(inout) :: this

        deallocate(this%cbuffY, this%fdump)

    end subroutine


end module 
