module kspreprocessing
    use kind_parameters, only: rkind, clen
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral
    use gaussianstuff, only: gaussian
    use decomp_2d
    use decomp_2d_io

    implicit none

    private
    public :: ksprep


    type :: ksprep
        private
        real(rkind), dimension(:,:,:), allocatable :: ffil,fdump, fxDnX, fxDnY, fxDnyDnY, fxDnyDnZ, fallDnZ, ztmp
        type(spectral), pointer :: spectE
        type(decomp_info), pointer :: gpE, sp_gpE
        type(decomp_info) :: gp_xDn, gp_xDnyDn, gp_allDn
        integer, dimension(:), allocatable :: planes2dumpC, planes2dumpF
        integer :: nxF, nyF, nzF
        character(len=clen) :: outputDir
        integer :: RunID
        complex(rkind), dimension(:,:,:), allocatable :: cbuffY
        type(gaussian) :: zfil1, zfil2

        contains 
            procedure :: init
            procedure :: destroy
            procedure :: LES_for_KS
    end type


contains 

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



    subroutine init(this, nx, ny, nz, spectE, gpE, outputdir, RunID, planes2DumpC, planes2DumpF)
        class(ksprep), intent(inout) :: this
        integer, intent(in) :: runID, nx, ny, nz
        type(decomp_info), target, intent(in) :: gpE 
        type(spectral), target, intent(in) :: spectE
        integer, dimension(:), intent(in) :: planes2dumpC, planes2dumpF
        character(len=*), intent(in) :: outputdir
        integer :: ierr

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
            call decomp_2d_write_plane(1,uE,dirid, pid, fname, this%gpE)
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
            call decomp_2d_write_plane(1,uE,dirid, pid, fname, this%gpE)
        end do  
     

    end subroutine 


    subroutine destroy(this)
        class(ksprep), intent(inout) :: this

        deallocate(this%cbuffY, this%fdump)

    end subroutine


end module 
