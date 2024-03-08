module IOsgridMod

    use mpi
    use kind_parameters, only : rkind,clen
    use decomp_2d,       only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                          transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y,&
                          update_halo, nrank, nproc
    use Lib_VTK_IO
    use IR_Precision
    use exits,           only: GracefulExit, message
    use SolidMixtureMod, only: solid_mixture
    implicit none

    type :: IOsgrid

        integer :: nprimary
        integer :: vizcount
        character(len=clen) :: file_prefix
        character(len=clen) :: vizdir
        character(len=clen), dimension(:), allocatable :: primary_names

    contains

        procedure :: init
        procedure :: destroy

        procedure :: WriteViz
        procedure :: SetVizcount
    
    end type

contains

    subroutine init(this, vizdir_, file_prefix_, nprimary_, primary_names_)
        class(IOsgrid),                         intent(inout) :: this
        character(len=*),                       intent(in)    :: vizdir_
        character(len=*),                       intent(in)    :: file_prefix_
        integer,                                intent(in)    :: nprimary_
        character(len=*), dimension(nprimary_), intent(in)    :: primary_names_

        integer :: i

        this%vizcount = 0
        this%vizdir = vizdir_

        ! Create vizdir if it does not exist
        ! call execute_command_line('mkdir -p ' // adjustl(trim(this%vizdir)))
        call system('mkdir -p ' // adjustl(trim(this%vizdir)))

        this%file_prefix = ''
        if (trim(file_prefix_) .NE. '') then
            this%file_prefix = trim(file_prefix_) // '_'
        end if
        
        this%nprimary = nprimary_

        if(size(primary_names_,1) .ne. this%nprimary) then
            call GracefulExit("Incompatible number of variables and number of variable names in IOsgrid",981)
        end if

        if (allocated(this%primary_names)) deallocate(this%primary_names)
        allocate(this%primary_names(this%nprimary))

        do i=1,this%nprimary
            this%primary_names(i) = trim(primary_names_(i))
        end do

    end subroutine

    subroutine destroy(this)
        class(IOsgrid), intent(inout) :: this

        this%vizcount = 0
        this%vizdir = ''
        this%file_prefix = ''
        this%nprimary = 0
        if (allocated(this%primary_names)) deallocate(this%primary_names)

    end subroutine

    subroutine WriteViz(this, gp, mesh, primary, mix, tsim)
        class(IOsgrid),                                                      intent(inout) :: this
        class(decomp_info),                                                  intent(in)    :: gp
        real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3),3),             intent(in)    :: mesh
        real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3),this%nprimary), intent(in)    :: primary
        class(solid_mixture),                                                intent(in)    :: mix
        real(rkind), optional,                                               intent(in)    :: tsim

        real(rkind), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3
        integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn
        integer :: nx,ny,nz
        integer, dimension(MPI_STATUS_SIZE) :: mpistatus
        integer :: i,ierr,E_IO,imat

        character(len=clen) :: dummy

        call system('mkdir -p ' // adjustl(trim(this%vizdir)//'/'//trim(strz(4,this%vizcount))))
        write(dummy,'(I4)') this%vizcount
        call message("Writing viz dump "//trim(dummy)//" to " //trim(this%vizdir)//'/'//trim(this%file_prefix)//trim(strz(4,this%vizcount))//'.pvts')

        nx = gp%xsz(1); ny = gp%ysz(2); nz = gp%zsz(3)

        nx1 = gp%yst(1)
        nx2 = gp%yen(1)+1
        ny1 = gp%yst(2)
        ny2 = gp%yen(2)      ! no +1 here since we're in the y decomposition
        nz1 = gp%yst(3)
        nz2 = gp%yen(3)+1

        ! No overlap point for boundary processors
        if ( (gp%yen(1) == nx) ) then
            nx2 = nx
        end if
        if ( (gp%yen(3) == nz) ) then
            nz2 = nz
        end if

        nn = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
    
        E_IO = VTK_INI_XML_WRITE(fformat='binary', &
                           filename=trim(this%vizdir)//'/'//trim(strz(4,this%vizcount))//'/'//trim(this%file_prefix)//trim(strz(4,this%vizcount))//'_'//trim(strz(6,nrank))//'.vts', &
                           mesh_topology='StructuredGrid', nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)
        
        E_IO = VTK_FLD_XML(fld_action='open')
        if (present(tsim)) then
            E_IO = VTK_FLD_XML(fld=tsim,fname='TIME')
        end if
        E_IO = VTK_FLD_XML(fld=this%vizcount,fname='CYCLE')
        E_IO = VTK_FLD_XML(fld_action='close')

        ! Halo update for x, y and z
        call update_halo(mesh(:,:,:,1),tmp1,1,gp,.FALSE.)
        call update_halo(mesh(:,:,:,2),tmp2,1,gp,.FALSE.)
        call update_halo(mesh(:,:,:,3),tmp3,1,gp,.FALSE.)

        E_IO = VTK_GEO_XML_WRITE(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2,NN=nn,&
                                 X=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1),          &
                                 Y=tmp2(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1),          &
                                 Z=tmp3(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1)           )

        if( allocated(tmp1) ) deallocate(tmp1)
        if( allocated(tmp2) ) deallocate(tmp2)
        if( allocated(tmp3) ) deallocate(tmp3)

        E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
        
        do i=1,this%nprimary
            call update_halo(primary(:,:,:,i),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname=trim(this%primary_names(i)),var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            ! E_IO = VTK_VAR_XML(NC_NN=nn,varname="vector",varX=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1),
            !                                              varY=tmp2(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1),
            !                                              varZ=tmp3(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
        end do

        do imat=1,mix%ns
            call update_halo(mix%material(imat)%g11,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g11',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g12,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g12',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g13,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g13',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g21,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g21',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g22,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g22',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g23,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g23',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g31,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g31',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g32,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g32',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%g33,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_g33',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%Ys,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_Ys',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%VF,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_VF',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%eh,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_eh',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%eel,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_eel',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%kap,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_kap',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%diff,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_diff',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%diff_g,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_diff_g',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%elastic%yield,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_yield',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%elastic%mu,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_mu',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%e_p,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_ep',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%e_pp,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_epp',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            !call update_halo((mix%material(imat)%e_pp-mix%material(imat)%e_p)/max(mix%material(imat)%e_p,1.0e-16),tmp1,1,gp,.FALSE.)
            call update_halo(mix%material(imat)%e_pp-mix%material(imat)%e_p,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_epdiff',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            !call update_halo(mix%material(imat)%pe/primary(:,:,:,1),tmp1,1,gp,.FALSE.) !pe/rho
            call update_halo(mix%material(imat)%pe,tmp1,1,gp,.FALSE.) !pe
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_pe',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%det_e,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_det_e',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%det_t,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_det_t',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%det_p,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_det_p',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%curl_e,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_curl_e',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%curl_t,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_curl_t',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%curl_p,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_curl_p',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt11,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt11',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt12,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt12',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt13,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt13',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt21,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt21',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt22,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt22',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt23,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt23',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt31,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt31',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt32,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt32',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gt33,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gt33',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%diff_gt,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_diff_gt',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp11,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp11',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp12,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp12',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp13,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp13',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp21,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp21',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp22,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp22',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp23,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp23',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp31,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp31',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp32,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp32',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%gp33,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_gp33',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%diff_gp,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_diff_gp',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%diff_pe,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_diff_pe',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
         
            call update_halo(mix%material(imat)%rhom,tmp1,1,gp,.FALSE.)
              E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_rhom',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
     
            call update_halo(mix%material(imat)%s,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_s',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
               
            call update_halo(mix%xi(:,:,:,imat),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_xi',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%VF_int(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_VFintx',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
           
            call update_halo(mix%material(imat)%VF_int(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_VFinty',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%VF_int(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_VFintz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%Ys_int(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_Ysintx',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%Ys_int(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_Ysinty',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%Ys_int(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_Ysintz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%rho_int(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_rhointx',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%rho_int(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_rhointy',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%rho_int(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_rhointz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
            
            call update_halo(mix%material(imat)%u_int,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_uint',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%v_int,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_vint',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)


            call update_halo(mix%material(imat)%w_int,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_wint',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%fluxYs,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_fluxYs',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%intSharp_aFV,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_aFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%intSharp_RFV,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_RFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%intSharp_a(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_a',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%adiff,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_aDiff',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
         
            call update_halo(mix%material(imat)%rhodiff,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_rhoDiffJ',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%YsLAD,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_YsLAD',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%vfLAD,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_vfLAD',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%material(imat)%fd,tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_fd',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

         !   call update_halo(mix%VF_fil(:,:,:,imat),tmp1,1,gp,.FALSE.)
         !   E_IO=VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_VFfil',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
         !   if( allocated(tmp1) ) deallocate(tmp1)

               
         end do

         call update_halo(mix%DerX,tmp1,1,gp,.FALSE.)
         E_IO = VTK_VAR_XML(NC_NN=nn,varname='DerX',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
         if( allocated(tmp1) ) deallocate(tmp1)

         call update_halo(mix%DerY,tmp1,1,gp,.FALSE.)
         E_IO = VTK_VAR_XML(NC_NN=nn,varname='DerY',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
         if( allocated(tmp1) ) deallocate(tmp1)
        
          call update_halo(mix%DerYstagg,tmp1,1,gp,.FALSE.)
         E_IO = VTK_VAR_XML(NC_NN=nn,varname='DerYstagg', var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
         if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intX_error,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intX_error',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intY_error,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intY_error',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%VF_intx,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='VFintX',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%VF_inty,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='VFinty',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%derX_error,tmp1,1,gp,.FALSE.)
          E_IO =VTK_VAR_XML(NC_NN=nn,varname='derX_error',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%derY_error,tmp1,1,gp,.FALSE.)
          E_IO =VTK_VAR_XML(NC_NN=nn,varname='derY_error',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%ddx_exact,tmp1,1,gp,.FALSE.)
          E_IO=VTK_VAR_XML(NC_NN=nn,varname='derX_exact',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%ddy_exact,tmp1,1,gp,.FALSE.)
          E_IO =VTK_VAR_XML(NC_NN=nn,varname='derY_exact',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%ddystagg_exact,tmp1,1,gp,.FALSE.)
          E_IO =VTK_VAR_XML(NC_NN=nn,varname='derYstagg_exact',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)
 
          call update_halo(mix%intx_exact,tmp1,1,gp,.FALSE.)
          E_IO=VTK_VAR_XML(NC_NN=nn,varname='intX_exact',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%inty_exact,tmp1,1,gp,.FALSE.)
          E_IO=VTK_VAR_XML(NC_NN=nn,varname='intY_exact',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%VF_intz,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='VFintz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%DerZ,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='DerZ',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%DivTest,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='DivTest',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%lapTest,tmp1,1,gp,.FALSE.)
          E_IO =  VTK_VAR_XML(NC_NN=nn,varname='lapTest',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%lap_error,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='lap_error',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%div_error,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='div_error',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intSharp_fFV(:,:,:,1),tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_uFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intSharp_fFV(:,:,:,2),tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_vFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intSharp_fFV(:,:,:,3),tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_wFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intSharp_hFV,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_hFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intSharp_kFV,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_kFV',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%intdiff,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_diff',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%antidiff,tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='intSharp_antidiff',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%Pmix,tmp1,1,gp,.FALSE.) 
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='Pmix',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%gradp(:,:,:,2),tmp1,1,gp,.FALSE.)
          E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradp_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%gradp(:,:,:,3),tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradp_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)

          call update_halo(mix%gradp(:,:,:,1),tmp1,1,gp,.FALSE.)
          E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradp_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
          if( allocated(tmp1) ) deallocate(tmp1)


         if((mix%use_surfaceTension)  ) then
            call update_halo(mix%kappa,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='kappa',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%maskKappa,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='maskKappa',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

           ! call update_halo(mix%phi,tmp1,1,gp,.FALSE.)
           ! E_IO = VTK_VAR_XML(NC_NN=nn,varname='phi',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
           ! if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%fmask,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='fmask',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%norm(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='norm_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
            
            call update_halo(mix%norm(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='norm_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%norm(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='norm_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%normFV(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='normFV_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%normFV(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='normFV_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%normFV(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='normFV_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradp(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradp_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradp(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradp_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradp(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradp_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradVF(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradVF_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradVF(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradVF_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradVF(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='gradVF_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)
           
            call update_halo(mix%gradVF_FV(:,:,:,1,1),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradVFFV_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradVF_FV(:,:,:,2,2),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradVFFV_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradVF_FV(:,:,:,3,3),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradVFFV_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradVF_FV(:,:,:,2,1),tmp1,1,gp,.FALSE.)
            E_IO=VTK_VAR_XML(NC_NN=nn,varname='gradVFFV_xy',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradxi(:,:,:,1),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradxi_x',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradxi(:,:,:,2),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradxi_y',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%gradxi(:,:,:,3),tmp1,1,gp,.FALSE.)
            E_IO =VTK_VAR_XML(NC_NN=nn,varname='gradxi_z',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

  
            call update_halo(mix%surfaceTension_f(:, :, :, 1),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fx',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)


            call update_halo(mix%surfaceTension_f(:, :, :, 2),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fy',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            
            call update_halo(mix%surfaceTension_f(:, :, :, 3),tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            
            call update_halo(mix%surfaceTension_e,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_e',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_pe,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_pe',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_fxx,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fxx',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_fyy,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fyy',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_fzz,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fzz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_fxz,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fxz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_fxy,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fxy',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

            call update_halo(mix%surfaceTension_fyz,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='surfaceTension_fyz',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
            if( allocated(tmp1) ) deallocate(tmp1)

       
       endif




      ! call update_halo(mix%xi,tmp1,1,gp,.FALSE.)
      ! E_IO = VTK_VAR_XML(NC_NN=nn,varname='xi',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
      ! if( allocated(tmp1) ) deallocate(tmp1)



        E_IO = VTK_DAT_XML(var_location='node',var_block_action='close')
        E_IO = VTK_GEO_XML_WRITE()
        
        E_IO = VTK_END_XML()

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (nrank == 0) then
            ! First process saves also the composite .pvts file
            E_IO = PVTK_INI_XML(filename = trim(this%vizdir)//'/'//trim(this%file_prefix)//trim(strz(4,this%vizcount))//'.pvts', mesh_topology = 'PStructuredGrid',&
                                nx1=1, nx2=nx, ny1=1, ny2=ny, nz1=1, nz2=nz, tp='Float64')
            do i=0,nproc-1
                if (i .NE. 0) then
                    call MPI_RECV(nx1,1,MPI_INTEGER,i,i        ,MPI_COMM_WORLD,mpistatus,ierr)
                    call MPI_RECV(nx2,1,MPI_INTEGER,i,i+  nproc,MPI_COMM_WORLD,mpistatus,ierr)
                    call MPI_RECV(ny1,1,MPI_INTEGER,i,i+2*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                    call MPI_RECV(ny2,1,MPI_INTEGER,i,i+3*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                    call MPI_RECV(nz1,1,MPI_INTEGER,i,i+4*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                    call MPI_RECV(nz2,1,MPI_INTEGER,i,i+5*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                end if
                E_IO = PVTK_GEO_XML(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2,&
                                    source=trim(strz(4,this%vizcount))//'/'//trim(this%file_prefix)//trim(strz(4,this%vizcount))//'_'//trim(strz(6,i))//'.vts')
            end do

            E_IO = PVTK_DAT_XML(var_location='node',var_block_action='open')
            
            do i=1,this%nprimary
                E_IO = PVTK_VAR_XML(varname=trim(this%primary_names(i)),tp='Float64')
            end do
            
            do imat=1,mix%ns
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g11', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g12', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g13', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g21', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g22', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g23', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g31', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g32', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_g33', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_Ys',  tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_VF',  tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_eh',  tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_eel', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_kap', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_diff',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_diff_g',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_yield',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_mu',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_ep',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_epp',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_epdiff',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_pe',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_det_e',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_det_t',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_det_p',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_curl_e',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_curl_t',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_curl_p',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt11', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt12', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt13', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt21', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt22', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt23', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt31', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt32', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gt33', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_diff_gt', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp11', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp12', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp13', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp21', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp22', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp23', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp31', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp32', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_gp33', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_diff_gp', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_diff_pe', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_rhom', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_s',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_xi', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_VFintx',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_VFinty',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_VFintz',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_Ysintx',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_Ysinty',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_Ysintz',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_rhointx',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_rhointy',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_rhointz',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_uint', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_vint', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_wint',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_fluxYs',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_aFV',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_RFV',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_a',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_aDiff',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_rhoDiffJ',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_YsLAD',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_vfLAD',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_fd',tp='Float64')
            end do

            E_IO = PVTK_VAR_XML(varname='DerX', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='DerY', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='DerYstagg', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intX_error', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intY_error', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='VF_intx', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='VF_inty', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='derX_error', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='derY_error', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='derX_exact', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='derY_exact', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='derYstagg_exact', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intX_exact', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intY_exact', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='VF_intz', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='DerZ', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='DivTest', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='lapTest', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='lap_error', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='div_error', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intSharp_uFV', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intSharp_vFV', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intSharp_wFV', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intSharp_hFV', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intSharp_diff', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='intSharp_antidiff', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='Pmix', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='normFV_z', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='gradp_x', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='gradp_y', tp='Float64')
            E_IO = PVTK_VAR_XML(varname='gradp_z', tp='Float64')




            if((mix%use_surfaceTension) ) then
                E_IO = PVTK_VAR_XML(varname='kappa', tp='Float64') 
                E_IO = PVTK_VAR_XML(varname='maskKappa', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='fmask', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='norm_x', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='norm_y', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='norm_z', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='normFV_x', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='normFV_y', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='normFV_z', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradp_x', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradp_y', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradp_z', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradVF_x', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradVF_y', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradVF_z', tp='Float64')                
                E_IO = PVTK_VAR_XML(varname='gradVFFV_x', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradVFFV_y', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradVFFV_z', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradVFFV_xy', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradxi_x', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradxi_y', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradxi_z', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fx', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fy', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fz', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_e', tp='Float64')       
                E_IO = PVTK_VAR_XML(varname='surfaceTension_pe', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fxx', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fyy', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fzz', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fxz', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fxy', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='surfaceTension_fyz', tp='Float64')
           endif


            E_IO = PVTK_DAT_XML(var_location='node',var_block_action='close')
            
            E_IO = PVTK_END_XML()
        else
            call MPI_SEND(nx1,1,MPI_INTEGER,0,nrank        ,MPI_COMM_WORLD,ierr)
            call MPI_SEND(nx2,1,MPI_INTEGER,0,nrank+  nproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(ny1,1,MPI_INTEGER,0,nrank+2*nproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(ny2,1,MPI_INTEGER,0,nrank+3*nproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(nz1,1,MPI_INTEGER,0,nrank+4*nproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(nz2,1,MPI_INTEGER,0,nrank+5*nproc,MPI_COMM_WORLD,ierr)
        end if

        ! Update vizcount
        this%vizcount = this%vizcount + 1

    end subroutine

    subroutine SetVizcount(this,step)
        class(IOsgrid), intent(inout) :: this
        integer, intent(in) :: step

        this%vizcount = step

    end subroutine

end module
