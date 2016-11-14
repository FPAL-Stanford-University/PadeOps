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

            call update_halo(mix%material(imat)%energy,tmp1,1,gp,.FALSE.)
            E_IO = VTK_VAR_XML(NC_NN=nn,varname='mat'//trim(strz(2,imat))//'_energy',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
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
        end do

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
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_energy',  tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_eel', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_kap', tp='Float64')
                E_IO = PVTK_VAR_XML(varname='mat'//trim(strz(2,imat))//'_diff',tp='Float64')
            end do

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
