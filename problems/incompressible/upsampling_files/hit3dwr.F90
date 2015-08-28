module hit3dwr
    use kind_parameters, only: rkind, clen 
    use constants, only: zero 
    implicit none

contains

    subroutine write_upsampled_file(dir,u, v, w)
        character(len=*), intent(in) :: dir
        real(rkind), dimension(:,:,:), intent(in) :: u, v, w
        character(len=8),parameter                        :: uFile = "U.000000"
        character(len=8),parameter                        :: vFile = "V.000000"
        character(len=8),parameter                        :: wFile = "W.000000"
        character(len=clen) :: fname 
        integer :: Nx, Ny, Nz
        integer, dimension(3) :: sizes
        integer :: fid = 123
        character(len=clen) :: newdir, command
        integer :: ierr, system, i, j, k 

        Nx = size(u,1)
        Ny = size(u,2)
        Nz = size(u,3)

        sizes(1) = Nx
        sizes(2) = Ny
        sizes(3) = Nz

        newdir = dir(:len_trim(dir))//'/'//"upsampled"
        command = "mkdir "//trim(newdir)
        ierr = system(trim(command))

        ! U file
        fname = newdir(:len_trim(newdir))//'/'//uFile
        open(fid,file=fname,form='unformatted', access='stream')  
        write(fid) sizes(1:3)
        write(fid) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        close(fid)

        ! V file
        fname = newdir(:len_trim(newdir))//'/'//vFile
        open(fid,file=fname,form='unformatted', access='stream')  
        write(fid) sizes(1:3)
        write(fid) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        ! W file
        fname = newdir(:len_trim(newdir))//'/'//wFile
        open(fid,file=fname,form='unformatted', access='stream')  
        write(fid) sizes(1:3)
        write(fid) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)

    end subroutine

    subroutine zeroPad(fin,fout)
        complex(rkind), dimension(:,:,:), intent(in) :: fin
        complex(rkind), dimension(:,:,:), intent(out) :: fout
        integer :: Nx, Ny, Nz
        integer :: NxL, NyL, NzL

        Nx = size(fin,1)
        Ny = size(fin,2)
        Nz = size(fin,3)

        NxL = size(fout,1)
        NyL = size(fout,2)
        NzL = size(fout,3)

        fout = zero

        fout(1:Nx/2,1:Ny/2,1:Nz/2) = fin(1:Nx/2,1:Ny/2,1:Nz/2);
        fout(Nx+1:NxL,1:Ny/2,1:Nz/2) = fin(Nx/2+1:Nx,1:Ny/2,1:Nz/2);
        fout(1:Nx/2,Ny+1:NyL,1:Nz/2) = fin(1:Nx/2,Ny/2+1:Ny,1:Nz/2);
        fout(Nx+1:NxL,Ny+1:NyL,1:Nz/2) = fin(Nx/2+1:Nx,Ny/2+1:Ny,1:Nz/2);
        fout(1:Nx/2,1:Ny/2,Nz+1:NzL) = fin(1:Nx/2,1:Ny/2,Nz/2+1:Nz);
        fout(Nx+1:NxL,1:Ny/2,Nz+1:NzL) = fin(Nx/2+1:Nx,1:Ny/2,Nz/2+1:Nz);
        fout(1:Nx/2,Ny+1:NyL,Nz+1:NzL) = fin(1:Nx/2,Ny/2+1:Ny,Nz/2+1:Nz);
        fout(Nx+1:NxL,Ny+1:NyL,Nz+1:NzL) = fin(Nx/2+1:Nx,Ny/2+1:Ny,Nz/2+1:Nz);

        fout = ((real(NxL)/real(Nx))*(real(NyL)/real(Ny))*(real(NzL)/real(Nz)))*fout
    end subroutine
end module 
