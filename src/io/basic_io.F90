module basic_io
    use kind_parameters, only: rkind, clen
    implicit none 
interface write_binary
    module procedure write_2D_binary, write_3D_binary
end interface

contains

    subroutine write_2d_ascii(raw_data,filename)
        real(kind=rkind), intent(in) :: raw_data(:,:)
        character(len=*), intent(in) :: filename
        character(len=30) :: rowfmt
        integer :: n1, n2
        integer :: i, j

        n1 = size(raw_data,1)
        n2 = size(raw_data,2)
        write(rowfmt,'(A,I4,A)') '(',n1,'(es25.17E3,1x))'
        OPEN(UNIT=10, FILE=trim(filename))

        do i=1,n1
            write(10,FMT=rowfmt) ( raw_data(i,j), j=1,n2 )
        enddo
        CLOSE(10)
    end subroutine   
 
   subroutine read_2d_ascii(data2read,filename)
        implicit none 
        real(rkind), intent(out), dimension(:,:), allocatable :: data2read
        character(len=*), intent(in) :: filename
        character(len=100000) :: columncount
        integer :: nc,nr, ierr, i, j 


        open(unit=10,file=filename,access='sequential',form='formatted')
        read(10, '(a)') columncount
        rewind(10)

        nc = 0
        do i = 1,len(columncount)
            if (columncount(i:i) == 'e') nc = nc +1
        end do 
 
        nr = 0
        do 
            read(10,*, iostat = ierr)
            if (ierr == -1) exit       
            nr = nr + 1
        end do 
        rewind(10)
        
        allocate(data2read(nr,nc))

        do i=1, nr
            read (10, *)  (data2read (i, j), j=1, nc)
        end do     

        close(10)

    end subroutine
    
    subroutine write_2D_binary(data2write,filename)
        implicit none
        real(rkind), intent(in), dimension(:,:) :: data2write
        character(len=*), intent(in) :: filename

        open(unit=8,file=trim(filename),form="unformatted")
        write(8)data2write(:,:)
        close(8)

    end subroutine

    subroutine write_3D_binary(data2write,filename)
        implicit none
        real(rkind), intent(in), dimension(:,:,:) :: data2write
        character(len=*), intent(in) :: filename

        open(unit=8,file=trim(filename),form="unformatted")
        write(8)data2write(:,:,:)
        close(8)

    end subroutine


end module 
