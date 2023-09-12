module basic_io
    use kind_parameters, only: rkind, clen
    implicit none 
interface write_binary
    module procedure write_2D_binary, write_3D_binary
end interface

contains
    subroutine write_0d_ascii(raw_data,filename)
        real(kind=rkind), intent(in) :: raw_data
        character(len=*), intent(in) :: filename
        character(len=30) :: rowfmt

        write(rowfmt,'(A,I4,A)') '(',1,'(es25.17E3,1x))'
        OPEN(UNIT=10, FILE=trim(filename))

        write(10,FMT=rowfmt) raw_data
        CLOSE(10)
    end subroutine   

    subroutine write_1d_ascii(raw_data,filename)
        real(kind=rkind), intent(in) :: raw_data(:)
        character(len=*), intent(in) :: filename
        character(len=30) :: rowfmt
        integer :: n1
        integer :: i

        n1 = size(raw_data)
        write(rowfmt,'(A,I4,A)') '(',1,'(es25.17E3,1x))'
        OPEN(UNIT=10, FILE=trim(filename))

        do i=1,n1
            write(10,FMT=rowfmt) raw_data(i)
        enddo
        CLOSE(10)
    end subroutine   
 
    subroutine write_2d_ascii(raw_data,filename)
        real(kind=rkind), intent(in) :: raw_data(:,:)
        character(len=*), intent(in) :: filename
        character(len=30) :: rowfmt
        integer :: n1, n2
        integer :: i, j

        n1 = size(raw_data,1)
        n2 = size(raw_data,2)
        write(rowfmt,'(A,I4,A)') '(',n2,'(es25.17E3,1x))'
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
        character(len=1000000) :: columncount
        integer :: nc,nr, ierr, i, j 
        logical :: fileExists = .true.

        inquire(file=filename,exist=fileExists)
        call assert(fileExists,"The file does not exist -- basic_io.F90")


        open(unit=10,file=filename,access='sequential',form='formatted')
        read(10, '(a)') columncount
        rewind(10)
        
        nc = 0
        do i = 1,len(columncount)
            if (columncount(i:i) == 'E') nc = nc +1
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
   
    subroutine read_ascii_2d_2rows(data2read,filename,nrows)
        implicit none 
        real(rkind), intent(out), dimension(nrows,2) :: data2read
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nrows
        integer :: i
        logical :: fileExists = .true.

        inquire(file=filename,exist=fileExists)
        call assert(fileExists,"The file does not exist -- basic_io.F90")

        open(12, file=filename)
        do i = 1,nrows
            read(12,*) data2read(i,:)
        end do 
        close(12)



    end subroutine 
    
   subroutine read_1d_ascii(data2read,filename)
        implicit none 
        real(rkind), intent(out), dimension(:), allocatable :: data2read
        character(len=*), intent(in) :: filename
        integer :: nr, ierr, i
        logical :: fileExists = .true.

        inquire(file=filename,exist=fileExists)
        call assert(fileExists,"The file does not exist -- basic_io.F90")
        
        open(unit=10,file=filename,access='sequential',form='formatted')
        nr = 0
        do 
            read(10,*, iostat = ierr)
            if (ierr == -1) exit       
            nr = nr + 1
        end do 
        rewind(10)
       
        allocate(data2read(nr))

        do i=1, nr
            read (10, *)  data2read (i)
        end do     

        close(10)

    end subroutine
   
    subroutine read_0d_ascii(data2read,filename)
        implicit none 
        real(rkind), intent(out)  :: data2read
        character(len=*), intent(in) :: filename
        logical :: fileExists = .true.

        inquire(file=filename,exist=fileExists)
        call assert(fileExists,"The file does not exist -- basic_io.F90")

        open(unit=10,file=filename,access='sequential',form='formatted')
        read (10, *)  data2read

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

    subroutine assert(statement,mssg) 
        logical, intent(in) :: statement
        character(len=*), intent(in) :: mssg
        if (.not. statement) then
           print*, trim(mssg)
           stop
        end if 
    end subroutine

end module 
