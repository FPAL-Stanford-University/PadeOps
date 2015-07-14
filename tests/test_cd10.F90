program test_cd10

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use timer,           only: tic,toc
    use cd10stuff,       only: cd10
    implicit none

    integer omp_get_num_threads, omp_get_thread_num
    double precision omp_get_wtime
    integer tid,nt,zblock

    integer:: nx = 4096, ny=64, nz=512
    logical, parameter :: periodic = .TRUE.

    type( cd10 ) :: mycd10
    real(rkind), dimension(:,:,:), allocatable :: x,f,df,df_exact
    real(rkind) :: dx,t0,t1

    integer :: y1,yn,z1,zn
    integer :: i,j,k,ierr

    allocate( x(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( df_exact(nx,ny,nz) )

    dx = two*pi/real(nx,rkind)

    do i=1,nx
        x(i,:,:) = real(i-1,rkind)*two*pi/real(nx,rkind)
        f(i,:,:) = sin(4._rkind * x(i,:,:))
        df_exact(i,:,:) = 4._rkind * cos( 4._rkind * x(i,:,:))
    end do

    ierr = mycd10%init( nx, dx, periodic, 0, 0)

    !$OMP PARALLEL PRIVATE(t0,t1,tid,nt,zblock,y1,yn,z1,zn) SHARED(f,df,mycd10,nz,ny)
    t0 = omp_get_wtime()
    nt = omp_get_num_threads()
    zblock = nz/nt
    tid = omp_get_thread_num()
    y1=1;yn=ny; z1=zblock*tid+1; zn=zblock*(tid+1)
    call mycd10 % cd10der1(f,df,ny,nz,y1,yn,z1,zn)
    t1 = omp_get_wtime()
    print*, "Elapsed time: ",t1-t0
    !$OMP END PARALLEL
    
    !y1=1;yn=512; z1=1; zn=128
    !call tic()
    !call mycd10 % cd10der1(f,df,ny,nz,y1,yn,z1,zn)
    !call toc("Time to get x derivatives")

    !y1=1;yn=512; z1=129; zn=256
    !call tic()
    !call mycd10 % cd10der1(f,df,ny,nz,y1,yn,z1,zn)
    !call toc("Time to get x derivatives")

    !y1=1;yn=512; z1=257; zn=256+128
    !call tic()
    !call mycd10 % cd10der1(f,df,ny,nz,y1,yn,z1,zn)
    !call toc("Time to get x derivatives")

    !y1=1;yn=512; z1=256+128+1; zn=512
    !call tic()
    !call mycd10 % cd10der1(f,df,ny,nz,y1,yn,z1,zn)
    !call toc("Time to get x derivatives")

    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( df_exact )

end program
