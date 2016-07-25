program test_dct1d
    use kind_parameters, only: rkind
    use dct_1d_stuff, only: dct1d
    use timer, only: tic, toc
    implicit none
    real(rkind), dimension(:,:,:), allocatable :: input, output, inputn
    complex(rkind), dimension(:,:,:), allocatable :: cinput, coutput, cinputn
    logical :: useOrthogonalNorm = .false. 
    integer(kind=8) :: c2cplanfwd, c2cplanbwd
    real(rkind) :: normfact
    type(dct1d) :: mydct
    integer :: nx = 1024, ny = 32, nz = 32
    integer :: i 
    include "fftw3.f"
    
    allocate(input(nx, ny, nz), output(nx, ny, nz), inputn(nx, ny, nz))

    do i = 1,nx
        input(i,:,:) = real(i,rkind)
    end do 
    
    call mydct%init(nx,ny,nz, useOrthogonalNorm)

    print*, "TEST 1: Out-of-place DCTs" 
    print*, "--------------------------"
    call tic() 
    call mydct%dctx(input,output)
    call toc()
    call tic()
    call mydct%idctx(output,inputn)
    call toc()
    print*, "Error:", maxval(abs(input - inputn))
    print*, "=========================="

    inputn = input
    print*, "TEST 2: In-place DCTs" 
    print*, "--------------------------"
    call tic()
    call mydct%dctx(inputn)
    call toc()
    call tic()
    call mydct%idctx(inputn)
    call toc()
    print*, "Error:", maxval(abs(input - inputn))
    print*, "=========================="

    call mydct%destroy()

    deallocate(input, output, inputn)

    allocate(cinput(2*nx, ny, nz), coutput(2*nx,ny,nz), cinputn(2*nx,ny,nz))
    
    do i = 1,2*nx
        cinput(i,:,:) = real(i,rkind)
    end do 
    normfact = 1.d0/(2.d0*real(nx,rkind))
    call dfftw_plan_many_dft(c2cplanfwd, 1, 2*nx, ny*nz, cinput,2*nx, 1, 2*nx, &
                            coutput, 2*nx, 1, 2*nx, FFTW_FORWARD, FFTW_EXHAUSTIVE)
    call dfftw_plan_many_dft(c2cplanbwd, 1, 2*nx, ny*nz, cinput,2*nx, 1, 2*nx, &
                            coutput, 2*nx, 1, 2*nx, FFTW_BACKWARD, FFTW_EXHAUSTIVE)

    print*, "TEST 3: Out-of-place FFT with 2x extension" 
    print*, "--------------------------"
    call tic()
    call dfftw_execute_dft( c2cplanfwd, cinput, coutput)  
    call toc()
    call tic()
    call dfftw_execute_dft( c2cplanbwd, coutput, cinputn)
    cinputn = cinputn*normfact
    call toc()
    print*, "Error:", maxval(abs(cinput - cinputn))
    print*, "=========================="
    call dfftw_destroy_plan ( c2cplanfwd )
    call dfftw_destroy_plan ( c2cplanbwd )
   
    call dfftw_plan_many_dft(c2cplanfwd, 1, 2*nx, ny*nz, cinputn,2*nx, 1, 2*nx, &
                            cinputn, 2*nx, 1, 2*nx, FFTW_FORWARD, FFTW_EXHAUSTIVE)
    call dfftw_plan_many_dft(c2cplanbwd, 1, 2*nx, ny*nz, cinputn,2*nx, 1, 2*nx, &
                            cinputn, 2*nx, 1, 2*nx, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
    cinputn = cinput
        
    print*, "TEST 4: In-place FFT with 2x extension" 
    print*, "--------------------------"
    call tic()
    call dfftw_execute_dft( c2cplanfwd, cinputn, cinputn)
    call toc()
    call tic()
    call dfftw_execute_dft( c2cplanbwd, cinputn, cinputn)
    cinputn = cinputn*normfact
    call toc()
    print*, "Error:", maxval(abs(cinput - cinputn))
    print*, "=========================="
    call dfftw_destroy_plan ( c2cplanfwd )
    call dfftw_destroy_plan ( c2cplanbwd )
   
    deallocate(cinput, cinputn, coutput)

end program 
