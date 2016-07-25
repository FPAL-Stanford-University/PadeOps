program test_dct1d
    use kind_parameters, only: rkind
    use dct_1d_stuff, only: dct1d
    
    real(rkind), dimension(:), allocatable :: input, output
    type(dct1d) :: mydct
    integer :: n = 8
    integer :: i 

    allocate(input(n), output(n))

    do i = 1,n
        input(i) = real(i,rkind)
    end do 
    
    call mydct%init(n,1,1)

    call mydct%dctx(input,output)
    print*, output 

    call mydct%idctx(output,input)
    print*, input


    call mydct%destroy()

    deallocate(input, output)

end program 
