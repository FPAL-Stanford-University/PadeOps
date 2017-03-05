    !!!!!!!!!!!!!! FRINGE METHOD STUFF !!!!!!!!!!!!!!!!!
    subroutine setupFringe(this,Fringe_xst, Fringe_xen, Fringe_delta_st, Fringe_delta_en)
        class(igrid), intent(inout) :: this
        real(rkind), intent(in) :: Fringe_xst, Fringe_xen, Fringe_delta_st, Fringe_delta_en
        real(rkind) :: Lx
        integer :: j, k
        real(rkind), dimension(:), allocatable :: x1, x2, Fringe_func, S1, S2


        ! Get the domain length
        Lx = p_maxval(maxval(this%mesh(:,:,:,1))) + this%dx
        allocate(this%Fringe_kernel_cells_x(this%nx, this%gpC%xsz(2), this%gpC%xsz(3)))
        allocate(this%Fringe_kernel_edges_x(this%nx, this%gpE%xsz(2), this%gpE%xsz(3)))

        allocate(x1         (this%nx))
        allocate(x2         (this%nx))
        allocate(S1         (this%nx))
        allocate(S2         (this%nx))
        allocate(Fringe_func(this%nx))
        
        x1 = ((this%mesh(:,1,1,1) -  Fringe_xst)/Fringe_delta_xst)
        x2 = ((this%mesh(:,1,1,1) -  Fringe_xen)/Fringe_delta_xen) + 1.d0
        call S_fringe(x1, S1)
        call S_fringe(x2, S2)
        Fringe_func = S1 - S2

        do k = 1,this%gpC%xsz(3)
           do j = 1,this%gpC%xsz(2)
               this%Fringe_kernel_cells_x(:,j,k) = Fringe_func    
           end do 
        end do

        do k = 1,this%gpE%xsz(3)
           do j = 1,this%gpE%xsz(2)
               this%Fringe_kernel_edges_x(:,j,k) = Fringe_func    
           end do 
        end do


    end subroutine 

    pure subroutine S_fringe(x, output)
       real(rkind), dimension(:), intent(in)    :: x
       real(rkind), dimension(:), intent(out)   :: output
       integer :: i

       do i = 1,size(x)
         if (x(i) < 0.d0) then
            output(i) = 0.d0
         else if (x(i) > 1.d0) then
            output(i) = 1.d0
         else
            output(i) = 1.d0/(1.d0 + exp((1.d0/(x(i) - 1.d0)) + (1.d0/(x(i)))))
         end if
       end do

    end subroutine
