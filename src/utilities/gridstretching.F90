module gridstretching
    use kind_parameters, only: rkind
    use constants, only: one, two, half
    use decomp_2d, only: decomp_info
    use exits,     only: GracefulExit
    implicit none

    !interface destroy_buffs
    !    module procedure destroy_buffs_real, destroy_buffs_complex 
    !end interface
    !
    !interface alloc_buffs
    !    module procedure alloc_buffs_real, alloc_buffs_complex 
    !end interface

contains
    
    subroutine init_stretched_grid(decomp, y, eta, ymetric, ymetric_flag, param1, param2, param3, param4)
    use constants,        only: zero, half, one
    use exits,            only: GracefulExit, message, nancheck

    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:),   intent(inout) :: y
    real(rkind), dimension(:,:,:  ), intent(inout) :: eta
    logical,                         intent(in   ) :: ymetric
    integer,                         intent(in)    :: ymetric_flag
    real(rkind),                     intent(in)    :: param1, param2, param3, param4
    integer     :: i,j,k
    real(rkind) :: yfocus, ytau, ystart, yh, alpha, beta, top,bot
    real(rkind) :: yfocus_adj, num, den, BB, yuniform_adj, BB2

    ! concentrate towards the center -- Pletcher, Tannehill, Anderson
    ! (Section 5.6, Transformation 3, pg. 332) 
    if(ymetric_flag==1) then
       ! Concentrate towards centre
       yfocus = param1; ytau = param2; ystart = param3; yh = param4
       yfocus_adj = yfocus - ystart
       num = one + (yfocus_adj/yh) * (exp( ytau) - one)
       den = one + (yfocus_adj/yh) * (exp(-ytau) - one)
       BB  = half/ytau*log(num/den)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = (eta(i,j,k) - ystart) !/ yh
                  num = sinh(ytau*BB)
                  y(i,j,k) = yfocus_adj * (one + sinh(ytau * (yuniform_adj/yh-BB))/num) + ystart
             end do
          end do
       end do
    elseif(ymetric_flag==2) then
       ! concentrate towards the two ends
       alpha = param1; beta = param2; ystart = param3; yh = param4
       BB   = (beta + 1) / (beta - 1)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = (eta(i,j,k) - ystart) / yh
                  BB2 = BB ** ( (yuniform_adj-alpha) / (1-alpha) )
                  num = ((beta+2*alpha)*BB2 - beta + 2*alpha ) * yh
                  y(i,j,k) = num/( (2*alpha+1)*(1+BB2) )   + ystart
             end do
          end do
       end do
    elseif(ymetric_flag==3) then
       ! concentrate at one end (at the starting/left/bottom side)
       beta = param2; ystart = param3; yh = param4
       BB   = (beta + 1) / (beta - 1)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = (eta(i,j,k) - ystart)/yh
                  BB2 = BB ** (1 - yuniform_adj)
                  num = (beta+1)-((beta-1)*BB2)
                  den = one + BB2
                  y(i,j,k) = yh * (num/den) + ystart
             end do
          end do
       end do
    elseif(ymetric_flag==4) then
       ! concentrate at arbitrary point (given by yfocus)
       beta = param2; ystart = param3; yh = param4; yfocus = param1 + abs(ystart)
       top = 1+ ((yfocus/yh)*(exp(beta)-1))
       bot = 1+ ((yfocus/yh)*(exp(-beta)-1))
       BB = (log(top/bot))/(2*beta)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = (eta(i,j,k) - ystart)/yh
                  num = sinh(beta*(yuniform_adj-BB))
                  den = sinh(beta*BB)
                  y(i,j,k) = yfocus * (1+(num/den))
             end do
          end do
       end do
       
       !call GracefulExit("flag = 3 (concentrate at arbitrary point) is incomplete right now",21)

    elseif(ymetric_flag==10) then
       ! finite-difference evaluation of metrics (reduces order of accuracy)
       call GracefulExit("flag = 4 (finite-difference evaluation of metrics) is incomplete right now",21)
    endif
    end subroutine
    subroutine init_metrics_der(nxp, nyp, nzp, xstretch, xuniform, flag, params, dxudxs, dxudxs_sq, d2xudxs2) 
        integer,                       intent(in)  :: nxp, nyp, nzp, flag
        real(rkind), dimension(5),     intent(in)  :: params
        real(rkind), dimension(:,:,:), intent(in)  :: xstretch, xuniform
        real(rkind), dimension(:,:,:), intent(out) :: dxudxs, dxudxs_sq, d2xudxs2

        real(rkind) :: xfocus, xfocus_adj, xtau, hh, xstart, num, num1, den, BB, BB2, xstretch_loc
        real(rkind) :: xbyxfocm1, xuniform_adj, alpha, beta, top, bot
        integer     :: i, j, k

        if(flag==1) then
           ! concentrate towards the center -- Pletcher, Tannehill, Anderson
           ! (Section 5.6, Transformation 3, pg. 332) 
           xfocus = params(1);  xtau   = params(2);  xstart = params(3); hh = params(4)
           !print '(5(e19.12,1x))', params(:)
           !print '(5(e19.12,1x))', xfocus, xtau, xstart, hh
           xfocus_adj = xfocus - xstart
           num = one + (xfocus_adj/hh) * (exp( xtau) - one)
           den = one + (xfocus_adj/hh) * (exp(-xtau) - one)
           BB = half/xtau*log(num/den)
           do k = 1, nzp
              do j = 1, nyp
                 do i = 1, nxp
                    ! adjust for starting point
                    xuniform_adj = (xuniform(i,j,k) - xstart) !/ hh

                    ! stretched location
                    num = sinh(xtau*BB)
                    xstretch_loc = xfocus_adj * (one + sinh(xtau * (xuniform_adj/hh-BB))/num) + xstart
                    xbyxfocm1 = (xstretch_loc-xstart)/xfocus_adj - one

                    ! metric for first derivative
                    dxudxs(i,j,k) = num * hh / (xtau * xfocus_adj * sqrt(one + (xbyxfocm1 * num)**2))

                    ! square of the metric for first derivative
                    dxudxs_sq(i,j,k) = dxudxs(i,j,k)**2

                    ! metric for second derivative
                    d2xudxs2(i,j,k) = -hh * num**3 * xbyxfocm1 / (xtau*xfocus_adj*xfocus_adj * (one + (num*xbyxfocm1)**2)**1.5d0)

                    !print '(i5,1x,4(e19.12,1x))', j, xstretch(1,j,1), xuniform(1,j,1), dxudxs(1,j,1), d2xudxs2(1,j,1)
                    ! compare with xstretch specified in meshgen
                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        !print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 1; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
        elseif(flag==2) then
           ! concentrate towards the two ends(alpha = 0.5) or one side at xend(alpha = 0)
           ! Pletcher, Tannehill, Anderson (Section 5.6, Transformation 2, pg. 335) 
           alpha  = params(1);  beta   = params(2);  xstart = params(3); hh = params(4)
           !print '(5(e19.12,1x))', params(:)
           !print '(5(e19.12,1x))', alpha, beta, xstart, hh
           BB = (beta + one) / (beta - one)
           num1 = two*beta * (two*alpha+one) * (one-alpha) / log(BB)
           do k = 1, nzp
              do j = 1, nyp
                 do i = 1, nxp
                    ! adjust for starting point
                    xuniform_adj = (xuniform(i,j,k) - xstart) / hh
                    ! stretched location
                    BB2 = BB ** ( (xuniform_adj-alpha) / (one-alpha) )
                    num = ((beta+two*alpha)*BB2 - beta + two*alpha ) * hh
                    xstretch_loc = num/( (two*alpha+1)*(one+BB2) )   + xstart !same strectching but done in locally

                    ! metric for first derivative
                    BB2 = two*alpha - ((xstretch_loc-xstart)/hh) * (two*alpha+one)
                    dxudxs(i,j,k) = num1 / (beta**2 - BB2**2)

                    ! square of the metric for first derivative
                    dxudxs_sq(i,j,k) = dxudxs(i,j,k)**2

                    ! metric for second derivative
                    d2xudxs2(i,j,k) = - ((num1/hh) * two * (two*alpha+one) * BB2) / ((beta**2 - BB2**2)**2)
                    !print '(i5,1x,4(e19.12,1x))', j, xstretch(1,j,1), xuniform(1,j,1), dxudxs(1,j,1), d2xudxs2(1,j,1)
                    ! compare with xstretch specified in meshgen
                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 2; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
        elseif(flag==3) then
           ! concentrate towards the start at xstart
           ! (Section 5.6, Transformation 1, pg. 335) 
           alpha  = params(1);  beta   = params(2);  xstart = params(3); hh = params(4)
           !print '(5(e19.12,1x))', params(:)
           !print '(5(e19.12,1x))', alpha, beta, xstart, hh
           BB = (beta + one) / (beta - one)
           num1 = two*beta / log(BB)
           do k = 1, nzp
              do j = 1, nyp
                 do i = 1, nxp
                    !)) adjust for starting point
                    xuniform_adj = (xuniform(i,j,k) - xstart) / hh
                    ! stretched location
                    BB2 = BB ** (one - xuniform_adj)
                    num = (beta+one) - ((beta - one) * BB2)
                    xstretch_loc = hh * (num / (BB2+one)) + xstart !same strectching but done in locally

                    ! metric for first derivative y=(xstretch_loc-xstart)
                    BB2 = (one - ((xstretch_loc-xstart)/hh)) !(1-y/h) 
                    dxudxs(i,j,k) = num1 / ((beta**2 - BB2**2)*hh)

                    ! square of the metric for first derivative
                    dxudxs_sq(i,j,k) = dxudxs(i,j,k)**2

                    ! metric for second derivative
                    d2xudxs2(i,j,k) = - ((num1/hh**2) * two * BB2) / ((beta**2 - BB2**2)**2)
                    !print '(i5,1x,4(e19.12,1x))', j, xstretch(1,j,1), xuniform(1,j,1), dxudxs(1,j,1), d2xudxs2(1,j,1)

                    ! compare with xstretch specified in meshgen
                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 3; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
        elseif(flag==4) then
           ! concentrate at any arbitary point
           ! (Section 5.6, Transformation 3, pg. 334)
           beta = params(2); xstart = params(3); hh = params(4); xfocus = params(1) + abs(xstart)
           top = one+ ((xfocus/hh)*(exp(beta)-one))
           bot = one+ ((xfocus/hh)*(exp(-beta)-one))
           BB = (log(top/bot))/(two*beta)
           do k = 1, nzp
              do j = 1, nyp
                 do i = 1, nxp
                    ! adjust for starting point
                    xuniform_adj = (xuniform(i,j,k) - xstart)/hh
                    num = sinh(beta*(xuniform_adj-BB))
                    den = sinh(beta*BB)
                    xstretch_loc = xfocus*(one+(num/den))

                    ! metric for first derivative
                    BB2 = one + (((xstretch_loc/xfocus)-1)*den)**2
                    dxudxs(i,j,k) = den/(beta*xfocus*dsqrt(BB))

                    ! square of the metric for first derivative
                    dxudxs_sq(i,j,k) = dxudxs(i,j,k)**2

                    ! metric for second derivative
                    d2xudxs2(i,j,k) = -((den**3)*((xstretch_loc/xfocus)-one))/(beta*xfocus*(BB**1.5))

                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 3; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
           ! stretching to lower side only (Formula .... on page ... of .....)
           ! call GracefulExit("flag = 4 (stretching to the lower side) is incomplete right now",21)
        elseif(flag==10) then
           ! finite-difference evaluation of metrics (reduces order of accuracy)
            call GracefulExit("flag = 10 (finite-difference evaluation of metrics) is incomplete right now",21)
        else
           ! incorrect input for flag
            call GracefulExit("flag must be 1, 2, 3, 4 or 10",21)
        endif

    end subroutine 

   !pure function mytrapz(x,f) result(intf)
   !     real(rkind), dimension(:), intent(in) :: x, f
   !     real(rkind) :: intf
   !     real(rkind) :: dx

   !     dx = x(2) - x(1)
   !     intf = sum(f)*dx
   ! end function 

   ! pure function logspace(x_min,x_max,n) result(x)
   !     integer, intent(in) :: n                    ! Desired size
   !     real(rkind), dimension(n) :: x ! Output array
   !     real(rkind), intent(in) :: x_min, x_max     ! Left and Right bounds in POWERS OF 10

   !     real(rkind) :: step
   !     integer :: i
   !     real(rkind), dimension(n) :: xtemp

   !     if (n .le. 1) then
   !         x = ten**(x_min)
   !         return 
   !     end if 
   !     step = (x_max - x_min)/(real(n,rkind)- one)
   !     xtemp(1) = x_min

   !     do i = 2,n
   !         xtemp(i) = xtemp(i-1) + step
   !     end do

   !      x = (ten)**(xtemp)

   ! end function

   ! pure function linspace(x_min,x_max,n) result(x)
   !     integer, intent(in) :: n
   !     real(rkind), dimension(n) :: x
   !     real(rkind), intent(in) :: x_min, x_max
   !     
   !     real(rkind) :: step
   !     integer :: i

   !     if (n .le. 1) then
   !         x = x_min
   !         return
   !     end if 
   !     step = (x_max - x_min)/(real(n,rkind) - one)
   !     
   !     x(1) = x_min
   !     do i = 2,n
   !         x(i) = x(i-1) + step
   !     end do 

   ! end function
   ! 
   ! subroutine destroy_buffs_real(buff)
   !     real(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: buff

   !     if (allocated(buff)) deallocate(buff)
   ! end subroutine

   ! subroutine destroy_buffs_complex(buff)
   !     complex(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: buff

   !     if (allocated(buff)) deallocate(buff)
   ! end subroutine

   ! subroutine alloc_buffs_real(buff,vars,dir,decomp)
   !     character(len=1), intent(in) :: dir
   !     class(decomp_info), intent(in) :: decomp
   !     integer, intent(in) :: vars
   !     real(rkind), dimension(:,:,:,:), allocatable, intent(out) :: buff

   !     if (allocated(buff)) deallocate(buff)

   !     select case (dir)
   !     case ("x")
   !         allocate(buff(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),vars))
   !     case("y")
   !         allocate(buff(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),vars))
   !     case("z")
   !         allocate(buff(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),vars))
   !     case default
   !         call GracefulExit("Incorrect direction selected in ALLOC_BUFFS_REAL subroutine", 13123)
   !     end select

   ! end subroutine

   ! subroutine alloc_buffs_complex(buff,vars,dir,decomp)
   !     character(len=1), intent(in) :: dir
   !     class(decomp_info), intent(in) :: decomp
   !     integer, intent(in) :: vars
   !     complex(rkind), dimension(:,:,:,:), allocatable, intent(out) :: buff

   !     if (allocated(buff)) deallocate(buff)

   !     select case (dir)
   !     case ("x")
   !         allocate(buff(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),vars))
   !     case("y")
   !         allocate(buff(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),vars))
   !     case("z")
   !         allocate(buff(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),vars))
   !     case default
   !         call GracefulExit("Incorrect direction selected in ALLOC_BUFFS_REAL subroutine", 13123)
   !     end select

   ! end subroutine

   ! pure subroutine upsample_Periodic_1d(vecin, vecout, nx)
   !     integer, intent(in) :: nx
   !     real(rkind), intent(in) , dimension(nx) :: vecin
   !     real(rkind), intent(out), dimension(2*nx) :: vecout
   !     integer :: i

   !     do i = 1,nx
   !         vecout(2*i - 1) = vecin(i)
   !     end do 

   !     do i = 2,2*nx-2,2
   !         vecout(i) = 0.5d0*(vecout(i-1) + vecout(i+1))
   !     end do 
   !     vecout(2*nx) = 0.5d0*(vecout(2*nx-1) + vecout(1))
   ! end subroutine
end module 
