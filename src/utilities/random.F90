#ifdef __GFORTRAN__
module ifport
end module ifport
#endif

module random
    use kind_parameters, only: rkind
    use constants, only: zero, one, two, pi
   
    implicit none
    private 
    public :: uniform_random, gaussian_random, randperm

    interface uniform_random
        module procedure unrand3R, unrand2R, unrand1R, unrand0R
    end interface 
    interface gaussian_random
        module procedure grand3R, grand2R, grand1R
    end interface
contains

    subroutine grand3R(array,mu,sigma,seed)
        real(rkind), dimension(:,:,:), intent(inout) :: array
        real(rkind), intent(in) :: mu, sigma
        integer, intent(in), optional :: seed
        real(rkind), dimension(:,:,:), allocatable, target :: uarr
        real(rkind), dimension(:,:,:), pointer :: uarr1, uarr2      
 
        allocate(uarr(size(array,1),size(array,2),2*size(array,3))) 
        
        if (present(seed)) then
            call uniform_random(uarr,zero,one,seed)
        else
            call uniform_random(uarr,zero,one)
        end if 
        uarr1 => uarr(:,:,1:size(array,3)) 
        uarr2 => uarr(:,:,size(array,3)+1:2*size(array,3)) 
       
        ! Now just use the Box-Muller algorithm 
        array = sqrt(-two*log(uarr1))*cos(two*pi*uarr2)        

        array = mu + sigma*array
        nullify(uarr1, uarr2)
        deallocate(uarr)
    end subroutine


    subroutine grand2R(array,mu,sigma,seed)
        real(rkind), dimension(:,:), intent(inout) :: array
        real(rkind), intent(in) :: mu, sigma
        integer, intent(in), optional :: seed
        
        real(rkind), dimension(:,:), allocatable, target :: uarr
        real(rkind), dimension(:,:), pointer :: uarr1, uarr2      
 
        allocate(uarr(size(array,1),2*size(array,2))) 
        
        if (present(seed)) then
            call uniform_random(uarr,zero,one,seed)
        else
            call uniform_random(uarr,zero,one)
        end if 
        
        
        uarr1 => uarr(:,1:size(array,2)) 
        uarr2 => uarr(:,size(array,2)+1:2*size(array,2)) 
       
        ! Now just use the Box-Muller algorithm 
        array = sqrt(-two*log(uarr1))*cos(two*pi*uarr2)        

        array = mu + sigma*array
        nullify(uarr1, uarr2)
        nullify(uarr1, uarr2)
        deallocate(uarr)
    end subroutine

    subroutine grand1R(array,mu,sigma,seed)
        real(rkind), dimension(:), intent(inout) :: array
        real(rkind), intent(in) :: mu, sigma
        integer, intent(in), optional :: seed
        
        real(rkind), dimension(:), allocatable, target :: uarr
        real(rkind), dimension(:), pointer :: uarr1, uarr2      
 
        allocate(uarr(2*size(array,1))) 
        
        if (present(seed)) then
            call uniform_random(uarr,zero,one,seed)
        else
            call uniform_random(uarr,zero,one)
        end if 
        
        uarr1 => uarr(1:size(array,1)) 
        uarr2 => uarr(size(array,1)+1:2*size(array,1)) 
       
        ! Now just use the Box-Muller algorithm 
        array = sqrt(-two*log(uarr1))*cos(two*pi*uarr2)        

        array = mu + sigma*array
        nullify(uarr1, uarr2)
        nullify(uarr1, uarr2)
        deallocate(uarr)
    end subroutine
    
    subroutine unrand3R(array,left, right, seed)
        real(rkind), dimension(:,:,:), intent(inout) :: array
        real(rkind), intent(in) :: left, right
        integer, optional, intent(in) :: seed
        real(rkind) :: diff
        integer, allocatable :: iseed(:)
        integer :: n

        diff = right - left
   
        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
        call random_number(array)
        array = diff*array
        array = array + left
    end subroutine

    subroutine unrand2R(array,left, right, seed)
        real(rkind), dimension(:,:), intent(inout) :: array
        real(rkind), intent(in) :: left, right
        real(rkind) :: diff
        integer, optional, intent(in) :: seed
        integer, allocatable :: iseed(:)
        integer :: n

        diff = right - left
    
        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
        call random_number(array)
        array = diff*array
        array = array + left
    end subroutine

    subroutine unrand1R(array,left, right, seed)
        real(rkind), dimension(:), intent(inout) :: array
        real(rkind), intent(in) :: left, right
        real(rkind) :: diff
        integer, optional, intent(in) :: seed
        integer, allocatable :: iseed(:)
        integer :: n

        diff = right - left
    
        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
        call random_number(array)
        array = diff*array
        array = array + left
    end subroutine

    subroutine unrand0R(val,left, right, seed)
        real(rkind), intent(inout) :: val
        real(rkind), intent(in) :: left, right
        integer, intent(in), optional :: seed
        integer, allocatable :: iseed(:)
        real(rkind) :: diff
        integer :: n

        diff = right - left

        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
    
        call random_number(val)
        val = diff*val
        val = val + left
    end subroutine

    subroutine randperm(N,seed,arrout)
      use fortran_assert, only: assert
      use arrayTools,     only: swap
      ! Return an array of size N with elements ranging 1:N without
      ! repitition and randomly ordered
      integer, intent(in) :: N, seed
      integer, dimension(N), intent(out) :: arrout
      integer :: i, a, newseed
      real(rkind) :: val

      do i = 1,N
        arrout(i) = i
      end do
      
      newseed = seed
      do i = 1,N-1
        ! Ensure seed isn't close to max possible integer value
        call assert(abs(huge(0) - newseed) > 1000000,&
          'Seed value exceeds maximum representable integer -- random.F90')
        
        call uniform_random(val,real(i,rkind),real(N,rkind),newseed)
        a = nint(val)

        ! Swap arrout(i) with arrout(a)
        call swap(arrout,a,i)

        ! Update seed
        newseed = seed + 125932
      end do

    end subroutine

    subroutine init_random_seed()
        ! Taken from GNU 
        use iso_fortran_env, only: int64
        use ifport
        implicit none
        integer, allocatable :: iseed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
        
        call random_seed(size = n)
        allocate(iseed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
             form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
           read(un) iseed
           close(un)
        else
           ! Fallback to XOR:ing the current time and pid. The PID is
           ! useful in case one launches multiple instances of the same
           ! program in parallel.
           call system_clock(t)
           if (t == 0) then
              call date_and_time(values=dt)
              t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                   + dt(3) * 24_int64 * 60 * 60 * 1000 &
                   + dt(5) * 60 * 60 * 1000 &
                   + dt(6) * 60 * 1000 + dt(7) * 1000 &
                   + dt(8)
           end if
           pid = getpid()
           t = ieor(t, int(pid, kind(t)))
           do i = 1, n
              iseed(i) = lcg(t)
           end do
        end if
        call random_seed(put=iseed)
    contains
        ! This simple PRNG might not be good enough for real work, but is
        ! sufficient for seeding a better PRNG.
        function lcg(s)
          integer :: lcg
          integer(int64) :: s
          if (s == 0) then
             s = 104729
          else
             s = mod(s, 4294967296_int64)
          end if
          s = mod(s * 279470273_int64, 4294967291_int64)
          lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
    end subroutine init_random_seed

end module 
