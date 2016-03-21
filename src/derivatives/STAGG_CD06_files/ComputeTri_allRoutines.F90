    subroutine ComputeTriD1_E2C(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 9._rkind/62._rkind

        n = this%n

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dup(n) = zero; dg(n) = (one - alpha)
        else
            dg(n) = (one + alpha); dup(n) = zero
        end if 

        if (this%isBotEven) then
            ddn(1) = zero; dg(1) = (one - alpha)
        else
            ddn(1) = zero; dg(1) = (one + alpha)
        end if 

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriD1_E2C(:,1) = ddn*den; this%TriD1_E2C(:,2) = den; this%TriD1_E2C(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine
    
    subroutine ComputeTriD1_C2E(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 9._rkind/62._rkind

        n = this%nE

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dup(n) = zero; ddn(n) = zero
        else
            ddn(n) = two*alpha; dup(n) = zero 
        end if 

        if (this%isBotEven) then
            ddn(1) = zero; dup(1) = zero
        else
            ddn(1) = zero; dup(1) = two*alpha
        end if 

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriD1_C2E(:,1) = ddn*den; this%TriD1_C2E(:,2) = den; this%TriD1_C2E(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine
    
    subroutine ComputeTriD1_C2C(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 1._rkind/3._rkind

        n = this%n

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dg(n) = (one - alpha)
        else
            dg(n) = (one + alpha)
        end if 

        if (this%isBotEven) then
            dg(1) = (one - alpha)
        else
            dg(1) = (one + alpha)
        end if

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriD1_C2C(:,1) = ddn*den; this%TriD1_C2C(:,2) = den; this%TriD1_C2C(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine

    subroutine ComputeTriD1_E2E(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 1._rkind/3._rkind

        n = this%nE

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dg(n) =  one; ddn(n) = zero
        else
            dg(n) =  one; ddn(n) = two*alpha
        end if 

        if (this%isBotEven) then
            dg(1) = one; dup(1) = zero
        else
            dg(1) = one; dup(1) = two*alpha 
        end if

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriD1_E2E(:,1) = ddn*den; this%TriD1_E2E(:,2) = den; this%TriD1_E2E(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine

    subroutine ComputeTriInterp_C2E(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 3._rkind/10._rkind

        n = this%nE

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            ddn(n) = two*alpha; dg(n) = one
        else
            ddn(n) = zero; dg(n) = one 
        end if 

        if (this%isBotEven) then
            dup(1) = two*alpha; dg(1) = one
        else
            dup(1) = zero; dg(1) = one
        end if 

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriInterp_C2E(:,1) = ddn*den; this%TriInterp_C2E(:,2) = den; this%TriInterp_C2E(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine


    subroutine ComputeTriInterp_E2C(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 3._rkind/10._rkind

        n = this%n

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dg(n) = (one + alpha)
        else
            dg(n) = (one - alpha)
        end if 

        if (this%isBotEven) then
            dg(1) = (one + alpha)
        else
            dg(1) = (one - alpha)
        end if 

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriInterp_E2C(:,1) = ddn*den; this%TriInterp_E2C(:,2) = den; this%TriInterp_E2C(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine


    subroutine ComputeTriD2_E2E(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 2._rkind/11._rkind

        n = this%nE

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dg(n) =  one; ddn(n) = two*alpha
        else
            dg(n) =  one; ddn(n) = zero
        end if 

        if (this%isBotEven) then
            dg(1) = one; dup(1) = two*alpha
        else
            dg(1) = one; dup(1) = zero 
        end if

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriD2_E2E(:,1) = ddn*den; this%TriD2_E2E(:,2) = den; this%TriD2_E2E(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine

    subroutine ComputeTriD2_C2C(this)
        class (cd06stagg), intent(inout) :: this
        integer             :: i, n 
        real(rkind), dimension(:), allocatable :: ddn, dg, dup, cp, den
        real(rkind), parameter :: alpha = 2._rkind/11._rkind

        n = this%n

        allocate(ddn(n));  allocate(dg(n)); allocate(dup(n)); allocate(cp(n)); allocate(den(n))
        ddn  = alpha; dg = one; dup = alpha 

        if (this%isTopEven) then
            dg(n) =  one + alpha
        else
            dg(n) =  one - alpha
        end if 

        if (this%isBotEven) then
            dg(1) = one + alpha
        else
            dg(1) = one - alpha 
        end if

        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))

        this%TriD2_C2C(:,1) = ddn*den; this%TriD2_C2C(:,2) = den; this%TriD2_C2C(:,3) = cp

        deallocate(ddn, dg, dup, den, cp)

    end subroutine
