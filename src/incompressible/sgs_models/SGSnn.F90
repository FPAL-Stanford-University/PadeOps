subroutine init_SGSnn(this,NNtype)
   use fortran_assert, only: assert
   class(sgs_igrid), intent(inout) :: this
   integer, intent(in) :: NNtype

   this%isEddyViscosityModel = .false.
     
   allocate(this%strain2(       1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))
   allocate(this%strain2_m(     1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))
   allocate(this%rot2(          1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))
   allocate(this%strainrot_m(   1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))
   allocate(this%rotstrainrot_m(1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))
   allocate(this%strain2rot_m(  1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))
   allocate(this%rotstrain2_m(  1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3))

   select case (NNtype)
   case (1) ! Unet
     allocate(this%invariants(1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),7))
     allocate(this%delta(1,8,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
   case (2) ! LSTM
     allocate(this%invariantsLSTM(1,10,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),7))
     allocate(this%deltaLSTM(     1,10,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),8))
   case default
       call assert(.false.,'NNtype must be 1 or 2')
   end select

   call message(1,"NN model initialized")
end subroutine

subroutine getSijRijforNNmod(duidxj, Sij, Rij, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in) :: duidxj
   real(rkind), dimension(nxL,nyL,nzL,9), intent(out) :: Sij, Rij
   integer :: i, j, k

   do k = 1,nzL
      do j = 1,nyL
         !$omp simd
         do i = 1,nxL
            Sij(i,j,k,1) = duidxj(i,j,k,1) ! S11 = dudx
            Sij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) + duidxj(i,j,k,4)) ! S12 = 0.5*(dudy + dvdx)
            Sij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) + duidxj(i,j,k,7)) ! S13 = 0.5*(dudz + dwdx)
           
            Sij(i,j,k,5) = duidxj(i,j,k,5) ! S22 = dvdy
            Sij(i,j,k,6) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
            
            Sij(i,j,k,9) = duidxj(i,j,k,9) ! S33 = dwdz

            Rij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) - duidxj(i,j,k,4)) ! R12 = 0.5*(dudy - dvdx)
            Rij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) - duidxj(i,j,k,7)) ! R13 = 0.5*(dudz - dwdx)
            Rij(i,j,k,6) = 0.5d0*(duidxj(i,j,k,6) - duidxj(i,j,k,8)) ! R23 = 0.5*(dvdz - dwdy)
         end do 
      end do 
   end do 
   Sij(:,:,:,4) = Sij(:,:,:,2) ! S21 = S12
   Sij(:,:,:,7) = Sij(:,:,:,3) ! S31 = S13
   Sij(:,:,:,8) = Sij(:,:,:,6) ! S32 = S23

   Rij(:,:,:,1) = 0.d0
   Rij(:,:,:,5) = 0.d0
   Rij(:,:,:,9) = 0.d0

   Rij(:,:,:,4) = -Rij(:,:,:,2) ! R21 = -R12
   Rij(:,:,:,7) = -Rij(:,:,:,3) ! R31 = -R13
   Rij(:,:,:,8) = -Rij(:,:,:,6) ! R32 = -R23
end subroutine

subroutine readInvariants(tid,runID,datadir,gp,dat)
    integer, intent(in) :: tid, runID
    character(len=*), intent(in) :: datadir
    class(decomp_info), intent(in) :: gp
    real(rkind), dimension(:,:,:,:), intent(out) :: dat
    character(len=clen) :: fname
    integer :: n
   
    do n = 1,size(dat,4) 
        write(fname,'(A,I2.2,A4,I1,A2,I6.6,A4)') trim(datadir)//'/Run',runID,'_inv',n,'_t',tid,'.out'
        call decomp_2d_read_one(1,dat(:,:,:,n),trim(fname),gp)
    end do
end subroutine

subroutine writeInvariants(tid,runID,datadir,gp,dat)
    integer, intent(in) :: tid, runID
    character(len=*), intent(in) :: datadir
    class(decomp_info), intent(in) :: gp
    real(rkind), dimension(:,:,:,:), intent(out) :: dat
    character(len=clen) :: fname
    integer :: n
   
    do n = 1,size(dat,4) 
        write(fname,'(A,I2.2,A4,I1,A2,I6.6,A4)') trim(datadir)//'/Run',runID,'_inv',n,'_t',tid,'.out'
        call decomp_2d_write_one(1,dat(:,:,:,n),trim(fname),gp)
    end do
end subroutine

subroutine compute_tauij_NN(this,tidNow)
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr) :: c_p
  real(c_float), pointer :: f_p(:,:,:,:,:)

  class(sgs_igrid), intent(inout) :: this
  integer, intent(in) :: tidNow
  integer :: tidPast, n
  integer :: i, j, k
  real, dimension(3, 3) :: arr_temp1
  real, dimension(3, 3) :: arr_temp2
  real, dimension(3, 3) :: arr_temp3
  real, dimension(3, 3) :: arr_temp4
  real, dimension(3, 3) :: arr_temp5
  real, dimension(3, 3) :: arr_temp6
  real, dimension(9) :: S_ij_a
  real, dimension(9) :: R_ij_a
  integer, dimension(9) :: timelist
  real :: del

  interface
      function loadtensorflow(delta, train_load, strain2, strain2_m, rot2, strainrot_m, rotstrainrot_m, strain2rot_m, &
        rotstrain2_m) bind(c)
        import :: c_ptr
        type(c_ptr) :: loadtensorflow
          real, dimension(2) :: train_load
          real, dimension(2) :: delta
          real, dimension(2) :: strain2
          real, dimension(2) :: strain2_m
          real, dimension(2) :: rot2
          real, dimension(2) :: strainrot_m
          real, dimension(2) :: rotstrainrot_m
          real, dimension(2) :: strain2rot_m
          real, dimension(2) :: rotstrain2_m
      end function loadtensorflow
  end interface

  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*7) :: train_loadin
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*8) :: deltain
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*7*10) :: train_loadinLSTM

  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2in
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2_min
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rot2in
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strainrot_min
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rotstrainrot_min
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2rot_min
  real, dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rotstrain2_min

 
  del = 2*3.141592/1024*32
  ! TODO: Andy, compute your invariants, call your NN, compute tauij
  do k = 1, this%gpC%xsz(3)
    do j = 1, this%gpC%xsz(2)
      do i = 1, this%gpC%xsz(1)
        S_ij_a(1) = this%S_ij_C(i, j, k, 1)
        S_ij_a(2) = this%S_ij_C(i, j, k, 2)
        S_ij_a(3) = this%S_ij_E(i, j, k, 3)
        S_ij_a(4) = this%S_ij_C(i, j, k, 2)
        S_ij_a(5) = this%S_ij_C(i, j, k, 4)
        S_ij_a(6) = this%S_ij_E(i, j, k, 5)
        S_ij_a(7) = this%S_ij_E(i, j, k, 3)
        S_ij_a(8) = this%S_ij_E(i, j, k, 5)
        S_ij_a(9) = this%S_ij_C(i, j, k, 6)

        R_ij_a(1) = this%R_ij_C(i, j, k, 1)
        R_ij_a(2) = this%R_ij_C(i, j, k, 2)
        R_ij_a(3) = this%R_ij_E(i, j, k, 3)
        R_ij_a(4) = this%R_ij_C(i, j, k, 2)
        R_ij_a(5) = this%R_ij_C(i, j, k, 4)
        R_ij_a(6) = this%R_ij_E(i, j, k, 5)
        R_ij_a(7) = this%R_ij_E(i, j, k, 3)
        R_ij_a(8) = this%R_ij_E(i, j, k, 5)
        R_ij_a(9) = this%R_ij_C(i, j, k, 6)

        arr_temp1 = matmul(reshape(S_ij_a(:),[3,3]),reshape(S_ij_a(:),[3, 3]))
        arr_temp3 = matmul(reshape(S_ij_a(:),[3,3]),arr_temp1)
        arr_temp2 = matmul(reshape(R_ij_a(:),[3,3]),reshape(R_ij_a(:),[3, 3]))
        arr_temp4 = matmul(reshape(S_ij_a(:),[3,3]),arr_temp2)
        arr_temp5 = matmul(arr_temp1,arr_temp2)
        arr_temp6 = matmul(matmul(arr_temp5,reshape(S_ij_a(:),[3,3])),reshape(R_ij_a(:),[3,3]))
        if (this%NNtype == 1) then
          this%invariants(1, i, j, k, 1) = S_ij_a(1)+S_ij_a(5)+S_ij_a(9)
          this%invariants(1, i, j, k, 2) = arr_temp1(1,1)+arr_temp1(2,2)+arr_temp1(3,3)
          this%invariants(1, i, j, k, 3) = arr_temp2(1,1)+arr_temp2(2,2)+arr_temp2(3,3)
          this%invariants(1, i, j, k, 4) = arr_temp3(1,1)+arr_temp3(2,2)+arr_temp3(3,3)
          this%invariants(1, i, j, k, 5) = arr_temp4(1,1)+arr_temp4(2,2)+arr_temp4(3,3)
          this%invariants(1, i, j, k, 6) = arr_temp5(1,1)+arr_temp5(2,2)+arr_temp5(3,3)
          this%invariants(1, i, j, k, 7) = arr_temp6(1,1)+arr_temp6(2,2)+arr_temp6(3,3)
          this%strain2(1, i, j, k, :, :) = reshape(S_ij_a(:),[3, 3])
          this%strain2_m(1, i, j, k, :, :) = matmul(reshape(S_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3]))
          this%rot2(1, i, j, k, :, :) =  matmul(reshape(R_ij_a(:),[3, 3]),reshape(R_ij_a(:),[3, 3]))
          this%strainrot_m(1, i, j, k, :, :) = matmul(reshape(S_ij_a(:),[3, 3]),reshape(R_ij_a(:),[3, 3]))
          this%rotstrainrot_m(1, i, j, k, :, :) = matmul(matmul(reshape(R_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3])),reshape(R_ij_a(:),[3, 3]))
          this%strain2rot_m(1, i, j, k, :, :) = matmul(matmul(reshape(S_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3])),reshape(R_ij_a(:),[3, 3]))
          this%rotstrain2_m(1, i, j, k, :, :) = matmul(matmul(reshape(R_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3])),reshape(S_ij_a(:),[3, 3]))
          this%delta(1, :, i, j, k) = (/del, del, del, del, del, del, del, del/)
        end if
        if (this%NNtype == 2) then
          this%invariantsLSTM(1, 10, i, j, k, 1) = S_ij_a(1)+S_ij_a(5)+S_ij_a(9)
          this%invariantsLSTM(1, 10, i, j, k, 2) = arr_temp1(1,1)+arr_temp1(2,2)+arr_temp1(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 3) = arr_temp2(1,1)+arr_temp2(2,2)+arr_temp2(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 4) = arr_temp3(1,1)+arr_temp3(2,2)+arr_temp3(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 5) = arr_temp4(1,1)+arr_temp4(2,2)+arr_temp4(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 6) = arr_temp5(1,1)+arr_temp5(2,2)+arr_temp5(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 7) = arr_temp6(1,1)+arr_temp6(2,2)+arr_temp6(3,3)
          this%strain2(1, i, j, k, :, :) = reshape(S_ij_a(:),[3, 3])
          this%strain2_m(1, i, j, k, :, :) = matmul(reshape(S_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3]))
          this%rot2(1, i, j, k, :, :) =  matmul(reshape(R_ij_a(:),[3, 3]),reshape(R_ij_a(:),[3, 3]))
          this%strainrot_m(1, i, j, k, :, :) = matmul(reshape(S_ij_a(:),[3, 3]),reshape(R_ij_a(:),[3, 3]))
          this%rotstrainrot_m(1, i, j, k, :, :) = matmul(matmul(reshape(R_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3])),reshape(R_ij_a(:),[3, 3]))
          this%strain2rot_m(1, i, j, k, :, :) = matmul(matmul(reshape(S_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3])),reshape(R_ij_a(:),[3, 3]))
          this%rotstrain2_m(1, i, j, k, :, :) = matmul(matmul(reshape(R_ij_a(:),[3, 3]),reshape(S_ij_a(:),[3, 3])),reshape(S_ij_a(:),[3, 3]))
        end if
      end do
    end do
  end do

  if (this%NNtype == 2) then
      if (tidNow < 512) then
        ! call Unetand onlu feed in LSTM at time step "10"
      end if
      if (tidNow >= 512) then
         timelist(1) = tidNow-511
         timelist(2) = tidNow-256
         timelist(3) = tidNow-64
         timelist(4) = tidNow-32
         timelist(5) = tidNow-16
         timelist(6) = tidNow-8
         timelist(7) = tidNow-4
         timelist(8) = tidNow-2
         timelist(9) = tidNow-1
         do n = 1,size(this%invariantsLSTM,2)-1
           call readInvariants(timelist(n),this%runID,this%datadir,this%gpC,this%invariantsLSTM(1,n,:,:,:,:))
         end do
      end if
      call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%invariantsLSTM(1,10,:,:,:,:))
  end if

  if (this%NNtype == 1) then
    train_loadin = reshape(this%invariants,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*7])
    deltain = reshape(this%delta,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*8])
    strain2in = reshape(this%strain2,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strain2_min = reshape(this%strain2_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rot2in = reshape(this%rot2,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strainrot_min = reshape(this%strainrot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rotstrainrot_min = reshape(this%rotstrainrot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strain2rot_min = reshape(this%strain2rot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rotstrain2_min = reshape(this%rotstrain2_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    c_p = loadtensorflow(deltain, train_loadin, strain2in, strain2_min, rot2in, strainrot_min, &
        rotstrainrot_min, strain2rot_min, rotstrain2_min)
    call c_f_pointer(c_p, f_p, [1, 6, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)])
    f_p = reshape(f_p, [1, 6, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)], order=[5, 4, 3, 2, 1])
  end if


  if (this%NNtype == 2) then
    train_loadin = reshape(this%invariantsLSTM,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*7*10])
    deltain = reshape(this%delta,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*8])
    strain2in = reshape(this%strain2,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strain2_min = reshape(this%strain2_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rot2in = reshape(this%rot2,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strainrot_min = reshape(this%strainrot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rotstrainrot_min = reshape(this%rotstrainrot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strain2rot_min = reshape(this%strain2rot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rotstrain2_min = reshape(this%rotstrain2_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    c_p = loadtensorflow(deltain, train_loadin, strain2in, strain2_min, rot2in, strainrot_min, &
        rotstrainrot_min, strain2rot_min, rotstrain2_min)
    call c_f_pointer(c_p, f_p, [1, 6, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)])
    f_p = reshape(f_p, [1, 6, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)], order=[5, 4, 3, 2, 1])
  end if

end subroutine
