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
     allocate(this%delta(1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),8))
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
  class(sgs_igrid), intent(inout) :: this
  integer, intent(in) :: tidNow
  integer :: tidPast, n

  ! TODO: Andy, compute your invariants, call your NN, compute tauij

  if (this%NNtype == 2) then
      do n = 1,size(this%invariantsLSTM,2)-1
          ! TODO: Andy, compute tidPast based on simulation tid
          call readInvariants(tidPast,this%runID,this%datadir,this%gpC,this%invariantsLSTM(1,n,:,:,:,:))
      end do
      call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%invariantsLSTM(1,10,:,:,:,:))
  end if
end subroutine
