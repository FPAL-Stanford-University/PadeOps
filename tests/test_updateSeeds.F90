program test_updateSeeds
  use kind_parameters, only: rkind
  use GaborModeRoutines, only: initializeSeeds, updateSeeds, getIdx
  implicit none

  real(rkind), dimension(7) :: seed
  integer :: nx, ny, nz
  integer :: ist, jst, kst, i, j, k
  integer :: idxOld, idxCurrent
  character(len=2) :: istStr, jstStr, kstStr
  
  call GETARG(1,istStr)
  call GETARG(2,jstStr)
  call GETARG(3,kstStr)

  read(istStr,'(I2)') ist
  read(jstStr,'(I2)') jst
  read(kstStr,'(I2)') kst
 
  nx = 32
  ny = 32
  nz = 32
   
  call initializeSeeds(seed,ist,jst,kst,nx,ny)

  idxOld = getIdx(ist,jst,kst,nx,ny) - 1
  do k = kst,nz
    do j = jst,ny
      do i = ist,nx
        idxCurrent = getIdx(i,j,k,nx,ny)
        call updateSeeds(seed,idxOld,idxCurrent)
        idxOld = idxCurrent

        if (i == 13 .and. j == 7 .and. k == 24) then
          print*, nint(seed(1))
        end if
      end do
    end do
  end do  

end program test_updateSeeds
