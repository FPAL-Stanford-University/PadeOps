program test_QHmesh
  use kind_parameters,    only: rkind, clen
  use QHmeshMod,          only: QHmesh
  use IncompressibleGrid, only: igrid
  use mpi
  use decomp_2d,          only: nproc, nrank

  type(igrid)  :: LES
  type(QHmesh) :: QHgrid
  character(len=clen) :: inputLES, inputQH
  integer :: ierr, k, n

  call MPI_Init(ierr)

  inputLES = '/work2/06632/ryanhass/stampede2/problems/Enrichment/tests/'&
    //'QHmesh/inputLS.dat'
  inputQH  = '/work2/06632/ryanhass/stampede2/problems/Enrichment/tests/'&
    //'QHmesh/inputQHmesh.dat'

  call LES%init(inputLES,.true.)
  call QHgrid%init(inputQH,LES)

  do n = 1,nproc
    if (nrank == n-1) then
      print*,"===================================================="
      print*, " "
      print*, "PE", nrank
      do k = 1,QHgrid%gpC%xsz(3)
        print*, "zmin:", QHgrid%zE(k)
      end do
      print*, "zmin:", QHgrid%zE(QHgrid%gpC%xsz(3)+1)

      print*, "------------------------------------------------------"
      print*, " "
      
      do k = 1,QHgrid%gpC%xsz(2)
        print*, "ymin:", QHgrid%yE(k)
      end do
      print*, "ymin:", QHgrid%yE(QHgrid%gpC%xsz(2)+1)
      print*, "------------------------------------------------------"
      print*, " "
      
      do k = 1,QHgrid%gpC%xsz(1)
        print*, "xmin:", QHgrid%xE(k)
      end do
      print*, "xmin:", QHgrid%xE(QHgrid%gpC%xsz(1)+1)
    end if

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end do

  call MPI_Finalize(ierr)
end program

subroutine meshgen_WallM(decomp, dx, dy, dz, mesh, inputfile)
    use kind_parameters,  only: rkind, clen
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    integer :: nxg, nyg, nzg
    real(rkind)  :: Lx = two*pi, Ly = two*pi, Lz = two*pi
    namelist /HIT_PeriodicINPUT/ Lx, Ly, Lz 

    !Lx = two*pi; Ly = two*pi; Lz = one
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_PeriodicINPUT)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = two*pi
    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nxg,rkind)
        dy = Ly/real(nyg,rkind)
        dz = Lz/real(nzg,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
                end do
            end do
        end do

        ! Shift everything to the origin 
        x = x - dx
        y = y - dy
        z = z - dz 

    end associate

end subroutine 

subroutine initfields_wallM(decompC, decompE, inpDirectory, mesh, fieldsC, fieldsE)
    use decomp_2d, only: decomp_info
    use kind_parameters, only : rkind
    type(decomp_info), intent(in) :: decompC
    type(decomp_info), intent(in) :: decompE
    character(len=*), intent(in) :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsE

end subroutine 

subroutine hook_source(tsim,mesh,Re,urhs,vrhs,wrhs)
    use kind_parameters, only : rkind
    real(rkind),                     intent(in)    :: tsim, Re
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:),   intent(inout) :: urhs, vrhs, wrhs
end subroutine
        
subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters, only : rkind
    character(len=*), intent(in) :: inputfile
    real(rkind), intent(out) :: Tref

end subroutine 

subroutine setInhomogeneousNeumannBC_Temp(inpDirectory, wTh_surf)
    use kind_parameters, only: rkind
    character(len=*), intent(in) :: inpDirectory
    real(rkind), intent(inout) :: wTh_surf

end subroutine 

subroutine setDirichletBC_Temp(inpDirectory, Tsurf, dTsurfdt)
    use kind_parameters, only: rkind
    character(len=*), intent(in) :: inpDirectory
    real(rkind), intent(out) :: Tsurf, dTsurfdt

end subroutine 
subroutine set_planes_io(xplanes, yplanes, zplanes)
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters, only: rkind
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs

end subroutine

subroutine set_KS_planes_io(planesCourseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCourseGrid
end subroutine
subroutine initScalar(decompC, inpDirectory, mesh, scalar_id, scalarField)
    use decomp_2d, only: decomp_info
    use kind_parameters, only : rkind
    type(decomp_info), intent(in) :: decompC
    character(len=*), intent(in) :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:), intent(out) :: scalarField
    integer, intent(in) :: scalar_id
end subroutine 
subroutine SetScalar_source(decompC, inpDirectory, mesh, scalar_id, scalarSource)
    use decomp_2d, only: decomp_info
    use kind_parameters, only : rkind
    type(decomp_info), intent(in) :: decompC
    character(len=*), intent(in) :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:), intent(out) :: scalarSource
    integer, intent(in) :: scalar_id
end subroutine 
