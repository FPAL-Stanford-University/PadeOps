#include "hitModules/hitCD10_RHS.F90"
#include "hitModules/hitCD10_timestep.F90"
#include "hitModules/hitCD10_IO.F90"

program hitCD10

    use kind_parameters,  only: rkind,clen
    use constants,        only: two,pi
    use exits,            only: message,newline
    use reductions,       only: P_MAXVAL
    use decomp_2d,        only: decomp_2d_init, decomp_2d_finalize, get_decomp_info, decomp_info
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters
    use hitCD10_timestep, only: timestep_x_to_z, timestep_z_to_x, get_divergence_x
    use hitCD10_IO,       only: GetHit3d_uvw

    implicit none

    type( decomp_info ) :: gp
    type( derivatives ) :: der
    type( filters     ) :: fil

    real(rkind), dimension(:,:,:,:), allocatable, target :: xfields,yfields,zfields,xRHS,yRHS,zRHS
    real(rkind), dimension(:,:,:), pointer               :: ux,vx,wx,uy,vy,wy,uz,vz,wz


    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow=0, pcol=0
    integer :: ierr

    real(rkind) :: dx,dy,dz

    real(rkind) :: Re = 50._rkind
    character(len=*), parameter :: dtype = "cd10"
    character(len=*), parameter :: ftype = "cf90"
    character(len=clen) :: dir = "/afs/ir/users/a/k/akshays/hit3d_data/"

    ! Initialize MPI and 2DECOMP
    call MPI_INIT(ierr)
    call decomp_2d_init(nx,ny,nz,prow,pcol)
    call get_decomp_info(gp)

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    ! Initialize the derivatives objects
    call der%init(    gp,               &
                      dx,    dy,    dz, &
                  .TRUE.,.TRUE.,.TRUE., &
                   dtype, dtype, dtype  )

    ! Initialize the filters objects
    call fil%init(    gp,               &
                  .TRUE.,.TRUE.,.TRUE., &
                   ftype, ftype, ftype  )

    ! Allocate arrays
    allocate(  xfields( gp%xsz(1), gp%xsz(2), gp%xsz(3), 3 ) )
    allocate(  xRHS   ( gp%xsz(1), gp%xsz(2), gp%xsz(3), 3 ) )
    
    allocate(  yfields( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate(  yRHS   ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    
    allocate(  zfields( gp%zsz(1), gp%zsz(2), gp%zsz(3), 3 ) )
    allocate(  zRHS   ( gp%zsz(1), gp%zsz(2), gp%zsz(3), 3 ) )

    ux => xfields(:,:,:,1); vx => xfields(:,:,:,2); wx => xfields(:,:,:,3)
    uy => yfields(:,:,:,1); vy => yfields(:,:,:,2); wy => yfields(:,:,:,3)
    uz => zfields(:,:,:,1); vz => zfields(:,:,:,2); wz => zfields(:,:,:,3)

    ! Initialize the fields in the X decomposition
    call getHit3d_uvw(nx,ny,nz,xfields,gp,dir)
    
    ! Check if initial condition is divergence free
    call get_divergence_x(ux,vx,wx,vy,wy,wz,der,gp,uz)

    call newline()
    call message(1,"Maximum absolute divergence",P_MAXVAL(ABS(uz)))
  
    call timestep_x_to_z(ux,vx,wx,uy,vy,wy,uz,vz,wz,xRHS,yRHS,zRHS,Re,der,gp)
    call timestep_z_to_x(ux,vx,wx,uy,vy,wy,uz,vz,wz,xRHS,yRHS,zRHS,Re,der,gp)
    
    ! Dellocate arrays
    deallocate(  xfields )
    deallocate(  xRHS    )
      
    deallocate(  yfields )
    deallocate(  yRHS    )
      
    deallocate(  zfields )
    deallocate(  zRHS    )

    call decomp_2d_finalize()
    call MPI_Finalize(ierr)

end program
