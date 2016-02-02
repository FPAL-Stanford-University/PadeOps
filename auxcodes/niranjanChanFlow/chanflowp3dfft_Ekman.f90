MODULE SUBROUTINES

  !USE mpi
  USE p3dfft
  IMPLICIT NONE
  INCLUDE "fftw3.f"
  INCLUDE "mpif.h"

  INTEGER, PARAMETER :: DP = KIND(1.0d0)
  COMPLEX(DP), PARAMETER :: iu = (0.0d0,1.0d0)
  REAL(DP), PARAMETER :: zero = 0.0d0
  REAL(DP), PARAMETER :: two = 2.0d0
  REAL(DP), PARAMETER :: twothird = 2.0d0/3.0d0
  REAL(DP), PARAMETER :: third = 1.0d0/3.0d0
  REAL(DP), PARAMETER :: ninth = 1.0d0/9.0d0
  REAL(DP), PARAMETER :: sixth = 1.0d0/6.0d0
  REAL(DP), PARAMETER :: half = 1.0d0/2.0d0
  REAL(DP), PARAMETER :: tseventh = 1.0d0/27.0d0
  REAL(DP), PARAMETER :: pi = 3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: twopi = 2.0d0*pi
  REAL(DP), PARAMETER :: piby3 = pi/3.0d0

  REAL(DP) :: Re, Pr, Ri, dt, tmax, aa, bb, tstatstart, tsnap, t, tsc, Rmax, Rmin
  REAL(DP) :: utaul1, utaul2, utau1, utau2
  INTEGER :: nxg, nyg, nzg, istretch, irestart, nprocz,nprocy, nvar, iitt, iwpov, nxgrst, nygrst, nzgrst

  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: tausgs, tautemp, lterm, mterm
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: u, dudx, dudy, dudz, conv, eddyvisc, strainmag, utemp, uxtemp
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: p, wp, rhow, ubufzm1, ubufzp1, secondEigv, sendbuf, recvbuf, dwpdx, dwpdy, dwpdz, tmp, divu, rhsp, dpdx, dpdy, dpdz, uuavg, savg, uin, hp, hprho, numg, deng, numrhog, denrhog, prodg, prodrhog,tauavg,spec!, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, drhodx, drhody, drhodz
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ubar, vbar, uavg, tavg, sendbuffilt, recvbuffilt
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: h, dz, delz2, deltabar2, x, y, z, zetaz, zetazz, zz, zzetaz, zzetazz, amg, bmg, cmg, amt, bmt, cmt, amn, bmn, cmn, am, bm, cm, amr, bmr, cmr, cp, senddiag, recvdiag, ssavg_d, ttavg, kg, omg, nutg, coef, coefrho, num, den, numrho, denrho, numl, denl, numrhol, denrhol, vol, epssgs, epssgs_rho, ssavg_rhod, tmp1d

  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: uh, sgstermh, sgstemph, viscnh, viscth, rhsh, rhs1h, convh, d2udzh, buoyh, uhtemp, fCorriolh
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:,:) :: ph, wph, dudzh, dvdzh, dwdzh, drhodzh, rhsph, dpdzh, tmph, dt_pgradh, tprofh, rhowh, uprofh, vprofh, onetr
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: dph, dsgsdzh, rhomhl, rhomh, tmp1dh

  !INTEGER, ALLOCATABLE, DIMENSION(:) :: k1, k2
  !INTEGER :: Nxdealias, Nydealias
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: k1, k2 ! REAL so as to handle domain sizes which are non-integral multiples of twopi
  REAL(DP) :: Nxdealias, Nydealias ! REAL so as to handle domain sizes which are non-integral multiples of twopi

  INTEGER :: MPIRP, MPICP, MPICOMM2D, MPICOORDX, MPICOORDY, MPICOORDZ, MYRANK, MPIFAIL, MPISRC, MPIDEST, mpistat(MPI_STATUS_SIZE), outid
  INTEGER :: nx, ny, nz, sx, sy, sz, ex, ey, ez, nxf, nyf, nzf, sxf, syf, szf, exf, eyf, ezf, itsnap, itstat, statindex, ievm, nxc, nyc, nzc
  REAL(DP) :: itfac, a, b, c, d, e, f, lx, ly, lz, dx, dy, delx2, dely2, pgr, maxdivu, mindivu, maxcfl, maxu1, maxu2, maxu3, maxu4, minu1, minu2, minu3, minu4, cs, Prt, mindz, mfr
  CHARACTER(LEN=3) :: fstr, bstr

  INTEGER*8 :: plan_f

  ! Sigma Model
  REAL(DP) :: g(3,3),gt(3,3),gtg(3,3),isgma(3),sigma1,sigma2,sigma3,sigma(3),dummy(3)

  ! Frequency of application of dynamic procedure
  INTEGER :: idyn, DYNFREQ
  REAL(DP) :: globcoef, globcoefrho

  CONTAINS

  SUBROUTINE startmpi

    IMPLICIT NONE
    CHARACTER(*), PARAMETER :: subname='startmpi'
    INTEGER :: mpierr, numprocs, ndims, dims(0:1), direc, displ
    LOGICAL :: periods(0:1), reorder, coords(0:1)

    call MPI_INIT(mpierr) ! Initialize the MPI environment

    if(mpierr /= MPI_SUCCESS) then
      print *,'MPI initialization error'
      stop
    endif

    MPIRP = MPI_DOUBLE_PRECISION
    MPICP = MPI_DOUBLE_COMPLEX
      
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpierr)
      
    if(numprocs /= nprocz*nprocy) then
      print*,'You must start MPI with ',nprocz*nprocy,'processors'
      return
    end if

    ! Set up a 2D Cartesian topology
    ndims = 2
    dims(0) = nprocz; dims(1) = nprocy
    periods(0) = .FALSE.; periods(1) = .TRUE.
    reorder = .TRUE.

    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods,reorder, MPICOMM2D, mpierr)
    !call checkMPI(mpierr, 'MPI_CART_CREATE', subname, .true.)
    call MPI_ERRHANDLER_SET(MPICOMM2D, MPI_ERRORS_RETURN, mpierr)
    !call checkMPI(mpierr, 'MPI_ERRHANDLER_SET', subname, .true.)

    ! Query where I am in the Cartesian topology 
    call MPI_COMM_RANK(MPICOMM2D, MYRANK, mpierr)
    !call checkMPI(mpierr, 'MPI_COMM_RANK', subname, .true.)

    call MPI_CART_COORDS(MPICOMM2D, MYRANK, ndims, coords, mpierr)
    !call checkMPI(mpierr, 'MPI_CART_COORDS', subname, .true.)
    MPICOORDX = 0; MPICOORDY = coords(0); MPICOORDZ = coords(1)

    direc = 0; displ = 1
    call MPI_CART_SHIFT(MPICOMM2D, direc, displ, MPISRC, MPIDEST, mpierr)
    !11111write(*,*) mpierr

    write(*,1), 'proc', MYRANK, MPICOORDX, MPICOORDY, MPICOORDZ, MPISRC, MPIDEST

1 FORMAT(a,7(1x,i4))

  END SUBROUTINE startmpi

  SUBROUTINE stopmpi

    IMPLICIT NONE
    INTEGER :: mpierr

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    call MPI_FINALIZE(mpierr)

  END SUBROUTINE stopmpi

  SUBROUTINE startp3dfft

    IMPLICIT NONE
    INTEGER :: dims(2), isize(3), istart(3), iend(3), fsize(3), fstart(3), fend(3)
    CHARACTER(LEN=13) :: fname

    nxc = nx; nyc = ny; nzc = nz

    dims(1) = nprocy; dims(2) = nprocz
    call p3dfft_setup(dims,nxg,nyg,nzg,MPICOMM2D)
    !call p3dfft_setup(dims,nxg,nyg,nzg,MPICOMM2D,nxc,nyc,nzc,.TRUE.)
    call p3dfft_get_dims(istart,iend,isize,1)
    call p3dfft_get_dims(fstart,fend,fsize,2)

    nx = isize(1); ny = isize(2); nz = isize(3)
    sx = istart(1); sy = istart(2); sz = istart(3)
    ex = iend(1); ey = iend(2); ez = iend(3)

    nxf = fsize(1);  nyf = fsize(2);  nzf = fsize(3)
    sxf = fstart(1); syf = fstart(2); szf = fstart(3)
    exf = fend(1);   eyf = fend(2);   ezf = fend(3)

    itfac = 1.0d0/dble(nxg*nyg)
    fstr = 'ffn'; bstr = 'nff'

    WRITE(fname,1) 'proc_',MYRANK,'.dat'
    OPEN(11,FILE=fname,STATUS='replace')
    WRITE(11,2) MYRANK, MPICOORDX, MPICOORDY, MPICOORDZ
    WRITE(11,2) nx, ny, nz
    WRITE(11,2) sx, sy, sz
    WRITE(11,2) ex, ey, ez
    WRITE(11,2) nxf, nyf, nzf
    WRITE(11,2) sxf, syf, szf
    WRITE(11,2) exf, eyf, ezf
    CLOSE(11)


1 FORMAT(a,i4.4,a)
2 FORMAT(4(1x,i4.4))

  END SUBROUTINE startp3dfft

  SUBROUTINE stopp3dfft

    IMPLICIT NONE

    CALL dfftw_destroy_plan(plan_f)

  END SUBROUTINE stopp3dfft

  SUBROUTINE memvars

    IMPLICIT NONE

    ALLOCATE(k1(1:nxg/2+1),k2(1:nyg))
    ALLOCATE(uh(sxf:exf,syf:eyf,szf:ezf,nvar),ph(sxf:exf,syf:eyf,szf:ezf),wph(sxf:exf,syf:eyf,szf:ezf))
    ALLOCATE(d2udzh(sxf:exf,syf:eyf,szf:ezf,nvar),dudzh(sxf:exf,syf:eyf,szf:ezf))
    ALLOCATE(dwdzh(sxf:exf,syf:eyf,szf:ezf),dvdzh(sxf:exf,syf:eyf,szf:ezf))
    ALLOCATE(drhodzh(sxf:exf,syf:eyf,szf:ezf))!,drhodzh(sxf:exf,syf:eyf,szf:ezf))
    ALLOCATE(rhsph(sxf:exf,syf:eyf,szf:ezf),dpdzh(sxf:exf,syf:eyf,szf:ezf))
    ALLOCATE(tmph(sxf:exf,syf:eyf,szf:ezf),dt_pgradh(sxf:exf,syf:eyf,szf:ezf),tprofh(sxf:exf,syf:eyf,2),uprofh(sxf:exf,syf:eyf,2),vprofh(sxf:exf,syf:eyf,2))

    ALLOCATE(u(sx:ex,sy:ey,sz:ez,nvar),conv(sx:ex,sy:ey,sz:ez,nvar),wp(sx:ex,sy:ey,sz:ez),p(sx:ex,sy:ey,sz:ez))
    ALLOCATE(dudx(sx:ex,sy:ey,sz:ez,nvar),dudy(sx:ex,sy:ey,sz:ez,nvar),dudz(sx:ex,sy:ey,sz:ez,nvar))
    ALLOCATE(dwpdx(sx:ex,sy:ey,sz:ez),dwpdy(sx:ex,sy:ey,sz:ez),dwpdz(sx:ex,sy:ey,sz:ez))
    ALLOCATE(ubufzm1(sx:ex,sy:ey,nvar+1),ubufzp1(sx:ex,sy:ey,nvar+1))
    ALLOCATE(sendbuf(sx:ex,sy:ey,nvar+1),recvbuf(sx:ex,sy:ey,nvar+1))
    ALLOCATE(divu(sx:ex,sy:ey,sz:ez))
    !ALLOCATE(dvdx(sx:ex,sy:ey,sz:ez),dvdy(sx:ex,sy:ey,sz:ez),dvdz(sx:ex,sy:ey,sz:ez))
    !ALLOCATE(dwdx(sx:ex,sy:ey,sz:ez),dwdy(sx:ex,sy:ey,sz:ez),dwdz(sx:ex,sy:ey,sz:ez))
    !ALLOCATE(drhodx(sx:ex,sy:ey,sz:ez),drhody(sx:ex,sy:ey,sz:ez),drhodz(sx:ex,sy:ey,sz:ez))
    ALLOCATE(x(nxg),y(nyg),z(nzg))
    ALLOCATE(zetaz(nzg),zetazz(nzg))
    ALLOCATE(zz(nzg-1),zzetaz(nzg-1),zzetazz(nzg-1))
    ALLOCATE(h(3),dz(nzg))
    ALLOCATE(amg(nzg),bmg(nzg),cmg(nzg))
    ALLOCATE(amt(nzg),bmt(nzg),cmt(nzg))
    ALLOCATE(amn(nzg),bmn(nzg),cmn(nzg))
    ALLOCATE(amr(nzg),bmr(nzg),cmr(nzg))
    ALLOCATE(am(nzg),bm(nzg),cm(nzg),cp(nzg),dph(nzg))
    ALLOCATE(senddiag(12),recvdiag(12))
    ALLOCATE(rhow(sx:ex,sy:ey,sz:ez))
    ALLOCATE(deltabar2(nzg),delz2(nzg),kg(nzg),omg(nzg),nutg(nzg),vol(nzg))
    ALLOCATE(eddyvisc(sx:ex,sy:ey,sz:ez,nvar-2),tausgs(sx:ex,sy:ey,sz:ez,nvar,3))
    ALLOCATE(strainmag(sx:ex,sy:ey,sz:ez,nvar-2))
    ALLOCATE(sgstermh(sxf:exf,syf:eyf,szf:ezf,nvar),sgstemph(sxf:exf,syf:eyf,szf:ezf,nvar),dsgsdzh(szf:ezf))
    ALLOCATE(viscnh(sxf:exf,syf:eyf,szf:ezf,nvar),viscth(sxf:exf,syf:eyf,szf:ezf,nvar))
    ALLOCATE(rhsh(sxf:exf,syf:eyf,szf:ezf,nvar),rhs1h(sxf:exf,syf:eyf,szf:ezf,nvar))
    ALLOCATE(convh(sxf:exf,syf:eyf,szf:ezf,nvar))
    ALLOCATE(buoyh(sxf:exf,syf:eyf,szf:ezf,nvar),rhowh(sxf:exf,syf:eyf,szf:ezf),rhomhl(szf:ezf),rhomh(szf:ezf))
    ALLOCATE(secondEigv(sx:ex,sy:ey,sz:ez),ubar(sx:ex,sy:ey),vbar(sx:ex,sy:ey))
    ALLOCATE(uxtemp(sx:ex,sy:ey,sz:ez,3))
    ALLOCATE(fCorriolh(sxf:exf,syf:eyf,szf:ezf,nvar),onetr(sxf:exf,syf:eyf,szf:ezf))

    !! dynamic gradmodel
    !ALLOCATE(sendbuffilt(sx:ex,sy:ey),recvbuffilt(sx:ex,sy:ey))
    !ALLOCATE(uin(sx:ex,sy:ey,sz:ez),hp(sx:ex,sy:ey,sz:ez),hprho(sx:ex,sy:ey,sz:ez))
    !ALLOCATE(lterm(sx:ex,sy:ey,sz:ez,3,3),mterm(sx:ex,sy:ey,sz:ez,3,3))
    !ALLOCATE(tautemp(sx:ex,sy:ey,sz:ez,3,3),utemp(sx:ex,sy:ey,sz:ez,nvar),uhtemp(sxf:exf,syf:eyf,szf:ezf,nvar))
    !ALLOCATE(coef(sz:ez),num(sz:ez),den(sz:ez),numl(sz:ez),denl(sz:ez))
    !ALLOCATE(coefrho(sz:ez),numrho(sz:ez),denrho(sz:ez),numrhol(sz:ez),denrhol(sz:ez))

    ! global dynamic gradmodel
    ALLOCATE(sendbuffilt(sx:ex,sy:ey),recvbuffilt(sx:ex,sy:ey))
    ALLOCATE(uin(sx:ex,sy:ey,sz:ez),hp(sx:ex,sy:ey,sz:ez),hprho(sx:ex,sy:ey,sz:ez))
    ALLOCATE(numg(sx:ex,sy:ey,sz:ez),deng(sx:ex,sy:ey,sz:ez),numrhog(sx:ex,sy:ey,sz:ez),denrhog(sx:ex,sy:ey,sz:ez))
    ALLOCATE(prodg(sx:ex,sy:ey,sz:ez),prodrhog(sx:ex,sy:ey,sz:ez))
    ALLOCATE(tautemp(sx:ex,sy:ey,sz:ez,3,3),utemp(sx:ex,sy:ey,sz:ez,nvar),uhtemp(sxf:exf,syf:eyf,szf:ezf,nvar))

    !! debug
    ALLOCATE(rhsp(sx:ex,sy:ey,sz:ez),dpdx(sx:ex,sy:ey,sz:ez),dpdy(sx:ex,sy:ey,sz:ez),dpdz(sx:ex,sy:ey,sz:ez))

    ! stats
    ALLOCATE(uavg(sz:ez,nvar),uuavg(sz:ez,3,3),ttavg(sz:ez),ssavg_d(sz:ez),savg(sz:ez,3,3),tavg(sz:ez,3))
    ALLOCATE(epssgs(sz:ez),epssgs_rho(sz:ez),tauavg(sz:ez,3,3),ssavg_rhod(sz:ez))
    ALLOCATE(tmp1d(1:nxg),tmp1dh(1:nxg/2+1),spec(1:nxg/2+1,sz:ez,4))

    ! init fftw for spectrum statistics
    call dfftw_plan_dft_r2c_1d(plan_f,nxg,tmp1d,tmp1dh,FFTW_ESTIMATE)

  END SUBROUTINE memvars

  SUBROUTINE memvarsdeal

    IMPLICIT NONE

    DEALLOCATE(k1,k2)
    DEALLOCATE(uh,ph,wph)
    DEALLOCATE(d2udzh,dudzh)
    DEALLOCATE(dwdzh,dvdzh)
    DEALLOCATE(drhodzh)!,drhodzh(sxf:exf,syf:eyf,szf:ezf))
    DEALLOCATE(rhsph,dpdzh)
    DEALLOCATE(tmph,dt_pgradh,tprofh,uprofh,vprofh)

    DEALLOCATE(u,conv,wp,p)
    DEALLOCATE(dudx,dudy,dudz)
    DEALLOCATE(dwpdx,dwpdy,dwpdz)
    DEALLOCATE(ubufzm1,ubufzp1)
    DEALLOCATE(sendbuf,recvbuf)
    DEALLOCATE(divu)
    !DEALLOCATE(dvdx(sx:ex,sy:ey,sz:ez),dvdy(sx:ex,sy:ey,sz:ez),dvdz(sx:ex,sy:ey,sz:ez))
    !DEALLOCATE(dwdx(sx:ex,sy:ey,sz:ez),dwdy(sx:ex,sy:ey,sz:ez),dwdz(sx:ex,sy:ey,sz:ez))
    !DEALLOCATE(drhodx(sx:ex,sy:ey,sz:ez),drhody(sx:ex,sy:ey,sz:ez),drhodz(sx:ex,sy:ey,sz:ez))
    DEALLOCATE(x,y,z)
    DEALLOCATE(zetaz,zetazz)
    DEALLOCATE(zz,zzetaz,zzetazz)
    DEALLOCATE(h,dz)
    DEALLOCATE(amg,bmg,cmg)
    DEALLOCATE(amt,bmt,cmt)
    DEALLOCATE(amn,bmn,cmn)
    DEALLOCATE(amr,bmr,cmr)
    DEALLOCATE(am,bm,cm,cp,dph)
    DEALLOCATE(senddiag,recvdiag)
    DEALLOCATE(rhow)
    DEALLOCATE(deltabar2,delz2,kg,nutg,omg,vol)
    DEALLOCATE(eddyvisc,tausgs)
    DEALLOCATE(strainmag)
    DEALLOCATE(sgstermh,sgstemph,dsgsdzh)
    DEALLOCATE(viscnh,viscth)
    DEALLOCATE(rhsh,rhs1h)
    DEALLOCATE(convh)
    DEALLOCATE(buoyh,rhowh,rhomh,rhomhl)
    DEALLOCATE(secondEigv,ubar,vbar)
    DEALLOCATE(uxtemp)
    DEALLOCATE(fCorriolh,onetr)

    !! dynamic gradmodel
    !DEALLOCATE(sendbuffilt,recvbuffilt)
    !DEALLOCATE(uin,hp,hprho)
    !DEALLOCATE(lterm,mterm)
    !DEALLOCATE(tautemp,utemp,uhtemp)
    !DEALLOCATE(coef,num,den,numl,denl)
    !DEALLOCATE(coefrho,numrho,denrho,numrhol,denrhol)

    ! global dynamic gradmodel
    DEALLOCATE(sendbuffilt,recvbuffilt)
    DEALLOCATE(uin,hp,hprho)
    DEALLOCATE(numg,deng,numrhog,denrhog)
    DEALLOCATE(prodg,prodrhog)
    DEALLOCATE(tautemp,utemp,uhtemp)

    !! debug
    DEALLOCATE(rhsp,dpdx,dpdy,dpdz)

    ! stats
    DEALLOCATE(uavg,uuavg,ttavg,ssavg_d,savg,tavg,tauavg,epssgs,epssgs_rho,ssavg_rhod)
    DEALLOCATE(tmp1d,tmp1dh,spec)

  END SUBROUTINE memvarsdeal

  SUBROUTINE arrayinit

    IMPLICIT NONE

    statindex = 0; itstat = -1
    uavg = 0.0d0; uuavg = 0.0d0; ttavg = 0.0d0; ssavg_d = 0.0d0
    savg = 0.0d0; tavg = 0.0d0
    tauavg = 0.0d0; epssgs = 0.0d0; epssgs_rho = 0.0d0; ssavg_rhod = 0.0d0
    spec = 0.0d0
    buoyh = (0.0d0,0.0d0)
    fCorriolh = (0.0d0,0.0d0)

  END SUBROUTINE arrayinit

  SUBROUTINE grid

    IMPLICIT NONE
    INTEGER :: i,j,k

    a = 0.0d0; b = twopi
    c = 0.0d0; d = 0.75d0*pi
    e = 0.0d0; f = 1.0d0

    lx = b-a; h(1) = lx/dble(nxg-1)
    do i = 1, nxg
      x(i) = a + h(1)*dble(i-1)
    enddo
    dx = h(1)

    ly = d-c; h(2) = ly/dble(nyg-1)
    do j = 1, nyg
      y(j) = c + h(2)*dble(j-1)
    enddo
    dy = h(2)

    lz = f-e; h(3) = lz/dble(nzg-1)
    do k = 1, nzg
      z(k) = e + h(3)*dble(k-1)
    enddo
    do k = 1, nzg-1
      zz(k) = 0.5d0*(z(k)+z(k+1))
    enddo
    dz = z(2) - z(1)

    do i = 1, nxg/2+1
      k1(i) = dble(i-1)
    enddo

    do j = 1, nyg/2+1
      k2(j) = dble(j-1)
    enddo
    do j = nyg/2+2, nyg
      k2(j) = dble(j-1-nyg)
    enddo

    Nxdealias = dble(nxg/3); Nydealias = dble(nyg/3)

    ! domains non-integral multiples of twopi
    k1 = k1*twopi/lx; k2 = k2*twopi/ly
    Nxdealias = Nxdealias*twopi/lx; Nydealias = Nydealias*twopi/ly

  END SUBROUTINE grid

  SUBROUTINE metrics

    IMPLICIT NONE
    INTEGER :: k
    REAL(DP) :: bterm, pow, numer, denom, arg
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: zeta

    IF(istretch==0) then
      zetaz = 1.0d0; zzetaz = 1.0d0
      zetazz = 0.0d0; zzetazz = 0.0d0
      mindz = dz(1)
    ELSE
      ALLOCATE(zeta(nzg))
      bterm = (bb+1.0d0)/(bb-1.0d0)
      do k = 1, nzg
        zeta(k) = (z(nzg-k+1)-e)/(f-e)
      enddo
      do k = 1, nzg
        ! Transformation
        pow = zeta(k)
        numer = bb*(bterm**pow - 1.0d0)
        denom = 1.0d0+bterm**pow
        z(k) = lz*(1.0d0-numer/denom) + e
        ! Metric derivatives
        numer = 2.0d0*bb
        arg = - bb**2 + (1.0d0 - (z(k)-e)/lz)**2
        denom = lz*arg*dlog(bterm)
        zetaz(k) = numer/denom
        zetazz(k) = numer/(denom*arg)* &
                   2.0d0*(1.0d0-(z(k)-e)/lz)/lz
      end do
      do k = 1,nzg-1
        pow = zeta(k)+(zeta(k+1)-zeta(k))*0.5d0
        numer = bb*(bterm**pow - 1.0d0)
        denom = 1.0d0+bterm**pow
        zz(k) = lz*(1.0d0-numer/denom) + e
        numer = 2.0*bb
        arg = -bb**2 + (1.0d0 - (zz(k)-e)/lz)**2
        denom = lz*arg*dlog(bterm)
        zzetaz(k) = numer/denom
        zzetazz(k) = numer/(denom*arg)* &
                   2.0d0*(1.0d0-(zz(k)-e)/lz)/lz
      end do

      ! new definition of dz
      dz(1) = (z(2) - z(1))
      do k = 2, nzg-1
        dz(k) = 0.5d0*(z(k+1) - z(k-1))
      end do
      dz(nzg) = (z(nzg) - z(nzg-1))
      mindz = minval(dz)

      !! one definition of dz
      !dz(1) = 2.0d0*(zz(1) - z(1))
      !do k = 2, nzg-1
      !  dz(k) = 0.5d0*(z(k) - z(k-1))
      !end do
      !dz(nzg) = 2.0d0*(z(nzg) - zz(nzg-1))

      DEALLOCATE(zeta) 
    ENDIF

    delx2 = dx**2; dely2 = dy**2
    do k = 1, nzg
      delz2(k) = dz(k)**2
      deltabar2(k) = (dx*dy*dz(k))**twothird
      vol(k) = dx*dy*dz(k)
    enddo

  END SUBROUTINE metrics

  SUBROUTINE rhsbc

    IMPLICIT NONE

    ! szf=1; ezf=nzg
    rhsh(:,:,szf,1:2) = 0.0d0
    !rhsh(:,:,ezf,1:2) = 0.0d0

    !rhsh(:,:,ezf,1) = uprofh(:,:,2)
    !rhsh(:,:,ezf,2) = vprofh(:,:,2)

    rhsh(:,:,szf,4) = tprofh(:,:,1)
    !rhsh(:,:,ezf,4) = tprofh(:,:,2)

  END SUBROUTINE rhsbc

  SUBROUTINE bch

    IMPLICIT NONE

    ! szf=1; ezf=nzg
    uh(:,:,szf,1:2) = 0.0d0

    !uh(:,:,ezf,1:2) = 0.0d0
    uh(:,:,ezf,1) = uh(:,:,ezf-1,1)
    uh(:,:,ezf,2) = uh(:,:,ezf-1,2)

  END SUBROUTINE bch

  SUBROUTINE bcp

  IF(sz==1) then
    wp(:,:,sz) = 0.0d0
  ENDIF
  IF(ez==nzg) then
    wp(:,:,ez) = 0.0d0
  ENDIF

  END SUBROUTINE bcp

  SUBROUTINE bc

  IF(sz==1) then
    u(:,:,sz,1:2) = 0.0d0
    u(:,:,sz,4) = 0.0d0
  ENDIF
  IF(ez==nzg) then
    u(:,:,ez,1) = u(:,:,ez-1,1)
    u(:,:,ez,2) = u(:,:,ez-1,2)
  !  u(:,:,ez,4) = -0.5d0
  ENDIF

  call bcp

  END SUBROUTINE bc

  SUBROUTINE dealiasuh

    IMPLICIT NONE
    INTEGER :: i, j, k, var

    do j = syf, eyf
     do i = sxf, exf
       if(k1(i) >= Nxdealias) then
         uh(i,j,:,:) = 0.0d0
         wph(i,j,:) = 0.0d0
       elseif(k2(j) >= Nydealias .or. k2(j) <= -Nydealias) then
         uh(i,j,:,:) = 0.0d0
         wph(i,j,:) = 0.0d0
       endif
     enddo
    enddo

  END SUBROUTINE dealiasuh

  SUBROUTINE init

    IMPLICIT NONE
    INTEGER :: i, j, k, var, rr, mpierr
    REAL(DP) :: epsfac, Uperiods, Vperiods, zpeak

      !do k = sz, ez
      !  write(100+MYRANK,*) k, z(k)
      !enddo

    !iitt = 0
    epsfac = 0.5d0
    Uperiods = 4.0d0; Vperiods = 4.0d0; zpeak = 0.3d0
    u = 0.0d0; p = 0.0d0
    do i = sx, ex
     do j = sy, ey
      do k = sz, ez
        !u(i,j,k,1) = (1.0d0-(z(k)-1.0d0)**8.0d0)+epsfac*0.5d0*lx*dcos(0.5d0*pi*(z(k)-1.0d0))*dcos(2.0d0*twopi/lx*x(i))*dsin(twopi/ly*y(j))
        !u(i,j,k,2) = -epsfac*0.5d0*ly*dcos(0.5d0*pi*(z(k)-1.0d0))*dsin(2.0d0*twopi/lx*x(i))*dcos(twopi/ly*y(j))
        !if(k<nzg) u(i,j,k,3) = -two*epsfac*(0.0d0-dsin(0.5d0*pi*(zz(k)-1.0d0)))*dsin(2.0d0*twopi/lx*x(i))*dsin(twopi/ly*y(j))

        !u(i,j,k,1) = (1.0d0-((z(k)-f)/Lz)**8.0d0)+epsfac*dexp(0.5d0)*(z(k)/Lz)*dcos(Uperiods*2.0d0**pi*y(j)/Ly)*dexp(-0.5d0*(z(k)/zpeak/Lz)**2.0d0)
        !u(i,j,k,2) =                              epsfac*dexp(0.5d0)*(z(k)/Lz)*dcos(Vperiods*2.0d0**pi*x(i)/Lx)*dexp(-0.5d0*(z(k)/zpeak/Lz)**2.0d0)

        u(i,j,k,1) = 1.0d0-dexp(-z(k)*20.0d0)*dcos(z(k)*20.0d0)+epsfac*dexp(0.5d0)*(z(k)/Lz)*dcos(Uperiods*2.0d0**pi*y(j)/Ly)*dexp(-0.5d0*(z(k)/zpeak/Lz)**2.0d0)
        u(i,j,k,2) =       dexp(-z(k)*20.0d0)*dsin(z(k)*20.0d0)+epsfac*dexp(0.5d0)*(z(k)/Lz)*dcos(Vperiods*2.0d0**pi*x(i)/Lx)*dexp(-0.5d0*(z(k)/zpeak/Lz)**2.0d0)

        u(i,j,k,4) = (0.0d0+z(k))!*0.5d0
        !if(i==13 .and. j==13) then
        !  write(*,*) k, z(k), dcos(0.5d0*pi*(z(k)-1.0d0))
        !endif
      enddo
     enddo
    enddo

    call massflowrate
    IF(MYRANK==outid) write(*,*) 'Mass Flow Rate = ', mfr

    ! set pressure gradient
    pgr = 9.812d0/Re
    ALLOCATE(tmp(sx:ex,sy:ey,sz:ez))
    tmp(:,:,:) = itfac*dt*pgr
    call p3dfft_ftran_r2c(tmp(:,:,:), dt_pgradh(:,:,:),fstr)
    DEALLOCATE(tmp)

    ! set onetr for Corriolis force
    ALLOCATE(tmp(sx:ex,sy:ey,sz:ez))
    tmp(:,:,:) = itfac*1.0d0
    call p3dfft_ftran_r2c(tmp(:,:,:), onetr(:,:,:),fstr)
    DEALLOCATE(tmp)

    ! set bounds for rho
    Rmax = 1.0d0
    Rmin = 0.0d0

    call bc

    do k = sz, ez
      write(myrank+100,*) k, z(k), u(12,12,k,1)
    enddo
    ! debug
    !do rr = 0, nprocz*nprocy-1
    !  IF(rr==myrank) then
    !    OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !    write(110,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1)) 
    !    write(110,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2)) 
    !    write(110,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3)) 
    !    write(110,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4)) 
    !    write(110,*) myrank
    !    CLOSE(110)
    !  ENDIF
    !  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    !enddo
    !IF(myrank==0) then
    !  OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !  write(110,*) '---'
    !  CLOSE(110)
    !ENDIF
    !write(myrank+100,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+100,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+100,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+100,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))

    do var = 1, nvar
      call p3dfft_ftran_r2c(u(:,:,:,var),uh(:,:,:,var),fstr)
    end do
    uh = itfac*uh
    ph = 0.0d0
    call vinterp
    call vinterph

    ! szf = 1; ezf = nzg
    uprofh(:,:,2) = uh(:,:,ezf,1)
    vprofh(:,:,2) = uh(:,:,ezf,2)

    tprofh(:,:,1) = uh(:,:,szf,4)
    tprofh(:,:,2) = uh(:,:,ezf,4)

    !! debug
    !call decode
    !do rr = 0, nprocz*nprocy-1
    !  IF(rr==myrank) then
    !    OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !    write(110,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1)) 
    !    write(110,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2)) 
    !    write(110,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3)) 
    !    write(110,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4)) 
    !    write(110,*) myrank
    !    CLOSE(110)
    !  ENDIF
    !  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    !enddo
    !IF(myrank==0) then
    !  OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !  write(110,*) '---'
    !  CLOSE(110)
    !ENDIF
    !write(myrank+100,*) '---'
    !write(myrank+100,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+100,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+100,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+100,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))

    call divergence
    IF(MYRANK==outid) then
      write(*,35) 'Init: ',t, maxdivu, maxu1, maxu2, maxcfl
    ENDIF

    call correct
    call bch
    call dealiasuh
    call decode
    call divergence
    IF(MYRANK==outid) then
      write(*,35) 'Init: ',t, maxdivu, maxu1, maxu2, maxcfl
    ENDIF

    !! debug
    !do rr = 0, nprocz*nprocy-1
    !  IF(rr==myrank) then
    !    OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !    write(110,'(i3.3,2(1x,e19.12))') rr,maxval(p(:,:,:)), minval(p(:,:,:)) 
    !    if(rr==nprocz*nprocy-1) write(110,*) '---'
    !    CLOSE(110)
    !  ENDIF
    !  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    !enddo

    do var = 1, nvar
      call p3dfft_ftran_r2c(u(:,:,:,var),uh(:,:,:,var),fstr)
    end do
    uh = itfac*uh
    call vinterp
    call vinterph
    call bch
    call dealiasuh
    call decode
    call divergence
    IF(MYRANK==outid) then
      write(*,35) 'Init: ',t, maxdivu, maxu1, maxu2, maxcfl
    ENDIF

    !! debug
    !do rr = 0, nprocz*nprocy-1
    !  IF(rr==myrank) then
    !    OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !    write(110,'(i3.3,2(1x,e19.12))') rr,maxval(p(:,:,:)), minval(p(:,:,:)) 
    !    if(rr==nprocz*nprocy-1) write(110,*) '---'
    !    CLOSE(110)
    !  ENDIF
    !  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    !enddo

    !call correct
    !call divergence
    !IF(MYRANK==outid) then
    !  write(*,35) 'Init: ',t, maxdivu, maxu1, maxu2, maxcfl
    !ENDIF
35 format(a,7(1x,e12.5),1x,i5,1x,i9)

  END SUBROUTINE init

  SUBROUTINE restartinit

    IMPLICIT NONE
    INTEGER :: var, k
    CHARACTER(LEN=15) :: fname
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: urst
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: prst
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: xrst, yrst, zrsy, xgrst, ygrst, zgrst
    INTEGER :: nxrst, nyrst, nzrst, sxrst, syrst, szrst, exrst, eyrst, ezrst


    IF(irestart==1) then
      ! irestart = 1 :: continue old run; same resolution as previous run. No
      ! interpolation needed. rhs1h needs to be read in.
      WRITE(fname,1) 'restart_',MYRANK, '.in'
      OPEN(65,file=fname,status='old',action='read')
      READ(65,*) t,u,p,rhs1h
      CLOSE(65)
    ELSEIF(irestart==2) then
!      ! irestart = 2 :: new run with initialization using previous data, interpolation needed. Do not read in rhs1h
!
!      ! read nxrst, nyrst, nzrst
!      WRITE(fname,1) 'proc_',MYRANK,'.dat'
!      OPEN(11,FILE=fname,STATUS='old',ACTION='read')
!      READ(11,2)
!      READ(11,2) nxrst, nyrst, nzrst
!      READ(11,2) sxrst, syrst, szrst
!      READ(11,2) exrst, eyrst, ezrst
!      CLOSE(11)
!      allocate(xrst(nxrst),yrst(nyrst),zrst(nzrst),urst(nxrst,nyrst,nzrst,nvar),prst(nxrst,nyrst,nzrst))
!      allocate(xgrst(nxgrst),ygrst(nygrst),zgrst(nzgrst))
!
!      ! read xrst, yrst, zrst, xgrst, ygrst, zgrst
!      WRITE(fname,1) 'grid.dat'
!      OPEN(11,FILE=fname,STATUS='old',ACTION='read')
!      ! read xrst, xgrst
!      read(11,*)
!      do i = 1, nxgrst
!        read(11,*) xgrst(i)
!      enddo
!      !xrst(:) = xgrst(sxrst:exrst)
!      read(11,*)
!      ! read yrst, ygrst
!      read(11,*)
!      do j = 1, nygrst
!        read(11,*) ygrst(j)
!      enddo
!      !yrst(:) = ygrst(syrst:eyrst)
!      read(11,*)
!      ! read zrst, zgrst
!      read(11,*)
!      do k = 1, nzgrst
!        read(11,*) zgrst(k)
!      enddo
!      read(11,*)
!      !zrst(:) = zgrst(szrst:ezrst)
!      CLOSE(11)
!
!      ! read u and p
!      WRITE(fname,1) 'restart_',MYRANK, '.in'
!      OPEN(65,file=fname,status='old',action='read')
!      READ(65,*) t,urst,prst
!      CLOSE(65)
!
!      ! reset t to 0 because this is a new simulation
!      t = 0.0d0
!
!      ! only z decomposition for now
!      ! determine extents of zrst needed for interpolation
!      ! lower end
!      kst = -1
!      do k = 1, nzgrst-1
!        if(zgrst(k)<z(sz) .and. zgrst(k+1)>z(sz)) then
!          kst = k
!          exit
!        endif
!      enddo
!      if(kst==-1) then
!        write(*,*) 'z(sz) is not bracketed by zgrst. Check details'
!        write(*,*) z(sz), zgrst(1), zgrst(nzgrst) 
!      endif
!      ! lower end
!      kend = -1
!      do k = 1, nzgrst-1
!        if(zgrst(k)<z(ez) .and. zgrst(k+1)>z(ez)) then
!          kend = k
!          exit
!        endif
!      enddo
!      if(kend==-1) then
!        write(*,*) 'z(ez) is not bracketed by zgrst. Check details'
!        write(*,*) z(ez), zgrst(1), zgrst(nzgrst) 
!      endif
!
!      if(kst<szrst) then
!        ! communication needed
!      endif
!      if(kend>ezrst) then
!        ! communication needed
!      endif
!
!      do k = sz, ez
!        do krst = kst, kend
!          if(zgrst(krst)<z(k) .and. zgrst(krst+1)>z(k)) exit
!        enddo
!        ! interpolate
!        if(krst<szrst) then
!        elseif(krst>ezrst) then
!        else
!          xl = zgrst(krst); xr = zgrst(krst+1); xm = z(k); fac = (xm-xl)/(xr-xl)
!          u(:,:,k,:) = urst(:,:,krst,:) + fac*(u(:,:,krst+1,:) - u(:,:,krst,:))
!        endif
!      enddo
!      kst = szrst
!      if(zgrst(szrst)>z(sz)) then
!        
!      endif
!
!      deallocate(xrst,yrst,zrst,urst,prst)
    ELSE
      ! irestart = 0 :: start new run; you should not be here
      write(*,*) 'Wrong irestart, ', irestart
      stop
    ENDIF

    !t = 0.0d0
    !do k = sz, ez
    !  u(:,:,k,4) = 0.5d0*(1.0d0+z(k))
    !enddo
    !u(:,:,:,4) = u(:,:,:,4)-0.5d0

    ! set pressure gradient
    pgr = 9.812d0/Re
    ALLOCATE(tmp(sx:ex,sy:ey,sz:ez))
    tmp(:,:,:) = itfac*dt*pgr
    call p3dfft_ftran_r2c(tmp(:,:,:), dt_pgradh(:,:,:),fstr)
    DEALLOCATE(tmp)

    ! set onetr for Corriolis force
    ALLOCATE(tmp(sx:ex,sy:ey,sz:ez))
    tmp(:,:,:) = itfac*1.0d0
    call p3dfft_ftran_r2c(tmp(:,:,:), onetr(:,:,:),fstr)
    DEALLOCATE(tmp)

    call bc

    ! set bounds for rho
    Rmax = 0.5d0
    Rmin = -0.5d0

    do var = 1, nvar
      call p3dfft_ftran_r2c(u(:,:,:,var),uh(:,:,:,var),fstr)
    end do
    uh = itfac*uh
    ph = 0.0d0

    ! szf = 1; ezf = nzg
    uprofh(:,:,2) = uh(:,:,ezf,1)
    vprofh(:,:,2) = uh(:,:,ezf,2)

    tprofh(:,:,1) = uh(:,:,szf,4)
    tprofh(:,:,2) = uh(:,:,ezf,4)

    call vinterp
    call vinterph
    call bch
    call dealiasuh
    call decode

    call divergence
    IF(MYRANK==outid) then
      write(*,35) 'Init: ',t, maxdivu, maxu1, maxu2, maxcfl
    ENDIF

1 FORMAT(a,i4.4,a)
35 format(a,7(1x,e12.5),1x,i5,1x,i9)

  END SUBROUTINE restartinit

  SUBROUTINE decode

    IMPLICIT NONE
    INTEGER :: var

    ! backward transform velocities and density
    do var = 1, nvar
      tmph(:,:,:) = uh(:,:,:,var)
      call p3dfft_btran_c2r(tmph(:,:,:),u(:,:,:,var),bstr)
    enddo

    ! backward transform pressure
    tmph(:,:,:) = ph(:,:,:)
    call p3dfft_btran_c2r(tmph(:,:,:),p(:,:,:),bstr)

  END SUBROUTINE decode

  SUBROUTINE exchangez(iwpov)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iwpov ! 1 :: (overwrite wp), 0 :: (no overwrite)
    INTEGER :: var, sendcnt, recvcnt, mpierr
    REAL(DP) :: xl, xr, xm, fac

    ! send zm1
    do var = 1, nvar
      sendbuf(:,:,var) = u(:,:,ez,var)
    enddo
    sendbuf(:,:,nvar+1) = wp(:,:,ez)
    sendcnt = nx*ny*(nvar+1); recvcnt = sendcnt
    call MPI_SENDRECV(sendbuf,sendcnt,MPIRP,MPIDEST,0,recvbuf,recvcnt,MPIRP,MPISRC,0,MPICOMM2D,mpistat,mpierr)
    ubufzm1(:,:,:) = recvbuf(:,:,:)

    ! send zp1
    sendbuf(:,:,1:nvar) = u(:,:,sz,1:nvar)
    sendbuf(:,:,nvar+1) = wp(:,:,sz)
    sendcnt = nx*ny*(nvar+1); recvcnt = sendcnt
    call MPI_SENDRECV(sendbuf,sendcnt,MPIRP,MPISRC,0,recvbuf,recvcnt,MPIRP,MPIDEST,0,MPICOMM2D,mpistat,mpierr)
    ubufzp1(:,:,1:nvar+1) = recvbuf(:,:,1:nvar+1)
    IF(iwpov==1) then ! overwrite wp
      IF(istretch==0) then
        ubufzp1(:,:,nvar+1) = 0.5d0*(u(:,:,ez,3) + ubufzp1(:,:,3))
      ELSE
        IF(ez/=nzg) then
          xl = zz(ez); xr = zz(ez+1); xm = z(ez+1); fac = (xm-xl)/(xr-xl)
          ubufzp1(:,:,nvar+1) = u(:,:,ez,3) + fac*(ubufzp1(:,:,3) - u(:,:,ez,3))
        ENDIF
      ENDIF
    ENDIF

    ! boundary conditions
    IF(sz==1) then
      ubufzm1(:,:,1:2) = -u(:,:,2,1:2)  ! Not staggered
      ubufzm1(:,:,3) = -u(:,:,1,3)      ! Staggered
      ubufzm1(:,:,4) = 2.0d0*u(:,:,1,4)-u(:,:,2,4)      ! Not staggered
      ubufzm1(:,:,nvar+1) = 2.0d0*wp(:,:,1)-wp(:,:,2)   ! Not staggered
    ENDIF
    IF(ez==nzg) then
      ubufzp1(:,:,1:2) = u(:,:,nzg-1,1:2)  ! Not staggered
      ubufzp1(:,:,3) = -u(:,:,nzg-1,3)      ! Staggered
      ubufzp1(:,:,4) = u(:,:,nzg-1,4)      ! Not staggered
      ubufzp1(:,:,nvar+1) = 2.0d0*wp(:,:,nzg)-wp(:,:,nzg-1)   ! Not staggered
    ENDIF

  END SUBROUTINE exchangez

  SUBROUTINE vinterph

    IMPLICIT NONE
    INTEGER :: k
    REAL(DP) :: xl, xr, xm, fac

    wph = 0.0d0
    ! interior points
    if(istretch==0)then
      do k = szf+1, ezf-1   ! can also be 2, nzg-1
        wph(:,:,k) = 0.5d0*(uh(:,:,k,3) + uh(:,:,k-1,3))
      end do
    else
      do k = szf+1, ezf-1   ! can also be 2, nzg-1
        xl=zz(k-1); xr=zz(k); xm=z(k); fac = (xm-xl)/(xr-xl)
        wph(:,:,k)=uh(:,:,k-1,3)+fac*(uh(:,:,k,3) - uh(:,:,k-1,3))
      end do
    endif

    ! boundary conditions - can also be 1, nzg
    wph(:,:,szf) = 0.0d0
    wph(:,:,ezf) = 0.0d0

  END SUBROUTINE vinterph

  SUBROUTINE vinterp

    IMPLICIT NONE
    INTEGER :: k
    REAL(DP) :: xl, xr, xm, fac

    wp = 0.0d0
    ! interior points
    if(istretch==0)then
      do k = sz+1, ez-1
        wp(:,:,k) = 0.5d0*(u(:,:,k,3) + u(:,:,k-1,3))
      end do
    else
      do k = sz+1, ez-1
        xl=zz(k-1); xr=zz(k); xm=z(k); fac = (xm-xl)/(xr-xl)
        wp(:,:,k)=u(:,:,k-1,3)+fac*(u(:,:,k,3) - u(:,:,k-1,3))
      end do
    endif

    if(sz.ne.1) then
      if(istretch==0) then
        wp(:,:,sz) = 0.5d0*(u(:,:,sz,3) + ubufzm1(:,:,3))
      else
        xl=zz(sz-1); xr=zz(sz); xm=z(sz); fac = (xm-xl)/(xr-xl)
        wp(:,:,sz)=ubufzm1(:,:,3)+fac*(u(:,:,sz,3) - ubufzm1(:,:,3))
      endif
    endif

    if(ez.ne.nzg) then
      if(istretch==0) then
        wp(:,:,ez) = 0.5d0*(u(:,:,ez,3) + u(:,:,ez-1,3))
      else
        xl=zz(ez-1); xr=zz(ez); xm=z(ez); fac = (xm-xl)/(xr-xl)
        wp(:,:,ez)=u(:,:,ez-1,3)+fac*(u(:,:,ez,3) - u(:,:,ez-1,3))
      endif
    endif

    call bcp

  END SUBROUTINE vinterp

  SUBROUTINE xyderivs

    IMPLICIT NONE
    INTEGER :: var, i, j

    ! x derivatives
    do var = 1, nvar
     do i = sxf, exf
       tmph(i,:,:) = iu*k1(i)*uh(i,:,:,var)
     enddo
     call p3dfft_btran_c2r(tmph,dudx(:,:,:,var),bstr)
    enddo

    ! dwpdx
    do i = sxf, exf
      tmph(i,:,:) = iu*k1(i)*wph(i,:,:)
    enddo
    call p3dfft_btran_c2r(tmph,dwpdx(:,:,:),bstr)

    ! y derivatives
    do var = 1, nvar
     do j = syf, eyf
       tmph(:,j,:) = iu*k2(j)*uh(:,j,:,var)
     enddo
     call p3dfft_btran_c2r(tmph,dudy(:,:,:,var),bstr)
    enddo

    ! dwpdy
    do j = syf, eyf
      tmph(:,j,:) = iu*k2(j)*wph(:,j,:)
    enddo
    call p3dfft_btran_c2r(tmph,dwpdy(:,:,:),bstr)

  END SUBROUTINE xyderivs

  SUBROUTINE zderivs

    IMPLICIT NONE
    INTEGER :: i,j

    do j = sy, ey
     do i = sx, ex
       call deriv1cen(dudz(i,j,:,1),u(i,j,:,1),ubufzm1(i,j,1),ubufzp1(i,j,1),0)
       call deriv1cen(dudz(i,j,:,2),u(i,j,:,2),ubufzm1(i,j,2),ubufzp1(i,j,2),0)
       call deriv1cen(dudz(i,j,:,3),u(i,j,:,3),ubufzm1(i,j,3),ubufzp1(i,j,3),1) ! staggered: BC handled via ubuf
       call deriv1cen(dudz(i,j,:,4),u(i,j,:,4),ubufzm1(i,j,4),ubufzp1(i,j,4),0)
       call deriv1cen(dwpdz(i,j,:), wp(i,j,:), ubufzm1(i,j,5),ubufzp1(i,j,5),0)
     enddo
    enddo

  END SUBROUTINE zderivs

  SUBROUTINE subzderivs

    IMPLICIT NONE

    dudx(:,:,:,3) = dwpdx
    dudy(:,:,:,3) = dwpdy
    dudz(:,:,:,3) = dwpdz

  END SUBROUTINE subzderivs


  ! 1st derivative centered differences complex
  SUBROUTINE deriv1cenh(phider,phi,ibclow,ibcupp)

    IMPLICIT NONE
    COMPLEX(DP), DIMENSION(szf:ezf), INTENT(OUT) :: phider
    COMPLEX(DP), DIMENSION(szf:ezf), INTENT(IN) :: phi
    INTEGER, INTENT(IN) :: ibclow, ibcupp
    INTEGER :: k
    REAL(DP) :: fac1, fac2

    if(ibclow==1) then
      ! explicitly set derivative to zero
      k = szf
      phider(k) = zero
    else
      ! one-sided derivatives at szf and ezf (or 1, nzg)
      k = szf
      phider(k) = (phi(k+1)-phi(k))/(z(k+1)-z(k))
    endif

    if(ibcupp==1) then
      k = ezf
      phider(k) = zero
    else
      k = ezf
      phider(k) = (phi(k)-phi(k-1))/(z(k)-z(k-1))
    endif

    ! interior
    do k = szf+1, ezf-1
      fac1 = (z(k)-z(k-1))/(z(k+1)-z(k)); fac2 = z(k+1)-z(k-1)
      phider(k) = (fac1*(phi(k+1)-phi(k)) + (phi(k)-phi(k-1))/fac1) / fac2
    enddo

  END SUBROUTINE deriv1cenh

  ! 1st derivative centered differences
  SUBROUTINE deriv1cen(phider,phi,philow,phiupp,istag)

    IMPLICIT NONE
    REAL(DP), DIMENSION(sz:ez), INTENT(OUT) :: phider
    REAL(DP), DIMENSION(sz:ez), INTENT(IN) :: phi
    REAL(DP), INTENT(IN) :: philow, phiupp
    INTEGER, INTENT(IN) :: istag
    INTEGER :: k, klow, kupp
    REAL(DP) :: fac1, fac2

    IF(istag==0) then ! non-staggered

      IF(sz==1) then
        k = sz
        phider(k) = 0.5d0*(phi(k+1)-philow)/(z(k+1)-z(k))
      ELSE
        k = sz; fac1 = (z(k)-z(k-1))/(z(k+1)-z(k)); fac2 = z(k+1)-z(k-1)
        phider(k) = (fac1*(phi(k+1)-phi(k)) + (phi(k)-philow)/fac1) / fac2
      ENDIF
      IF(ez==nzg) then
        k = ez
        phider(k) = 0.5d0*(phiupp-phi(k-1))/(z(k)-z(k-1))
      ELSE
        k = ez; fac1 = (z(k)-z(k-1))/(z(k+1)-z(k)); fac2 = z(k+1)-z(k-1)
        phider(k) = (fac1*(phiupp-phi(k)) + (phi(k)-phi(k-1))/fac1) / fac2
      ENDIF
      do k = sz+1, ez-1
        fac1 = (z(k)-z(k-1))/(z(k+1)-z(k)); fac2 = z(k+1)-z(k-1)
        phider(k) = (fac1*(phi(k+1)-phi(k)) + (phi(k)-phi(k-1))/fac1) / fac2
      enddo

    ELSE ! staggered

      klow = sz+1; kupp = ez-1

      IF(sz==1) then
        k = sz
        phider(k) = 0.5d0*(phi(k+1)-philow)/(zz(k+1)-zz(k))
      ELSE
        k = sz; fac1 = (zz(k)-zz(k-1))/(zz(k+1)-zz(k)); fac2 = zz(k+1)-zz(k-1)
        phider(k) = (fac1*(phi(k+1)-phi(k)) + (phi(k)-philow)/fac1) / fac2
      ENDIF
      IF(ez==nzg) then
        kupp = ez-2
        phider(nzg-1) = 0.5d0*(phiupp-phi(nzg-2))/(zz(nzg-1)-zz(nzg-2))
      ELSE
        k = ez; fac1 = (zz(k)-zz(k-1))/(zz(k+1)-zz(k)); fac2 = zz(k+1)-zz(k-1)
        phider(k) = (fac1*(phiupp-phi(k)) + (phi(k)-phi(k-1))/fac1) / fac2
      ENDIF
      do k = klow, kupp
        fac1 = (zz(k)-zz(k-1))/(zz(k+1)-zz(k)); fac2 = zz(k+1)-zz(k-1)
        phider(k) = (fac1*(phi(k+1)-phi(k)) + (phi(k)-phi(k-1))/fac1) / fac2
      enddo

    ENDIF

  END SUBROUTINE deriv1cen

  ! 2nd derivative centered differences
  SUBROUTINE deriv2cen(phider,phi,ibclow,ibcupp,istag)

    IMPLICIT NONE
    COMPLEX(DP), DIMENSION(szf:ezf), INTENT(OUT) :: phider
    COMPLEX(DP), DIMENSION(szf:ezf), INTENT(IN) :: phi
    INTEGER, INTENT(IN) :: ibclow, ibcupp, istag
    INTEGER :: k

    IF(istag==0) then ! non-staggered

      IF(ibclow==0) then ! Dirichlet
        phider(szf) = 0.0d0
      ELSEIF(ibclow==1) then ! Neumann
        phider(szf) = 2.0d0*(phi(2)-phi(1))/(z(2)-z(1))**2
      ENDIF

      IF(ibcupp==0) then ! Dirichlet
        phider(ezf) = 0.0d0 
      ELSEIF(ibcupp==1) then ! Neumann
        phider(ezf) = 2.0d0*(phi(nzg-1)-phi(nzg))/(z(nzg)-z(nzg-1))**2
      ENDIF

      do k = szf+1, ezf-1 ! can also be 2, nzg-1
        phider(k) = 2.0d0*phi(k+1)/((z(k+1)-z(k-1))*(z(k+1)-z(k))) + &
                    2.0d0*phi(k-1)/((z(k+1)-z(k-1))*(z(k)-z(k-1))) - &
                    2.0d0*phi(k)/((z(k+1)-z(k))*(z(k)-z(k-1)))
      enddo

    ELSE ! staggered

      IF(ibclow==0) then ! Dirichlet
        phider(szf) = 2.0d0/((zz(2)+zz(1)-2.0d0*z(1))*(zz(2)-zz(1)))*(phi(2)-phi(1)*(zz(2)-z(1))/(zz(1)-z(1)))
      ELSEIF(ibclow==1) then ! Neumann
        phider(szf) = 2.0d0/((zz(2)+zz(1)-2.0d0*z(1))*(zz(2)-zz(1)))*(phi(2)-phi(1))
      ENDIF

      IF(ibcupp==0) then ! Dirichlet
        phider(ezf-1) = 2.0d0/((2.0d0*z(nzg)-zz(nzg-1)-zz(nzg-2))*(zz(nzg-1)-zz(nzg-2)))*(phi(nzg-2)-phi(nzg-1)*(z(nzg)-zz(nzg-2))/(z(nzg)-zz(nzg-1)))
      ELSEIF(ibcupp==1) then ! Neumann
        phider(ezf-1) = 2.0d0/((2.0d0*z(nzg)-zz(nzg-1)-zz(nzg-2))*(zz(nzg-1)-zz(nzg-2)))*(phi(nzg-2)-phi(nzg-1))
      ENDIF

      do k = szf+1, ezf-2 ! can also be 2, nzg-2
        phider(k) = 2.0d0*phi(k+1)/((zz(k+1)-zz(k-1))*(zz(k+1)-zz(k))) + &
                    2.0d0*phi(k-1)/((zz(k+1)-zz(k-1))*(zz(k)-zz(k-1))) - &
                    2.0d0*phi(k)/((zz(k+1)-zz(k))*(zz(k)-zz(k-1)))
      enddo

    ENDIF

  END SUBROUTINE deriv2cen

  SUBROUTINE process

    IMPLICIT NONE
    integer::i,j,k
    INTEGER :: iwork(5*3),ifail(3),eigfound,info
    REAL(DP)::eigv_tol
    REAL(DP) :: eigvals(3),eigvects(3,3),work(105),Qij(3,3)


    !call subzderivs

    eigv_tol=1.0d-30;

    do k = sz, ez
       do j = sy, ey
          do i = sx, ex


             ! not calculating inv2gradvel
!!$             inv2gradvel(i,j,k)=(dudx(i,j,k)*dvdy(i,j,k) + &
!!$                  dudx(i,j,k)*dwdz(i,j,k) + &
!!$                  dvdy(i,j,k)*dwdz(i,j,k)) - &
!!$                  (dudy(i,j,k)*dvdx(i,j,k) + &
!!$                  dudz(i,j,k)*dwdx(i,j,k) + &
!!$                  dvdz(i,j,k)*dwdy(i,j,k))

             Qij(1,1)=dudx(i,j,k,1)**2+ &
                  dudy(i,j,k,1)*dudx(i,j,k,2)+dudz(i,j,k,1)*dwpdx(i,j,k)
             Qij(2,2)=dudy(i,j,k,2)**2+ &
                  dudx(i,j,k,2)*dudy(i,j,k,1)+dudz(i,j,k,2)*dwpdy(i,j,k)
             Qij(3,3)=dwpdz(i,j,k)**2+ &
                  dwpdx(i,j,k)*dudz(i,j,k,1)+dwpdy(i,j,k)*dudz(i,j,k,2)
             Qij(2,1)=0.50*((dudx(i,j,k,1)+dudy(i,j,k,2))* &
                  (dudx(i,j,k,2)+dudy(i,j,k,1)) &
                  + dwpdx(i,j,k)*dudz(i,j,k,2) &
                  + dwpdy(i,j,k)*dudz(i,j,k,1))
             Qij(3,1)=0.50*((dudx(i,j,k,1)+dwpdz(i,j,k))* &
                  (dwpdx(i,j,k)+dudz(i,j,k,1)) &
                  + dudx(i,j,k,2)*dwpdy(i,j,k) &
                  + dudz(i,j,k,2)*dudy(i,j,k,1))
             Qij(3,2)=0.50*((dudy(i,j,k,2)+dwpdz(i,j,k))* &
                  (dwpdy(i,j,k)+dudz(i,j,k,2)) &
                  + dudy(i,j,k,1)*dwpdx(i,j,k) &
                  + dudz(i,j,k,1)*dudx(i,j,k,2))

             ! Now that we've filled the Q matrix for this point,
             ! calculate the eigenvalues.  Select the second eigenvalue
             ! and place it in the result matrix for later use.

             call DSYEVX('N',    &     ! JOBZ: Eigenvalues only
                  'I',      &    ! RANGE: Specify eigenvalues by index
                  'L',      &    ! UPLO: Lower triangular matrix supplied
                  3,        &     ! N: Order of matrix
                  Qij,      &     ! A: Data matrix (destroyed on return)
                  3,        &     ! LDA: Leading dimension of Qij
                  0.d0,     &    ! VL: Ignored since RANGE='I'
                  0.d0,     &    ! VU: Ignored since RANGE='I'
                  2,        &    ! IL: Just calculate second eigenvalue
                  2,        &     ! IU: Just calculate second eigenvalue
                  eigv_tol, &   ! ABSTOL: Convergence criterion
                  eigfound, &   ! M: number of eigenvalues found
                  eigvals,  &   ! W: Eigenvalue vector
                  eigvects, &   ! Z: Eigenvector array (not referenced)
                  3,        &    ! LDZ: Leading dimension of Z
                  work,     &     ! WORK: Work array (work(1)=opt. LWORK)
                  105,      &     ! LWORK: Length of LWORK
                  iwork,    &     ! IWORK: Integer workspace
                  ifail,    &     ! IFAIL: not referenced (JOBZ='N')
                  info)        ! INFO: =0: success, <0: -info arg has
             !       illegal val, >0: info eigvects
             !       failed to converge (indices in
             !       ifail)
             secondEigv(i,j,k)=eigvals(1) ! since only calculating one, stored

             ! in first index of eigvals()
!!$ vortmag(i,j,k)=sqrt(vort(i,j,k,1)**2+vort(i,j,k,2)**2+vort(i,j,k,3)**2)

          enddo
       enddo
    enddo

  END SUBROUTINE process

  SUBROUTINE writegrid

    IMPLICIT NONE
    INTEGER :: i, j, k

    IF(MYRANK==0) then
      OPEN(10,FILE='grid.dat',STATUS='replace')
      WRITE(10,*) 'x'
      DO i = 1, nxg
        WRITE(10,*) x(i)
      ENDDO
      WRITE(10,*) 'Done x'
      WRITE(10,*) 'y'
      DO j = 1, nyg
        WRITE(10,*) y(j)
      ENDDO
      WRITE(10,*) 'Done y'
      WRITE(10,*) 'z'
      DO k = 1, nzg
        WRITE(10,*) z(k)
      ENDDO
      WRITE(10,*) 'Done z'
      CLOSE(10)
    ENDIF

  END SUBROUTINE writegrid

  SUBROUTINE writedata

    IMPLICIT NONE
    INTEGER :: i, j, k
    CHARACTER(LEN=17) :: fname

    itsnap = itsnap + 1
    WRITE(fname,1) 'sol_',MYRANK,'_',itsnap,'.dat'
    OPEN(10,FILE=fname,STATUS='replace')
    WRITE(10,2) sx, ex, sy, ey, sz, ez, t
    DO k = sz, ez
     DO j = sy, ey
      DO i = sx, ex
        WRITE(10,3) u(i,j,k,1), u(i,j,k,2), wp(i,j,k), u(i,j,k,4), p(i,j,k), secondEigv(i,j,k)
      ENDDO
     ENDDO
    ENDDO
    CLOSE(10)

1 FORMAT(a,i4.4,a,i4.4,a)
2 FORMAT(6(i4.4,1x),e19.12)
3 FORMAT(6(e19.12,1x))

  END SUBROUTINE writedata

  SUBROUTINE outputIC

    IMPLICIT NONE

    call decode
    iwpov = 1; call exchangez(iwpov)
    call vinterp
    call vinterph
    !call bc
    call xyderivs
    call zderivs
    call process
    call writegrid
    itsnap = -1
    call writedata

  END SUBROUTINE outputIC

  SUBROUTINE init_tridiag

    IMPLICIT NONE
    INTEGER :: k

    amg = 0.0d0; bmg = 0.0d0; cmg = 0.0d0
    DO k = szf, ezf  ! can also be 1, nzg
      IF(k==szf) then
        amg(k) = -2.0d0/(z(2)-z(1))**2
        cmg(k) = 2.0d0/(z(2)-z(1))**2
      ELSEIF(k==ezf) then
        amg(k) = -2.0d0/(z(nzg)-z(nzg-1))**2
        bmg(k) = 2.0d0/(z(nzg)-z(nzg-1))**2
      ELSE
        amg(k) = -2.0d0/((z(k+1)-z(k))*(z(k)-z(k-1)))
        bmg(k) = 2.0d0/((z(k+1)-z(k-1))*(z(k)-z(k-1)))
        cmg(k) = 2.0d0/((z(k+1)-z(k-1))*(z(k+1)-z(k)))
      ENDIF
    ENDDO

  END SUBROUTINE init_tridiag

  SUBROUTINE init_tridiag_vel

    IMPLICIT NONE
    INTEGER :: k

    ! tangential: for u1, u2, with fixed values (lower) and zero gradient
    ! (upper) BC inbuilt
    amt = 0.0d0; bmt = 0.0d0; cmt = 0.0d0
    DO k = szf, ezf  ! can also be 1, nzg
      IF(k==szf) then
        amt(k) = 1.0d0
      ELSEIF(k==ezf) then
        amt(k) = 1.0d0+dt/(2.0d0*Re)*2.0d0/((z(k)-z(k-1))**2.0d0)
        cmt(k) = -dt/(2.0d0*Re)*2.0d0/((z(k)-z(k-1))**2.0d0)
      ELSE
        amt(k) = 1.0d0+dt/(2.0d0*Re)*2.0d0/((z(k+1)-z(k))*(z(k)-z(k-1)))
        bmt(k) = -dt/(2.0d0*Re)*2.0d0/((z(k+1)-z(k-1))*(z(k)-z(k-1)))
        cmt(k) = -dt/(2.0d0*Re)*2.0d0/((z(k+1)-z(k-1))*(z(k+1)-z(k)))
      ENDIF
    ENDDO

    ! normal: for u3, with zero gradient values BC inbuilt
    amn = 0.0d0; bmn = 0.0d0; cmn = 0.0d0
    DO k = szf, ezf  ! can also be 1, nzg
      IF(k==szf) then
        amn(k) = 1.0d0+dt/(2.0d0*Re)*2.0d0/(zz(2)-zz(1))**2
        cmn(k) = -dt/(2.0d0*Re)*2.0d0/(zz(2)-zz(1))**2
      ELSEIF(k==ezf-1) then
        amn(k) = 1.0d0+dt/(2.0d0*Re)*2.0d0/(zz(nzg-1)-zz(nzg-2))**2
        bmn(k) = -dt/(2.0d0*Re)*2.0d0/(zz(nzg-1)-zz(nzg-2))**2
      ELSEIF(k==ezf) then
      ELSE
        amn(k) = 1.0d0+dt/(2.0d0*Re)*2.0d0/((zz(k+1)-zz(k))*(zz(k)-zz(k-1)))
        bmn(k) = -dt/(2.0d0*Re)*2.0d0/((zz(k+1)-zz(k-1))*(zz(k)-zz(k-1)))
        cmn(k) = -dt/(2.0d0*Re)*2.0d0/((zz(k+1)-zz(k-1))*(zz(k+1)-zz(k)))
      ENDIF
    ENDDO

    ! tangential: for rho, with fixed values BC inbuilt
    amr = 0.0d0; bmr = 0.0d0; cmr = 0.0d0
    DO k = szf, ezf  ! can also be 1, nzg
      IF(k==szf) then
        amr(k) = 1.0d0
      ELSEIF(k==ezf) then
        amr(k) = 1.0d0+dt/(2.0d0*Re*Pr)*2.0d0/((z(nzg)-z(nzg-1))**two)
        cmr(k) = -dt/(2.0d0*Re*Pr)*2.0d0/((z(nzg)-z(nzg-1))**two)
      ELSE
        amr(k) = 1.0d0+dt/(2.0d0*Re*Pr)*2.0d0/((z(k+1)-z(k))*(z(k)-z(k-1)))
        bmr(k) = -dt/(2.0d0*Re*Pr)*2.0d0/((z(k+1)-z(k-1))*(z(k)-z(k-1)))
        cmr(k) = -dt/(2.0d0*Re*Pr)*2.0d0/((z(k+1)-z(k-1))*(z(k+1)-z(k)))
      ENDIF
    ENDDO

  END SUBROUTINE init_tridiag_vel

  SUBROUTINE convect

    IMPLICIT NONE
    INTEGER :: i, j, k, var
    REAL(DP) :: xl, xr, xm, fac

    ! non staggered
    conv(:,:,:,1) = u(:,:,:,1)*dudx(:,:,:,1) + u(:,:,:,2)*dudy(:,:,:,1) + wp(:,:,:)*dudz(:,:,:,1)
    conv(:,:,:,2) = u(:,:,:,1)*dudx(:,:,:,2) + u(:,:,:,2)*dudy(:,:,:,2) + wp(:,:,:)*dudz(:,:,:,2)
    conv(:,:,:,4) = u(:,:,:,1)*dudx(:,:,:,4) + u(:,:,:,2)*dudy(:,:,:,4) + wp(:,:,:)*dudz(:,:,:,4)

    ! staggered
    do k = sz, ez-1
     xl = z(k); xr = z(k+1); xm = zz(k); fac = (xm-xl)/(xr-xl)
     ubar(:,:) = u(:,:,k,1) + fac*(u(:,:,k+1,1)-u(:,:,k,1))
     vbar(:,:) = u(:,:,k,2) + fac*(u(:,:,k+1,2)-u(:,:,k,2))
     conv(:,:,k,3) = ubar(:,:)*dudx(:,:,k,3) + vbar(:,:)*dudy(:,:,k,3) + u(:,:,k,3)*dudz(:,:,k,3)
    enddo
    k = ez
    if(ez==nzg) then
      conv(:,:,k,3) = 0.0d0
    else
      xl = z(k); xr = z(k+1); xm = zz(k); fac = (xm-xl)/(xr-xl)
      ubar(:,:) = u(:,:,k,1) + fac*(ubufzp1(:,:,1)-u(:,:,k,1))
      vbar(:,:) = u(:,:,k,2) + fac*(ubufzp1(:,:,2)-u(:,:,k,2))
      conv(:,:,k,3) = ubar(:,:)*dudx(:,:,k,3) + vbar(:,:)*dudy(:,:,k,3) + u(:,:,k,3)*dudz(:,:,k,3)
    endif

    do var = 1, nvar
      call p3dfft_ftran_r2c(conv(:,:,:,var),convh(:,:,:,var),fstr)
    end do
    convh = itfac*convh

  END SUBROUTINE convect

  SUBROUTINE viscous

    IMPLICIT NONE
    INTEGER :: i, j, var

    do j = syf, eyf
     do i = sxf, exf
       call deriv2cen(d2udzh(i,j,:,1),uh(i,j,:,1),0,1,0)
       call deriv2cen(d2udzh(i,j,:,2),uh(i,j,:,2),0,1,0)
       call deriv2cen(d2udzh(i,j,:,3),uh(i,j,:,3),0,0,1)
       call deriv2cen(d2udzh(i,j,:,4),uh(i,j,:,4),0,1,0)
     enddo
    enddo

    do j = syf, eyf
     do i = sxf, exf
      do var = 1, 3
        viscth(i,j,:,var) = -(k1(i)**2+k2(j)**2)*uh(i,j,:,var)/Re
        viscnh(i,j,:,var) = 0.5d0*d2udzh(i,j,:,var)/Re
      enddo
      ! rho
      var = 4
      viscth(i,j,:,var) = -(k1(i)**2+k2(j)**2)*uh(i,j,:,var)/(Re*Pr)
      viscnh(i,j,:,var) = 0.5d0*d2udzh(i,j,:,var)/(Re*Pr)
     enddo
    enddo

  END SUBROUTINE viscous

  SUBROUTINE Corriolis

    IMPLICIT NONE
    REAL(DP) :: fac

    !fac = 2.0d0/Re
    !fac = 1.0d0/9.9d0
    fac = 1.0d0/10.0d0

    ! Corriolis forcing for Ekman layer
    fCorriolh(:,:,:,1) = fac*uh(:,:,:,2)                    ! x direction
    fCorriolh(:,:,:,2) = fac*(onetr(:,:,:) - uh(:,:,:,1))   ! y direction

  END SUBROUTINE Corriolis


  SUBROUTINE buoyancy

    IMPLICIT NONE
    INTEGER :: k, mpierr
    REAL(DP) :: xl, xr, xm, fac

    if(istretch==0)then
      do k = szf, ezf-1     ! 1, nzg-1
        rhowh(:,:,k) = 0.5d0*(uh(:,:,k,4) + uh(:,:,k+1,4))
      end do
    else
      do k = szf, ezf-1    ! 1, nzg-1
        xl=z(k); xr=z(k+1); xm=zz(k); fac = (xm-xl)/(xr-xl)
        rhowh(:,:,k)=uh(:,:,k,4)+fac*(uh(:,:,k+1,4) - uh(:,:,k,4))
      end do
    endif

    rhomhl = (0.0d0,0.0d0); rhomh = (0.0d0,0.0d0)
    do k = szf, ezf
      rhomhl(k) = sum(rhowh(:,:,k))/dble(nxg*nyg)
    enddo
    call MPI_ALLREDUCE(rhomhl,rhomh,nzg,MPICP,MPI_SUM,MPICOMM2D,mpierr)
    do k = szf, ezf
      buoyh(:,:,k,3) = -Ri*(rhowh(:,:,k)-rhomh(k))
    enddo

  END SUBROUTINE buoyancy

      SUBROUTINE gradmodel

      IMPLICIT NONE
      INTEGER :: i,j,k, ind,jnd
      REAL(DP) :: gbar(3,4),gbarden,gbartden,pamodelv,pamodelt,pamodel,pamodel2,pgrmod

      ievm = 0
      call magstrainrate

      ! store w derivatives
      uxtemp(:,:,:,1) = dudx(:,:,:,3)
      uxtemp(:,:,:,2) = dudy(:,:,:,3)
      uxtemp(:,:,:,3) = dudz(:,:,:,3)
      ! copy wp derivatives
      dudx(:,:,:,3) = dwpdx(:,:,:)
      dudy(:,:,:,3) = dwpdy(:,:,:)
      dudz(:,:,:,3) = dwpdz(:,:,:)

      do k = sz, ez
       do j = sy, ey
        do i = sx, ex
          ! Porte-Agel model
          do ind=1,3
           do jnd=1,4
            gbar(ind,jnd) = delx2*dudx(i,j,k,ind)*dudx(i,j,k,jnd) + &
                            dely2*dudy(i,j,k,ind)*dudy(i,j,k,jnd) + &
                            delz2(k)*dudz(i,j,k,ind)*dudz(i,j,k,jnd)
           enddo
          enddo
          gbarden = gbar(1,1)+gbar(2,2)+gbar(3,3)
          gbartden = dsqrt(sum(gbar(:,4)**2))

          pamodelv = -(gbar(1,1)*dudx(i,j,k,1) + gbar(2,2)*dudy(i,j,k,2) + &
                  gbar(3,3)*dudz(i,j,k,3) +                                &
            0.5d0*(gbar(1,2)+gbar(2,1))*(dudx(i,j,k,2)+dudy(i,j,k,1)) +    &
            0.5d0*(gbar(1,3)+gbar(3,1))*(dudx(i,j,k,3)+dudz(i,j,k,1)) +    &
            0.5d0*(gbar(2,3)+gbar(3,2))*(dudy(i,j,k,3)+dudz(i,j,k,2)))/gbarden 

          pamodelt = -(gbar(1,4)*dudx(i,j,k,4)+gbar(2,4)*dudy(i,j,k,4) + &
                       gbar(3,4)*dudz(i,j,k,4))/gbartden

          if(pamodelv<0.0d0) then
            pamodel2 = 0.0d0
          else
            pamodel2 = pamodelv*pamodelv
          endif
          if(gbarden<1.0d-12) then
            tausgs(i,j,k,1:3,1:3) = 0.0d0
          else
            tausgs(i,j,k,1:3,1:3) = 8.0d0*deltabar2(k)*pamodel2*gbar(1:3,1:3)/gbarden
          endif

          ! MGM
          if(pamodelv<0.0d0 .or. pamodelt<0.0d0) then
            pamodel = 0.0d0
          else
            pamodel = pamodelv*pamodelt
          endif
          if(gbartden<1.0d-12) then
            tausgs(i,j,k,4,1:3) = 0.0d0
          else
            tausgs(i,j,k,4,1:3) = 2.0d0*deltabar2(k)/Pr*pamodel*gbar(1:3,4)/gbartden
          endif

          !! Smag grad model
          !if(pamodelt<0.0d0 .or. gbartden<1.0d-12) then
          !  tausgs(i,j,k,4,1:3) = 0.0d0
          !else
          !  tausgs(i,j,k,4,1:3) = 0.075d0*deltabar2(k)*strainmag(i,j,k,1)*strainmag(i,j,k,2)*gbar(1:3,4)/gbartden
          !endif
        enddo
       enddo
      enddo

    ! debug
    !WRITE(myrank+100,*) 'gm-1', maxval(dabs(tausgs(:,:,:,1,1))), minval(dabs(tausgs(:,:,:,1,1)))
    !WRITE(myrank+100,*) 'gm-2', maxval(dabs(tausgs(:,:,:,1,2))), minval(dabs(tausgs(:,:,:,1,2)))
    !WRITE(myrank+100,*) 'gm-3', maxval(dabs(tausgs(:,:,:,1,3))), minval(dabs(tausgs(:,:,:,1,3)))
    !WRITE(myrank+100,*) 'mgm1', maxval(dabs(tausgs(:,:,:,4,1))), minval(dabs(tausgs(:,:,:,4,1)))
    !WRITE(myrank+100,*) 'mgm2', maxval(dabs(tausgs(:,:,:,4,2))), minval(dabs(tausgs(:,:,:,4,2)))
    !WRITE(myrank+100,*) 'mgm3', maxval(dabs(tausgs(:,:,:,4,3))), minval(dabs(tausgs(:,:,:,4,3)))
    !WRITE(myrank+100,*) '---'

      ! restore w derivatives
      dudx(:,:,:,3) = uxtemp(:,:,:,1)
      dudy(:,:,:,3) = uxtemp(:,:,:,2)
      dudz(:,:,:,3) = uxtemp(:,:,:,3)

  END SUBROUTINE gradmodel

  SUBROUTINE gradkernel(ifilt)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifilt
    INTEGER :: i,j,k, ind,jnd
    REAL(DP) :: gbar(3,4),gbarden,gbartden,pamodelv,pamodelt,pamodel,pamodel2,pgrmod

    ievm = 0
    !call magstrainrate

    ! store w derivatives
    uxtemp(:,:,:,1) = dudx(:,:,:,3)
    uxtemp(:,:,:,2) = dudy(:,:,:,3)
    uxtemp(:,:,:,3) = dudz(:,:,:,3)
    ! copy wp derivatives
    dudx(:,:,:,3) = dwpdx(:,:,:)
    dudy(:,:,:,3) = dwpdy(:,:,:)
    dudz(:,:,:,3) = dwpdz(:,:,:)

    do k = sz, ez
     do j = sy, ey
      do i = sx, ex
        ! Porte-Agel model
        do ind=1,3
         do jnd=1,4
          gbar(ind,jnd) = delx2*dudx(i,j,k,ind)*dudx(i,j,k,jnd) + &
                          dely2*dudy(i,j,k,ind)*dudy(i,j,k,jnd) + &
                          delz2(k)*dudz(i,j,k,ind)*dudz(i,j,k,jnd)
         enddo
        enddo
        gbarden = gbar(1,1)+gbar(2,2)+gbar(3,3)
        gbartden = dsqrt(sum(gbar(:,4)**2))

        pamodelv = -(gbar(1,1)*dudx(i,j,k,1) + gbar(2,2)*dudy(i,j,k,2) + &
                gbar(3,3)*dudz(i,j,k,3) +                                &
          0.5d0*(gbar(1,2)+gbar(2,1))*(dudx(i,j,k,2)+dudy(i,j,k,1)) +    &
          0.5d0*(gbar(1,3)+gbar(3,1))*(dudx(i,j,k,3)+dudz(i,j,k,1)) +    &
          0.5d0*(gbar(2,3)+gbar(3,2))*(dudy(i,j,k,3)+dudz(i,j,k,2)))/gbarden 

        pamodelt = -(gbar(1,4)*dudx(i,j,k,4)+gbar(2,4)*dudy(i,j,k,4) + &
                     gbar(3,4)*dudz(i,j,k,4))/gbartden

        ! clipping
        !if(pamodelv<0.0d0) then
        !  pamodel2 = 0.0d0
        !else
        !  pamodel2 = pamodelv*pamodelv
        !endif
        ! deferred clipping
        pamodel2 = pamodelv*pamodelv
        if(ifilt==0) then     ! update hp only for grid-filtered quantities
          if(pamodelv<0.0d0) then
            hp(i,j,k) = 0
          else
            hp(i,j,k) = 1
          endif
        endif
        if(gbarden<1.0d-12) then
          tausgs(i,j,k,1:3,1:3) = 0.0d0
        else
          tausgs(i,j,k,1:3,1:3) = 8.0d0*deltabar2(k)*pamodel2*gbar(1:3,1:3)/gbarden
        endif

        !! MGM
        !! clipping
        !if(pamodelv<0.0d0 .or. pamodelt<0.0d0) then
        !  pamodel = 0.0d0
        !else
        !  pamodel = pamodelv*pamodelt
        !endif
        ! deferred clipping
        pamodel = pamodelv*pamodelt
        if(ifilt==0) then     ! update hp only for grid-filtered quantities
          if(pamodelv<0.0d0 .or. pamodelt<0.0d0) then
            hprho(i,j,k) = 0
          else
            hprho(i,j,k) = 1
          endif
        endif
        if(gbartden<1.0d-12) then
          tausgs(i,j,k,4,1:3) = 0.0d0
        else
          tausgs(i,j,k,4,1:3) = 2.0d0*deltabar2(k)/Pr*pamodel*gbar(1:3,4)/gbartden
        endif

        !! Smag grad model
        !if(pamodelt<0.0d0 .or. gbartden<1.0d-12) then
        !  tausgs(i,j,k,4,1:3) = 0.0d0
        !else
        !  tausgs(i,j,k,4,1:3) = 0.075d0*deltabar2(k)*strainmag(i,j,k,1)*strainmag(i,j,k,2)*gbar(1:3,4)/gbartden
        !endif

        ! for global dynamic procedure
        prodg(i,j,k) = 8.0d0*deltabar2(k)*pamodel2*(-pamodelv)
        prodrhog(i,j,k) = 2.0d0*deltabar2(k)/Pr*pamodel*(-pamodelt)
      enddo
     enddo
    enddo

    ! restore w derivatives
    dudx(:,:,:,3) = uxtemp(:,:,:,1)
    dudy(:,:,:,3) = uxtemp(:,:,:,2)
    dudz(:,:,:,3) = uxtemp(:,:,:,3)

  END SUBROUTINE gradkernel

  SUBROUTINE ftran_filtered

    IMPLICIT NONE
    INTEGER :: var

    do var = 1, nvar
      IF(var==3) then
        call p3dfft_ftran_r2c(wp(:,:,:),wph(:,:,:),fstr)
        wph(:,:,:) = itfac*wph(:,:,:)
      ELSE
        call p3dfft_ftran_r2c(u(:,:,:,var),uh(:,:,:,var),fstr)
        uh(:,:,:,var) = itfac*uh(:,:,:,var)
      ENDIF
    enddo

    call dealiasuh

  END SUBROUTINE ftran_filtered

  SUBROUTINE filter3d(uin,uout)

    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(sx:ex,sy:ey,sz:ez) :: uin
    REAL(DP), INTENT(OUT), DIMENSION(sx:ex,sy:ey,sz:ez) :: uout
    INTEGER :: i, j, k, sendcnt, recvcnt, mpierr, klow, kupp
    REAL(DP) :: df1, df2, fac1, fac2, fac3

    ! filter in x direction
    ! no decomposition in x direction
    do i = sx, ex
     IF(i==sx) then
       uout(i,:,:) = sixth*uin(nxg-1,:,:)+twothird*uin(i,:,:)+sixth*uin(i+1,:,:)
     ELSEIF(i==ex) then
       uout(i,:,:) = sixth*uin(i-1,:,:)+twothird*uin(i,:,:)+sixth*uin(2,:,:)
     ELSE
       uout(i,:,:) = sixth*uin(i-1,:,:)+twothird*uin(i,:,:)+sixth*uin(i+1,:,:)
     ENDIF
    enddo

    ! filter in y direction
    ! no decomposition in y direction - for the moment
    do j = sy, ey
     IF(j==sy) then
       uin(:,j,:) = sixth*uout(:,nyg-1,:)+twothird*uout(:,j,:)+sixth*uout(:,j+1,:)
     ELSEIF(j==ey) then
       uin(:,j,:) = sixth*uout(:,j-1,:)+twothird*uout(:,j,:)+sixth*uout(:,2,:)
     ELSE
       uin(:,j,:) = sixth*uout(:,j-1,:)+twothird*uout(:,j,:)+sixth*uout(:,j+1,:)
     ENDIF
    enddo

    ! filter in z direction
    ! send zm1
    sendbuffilt(:,:) = uin(:,:,ez)
    sendcnt = nx*ny; recvcnt = sendcnt
    call MPI_SENDRECV(sendbuffilt,sendcnt,MPIRP,MPIDEST,0,recvbuffilt,recvcnt,MPIRP,MPISRC,0,MPICOMM2D,mpistat,mpierr)
    ubufzm1(:,:,1) = recvbuffilt(:,:)

    ! send zp1
    sendbuffilt(:,:) = uin(:,:,sz)
    sendcnt = nx*ny; recvcnt = sendcnt
    call MPI_SENDRECV(sendbuffilt,sendcnt,MPIRP,MPISRC,0,recvbuffilt,recvcnt,MPIRP,MPIDEST,0,MPICOMM2D,mpistat,mpierr)
    ubufzp1(:,:,1) = recvbuffilt(:,:)

    ! boundary conditions
    IF(sz==1) then
      ubufzm1(:,:,1) = 2.0d0*uin(:,:,1)-uin(:,:,2)      ! Not staggered
    ENDIF
    IF(ez==nzg) then
      ubufzp1(:,:,1) = 2.0d0*uin(:,:,nzg)-uin(:,:,nzg-1)      ! Not staggered
    ENDIF

    IF(sz==1) then
      uout(:,:,sz) = uin(:,:,sz)
      klow = sz+1
    ELSE
      klow = sz
    ENDIF
    IF(ez==nzg) then
      uout(:,:,ez) = uin(:,:,ez)
      kupp = ez-1
    ELSE
      kupp = ez
    ENDIF

    do k = klow, kupp
     df1 = z(k)-z(k-1); df2 = z(k+1)-z(k)
     IF(df1<df2) then
       fac1 = sixth*(1.0d0-(df2-df1)/(df2+df1))
       fac2 = twothird+third*(1.0d0-df1/df2)
       fac3 = third*df1*df1/(df2*(df1+df2))
     ELSE
       fac3 = sixth*(1.0d0-(df1-df2)/(df2+df1))
       fac2 = twothird+third*(1.0d0-df2/df1)
       fac1 = third*df2*df2/(df1*(df1+df2))
     ENDIF
     IF(k==sz) then
       uout(:,:,k) = fac1*ubufzm1(:,:,1)+fac2*uin(:,:,k)+fac3*uin(:,:,k+1)
     ELSEIF(k==ez) then
       uout(:,:,k) = fac1*uin(:,:,k-1)+fac2*uin(:,:,k)+fac3*ubufzp1(:,:,1)
     ELSE
       uout(:,:,k) = fac1*uin(:,:,k-1)+fac2*uin(:,:,k)+fac3*uin(:,:,k+1)
     ENDIF
    enddo

  END SUBROUTINE filter3d

  SUBROUTINE globdynamic_gradmodel

    IMPLICIT NONE
    INTEGER :: i, j, k, mpierr, ifilt
    REAL(DP) :: term2,coef,coefrho,coeffacl(4),coeffac(4)

    ifilt = 0  ! grid-filtered quantities
    call gradkernel(ifilt)

    idyn = idyn + 1
    if(idyn==DYNFREQ) then

      idyn = 0

      ! store tausgs
      tautemp(:,:,:,1,1) = tausgs(:,:,:,1,1)
      tautemp(:,:,:,1,2) = tausgs(:,:,:,1,2)
      tautemp(:,:,:,1,3) = tausgs(:,:,:,1,3)
      tautemp(:,:,:,2,2) = tausgs(:,:,:,2,2)
      tautemp(:,:,:,2,3) = tausgs(:,:,:,2,3)
      tautemp(:,:,:,3,3) = tausgs(:,:,:,3,3)
      tautemp(:,:,:,2,1) = tausgs(:,:,:,4,1)
      tautemp(:,:,:,3,1) = tausgs(:,:,:,4,2)
      tautemp(:,:,:,3,2) = tausgs(:,:,:,4,3)
  
      ! first parts of numerator and denominator terms
      uin(:,:,:) = dudx(:,:,:,1)**2 +  dudx(:,:,:,2)**2 +  dwpdx(:,:,:)**2 & 
                 + dudy(:,:,:,1)**2 +  dudy(:,:,:,2)**2 +  dwpdy(:,:,:)**2 & 
                 + dudz(:,:,:,1)**2 +  dudz(:,:,:,2)**2 +  dwpdz(:,:,:)**2
      call filter3d(uin(:,:,:),numg(:,:,:))
      call filter3d(prodg(:,:,:),deng(:,:,:))
      uin(:,:,:) = wp(:,:,:)*u(:,:,:,4); call filter3d(uin(:,:,:), utemp(:,:,:,1))
      numg(:,:,:) = numg(:,:,:)/Re - Ri*utemp(:,:,:,1)
  
      uin(:,:,:) = dudx(:,:,:,4)**2 +  dudy(:,:,:,4)**2 +  dudz(:,:,:,4)**2
      call filter3d(uin(:,:,:),numrhog(:,:,:))
      call filter3d(prodrhog(:,:,:),denrhog(:,:,:))
  
      ! store u and uh
      utemp(:,:,:,1) = u(:,:,:,1)
      utemp(:,:,:,2) = u(:,:,:,2)
      utemp(:,:,:,3) = wp(:,:,:)
      utemp(:,:,:,4) = u(:,:,:,4)
      uhtemp(:,:,:,1) = uh(:,:,:,1)
      uhtemp(:,:,:,2) = uh(:,:,:,2)
      uhtemp(:,:,:,3) = wph(:,:,:)
      uhtemp(:,:,:,4) = uh(:,:,:,4)
  
      ! test-filter individual velocities and scalar
      uin(:,:,:) = u(:,:,:,1); call filter3d(uin(:,:,:),u(:,:,:,1))
      uin(:,:,:) = u(:,:,:,2); call filter3d(uin(:,:,:),u(:,:,:,2))
      uin(:,:,:) = wp(:,:,:);  call filter3d(uin(:,:,:),wp(:,:,:))
      uin(:,:,:) = u(:,:,:,4); call filter3d(uin(:,:,:),u(:,:,:,4))
  
      ! compute derivatives
      call ftran_filtered
      iwpov = 0; call exchangez(iwpov)
      !call vinterp
      !call vinterph
      call bc
      call xyderivs
      call zderivs
  
      ifilt = 1  ! test-filtered quantities
      call gradkernel(ifilt)
  
      ! second parts of numerator and denominator terms; global coefficients
      coeffacl = 0.0d0
      DO k = sz, ez
       DO j = sy, ey
        DO i = sx, ex
          term2 = dudx(i,j,k,1)**2 +  dudx(i,j,k,2)**2 +  dwpdx(i,j,k)**2 & 
                + dudy(i,j,k,1)**2 +  dudy(i,j,k,2)**2 +  dwpdy(i,j,k)**2 & 
                + dudz(i,j,k,1)**2 +  dudz(i,j,k,2)**2 +  dwpdz(i,j,k)**2
          numg(i,j,k) = numg(i,j,k) - term2/Re + Ri*wp(i,j,k)*u(i,j,k,4)
          deng(i,j,k) = deng(i,j,k) - 4.0d0*prodg(i,j,k)
  
          term2 = dudx(i,j,k,4)**2 +  dudy(i,j,k,4)**2 +  dudz(i,j,k,4)**2 
          numrhog(i,j,k) = (numrhog(i,j,k) - term2)/(Re*Pr)
          denrhog(i,j,k) = denrhog(i,j,k) - 4.0d0*prodrhog(i,j,k)
  
          coeffacl(1) = coeffacl(1) + numg(i,j,k)*vol(k)   ! num
          coeffacl(2) = coeffacl(2) + deng(i,j,k)*vol(k)   ! den
          coeffacl(3) = coeffacl(3) + numrhog(i,j,k)*vol(k)    ! numrho
          coeffacl(4) = coeffacl(4) + denrhog(i,j,k)*vol(k)    ! denrho
        ENDDO
       ENDDO
      ENDDO
  
      ! compute global coefficients
      call MPI_ALLREDUCE(coeffacl,coeffac,4,MPIRP,MPI_SUM,MPICOMM2D,mpierr)
      coef = 0.0d0; coefrho = 0.0d0
      IF(dabs(coeffac(2))>1.0d-12) coef = coeffac(1)/coeffac(2)
      IF(coef<1.0d-12) coef = 1.0d0
      IF(dabs(coeffac(4))>1.0d-12) coefrho = coeffac(3)/coeffac(4)
      IF(coefrho<1.0d-12) coefrho = 1.0d0
  
      ! restore u and uh
      u(:,:,:,1) = utemp(:,:,:,1)
      u(:,:,:,2) = utemp(:,:,:,2)
      wp(:,:,:) = utemp(:,:,:,3)
      u(:,:,:,4) = utemp(:,:,:,4)
      uh(:,:,:,1) = uhtemp(:,:,:,1)
      uh(:,:,:,2) = uhtemp(:,:,:,2)
      wph(:,:,:) = uhtemp(:,:,:,3)
      uh(:,:,:,4) = uhtemp(:,:,:,4)
  
      globcoef = coef; globcoefrho = coefrho

    endif

    write(223,*) t, globcoef, globcoefrho

    ! copmute tausgs
    tausgs(:,:,:,1,1) = globcoef*tautemp(:,:,:,1,1)*hp(:,:,:)
    tausgs(:,:,:,1,2) = globcoef*tautemp(:,:,:,1,2)*hp(:,:,:)
    tausgs(:,:,:,1,3) = globcoef*tautemp(:,:,:,1,3)*hp(:,:,:)
    tausgs(:,:,:,2,2) = globcoef*tautemp(:,:,:,2,2)*hp(:,:,:)
    tausgs(:,:,:,2,3) = globcoef*tautemp(:,:,:,2,3)*hp(:,:,:)
    tausgs(:,:,:,3,3) = globcoef*tautemp(:,:,:,3,3)*hp(:,:,:)
    tausgs(:,:,:,2,1) = tausgs(:,:,:,1,2)
    tausgs(:,:,:,3,1) = tausgs(:,:,:,1,3)
    tausgs(:,:,:,3,2) = tausgs(:,:,:,2,3)

    tausgs(:,:,:,4,1) = globcoefrho*tautemp(:,:,:,2,1)*hprho(:,:,:)
    tausgs(:,:,:,4,2) = globcoefrho*tautemp(:,:,:,3,1)*hprho(:,:,:)
    tausgs(:,:,:,4,3) = globcoefrho*tautemp(:,:,:,3,2)*hprho(:,:,:)
    
  END SUBROUTINE globdynamic_gradmodel

  SUBROUTINE dynamic_gradmodel

    IMPLICIT NONE
    INTEGER :: i, j, k, mpierr, ifilt
    REAL(DP) :: lm, mm, lmrho, mmrho

    ifilt = 0  ! grid-filtered quantities
    call gradkernel(ifilt)
    ! store tausgs, u and uh
    tautemp(:,:,:,1,1) = tausgs(:,:,:,1,1)
    tautemp(:,:,:,1,2) = tausgs(:,:,:,1,2)
    tautemp(:,:,:,1,3) = tausgs(:,:,:,1,3)
    tautemp(:,:,:,2,2) = tausgs(:,:,:,2,2)
    tautemp(:,:,:,2,3) = tausgs(:,:,:,2,3)
    tautemp(:,:,:,3,3) = tausgs(:,:,:,3,3)
    tautemp(:,:,:,2,1) = tausgs(:,:,:,4,1)
    tautemp(:,:,:,3,1) = tausgs(:,:,:,4,2)
    tautemp(:,:,:,3,2) = tausgs(:,:,:,4,3)

    idyn = idyn+1
    IF(idyn==1) then   ! compute coefficients every 10 time steps

    idyn = 0

    ! store tausgs, u and uh
    utemp(:,:,:,1) = u(:,:,:,1)
    utemp(:,:,:,2) = u(:,:,:,2)
    utemp(:,:,:,3) = wp(:,:,:)
    utemp(:,:,:,4) = u(:,:,:,4)
    uhtemp(:,:,:,1) = uh(:,:,:,1)
    uhtemp(:,:,:,2) = uh(:,:,:,2)
    uhtemp(:,:,:,3) = wph(:,:,:)
    uhtemp(:,:,:,4) = uh(:,:,:,4)

    uin(:,:,:) = u(:,:,:,1)*u(:,:,:,1); call filter3d(uin(:,:,:),lterm(:,:,:,1,1))
    uin(:,:,:) = u(:,:,:,1)*u(:,:,:,2); call filter3d(uin(:,:,:),lterm(:,:,:,1,2))
    uin(:,:,:) = u(:,:,:,1)*wp(:,:,:);  call filter3d(uin(:,:,:),lterm(:,:,:,1,3))
    uin(:,:,:) = u(:,:,:,2)*u(:,:,:,2); call filter3d(uin(:,:,:),lterm(:,:,:,2,2))
    uin(:,:,:) = u(:,:,:,2)*wp(:,:,:);  call filter3d(uin(:,:,:),lterm(:,:,:,2,3))
    uin(:,:,:) = wp(:,:,:)*wp(:,:,:);   call filter3d(uin(:,:,:),lterm(:,:,:,3,3))
    uin(:,:,:) = u(:,:,:,1)*u(:,:,:,4); call filter3d(uin(:,:,:),lterm(:,:,:,2,1))
    uin(:,:,:) = u(:,:,:,2)*u(:,:,:,4); call filter3d(uin(:,:,:),lterm(:,:,:,3,1))
    uin(:,:,:) = wp(:,:,:)*u(:,:,:,4);  call filter3d(uin(:,:,:),lterm(:,:,:,3,2))

    uin(:,:,:) = tausgs(:,:,:,1,1); call filter3d(uin(:,:,:),mterm(:,:,:,1,1))
    uin(:,:,:) = tausgs(:,:,:,1,2); call filter3d(uin(:,:,:),mterm(:,:,:,1,2))
    uin(:,:,:) = tausgs(:,:,:,1,3); call filter3d(uin(:,:,:),mterm(:,:,:,1,3))
    uin(:,:,:) = tausgs(:,:,:,2,2); call filter3d(uin(:,:,:),mterm(:,:,:,2,2))
    uin(:,:,:) = tausgs(:,:,:,2,3); call filter3d(uin(:,:,:),mterm(:,:,:,2,3))
    uin(:,:,:) = tausgs(:,:,:,3,3); call filter3d(uin(:,:,:),mterm(:,:,:,3,3))
    uin(:,:,:) = tausgs(:,:,:,4,1); call filter3d(uin(:,:,:),mterm(:,:,:,2,1))
    uin(:,:,:) = tausgs(:,:,:,4,2); call filter3d(uin(:,:,:),mterm(:,:,:,3,1))
    uin(:,:,:) = tausgs(:,:,:,4,3); call filter3d(uin(:,:,:),mterm(:,:,:,3,2))

    uin(:,:,:) = u(:,:,:,1); call filter3d(uin(:,:,:),u(:,:,:,1))
    uin(:,:,:) = u(:,:,:,2); call filter3d(uin(:,:,:),u(:,:,:,2))
    uin(:,:,:) = wp(:,:,:);  call filter3d(uin(:,:,:),wp(:,:,:))
    uin(:,:,:) = u(:,:,:,4); call filter3d(uin(:,:,:),u(:,:,:,4))

    call ftran_filtered

    ! prepare for stats or next time step
    iwpov = 0; call exchangez(iwpov)
    !call vinterp
    !call vinterph
    call bc
    call xyderivs
    call zderivs
    ifilt = 1  ! test-filtered quantities
    call gradkernel(ifilt)

    lterm(:,:,:,1,1) = lterm(:,:,:,1,1) - u(:,:,:,1)*u(:,:,:,1)
    lterm(:,:,:,1,2) = lterm(:,:,:,1,2) - u(:,:,:,1)*u(:,:,:,2)
    lterm(:,:,:,1,3) = lterm(:,:,:,1,3) - u(:,:,:,1)*wp(:,:,:)
    lterm(:,:,:,2,2) = lterm(:,:,:,2,2) - u(:,:,:,2)*u(:,:,:,2)
    lterm(:,:,:,2,3) = lterm(:,:,:,2,3) - u(:,:,:,2)*wp(:,:,:)
    lterm(:,:,:,3,3) = lterm(:,:,:,3,3) - wp(:,:,:)*wp(:,:,:)
    lterm(:,:,:,2,1) = lterm(:,:,:,2,1) - u(:,:,:,1)*u(:,:,:,4)
    lterm(:,:,:,3,1) = lterm(:,:,:,3,1) - u(:,:,:,2)*u(:,:,:,4)
    lterm(:,:,:,3,2) = lterm(:,:,:,3,2) - wp(:,:,:)*u(:,:,:,4)

    mterm(:,:,:,1,1) = 4.0d0*tausgs(:,:,:,1,1) - mterm(:,:,:,1,1)
    mterm(:,:,:,1,2) = 4.0d0*tausgs(:,:,:,1,2) - mterm(:,:,:,1,2)
    mterm(:,:,:,1,3) = 4.0d0*tausgs(:,:,:,1,3) - mterm(:,:,:,1,3)
    mterm(:,:,:,2,2) = 4.0d0*tausgs(:,:,:,2,2) - mterm(:,:,:,2,2)
    mterm(:,:,:,2,3) = 4.0d0*tausgs(:,:,:,2,3) - mterm(:,:,:,2,3)
    mterm(:,:,:,3,3) = 4.0d0*tausgs(:,:,:,3,3) - mterm(:,:,:,3,3)
    mterm(:,:,:,2,1) = 4.0d0*tausgs(:,:,:,4,1) - mterm(:,:,:,2,1)
    mterm(:,:,:,3,1) = 4.0d0*tausgs(:,:,:,4,2) - mterm(:,:,:,3,1)
    mterm(:,:,:,3,2) = 4.0d0*tausgs(:,:,:,4,3) - mterm(:,:,:,3,2)

    numl = 0.0d0; denl = 0.0d0; numrhol = 0.0d0; denrhol = 0.0d0
    DO k = sz, ez
     DO j = sy, ey
      DO i = sx, ex
        lm = lterm(i,j,k,1,1)*mterm(i,j,k,1,1) + lterm(i,j,k,2,2)*mterm(i,j,k,2,2) + lterm(i,j,k,3,3)*mterm(i,j,k,3,3) &
           + 2.0d0*(lterm(i,j,k,1,2)*mterm(i,j,k,1,2) + lterm(i,j,k,1,3)*mterm(i,j,k,1,3) + lterm(i,j,k,2,3)*mterm(i,j,k,2,3))

        mm = mterm(i,j,k,1,1)**2 + mterm(i,j,k,2,2)**2 + mterm(i,j,k,3,3)**2 &
           + 2.0d0*(mterm(i,j,k,1,2)**2 + mterm(i,j,k,1,3)**2 + mterm(i,j,k,2,3)**2)

        lmrho = lterm(i,j,k,2,1)*mterm(i,j,k,2,1) + lterm(i,j,k,3,1)*mterm(i,j,k,3,1) + lterm(i,j,k,3,2)*mterm(i,j,k,3,2)
        mmrho = mterm(i,j,k,2,1)**2 + mterm(i,j,k,3,1)**2 + mterm(i,j,k,3,2)**2

        numl(k) = numl(k) + lm
        denl(k) = denl(k) + mm
        numrhol(k) = numrhol(k) + lmrho
        denrhol(k) = denrhol(k) + mmrho
      ENDDO
     ENDDO
    ENDDO

    num = numl; den = denl; numrho = numrhol; denrho = denrhol
    !call MPI_ALLREDUCE(numl,num,nzg,MPIRP,MPI_SUM,MPICOMM2D,mpierr)
    !call MPI_ALLREDUCE(denl,den,nzg,MPIRP,MPI_SUM,MPICOMM2D,mpierr)
    !call MPI_ALLREDUCE(numrhol,numrho,nzg,MPIRP,MPI_SUM,MPICOMM2D,mpierr)
    !call MPI_ALLREDUCE(denrhol,denrho,nzg,MPIRP,MPI_SUM,MPICOMM2D,mpierr)

    !write(MYRANK+100,'(8(e12.5,1x))') maxval(num), minval(num), maxval(den), minval(den), maxval(numrho), minval(numrho), maxval(denrho), minval(denrho)

    DO k = sz, ez
      if(dabs(den(k))>1.0d-12) then
        coef(k) = num(k)/den(k)
      else
        coef(k) = 1.0d0
      endif
      if(coef(k)<1.0d-12) coef(k) = 1.0d0
      if(dabs(denrho(k))>1.0d-12) then
        coefrho(k) = numrho(k)/denrho(k)
      endif
      if(coefrho(k)<1.0d-12) coefrho(k) = 1.0d0
    ENDDO

    ENDIF

    DO k = sz, ez
      tausgs(:,:,k,1,1) = coef(k)*tautemp(:,:,k,1,1)*hp(:,:,k)
      tausgs(:,:,k,1,2) = coef(k)*tautemp(:,:,k,1,2)*hp(:,:,k)
      tausgs(:,:,k,1,3) = coef(k)*tautemp(:,:,k,1,3)*hp(:,:,k)
      tausgs(:,:,k,2,2) = coef(k)*tautemp(:,:,k,2,2)*hp(:,:,k)
      tausgs(:,:,k,2,3) = coef(k)*tautemp(:,:,k,2,3)*hp(:,:,k)
      tausgs(:,:,k,3,3) = coef(k)*tautemp(:,:,k,3,3)*hp(:,:,k)
      tausgs(:,:,k,2,1) = tausgs(:,:,k,1,2)
      tausgs(:,:,k,3,1) = tausgs(:,:,k,1,3)
      tausgs(:,:,k,3,2) = tausgs(:,:,k,2,3)

      tausgs(:,:,k,4,1) = coefrho(k)*tautemp(:,:,k,2,1)*hprho(:,:,k)
      tausgs(:,:,k,4,2) = coefrho(k)*tautemp(:,:,k,3,1)*hprho(:,:,k)
      tausgs(:,:,k,4,3) = coefrho(k)*tautemp(:,:,k,3,2)*hprho(:,:,k)
    ENDDO
    
    ! restore tausgs, u and uh
    u(:,:,:,1) = utemp(:,:,:,1)
    u(:,:,:,2) = utemp(:,:,:,2)
    wp(:,:,:) = utemp(:,:,:,3)
    u(:,:,:,4) = utemp(:,:,:,4)
    uh(:,:,:,1) = uhtemp(:,:,:,1)
    uh(:,:,:,2) = uhtemp(:,:,:,2)
    wph(:,:,:) = uhtemp(:,:,:,3)
    uh(:,:,:,4) = uhtemp(:,:,:,4)

  END SUBROUTINE dynamic_gradmodel

!  SUBROUTINE dynamic_smagorinsky
!
!    IMPLICIT NONE
!
!    numerl = 0.0d0; denoml = 0.0d0
!    do k=sz,ez
!     do j=sy,ey
!      do i=sx,ex
!        numerl(k,1) = numerl(k,1) + matmul(Lij-Nij,Mij)
!        denoml(k,1) = denoml(k,1) + matmul(Mij,Mij)
!        numerl(k,2) = numerl(k,2) + matmul(LTi,MTi)
!        denoml(k,2) = denoml(k,2) + matmul(MTi,MTi)
!      enddo
!     enddo
!    enddo
!    call MPI_ALLREDUCE(numerl,numer,2*nzg,MPIRP,MPI_SUM,MPICOMM2D,mpierr)
!    call MPI_ALLREDUCE(denoml,denom,2*nzg,MPIRP,MPI_SUM,MPICOMM2D,mpierr)
!    do k=sz,ez
!      if(dabs(denom(k,1))>1.0d-12) then
!        cs(k,1) = -0.5d0*numer(k,1)/denom(k,1)
!      else
!        cs(k,1) = 0.0d0
!      endif
!      if(dabs(denom(k,2))>1.0d-12) then
!        cs(k,2) = -0.5d0*numer(k,2)/denom(k,2)
!      else
!        cs(k,2) = 0.0d0
!      endif
!
!      eddyvisc(:,:,k,1) = cs(k,1)*deltabar2(k)*strainmag(:,:,k,1)
!      eddyvisc(:,:,k,2) = cs(k,2)*deltabar2(k)*strainmag(:,:,k,1)
!
!    do j=sy,ey
!      do i=sx,ex
!      enddo
!     enddo
!    enddo
!
!  END SUBROUTINE dynamic_smagorinsky

  SUBROUTINE sigmamodel_methodb

    IMPLICIT NONE
    INTEGER :: i,j,k, numneg1, numneg2, numneg3
    REAL(DP) :: invar1, invar2, invar3, alpha1, alpha2, alpha3, trgtg2, fac1, fac2, lam1, lam2, lam3

    ievm = 1
    cs = 1.35d0; Prt = 0.7d0
    numneg1 = 0; numneg2 = 0; numneg3 = 0

    !call subzderivs

    do k = sz, ez
     do j = sy, ey
      do i = sx, ex
        g(1,1)=dudx(i,j,k,1); g(1,2)=dudy(i,j,k,1); g(1,3)=dudz(i,j,k,1)!*zetaz(k)
        g(2,1)=dudx(i,j,k,2); g(2,2)=dudy(i,j,k,2); g(2,3)=dudz(i,j,k,2)!*zetaz(k)
        g(3,1)=dwpdx(i,j,k);  g(3,2)=dwpdy(i,j,k);  g(3,3)=dwpdz(i,j,k)!*zetaz(k)

        !g(1,1)=1.0d0; g(1,2)=2.0d0; g(1,3)=3.0d0
        !g(2,1)=4.0d0; g(2,2)=5.0d0; g(2,3)=6.0d0
        !g(3,1)=7.0d0; g(3,2)=8.0d0; g(3,3)=10.0d0

        gt(1,1)=g(1,1); gt(1,2)=g(2,1); gt(1,3)=g(3,1)
        gt(2,1)=g(1,2); gt(2,2)=g(2,2); gt(2,3)=g(3,2)
        gt(3,1)=g(1,3); gt(3,2)=g(2,3); gt(3,3)=g(3,3)
        gtg = matmul(gt,g)

        !print *, 'g: ', g
        !print *, 'gt: ', gt
        !print *, 'gtg: ', gtg

        invar1 = gtg(1,1)+gtg(2,2)+gtg(3,3)
        trgtg2 = gtg(1,1)**2+gtg(2,2)**2+gtg(3,3)**2+2.0d0*(gtg(1,2)*gtg(2,1)+gtg(1,3)*gtg(3,1)+gtg(2,3)*gtg(3,2))
        !invar2 = 0.5d0*invar1**2-trgtg2
        invar2 = gtg(1,1)*gtg(2,2)+gtg(1,1)*gtg(3,3)+gtg(2,2)*gtg(3,3)-gtg(1,2)**2-gtg(1,3)**2-gtg(2,3)**2
        invar3 = gtg(1,1)*(gtg(2,2)*gtg(3,3)-gtg(2,3)*gtg(3,2)) &
               - gtg(2,1)*(gtg(1,2)*gtg(3,3)-gtg(3,2)*gtg(1,3)) &
               + gtg(3,1)*(gtg(1,2)*gtg(2,3)-gtg(2,2)*gtg(1,3))

        !print *, 'invars: ', invar1, invar2, invar3

        alpha1 = ninth*invar1**2-third*invar2
        alpha2 = tseventh*invar1**3-sixth*invar1*invar2+half*invar3
        fac1 = alpha1**1.5; fac2 = alpha2/fac1
        IF(fac1 < 1.0d-10 .or. dabs(fac2)>1.0d0) then
          alpha3=0.0d0
        ELSE
          alpha3 = third*dacos(fac2)
        ENDIF
        ! first check: 0 < alpha3 < pi/3
        IF(alpha3<0.0d0 .or. alpha3>piby3) then
          print *, 'alpha3 outside [0, pi/3]', i, j, k, alpha3
          stop
        ENDIF

        IF(alpha1<0.0d0) then
          write(*,*) 'alpha1: ', alpha1
          alpha1 = 0.0d0
        ELSE
          alpha1 = 2.0d0*dsqrt(alpha1)
        ENDIF

        lam1 = third*invar1+alpha1*dcos(alpha3)
        lam2 = third*invar1-alpha1*dcos(piby3+alpha3)
        lam3 = third*invar1-alpha1*dcos(piby3-alpha3)

        IF(dabs(lam1)<1.0d-7) lam1=0.0d0
        IF(dabs(lam2)<1.0d-4*dabs(lam1)) lam2=0.0d0
        IF(dabs(lam3)<1.0d-4*dabs(lam1)) lam3=0.0d0

        IF(lam1<0.0d0) then
          !write(*,*) 'sigma1: ', sigma1
          numneg1=numneg1+1
          sigma1 = 0.0d0
        ELSE
          sigma1 = dsqrt(lam1)
        ENDIF

        IF(lam2<0.0d0) then
          !write(*,*) 'sigma2: ', sigma2
          numneg2=numneg2+1
          sigma2 = 0.0d0
        ELSE
          sigma2 = dsqrt(lam2)
        ENDIF
        !sigma3 = invar1-sigma1**2-sigma2**2
        IF(lam3<0.0d0) then
          !write(*,2) 'sigma3: ', sigma3, invar1, alpha1, alpha3, piby3-alpha3,
          !dcos(piby3-alpha3)
          !write(*,1) sigma1, sigma2, sigma3
          !write(*,*) 'g: ', g
          !write(*,*) 'sigma3: ', sigma3
          numneg3=numneg3+1
          sigma3 = 0.0d0
          !stop
        ELSE
          sigma3 = dsqrt(lam3)
        ENDIF
        !write(*,1) invar1, invar2,invar3!alpha1,sigma1,sigma2,sigma3
        !write(*,1) gtg, gt, g
        !write(*,1) sigma1, sigma2, sigma3
        IF(sigma1>0.0d0) then
          eddyvisc(i,j,k,1) = cs**2*deltabar2(k)*(sigma3*(sigma1-sigma2)*(sigma2-sigma3)/sigma1**2)
        ELSE
          eddyvisc(i,j,k,1) = 0.0d0
        ENDIF
        eddyvisc(i,j,k,2) = eddyvisc(i,j,k,1)/Prt        
      enddo
     enddo
    enddo

    1 FORMAT(6(1x,e12.5))
    2 FORMAT(a,6(1x,e12.5))

  END SUBROUTINE sigmamodel_methodb

  SUBROUTINE magstrainrate

    IMPLICIT NONE
    INTEGER :: i,j,k

    do k = sz, ez
     do j = sy, ey
       do i = sx, ex
         !strain rate
         tausgs(i,j,k,1,1) = dudx(i,j,k,1)
         tausgs(i,j,k,2,2) = dudy(i,j,k,2)
         tausgs(i,j,k,3,3) = dwpdz(i,j,k)
         tausgs(i,j,k,1,2) = 0.5d0*(dudy(i,j,k,1) + dudx(i,j,k,2))
         tausgs(i,j,k,2,1) = tausgs(i,j,k,1,2)
         tausgs(i,j,k,1,3) = 0.5d0*(dudz(i,j,k,1) + dwpdx(i,j,k))
         tausgs(i,j,k,3,1) = tausgs(i,j,k,1,3)
         tausgs(i,j,k,2,3) = 0.5d0*(dudz(i,j,k,2) + dwpdy(i,j,k))
         tausgs(i,j,k,3,2) = tausgs(i,j,k,2,3)

         ! magnitude of strain rate
         strainmag(i,j,k,1) = dsqrt(2.0d0*(tausgs(i,j,k,1,1)**2 + &
                            tausgs(i,j,k,2,2)**2 + &
                            tausgs(i,j,k,3,3)**2 + &
                         2.0d0*tausgs(i,j,k,1,2)**2 + &
                         2.0d0*tausgs(i,j,k,1,3)**2 + &
                         2.0d0*tausgs(i,j,k,2,3)**2) )

         ! scalar strain
         strainmag(i,j,k,2) = dsqrt(dudx(i,j,k,4)**2.0d0+dudy(i,j,k,4)**2.0d0+dudz(i,j,k,4)**2.0d0)
        end do
      end do
    end do
  !$OMP END PARALLEL DO


  END SUBROUTINE magstrainrate

  SUBROUTINE sgstress

    IMPLICIT NONE
    INTEGER :: i,j,k

    IF(ievm==1) then
     call magstrainrate
     do k = sz, ez
      do j = sy, ey
       do i = sx, ex
         tausgs(i,j,k,1:3,:) = -2.0d0*eddyvisc(i,j,k,1)*tausgs(i,j,k,1:3,:)
         tausgs(i,j,k,4,1) = -eddyvisc(i,j,k,2)*dudx(i,j,k,4)
         tausgs(i,j,k,4,2) = -eddyvisc(i,j,k,2)*dudy(i,j,k,4)
         tausgs(i,j,k,4,3) = -eddyvisc(i,j,k,2)*dudz(i,j,k,4)
       enddo
      enddo
     enddo
    ENDIF
  
    !write(222,*) t, maxval(tausgs), minval(tausgs)
  
    ! Take div of tau for sgsterm
    call taudiv

  END SUBROUTINE sgstress

  SUBROUTINE taudiv

    IMPLICIT NONE

    INTEGER :: i,j,k,var
    REAL(DP) :: xl,xr,xm,fac

    ! dtauds x derivatives
    do var = 1, nvar
      call p3dfft_ftran_r2c(tausgs(:,:,:,var,1),tmph(:,:,:),fstr)
      do i = sxf, exf
        sgstermh(i,:,:,var) = iu*k1(i)*itfac*tmph(i,:,:)
      enddo
    enddo

    ! dtauds y derivatives
    do var = 1, nvar
      call p3dfft_ftran_r2c(tausgs(:,:,:,var,2),tmph(:,:,:),fstr)
      do j = syf, eyf
        sgstermh(:,j,:,var) = sgstermh(:,j,:,var) + iu*k2(j)*itfac*tmph(:,j,:)
      enddo
    enddo

    ! store sgstemph for z derivs to be computed later
    do var = 1, nvar
      call p3dfft_ftran_r2c(tausgs(:,:,:,var,3),tmph(:,:,:),fstr)
      do j = syf, eyf
       do i = sxf, exf
         call deriv1cenh(dsgsdzh(szf:ezf),tmph(i,j,szf:ezf),0,1)
         sgstermh(i,j,:,var) = sgstermh(i,j,:,var) + itfac*dsgsdzh(:)
       enddo
      enddo
    enddo

    ! interpolate in z direction
    tmph(:,:,:) = sgstermh(:,:,:,3)
    DO k = szf, ezf-1
      xl = z(k); xr = z(k+1); xm = zz(k); fac = (xm-xl)/(xr-xl)
      sgstermh(:,:,k,3) = tmph(:,:,k) + fac*(tmph(:,:,k+1)-tmph(:,:,k))
    ENDDO

  END SUBROUTINE taudiv

  SUBROUTINE rhsmom_Euler

    IMPLICIT NONE

    !! already done at the end ofprevious time step
    !call exchangez
    !call vinterp
    !call vinterph
    !call xyderivs
    !call zderivs

    call convect
    call viscous

    call Corriolis
    call buoyancy

    call globdynamic_gradmodel
    call sgstress

    rhsh = uh + dt*(-convh + viscth + viscnh + buoyh - sgstermh - fCorriolh)
    !rhsh(:,:,:,1) = rhsh(:,:,:,1) + dt_pgradh   ! only for channel flow

    ! stored for next time ABCN step
    rhs1h = -convh + viscth + buoyh - sgstermh - fCorriolh

  END SUBROUTINE rhsmom_Euler

  SUBROUTINE solvetridiag(i,j,kst,kend)
  INTEGER, INTENT(IN) :: i,j,kst,kend
  INTEGER :: k
  REAL(DP) :: den
  
  !
  ! am(1) cm(1) 0 ----------0
  ! bm(1) am(2) cm(2) 0 -----0
  ! 0     bm(2) am(3) cm(3)---0
  !                am(n-1)cm(n-1) 
  !                bm(n-1)am(n) 
  !
  
    cp(kst) = cm(kst)/am(kst)
    dph(kst) = rhsph(i,j,kst)/am(kst)
    DO k = kst+1, kend
      den = am(k) - cp(k-1)*bm(k)
      cp(k) = cm(k)/den
      dph(k) = (rhsph(i,j,k)-dph(k-1)*bm(k))/den
    ENDDO
    ph(i,j,kend) = dph(kend)
    DO k=kend-1,kst,-1
      ph(i,j,k) = dph(k) - cp(k)*ph(i,j,k+1)
    ENDDO

  END SUBROUTINE solvetridiag

  SUBROUTINE solvetridiag_vel(i,j,var,kst,kend)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j,var,kst,kend
    INTEGER :: k
    REAL(DP) :: den

    cp(kst) = cm(kst)/am(kst)
    dph(kst) = rhsh(i,j,kst,var)/am(kst)
    DO k=kst+1,kend
      den = am(k) - cp(k-1)*bm(k)
      cp(k) = cm(k)/den
      dph(k) = (rhsh(i,j,k,var)-dph(k-1)*bm(k))/den
    ENDDO
    uh(i,j,kend,var) = dph(kend)
    DO k=kend-1,kst,-1
      uh(i,j,k,var) = dph(k) - cp(k)*uh(i,j,k+1,var)
    ENDDO

  END SUBROUTINE solvetridiag_vel

  SUBROUTINE solveuhat

    IMPLICIT NONE
    INTEGER :: i, j, var

    DO j=syf,eyf
     DO i=sxf,exf
       ! u :: tangential
       var = 1; am = amt; bm = bmt; cm = cmt
       call solvetridiag_vel(i,j,var,szf,ezf)      ! szf=1; ezf=nzg
       ! v :: tangential
       var = 2; am = amt; bm = bmt; cm = cmt
       call solvetridiag_vel(i,j,var,szf,ezf)      ! szf=1; ezf=nzg
       ! w :: normal
       var = 3; am = amn; bm = bmn; cm = cmn
       call solvetridiag_vel(i,j,var,szf,ezf-1)    ! szf=1; ezf=nzg
       ! rho :: tangential
       var = 4; am = amr; bm = bmr; cm = cmr
       call solvetridiag_vel(i,j,var,szf,ezf)      ! szf=1; ezf=nzg
     ENDDO
    ENDDO

  END SUBROUTINE solveuhat

  SUBROUTINE rhspres

    IMPLICIT NONE
    INTEGER :: i, j, k

    rhsph = 0.0d0
    do k = szf, ezf           ! can also be 1, nzg
     do j = syf, eyf
      do i = sxf, exf
       if(k==szf) then
         dwdzh(i,j,k) = uh(i,j,k,3)/(zz(k)-z(k))
       elseif(k==ezf) then
         dwdzh(i,j,k) = -uh(i,j,k-1,3)/(z(k)-zz(k-1))
       else
         dwdzh(i,j,k) = (uh(i,j,k,3)-uh(i,j,k-1,3))/(zz(k)-zz(k-1))
       endif
       rhsph(i,j,k) = (iu*(k1(i)*uh(i,j,k,1) + k2(j)*uh(i,j,k,2)) + dwdzh(i,j,k))/dt
      enddo
     enddo
    enddo

  END SUBROUTINE rhspres

  SUBROUTINE solvePoisson

    IMPLICIT NONE
    INTEGER :: i, j
    REAL(DP) :: kmag2

    !iitt = iitt+1000
    DO j = syf, eyf
     DO i = sxf, exf
      kmag2 = k1(i)**2+k2(j)**2
      am(:) = amg(:)-kmag2; bm = bmg; cm = cmg;
      IF(i==1 .and. j==1) then
        am(szf) = am(szf); bm(szf) = 0.0d0; cm(szf) = 0.0d0; rhsph(1,1,szf) = (0.0d0,0.0d0)
        am(ezf) = am(ezf); bm(ezf) = 0.0d0; cm(ezf) = 0.0d0; rhsph(1,1,ezf) = (0.0d0,0.0d0)
      ENDIF
      call solvetridiag(i,j,szf,ezf)
      !write(myrank+100+iitt,'(2(i3.3,1x),3(e19.12,1x))') i, j, maxval(cdabs(ph(i,j,:))), maxval(cdabs(rhsph(i,j,:))), kmag2
     ENDDO
    ENDDO
    !ph(nptsx/2+1,:,:) = 0.0d0   ! odd ball modes
    !write(myrank+100,*) '---'

  END SUBROUTINE solvePoisson

  SUBROUTINE projector

    IMPLICIT NONE
    INTEGER :: i, j, k
    REAL(DP) :: den

    do i = sxf, exf
      uh(i,:,:,1) = uh(i,:,:,1) - dt*iu*k1(i)*ph(i,:,:)
    enddo

    do j = syf, eyf
      uh(:,j,:,2) = uh(:,j,:,2) - dt*iu*k2(j)*ph(:,j,:)
    enddo

    do k = szf, ezf-1 ! or 1, nzg-1
      den = z(k+1)-z(k)
      dpdzh(:,:,k) = (ph(:,:,k+1) - ph(:,:,k))/den
      uh(:,:,k,3) = uh(:,:,k,3) - dt*dpdzh(:,:,k)
    enddo

  END SUBROUTINE projector

  SUBROUTINE correct

    IMPLICIT NONE

    call rhspres
    call solvePoisson
    call projector

  END SUBROUTINE correct

  SUBROUTINE clip

    IMPLICIT NONE
    INTEGER :: i, j, k, var

    do k = sz, ez
     do j = sy, ey
      do i = sx, ex
        u(i,j,k,4) = min(u(i,j,k,4),Rmax)
        u(i,j,k,4) = max(u(i,j,k,4),Rmin)
      enddo
     enddo
    enddo

    var = 4
    call p3dfft_ftran_r2c(u(:,:,:,var),uh(:,:,:,var),fstr)
    uh(:,:,:,var) = itfac*uh(:,:,:,var)

  END SUBROUTINE clip

  SUBROUTINE divergence

    IMPLICIT NONE
    INTEGER :: i, j, k, mpierr

    do k = szf, ezf           ! can also be 1, nzg
     do j = syf, eyf
      do i = sxf, exf
       if(k==szf) then
         dwdzh(i,j,k) = uh(i,j,k,3)/(zz(k)-z(k))
       elseif(k==ezf) then
         dwdzh(i,j,k) = -uh(i,j,k-1,3)/(z(k)-zz(k-1))
       else
         dwdzh(i,j,k) = (uh(i,j,k,3)-uh(i,j,k-1,3))/(zz(k)-zz(k-1))
       endif
       tmph(i,j,k) = iu*(k1(i)*uh(i,j,k,1) + k2(j)*uh(i,j,k,2)) + dwdzh(i,j,k)
      enddo
     enddo
    enddo

    call p3dfft_btran_c2r(tmph(:,:,:),divu(:,:,:),fstr)

    IF(sz==1) divu(:,:,sz) = 0.0d0
    IF(ez==nzg) divu(:,:,ez) = 0.0d0

    ! compute utau
    if(sz==1) then
      utaul1 = maxval((u(:,:,2,1)-u(:,:,1,1)))/(z(2)-z(1))
      utaul1 = dsqrt(utaul1**2.0d0 + (maxval((u(:,:,2,2)-u(:,:,1,2)))/(z(2)-z(1)))**2.0d0)
    else
      utaul1 = 0.0d0
    endif
    if(ez==nzg) then
      utaul2 = maxval((u(:,:,nzg,1)-u(:,:,nzg-1,1)))/(z(nzg)-z(nzg-1))
      utaul2 = dsqrt(utaul2**2.0d0 + (maxval((u(:,:,nzg,2)-u(:,:,nzg-1,2)))/(z(nzg)-z(nzg-1)))**2.0d0)
    else
      utaul2 = 0.0d0
    endif

    senddiag(1) = maxval(divu); senddiag(2) = minval(divu)
    senddiag(3) = maxval(u(:,:,:,1)); senddiag(4) = -minval(u(:,:,:,1))
    senddiag(5) = maxval(u(:,:,:,2)); senddiag(6) = -minval(u(:,:,:,2))
    senddiag(7) = maxval(u(:,:,:,3)); senddiag(8) = -minval(u(:,:,:,3))
    senddiag(9) = maxval(u(:,:,:,4)); senddiag(10) = -minval(u(:,:,:,4))
    senddiag(11) = utaul1; senddiag(12) = utaul2

    call MPI_REDUCE(senddiag,recvdiag,12,MPIRP,MPI_MAX,outid,MPICOMM2D,mpierr)

    maxdivu = recvdiag(1); mindivu = -recvdiag(2)
    maxu1 = recvdiag(3); minu1 = -recvdiag(4)
    maxu2 = recvdiag(5); minu2 = -recvdiag(6)
    maxu3 = recvdiag(7); minu3 = -recvdiag(8)
    maxu4 = recvdiag(9); minu4 = -recvdiag(10)
    utau1 = recvdiag(11); utau2 = recvdiag(12)

    maxcfl = dsqrt(maxu1**2+maxu2**2+maxu3**2)*dt/mindz

  END SUBROUTINE divergence

  SUBROUTINE rhsmom_ABCN

    IMPLICIT NONE

    !! already done at the end ofprevious time step
    !call exchangez
    !call vinterp
    !call vinterph
    !call xyderivs
    !call zderivs

    call convect
    call viscous

    call Corriolis
    call buoyancy

    call globdynamic_gradmodel
    call sgstress

    rhsh = uh + dt*(1.5d0*(-convh+viscth+buoyh-sgstermh-fCorriolh) -0.5d0*rhs1h + viscnh)
    !rhsh(:,:,:,1) = rhsh(:,:,:,1) + dt_pgradh   ! only for channel flow

    ! stored for next time ABCN step
    rhs1h = -convh + viscth + buoyh - sgstermh - fCorriolh


  END SUBROUTINE rhsmom_ABCN

  SUBROUTINE init_gather_write_line

    IMPLICIT NONE
    INTEGER :: k

    OPEN(15,FILE='linedat.dat',STATUS='replace')
    WRITE(15,*) 'VARIABLES="Z","k","omega","nut"'
    CLOSE(15)

  END SUBROUTINE init_gather_write_line

  SUBROUTINE gather_write_line

    IMPLICIT NONE
    INTEGER :: i, j, k, mpierr

    !i = nxg/3; j = nyg/3
    call MPI_GATHER(coef(sz:ez),nz,MPIRP,kg,nz,MPIRP,outid,MPI_COMM_WORLD,mpierr)
    call MPI_GATHER(coefrho(sz:ez),nz,MPIRP,omg,nz,MPIRP,outid,MPI_COMM_WORLD,mpierr)
    call MPI_GATHER(coef(sz:ez),nz,MPIRP,nutg,nz,MPIRP,outid,MPI_COMM_WORLD,mpierr)

    IF(MYRANK==outid) then
    OPEN(15,FILE='linedat.dat',STATUS='old',ACTION='write',POSITION='append')
    WRITE(15,'(a,e12.5,a,i4)') 'ZONE T="t=', t, '", F=POINT, I=', nzg
    DO k = 1, nzg
      WRITE(15,'(4(e12.5,1x))') z(k), kg(k), omg(k), nutg(k)
    ENDDO
    CLOSE(15)
    ENDIF

  END SUBROUTINE gather_write_line


  SUBROUTINE Euler

    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER :: rr, mpierr

    !write(myrank+100,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+100,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+100,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+100,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+100,*) '----'
    call rhsmom_Euler
    !! debug
    !call p3dfft_btran_c2r(rhsh(:,:,:,1),dpdx(:,:,:),bstr)
    !call p3dfft_btran_c2r(rhsh(:,:,:,2),dpdy(:,:,:),bstr)
    !call p3dfft_btran_c2r(rhsh(:,:,:,3),dpdz(:,:,:),bstr)
    !call p3dfft_btran_c2r(rhsh(:,:,:,4),rhsp(:,:,:),bstr)
    !write(myrank+100,*) 'rhs1', maxval(dpdx), minval(dpdx)
    !write(myrank+100,*) 'rhs2', maxval(dpdy), minval(dpdy)
    !write(myrank+100,*) 'rhs3', maxval(dpdz), minval(dpdz)
    !write(myrank+100,*) 'rhs4', maxval(rhsp), minval(rhsp)
    !write(myrank+100,*) '----'
    !write(myrank+100,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+100,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+100,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+100,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+100,*) '----'
    
    call rhsbc
    call solveuhat

    !! debug
    !WRITE(myrank+100,*) 'u1-1', maxval(cdabs(uh(:,:,szf,1))), minval(cdabs(uh(:,:,szf,1)))
    !WRITE(myrank+100,*) 'u1-2', maxval(cdabs(uh(:,:,ezf,1))), minval(cdabs(uh(:,:,ezf,1)))
    !WRITE(myrank+100,*) 'u2-1', maxval(cdabs(uh(:,:,szf,2))), minval(cdabs(uh(:,:,szf,2)))
    !WRITE(myrank+100,*) 'u2-2', maxval(cdabs(uh(:,:,ezf,2))), minval(cdabs(uh(:,:,ezf,2)))
    !!call decode
    !!write(myrank+100,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !!write(myrank+100,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !!write(myrank+100,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !!write(myrank+100,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+100,*) 'p ', maxval(p(:,:,:)), minval(p(:,:,:))
    !write(myrank+100,*) '----'
    !do j = syf, eyf
    !  rhsph(:,j,:) = iu*k2(j)*uh(:,j,:,2)
    !enddo
    !call p3dfft_btran_c2r(rhsph(:,:,:),dpdy(:,:,:),bstr)
    !do i = sxf, exf
    !  rhsph(i,:,:) = iu*k1(i)*uh(i,:,:,1)
    !enddo
    !call p3dfft_btran_c2r(rhsph(:,:,:),dpdx(:,:,:),bstr)
    !write(myrank+100,*) 'dudx', maxval(dpdx(:,:,:)), minval(dpdx(:,:,:))
    !write(myrank+100,*) 'dvdy', maxval(dpdy(:,:,:)), minval(dpdy(:,:,:))

    call correct
    !WRITE(myrank+100,*) 'u1-1', maxval(cdabs(uh(:,:,szf,1))), minval(cdabs(uh(:,:,szf,1)))
    !WRITE(myrank+100,*) 'u1-2', maxval(cdabs(uh(:,:,ezf,1))), minval(cdabs(uh(:,:,ezf,1)))
    !WRITE(myrank+100,*) 'u2-1', maxval(cdabs(uh(:,:,szf,2))), minval(cdabs(uh(:,:,szf,2)))
    !WRITE(myrank+100,*) 'u2-2', maxval(cdabs(uh(:,:,ezf,2))), minval(cdabs(uh(:,:,ezf,2)))
    !write(myrank+100,*) '----'
    call dealiasuh
    !WRITE(myrank+100,*) 'u1-1', maxval(cdabs(uh(:,:,szf,1))), minval(cdabs(uh(:,:,szf,1)))
    !WRITE(myrank+100,*) 'u1-2', maxval(cdabs(uh(:,:,ezf,1))), minval(cdabs(uh(:,:,ezf,1)))
    !WRITE(myrank+100,*) 'u2-1', maxval(cdabs(uh(:,:,szf,2))), minval(cdabs(uh(:,:,szf,2)))
    !WRITE(myrank+100,*) 'u2-2', maxval(cdabs(uh(:,:,ezf,2))), minval(cdabs(uh(:,:,ezf,2)))
    call bch

    call decode
    !! debug
    !do rr = 0, nprocz*nprocy-1
    !  IF(rr==myrank) then
    !    OPEN(110,FILE='fort110',STATUS='old',ACTION='write',POSITION='append')
    !    write(110,'(i3.3,2(1x,e19.12))') rr,maxval(p(:,:,:)), minval(p(:,:,:)) 
    !    if(rr==nprocz*nprocy-1) write(110,*) '---'
    !    CLOSE(110)
    !  ENDIF
    !  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    !enddo

    !call correct
    !call bc
    !! debug
    !call clip
    !call p3dfft_btran_c2r(dwdzh(:,:,:),dudz(:,:,:,3),bstr)
    !!write(myrank+100,*) 'dwdz', maxval(dudz(:,:,:,3)), minval(dudz(:,:,:,3))
    !call p3dfft_btran_c2r(rhsph(:,:,:),rhsp(:,:,:),bstr)
    !do i = sxf, exf
    !  rhsph(i,:,:) = iu*k1(i)*ph(i,:,:)
    !enddo
    !call p3dfft_btran_c2r(rhsph(:,:,:),dpdx(:,:,:),bstr)
    !do j = syf, eyf
    !  rhsph(:,j,:) = iu*k2(j)*ph(:,j,:)
    !enddo
    !call p3dfft_btran_c2r(rhsph(:,:,:),dpdy(:,:,:),bstr)
    !call p3dfft_btran_c2r(dpdzh(:,:,:),dpdz(:,:,:),str)
    !write(myrank+100,*) '----'
    !IF(sz==1) WRITE(myrank+100,*) 'u1-1', maxval(dabs(u(:,:,sz,1))), minval(dabs(u(:,:,sz,1)))
    !IF(ez==nzg) WRITE(myrank+100,*) 'u1-2', maxval(dabs(u(:,:,ez,1))), minval(dabs(u(:,:,ez,1)))
    !IF(sz==1) WRITE(myrank+100,*) 'u2-1', maxval(dabs(u(:,:,sz,2))), minval(dabs(u(:,:,sz,2)))
    !IF(ez==nzg) WRITE(myrank+100,*) 'u2-2', maxval(dabs(u(:,:,ez,2))), minval(dabs(u(:,:,ez,2)))
    !write(myrank+100,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+100,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+100,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+100,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+100,*) 'p ', maxval(p(:,:,:)), minval(p(:,:,:))
    !write(myrank+100,*) 'rhsp ', maxval(rhsp(:,:,:)), minval(rhsp(:,:,:))
    !write(myrank+100,*) 'dpdx ', maxval(dpdx(:,:,:)), minval(dpdx(:,:,:))
    !write(myrank+100,*) 'dpdy ', maxval(dpdy(:,:,:)), minval(dpdy(:,:,:))
    !write(myrank+100,*) 'dpdz ', maxval(dpdz(:,:,:)), minval(dpdz(:,:,:))
    ! end debug
    call clip
    call divergence

    IF(MYRANK==outid) then
      write(*,35) 'Eul: ',t, max(maxdivu, mindivu), maxu1, maxu2, maxcfl, utau1
    ENDIF

    ! prepare for stats or next time step
    iwpov = 1; call exchangez(iwpov)
    call vinterp
    call vinterph
    call bc
    call xyderivs
    call zderivs

35 format(a,7(1x,e12.5),1x,i5,1x,i9)

  END SUBROUTINE Euler

  SUBROUTINE ABCN

    IMPLICIT NONE
    INTEGER :: i, j

    !write(myrank+200,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+200,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+200,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+200,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+200,*) '----'
    call rhsmom_ABCN
    !! debug
    !call p3dfft_btran_c2r(rhsh(:,:,:,1),dpdx(:,:,:),bstr)
    !call p3dfft_btran_c2r(rhsh(:,:,:,2),dpdy(:,:,:),bstr)
    !call p3dfft_btran_c2r(rhsh(:,:,:,3),dpdz(:,:,:),bstr)
    !call p3dfft_btran_c2r(rhsh(:,:,:,4),rhsp(:,:,:),bstr)
    !write(myrank+200,*) 'rhs1', maxval(dpdx), minval(dpdx)
    !write(myrank+200,*) 'rhs2', maxval(dpdy), minval(dpdy)
    !write(myrank+200,*) 'rhs3', maxval(dpdz), minval(dpdz)
    !write(myrank+200,*) 'rhs4', maxval(rhsp), minval(rhsp)
    !write(myrank+200,*) '----'
    !write(myrank+200,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+200,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+200,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+200,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+200,*) '----'
    
    call rhsbc
    call solveuhat

    !! debug
    !call decode
    !write(myrank+200,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+200,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+200,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+200,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+200,*) 'p ', maxval(p(:,:,:)), minval(p(:,:,:))
    !write(myrank+200,*) '----'

    call correct
    call dealiasuh
    call bch

    call decode
    !call bc
    !! debug
    !call bc
    !call p3dfft_btran_c2r(dwdzh(:,:,:),dudz(:,:,:,3),bstr)
    !!write(myrank+200,*) 'dwdz', maxval(dudz(:,:,:,3)), minval(dudz(:,:,:,3))
    !call p3dfft_btran_c2r(rhsph(:,:,:),rhsp(:,:,:),bstr)
    !do i = sxf, exf
    !  rhsph(i,:,:) = iu*k1(i)*ph(i,:,:)
    !enddo
    !call p3dfft_btran_c2r(rhsph(:,:,:),dpdx(:,:,:),bstr)
    !do j = syf, eyf
    !  rhsph(:,j,:) = iu*k2(j)*ph(:,j,:)
    !enddo
    !call p3dfft_btran_c2r(rhsph(:,:,:),dpdy(:,:,:),bstr)
    !call p3dfft_btran_c2r(dpdzh(:,:,:),dpdz(:,:,:),bstr)
    !write(myrank+200,*) '----'
    !write(myrank+200,*) 'u1', maxval(u(:,:,:,1)), minval(u(:,:,:,1))
    !write(myrank+200,*) 'u2', maxval(u(:,:,:,2)), minval(u(:,:,:,2))
    !write(myrank+200,*) 'u3', maxval(u(:,:,:,3)), minval(u(:,:,:,3))
    !write(myrank+200,*) 'u4', maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    !write(myrank+200,*) 'p ', maxval(p(:,:,:)), minval(p(:,:,:))
    call clip
    call divergence

    IF(MYRANK==outid) then
      write(*,35) 'ABCN:',t, max(maxdivu, mindivu), maxu1, maxu2, maxcfl, utau1
    ENDIF

35 format(a,7(1x,e12.5),1x,i5,1x,i9)

    ! prepare for stats or next time step
    iwpov = 1; call exchangez(iwpov)
    call vinterp
    call vinterph
    call bc
    call xyderivs
    call zderivs

  END SUBROUTINE ABCN

  SUBROUTINE massflowrate

    IMPLICIT NONE
    INTEGER :: i,j,k,mpierr
    REAL(DP) :: mfr_l

    mfr_l = 0.0d0
    do k = sz, ez
     do j = sy, ey
      do i = sx, ex
        mfr_l = mfr_l + u(i,j,k,1)*dy*dz(k)
      enddo
     enddo
    enddo
    
    call MPI_REDUCE(mfr_l,mfr,1,MPIRP,MPI_SUM,outid,MPICOMM2D,mpierr)
    mfr = mfr/dble(nxg*ly*lz)
    IF(MYRANK==outid) WRITE(222,*) t, mfr

  END SUBROUTINE massflowrate

  SUBROUTINE setPresGrad

    IMPLICIT NONE

    !pgr = utau**2

  END SUBROUTINE setPresGrad

  SUBROUTINE stats

    IMPLICIT NONE
    INTEGER :: var,j,k
    REAL(DP) :: nxny

    !! already done at end of previous time step
    !call decode
    !call exchangez
    !call vinterp
    !call vinterph
    !call bc
    !call xyderivs
    !call zderivs
    !call subzderivs

    statindex = statindex + 1
    nxny = dble(nxg*nyg)

    do k = sz, ez

      uavg(k,1) = uavg(k,1) + sum(u(:,:,k,1))/nxny
      uavg(k,2) = uavg(k,2) + sum(u(:,:,k,2))/nxny
      uavg(k,3) = uavg(k,3) + sum(wp(:,:,k))/nxny
      uavg(k,4) = uavg(k,4) + sum(u(:,:,k,4))/nxny

      uuavg(k,1,1) = uuavg(k,1,1)+sum(u(:,:,k,1)*u(:,:,k,1))/nxny
      uuavg(k,1,2) = uuavg(k,1,2)+sum(u(:,:,k,1)*u(:,:,k,2))/nxny
      uuavg(k,1,3) = uuavg(k,1,3)+sum(u(:,:,k,1)*wp(:,:,k))/nxny
      uuavg(k,2,2) = uuavg(k,2,2)+sum(u(:,:,k,2)*u(:,:,k,2))/nxny
      uuavg(k,2,3) = uuavg(k,2,3)+sum(u(:,:,k,2)*wp(:,:,k))/nxny
      uuavg(k,3,3) = uuavg(k,3,3)+sum(wp(:,:,k)*wp(:,:,k))/nxny
      uuavg(k,2,1) = uuavg(k,2,1)+sum(u(:,:,k,1)*u(:,:,k,4))/nxny
      uuavg(k,3,1) = uuavg(k,3,1)+sum(u(:,:,k,2)*u(:,:,k,4))/nxny
      uuavg(k,3,2) = uuavg(k,3,2)+sum(wp(:,:,k)*u(:,:,k,4))/nxny
      ttavg(k) = ttavg(k)+sum(u(:,:,k,4)*u(:,:,k,4))/nxny

      ssavg_d(k) = ssavg_d(k) + sum(2.0d0*dudx(:,:,k,1)**2+dudx(:,:,k,2)**2+dwpdx(:,:,k)**2   &
                              + dudy(:,:,k,1)**2+2.0d0*dudy(:,:,k,2)**2+dwpdy(:,:,k)**2       &
                              + dudz(:,:,k,1)**2+dudz(:,:,k,2)**2+2.0d0*dwpdz(:,:,k)**2       &
                              + 2.0d0*dudy(:,:,k,1)*dudx(:,:,k,2) &
                              + 2.0d0*dudz(:,:,k,1)*dwpdx(:,:,k) &
                              + 2.0d0*dudz(:,:,k,2)*dwpdy(:,:,k))/nxny

      ssavg_rhod(k) = ssavg_rhod(k) + sum(dudx(:,:,k,4)**2+dudy(:,:,k,4)**2+dudz(:,:,k,4)**2)/nxny

      savg(k,1,1) = savg(k,1,1) + sum(dudx(:,:,k,1))/nxny
      savg(k,1,2) = savg(k,1,2) + sum(dudy(:,:,k,1))/nxny
      savg(k,1,3) = savg(k,1,3) + sum(dudz(:,:,k,1))/nxny
      savg(k,2,1) = savg(k,2,1) + sum(dudx(:,:,k,2))/nxny
      savg(k,2,2) = savg(k,2,2) + sum(dudy(:,:,k,2))/nxny
      savg(k,2,3) = savg(k,2,3) + sum(dudz(:,:,k,2))/nxny
      savg(k,3,1) = savg(k,3,1) + sum(dwpdx(:,:,k))/nxny
      savg(k,3,2) = savg(k,3,2) + sum(dwpdy(:,:,k))/nxny
      savg(k,3,3) = savg(k,3,3) + sum(dwpdz(:,:,k))/nxny

      tavg(k,1) = tavg(k,1) + sum(dudx(:,:,k,4))/nxny
      tavg(k,2) = tavg(k,2) + sum(dudy(:,:,k,4))/nxny
      tavg(k,3) = tavg(k,3) + sum(dudz(:,:,k,4))/nxny

      epssgs(k) = epssgs(k) - sum(tausgs(:,:,k,1,1)*dudx(:,:,k,1) + &
                                  tausgs(:,:,k,2,2)*dudy(:,:,k,2) + &
                                  tausgs(:,:,k,3,3)*dwpdz(:,:,k)  + &
                                  tausgs(:,:,k,1,2)*(dudy(:,:,k,1)+dudx(:,:,k,2)) + &
                                  tausgs(:,:,k,1,3)*(dudz(:,:,k,1)+dwpdx(:,:,k)) + &
                                  tausgs(:,:,k,2,3)*(dudz(:,:,k,2)+dwpdy(:,:,k)))/nxny

      epssgs_rho(k) = epssgs_rho(k) - sum(tausgs(:,:,k,4,1)*dudx(:,:,k,4) + &
                                          tausgs(:,:,k,4,2)*dudy(:,:,k,4) + &
                                          tausgs(:,:,k,4,3)*dudz(:,:,k,4))/nxny

      tauavg(k,1,1) = tauavg(k,1,1) + sum(tausgs(:,:,k,1,1))/nxny
      tauavg(k,1,2) = tauavg(k,1,2) + sum(tausgs(:,:,k,1,2))/nxny
      tauavg(k,1,3) = tauavg(k,1,3) + sum(tausgs(:,:,k,1,3))/nxny
      tauavg(k,2,2) = tauavg(k,2,2) + sum(tausgs(:,:,k,2,2))/nxny
      tauavg(k,2,3) = tauavg(k,2,3) + sum(tausgs(:,:,k,2,3))/nxny
      tauavg(k,3,3) = tauavg(k,3,3) + sum(tausgs(:,:,k,3,3))/nxny
      tauavg(k,2,1) = tauavg(k,2,1) + sum(tausgs(:,:,k,4,1))/nxny
      tauavg(k,3,1) = tauavg(k,3,1) + sum(tausgs(:,:,k,4,2))/nxny
      tauavg(k,3,2) = tauavg(k,3,2) + sum(tausgs(:,:,k,4,3))/nxny

    enddo

    do var = 1, nvar
     if(var==3) then
      do k=sz,ez
       do j=sy,ey
        tmp1d = wp(:,j,k)
        call dfftw_execute_dft_r2c(plan_f,tmp1d,tmp1dh)
        spec(:,k,var) = spec(:,k,var) + dconjg(tmp1dh)*tmp1dh/nxny
       enddo
      enddo
     else
      do k=sz,ez
       do j=sy,ey
        tmp1d = u(:,j,k,var)
        call dfftw_execute_dft_r2c(plan_f,tmp1d,tmp1dh)
        spec(:,k,var) = spec(:,k,var) + dconjg(tmp1dh)*tmp1dh/nxny
       enddo
      enddo
     endif
    enddo

  END SUBROUTINE stats

  SUBROUTINE writestats

    IMPLICIT NONE
    INTEGER :: i, k
    REAL(DP) :: dissp
    CHARACTER(LEN=19) :: fname

    itstat = itstat + 1
    WRITE(fname,1) 'stats_',MYRANK,'_',itstat,'.dat'
    OPEN(10,FILE=fname,STATUS='replace')
    WRITE(10,2) sx, ex, sy, ey, sz, ez, t, statindex
    DO k = sz, ez
      !dissp = ssavg_d(k) - (2.0d0*savg(k,1,1)**2 + savg(k,1,2)**2 + savg(k,1,3)**2 &
      !                      +  savg(k,2,1)**2 + 2.0d0*savg(k,2,2)**2 + savg(k,2,3)**2 &
      !                      +  savg(k,3,1)**2 + savg(k,3,2)**2 + 2.0d0*savg(k,3,3)**2 &
      !                      + 2.0d0*savg(k,1,2)*savg(k,2,1) &
      !                      + 2.0d0*savg(k,1,3)*savg(k,3,1) &
      !                      + 2.0d0*savg(k,2,3)*savg(k,3,2) )

      WRITE(10,3) uavg(k,1:4), uuavg(k,1,1:3), uuavg(k,2,2:3), uuavg(k,3,3), uuavg(k,2,1), uuavg(k,3,1:2), ttavg(k), ssavg_d(k), savg(k,1,1:3), savg(k,2,1:3), savg(k,3,1:3), tavg(k,:), tauavg(k,1,1:3), tauavg(k,2,2:3), tauavg(k,3,3), tauavg(k,2,1), tauavg(k,3,1:2), epssgs(k), epssgs_rho(k), ssavg_rhod(k)
    ENDDO
    CLOSE(10)

    WRITE(fname,1) 'spech_',MYRANK,'_',itstat,'.dat'
    OPEN(10,FILE=fname,STATUS='replace')
    WRITE(10,2) sx, ex, sy, ey, sz, ez, t, statindex
    DO k = sz, ez
      WRITE(10,4) spec(:,k,1)
      WRITE(10,4) spec(:,k,2)
      WRITE(10,4) spec(:,k,3)
      WRITE(10,4) spec(:,k,4)
    ENDDO
    CLOSE(10)

1 FORMAT(a,i4.4,a,i4.4,a)
2 FORMAT(6(i4.4,1x),e19.12,1x,i7.7)
3 FORMAT(39(e19.12,1x))
4 FORMAT(100(e19.12,1x))

  END SUBROUTINE writestats

  SUBROUTINE output

    IMPLICIT NONE

    call decode
    iwpov = 1; call exchangez(iwpov)
    call vinterp
    call vinterph
    !call bc
    call xyderivs
    call zderivs
    call process
    call writedata

  END SUBROUTINE output

END MODULE SUBROUTINES

PROGRAM cfp3dfft

  USE subroutines
  IMPLICIT NONE
  INTEGER :: numTimeSnaps, numStatSnaps, ierr 
  REAL(DP) :: endtime, starttime, tmfr, tg
  CHARACTER(LEN=16) :: fname

  OPEN(10,FILE='input.in',STATUS='old',ACTION='read')
  READ(10,*) Re, Pr, Ri, nvar
  READ(10,*) dt, tmax, nxg, nyg, nzg
  READ(10,*) nprocy, nprocz, outid
  READ(10,*) istretch, aa, bb
  READ(10,*) tstatstart, irestart, tsnap, numTimeSnaps, numStatSnaps
  READ(10,*) nxgrst, nygrst, nzgrst
  CLOSE(10)

  call startmpi
  call startp3dfft

  call memvars
  call arrayinit
  call grid
  call metrics

  call init_tridiag
  call init_tridiag_vel

  if(irestart==0) then
    call init
    t = 0.0d0
    tmfr = 0
  else
    call restartinit
    tmfr = 0
  endif
  call outputIC

  ! Frequency of application of dynamic procedure
  DYNFREQ = 1
  idyn = DYNFREQ-1 ! ensure global coefficient is computed at first time step
  ! for constant coefficient model
  DYNFREQ = nint((tmax-t)/dt) + 100000       ! an integer larger than the total number of time steps
  idyn = -1                                    ! ensure global coefficient is never computed
  globcoef = 1.0d0; globcoefrho = 1.0d0        ! directly set coefficient values fort first and all later time steps

  starttime = MPI_WTime()

  if(irestart==1) then
    call ABCN
  else
    call Euler
  endif

  t = t+dt
  tsc = tsc+dt

  !! debug - local dynamic procedure
  !call init_gather_write_line
  !call gather_write_line
  !tg = 0.0d0

  do while(t<tmax)
    call setPresGrad
    call ABCN
    t = t+dt
    tsc = tsc+dt

    !!debug - local dynamic procedure
    !tg = tg+dt
    !if(tg.gt.0.02d0*tsnap) then
    !  tg=0.0d0
    !  call gather_write_line
    !endif

    tmfr = tmfr+dt
    if(t>=tstatstart) call stats
    if(tsc.ge.tsnap) then
      call output
      tsc=0.0
      if(t>=tstatstart) call writestats
      WRITE(fname,1) 'restart_',MYRANK,'.out'
      open(66,file=fname,status='unknown')
      WRITE(66,*) t,u,p,rhs1h
      close(66)
      !print *, maxval(u(:,:,:,1)), minval(u(:,:,:,1))
      !print *, maxval(u(:,:,:,2)), minval(u(:,:,:,2))
      !print *, maxval(u(:,:,:,3)), minval(u(:,:,:,3))
      !print *, maxval(u(:,:,:,4)), minval(u(:,:,:,4))
    endif
    if(tmfr.ge.1.0d0) then
      tmfr = 0.0d0
      call massflowrate
    endif
    open(unit=20,file='exit',status='old',iostat=ierr)
    if(ierr==0) exit
  end do

  endtime = MPI_WTime()
  IF(MYRANK==outid) print *, 'Computation Time =', endtime-starttime

  call output
  if(t>=tstatstart) call writestats

  call memvarsdeal
  call stopp3dfft
  call stopmpi

1 FORMAT(a,i4.4,a)

END PROGRAM cfp3dfft
