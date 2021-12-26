    subroutine equilibratePressureTemperature(this,mixRho,mixE,mixP,mixT,isub)
        use reductions, only : P_MINVAL,P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(4)             :: iparams
        integer, intent(in)               :: isub

        integer :: i,j,k,imat,icount
        real(rkind) :: maxp, peqb,rhomin,pmin,pmax
        
        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp

           if (not in interface, use old method)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
           !! OLD pT eqb code
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

            ! set fparams
            fparams(1) = mixRho(i,j,k)*ehmix(i,j,k)

            ! set iparams
            iparams(1) = 2
            iparams(2) = i; iparams(3) = j; iparams(4) = k;

            maxp = zero; peqb = zero
            do imat=1,this%ns
              ! set initial guess
              peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k)
            end do

            ! solve non-linear equation
            call this%rootfind_nr_1d(peqb,fparams,iparams) !old

            !Assign equilibrium pressure to mixP output variable    
            mixP(i,j,k) = peqb !*maxp

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

           
           else (in interface, use new method)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
           !! NEW incompressible phase EOS code
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
          
           ! 1) calculate the species densities and volume fractions
           ! (conservation of mass)
           ! reference density is stored in this%material(imat)%elastic%rho0
           ! this%material(2)%rho = this%material(2)%elastic%rho0

           ! 2) iteratively solve for equilibrium temperature
           ! 2a) set initial guess for temperaures and energies based on
           ! previous timestep
           ! 2b) calculate partial derivatives (analytical based Stiff EOS)
           ! 2c) calculate update of T_eqb until it converges 
                 ! (also update species energies)       
                 ! (may be necessary to use a limiter)

           ! 3) calculate species pressures analytical based on
           ! densities and temperature
           ! 
           ! 4) calculate mixture pressure based on the volume fractions
           ! and species pressures
            maxp = zero; peqb = zero
            do imat=1,this%ns
              ! set initial guess
              peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k)
            end do
        
           endif     

          enddo
         enddo
        enddo

        IF in interface (OLD METHOD) 

        !Calculate mixture temperature based on mixture pressure                
        mixT = zero
        do i = 1, this%ns
          mixT = mixT + this%material(i)%Ys*this%material(i)%hydro%Cv * &
                (mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)/(mixP + this%material(i)%hydro%PInf)
        enddo
        mixT = ehmix/mixT

        do i = 1, this%ns
            this%material(i)%T = mixT !assign species temperatures to mixture T
            this%material(i)%p = mixP !assign species pressures    to mixture p
            this%material(i)%VF = mixRho*this%material(i)%Ys*(this%material(i)%hydro%gam-one)* &
                                         this%material(i)%hydro%Cv*mixT/(mixP + this%material(i)%hydro%PInf)

            this%material(i)%eh = this%material(i)%hydro%Cv*mixT*(mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf) / &
                                                                 (mixP + this%material(i)%hydro%PInf)
        end do

        ELSE
        
        !TODO: Make sure that all output arguments are assigned

        ENDIF

    end subroutine


!!! TODOs
! 1) create arrays for species densities in the structure this%material(:)%rho in SolidMod.F90
