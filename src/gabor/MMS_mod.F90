module MMS_mod
    use kind_parameters, only: rkind
    use constants, only: pi
    implicit none
    real(rkind) :: tMMS = 0.d0

    contains
      pure subroutine get_MMS_source(f,g)
          real(rkind), dimension(3), intent(inout) :: f, g
          
          f(1) = cos(2.d0*pi*tMMS)*(0.4d0*cos(2.d0*pi*tMMS)**2.d0 + 0.121d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.025d0*sin(4.d0*pi*tMMS)**2.d0) - &
                 ((8.d0*cos(2.d0*pi*tMMS)**2.d0)/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.25d0*sin(4.d0*pi*tMMS)**2.d0) - &
                 1.d0)*(1.1d0*cos(2.d0*pi*tMMS) + 1.3d0*cos(4.d0*pi*tMMS) + 1.2d0*sin(2.d0*pi*tMMS)) - &
                 2.d0*pi*sin(2.d0*pi*tMMS) - &
                 (4.4d0*cos(2.d0*pi*tMMS)*sin(2.d0*pi*tMMS)*(3.1d0*cos(2.d0*pi*tMMS) + 3.3d0*cos(4.d0*pi*tMMS) + &
                 3.2d0*sin(2.d0*pi*tMMS)))/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.25d0*sin(4.d0*pi*tMMS)**2.d0) - &
                 (2.d0*cos(2.d0*pi*tMMS)*sin(4.d0*pi*tMMS)*(2.1d0*cos(2.d0*pi*tMMS) + 2.3d0*cos(4.d0*pi*tMMS) + &
                 2.2d0*sin(2.d0*pi*tMMS)))/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.25d0*sin(4.d0*pi*tMMS)**2.d0)
          f(2) = sin(2.d0*pi*tMMS)*(0.4d0*cos(2.d0*pi*tMMS)**2.d0 + 0.121d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.025d0*sin(4.d0*pi*tMMS)**2.d0) - ((0.5d0*sin(4.d0*pi*tMMS)**2.d0)/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + &
                 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + 0.25d0*sin(4.d0*pi*tMMS)**2.d0) - 1)*(2.1d0*cos(2.d0*pi*tMMS)&
                 + 2.3d0*cos(4.d0*pi*tMMS) + 2.2d0*sin(2.d0*pi*tMMS)) + 2.d0*pi*cos(2.d0*pi*tMMS) - &
                 (1.1d0*sin(2.d0*pi*tMMS)*sin(4.d0*pi*tMMS)*(3.1d0*cos(2.d0*pi*tMMS) + &
                 3.3d0*cos(4.d0*pi*tMMS) + 3.2d0*sin(2.d0*pi*tMMS)))/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + &
                 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + 0.25d0*sin(4.d0*pi*tMMS)**2.d0) - &
                 (2.d0*cos(2.d0*pi*tMMS)*sin(4.d0*pi*tMMS)*(1.1d0*cos(2.d0*pi*tMMS) + 1.3d0*cos(4.d0*pi*tMMS) + &
                 1.2d0*sin(2.d0*pi*tMMS)))/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.25d0*sin(4.d0*pi*tMMS)**2.d0)
          f(3) = cos(4.d0*pi*tMMS)*(0.4d0*cos(2.d0*pi*tMMS)**2.d0 + 0.121d0*sin(2.d0*pi*tMMS)**2.d0 + &
                 0.025d0*sin(4.d0*pi*tMMS)**2.d0) - 4.d0*pi*sin(4.d0*pi*tMMS) - &
                 ((2.42d0*sin(2.d0*pi*tMMS)**2.d0)/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + &
                 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + 0.25d0*sin(4.d0*pi*tMMS)**2.d0) - &
                 1)*(3.1d0*cos(2.d0*pi*tMMS) + 3.3d0*cos(4.d0*pi*tMMS) + &
                 3.2d0*sin(2.d0*pi*tMMS)) - &
                 (1.1d0*sin(2.d0*pi*tMMS)*sin(4.d0*pi*tMMS)*(2.1d0*cos(2.d0*pi*tMMS) + &
                 2.3d0*cos(4.d0*pi*tMMS) + 2.2d0*sin(2.d0*pi*tMMS)))/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + &
                 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + 0.25d0*sin(4.d0*pi*tMMS)**2.d0) - &
                 (4.4d0*cos(2.d0*pi*tMMS)*sin(2.d0*pi*tMMS)*(1.1d0*cos(2.d0*pi*tMMS) + &
                 1.3d0*cos(4.d0*pi*tMMS) + 1.2d0*sin(2.d0*pi*tMMS)))/(4.d0*cos(2.d0*pi*tMMS)**2.d0 + &
                 1.21d0*sin(2.d0*pi*tMMS)**2.d0 + 0.25d0*sin(4.d0*pi*tMMS)**2.d0)
          g(1) = 2.2d0*cos(2.d0*pi*tMMS) + 3.41d0*sin(2.d0*pi*tMMS) + 1.05d0*sin(4.d0*pi*tMMS) &
                 - 4.d0*pi*sin(2.d0*pi*tMMS)
          g(2) = 2.4d0*cos(2.d0*pi*tMMS) + 3.52d0*sin(2.d0*pi*tMMS) + 1.1d0*sin(4.d0*pi*tMMS) &
                 + 2.d0*pi*cos(4.d0*pi*tMMS)
          g(3) = 2.6d0*cos(2.d0*pi*tMMS) + 3.63d0*sin(2.d0*pi*tMMS) + 1.15d0*sin(4.d0*pi*tMMS) &
                 + 2.2d0*pi*cos(2.d0*pi*tMMS)

      end subroutine
end module MMS_mod
