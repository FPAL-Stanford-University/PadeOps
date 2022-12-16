program test_seed
    use seedGen
    use random,    only: uniform_random
    use kind_parameters, only: rkind 
    use constants
    implicit none
    integer :: seed1, seed3, seed2
    real(rkind), dimension(5) :: theta

    seed1 = get_seed_from_location(-0.17500000000000002d0,-0.22500000000000001d0,-0.22500000000000001d0, 1.D18)
    seed2 = get_seed_from_location(-0.1750005000000002d0, -0.22500000000000001d0,-0.22500000000000001d0, 1.D18)
    seed3 = get_seed_from_location(-0.17500000000000002d0, 0.22500000000000001d0,-0.22500000000000001d0, 1.D18)

    print*, seed1, seed2, seed3
    call uniform_random(theta,0.d0,2.d0*pi,seed1)
    print*, "From seed1:", theta
    call uniform_random(theta,0.d0,2.d0*pi,seed2)
    print*, "From seed2:", theta
    call uniform_random(theta,0.d0,2.d0*pi,seed3)
    print*, "From seed3:", theta

end program 
