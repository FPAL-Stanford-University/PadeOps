program test_seed
    use seedGen
    implicit none
    integer :: seed1, seed3, seed2

    seed1 = get_seed_from_location(0.975000000000000d0, 0.225000000000000d0, 0.225000000000000d0, 1.D18)
    seed2 = get_seed_from_location(0.925000000000000d0, 0.225000000000000d0, 0.225000000000000d0, 1.D18)
    seed3 = get_seed_from_location(0.925000000000000d0, 0.225000000000000d0, 0.225000000000000d0, 2.D18)

    print*, seed1, seed2, seed3

end program 
