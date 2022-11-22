if (periodicBCs(1)) then
  call message(2,'Adding x-periodic contribution to velocity')
  if (lastY > 0 .and. lastZ > 0) then
    call this%renderLocalVelocity(&
      [this%x,     haloBuffY(:,1),  haloBuffZ(:,1) ] + Lx,  &
      [this%y,     haloBuffY(:,2),  haloBuffZ(:,2) ],  &
      [this%z,     haloBuffY(:,3),  haloBuffZ(:,3) ],  &
      [this%kx,    haloBuffY(:,4),  haloBuffZ(:,4) ],  &
      [this%ky,    haloBuffY(:,5),  haloBuffZ(:,5) ],  &
      [this%kz,    haloBuffY(:,6),  haloBuffZ(:,6) ],  &
      [this%uhatR, haloBuffY(:,7),  haloBuffZ(:,7) ],  &
      [this%uhatI, haloBuffY(:,8),  haloBuffZ(:,8) ],  &
      [this%vhatR, haloBuffY(:,9),  haloBuffZ(:,9) ],  &
      [this%vhatI, haloBuffY(:,10), haloBuffZ(:,10)], &
      [this%whatR, haloBuffY(:,11), haloBuffZ(:,11)], &
      [this%whatI, haloBuffY(:,12), haloBuffZ(:,12)])
    call this%renderLocalVelocity(&
      [this%x,     haloBuffY(:,1),  haloBuffZ(:,1) ] - Lx,  &
      [this%y,     haloBuffY(:,2),  haloBuffZ(:,2) ],  &
      [this%z,     haloBuffY(:,3),  haloBuffZ(:,3) ],  &
      [this%kx,    haloBuffY(:,4),  haloBuffZ(:,4) ],  &
      [this%ky,    haloBuffY(:,5),  haloBuffZ(:,5) ],  &
      [this%kz,    haloBuffY(:,6),  haloBuffZ(:,6) ],  &
      [this%uhatR, haloBuffY(:,7),  haloBuffZ(:,7) ],  &
      [this%uhatI, haloBuffY(:,8),  haloBuffZ(:,8) ],  &
      [this%vhatR, haloBuffY(:,9),  haloBuffZ(:,9) ],  &
      [this%vhatI, haloBuffY(:,10), haloBuffZ(:,10)], &
      [this%whatR, haloBuffY(:,11), haloBuffZ(:,11)], &
      [this%whatI, haloBuffY(:,12), haloBuffZ(:,12)])
  else if (lastY > 0) then
    call this%renderLocalVelocity([this%x,     haloBuffY(:,1)] + Lx,  &
                                  [this%y,     haloBuffY(:,2)],  &
                                  [this%z,     haloBuffY(:,3)],  &
                                  [this%kx,    haloBuffY(:,4)],  &
                                  [this%ky,    haloBuffY(:,5)],  &
                                  [this%kz,    haloBuffY(:,6)],  &
                                  [this%uhatR, haloBuffY(:,7)],  &
                                  [this%uhatI, haloBuffY(:,8)],  &
                                  [this%vhatR, haloBuffY(:,9)],  &
                                  [this%vhatI, haloBuffY(:,10)], &
                                  [this%whatR, haloBuffY(:,11)], &
                                  [this%whatI, haloBuffY(:,12)])
    call this%renderLocalVelocity([this%x,     haloBuffY(:,1)] - Lx,  &
                                  [this%y,     haloBuffY(:,2)],  &
                                  [this%z,     haloBuffY(:,3)],  &
                                  [this%kx,    haloBuffY(:,4)],  &
                                  [this%ky,    haloBuffY(:,5)],  &
                                  [this%kz,    haloBuffY(:,6)],  &
                                  [this%uhatR, haloBuffY(:,7)],  &
                                  [this%uhatI, haloBuffY(:,8)],  &
                                  [this%vhatR, haloBuffY(:,9)],  &
                                  [this%vhatI, haloBuffY(:,10)], &
                                  [this%whatR, haloBuffY(:,11)], &
                                  [this%whatI, haloBuffY(:,12)])
  else if (lastZ > 0) then
    call this%renderLocalVelocity([this%x,     haloBuffZ(:,1)] + Lx,  &
                                  [this%y,     haloBuffZ(:,2)],  &
                                  [this%z,     haloBuffZ(:,3)],  &
                                  [this%kx,    haloBuffZ(:,4)],  &
                                  [this%ky,    haloBuffZ(:,5)],  &
                                  [this%kz,    haloBuffZ(:,6)],  &
                                  [this%uhatR, haloBuffZ(:,7)],  &
                                  [this%uhatI, haloBuffZ(:,8)],  &
                                  [this%vhatR, haloBuffZ(:,9)],  &
                                  [this%vhatI, haloBuffZ(:,10)], &
                                  [this%whatR, haloBuffZ(:,11)], &
                                  [this%whatI, haloBuffZ(:,12)])
    call this%renderLocalVelocity([this%x,     haloBuffZ(:,1)] - Lx,  &
                                  [this%y,     haloBuffZ(:,2)],  &
                                  [this%z,     haloBuffZ(:,3)],  &
                                  [this%kx,    haloBuffZ(:,4)],  &
                                  [this%ky,    haloBuffZ(:,5)],  &
                                  [this%kz,    haloBuffZ(:,6)],  &
                                  [this%uhatR, haloBuffZ(:,7)],  &
                                  [this%uhatI, haloBuffZ(:,8)],  &
                                  [this%vhatR, haloBuffZ(:,9)],  &
                                  [this%vhatI, haloBuffZ(:,10)], &
                                  [this%whatR, haloBuffZ(:,11)], &
                                  [this%whatI, haloBuffZ(:,12)])
  else
    call this%renderLocalVelocity(this%x + Lx, this%y, this%z, this%kx, &
      this%ky, this%kz, this%uhatR, this%uhatI, this%vhatR, this%vhatI, &
      this%whatR, this%whatI)
    call this%renderLocalVelocity(this%x - Lx, this%y, this%z, this%kx, &
      this%ky, this%kz, this%uhatR, this%uhatI, this%vhatR, this%vhatI, &
      this%whatR, this%whatI)
  end if
end if
