if (lastY > 0 .and. lastZ > 0) then
  call message(2,'Adding halo mode contribution to velocity')
  call this%renderLocalVelocity([haloBuffY(:,1),  haloBuffZ(:,1)], &
                                [haloBuffY(:,2),  haloBuffZ(:,2)], &
                                [haloBuffY(:,3),  haloBuffZ(:,3)], &
                                [haloBuffY(:,4),  haloBuffZ(:,4)], &
                                [haloBuffY(:,5),  haloBuffZ(:,5)], &
                                [haloBuffY(:,6),  haloBuffZ(:,6)], &
                                [haloBuffY(:,7),  haloBuffZ(:,7)], &
                                [haloBuffY(:,8),  haloBuffZ(:,8)], &
                                [haloBuffY(:,9),  haloBuffZ(:,9)], &
                                [haloBuffY(:,10), haloBuffZ(:,10)], &
                                [haloBuffY(:,11), haloBuffZ(:,11)], &
                                [haloBuffY(:,12), haloBuffZ(:,12)])
else if (lastY > 0) then
  call message(2,'Adding halo mode contribution to velocity')
  haloBuffY => this%haloBuffY(1:lastY,:)
  call this%renderLocalVelocity(haloBuffY(:,1), haloBuffY(:,2), &
    haloBuffY(:,3), haloBuffY(:,4), haloBuffY(:,5), &
    haloBuffY(:,6), haloBuffY(:,7), haloBuffY(:,8), &
    haloBuffY(:,9), haloBuffY(:,10), haloBuffY(:,11), &
    haloBuffY(:,12))
else if (lastZ > 0) then
  call message(2,'Adding halo mode contribution to velocity')
  haloBuffZ => this%haloBuffZ(1:lastZ,:)
  call this%renderLocalVelocity(haloBuffZ(:,1), haloBuffZ(:,2), &
    haloBuffZ(:,3), haloBuffZ(:,4), haloBuffZ(:,5), &
    haloBuffZ(:,6), haloBuffZ(:,7), haloBuffZ(:,8), &
    haloBuffZ(:,9), haloBuffZ(:,10), haloBuffZ(:,11), &
    haloBuffZ(:,12))
end if
