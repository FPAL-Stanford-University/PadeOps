  if (periodicBCs(1)) then
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if
  
  if (periodicBCs(2)) then
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if
  
  if (periodicBCs(3)) then
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2),  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2),  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if

  if (periodicBCs(1) .and. periodicBCs(2)) then
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if

  if (periodicBCs(1) .and. periodicBCs(3)) then
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if

  if (periodicBCs(2) .and. periodicBCs(3)) then
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if

  if (periodicBCs(1) .and. periodicBCs(2) .and. periodicBCs(3)) then
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2) + Ly,  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3) + Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2) - Ly,  this%renderModeData(:,3) - Lz, &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if
