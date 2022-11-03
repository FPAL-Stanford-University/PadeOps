subroutine runTests(this)
  class(enrichmentOperator), intent(inout) :: this
  call testMPIsendRecv()
end subroutine
