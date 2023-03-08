subroutine togglePointer(this)
  class(enrichmentOperator), intent(inout), target :: this  
 
  this%activeIndex = mod(this%activeIndex + 1,2)
  
  ! Link pointers
  this%modeData => this%rawModeData(1:this%nmodes,  :,this%activeIndex)
  this%x        => this%rawModeData(1:this%nmodes,  1, this%activeIndex)
  this%y        => this%rawModeData(1:this%nmodes,  2, this%activeIndex)
  this%z        => this%rawModeData(1:this%nmodes,  3, this%activeIndex)
  this%kx       => this%rawModeData(1:this%nmodes,  4, this%activeIndex)
  this%ky       => this%rawModeData(1:this%nmodes,  5, this%activeIndex)
  this%kz       => this%rawModeData(1:this%nmodes,  6, this%activeIndex)
  this%uhatR    => this%rawModeData(1:this%nmodes,  7, this%activeIndex)
  this%uhatI    => this%rawModeData(1:this%nmodes,  8, this%activeIndex)
  this%vhatR    => this%rawModeData(1:this%nmodes,  9, this%activeIndex)
  this%vhatI    => this%rawModeData(1:this%nmodes, 10, this%activeIndex)
  this%whatR    => this%rawModeData(1:this%nmodes, 11, this%activeIndex)
  this%whatI    => this%rawModeData(1:this%nmodes, 12, this%activeIndex)
  if (this%isStratified) this%T => this%rawModeData(1:this%nmodes, 13, this%activeIndex)
end subroutine
