  do j=1,n2
     do i=1,n1
        call dgbtrs('N',this%n,7,7,1,d15,this%n,ipiv,y(i,j,:),this%n,info)
     enddo
  enddo
