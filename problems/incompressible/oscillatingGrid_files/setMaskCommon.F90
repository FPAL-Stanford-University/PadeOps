      ist = max(ceiling(xmin/dx),gp%xst(1))
      ien = min(floor(  xmax/dx),gp%xen(1))
      ist = ist - gp%xst(1) + 1
      ien = ien - gp%xst(1) + 1
      
      jst = max(ceiling(ymin/dy),gp%xst(2))
      jen = min(floor(  ymax/dy),gp%xen(2))
      jst = jst - gp%xst(2) + 1
      jen = jen - gp%xst(2) + 1
      
      do k = kst, ken
        do j = jst, jen
          do i = ist, ien
            mask(i,j,k) = one 
          end do
        end do
      end do
