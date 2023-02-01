      call getStEnIndices(xmin,xmax,dx,gp%xst(1),gp%xen(1),ist,ien)
      call getStEnIndices(ymin,ymax,dy,gp%xst(2),gp%xen(2),jst,jen)
      call setMask(ist,ien,jst,jen,kst,ken,gp,mask,one) 
