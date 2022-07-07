      subroutine shpMeasure (x,     ien,        shp, shgl,
     &                       meshq, tempVolume, errorcount)
c---------------------------------------------------------
c This subroutine calculate the shape quality
c based on the mean ratio method.
c
c---------------------------------------------------------
c
        include "common.h"
c
        real*8    meshq(npro)
        dimension x(numnp, nsd),        ien(npro, nshl),
     &            xl(npro, nenl, nsd)
c
        dimension tetEdge(npro, 6, nsd),
     &            tetSumEdge(npro)
c
        dimension shp(nshl, ngauss),      sgn(npro,nshl),
     &            shgl(nsd,nshl,ngauss),  shg(npro,nshl,nsd),
     &            shdrv(npro,nsd,nshl),   shape(npro,nshl),
     &            tempVolume(npro)
c
        real*8, dimension(:),     pointer :: WdetJ
        real*8, dimension(:,:,:), pointer :: dxidx
c
        integer errorcount(2)
c
c---------------------------------------------------------
c
       call localx(x,      xl,     ien,    nsd,    'gather  ')
c
       tetEdge(:,1,:) = xl(:,1,:) - xl(:,4,:)
       tetEdge(:,2,:) = xl(:,2,:) - xl(:,4,:)
       tetEdge(:,3,:) = xl(:,3,:) - xl(:,4,:)
       tetEdge(:,4,:) = xl(:,1,:) - xl(:,3,:)
       tetEdge(:,5,:) = xl(:,2,:) - xl(:,3,:)
       tetEdge(:,6,:) = xl(:,1,:) - xl(:,2,:)
c
c... try to use Jacobian determinant to calculate volume
c
       allocate (dxidx(npro,nsd,nsd))
       allocate (WdetJ(npro))


       tempVolume = zero
       if (ipord .gt. 1) then
         call getsgn(ien,sgn)
       endif
       do intp = 1, ngauss
         if (Qwt(lcsyst,intp) .eq. zero) cycle    ! precaution
         call getshp(shp, shgl, sgn, shape, shdrv)
         call e3metric(xl, shdrv, dxidx, shg, WdetJ)
         tempVolume = tempVolume + WdetJ
       enddo
       deallocate ( dxidx )
       deallocate ( WdetJ )

       
       
c
c... calculate the sum of length of all edges for tet only
c
       tetSumEdge(:) = tetEdge(:,1,1) * tetEdge(:,1,1)
     &               + tetEdge(:,1,2) * tetEdge(:,1,2)
     &               + tetEdge(:,1,3) * tetEdge(:,1,3)
     &               + tetEdge(:,2,1) * tetEdge(:,2,1)
     &               + tetEdge(:,2,2) * tetEdge(:,2,2)
     &               + tetEdge(:,2,3) * tetEdge(:,2,3)
     &               + tetEdge(:,3,1) * tetEdge(:,3,1)
     &               + tetEdge(:,3,2) * tetEdge(:,3,2)
     &               + tetEdge(:,3,3) * tetEdge(:,3,3)
     &               + tetEdge(:,4,1) * tetEdge(:,4,1)
     &               + tetEdge(:,4,2) * tetEdge(:,4,2)
     &               + tetEdge(:,4,3) * tetEdge(:,4,3)
     &               + tetEdge(:,5,1) * tetEdge(:,5,1)
     &               + tetEdge(:,5,2) * tetEdge(:,5,2)
     &               + tetEdge(:,5,3) * tetEdge(:,5,3)
     &               + tetEdge(:,6,1) * tetEdge(:,6,1)
     &               + tetEdge(:,6,2) * tetEdge(:,6,2)
     &               + tetEdge(:,6,3) * tetEdge(:,6,3)
c
       meshq(:) = 15552.0 * tempVolume(:) * abs(tempVolume(:))
     &          / ( tetSumEdge(:) * tetSumEdge(:) * tetSumEdge(:) )
c
       do i = 1, npro
          if (meshq(i) .gt. 1.0000001) then
             errorcount(1) = errorcount(1) + 1
          endif
          if (meshq(i) .lt. 0.0) then
             errorcount(2) = errorcount(2) + 1
          endif
       enddo
c
c.... return
c
       return
       end
