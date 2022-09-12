      subroutine timeDependBCElas(x, iBC, BC_flow)
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      use m2gfields ! read m2g fields
c
      include "common.h"
      include "mpif.h"
c
c
      real*8    x(numnp,nsd)
      dimension iBC(nshg), BC_flow(nshg,4)
      integer   casenumber

c.... precribed BC for slip surfaces 
      real*8    ang, x_center(nshg,3), rotx(nshg,3)
      integer   ftag1, ftag2, ftag3, ftag4, ftag5
      integer   etag1, etag2, etag3, etag4
      integer   maxDir(1)
      real*8    pnormal(nshg,3), inormal(nshg,3)
      real*8    m2gD(nshg,3)
      real*8    cos_theta, sin_theta
      integer   counter

      casenumber = 0

      if (tdbcflowcase .gt. 0) then
        if (myrank .eq. master) then
          write(*,*) "use time dependent flow bc case :", tdbcflowcase
        endif
        casenumber = tdbcflowcase
      endif
c
c.... Update BC value based on geom and iBC
c

c
c.... case 1
c
      if ( casenumber .eq. 1 ) then

c... Tags for bullet                

        ftag1 = 145 ! nose
        ftag2 = 149 ! mid1
        ftag3 = 154 ! mid2
        ftag4 = 159 ! mid3
        ftag5 = 164 ! back

        etag1 = 137 ! f12
        etag2 = 139 ! f23
        etag3 = 141 ! f34
        etag4 = 143 ! f45

          do i = 1, numnp

            x_center(i,1) = x(i,1)
            x_center(i,2) = 0d0
            x_center(i,3) = 0d0

              call core_is_in_closure(m2gClsfcn(i,1), m2gClsfcn(i,2),
     &                            1,              etag4,
     &                            answer)

c Prescribed BC on the faces and edges of the bullet
            if( (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag1) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag2) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag3) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag4) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag5) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag1) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag2) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag3) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag4) .or.
               answer.ne.0)
              then


c ... calc normal call to get this value
          
c...  For back face (ftag5), at 0 deg AoA, normal will be x aligned
c..   Normal calculations will need to change for non-zero AoA

              if( m2gClsfcn(i,2) .eq. ftag5) then
                pnormal(i,1) = 1
                pnormal(i,2) = 0d0
                pnormal(i,3) = 0d0
              else
                pnormal(i,:) = x(i,:) - x_center(i,)   
              endif

              tmp = sqrt( pnormal(i,1)*pnormal(i,1)
          &             + pnormal(i,2)*pnormal(i,2)
          &             + pnormal(i,3)*pnormal(i,3) )
              inormal(i,:) = pnormal(i,:) / (tmp)

              maxDir = maxloc(abs(inormal(i,:)))
              select case(maxDir(1))
c ... if x of normal is the max
              case(1)
              iBC(i) = ibset(iBC(i),3)
              iBC(i) = ibclr(iBC(i),4)
              iBC(i) = ibclr(iBC(i),5)
              BC_flow(i,1) = 0d0 / pnormal(i,1)
              BC_flow(i,2) = pnormal(i,2) / pnormal(i,1)
              BC_flow(i,3) = pnormal(i,3) / pnormal(i,1)
c ... if y of normal is the max
              case(2)
              iBC(i) = ibclr(iBC(i),3)
              iBC(i) = ibset(iBC(i),4)
              iBC(i) = ibclr(iBC(i),5)
              BC_flow(i,1) = 0d0 / pnormal(i,2)
              BC_flow(i,2) = pnormal(i,1) / pnormal(i,2)
              BC_flow(i,3) = pnormal(i,3) / pnormal(i,2)
c ... if z of normal is the max
              case(3)
              iBC(i) = ibclr(iBC(i),3)
              iBC(i) = ibclr(iBC(i),4)
              iBC(i) = ibset(iBC(i),5)
              BC_flow(i,1) = 0d0 / pnormal(i,3)
              BC_flow(i,2) = pnormal(i,1) / pnormal(i,3)
              BC_flow(i,3) = pnormal(i,2) / pnormal(i,3)
              end select
            endif
          enddo
        endif



c
c.... end case 1
c

      return
      end
c

