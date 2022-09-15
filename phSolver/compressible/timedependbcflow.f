      subroutine timeDependBCFlow(x, iBC, BC_flow)
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
      integer   ftag1, ftag2, ftag3, ftag4, ftag5
      integer   etag1, etag2, etag3, etag4
      integer   maxDir(1)
      real*8    inormal(3)
      integer   answer
      real*8    n1, n2, n3

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
        answer = 0

        do i = 1, numnp
        
            if(m2gClsfcn(i,1).eq.0)
              call error('timedependbcflow, vertex slipBC not supported')

c Prescribed BC on the faces and edges of the bullet
            if( (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag1) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag2) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag3) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag4) .or.
          &   (m2gClsfcn(i,1).eq.2.and.m2gclsfcn(i,2).eq.ftag5) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag1) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag2) .or.
          &   (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag3))
              then

              if(m2gClsfcn(i,1).eq.1) then
                e_or_v_tag = m2gClsfcn(i,2) 
                if (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag1) then
                  f_tag = ftag1 
                endif
                if (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag2) then
                  f_tag = ftag2
                endif 
                if (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag3) then
                  f_tag = ftag3
                endif
              else
                f_tag = m2gClsfcn(i,2)  
              endif

c ... get unit normal
              call core_get_surf_normal(m2gClsfcn(i,1), f_tag, e_or_v_tag,
     &            m2gParCoord(i,1), m2gParCoord(i,2), n1, n2, n3)

              inormal(1) = n1
              inormal(2) = n2
              inormal(3) = n3

              maxDir = maxloc(abs(inormal(:)))
              select case(maxDir(1))
c ... if x of normal is the max
              case(1)
              iBC(i) = ibset(iBC(i),3)
              iBC(i) = ibclr(iBC(i),4)
              iBC(i) = ibclr(iBC(i),5)
              BC_flow(i,1) = 0d0 / inormal(1)
              BC_flow(i,2) = inormal(2) / inormal(1)
              BC_flow(i,3) = inormal(3) / inormal(1)
c ... if y of normal is the max
              case(2)
              iBC(i) = ibclr(iBC(i),3)
              iBC(i) = ibset(iBC(i),4)
              iBC(i) = ibclr(iBC(i),5)
              BC_flow(i,1) = 0d0 / inormal(2)
              BC_flow(i,2) = inormal(1) / inormal(2)
              BC_flow(i,3) = inormal(3) / inormal(2)
c ... if z of normal is the max
              case(3)
              iBC(i) = ibclr(iBC(i),3)
              iBC(i) = ibclr(iBC(i),4)
              iBC(i) = ibset(iBC(i),5)
              BC_flow(i,1) = 0d0 / inormal(3)
              BC_flow(i,2) = inormal(1) / inormal(3)
              BC_flow(i,3) = inormal(2) / inormal(3)
              end select
            else if (m2gClsfcn(i,1).eq.1.and.m2gclsfcn(i,2).eq.etag4) then
              f_tag = ftag4
              call core_get_surf_normal(m2gClsfcn(i,1), f_tag, etag4,
     &            m2gParCoord(i,1), m2gParCoord(i,2), n1, n2, n3)
              inormal(1) = n1
              inormal(2) = n2
              inormal(3) = n3

              maxDir = maxloc(abs(inormal(:)))
              select case(maxDir(1))
c... if x of normal is the max
              case (1)
              call error('timedepenfbcflow, unit normal 
     &            on back edge nb1 is max, which is not allowed')
c... if y of normal is the max
              case (2)
              iBC(i) = ibset(iBC(i), 3)
              iBC(i) = ibset(iBC(i), 4)
              iBC(i) = ibclr(iBC(i), 5)
              BC_flow(i,1) = 0d0
              BC_flow(i,2) = 0d0
              BC_flow(i,3) = 0d0 / inormal(2)
              BC_flow(i,4) = inormal(3) / inormal(2)
c... if z of normal is the max
              case (3)
              iBC(i) = ibset(iBC(i), 3)
              iBC(i) = ibclr(iBC(i), 4)
              iBC(i) = ibset(iBC(i), 5)
              BC_flow(i,1) = 0d0
              BC_flow(i,2) = 0d0
              BC_flow(i,3) = 0d0 / inormal(3)
              BC_flow(i,4) = inormal(2) / inormal(3)
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

