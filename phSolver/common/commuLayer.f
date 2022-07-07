c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine layerCommuAssembly(global, rtemp, ilwork, n)
c
         include "common.h"
         include "auxmpi.h"
c
         dimension global(nshg,n),
     &             rtemp(maxfront*n,maxtask),
     &             ilwork(nlwork)
c
         real*8  locFlag, revFlag
c
         integer itkbeg, jdl, iacc, numseg, isgbeg, icsg, icid,
     &           lenseg, itemp, lfront, j, idof
c
         numtask = ilwork(1)
         itkbeg=1 ! slot in ilwork
         jdl=0    ! slot in rtemp
c
         do j=1,numtask        ! loop over all tasks
c
c.... total number of nodes involved in this task (lfront)
c
            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)
            lfront = 0
            do is = 1,numseg
              lenseg = ilwork (itkbeg + 4 + 2*is)
              lfront = lfront + lenseg
            enddo
c
            if(iacc.eq.1) then
               jdl=jdl+1  ! keep track of order of rtemp's
c
c.... only add the data from growth curve once
c
               itemp = 1
               do is = 1,numseg
                 isgbeg = ilwork (itkbeg + 3 + 2*is)
                 lenseg = ilwork (itkbeg + 4 + 2*is)
                 do icsg = 1,lenseg
                   icid  = isgbeg + icsg - 1
                   locFlag = global(icid,4) ! if it has been updated
                   if(lfront .gt. 0) then
                     revFlag = rtemp(3*lfront+itemp,jdl)
                   else
                     revFlag = 0.0;
                   endif
c
c.... if received from a growth curve and it's not been updated
c
                   if(revFlag .gt. 0.5 .and. locFlag .lt. 0.5) then
                     global(icid,1) = rtemp (itemp,jdl)          ! x
                     global(icid,2) = rtemp (lfront+itemp,jdl)   ! y
                     global(icid,3) = rtemp (2*lfront+itemp,jdl) ! z
                     global(icid,4) = revFlag
                   endif
                   itemp = itemp + 1
                 enddo ! end within each segment
               enddo ! end segments
            endif ! end of receive (iacc=1)
            itkbeg = itkbeg + 4 + 2*numseg
         enddo ! end tasks
c
c.... end
c
        return
        end
c
