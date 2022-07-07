      subroutine bc3ResElas (iBC,  BC,  elasres, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the residual vector for 3D elements.
c
c input: 
c  iBC   (nshg)            : Boundary Condition Code
c  BC    (nshg,ndofBC)     : the boundary condition constraint parameters
c  elasres   (nshg,nelas)  : residual before BC is applied
c
c output:
c  elasres   (nshg,nelas)  : residual after satisfaction of BC
c
c----------------------------------------------------------------------
c
      include "common.h"
c     
      dimension iBC(nshg),             BC(nshg,4),  
     &          elasres(nshg,nelas),   ilwork(nlwork),
     &          iper(nshg)
c
c.... displacement
c
c.... x1-displacement
c
        where (btest(iBC,14)) 
          elasres(:,1) = zero
        endwhere
c
c.... x2-displacement
c
        where (btest(iBC,15)) 
          elasres(:,2) = zero
        endwhere
c
c.... x3-displacement
c
        where (btest(iBC,16)) 
          elasres(:,3) = zero
        endwhere
c
c.... parallel communication
c
      if(numpe.gt.1) then
c
c.... nodes treated on another processor are eliminated
c
        numtask = ilwork(1)
        itkbeg = 1

        do itask = 1, numtask

          iacc   = ilwork (itkbeg + 2)
          numseg = ilwork (itkbeg + 4)

          if (iacc .eq. 0) then
            do is = 1,numseg
              isgbeg = ilwork (itkbeg + 3 + 2*is)
              lenseg = ilwork (itkbeg + 4 + 2*is)
              isgend = isgbeg + lenseg - 1
              elasres(isgbeg:isgend,:) = zero
            enddo
          endif
          
          itkbeg = itkbeg + 4 + 2*numseg

        enddo
        endif
c
c.... return
c
        return
        end
c
