        subroutine bc3BDgElas (iBC,  BC, elasBDiag, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the block-diagonal preconditioning
c   matrix for 3D elements.
c
c input:
c  iBC    (nshg)        : boundary condition code
c  BC     (nshg,ndofBC) : Dirichlet BC constraint parameters
c  elasBDiag   (nshg,nelas,nelas) : preconditionning matrix before BC
c                                   (only upper part)
c
c output:
c  elasBDiag   (nshg,nelas,nelas) : preconditionning matrix after BC
c                                   is satisfied
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),                     BC(nshg,4),
     &            elasBDiag(nshg,nelas,nelas),   ilwork(nlwork),
     &            iper(nshg)
c
c.... displacement
c
c.... x1-displacement  
c
        where (btest(iBC,14))
          elasBDiag(:,1,1) = one
          elasBDiag(:,1,2) = zero
          elasBDiag(:,1,3) = zero
          elasBDiag(:,2,1) = zero
          elasBDiag(:,3,1) = zero
        endwhere
c
c.... x2-displacement  
c
        where (btest(iBC,15))
          elasBDiag(:,2,1) = zero
          elasBDiag(:,2,2) = one
          elasBDiag(:,2,3) = zero
          elasBDiag(:,1,2) = zero
          elasBDiag(:,3,2) = zero
        endwhere
c
c.... x3-displacement  
c
        where (btest(iBC,16))
          elasBDiag(:,3,1) = zero
          elasBDiag(:,3,2) = zero
          elasBDiag(:,3,3) = one
          elasBDiag(:,1,3) = zero
          elasBDiag(:,2,3) = zero
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
              elasBDiag(isgbeg:isgend,:,:) = zero
              elasBDiag(isgbeg:isgend,1,1) = one
              elasBDiag(isgbeg:isgend,2,2) = one
              elasBDiag(isgbeg:isgend,3,3) = one
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
