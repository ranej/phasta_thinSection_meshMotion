        subroutine itrBCElas (umesh, disp, iBC, BC, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This program satisfies the boundary conditions on the mesh-elastic BCs.
c
c input:
c  disp   (nshg,nelas)  : mesh-elastic displacement 
c  iBC    (nshg)        : Boundary Condition Code
c  BC     (nshg,4)      : boundary condition constraint parameters
c                         pass only 4 componenets to this subroutine
c
c output:
c  disp   (nshg,nelas)  : Adjusted mesh displacement value(s) corresponding
c                          to a constraint d.o.f.
c  
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension umesh(nshg,nelas),        disp(nshg,nelas),
     &            iBC(nshg),                BC(nshg,4)
c
        dimension ilwork(nlwork),           iper(nshg)
c
c.... -----------------> Mesh Elastic Velocity <-----------------------
c.... 3D
c
c.... x1-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 1)
            disp(:,1)  = BC(:,1) - BC(:,2) * disp(:,2)
     &                           - BC(:,3) * disp(:,3)
            umesh(:,1) = disp(:,1) / Delt(1)
          endwhere
c
c.... x2-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 2)
            disp(:,2)  = BC(:,1) - BC(:,2) * disp(:,1)
     &                           - BC(:,3) * disp(:,3)
            umesh(:,2) = disp(:,2) / Delt(1)
          endwhere
c
c.... x1-velocity and x2-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 3)
            disp(:,1)  = BC(:,1) - BC(:,2) * disp(:,3)
            disp(:,2)  = BC(:,3) - BC(:,4) * disp(:,3)
            umesh(:,1) = disp(:,1) / Delt(1)
            umesh(:,2) = disp(:,2) / Delt(1)
          endwhere
c
c.... x3-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 4)
            disp(:,3)  = BC(:,1) - BC(:,2) * disp(:,1)
     &                           - BC(:,3) * disp(:,2)
            umesh(:,3) = disp(:,3) / Delt(1)
          endwhere
c
c.... x1-velocity and x3-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 5)
            disp(:,1)  = BC(:,1) - BC(:,2) * disp(:,2)
            disp(:,3)  = BC(:,3) - BC(:,4) * disp(:,2)
            umesh(:,1) = disp(:,1) / Delt(1)
            umesh(:,3) = disp(:,3) / Delt(1)
          endwhere
c
c.... x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 6)
            disp(:,2)  = BC(:,1) - BC(:,2) * disp(:,1)
            disp(:,3)  = BC(:,3) - BC(:,4) * disp(:,1)
            umesh(:,2) = disp(:,2) / Delt(1)
            umesh(:,3) = disp(:,3) / Delt(1)
          endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,14,3) .eq. 7)
            disp(:,1)  = BC(:,1)
            disp(:,2)  = BC(:,2)
            disp(:,3)  = BC(:,3)
            umesh(:,1) = disp(:,1) / Delt(1)
            umesh(:,2) = disp(:,2) / Delt(1)
            umesh(:,3) = disp(:,3) / Delt(1)
          endwhere
c
c       endif
c
c.... end of velocity
c
c.... communications
c 
        if (numpe > 1) then
           call commu (disp, ilwork, nelas, 'out')
        endif
c
c       slave has masters value, for abc we need to rotate it
c
        if(iabc==1) then        !are there any axisym bc's
           call rotabc(disp, iBC, 'out')
        endif
c
c.... return
c
        return
        end

