c----------------------------------------------------------------------
c
c This subroutine satisfies the flow BCs on the actual_vi array. 
c actual_vi contatins the nodal values of the interface velocity. The 
c structure `of this code is derived from itrbc.f. This subroutine is 
c called from if_velocity.f. After the modification in this subroutine, 
c the actual_vi  array is used to set the mesh velocity at the interface 
c in if_velocity.f.
c 
c INPUT:
c  actual_vi	(nshg,nsd)   : Nodal velocity at the DG interface 
c  iBC	(nshg)        : Boundary Condition Code
c  BC	(nshg,ndofBC) : boundary condition constraint parameters
c
c
c OUTPUT:
c  actual_vi	(nshg,nsd)   : Adjusted actual_vi values based on flow BCs 
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
        subroutine itrBCvi (actual_vi, iBC, BC)
c
        include "common.h"
c
        dimension actual_vi(nshg,nsd),             iBC(nshg),
     &            BC(nshg,ndofBC)
          
c.... ------------------------->  Velocity  <--------------------------
c.... 3D
c
c.... x1-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 1)
            actual_vi(:,1) =  BC(:,3)  - BC(:,4) * actual_vi(:,2)
     &                         - BC(:,5) * actual_vi(:,3)
          endwhere
c
c.... x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 2)
            actual_vi(:,2) = BC(:,3)  - BC(:,4) * actual_vi(:,1)
     &                        - BC(:,5) * actual_vi(:,3)
          endwhere
c
c.... x1-velocity and x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 3)
            actual_vi(:,1) =  BC(:,3)  - BC(:,4) * actual_vi(:,3)
            actual_vi(:,2) =  BC(:,5)  - BC(:,6) * actual_vi(:,3)
          endwhere
c
c.... x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 4)
            actual_vi(:,3) = BC(:,3) - BC(:,4) * actual_vi(:,1)
     &                       - BC(:,5) * actual_vi(:,2)
          endwhere
c
c.... x1-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 5)
            actual_vi(:,1) = BC(:,3) - BC(:,4) * actual_vi(:,2)
            actual_vi(:,3) = BC(:,5) - BC(:,6) * actual_vi(:,2)
          endwhere
c
c.... x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 6)
            actual_vi(:,2) = BC(:,3)  - BC(:,4) * actual_vi(:,1)
            actual_vi(:,3) = BC(:,5)  - BC(:,6) * actual_vi(:,1)
          endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 7)
            actual_vi(:,1) =  BC(:,3)
            actual_vi(:,2) =  BC(:,4)
            actual_vi(:,3) =  BC(:,5) 
          endwhere
c
c
c.... end of velocity
c
        return
        end subroutine itrBCvi

