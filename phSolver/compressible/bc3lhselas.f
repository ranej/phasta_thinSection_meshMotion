      subroutine bc3LHSElas (iBC,  BC,  iens,  Estiff )
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of LHS Estiff matrix for a single 
c element.
c
c input:
c  iBC   (nshg)         : boundary condition code
c  BC    (nshg,nelas)      : Dirichlet BC constraint parameters
c  ien   (npro,nshl)    : ien array for this element
c  Estiff(npro,nshl*nelas,nshl*nelas) : element stiffness matrix (LHS) before BC
c  Eforce(npro,nshl,nelas)            : element force array (RHS) before BC
c
c output:
c  Estiff(npro,nshl*nelas,nshl*nelas) : LHS stiffness matrix after BC is satisfied
c  Eforce(npro,nshl,nelas)            : RHS force array after BC is satisfied
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),       ien(npro,nshl),
     &            BC(nshg,4),  Estiff(npro,ndofelas,ndofelas)
c
        integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
        ien=abs(iens)
c
c.... loop over elements
c
        do iel = 1, npro
c
c.... loop over number of shape functions for this element
c
        do inod = 1, nshl
c
c.... set up parameters
c
        in  = ien(iel,inod)
        if (ibits(iBC(in),14,3) .eq. 0)  cycle
        ioff = (inod - 1) * nelas
        i1 = ioff + 1
        i2 = ioff + 2
        i3 = ioff + 3
c
c.... displacement
c
c
c.... x1-displacement
c
        if ( btest(iBC(in),14) ) then
          Estiff(iel,i1,:) = zero
          Estiff(iel,:,i1) = zero
          Estiff(iel,i1,i1) = one
        endif
c
c.... x2-displacement
c
        if ( btest(iBC(in),15) ) then
          Estiff(iel,i2,:) = zero
          Estiff(iel,:,i2) = zero
          Estiff(iel,i2,i2) = one
        endif
c
c.... x3-displacement
c
        if ( btest(iBC(in),16) ) then
          Estiff(iel,i3,:) = zero
          Estiff(iel,:,i3) = zero
          Estiff(iel,i3,i3) = one
        endif
c
c.... end loop over shape functions (nodes)
c
      enddo
c
c.... end loop over elements
c
      enddo
c
c.... return
c
      return
      end
