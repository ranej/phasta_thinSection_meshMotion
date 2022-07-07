        subroutine e3Elas (xl,        displ,   meshq,    meshV,
     &                     shp,       shgl,    sgn,
     &                     Estiff,    Eforce,  elasBDiagl )
c
c----------------------------------------------------------------------
c
c This routine compute the stiffness of the mesh-elastic equation:
c    stiffness * displacement = force
c
c input:
c    xl(npro,nenl,nsd):     local coordinates of the vertices
c    elasCo(2):             two coefficients used to modify elastic parameters
c                           to adjust the stiffness of this element
c    displ(npro,nenl,nsd):  local displacement of the vertices
c    shgl(nsd,nshl,ngauss): element local grad-shape-function N_{a,xi}
c
c output:
c    Estiff(npro,nshl*nelas,nshl*nelas): element stiffness matrix (LHS)
c    Eforce(npro,nshl,nelas):            element force matrix (RHS)
c                                        the predicted solustion is zero
c    elasBDiagl(npro,nshl,nelas):        block-preconditioner (Jacobi)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension xl(npro,nenl,nsd),       displ(npro,nenl,nsd),
     &            shp(nshl, ngauss),       sgn(npro,nshl),
     &            shgl(nsd,nshl,ngauss),   shg(npro,nshl,nsd),
     &            shdrv(npro,nsd,nshl),    shape(npro,nshl)
c
C        dimension dxidxelas(npro,nsd,nsd), WdetJelas(npro)
c
        dimension Estiff(npro,ndofelas,ndofelas),
     &            Eforce(npro,nshl,nelas),
     &            elasBDiagl(npro,nshl,nelas,nelas)
c
        real*8, dimension(:),     pointer :: WdetJelas
        real*8, dimension(:,:,:), pointer :: Kelas, dxidxelas
        real*8  lamda(npro),  mu(npro), meshq(npro), meshV(npro)
        real*8  youngMod(npro), poisnRat(npro)
c
c
c----------------------------------------------------------------------
c
c        if (associated(dxidxelas)) deallocate(dxidxelas)
        allocate (dxidxelas(npro,nsd,nsd))
        allocate (WdetJelas(npro))
        allocate (Kelas(npro,nelas,nelas))
c
c.... initial setup
c
c.... modify Young's Modulus and Poisson Ratio
c
        if(datelas_volume_YM .eq. 1) then
          youngMod(:) = datelas(1,1) / meshV(:)
        else
          youngMod(:) = datelas(1,1)
        endif
        poisnRat(:) = datelas(1,2)
c
c.... enlarge the stiffness for wedge elements
c
        if (lcsyst .eq. 3) then ! wedge
          youngMod(:) = blfactor * youngMod(:)
        endif
c
        lamda(:) = youngMod(:) * poisnRat(:) /
     &         (1.0+poisnRat(:)) / (1.0-2.0*poisnRat(:))
        mu(:)    = youngMod(:) / 2.0 / (1.0+poisnRat(:))
c
c.... loop through the integration points
c
        do intp = 1, ngauss
c
c.... if Det. .eq. 0, do not include this point
c
        if (Qwt(lcsyst,intp) .eq. zero) cycle    ! precaution
c
        call getshp(shp,          shgl,      sgn,
     &              shape,        shdrv)
c
        call e3metric (xl, shdrv, dxidxelas, shg, WdetJelas)
c
c.... loop over vertices in the element
c
        do i = 1,nenl
          do j = 1,nenl
            Kelas = zero
c
c.... compute BDB
c
            Kelas(:,1,1) = WdetJelas(:)
     &               * ( shg(:,i,2) * shg(:,j,2) * mu(:)
     &               +   shg(:,i,3) * shg(:,j,3) * mu(:)
     &               +   shg(:,i,1) * shg(:,j,1) * (lamda(:) + 2*mu(:)))
            Kelas(:,1,2) = WdetJelas(:)
     &               * ( shg(:,i,1) * shg(:,j,2) * lamda(:)
     &               +   shg(:,i,2) * shg(:,j,1) * mu(:)    )
            Kelas(:,1,3) = WdetJelas(:)
     &               * ( shg(:,i,1) * shg(:,j,3) * lamda(:)
     &               +   shg(:,i,3) * shg(:,j,1) * mu(:)    )
            Kelas(:,2,1) = WdetJelas(:)
     &               * ( shg(:,i,2) * shg(:,j,1) * lamda(:)
     &               +   shg(:,i,1) * shg(:,j,2) * mu(:)    )
            Kelas(:,2,2) = WdetJelas(:)
     &               * ( shg(:,i,1) * shg(:,j,1) * mu(:)
     &               +   shg(:,i,3) * shg(:,j,3) * mu(:)
     &               +   shg(:,i,2) * shg(:,j,2) * (lamda(:) + 2*mu(:)))
            Kelas(:,2,3) = WdetJelas(:)
     &               * ( shg(:,i,2) * shg(:,j,3) * lamda(:)
     &               +   shg(:,i,3) * shg(:,j,2) * mu(:)    )
            Kelas(:,3,1) = WdetJelas(:)
     &               * ( shg(:,i,3) * shg(:,j,1) * lamda(:)
     &               +   shg(:,i,1) * shg(:,j,3) * mu(:)    )
            Kelas(:,3,2) = WdetJelas(:)
     &               * ( shg(:,i,3) * shg(:,j,2) * lamda(:)
     &               +   shg(:,i,2) * shg(:,j,3) * mu(:)    )
            Kelas(:,3,3) = WdetJelas(:)
     &               * ( shg(:,i,1) * shg(:,j,1) * mu(:)
     &               +   shg(:,i,2) * shg(:,j,2) * mu(:)
     &               +   shg(:,i,3) * shg(:,j,3) * (lamda(:) + 2*mu(:)))
c
c... Form the element stiffness matrix
c
            Estiff(:,3*(i-1)+1,3*(j-1)+1) =
     &      Estiff(:,3*(i-1)+1,3*(j-1)+1) + Kelas(:,1,1)
            Estiff(:,3*(i-1)+1,3*(j-1)+2) =
     &      Estiff(:,3*(i-1)+1,3*(j-1)+2) + Kelas(:,1,2)
            Estiff(:,3*(i-1)+1,3*(j-1)+3) =
     &      Estiff(:,3*(i-1)+1,3*(j-1)+3) + Kelas(:,1,3)
            Estiff(:,3*(i-1)+2,3*(j-1)+1) =
     &      Estiff(:,3*(i-1)+2,3*(j-1)+1) + Kelas(:,2,1)
            Estiff(:,3*(i-1)+2,3*(j-1)+2) =
     &      Estiff(:,3*(i-1)+2,3*(j-1)+2) + Kelas(:,2,2)
            Estiff(:,3*(i-1)+2,3*(j-1)+3) =
     &      Estiff(:,3*(i-1)+2,3*(j-1)+3) + Kelas(:,2,3)
            Estiff(:,3*(i-1)+3,3*(j-1)+1) =
     &      Estiff(:,3*(i-1)+3,3*(j-1)+1) + Kelas(:,3,1)
            Estiff(:,3*(i-1)+3,3*(j-1)+2) =
     &      Estiff(:,3*(i-1)+3,3*(j-1)+2) + Kelas(:,3,2)
            Estiff(:,3*(i-1)+3,3*(j-1)+3) =
     &      Estiff(:,3*(i-1)+3,3*(j-1)+3) + Kelas(:,3,3)
c
c.... Form the element force vector
c
            Eforce(:,i,1) = Eforce(:,i,1)
     &                     - Kelas(:,1,1)*displ(:,j,1)
     &                     - Kelas(:,1,2)*displ(:,j,2)
     &                     - Kelas(:,1,3)*displ(:,j,3)
            Eforce(:,i,2) = Eforce(:,i,2)
     &                     - Kelas(:,2,1)*displ(:,j,1)
     &                     - Kelas(:,2,2)*displ(:,j,2)
     &                     - Kelas(:,2,3)*displ(:,j,3)
            Eforce(:,i,3) = Eforce(:,i,3)
     &                     - Kelas(:,3,1)*displ(:,j,1)
     &                     - Kelas(:,3,2)*displ(:,j,2)
     &                     - Kelas(:,3,3)*displ(:,j,3)
c
          enddo
        enddo
c
        enddo ! end loop integration points
c
        deallocate ( dxidxelas )
        deallocate ( WdetJelas )
        deallocate ( Kelas )
c
c.... return
c
      return
      end
