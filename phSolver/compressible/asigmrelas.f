        subroutine AsIGMRElas (x,     disp,    shp,       shgl, 
     &                         ien,   elasres, elasBDiag, Estiff,
     &                         meshq, meshV )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension x(numnp,nsd),             disp(numnp,nsd), 
     &            shp(nshl,MAXQPT),         sgn(npro,nshl), 
     &            shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),           elasres(nshg,nelas),
     &            elasBDiag(nshg,nelas,nelas)
c
        dimension xl(npro,nenl,nsd),        displ(npro,nenl,nsd),
     &            Eforce(npro,nshl,nelas),      
     &            elasBDiagl(npro,nshl,nelas,nelas),
     &            Estiff(npro,ndofelas,ndofelas)
c   
        real*8    meshq(npro),              meshV(npro)
c       
c.... create the matrix of mode signs for the hierarchic basis 
c     function
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call localx(disp,   displ,  ien,    nsd,    'gather  ')
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        Eforce     = zero
        elasBDiagl = zero
        elasCo     = zero
c
        call e3Elas  (xl,     displ,  meshq,     meshV,  
     &                shp,    shgl,   sgn,
     &                Estiff, Eforce, elasBDiagl ) 
c
c.... assemble the residual and modified residual
c
        call local (elasres,    Eforce,     ien,    nelas,  'scatter ')
c
c.... extract and assemble the Block-Diagonal (see note in elmgmr, line 280)
c
        if (iprec .ne. 0) then 
           do i = 1, nshl
              do j = 1, nelas
                 i0 = (i - 1) * nelas + j
                 do k = 1, nelas
                    j0 = (i - 1) * nelas + k
                    elasBDiagl(:,i,j,k) = Estiff(:,i0,j0)
                 enddo
              enddo
           enddo
           call local (elasBDiag,  elasBDiagl,
     &                 ien,        nelas*nelas, 'scatter ')
        endif
c        
c.... end
c
        return
        end
