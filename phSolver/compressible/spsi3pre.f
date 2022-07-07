      subroutine Spsi3pre (BDiag,  lhsK,  nvar, col, row)
c
c------------------------------------------------------------------------------
c This is the initialization routine for the Sparse-GMRES solver.
c It pre-preconditions the LHS mass matrix and sets up the 
c sparse preconditioners. (pre-preconditioning is block diagonal
c scaling). 
c
c input:
c     BDiag  (nshg,nvar,nvar)    : block diagonal scaling matrix 
c                                  which is already LU factored.
c     lhsK(nvar*nvar,nnz_tot) : sparse LHS mass matrix
c
c output:
c     lhsK(nvar*nvar,nnz_tot)  : pre-preconditioned (scaled) mass matrix
c
c Nahid Razmara, Spring 2000.	(Sparse Matrix)
c------------------------------------------------------------------------------
c
      use pointer_data

      include "common.h"
c
        integer col(nshg+1),  row(nnz*nshg)
        integer sparseloc, c, c1
        real*8 lhsK(nvar*nvar,nnz_tot)
c
      dimension  BDiag(nshg,nvar,nvar)

      
c
c.... block-diagonal pre-precondition LHS 
c
c     
c.... reduce by columns, (left block diagonal preconditioning)
c
c     lhsK  <-- inverse(L^tilde) lhsK     
c     
c
      if(nflow.ne.5) then
         write(*,*)' spsi3pre.f assumed nflow=5'
         stop
      endif
c
      do i = 1, nshg 
          do k = col(i), col(i+1)-1
c             
            do ic = 2, nvar    
              do j = 1, nvar
                id = (j-1) * nvar + ic
                do il = 1, ic-1
                  lhsK(id,k) = lhsK(id,k) 
     &                       - BDiag(i,ic,ic-il) * lhsK(id-il,k)
                enddo 
              enddo
            enddo
c            
          enddo
        enddo
            
        do i = 1, nshg
c
            do k = col(i), col(i+1)-1
                j = row(k)
c     
c.... reduce by rows, (right block diagonal preconditioning)
c
c     lhsK   <-- lhsK  inverse(U^tilde)
c     
                  do jj = 1, nvar
                       lhsK(jj,k)  = BDiag(j,1,1)*lhsK(jj,k) 
                  enddo
c              
                  do ic = 2, nvar    
                    do jj = 1, nvar
                      id = (ic-1) * nvar + jj
                      do il = 1, ic-1
                        iid = id - nvar * il
                        lhsK(id,k) = lhsK(id,k) 
     &                             - BDiag(j,ic-il,ic) * lhsK(iid,k) 
                      enddo
                      lhsK(id,k) = BDiag(j,ic,ic) * lhsK(id,k) 
                    enddo
                  enddo
c 
             enddo
        enddo
c     
      return

      end

      
