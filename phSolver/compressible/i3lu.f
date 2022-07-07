        subroutine i3LU (Diag, r, nvar, code)
c
c----------------------------------------------------------------------
c
c This routine performs a LU factorization/solve of a set of matrices,
c used for block diagonal preconditioning in the iterative driver.
 
c input:
c  Diag (nshg,nflow,nflow)  : block diagonal (symmetric storage)
c  r    (nshg,nflow)    : residual
c  code                 : operation code
c                           .eq. 'LU_Fact ', Cholesky Factor
c                           .eq. 'forward ', forward reduction
c                           .eq. 'backward', backward substitution
c                           .eq. 'product ', product Diag.r
c
c output:
c  Diag (nshg,nvar,nvar)  : LU decomp. of block diagonal
c  r    (nshg,nvar)    : reduced residual
c
c
c Note: the formulation used here is taken from Golub's "Matrix
c       Computations" Book (1984), pages 82-85 algorithm P5.1-3.
c
c       L and U overwrite Diag
c
c       The diagonal terms (i,i) of U are stored in inverted form.
c
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension Diag(nshg,nvar,nvar),        r(nshg,nvar)
c
        character*8 code
        
c
c.... perform LU decomposition with the Diagonal terms inverted
c
        if (code .eq. 'LU_Fact ') then
c       
        do i = 1, nvar
          if (i .gt. 1) then
            do j = i, nvar
              do k = 1, i-1
                Diag(:,i,j) = Diag(:,i,j)
     &                      - Diag(:,i,k) * Diag(:,k,j)
              enddo
            enddo
          endif
c      
          Diag(:,i,i) = one / Diag(:,i,i)
c
          if (i .lt. nvar) then
            do j = i+1, nvar
              if (i .gt. 1) then 
                do k = 1, i-1
                  Diag(:,j,i) = Diag(:,j,i) 
     &                        - Diag(:,j,k) * Diag(:,k,i)
                enddo
              endif
              Diag(:,j,i) = Diag(:,i,i) * Diag(:,j,i)
            enddo
          endif
        enddo
c
c  INACCURATE NOW     !      flops = flops + 110*nshg
c
          return
        endif  ! end LU_fact
c
c.... perform forward reduction
c
        if (code .eq. 'forward ') then
c
        do i = 2, nvar
           do j = 1, i-1
              r(:,i) = r(:,i) - Diag(:,i,j) * r(:,j)
           enddo
        enddo
c
c.... flop count
c
c  INACCURATE      !      flops = flops + 25*nshg
c
          return
        endif ! end forward
c
c.... perform backward substitution
c
        if (code .eq. 'backward') then
c
        r(:,nvar) = Diag(:,nvar,nvar) * r(:,nvar)
        
        do i = nvar-1, 1, -1
           do j = i+1, nvar
              r(:,i) = r(:,i) - r(:,j) * Diag(:,i,j)
           enddo
           r(:,i) = Diag(:,i,i) * r(:,i)
        enddo
c
c.... flop count
c
!      flops = flops + 25*nshg

          return
        endif ! end backward
c
c.... perform product U.r
c
        if (code .eq. 'product ') then
c
        do i = 1, nvar
           r(:,i) = r(:,i) / Diag(:,i,i)
           if (i .lt. nvar) then
              do j = i+1, nvar
                 r(:,i) = r(:,i) + r(:,j) * Diag(:,i,j)
              enddo
           endif
        enddo
c
c.... flop count
c
!      flops = flops + 40*nshg

          return
        endif ! end product
c
        call error ('i3LU    ', code, 0)
c
c.... return
c
c
111     format(5(e14.7,2x))
        return
        end
