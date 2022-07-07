	subroutine fillsparseI(	iens, xKebe,	lhsK,
     &                               xGoC,      lhsP,
     1			             row,	col)
c
c
c
	include "common.h"
	real*8	xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
	integer	ien(npro,nshl),	col(nshg+1), row(nshg*nnz)
	real*8	lhsK(9,nnz_tot),	lhsP(4,nnz_tot)
c
	integer	aa,	b,	c,	e,	i,	k,	n
c
	integer sparseloc

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c       
c.... Accumulate the lhs
c
	do e = 1, npro ! loop over the elements
	    do aa = 1, nshl ! loop over the local equation numbers
		i = ien(e,aa) ! finds the global equation number or
			      ! block-row of our matrix
		c = col(i)    ! starting point to look for the matching column
		n = col(i+1) - c  !length of the list of entries in rowp
		do b = 1, nshl ! local variable number tangent respect
			       ! to
c function that searches row until it finds the match that gives the
c		   global equation number

		    k = sparseloc( row(c), n, ien(e,b) ) + c-1
c
c                                             *         *
c                   dimension egmass(npro,ndof,nenl,ndof,nenl)
c
c compressible      lhsT(1:5,1:5,k)=lhsT(1:5,1:5,k)+egmass(e,1:5,aa,1:5,b)
c
		    lhsK(1,k) = lhsK(1,k) + xKebe(e,1,aa,b)
		    lhsK(2,k) = lhsK(2,k) + xKebe(e,2,aa,b)
		    lhsK(3,k) = lhsK(3,k) + xKebe(e,3,aa,b)
		    lhsK(4,k) = lhsK(4,k) + xKebe(e,4,aa,b)
		    lhsK(5,k) = lhsK(5,k) + xKebe(e,5,aa,b)
		    lhsK(6,k) = lhsK(6,k) + xKebe(e,6,aa,b)
		    lhsK(7,k) = lhsK(7,k) + xKebe(e,7,aa,b)
		    lhsK(8,k) = lhsK(8,k) + xKebe(e,8,aa,b)
		    lhsK(9,k) = lhsK(9,k) + xKebe(e,9,aa,b)
c
		    lhsP(1,k) = lhsP(1,k) + xGoC(e,1,aa,b)
		    lhsP(2,k) = lhsP(2,k) + xGoC(e,2,aa,b)
		    lhsP(3,k) = lhsP(3,k) + xGoC(e,3,aa,b)
		    lhsP(4,k) = lhsP(4,k) + xGoC(e,4,aa,b)
		enddo
	    enddo
	enddo
c
c.... end
c
	return
	end


	subroutine fillsparseC(	iens, EGmass,	lhsK,
     1			             row,	col)
c
c-----------------------------------------------------------
c This routine fills up the spasely stored LHS mass matrix
c
c Nahid Razmra, Spring 2000. (Sparse Matrix)
c-----------------------------------------------------------
c
c

	include "common.h"

 	real*8	EGmass(npro,nedof,nedof)
        integer	ien(npro,nshl),	col(nshg+1), row(nnz*nshg)
        real*8 lhsK(nflow*nflow,nnz_tot)

c
        integer aa, b, c, e, i, k, n, n1
        integer f, g, r, s, t
c
	integer sparseloc

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c
c.... Accumulate the lhsK
c
      do e = 1, npro
        do aa = 1, nshl  !loop over matrix block column index
          i = ien(e,aa)  !get mapping from aath node of element e to 
                         ! the ith row of lhsk

          c = col(i)     !Get the mapping from the ith row of lhsk to
                         ! to the corresponding location in row
          n = col(i+1) - c   !number of nonzero blocks in the row
          r = (aa-1)*nflow   !starting index of the ath node in EGmass
          do b = 1, nshl     
            s = (b-1)*nflow  !starting index of the bth node's 
                             !contribution to node aa. 
            k = sparseloc( row(c), n, ien(e,b) ) + c-1 
                !Find the index of row which corresponds to node b of 
                !element e. This is where contributions from EGmass 
                !will actually get stored. 
            
            do g = 1, nflow    !loop over columns in a block
               t = (g-1)*nflow 
               do f = 1, nflow !loop over rows in a block
                 lhsK(t+f,k) = lhsK(t+f,k) + EGmass(e,r+f,s+g)
               enddo
            enddo
          enddo
        enddo
      enddo

      return
      end
 
      subroutine fillsparseC_BC(    iens, EGmass,
     &                              lhsk, row,    col)
!
!-----------------------------------------------------------
! This routine adds contributions from heat flux BCs to the 
! spasely stored LHS mass matrix. This routine is modified 
! from fillsparseC. 
!
! Nicholas Mati, Summer 2014. (Sparse Matrix)
!-----------------------------------------------------------
!
      include "common.h"

      real*8     EGmass(npro, nshl, nshl)  !only contains term (5,5)
      integer    ien(npro,nshl),    col(nshg+1), row(nnz*nshg)
      real*8     lhsK(nflow*nflow,nnz_tot)

      integer aa, b, c, e, i, k, n, n1
      integer f, g, r, s, t
      
      integer sparseloc
      integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
      ien=abs(iens)

      !Accumulate the lhsK
      do e = 1, npro
        do aa = 1, nshl  !loop over matrix block column index
          i = ien(e,aa)  !get mapping from aath node of element e to 
                         ! the ith row of lhsk

          c = col(i)     !Get the mapping from the ith row of lhsk to
                         ! to the corresponding location in row
          n = col(i+1) - c   !number of nonzero blocks in the row
!          r = (aa-1)*nflow   !starting index of the ath node in EGmass
          do b = 1, nshl     
!            s = (b-1)*nflow  !starting index of the bth node's 
                             !contribution to node aa. 
            k = sparseloc( row(c), n, ien(e,b) ) + c-1 
                !Find the index of row which corresponds to node b of 
                !element e. This is where contributions from EGmass 
                !will actually get stored. 
            
            lhsk(25, k) = lhsk(25, k) + EGmass(e, aa, b)
!            do g = 1, nflow    !loop over columns in a block
!               t = (g-1)*nflow 
!               do f = 1, nflow !loop over rows in a block
!                 lhsK(t+f,k) = lhsK(t+f,k) + EGmass(e,r+f,s+g)
!               enddo
!            enddo
          enddo !loop over node b
        enddo !loop over node a
      enddo !loop over elements

      !end
      return
      end


      subroutine fillsparseSclr(   iens,      xSebe,    lhsS,
     1                         row,    col)
      
      include "common.h"
      real*8    xSebe(npro,nshl,nshl)
      integer    ien(npro,nshl),    col(nshg+1), row(nshg*nnz)
      real*8    lhsS(nnz_tot)    

      integer    aa,    b,    c,    e,    i,    k,    n

      integer sparseloc

      integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
      ien=abs(iens)
c
c.... Accumulate the lhs
c
      do e = 1, npro
        do aa = 1, nshl
        i = ien(e,aa)
        c = col(i)
        n = col(i+1) - c
        do b = 1, nshl
            k = sparseloc( row(c), n, ien(e,b) ) + c-1
c
            lhsS(k) = lhsS(k) + xSebe(e,aa,b)
        enddo
        enddo
      enddo
c
c.... end
c
      return
      end

      subroutine fillsparse_if
     &  ( lhsk,
     &    ienif0,ienif1,
     &    col,row,
     &    egmass,
     &    npro,
     &    nedof0, nedof1,
     &    nflow,nshg,nnz,nnz_tot)
c
        use workfc_m
c
        implicit none
c
        real*8, intent(inout) :: lhsK(nflow*nflow,nnz_tot)
        integer, dimension(:,:), pointer, intent(in) :: ienif0,ienif1
        integer, intent(in) :: col(nshg+1), row(nnz*nshg)
        real*8, intent(in) :: egmass(npro,nedof0,nedof1)
        integer, intent(in) :: nflow,nshg,nnz,nnz_tot, npro, nedof0, nedof1
c
        integer :: e,a,b,i,c,n,r,s,k,g,t,f
        integer :: nshl0,nshl1
        integer :: sparseloc
c
        nshl0 = size(ienif0,2)
        nshl1 = size(ienif1,2)
c
        do e = 1,npro
          do a = 1,nshl0
            i = ienif0(e,a)
            c = col(i)
            n = col(i+1) - c
            r = (a-1)*nflow
            do b = 1,nshl1
              s = (b-1)*nflow
              k = sparseloc( row(c), n, ienif1(e,b) ) + c-1
              do g = 1,nflow
                t = (g-1)*nflow
                do f = 1,nflow
c
                  lhsK(t+f,k) = lhsK(t+f,k) + EGmass(e,r+f,s+g)
c
                enddo
              enddo
            enddo
          enddo
        enddo        
c
      end subroutine fillsparse_if
c
        subroutine fillsparseElas(   iens,   Estiff,	lhsK,
     &			             row,    col)
c
c-----------------------------------------------------------
c
c This routine fills up the spasely stored LHS mesh-elastic
c stiffness matrix. 
c
c-----------------------------------------------------------
c
	include "common.h"
c
 	real*8	Estiff(npro,ndofelas,ndofelas)
        integer	ien(npro,nshl),	col(nshg+1), row(nnz*nshg)
        real*8 lhsK(nelas*nelas,nnz_tot)
c
        integer aa, b, c, e, i, k, n, n1
        integer f, g, r, s, t
c
	integer sparseloc
c
	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c
c.... Accumulate the lhsK
c
	do e = 1, npro
	    do aa = 1, nshl
		i = ien(e,aa)
		c = col(i)
		n = col(i+1) - c
		r = (aa-1)*nelas
		do b = 1, nshl
		    s = (b-1)*nelas
                    k = sparseloc( row(c), n, ien(e,b) ) + c-1
c
		    do g = 1, nelas
		       t = (g-1)*nelas
		       do f = 1, nelas
c
			  lhsK(t+f,k) = lhsK(t+f,k) + Estiff(e,r+f,s+g)
c
                       enddo
                   enddo
		enddo
	    enddo
	enddo
c
c.... end
c
	return
	end subroutine fillsparseElas
c
      integer    function sparseloc( list, n, target )

c-----------------------------------------------------------
c This function finds the location of the non-zero elements
c of the LHS matrix in the sparsely stored matrix 
c lhsK(nflow*nflow,nnz*numnp)
c
c Nahid Razmara, Spring 2000.     (Sparse Matrix)
c-----------------------------------------------------------

      integer    list(n),    n,    target
      integer    rowvl,    rowvh,    rowv

c
c.... Initialize
c
      rowvl = 1
      rowvh = n + 1
c
c.... do a binary search
c
100    if ( rowvh-rowvl .gt. 1 ) then
        rowv = ( rowvh + rowvl ) / 2
        if ( list(rowv) .gt. target ) then
        rowvh = rowv
        else
        rowvl = rowv
        endif
        goto 100
      endif
c
c.... return
c
      sparseloc = rowvl
c
      return
      end

