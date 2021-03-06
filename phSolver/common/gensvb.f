        subroutine gensvb (ienb,   ibcb,    bcb, im2gb, rotBandIndex)
c
          use mattype_m
c----------------------------------------------------------------------
c
c  This routine saves the boundary element block.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        integer, dimension(npro,nshl),        intent(out) :: ienb
        integer, dimension(npro,ndibcb),      intent(out) :: ibcb
        real*8,  dimension(npro,nshlb,ndbcb), intent(out) :: bcb
        integer, dimension(npro,3),           intent(out) :: im2gb
        integer, dimension(npro),             intent(out) :: rotBandIndex

c
c.... generate the boundary element mapping
c
        do i = 1, nshl
          ienb(1:npro,i) = ientmp(1:npro,i)
        enddo
c
c.... save the boundary element data
c
        iBCB   = iBCBtmp
        im2gb  = im2gbtmp
        rotBandIndex = rotBandIndextmp
        do i=1, npro
           write(*,*) "rotBandIdex gensvb:" , rotBandIndex(i)
        end do

        do i = 1, nenbl ! This is NOT NSHLB as we are just copying the
                        ! piecewise constant data given by NSpre and
                        ! higher order coefficients must be zero
           do j = 1, ndBCB
              BCB(:,i,j)   = BCBtmp(:,j)
           end do
        end do
        do i = nenbl+1, nshlb
           do j = 1, ndBCB
              BCB(:,i,j)   = zero
           end do
        end do
c
c.... return
c
        return
        end
