      module ifbc_def_m
c
        use global_const_m
        use conpar_m
        use if_velocity_m
c
        implicit none
c
        integer, parameter :: nifbc = 1
        integer, parameter :: ivapor_frac = 1
        integer, parameter :: iwt = 2
c
        real*8, pointer :: if_normal(:,:)
        real*8, pointer :: ifbc(:,:)
c
      end module ifbc_def_m
c
      module ifbc_m
c
        use ifbc_def_m
        use workfc_m
        use number_def_m
        use pointer_data
c
        implicit none
c
        contains
c
        subroutine ifbc_malloc
          allocate(if_normal(nshg,nsd))
          allocate(ifbc(nshg,nifbc+1))
          ifbc(:,ivapor_frac) = zero
          ifbc(:,iwt) = one
        end subroutine ifbc_malloc
c
        subroutine ifbc_mfree
          deallocate (if_normal)
          deallocate (ifbc)
        end subroutine ifbc_mfree
c
        subroutine ifbc_set
     &  (
     &    y
     &,   BC
     &,   iBC
     &,   ilwork
     &,   nlwork
     &  )
c
          use elmpar_m
          use blkdat_m
          use shpdat_m
          use genpar_m
          use matdat_def_m
c
          include "mpif.h"
c
          real*8, intent(inout) :: y(nshg,ndof)
          real*8,  intent(inout) ::  BC(nshg,ndofbc)
          integer, intent(inout) :: iBC(nshg)
          integer, intent(in)    :: ilwork(nlwork)
          integer, intent(in)    :: nlwork
c
          integer :: i0,i1,n,iel,iblk,materif0,materif1, npro, ierr, isclr
          integer, pointer :: ienif0(:,:),ienif1(:,:)
          logical :: is_ideal_gas,is_ideal_gas_mixture
c
          if (numpe > 1) then
            call commu (ifbc(:,ivapor_frac), ilwork, 1, 'in ')
            call commu (ifbc(:,iwt), ilwork, 1, 'in ')
            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
            call commu (ifbc(:,ivapor_frac), ilwork, 1, 'out')
            call commu (ifbc(:,iwt), ilwork, 1, 'out ')
            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          endif
c
          where(ifbc(:,ivapor_frac) .ne. zero)
            BC(:,6+isclr) = ifbc(:,ivapor_frac)/ifbc(:,iwt)
          endwhere
c
          return            
c
          do iblk = 1,nelblif
c
            npro  = lcblkif(1,iblk+1) - lcblkif(1,iblk)
c
            ienif0 => mienif0(iblk)%p
            ienif1 => mienif1(iblk)%p
c
            do iel = 1,npro
c              do n = 1,nshl0
              do n = 1,3  ! only triangles on the surface
                i0 = ienif0(iel,n)
                i1 = ienif1(iel,n)
c
c... set vapor fraction
c
      isclr = 1
                BC(i0,6+isclr) = ifbc(i0,ivapor_frac)/ifbc(i0,iwt)
C                BC(i0,6+isclr) = zero
                BC(i1,6+isclr) = one
                y(i0,5+isclr) = BC(i0,6+isclr)
                y(i1,5+isclr) = BC(i1,6+isclr)
c      write(*,*) 'i0, i1: ',i0,i1,iBC(i0),iBC(i1),btest(ibc(i0),6),btest(ibc(i1),6)
C      write(*,*) '[',myrank,'] i0: ',i0,ifbc(i0,ivapor_frac),ifbc(i0,iwt),BC(i0,6+isclr),btest(ibc(i0),6)
c      write(*,*) 'i1: ',i1,ifbc(i1,ivapor_frac),ifbc(i1,iwt),BC(i1,6+isclr),btest(ibc(i1),6)
c
              enddo
            enddo
          enddo
c
        end subroutine ifbc_set
c
      end module ifbc_m
