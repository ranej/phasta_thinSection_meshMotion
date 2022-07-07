      module if_velocity_m
c
c----------------------------------------
c    aims to calculate the interface
c    velocity at the global level
c----------------------------------------
        use workfc_m
        use number_def_m
        use pointer_data
        use blkdat_m
        use interfaceflag
c	use bc_on_vi_m
c
        implicit none
c
        real*8, pointer :: sum_vi_area(:,:)    ! interface velocity weighted by interface area
c
      contains
c
      subroutine init_sum_vi_area(nshg,nsd)
        integer, intent(in) :: nshg,nsd
        if (associated(sum_vi_area)) 
     &    deallocate (sum_vi_area)
        allocate (sum_vi_area(nshg,nsd+1))
        sum_vi_area(:,1:nsd) = zero
        sum_vi_area(:,1+nsd) = one
      end subroutine init_sum_vi_area
c
      subroutine destruct_sum_vi_area
        if (associated(sum_vi_area))
     &    deallocate(sum_vi_area)
      end subroutine destruct_sum_vi_area
c
      subroutine set_if_velocity 
     & (
     &  BC, iBC, umesh, disp, x, dt, ilwork,
     &  nshg, ndofBC, nsd, nelblif, nlwork, ndof
     & )
c
        include "mpif.h"
c
        real*8,  intent(inout) ::  BC(nshg,ndofBC)
        integer, intent(inout) :: iBC(nshg)
        real*8,  dimension(nshg,nsd), intent(inout)    :: umesh, disp
        real*8,  dimension(nshg,nsd), intent(in)    :: x
        real*8, intent(in) :: dt
        integer, intent(in)    :: ilwork(nlwork)
        integer, intent(in) :: nshg, ndofBC, nsd, nelblif, nlwork, ndof
c
        integer :: iblk, iel, npro,inode, i0, i1, n, ierr
        integer, pointer :: ienif0(:,:), ienif1(:,:)
	real*8, dimension(nshg,3) :: actual_vi

        if (numpe > 1) then
          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'in ')
          call commu (sum_vi_area(:,4), ilwork, 1, 'in ')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c
        if (numpe > 1) then
          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'out')
          call commu (sum_vi_area(:,4), ilwork, 1, 'out')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c
        actual_vi = zero
        do inode = 1, nshg
          if ( ifFlag(inode) .eq. 1 ) then
c            write(*,*) "rank",myrank,"i",inode,"x=",x(inode,:)
            actual_vi(inode,:) = sum_vi_area(inode,:) / sum_vi_area(inode,nsd+1)
          endif
        enddo
c
        call itrBCvi ( actual_vi ,iBC ,BC )
c
        do inode = 1, nshg
          if ( ifFlag(inode) .eq. 1 ) then
            umesh(inode,:) = actual_vi(inode,:)
c
c.... the following line is moved to solve mesh part
c.... since in restart case, we should update interface mesh BC
c.... before we solve the mesh
c            BC(inode,ndof+2:ndof+4) = umesh(inode,:) * dt
          endif
        enddo
c
c
100   format(a,'[',i2,'] ',i6,3f7.3,x,7e14.6)
200   format(a,'[',i2,'] ',i6,3e14.6)
      end subroutine set_if_velocity
c
      end module if_velocity_m
