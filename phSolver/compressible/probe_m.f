      module probe_def_m
        implicit none
        type probe_t
          real*8 :: mass, momentum(3), energy
        end type probe_t
        integer, parameter :: nprobe = 1
     &,                       iprobe_domain = 1
      end module probe_def_m
c
      module probe_m
c
c----------------------------------------
c    aims to monitor the conservation laws over
c    the entire domain and print out the values
c----------------------------------------
        use probe_def_m
        use global_const_m
        use number_def_m
        use workfc_m
        use block_m
        use e3_func_m
        use local_m
        implicit none
c
        type(probe_t) :: probe(nprobe)
c
      contains
c
      subroutine probe_conservation(y,x,shp_table,shgl_table)
        real*8, intent(in) :: y(nshg,ndof), x(numnp,nsd)
        real*8, intent(in) :: shp_table(MAXTOP,maxsh,MAXQPT), shgl_table(MAXTOP,nsd,maxsh,MAXQPT)
        integer :: iblk
        integer :: sgn(npro,nshl)
        real*8, pointer :: s(:,:)
c
        call init_probe(probe(1))
c
        set_block_ptr => set_interior_block
        do iblk = 1,nelblk
c
          call set_interior_block(iblk)
          call set_block_ptr(iblk)
c
          call e3_malloc_ptr
          allocate(yl(npro,nshl,ndof),xl(npro,nenl,nsd))
          allocate(s(npro,nflow))
c
          if (ipord .gt. 1) then
            call getsgn(ien,sgn)
          endif
c
          call localy(y,yl,ien,ndofl,'gather  ',nshg,nshl,npro,ipord)
          call localx(x,xl,ien,nsd,'gather  ',nshg,nshl,npro)

          call e3_int_conserv(s,shp_table(lcsyst,1:nshl,:),shgl_table(lcsyst,1:nsd,1:nshl,:),sgn)
c
c...apply any material condtion here..
c
          call set_probe(probe(1),s)
c
          call e3_mfree_ptr
          deallocate(yl,xl)
          deallocate(s)
c
        enddo
c
        call get_probe_global
c
      end subroutine probe_conservation
c
      subroutine init_probe(this)
        type(probe_t), intent(inout) :: this
        this%mass     = zero
        this%momentum = zero
        this%energy   = zero
      end subroutine init_probe
c
      subroutine set_probe(this,dui)
        type(probe_t), intent(inout) :: this
        real*8, dimension(npro,nflow), intent(in) :: dui
        this%mass        = this%mass        + sum(dui(:,1))
        this%momentum(1) = this%momentum(1) + sum(dui(:,2))
        this%momentum(2) = this%momentum(2) + sum(dui(:,3))
        this%momentum(3) = this%momentum(3) + sum(dui(:,4))
        this%energy      = this%energy      + sum(dui(:,5))
      end subroutine set_probe
c
      subroutine get_probe_global
        use timdat_m
        use mio_m
        implicit none
        include "mpif.h"
        real*8, dimension(nprobe*nflow) :: myprobe,globalprobe
        integer :: i, ierr
        do i = 1,nprobe
          myprobe(1+(i-1)*nflow) = probe(i)%mass
          myprobe(2+(i-1)*nflow) = probe(i)%momentum(1)
          myprobe(3+(i-1)*nflow) = probe(i)%momentum(2)
          myprobe(4+(i-1)*nflow) = probe(i)%momentum(3)
          myprobe(5+(i-1)*nflow) = probe(i)%energy
        enddo
        globalprobe = zero
        call mpi_allreduce(myprobe,globalprobe,nprobe*nflow,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
        if (myrank == 0) then
          write(*,200) lstep,globalprobe(1:5)
          write(iconserv,100) lstep,globalprobe(1:5)
        endif
100   format(1x,i6,5e24.16)
200   format(1p,'Conservation: ',i6,5(2x,e10.3))
      end subroutine get_probe_global
c
      end module probe_m
