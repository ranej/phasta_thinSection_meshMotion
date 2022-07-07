      module e3_func_m
c
c----------------------------------------
c    aims to calculate the conservative
c    and primitive variables for the
c    interior elements
c----------------------------------------
        use conpar_m
        use elmpar_m
        use genpar_m
        use e3_param_m
        use eqn_state_m
        use shpdat_m
        use intpt_m
        use e3metric_m
        implicit none
c
        real*8, dimension(:), pointer :: detJ
        real*8, dimension(:,:,:), pointer :: yl,xl
c
      contains
c
      subroutine e3_int_conserv(s,shp,shgl,sgn)
        implicit none
        real*8, dimension(:,:), pointer :: s
        real*8, dimension(nshl,ngauss), intent(in) :: shp
        real*8, dimension(nsd,nshl,ngauss), intent(in) :: shgl
        integer, dimension(npro,nshl), intent(in) :: sgn
c
        integer :: iflow
        real*8 :: shp_qpt(npro,nshl)
        real*8 :: dui(npro,nflow)
c        real*8 :: dxidx(npro,nsd,nsd)
        real*8, dimension(:), pointer :: WdetJ
        real*8, dimension(:,:,:), pointer :: shg_qpt, shgl_qpt
        real*8, dimension(:,:,:), pointer :: dxidx ! added
c
        allocate(WdetJ(npro))
        allocate(shg_qpt(npro,nshl,nsd))
        allocate(shgl_qpt(npro,nsd,nshl))
        allocate(dxidx(npro,nsd,nsd)) !added
c
        s = zero
c
        do intp = 1,ngauss
c
          call getshp(shp, shgl, sgn, shp_qpt, shgl_qpt)
          call calc_primitive(dui,shp_qpt)
          call getthm6_ptr
          call calc_conservative(dui)
          call e3metric(shg_qpt,dxidx,shgl_qpt,xl)
c
c... multiply by weight
c
          do iflow = 1,nflow
            s(:,iflow) = s(:,iflow) + wdetj(:)*dui(:,iflow)
          enddo
c
        enddo
c
        deallocate(WdetJ,shg_qpt,shgl_qpt)
c
      end subroutine e3_int_conserv
c
      subroutine calc_primitive(dui,shp)
        real*8, dimension(npro,nflow), intent(out) :: dui
        real*8, dimension(npro,nshl), intent(in) :: shp
        integer :: n
        dui = zero
        do n = 1, nshl
           dui(:,1) = dui(:,1) + shp(:,n) * yl(:,n,1) ! p
           dui(:,2) = dui(:,2) + shp(:,n) * yl(:,n,2) ! u1
           dui(:,3) = dui(:,3) + shp(:,n) * yl(:,n,3) ! u2
           dui(:,4) = dui(:,4) + shp(:,n) * yl(:,n,4) ! u3
           dui(:,5) = dui(:,5) + shp(:,n) * yl(:,n,5) ! T
        enddo
        pres = dui(:,1)
        T = dui(:,5)
      end subroutine calc_primitive
c
      subroutine calc_conservative(dui)
        real*8, dimension(npro,nflow), intent(inout) :: dui
        dui(:,1) = rho
        dui(:,2) = rho * dui(:,2)
        dui(:,3) = rho * dui(:,3)
        dui(:,4) = rho * dui(:,4)
        dui(:,5) = rho * (ei + rk)
      end subroutine calc_conservative
c
      end module e3_func_m
