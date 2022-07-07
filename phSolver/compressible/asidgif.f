      subroutine asidgif
     & (
     &   res,
     &   y,        x,       umesh,
     &   shpif0,   shpif1,  shgif0,  shgif1,
     &   qwtif, qwtif0,   qwtif1,
     &   ienif0,   ienif1,
     &   BDiag)  
c
c----------------------------------------
c    aims to localize, compute and assemble
c    the contribution of DG interface
c----------------------------------------
          use hierarchic_m
          use local_m
          use e3if_m
          use e3if_geom_m
          use if_global_m
          use conpar_m
          use weighted_normal_data_m, only:w_normal_l0, w_normal_l1, w_normal_global !for weighted normal on the interface
          use genpar_m, only: iprec
          use dgifinp_m, only: i_w_normal
c
          implicit none
c
          real*8, dimension(nshg,nflow), intent(inout) :: res
          real*8, dimension(nshg,ndof),  intent(in)    :: y
          real*8, dimension(nshg,nsd),   intent(in)    :: x
          real*8, dimension(nshg, nsd), intent(inout) :: umesh
          real*8, dimension(nshl0,nqpt),intent(in)   :: shpif0
          real*8, dimension(nshl1,nqpt),intent(in)   :: shpif1
          real*8, dimension(nsd,nshl0,nqpt), intent(in)  :: shgif0
          real*8, dimension(nsd,nshl1,nqpt), intent(in)  :: shgif1
          real*8, dimension(nqpt), intent(in) :: qwtif, qwtif0, qwtif1
          real*8, dimension(nshg,nflow,nflow), intent(inout) :: BDiag
          integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
          integer :: i0,i1,iel,n,i,npro_,imin(5),imax(5)
c
          integer :: j, k, j0, loc  
c
#define debug 0
#define debugifbc 1
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
          if (ipord .gt. 1) then
           call getsgn(ienif0,sgn0,nshl0,nenl0)
           call getsgn(ienif1,sgn1,nshl1,nenl1)
        endif
c
c... localize
c
        call localy(y, ycl0, ienif0, ndof, 'gather  ', nshg, nshl0, npro, ipord)
        call localy(y, ycl1, ienif1, ndof, 'gather  ', nshg, nshl1, npro, ipord)
c
        call localx(x, xl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
        call localx(x, xl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
c
        call localx(umesh, umeshl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
        call localx(umesh, umeshl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
c
        call localx(if_normal, if_normal_l0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
        call localx(if_normal, if_normal_l1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
c
        if (associated(if_kappa)) then
          call localx(if_kappa, if_kappa_l0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
          call localx(if_kappa, if_kappa_l1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
        endif
c... localize the weighted normal
        if (i_w_normal .eq. 1) then
          call localx(w_normal_global, w_normal_l0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
          call localx(w_normal_global, w_normal_l1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
        endif
c
        call e3if(shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
#if debug==1
      npro_ = npro
      imin = minloc(res,1)
      imax = maxloc(res,1)
c      i = 198
c      i = imax(1)
c      i = imin(1)
      i = 97
      do iel = 1,npro_
        do n = 1,nshl0
          if (ienif0(iel,n) .ne. i) cycle
          write(*,10) 'x0: ',iel,n,ienif0(iel,n),xl0(iel,n,:)
        enddo
      enddo
      do iel = 1,npro_
        do n = 1,nshl1
          if (ienif1(iel,n) .ne. i) cycle
          write(*,10) 'x1: ',iel,n,ienif1(iel,n),xl1(iel,n,:)
        enddo
      enddo
      do iel = 1,npro_
        do n = 1,nshl0
          if (ienif0(iel,n) .ne. i) cycle
          write(*,20) 'res before: ',ienif0(iel,n),res(ienif0(iel,n),:)
        enddo
      enddo
      do iel = 1,npro_
        do n = 1,nshl1
          if (ienif1(iel,n) .ne. i) cycle
          write(*,20) 'res before: ',ienif1(iel,n),res(ienif1(iel,n),:)
        enddo
      enddo
#endif
c
c.... assemble the local residual arrays
c
        call local (res, rl0, ienif0, nflow, 'scatter ', nshg,nshl0,npro,ipord)
        call local (res, rl1, ienif1, nflow, 'scatter ', nshg,nshl1,npro,ipord)
c
#if debug==1
      do iel = 1,npro_
        do n = 1,nshl0
          if (ienif0(iel,n) .ne. i) cycle
          write(*,10) 'rl0: ',iel,n,ienif0(iel,n),rl0(iel,n,:)
        enddo
      enddo
      do iel = 1,npro_
        do n = 1,nshl1
          if (ienif1(iel,n) .ne. i) cycle
          write(*,10) 'rl1: ',iel,n,ienif1(iel,n),rl1(iel,n,:)
        enddo
      enddo
      do iel = 1,npro_
        do n = 1,nshl0
          if (ienif0(iel,n) .ne. i) cycle
          write(*,20) 'res after: ',ienif0(iel,n),res(ienif0(iel,n),:)
        enddo
      enddo
      do iel = 1,npro_
        do n = 1,nshl1
          if (ienif1(iel,n) .ne. i) cycle
          write(*,20) 'res after: ',ienif1(iel,n),res(ienif1(iel,n),:)
        enddo
      enddo
#endif
c... get the element preconditioner
c... extract and assemble the Block-Diagonal (same as ASIGMR)
        BDiagl_00 = zero
        BDiagl_11 = zero
c        
        if (iprec .ne. 0) then 
          do i = 1, nshl0
              do j = 1, nflow
                 i0 = (i - 1) * nflow + j
                 do k = 1, nflow
                    j0 = (i - 1) * nflow + k
                    loc = (k-1)*nflow + j ! the index number when converting 5 by 5 into 25 by 1 vector
                    BDiagl_00(:,i,loc) = egmass00(:,i0,j0)
                 enddo
              enddo
          enddo
c
           do i = 1, nshl1
              do j = 1, nflow
                 i0 = (i - 1) * nflow + j
                 do k = 1, nflow
                    j0 = (i - 1) * nflow + k
                    loc = (k-1)*nflow + j
                    BDiagl_11(:,i,loc) = egmass11(:,i0,j0)
                 enddo
              enddo
           enddo
c           
           call local (BDiag,  BDiagl_00, ienif0, nflow*nflow, 'scatter ', nshg, nshl0,npro,ipord)
           call local (BDiag,  BDiagl_11, ienif1, nflow*nflow, 'scatter ', nshg, nshl1,npro,ipord)
        endif
c
c
        call local (sum_vi_area, sum_vi_area_l0, ienif0, nsd+1, 'scatter ', nshg, nshl0,npro,ipord)
        call local (sum_vi_area, sum_vi_area_l1, ienif1, nsd+1, 'scatter ', nshg, nshl1,npro,ipord)
c
        call local (ifbc, ifbc_l0, ienif0, nifbc+1, 'scatter ', nshg, nshl0, npro, ipord)
        call local (ifbc, ifbc_l1, ienif1, nifbc+1, 'scatter ', nshg, nshl1, npro, ipord)
c
10    format(a12,3i6,5e24.16)
20    format(a12,1i6,5e24.16)
30    format(a12,1i6,e24.16)
c
      end subroutine asidgif
