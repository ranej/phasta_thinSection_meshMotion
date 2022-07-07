      module e3_solid_m
c
        use e3_param_m
        use solid_data_m
c
        implicit none
c
        integer :: iblk_solid
c
        real*8, dimension(:,:), pointer :: d
        real*8, dimension(:), pointer :: det_d, det_baf
        real*8, dimension(:), pointer :: bulkMod, shearMod, Ja_def
c
        contains
c
        subroutine e3_malloc_solid
c
          call e3_malloc
c
          allocate(d(npro,b_size))
          allocate(det_d(npro),det_baf(npro))
          allocate(Ja_def(npro),bulkMod(npro),shearMod(npro))
c
        end subroutine e3_malloc_solid
c
        subroutine e3_mfree_solid
c
          call e3_mfree
c
          deallocate(d)
          deallocate(det_d,det_baf)
          deallocate(Ja_def,bulkMod,shearMod)
c
        end subroutine e3_mfree_solid
c
      end module e3_solid_m
c
      module e3_solid_func_m
c
        use e3_solid_m
        use number_def_m
        use global_const_m
        use intpt_m
        use inpdat_m
        use conpar_m
        use timdat_m
c
        implicit none
c
        real*8, dimension(:,:), pointer :: dudx, dudy, dudz
        real*8, dimension(:,:,:), pointer :: AS
c
      contains
c
      subroutine calc_solid
c
        implicit none
c
        allocate(AS(npro,6,6))
c
        call calc_as_matrix
        d(:,:) = almBi * b(iblk_solid)%p(:,intp,:)
     &+      alfBi * Delt(1) * (almBi - gamBi) 
     &*      b_dot(iblk_solid)%p(:,intp,:)
        call setB_af
        call get_det(b_af(iblk_solid)%p(:,intp,:),det_baf,npro)
        Ja_def= (det_baf)**0.5
        call get_det(d,det_d,npro)
c
        deallocate(AS)
c
      end subroutine calc_solid
c
      subroutine calc_solid_bdy
c
        implicit none
c
        allocate(AS(npro,6,6))
c
        call calc_as_matrix
        d(:,:) = almBi * bdy_b(iblk_solid)%p(:,intp,:)
     &+      alfBi * Delt(1) * (almBi - gamBi) 
     &*      bdy_b_dot(iblk_solid)%p(:,intp,:)
        call setB_af_bdy
        call get_det(bdy_b_af(iblk_solid)%p(:,intp,:),det_baf,npro)
        Ja_def= (det_baf)**0.5
        call get_det(d,det_d,npro)
c
        deallocate(AS)
c
      end subroutine calc_solid_bdy
c
      subroutine calc_as_matrix
c
        implicit none
c
c initialize the AS matrix
         AS = zero !add
c
         AS(:,1,1) = 2 * dudx(:,1)
         AS(:,1,5) = 2 * dudz(:,1)
         AS(:,1,6) = 2 * dudy(:,1)
c
         AS(:,2,2) = 2 * dudy(:,2)
         AS(:,2,4) = 2 * dudz(:,2)
         AS(:,2,6) = 2 * dudx(:,2)
c
         AS(:,3,3) = 2 * dudz(:,3)
         AS(:,3,4) = 2 * dudy(:,3)
         AS(:,3,5) = 2 * dudx(:,3)
c
         AS(:,4,2) = dudy(:,3)
         AS(:,4,3) = dudz(:,2)
         AS(:,4,4) = dudy(:,2) + dudz(:,3)
         AS(:,4,5) = dudx(:,2)
         AS(:,4,6) = dudx(:,3)
c
         AS(:,5,1) = dudx(:,3)
         AS(:,5,3) = dudz(:,1)
         AS(:,5,4) = dudy(:,1)
         AS(:,5,5) = dudx(:,1) + dudz(:,3)
         AS(:,5,6) = dudy(:,3)
c         
         AS(:,6,1) = dudx(:,2)
         AS(:,6,2) = dudy(:,1)
         AS(:,6,4) = dudz(:,1)
         AS(:,6,5) = dudz(:,2)
         AS(:,6,6) = dudx(:,1) + dudy(:,2)
c
c
      end subroutine calc_as_matrix
c
      subroutine setB_af
c
c... calculate the left Cauchy-green tensor at time step n+af
c
         implicit none 
c
         integer,parameter :: nsize = 6 
c
         real*8, dimension(6,6) :: ident
         real*8, dimension(6,6) :: temp_matrix
         integer :: i
c
         ident = zero
         do i = 1, nsize
           ident(i,i) = one
         enddo
c
         do i = 1, npro
          temp_matrix(:,:) = almBi * ident(:,:) + gamBi * Delt(1)
     &                     *alfBi * AS(i,:,:) !check here
          b_af(iblk_solid)%p(i,intp,:) = matmul(temp_matrix(:,:) , d(i,:))
         enddo
c
       end subroutine setB_af 
c
      subroutine setB_af_bdy
c
c... calculate the left Cauchy-green tensor at time step n+af
c
         implicit none 
c
         integer,parameter :: nsize = 6 
c
         real*8, dimension(6,6) :: ident
         real*8, dimension(6,6) :: temp_matrix
         integer :: i
c
         ident = zero
         do i = 1, nsize
           ident(i,i) = one
         enddo
c
         do i = 1, npro
          temp_matrix(:,:) = almBi * ident(:,:) + gamBi * Delt(1)
     &                     *alfBi * AS(i,:,:) !check here
          bdy_b_af(iblk_solid)%p(i,intp,:) = matmul(temp_matrix(:,:) , d(i,:))
         enddo
c
       end subroutine setB_af_bdy
c
       subroutine get_det(matr,det_matr,npro)
c
        implicit none 
c
        real*8, dimension(npro,6),intent(in) :: matr
        real*8, dimension(npro) :: det_matr
        integer, intent(in) :: npro
c
        det_matr(:) = matr(:,1) * matr(:,2) * matr(:,3)
     &+             matr(:,6) * matr(:,4) * matr(:,5)
     &+             matr(:,5) * matr(:,6) * matr(:,4) 
     &-             matr(:,1) * matr(:,4) * matr(:,4)
     &-             matr(:,5) * matr(:,2) * matr(:,5)
     &-             matr(:,6) * matr(:,6) * matr(:,3)
c
      end subroutine get_det
c
      end module e3_solid_func_m
