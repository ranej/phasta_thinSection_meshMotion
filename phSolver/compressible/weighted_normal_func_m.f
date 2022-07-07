      module weighted_normal_func_m
c
c------------------------------------------------------------------------------
c  calculating the weighted normal at quadrature points
c------------------------------------------------------------------------------
        implicit none
c
c
      contains
        subroutine get_weighted_normal(nv, w_normal, shp, nshl)
c-------------------------------------------------------------------------------
          use e3if_func_m, only:sum_qpt
          use propar_m, only: npro
          use global_const_m, only: nsd
          implicit none
c
          real*8, dimension(npro,nsd),intent(inout) :: nv
          real*8, dimension(npro,nshl,nsd),intent(in) :: w_normal
          real*8, dimension(npro,nshl), intent(in) :: shp
          integer, intent(in) :: nshl
c          
          real*8, dimension(npro,nsd) :: w_n_qpt
          real*8, dimension(npro) :: temp_len
          integer :: iel

c          
c... get the weighted normal at qudrature points
          do iel = 1,npro
            w_n_qpt(iel,1) = sum_qpt(nshl,w_normal(iel,:,1),shp(iel,:))
            w_n_qpt(iel,2) = sum_qpt(nshl,w_normal(iel,:,2),shp(iel,:))
            w_n_qpt(iel,3) = sum_qpt(nshl,w_normal(iel,:,3),shp(iel,:))

c... get the length of weighted normal at qudrature points
            temp_len(iel) = sqrt(w_n_qpt(iel,1)*w_n_qpt(iel,1)
     &                    + w_n_qpt(iel,2)*w_n_qpt(iel,2)
     &                    + w_n_qpt(iel,3)*w_n_qpt(iel,3))
c... normalization of the weighted normal
            nv(iel,1)  = w_n_qpt(iel,1) / temp_len(iel)
            nv(iel,2)  = w_n_qpt(iel,2) / temp_len(iel)
            nv(iel,3)  = w_n_qpt(iel,3) / temp_len(iel)
          enddo
c        
        end subroutine get_weighted_normal
      end module weighted_normal_func_m
