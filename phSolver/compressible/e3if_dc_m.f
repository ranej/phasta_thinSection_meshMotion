      module e3if_dc_m
c
c------------------------------------------------------------------------------
c  calculating the discontinuity capturing (DC) operator for the interface
c  macro-elements
c------------------------------------------------------------------------------
        implicit none
c
c
      contains
        subroutine calc_projector(proj, nv)
c..............................................................................
c  calculating the projector which project a quanlity into its tangential
c  components as :
c  \bm{P} =  \bm{I} - \bm{n} \otimes \bm{n}
c..............................................................................        
          use number_def_m
          use propar_m, only: npro
          use global_const_m, only: nsd
          implicit none
c          
          real*8,  dimension(npro,nsd), intent(in) :: nv ! normal vector
          real*8,  dimension(npro,nsd,nsd), intent(out) :: proj ! projector
c
          proj = zero
c
          proj(:,1,1) = one - nv(:,1)*nv(:,1)
          proj(:,1,2) =     - nv(:,1)*nv(:,2)
          proj(:,1,3) =     - nv(:,1)*nv(:,3)
c
          proj(:,2,1) =     - nv(:,2)*nv(:,1)
          proj(:,2,2) = one - nv(:,2)*nv(:,2)
          proj(:,2,3) =     - nv(:,2)*nv(:,3)
c
          proj(:,3,1) =     - nv(:,3)*nv(:,1)
          proj(:,3,2) =     - nv(:,3)*nv(:,2)
          proj(:,3,3) = one - nv(:,3)*nv(:,3)
  
c
        end subroutine calc_projector
c
        subroutine calc_ch(ch0, ch1, f_jump)
c..............................................................................
c  calculating the c^h of the DC operator for the interface, which is analogous 
c  to the nu^h for the volumn elements. 
c  Unit of c^h ~ [Length]/ [Time]
c..............................................................................        
          use number_def_m
          use propar_m, only: npro
          use conpar_m, only: nflow
          use dgifinp_m, only: if_e_dc
          implicit none
c
          real*8, dimension(npro),intent(out) :: ch0, ch1 
          real*8, dimension(npro,nflow),intent(in) :: f_jump ! flux jump
c
          real*8, dimension(npro,nflow) :: u_ref_0, u_ref_1 ! reference conservative
                                                            ! variables
          real*8, dimension(npro,nflow) :: temp0, temp1 ! local temporary array
          integer :: iel                                                  
c
          u_ref_0(:,1) = 1.000000000000000d-2 ! hacking, rho_ref, gas phase
          u_ref_0(:,2) = 1.000000000000000d0 ! hacking, rho_ref*v_ref
          u_ref_0(:,3) = 1.000000000000000d0 ! hacking, rho_ref*v_ref
          u_ref_0(:,4) = 1.000000000000000d0 ! hacking, rho_ref*v_ref
          u_ref_0(:,5) = 2.114165517241379d3 ! hacking, rho_ref*(ei + 0.5*v_ref^2)
c
          u_ref_1(:,1) = 1.000000000000000d0 ! hacking, rho_ref, liquid phase
          u_ref_1(:,2) = 1.000000000000000d0  ! hacking, rho_ref*v_ref
          u_ref_1(:,3) = 1.000000000000000d0 ! hacking, rho_ref*v_ref
          u_ref_1(:,4) = 1.000000000000000d0 ! hacking, rho_ref*v_ref
          u_ref_1(:,5) = 4.204805000000000d5! hacking, rho_ref*(ei + 0.5*v_ref^2)
c... diag(U1_ref ... U5_ref) * flux_jump          
          temp0(:,1) = ( one/u_ref_0 (:,1) )*f_jump(:,1)
          temp0(:,2) = ( one/u_ref_0 (:,2) )*f_jump(:,2)
          temp0(:,3) = ( one/u_ref_0 (:,3) )*f_jump(:,3)
          temp0(:,4) = ( one/u_ref_0 (:,4) )*f_jump(:,4)
          temp0(:,5) = ( one/u_ref_0 (:,5) )*f_jump(:,5)
c
          temp1(:,1) = ( one/u_ref_1 (:,1) )*f_jump(:,1)
          temp1(:,2) = ( one/u_ref_1 (:,2) )*f_jump(:,2)
          temp1(:,3) = ( one/u_ref_1 (:,3) )*f_jump(:,3)
          temp1(:,4) = ( one/u_ref_1 (:,4) )*f_jump(:,4)
          temp1(:,5) = ( one/u_ref_1 (:,5) )*f_jump(:,5)
c... c^h for each phase
          do iel = 1,npro
            ch0(iel) = if_e_dc 
     &               * sqrt(dot_product(temp0(iel,:),temp0(iel,:)))
            ch1(iel) = if_e_dc 
     &               * sqrt(dot_product(temp1(iel,:),temp1(iel,:)))
          enddo
c          
        end subroutine calc_ch
c
        subroutine e3if_dc_res(ri, var, A0, ch, pt_g_p)
c..............................................................................
c  calculation of the contribution of DC operator to the integrand of residual
c  for both phases
c  G_b = c^h g^{km}_{I} N_{a,k} A0 Y_{,m} d \Gama
c..............................................................................
          use propar_m, only: npro
          use conpar_m, only: nflow
          use global_const_m, only: nsd
          use e3if_param_m, only: var_t 
          implicit none
c
          real*8, dimension(npro,nflow*(nsd+1)), intent(inout) :: ri
          real*8, dimension(npro,nflow,nflow), intent(in) :: A0
          type(var_t), dimension(npro), intent(in) :: var
          real*8, dimension(npro), intent(in) :: ch
          real*8, dimension(npro,nsd,nsd),intent(in) :: pt_g_p
c
          real*8, dimension(npro,nflow,nsd) :: A0_grad_y
          integer :: iel, k, k0
c
          do iel = 1,npro
c... A0 Y_{,m}          
            A0_grad_y(iel,:,:) = matmul( A0(iel,:,:),
     &                                   transpose(var(iel)%grad_y(:,:)))
c... c^h g^{km}_{I} A0 Y_{,m}
            do k = 1,nsd
              k0 = (k-1)*nflow !k0+1, starting index in ri for kth direction
                               !k0+nflow, ending index in ri for kth direction
              ri(iel,k0+1:k0+nflow) = ri(iel,k0+1:k0+nflow)
     &        + ch(iel) 
     &        * matmul( A0_grad_y(iel,:,:), pt_g_p(iel,k,:) )
            enddo
          enddo 
c                
        end subroutine e3if_dc_res
c
        subroutine e3if_dc_egmass(egmass,nshl, shg,A0, ch,gI, WdetJ)  
c..............................................................................
c  calculation of the contribution of DC operator to the local stiffness matrix
c  for both phases
c  \partial G_b / \partial Y_{a} = 
c  c^h g^{km}_{I} N_{a,k} A0 N_{b,m} d \Gama
c..............................................................................
          use propar_m, only: npro
          use conpar_m, only: nflow
          use global_const_m, only: nsd
          implicit none
c
          real*8, dimension(:,:,:),   intent(inout) :: egmass
          real*8, dimension(:,:,:),   intent(in) :: shg
          real*8, dimension(npro,nflow,nflow), intent(in) :: A0
          real*8, dimension(npro), intent(in) :: ch
          real*8, dimension(npro,nsd,nsd),intent(in) :: gI !g^{ij}_{I} = 
                                                           ! proj^T g^{ij} proj
          real*8, dimension(npro),intent(in) :: WdetJ
          integer, intent(in)  :: nshl
c
          real*8, dimension(npro,nsd) :: shga_g
          real*8, dimension(npro) :: shga_g_shgb
          integer :: iel, ia, ib, a_row, b_col
c
          do iel = 1,npro
            do ia = 1,nshl
              a_row = (ia -1)*nflow ! a_row+1: starting row# for dof ia
                                    ! a_row+nflow: ending row# for dof ia
c... g^{km}_{I} N_{ia,k}                                    
              shga_g(iel,:) = matmul( shg(iel,ia,:), gI(iel,:,:) )
              do ib = 1,nshl
                b_col = (ib -1)*nflow ! b_col+1: starting column# for dof ib
                                      ! b_col+nflow: ending column# for dof ib
c... g^{km}_{I} N_{ia,k} N_{ib,m}                                                    
                shga_g_shgb(iel) = dot_product( shga_g(iel,:),shg(iel,ib,:))
c
                egmass(iel, a_row+1:a_row+nflow, b_col+1:b_col+nflow) = 
     &          egmass(iel, a_row+1:a_row+nflow, b_col+1:b_col+nflow)
     &        + ch(iel) * shga_g_shgb(iel) * A0(iel,:, :) * WdetJ(iel)
              enddo
            enddo 
          enddo
c         
        end subroutine e3if_dc_egmass
        
        subroutine e3if_dc
c..............................................................................
c  calculating the discontinuity capturing (DC) operator
c..............................................................................        
          use e3if_param_m
          use e3gij_m, only: e3giju
          implicit none
c
c
          real*8, dimension(npro) :: ch0, ch1 ! nu^h for DC on interface
          real*8, dimension(npro,nsd,nsd) :: proj  ! tangential projector
          real*8, dimension(npro,6) :: giju0, giju1 ! g^{ij} vector form
          real*8, dimension(npro,nsd,nsd) :: giju_f0, giju_f1 ! g^{ij} full
          real*8, dimension(npro,nsd,nsd) :: pt_g0, pt_g1 ! proj_{ik} g^ij
                                                          ! proj is symetric
          real*8, dimension(npro,nsd,nsd) :: pt_g_p0, pt_g_p1 ! proj_{ik} g^ij proj_{jm}
          integer :: iel
c... get the tangential projector          
          call calc_projector(proj, nv0) ! get the tangential projector
c... get the c^h
          call calc_ch(ch0, ch1, f_jump)
c          
          call e3giju (giju0, dxidx0, npro, nsd, lcsyst0) ! get g^{ij}, 0 side
          call e3giju (giju1, dxidx1, npro, nsd, lcsyst1) ! get g^{ij}, 1 side
c.... get the full g^ij matrix
          giju_f0(:,1,1) = giju0(:,1)
          giju_f0(:,1,2) = giju0(:,4)
          giju_f0(:,1,3) = giju0(:,5)
          giju_f0(:,2,1) = giju0(:,4)
          giju_f0(:,2,2) = giju0(:,2)
          giju_f0(:,2,3) = giju0(:,6)
          giju_f0(:,3,1) = giju0(:,5)
          giju_f0(:,3,2) = giju0(:,6)
          giju_f0(:,3,3) = giju0(:,3) 
c
          giju_f1(:,1,1) = giju1(:,1)
          giju_f1(:,1,2) = giju1(:,4)
          giju_f1(:,1,3) = giju1(:,5)
          giju_f1(:,2,1) = giju1(:,4)
          giju_f1(:,2,2) = giju1(:,2)
          giju_f1(:,2,3) = giju1(:,6)
          giju_f1(:,3,1) = giju1(:,5)
          giju_f1(:,3,2) = giju1(:,6)
          giju_f1(:,3,3) = giju1(:,3)
c... calculate proj_{ik} g^ij, notice proj is symetric
          do iel = 1,npro
            pt_g0(iel,:,:) = matmul(proj(iel,:,:), giju_f0(iel,:,:))
            pt_g1(iel,:,:) = matmul(proj(iel,:,:), giju_f1(iel,:,:))
c... calculate proj_{ik} g^ij proj_{jm}
            pt_g_p0(iel,:,:) = matmul(pt_g0(iel,:,:),proj(iel,:,:))
            pt_g_p1(iel,:,:) = matmul(pt_g1(iel,:,:),proj(iel,:,:))             
          enddo         
c... calculate the local residual
          call e3if_dc_res(ri0, var0, A0_0, ch0,pt_g_p0)
          call e3if_dc_res(ri1, var1, A0_1, ch1,pt_g_p1)
c... calculate the local stiffness matrix
          if (lhs_dg .eq. 1) then
            call e3if_dc_egmass(egmass00, nshl0, shg0,A0_0, ch0, pt_g_p0,
     &                          WdetJif0)
            call e3if_dc_egmass(egmass11, nshl1, shg1,A0_1, ch1, pt_g_p1,
     &                          WdetJif1)       
          endif          
c                                                
          
c
        end subroutine e3if_dc  
c
      end module e3if_dc_m
