      module e3if_param_m
c
c----------------------------------------
c    aims to allocate/deallocate the variables
c    for the macro elements for DG interface
c----------------------------------------
        use e3_param_m
        use number_def_m
        use pointer_data
        use elmpar_m
        use shpdat_m
        use conpar_m
        use genpar_m
        use propar_m
        use ifbc_def_m
c
        implicit none
c
        abstract interface
          subroutine calc_vi2
          end subroutine calc_vi2
          subroutine e3ifvar
          end subroutine e3ifvar
        end interface
c
        procedure(e3ifvar), pointer :: get_vap_frac0
c
        integer, dimension(:,:),   pointer :: sgn0, sgn1
c
        real*8,  dimension(:,:,:), pointer :: ycl0, ycl1
        real*8,  dimension(:,:,:), pointer :: acl0, acl1
        real*8,  dimension(:,:,:), pointer :: xl0, xl1      ! nodal coordinates
        real*8,  dimension(:,:,:), pointer :: umeshl0, umeshl1 ! mesh velocity
        real*8,  dimension(:,:),   pointer :: nv0, nv1      ! interface normal vectors
        real*8,  dimension(:,:),   pointer :: vi            ! interface velocity (at integration point)
        real*8,  dimension(:),     pointer :: area
        real*8,  dimension(:),     pointer :: kappa0, kappa1
        real*8,  dimension(:,:,:), pointer :: sum_vi_area_l0, sum_vi_area_l1
        real*8,  dimension(:,:,:), pointer :: if_normal_l0, if_normal_l1,if_kappa_l0,if_kappa_l1
        real*8,  dimension(:,:,:), pointer :: ifbc_l0,ifbc_l1
        real*8,  dimension(:,:,:), pointer :: cmtrx,ctc     ! kinematic continuity matrix C
        real*8,  dimension(:,:),   pointer :: shp0, shp1    ! element shape function at quadrature point
        real*8,  dimension(:,:,:), pointer :: shgl0, shgl1  ! element shape function gradient at a quadrature point
        real*8,  dimension(:,:,:), pointer :: shg0, shg1    ! shape function gradient at a quadrature point
        real*8,  dimension(:,:,:), pointer :: dxdxi0, dxdxi1  ! element deformation tensor
        real*8,  dimension(:),     pointer :: WdetJif0, WdetJif1
c
        real*8, dimension(:,:), pointer :: y0,y1
        real*8, dimension(:,:,:), pointer :: y_jump
        real*8,  dimension(:,:,:), pointer :: rl0, rl1      ! residual over the element
        real*8,  dimension(:,:),   pointer :: ri0, ri1      ! residual at the integration point
c
        real*8, dimension(:,:,:), pointer :: A0_0,A0_1,A0Na_0, A0Na_1
        real*8, dimension(:,:,:,:), pointer :: Ai0, Ai1, AiNa0, AiNa1, KijNaj0, KijNaj1, KijNajC0, KijNajC1
        real*8, dimension(:,:,:,:,:), pointer :: Kij0, Kij1
        real*8, dimension(:,:,:), pointer :: egmass00, egmass01, egmass10, egmass11
        real*8, dimension(:,:,:), pointer :: BDiagl_00, BDiagl_11  ! local Block-Diagonal preconditioner
c
c... evaluated parameters at the integration point
c
        real*8, dimension(:), pointer :: rho0,  rho1     ! density
     &,                                  pres0, pres1    ! pressure
     &,                                  T0,    T1       ! temperature
     &,                                  ei0,   ei1      ! internal energy
     &,                                  c0,    c1       ! speed of sound
     &,                                  vap_frac0       ! mass fraction of vaporized liquid
     &,                                  alpha_LF
        real*8, dimension(:,:), pointer :: u0, u1, um0, um1
c
        real*8, dimension(:), pointer :: rk0, h0, cp0, alfaP0, betaT0
        real*8, dimension(:), pointer :: rk1, h1, cp1, alfaP1, betaT1
c
        integer :: istep
        real*8 :: time,deltat
        real*8 :: s, e, length_h, emu, ek  
        real*8, pointer :: mu(:,:)
c
        integer :: lhs_dg
c
c... properties
c
        type prop_t
          integer :: mater
          real*8 :: rmu, rlm, rlm2mu, con
          real*8, dimension(15,15) :: stiff
        end type prop_t
c
        type qpt_t           ! type def for integration point
          integer :: nshl
          real*8, pointer :: shp(:), shg(:,:), shgl(:,:)  ! shape function and gradient at quadrature point
        end type qpt_t
c
        type var_t            ! type def for integration variables
          real*8 :: rho, u(3), p, T, ei, y(5), grad_y(3,5), grad_yl(3,5)
        end type var_t
c
        type element_t
          real*8 :: dxdxi(3,3)
          real*8 :: WdetJ
        end type element_t
c
        type(qpt_t), dimension(:), pointer :: qpt0, qpt1
        type(var_t),  dimension(:), pointer :: var0,  var1
        type(prop_t), dimension(:), pointer :: prop0, prop1
        type(element_t), dimension(:), pointer :: e0, e1
c
        procedure(getthm2), pointer :: getthmif0_ptr, getthmif1_ptr
        procedure(calc_vi2), pointer :: calc_vi_cavitation
c... added for interface DC operator
        real*8, dimension(:,:), pointer :: f_jump ! flux jump
        real*8, dimension(:,:,:), pointer :: dxidx0, dxidx1  ! inverse of element deformation tensor      
c
      contains
c
        subroutine e3if_malloc
c
          integer :: i
c
          allocate(ycl0(npro,nshl0,ndof))
          allocate(ycl1(npro,nshl1,ndof))
          allocate(umeshl0(npro,nenl0,nsd))
          allocate(umeshl1(npro,nenl1,nsd))
          allocate(vi(npro,nsd))
          allocate(sum_vi_area_l0(npro,nshl0,nsd+1))
          allocate(sum_vi_area_l1(npro,nshl1,nsd+1))
          allocate(ifbc_l0(npro,nshl0,nifbc+1))
          allocate(ifbc_l1(npro,nshl1,nifbc+1))
          allocate(cmtrx(npro,nflow,nflow))
          allocate(ctc(npro,nflow,nflow))
          allocate(acl0(npro,nshl0,ndof))
          allocate(acl1(npro,nshl1,ndof))
c
          allocate(shg0(npro,nshl0,nsd))
          allocate(shg1(npro,nshl1,nsd))
          allocate(dxdxi0(npro,nsd,nsd),dxdxi1(npro,nsd,nsd))
          allocate(dxidx0(npro,nsd,nsd),dxidx1(npro,nsd,nsd))
c
          allocate(rl0(npro,nshl0,nflow))
          allocate(rl1(npro,nshl1,nflow))
c
          allocate(y0(npro,nflow),y1(npro,nflow))
          allocate(y_jump(npro,nflow,nsd))
          allocate(f_jump(npro,nflow))
          allocate(ri0(npro,nflow*(nsd+1)))
          allocate(ri1(npro,nflow*(nsd+1)))
          allocate(rho0(npro),u0(npro,nsd),pres0(npro),T0(npro),ei0(npro),um0(npro,nsd))
          allocate(rho1(npro),u1(npro,nsd),pres1(npro),T1(npro),ei1(npro),um1(npro,nsd))
          allocate(rk0(npro),h0(npro),cp0(npro),alfaP0(npro),betaT0(npro))
          allocate(rk1(npro),h1(npro),cp1(npro),alfaP1(nprO),betaT1(npro))
          allocate(c0(npro), c1(npro))
          allocate(vap_frac0(npro))
          allocate(alpha_LF(npro))
c
          allocate(qpt0(npro),qpt1(npro))
          allocate(var0 (npro), var1 (npro))
          allocate(prop0(npro), prop1(npro))
          allocate(e0(npro),e1(npro))
c
          do i = 1,npro
            allocate(qpt0(i)%shp(nshl0))
            allocate(qpt1(i)%shp(nshl1))
            allocate(qpt0(i)%shg(nsd,nshl0))
            allocate(qpt1(i)%shg(nsd,nshl1))
            allocate(qpt0(i)%shgl(nsd,nshl0))
            allocate(qpt1(i)%shgl(nsd,nshl1))
          enddo
c
          allocate(A0_0(npro,nflow,nflow))
          allocate(A0_1(npro,nflow,nflow))
          allocate(A0Na_0(npro,nflow,nflow))
          allocate(A0Na_1(npro,nflow,nflow))
          allocate(Ai0(npro,nsd,nflow,nflow))
          allocate(Ai1(npro,nsd,nflow,nflow))
          allocate(AiNa0(npro,nsd,nflow,nflow))
          allocate(AiNa1(npro,nsd,nflow,nflow))
          allocate(Kij0(npro,nsd,nsd,nflow,nflow))
          allocate(Kij1(npro,nsd,nsd,nflow,nflow))
          allocate(KijNaj0(npro,nsd,nflow,nflow))
          allocate(KijNaj1(npro,nsd,nflow,nflow))
          allocate(KijNajC0(npro,nsd,nflow,nflow))
          allocate(KijNajC1(npro,nsd,nflow,nflow))
c
          allocate(kappa0(npro),kappa1(npro))
c
          allocate(mu(npro,nflow))
c
          allocate(BDiagl_00(npro,nshl0,nflow*nflow))
          allocate(BDiagl_11(npro,nshl1,nflow*nflow))          
c
        end subroutine e3if_malloc
c
        subroutine e3if_mfree
c
          integer :: i
c
          deallocate(ycl0,ycl1)
          deallocate(umeshl0,umeshl1)
          deallocate(vi)
          deallocate(sum_vi_area_l0,sum_vi_area_l1)
          deallocate(ifbc_l0,ifbc_l1)
          deallocate(cmtrx)
          deallocate(ctc)
          deallocate(acl0,acl1)
          deallocate(shg0,shg1)
          deallocate(dxdxi0,dxdxi1)
          deallocate(dxidx0,dxidx1)
          deallocate(y0,y1)
          deallocate(y_jump)
          deallocate(f_jump)
          deallocate(rl0,rl1)
          deallocate(ri0,ri1)
          deallocate(rho0,u0,pres0,T0,ei0,um0)
          deallocate(rho1,u1,pres1,T1,ei1,um1)
          deallocate(rk0,h0,cp0,alfaP0,betaT0)
          deallocate(rk1,h1,cp1,alfaP1,betaT1)
          deallocate(c0,c1)
          deallocate(vap_frac0)
          deallocate(alpha_LF)
c
          do i = 1,npro
            deallocate(qpt0(i)%shp,qpt1(i)%shp)
            deallocate(qpt0(i)%shg,qpt1(i)%shg)
            deallocate(qpt0(i)%shgl,qpt1(i)%shgl)
          enddo
c
          deallocate(qpt0,qpt1)
          deallocate(var0, var1)
          deallocate(prop0, prop1)
          deallocate(e0,e1)
c
          deallocate(Ai0,Ai1)
          deallocate(A0_0,A0_1)
          deallocate(A0Na_0,A0Na_1)
          deallocate(AiNa0,AiNa1)
          deallocate(Kij0,Kij1)
          deallocate(KijNaj0,KijNaj1)
          deallocate(KijNajC0,KijNajC1)
c
          nullify(egmass00)
          nullify(egmass01)
          nullify(egmass10)
          nullify(egmass11)
c
          deallocate(kappa0,kappa1)
c
          deallocate(mu)
c
          deallocate(BDiagl_00, BDiagl_11)
c
        end subroutine e3if_mfree
c
        subroutine e3if_empty
          return
        end subroutine e3if_empty
c
        subroutine e3if_setparam_0
          rho  => rho0
          ei   => ei0
          pres => pres0
          T    => T0
          h    => h0
          cp   => cp0
          alphaP => alfaP0
          betaT => betaT0
          c     => c0
          mater = mater0
          vap_frac => vap_frac0
        end subroutine e3if_setparam_0
c
        subroutine e3if_setparam_1
          rho  => rho1
          ei   => ei1
          pres => pres1
          T    => T1
          h    => h1
          cp   => cp1
          alphaP => alfaP1
          betaT => betaT1
          c     => c1
          mater = mater1
        end subroutine e3if_setparam_1
c
      end module e3if_param_m
