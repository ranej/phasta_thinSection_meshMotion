      module e3if_m
c
c----------------------------------------
c   this is the main subroutine in calculating
c   the element level contribution for the
c   interface. It calls many other subroutines
c   to collect the data for both RHS and LHS
c----------------------------------------
        use workfc_m
        use hierarchic_m
        use matdat_def_m
        use e3if_param_m
        use e3if_geom_m
        use e3if_func_m
        use e3if_diff_m
        use eqn_state_m
        use e3metric_m
        use e3if_lhs_m
        use e3if_vi_m
        use if_global_m
        use e3if_dc_m ! DC operator for interface
        use weighted_normal_func_m ! for weighted normal
        use weighted_normal_data_m ! for weighted normal
        use dgifinp_m, only: i_w_normal, i_if_dc
c
        implicit none
c
      contains

        subroutine e3if(shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
          real*8, dimension(nshl0,nqpt), intent(in) :: shpif0
          real*8, dimension(nshl1,nqpt), intent(in) :: shpif1
          real*8, dimension(nsd,nshl0,nqpt), intent(in) :: shgif0
          real*8, dimension(nsd,nshl1,nqpt), intent(in) :: shgif1
          real*8, dimension(nqpt), intent(in) :: qwtif0, qwtif1
c
          integer :: intp, intp0, intp1
      integer :: iel,isd,n
      real*8 ::sum0,sumg0
      real*8, allocatable :: tmpmu0(:,:),tmpmu1(:,:)
#define debug 0
c
c      write(*,*) 'In e3if...'
c
c      do n = 1,nshl0
c        write(*,33) n,shpif0(n,1),shgif0(1:3,n,1)
c      enddo
c      write(*,34) sum(shpif0(1:nshl0,1)),sum(shgif0(1,1:nshl0,1)),sum(shgif0(2,1:nshl0,1)),sum(shgif0(3,1:nshl0,1))
33    format('n,shpif0,shgif0',i4,4e24.16)
34    format('sum: shpif0,shgif0',4e24.16)
c
          rl0 = zero
          rl1 = zero
c
          sum_vi_area_l0 = zero
          sum_vi_area_l1 = zero
c
          ifbc_l0 = zero
c
c
c------------------------------------------------------------------------------------
c------------ THIS IMPLEMENTATION IS VALID ONLY FOR LINEAR TETS -----------------------
c------------------------------------------------------------------------------------
c
c
c------------------- Begin loop over quadrature points ----------------------------
          do intp = 1, nqpt
c		
            ri0 = zero
            ri1 = zero
            intp0 = intp
            intp1 = mod(2*intp+1,3)+1
c node #1 is matched at the chef level -- nodes 2 and 3 always unmatched at interface
c intp0 and intp1 intoduced to match physical location of integration points
c
c... do not include this quadrature point, if Det. .eq. 0
c
            if (qwtif0(intp0) .eq. zero) cycle
            if (qwtif1(intp1) .eq. zero) cycle
c
c... create a matrix of shape fucntions (and derivatives) for each
c    element at this quadrature point. These arrays will contain
c    the correct signs for the higher order hierarchic basis
c
            call  getshp(shp0, shgl0, shpif0, shgif0, sgn0, npro, nsd, nshl0, nqpt, nenl0, intp0, ipord)
            call  getshp(shp1, shgl1, shpif1, shgif1, sgn1, npro, nsd, nshl1, nqpt, nenl1, intp1, ipord)
c      call getshp_if(shp0,shp1,shgl0,shgl1,shpif0,shpif1,shgif0,shgif1,xl0,xl1,npro,nsd,nshl0,nqpt,nenl0,intp)
c
c... Element Metrics
c
            call e3metric(shg0,dxidx0,shgl0,xl0)
            call e3metric(shg1,dxidx1,shgl1,xl1)
c      print*, 'shg0: ',shg0(1,:,1)
c      print*, 'shg1: ',shg1(1,:,1)
c
c... Normal Vectors
c
            call calc_normal_vectors(nv0,area,WdetJif0,xl0,qwtif0,itpid,lcsyst0,intp0,npro)
            call calc_normal_vectors(nv1,area,WdetJif1,xl1,qwtif1,itpid,lcsyst1,intp1,npro)
c
c... replace the normal with the weighted normal if needed
            if (i_w_normal .eq. 1) then
              call get_weighted_normal(nv0, w_normal_l0, shp0, nshl0)                        
              call get_weighted_normal(nv1, w_normal_l1, shp1, nshl1)
            endif
c
c... calculate the integration varibles
c
            call e3if_var
c
            call e3if_mtrx
c
c... calculate the contribution of the jump in the fluxes, across the interface
c
      prop0%mater = mater0
      prop1%mater = mater1
c
            call e3var(y0, var0, ycl0, shp0, shgl0, shg0, nshl0) 
            call e3var(y1, var1, ycl1, shp1, shgl1, shg1, nshl1) 
c
            call calc_stiff(prop0, var0, mater0)
            call calc_stiff(prop1, var1, mater1)
c
c ... Interface flux 
c
            call e3if_flux
c
c            call flux_jump
c
            call calc_cmtrx
            call calc_y_jump
c
            call kinematic_condition(ri0,Kij0)
            call kinematic_condition(ri1,Kij1)
c
c ... kinematic condition term:
c     set the mu coeff to the max of the 0,1 materials 
c
      allocate(tmpmu0(npro,nflow),tmpmu1(npro,nflow))
c
      select case (mat_eos(mater0,1))
      case (ieos_ideal_gas,ieos_ideal_gas_mixture,ieos_liquid_1)
        tmpmu0(:,1) = zero
        tmpmu0(:,2) = prop0%stiff(3,3)
        tmpmu0(:,3) = prop0%stiff(3,3)
        tmpmu0(:,4) = prop0%stiff(3,3)
        tmpmu0(:,5) = prop0%stiff(5,5)
      case (ieos_solid_1)
C
C Yu please fill this out
c
      case default
        call error ('getthm  ', 'wrong material', mater)
      end select
c
      select case (mat_eos(mater1,1))
      case (ieos_ideal_gas,ieos_ideal_gas_mixture,ieos_liquid_1)
        tmpmu1(:,1) = zero
        tmpmu1(:,2) = prop1%stiff(3,3)
        tmpmu1(:,3) = prop1%stiff(3,3)
        tmpmu1(:,4) = prop1%stiff(3,3)
        tmpmu1(:,5) = prop1%stiff(5,5)
      case (ieos_solid_1)
C
C Yu please fill this out
c
      case default
        call error ('getthm  ', 'wrong material', mater)
      end select
c
      mu(:,1) = zero
      mu(:,2) = max(tmpmu0(:,2),tmpmu1(:,2))*emu
      mu(:,3) = max(tmpmu0(:,3),tmpmu1(:,3))*emu
      mu(:,4) = max(tmpmu0(:,4),tmpmu1(:,4))*emu
      mu(:,5) = max(tmpmu0(:,5),tmpmu1(:,5))*ek
c
      deallocate(tmpmu0,tmpmu1)
c
            call dg_penalty(ri0,y0,y1)
            call dg_penalty(ri1,y1,y0)
c... discontinuity capturing term for the interface
            if (i_if_dc .eq. 1) then
              call e3if_dc
            endif
c...LHS calculations...
c
            if (lhs_dg .eq. 1) then
c
c              call set_lhs_matrices ! Unnecessary now
c
c              call calc_egmass(egmass00,egmass01,
c     &                         A0_0, A0_1, Ai0, Ai1,
c     &                         Kij0, Kij1,
c     &                         AiNa0,AiNa1,KijNaj0,KijNaj1,KijNajC0,KijNajC1,
c     &                         shp0,nv0,nv1,WdetJif0,prop0,nshl0,nshl1)
c
c              call calc_egmass(egmass11,egmass10,
c     &                         A0_1, A0_0, Ai1, Ai0,
c     &                         Kij1, Kij0,
c     &                         AiNa1,AiNa0,KijNaj1,KijNaj0,KijNajC1,KijNajC0,
c     &                         shp1,nv1,nv0,WdetJif1,prop1,nshl1,nshl0)

c               call calc_egmass_(egmass00,shp0,shp0,shg0,shg0,Ai0,Kij0,Kij0,nv0,nv0,WdetJif0,nshl0,nshl0)
c               call calc_egmass_(egmass01,shp0,shp1,shg0,shg1,Ai1,Kij0,Kij1,nv0,nv1,WdetJif0,nshl0,nshl1)
c               call calc_egmass_(egmass10,shp1,shp0,shg1,shg0,Ai0,Kij1,Kij0,nv1,nv0,WdetJif1,nshl1,nshl0)
c               call calc_egmass_(egmass11,shp1,shp1,shg1,shg1,Ai1,Kij1,Kij1,nv1,nv1,WdetJif1,nshl1,nshl1)
               call calc_egmass_fix_sign(egmass00,shp0,shp0,shg0,shg0,Ai0,Kij0,Kij0,nv0,nv0,WdetJif0,nshl0,nshl0,'same')
               call calc_egmass_fix_sign(egmass01,shp0,shp1,shg0,shg1,Ai1,Kij0,Kij1,nv0,nv1,WdetJif0,nshl0,nshl1,'diff')
               call calc_egmass_fix_sign(egmass10,shp1,shp0,shg1,shg0,Ai0,Kij1,Kij0,nv1,nv0,WdetJif1,nshl1,nshl0,'diff')
               call calc_egmass_fix_sign(egmass11,shp1,shp1,shg1,shg1,Ai1,Kij1,Kij1,nv1,nv1,WdetJif1,nshl1,nshl1,'same')
c
            endif
c
            call e3if_wmlt(rl0, ri0, shp0, shg0, WdetJif0, nshl0)
            call e3if_wmlt(rl1, ri1, shp1, shg1, WdetJif1, nshl1)
c      do iel = 1,npro
c        write(*,10) 'E3IF: ri0: ',iel,ri0(iel,16:20)
c        write(*,10) 'E3IF: ri1: ',iel,ri1(iel,16:20)
c      enddo
c      do iel = 1,1
c        do n = 1,nshl0
c          write(*,10) 'E3IF: rl0: ',n,rl0(iel,n,:)
c        enddo
c      enddo
c      do iel = 1,1
c        do n = 1,nshl1
c          write(*,10) 'E3IF: rl1: ',n,rl1(iel,n,:)
c        enddo
c      enddo
c      do iel = 1,8
c        write(*,10) 'WdetJ0,1: ',iel, WdetJif0(iel),WdetJif1(iel)
c      enddo
c      write(*,10) 'qwtif0,1: ',iel, qwtif0(intp),qwtif1(intp)
10    format(a,i6,5e24.16)
c
c ... Interface velocity calculations...
c
            call calc_vi(pres0,u1)
c
            call calc_vi_area_node(sum_vi_area_l0,shp0,WdetJif0,nshl0)
            call calc_vi_area_node(sum_vi_area_l1,shp1,WdetJif1,nshl1)
c
c            call calc_vapor_frac_node
c
          enddo 
c------------------- End of Quadrature loop --------------------------------
        end subroutine e3if
c
        subroutine e3var(y,var,ycl,shp,shgl,shg,nshl)
c
          real*8, dimension(:,:), intent(out) :: y
          type(var_t), pointer, intent(out) :: var(:)
          real*8, pointer, intent(in) :: shp(:,:),shgl(:,:,:), shg(:,:,:)
          real*8, dimension(npro,nshl,nflow), intent(in) :: ycl
          integer, intent(in) :: nshl
c
          integer :: iel,ivar,isd,jsd,n
          real*8 :: grad_y(npro)
c
          do ivar = 1,nflow
            y(:,ivar) = zero
            var(:)%y(ivar) = zero
            var(:)%grad_y(1,ivar) = zero
            var(:)%grad_y(2,ivar) = zero
            var(:)%grad_y(3,ivar) = zero
            do n = 1,nshl
              y(:,ivar) = y(:,ivar) + ycl(:,n,ivar)*shp(:,n)
              var(:)%y(ivar) = var(:)%y(ivar) + ycl(:,n,ivar)*shp(:,n)
              var(:)%grad_y(1,ivar) = var(:)%grad_y(1,ivar) + ycl(:,n,ivar)*shg(:,n,1)
              var(:)%grad_y(2,ivar) = var(:)%grad_y(2,ivar) + ycl(:,n,ivar)*shg(:,n,2)
              var(:)%grad_y(3,ivar) = var(:)%grad_y(3,ivar) + ycl(:,n,ivar)*shg(:,n,3)
c      if (ivar==1) then
c        write(*,10) 'node, shg, ycl, grad',n,shg(1,n,1:3),ycl(1,n,ivar),var(1)%grad_y(1:3,ivar)
c      endif
            enddo
          enddo
10    format(a,i4,7e24.16)

          return
c
        end subroutine e3var
c
        subroutine     e3if_var
c
          implicit none 
c
          integer :: iel
          real*8, dimension(npro) :: tmp
c
c... compute variables at the integration point
c
          call get_var(npro,nshl0,pres0,u0,T0,ycl0,shp0,var0)
          call get_var(npro,nshl1,pres1,u1,T1,ycl1,shp1,var1)
c
          call get_vap_frac0
c
          call get_mesh_velocity(um0,umeshl0,shp0,npro,nshl0)
          call get_mesh_velocity(um1,umeshl1,shp1,npro,nshl1)
c
          call getthmif0
          call getthmif1
c
        end subroutine e3if_var
c
        subroutine e3if_flux
c
          integer :: iel,iflow,isd,jsd,jflow,n
          real*8, dimension(nsd,nflow) :: f0, f1, fconv0, fconv1, fdiff0, fdiff1
          real*8, dimension(nflow) :: f0n0, f0n1, f1n0, f1n1
c
          real*8 :: etot, diff_flux(nsd,nflow), dTdx0
          real*8 :: kappa0(nsd), kappa1(nsd), k0,k1 ! mean curvature
c
          real*8 :: alpha,jump_u(5),climit,jump_y(5),A0_jump_y(5)
c
          element_loop: do iel = 1,npro
c
            call calc_conv_flux(fconv0,rho0(iel),u0(iel,:),um0(iel,:),pres0(iel),ei0(iel))
            call calc_conv_flux(fconv1,rho1(iel),u1(iel,:),um1(iel,:),pres1(iel),ei1(iel))
c
            call calc_diff_flux(fdiff0,var0(iel),prop0(iel))
            call calc_diff_flux(fdiff1,var1(iel),prop1(iel))
c        write(*,500) myrank,iel,fconv0(1,:)
c        write(*,500) myrank,iel,fconv1(1,:)
c        write(*,500) myrank,iel,fdiff0(1,:)
c        write(*,500) myrank,iel,fdiff1(1,:)
c        write(*,500) myrank,iel,fdiff0(:,5)-fdiff1(:,5)
c        write(*,500) myrank,iel,var0(iel)%grad_y(:,5)-var1(iel)%grad_y(:,5)
c
c... calculate flux in normal direction...
c
            do iflow = 1,nflow
              f0(:,iflow) = fconv0(:,iflow) - fdiff0(:,iflow)
              f1(:,iflow) = fconv1(:,iflow) - fdiff1(:,iflow)
              f0n0(iflow) = dot_product(f0(:,iflow),nv0(iel,:))
              f0n1(iflow) = dot_product(f0(:,iflow),nv1(iel,:))
              f1n0(iflow) = dot_product(f1(:,iflow),nv0(iel,:))
              f1n1(iflow) = dot_product(f1(:,iflow),nv1(iel,:))
            enddo
c        write(*,500) myrank,iel,f0n0(:)
c
c
c      if (iel == 1) then
c        write(*,10) 'rho0, u0, p0, T0:', rho0(iel), u0(iel,1), pres0(iel), ei0(iel)
c        write(*,10) 'conv0: ',fconv0(1,:)
c        write(*,10) 'diff0: ',fdiff0(1,:)
c        write(*,11) 'grad_y u 0:', iel, var0(iel)%grad_y(:,2)
c        write(*,11) 'grad_y T 0:', iel, var0(iel)%grad_y(:,5)
c        write(*,10) 'f0: ',f0(1,:)
c        write(*,*) 'stiff0: ',prop0(1)%stiff(15,15)
c        write(*,10) 'rho1, u1, p1, ei1:', rho1(iel), u1(iel,1), pres1(iel), ei1(iel)
c        write(*,10) 'conv1: ',fconv1(1,:)
c        write(*,10) 'diff1: ',fdiff1(1,:)
c        write(*,11) 'grad_y u 1:', iel, var1(iel)%grad_y(:,2)
c        write(*,11) 'grad_y T 1:', iel, var1(iel)%grad_y(:,5)
c        write(*,10) 'f1: ',f1(1,:)
c        write(*,*) 'stiff1: ',prop1(1)%stiff(15,15)
c        write(*,11) 'f0n0:',iel,f0n0
c        write(*,11) 'f1n0:',iel,f1n0
c        write(*,11) 'f1n1:',iel,f1n1
c        write(*,11) 'f0n1:',iel,f0n1
c      endif
c... calculating the flux jump
            f_jump(iel,1:5) =  f1n1(1:5) + f0n0(1:5)
c
            ri0(iel,16:20) = ri0(iel,16:20) + pt50 * ( f0n0(1:5) + f1n0(1:5) )
            ri1(iel,16:20) = ri1(iel,16:20) + pt50 * ( f1n1(1:5) + f0n1(1:5) )
c
#if debug==1
      if (iel == 1) then
        write(*,11) 'ri0  : ',iel,ri0(iel,16:20)
        write(*,11) 'ri1  : ',iel,ri1(iel,16:20)
      endif
#endif
c
c...UPWIND????
c
C      ri0(iel,16:20) = ri0(iel,16:20) + f1n0(1:5)
C      ri1(iel,16:20) = ri1(iel,16:20) + f1n1(1:5)
c
c
c... Here is the additional stability terms from the Lax-Friedrichs flux calculations
c
            climit = zero
c            climit = one
CC            climit = 1.e-1
            alpha_LF(iel) = climit * max(abs(dot_product(u0(iel,:)-um0(iel,:),nv0(iel,:))-c0(iel)),
     &                  abs(dot_product(u1(iel,:)-um1(iel,:),nv1(iel,:))-c1(iel)))
            alpha = alpha_LF(iel)
c
            jump_u(1) = rho0(iel) - rho1(iel)
            jump_u(2) = rho0(iel)*u0(iel,1) - rho1(iel)*u1(iel,1)
            jump_u(3) = rho0(iel)*u0(iel,2) - rho1(iel)*u1(iel,2)
            jump_u(4) = rho0(iel)*u0(iel,3) - rho1(iel)*u1(iel,3)
            jump_u(5) = rho0(iel)*(ei0(iel)+pt50*dot_product(u0(iel,:),u0(iel,:))) 
     &                - rho1(iel)*(ei1(iel)+pt50*dot_product(u1(iel,:),u1(iel,:)))
c
            jump_y(1) = pres0(iel) - pres1(iel)
            jump_y(2) = u0(iel,1) - u1(iel,1)
            jump_y(3) = u0(iel,2) - u1(iel,2)
            jump_y(4) = u0(iel,3) - u1(iel,3)
            jump_y(5) = T0(iel) - T1(iel)
c
            A0_jump_y = zero
c
            do iflow = 1,nflow
              do jflow = 1,nflow
                A0_jump_y(iflow) = A0_jump_y(iflow) + pt50*(A0_0(iel,iflow,jflow)+A0_1(iel,iflow,jflow))*jump_y(jflow)
              enddo 
            enddo
c
CC      ri0(iel,16:20) = ri0(iel,16:20) + pt50*alpha*jump_u(1:5)
CC      ri1(iel,16:20) = ri1(iel,16:20) - pt50*alpha*jump_u(1:5)
c      write(*,11) 'ri0 a: ',iel,ri0(iel,16:20)
c      write(*,11) 'ri1 a: ',iel,ri1(iel,16:20)
      ri0(iel,16:20) = ri0(iel,16:20) + alpha*A0_jump_y(1:5)
      ri1(iel,16:20) = ri1(iel,16:20) - alpha*A0_jump_y(1:5)
c
C... Do we account for surface tension in jump?
c
            if (associated(if_kappa)) then
c
              kappa0 = zero
              kappa1 = zero
              do n = 1,3
                kappa0 = kappa0 + shp0(iel,n)*if_kappa_l0(iel,n,1:nsd)
                kappa1 = kappa1 + shp1(iel,n)*if_kappa_l1(iel,n,1:nsd)
              enddo
              k0 = dot_product(kappa0,nv0(iel,:))
              k1 = dot_product(kappa1,nv1(iel,:))
c
c              ri0(iel,17:19) = ri0(iel,17:19) + pt50 * surface_tension_coeff * k0 * nv0(iel,1:nsd)
c              ri1(iel,17:19) = ri1(iel,17:19) + pt50 * surface_tension_coeff * k1 * nv1(iel,1:nsd)
              ri0(iel,17:19) = ri0(iel,17:19) + pt50 * surface_tension_coeff * kappa0
              ri1(iel,17:19) = ri1(iel,17:19) + pt50 * surface_tension_coeff * kappa1
c
            endif
c
          enddo element_loop
c
10    format(a,5e24.16)
11    format(a,i6,5e24.16)
20    format(a,1e24.16)
500   format('[',i2,'] ',i3,x,5e24.16)
c
        end subroutine e3if_flux
c
c
c
c        subroutine flux_jump
c
c          real*8, dimension(npro) :: vi0,vi1,etot0,etot1
c
c          vi0 = + (vi(:,1)-um0(:,1))*nv0(:,1)
c    &          + (vi(:,2)-um0(:,2))*nv0(:,2)
c     &          + (vi(:,3)-um0(:,3))*nv0(:,3)
c
c          vi1 = + (vi(:,1)-um1(:,1))*nv1(:,1)
c     &          + (vi(:,2)-um1(:,2))*nv1(:,2)
c     &          + (vi(:,3)-um1(:,3))*nv1(:,3)
c
c      write(*,*) 'vi0:',vi0
c      write(*,*) 'vi1:',vi1
c
c          etot0 = ei0 + pt50 * (u0(:,1)*u0(:,1)+u0(:,2)*u0(:,2)+u0(:,3)*u0(:,3))
c          etot1 = ei1 + pt50 * (u1(:,1)*u1(:,1)+u1(:,2)*u1(:,2)+u1(:,3)*u1(:,3))
c
c          ri0(:,16) = ri0(:,16) + vi0*rho0
c          ri0(:,17) = ri0(:,17) + vi0*rho0*u0(:,1)
c          ri0(:,18) = ri0(:,18) + vi0*rho0*u0(:,2)
c          ri0(:,19) = ri0(:,19) + vi0*rho0*u0(:,3)
c          ri0(:,20) = ri0(:,20) + vi0*rho0*etot0
c
c          ri1(:,16) = ri1(:,16) + vi1*rho1
c          ri1(:,17) = ri1(:,17) + vi1*rho1*u1(:,1)
c          ri1(:,18) = ri1(:,18) + vi1*rho1*u1(:,2)
c          ri1(:,19) = ri1(:,19) + vi1*rho1*u1(:,3)
c          ri1(:,20) = ri1(:,20) + vi1*rho1*etot1
c
c        end subroutine flux_jump
c
c
c
        subroutine kinematic_condition(ri,Kij)
c
           real*8, dimension(:,:), intent(inout) :: ri
           real*8, dimension(:,:,:,:,:), intent(in) :: Kij
c
           integer :: iflow,jflow,kflow,isd,jsd
           real*8 :: this_kcy(npro)
           real*8, dimension(npro,nflow,nflow,nsd,nsd) :: CKij
           real*8,dimension(npro,nflow,nsd) :: cy_jump, kcy
c
           do iflow = 1,nflow
c
             cy_jump(:,iflow,:) = zero
c
             do jflow = 1,nflow
               cy_jump(:,iflow,1) = cy_jump(:,iflow,1) + cmtrx(:,iflow,jflow)*y_jump(:,jflow,1)
               cy_jump(:,iflow,2) = cy_jump(:,iflow,2) + cmtrx(:,iflow,jflow)*y_jump(:,jflow,2)
               cy_jump(:,iflow,3) = cy_jump(:,iflow,3) + cmtrx(:,iflow,jflow)*y_jump(:,jflow,3)
             enddo
c
           enddo
c
           do iflow = 1,nflow
             do isd = 1,nsd
c
             this_kcy = zero
c
             do jflow = 1,nflow
               do jsd = 1,nsd
!                 this_kcy = this_kcy + Kij(:,iflow,jflow,isd,jsd)*cy_jump(:,jflow,jsd)
		this_kcy = this_kcy + Kij(:,jsd,isd,iflow,jflow)*cy_jump(:,jflow,jsd) !DOUBLE CHECK THIS
               enddo
             enddo
c
             ri(:,nflow*(isd-1)+iflow) = ri(:,nflow*(isd-1)+iflow) + pt50 * s * this_kcy
c      write(*,*) 'KINEMATIC: ',ri(1,nflow*(isd-1)+iflow),this_kcy(1)
c
             enddo
           enddo
c
        end subroutine kinematic_condition
c
        subroutine dg_penalty(ri,y0,y1)
c
           real*8, dimension(:,:), intent(inout) :: ri
           real*8, dimension(:,:), intent(in) :: y0,y1
c
           integer :: iflow,jflow,isd
           real*8 :: this_sum(npro)
c
           do iflow = 1,nflow
c
               this_sum = zero
c
               do jflow = 1,nflow
                 this_sum = this_sum + ctc(:,iflow,jflow)*(y0(:,jflow)-y1(:,jflow))
               enddo
c
               ri(:,3*nflow+iflow) = ri(:,3*nflow+iflow) + e*mu(:,iflow)/length_h * this_sum
c
           enddo
c
        end subroutine dg_penalty
c
        subroutine calc_cmtrx
c
          integer :: p,q,r
c
          cmtrx = zero
c
c          cmtrx(:,1,1) = one
c
          cmtrx(:,2,2) = one - nv0(:,1)*nv0(:,1)
          cmtrx(:,2,3) =     - nv0(:,1)*nv0(:,2)
          cmtrx(:,2,4) =     - nv0(:,1)*nv0(:,3)
c
          cmtrx(:,3,2) =     - nv0(:,2)*nv0(:,1)
          cmtrx(:,3,3) = one - nv0(:,2)*nv0(:,2)
          cmtrx(:,3,4) =     - nv0(:,2)*nv0(:,3)
c
          cmtrx(:,4,2) =     - nv0(:,3)*nv0(:,1)
          cmtrx(:,4,3) =     - nv0(:,3)*nv0(:,2)
          cmtrx(:,4,4) = one - nv0(:,3)*nv0(:,3)
c
          cmtrx(:,5,5) = one
c
c...NOTE: need to derive an expression for it, instead of this:
c
          do q = 1,nflow
            do p = 1,nflow
              ctc(:,p,q) = zero
              do r = 1,nflow
                ctc(:,p,q) = ctc(:,p,q) + cmtrx(:,p,r)*cmtrx(:,r,q)
              enddo
            enddo
          enddo
c
        end subroutine calc_cmtrx
c
        subroutine calc_conv_flux(flux,rho,u,um,p,ei)
c
          real*8, dimension(nsd,nflow), intent(out) :: flux
          real*8, intent(in) :: rho,p,ei,u(nsd),um(nsd)
          real*8 :: delta(3,3)
          data delta /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
c
          integer :: i
          real*8  :: etot
c
          etot = ei + pt50 * dot_product(u,u)
c
          do i = 1,nsd
c
            flux(i,1) = rho*(u(i)-um(i))
            flux(i,2) = rho*(u(i)-um(i))*u(1) + p*delta(1,i)
            flux(i,3) = rho*(u(i)-um(i))*u(2) + p*delta(2,i)
            flux(i,4) = rho*(u(i)-um(i))*u(3) + p*delta(3,i)
            flux(i,5) = rho*(u(i)-um(i))*etot + p*u(i)
c
          enddo
c
        end subroutine calc_conv_flux
c
        subroutine calc_diff_flux(flux,var,prop)
c
          real*8, dimension(nsd,nflow), intent(out) :: flux
          type(var_t), intent(in) :: var
          type(prop_t), intent(in) :: prop
c
          integer :: isd,jsd,iflow,jflow
          real*8  :: this_sum
c
c... diffusive flux
c
            do iflow = 1,nflow
              do isd = 1,nsd
                this_sum = zero
                ! multiply matrix: K_ij * Y,j
                do jflow = 1,nflow
                  do jsd = 1,nsd
                    this_sum = this_sum +
     &                prop%stiff((isd-1)*nflow+iflow,(jsd-1)*nflow+jflow)*var%grad_y(jsd,jflow)
                  enddo
                enddo
                flux(isd,iflow) = this_sum
              enddo
            enddo
c
        end subroutine calc_diff_flux
c
        subroutine e3if_wmlt(rl,ri,shp,shg,WdetJ,nshl)
c
          integer, intent(in) :: nshl
          real*8, dimension(:,:,:), pointer, intent(out) :: rl
          real*8, dimension(:,:),   pointer, intent(in)  :: ri,shp
          real*8, dimension(:,:,:), pointer, intent(in)  :: shg
          real*8, dimension(:),     pointer, intent(in)  :: WdetJ
c
          integer :: n,iflow,isd
c
          do isd = 1,nsd
            do iflow = 1,nflow
              do n = 1, nshl
                rl(:,n,iflow) = rl(:,n,iflow) +
     &            WdetJ * shg(:,n,isd) * ri(:,nflow*(isd-1)+iflow)
              enddo
            enddo
          enddo
c
          do iflow = 1,nflow
            do n = 1,nshl
              rl(:,n,iflow) = rl(:,n,iflow) + shp(:,n)*WdetJ(:) * ri(:,3*nflow+iflow)
            enddo
          enddo
c
        end subroutine e3if_wmlt
c
      end module e3if_m
