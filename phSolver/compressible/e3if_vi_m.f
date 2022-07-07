      module e3if_vi_m
c
c----------------------------------------
c    aims to calculate the interface
c    velocity at element level
c----------------------------------------
        use dgifinp_m
        use e3if_param_m
        use matdat_def_m
c
        implicit none
c
      contains
c
      subroutine calc_vi(p,u)
c
        use sclrs_m
c
        real*8, pointer, intent(in) :: p(:), u(:,:)
        real*8 :: vi_max
c
        integer :: isd,i,n,iel
        real*8 :: vimag,v1,t1,c1
        real*8 :: pvap,cprod,cdest
c... Clausius-Clapeyron:
        real*8, dimension(npro) :: x_vapor     ! vapor mole fraction in gas/vapor mixture
     &,                            y_vapor     ! vapor mass fraction ...
     &,                            dydn        ! vapor mass fraction gradient in normal direction of liquid 1
     &,                            rho_mix     ! mixture density
     &,                            mw_mix      ! mixture molecular weight
     &,                            vap_rate
     &,                            un0,un1     ! velocity in normal direction
c
#define debug 0
        select case (vi_ramping)
        case (no_ramp)
          c1 = one
        case (linear_ramp)
          t1 = ramp_time
          v1 = vimag
c
          if (time <= t1) then
            c1 = time/t1
          else
            c1 = one
          endif
        case default
          call error ('ERROR in e3if_vi:',' vi_ramping is not supported,',vi_ramping)
        end select
c
        select case (phase_change_model)
        case (no_vi)
          vi = zero
        case (const_vi)
          vi(:,1) = c1 * (vi_mag * nv0(:,1) + u1(:,1))
          vi(:,2) = c1 * (vi_mag * nv0(:,2) + u1(:,2))
          vi(:,3) = c1 * (vi_mag * nv0(:,3) + u1(:,3))
c          vi(:,1) = vi_mag * nv0(:,1)
c          vi(:,2) = vi_mag * nv0(:,2)
c          vi(:,3) = vi_mag * nv0(:,3)
c      write(*,100) 'vi_mag: ',vi_mag
c      write(*,100) 'vi    : ',vi(:,1)
c      write(*,100) 'nv0   : ',nv0(:,1)
c      write(*,100) 'u1    : ',u1(:,1)
100   format(a,8e16.5)
          return
        case (vieilles_burning)
c
          vi(:,1) = burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * nv0(:,1)
          vi(:,2) = burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * nv0(:,2)
          vi(:,3) = burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * nv0(:,3)
c
        case (clausius_clapeyron)
c
          x_vapor = exp( - hfg_liquid / Ru * (one/T0 - one/T_boil_liquid) )
          mw_mix  = x_vapor*MW_liquid + (one-x_vapor)*mat_prop(mater0,iprop_ideal_gas_mw, 1)
c
          y_vapor = x_vapor*MW_liquid/MW_mix
          y_vapor = max(zero,y_vapor)
          y_vapor = min(one, y_vapor)
c
          vap_frac0 = y_vapor
c
          rho_mix = pres0 / (Ru/mw_mix*1.d3*T0)
c
          dydn = zero
c
          do n = 1,nshl0
            dydn = dydn + ycl0(:,n,ndof)*shg0(:,n,1)*nv1(:,1)
     &                  + ycl0(:,n,ndof)*shg0(:,n,2)*nv1(:,2)
     &                  + ycl0(:,n,ndof)*shg0(:,n,3)*nv1(:,3)
          enddo
c
          un0 = (u0(:,1)-um0(:,1))*nv0(:,1) + (u0(:,2)-um0(:,2))*nv0(:,2) + (u0(:,3)-um0(:,3))*nv0(:,3)
          un1 = (u1(:,1)-um1(:,1))*nv1(:,1) + (u1(:,2)-um1(:,2))*nv1(:,2) + (u1(:,3)-um1(:,3))*nv1(:,3)
c
          vap_rate = rho_mix/rho1*(-y_vapor*un0 - scdiff(1)*dydn)
c
          vi(:,1) = c1 * (vap_rate * nv0(:,1) + u1(:,1))
          vi(:,2) = c1 * (vap_rate * nv0(:,2) + u1(:,2))
          vi(:,3) = c1 * (vap_rate * nv0(:,3) + u1(:,3))
c          vi(:,1) = vap_rate * nv0(:,1)
c          vi(:,2) = vap_rate * nv0(:,2)
c          vi(:,3) = vap_rate * nv0(:,3)
c
#if debug == 1
      do iel = 1,npro
        write(*,10) iel, y_vapor(iel),vap_rate(iel),vi(iel,1)
      enddo
10    format('E3IF_VI: iel,y_vapor,vap_rate,vi(1):',i4,3e18.6)
#endif
C
c... NOTE: we don't add liquid velocity because it's already there...
c
      return
c
        case (cavitation)
c
c >>> NOTE: HARD CODED NUMBERS HERE
c
          pvap = 110.0d3
          cdest = 1.d-4
          cprod = 1.d-4
c
c          vi(:,1) = -nv0(:,1)*(cdest*min(zero,pres0(i)-pvap) + cprod*max(zero,pres0-pvap))
c          vi(:,2) = -nv0(:,2)*(cdest*min(zero,pres0(i)-pvap) + cprod*max(zero,pres0-pvap))
c          vi(:,3) = -nv0(:,3)*(cdest*min(zero,pres0(i)-pvap) + cprod*max(zero,pres0-pvap))
          vi(:,1) = nv1(:,1)*(cprod*max(zero,pres0-pvap))
          vi(:,2) = nv1(:,2)*(cprod*max(zero,pres0-pvap))
          vi(:,3) = nv1(:,3)*(cprod*max(zero,pres0-pvap))
c      vi_max = 10.0d0
c      vi(:,1) = nv1(:,1)*(min(vi_max,cprod*max(zero,pres0-pvap)))
c      vi(:,2) = nv1(:,2)*(min(vi_max,cprod*max(zero,pres0-pvap)))
c      vi(:,3) = nv1(:,3)*(min(vi_max,cprod*max(zero,pres0-pvap)))
c
        case default
          call error ('ERROR in e3if_vi:',' phase_change_model is not supported.',phase_change_model)
        end select
c
c... add flow velocity
c
        vi(:,1) = c1*(vi(:,1) + u(:,1))
        vi(:,2) = c1*(vi(:,2) + u(:,2))
        vi(:,3) = c1*(vi(:,3) + u(:,3))
c
      end subroutine calc_vi
c
      subroutine calc_vi_area_node(sum_vi_area_l,shp,WdetJif,nshl)
c
        real*8, dimension(:,:,:), intent(inout) :: sum_vi_area_l
        real*8, dimension(:,:),   intent(in)    :: shp
        real*8, dimension(:),     intent(in)    :: WdetJif
        integer, intent(in) :: nshl
c
        integer :: n,isd
c
        do n = 1,nshl
          do isd = 1,nsd
            sum_vi_area_l(:,n,isd) = sum_vi_area_l(:,n,isd) 
     &                             + shp(:,n)*vi(:,isd)*area(:)
c     &                             + shp(:,n)*vi(:,isd)*WdetJif(:)
          enddo
          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*area(:)
c          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*WdetJif(:)
        enddo
c
      end subroutine calc_vi_area_node
c
      subroutine calc_vapor_frac_node
        integer :: n
        do n = 1,nshl0
          ifbc_l0(:,n,ivapor_frac) = ifbc_l0(:,n,ivapor_frac) + shp0(:,n)*area*vap_frac0
          ifbc_l0(:,n,iwt)         = ifbc_l0(:,n,iwt)         + shp0(:,n)*area
          ifbc_l1(:,n,ivapor_frac) = ifbc_l1(:,n,ivapor_frac) + shp1(:,n)
          ifbc_l1(:,n,iwt)         = ifbc_l1(:,n,iwt)         + shp1(:,n)
        enddo
      end subroutine calc_vapor_frac_node
c
      end module e3if_vi_m
