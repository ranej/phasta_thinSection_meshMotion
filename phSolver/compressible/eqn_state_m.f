      module eqn_state_m
c
c----------------------------------------
c    aims to compute the equation of
c    state
c----------------------------------------
        use e3_param_m
        use number_def_m
        use matdat_def_m
        use mmatpar_m
        use dgifinp_m
c
        implicit none
c
        real*8 :: rho_ref, p_ref, T_ref, alpha_P, beta_T, cv_liq
        real*8 :: rho_ref_s, p_ref_s, T_ref_s, alpha_P_s, cv_s
        real*8 :: bulkMod_s, shearMod_s 
c
      contains
c
      subroutine getthm6_ideal_gas
c
        mw    = mat_prop(mater,iprop_ideal_gas_mw, 1)
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = Ru/mw*1.0d3
        gamma1 = gamma - one
c
        rho = pres / (Rgas*T)
        ei  = T * Rgas / gamma1
c
      end subroutine getthm6_ideal_gas
c
      subroutine getthm7_ideal_gas
c
        call getthm6_ideal_gas
c
        h   = T * Rgas / gamma1 * gamma
        cp  = Rgas*gamma / gamma1
        alphaP = one / T
        betaT  = one / pres
        if (associated(cv)) cv  = Rgas / gamma1
        if (associated(gamb)) gamb = gamma1
        if (associated(c)) c =  sqrt( gamma * Rgas * T )
c
      end subroutine getthm7_ideal_gas
c
      subroutine getthm6_ideal_gas_mixture
c
        integer :: iel
        real*8, dimension(npro) :: mw,rgas,y
c
        y = vap_frac
        y = max(zero,y)
        y = min(one, y)
c
        mw  = y*MW_liquid + (one-y)*mat_prop(mater,iprop_ideal_gas_mw, 1)
c
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = Ru/mw*1.0d3
        gamma1 = gamma - one
c
        rho = pres / (Rgas*T)
        ei  = T * Rgas / gamma1
c
      end subroutine getthm6_ideal_gas_mixture
c
      subroutine getthm7_ideal_gas_mixture
c
        real*8, dimension(npro) :: mw,rgas,y,gamma_mix,cp_gas,cv_gas,cv_mix,cp_mix
c
        y = vap_frac
        y = max(zero,y)
        y = min(one, y)
c
        mw  = y*MW_liquid + (one-y)*mat_prop(mater,iprop_ideal_gas_mw, 1)
c
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        gamma1 = gamma - one
        Rgas  = Ru/mw*1.0d3
        rho = pres / (Rgas*T)
        cp_gas  = Rgas*gamma / gamma1
        cv_gas  = Rgas / gamma1
        cv_liq  = mat_prop(2,iprop_liquid_1_cv,     1)
c
        cp_mix = y*cv_liq + (one-y)*cp_gas
        cv_mix = y*cv_liq + (one-y)*cv_gas
c
        ei  = T * cv_mix
        h   = T * cp_mix
c
        alphaP = one / T
        betaT  = one / pres
c
c        gamb = gamma1
        if (associated(cv))   cv = cv_mix
        if (associated(gamb)) gamb = cp_mix / cv_mix - one
        if (associated(C))    c =  sqrt( gamma * Rgas * T )
c
      end subroutine getthm7_ideal_gas_mixture
c
      function rho_ideal_gas(p,R,T) result(rho)
        implicit none
        real*8 p,R,T,rho
        rho = p / (R*T)
      end function rho_ideal_gas
c
      subroutine getthm6_liquid_1
c
        rho_ref = mat_prop(mater,iprop_liquid_1_rho_ref,1)
        p_ref   = mat_prop(mater,iprop_liquid_1_p_ref,  1)
        T_ref   = mat_prop(mater,iprop_liquid_1_T_ref,  1)
        cv_liq  = mat_prop(mater,iprop_liquid_1_cv,     1)
        alpha_P = mat_prop(mater,iprop_liquid_1_alphaP, 1)
        beta_T  = mat_prop(mater,iprop_liquid_1_betaT,  1)
c
        rho = rho_ref * (one - alpha_P*(T-T_ref) + beta_T*(pres-P_ref))
        ei  = cv_liq*T + if_reaction_heat
c
      end subroutine getthm6_liquid_1
c
      subroutine getthm7_liquid_1
c
        call getthm6_liquid_1
c
        h   = ei + pres/rho
        cp  = cv_liq
        alphaP = alpha_P
        betaT  = beta_T
        if (associated(cv)) cv  = cv_liq
        if (associated(gamb)) gamb = zero
c        c =  sqrt(one/(rho_ref*betaT))
        if (associated(c)) c =  sqrt(one/(rho*betaT))
c
      end subroutine getthm7_liquid_1
c
      subroutine getthm6_solid_1
c
        use e3_solid_m
        implicit none
c
        rho_ref_s = mat_prop(mater,iprop_solid_1_rho_ref,1)
        p_ref_s   = mat_prop(mater,iprop_solid_1_p_ref,  1)
        T_ref_s   = mat_prop(mater,iprop_solid_1_T_ref,  1)
        cv_s     = mat_prop(mater,iprop_solid_1_cv,     1)
        alpha_P_s = mat_prop(mater,iprop_solid_1_alphaP, 1)
        bulkMod_s = mat_prop(mater,iprop_solid_1_bulkMod, 1)
        shearMod_s = mat_prop(mater,iprop_solid_1_shearMod, 1)
c        beta_T  = mat_prop(mater,iprop_solid_1_betaT,  1)
c
        alphaP = alpha_P_s
        betaT  = one /(bulkMod_s * Ja_def) ! double check here
        rho = rho_ref_s * (one - alphaP*(T - T_ref_s) 
     &                    + betaT*(pres - p_ref_s))
c        ei  = ( cv_s - pres * alpha_P_s/rho)* T 
c    &        + (betaT * pres - alphaP * T)/rho * pres
        ei  = cv_s * T
c
      end subroutine getthm6_solid_1
c
      subroutine getthm7_solid_1
c
        use e3_solid_m
        implicit none
c
        call getthm6_solid_1
c
        h   = ei + pres/rho
        cv  = cv_s
        cp  = cv_s
        bulkMod = bulkMod_s
        shearMod = shearMod_s
c        c =  sqrt(one/(rho_ref*betaT))
C         c =  sqrt(one/(rho*betaT))
C         gamb = zero
c
c
      end subroutine getthm7_solid_1
c
      subroutine getthmif0
        use e3if_param_m
        call e3if_setparam_0
        call getthmif0_ptr
      end subroutine getthmif0
c
      subroutine getthmif1
        use e3if_param_m
        call e3if_setparam_1
        call getthmif1_ptr
      end subroutine getthmif1
c
      subroutine getthmif_solid_0
c
      end subroutine getthmif_solid_0
c
      subroutine getthmif_solid_1
c
      end subroutine getthmif_solid_1
c
      end module eqn_state_m
