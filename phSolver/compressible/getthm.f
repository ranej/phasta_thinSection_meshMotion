      subroutine getthm (rho_,ei_,p_,T_,npro_,mater_
     &,                  h_,  cv_,cp_,alphaP_,betaT_,gamb_,c_)
c
        use eqn_state_m
        use e3_solid_m
c
        implicit none
c
        real*8, dimension(npro_), target, intent(out) :: rho_,ei_,h_,cv_,cp_,alphaP_,betaT_,gamb_,c_
        real*8, dimension(npro_), target, intent(in) :: p_,T_
        integer, intent(in) :: npro_, mater_
c
        mater = mater_
c
        rho => rho_
        ei  => ei_
        p   => p_
        T   => T_
        h   => h_
        cv  => cv_
        cp  => cp_
        alphaP => alphaP_
        betaT  => betaT_
        gamb  => gamb_
        c => c_
c 
        select case (mat_eos(mater,1))
        case (ieos_ideal_gas,ieos_ideal_gas_2)
c
c          call getthm_ideal_gas
c
        case (ieos_liquid_1)
c
c          call getthm_liquid_1
c
        case (ieos_solid_1)
c
c          call getthm_solid_1
c
        case default
          call error ('getthm  ', 'wrong material', mater)
        end select
c
      end subroutine getthm
