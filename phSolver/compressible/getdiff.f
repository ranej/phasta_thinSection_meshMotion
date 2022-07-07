      subroutine getdiff(rmu, rlm, rlm2mu, con, npro, mater)
c
        use number_def_m
        use matdat_def_m
c
        implicit none
c
        integer, intent(in) :: npro, mater
        real*8, dimension(npro), intent(out) :: rmu,rlm,rlm2mu,con
c
        integer :: ivisc, icond
c
        select case (mat_eos(mater,1))
        case (ieos_ideal_gas,ieos_ideal_gas_2,ieos_ideal_gas_mixture)
          ivisc = iprop_gas_visc
          icond = iprop_gas_cond
        case (ieos_liquid_1)
          ivisc = iprop_liquid_1_visc
          icond = iprop_liquid_1_cond
        case (ieos_solid_1)
          icond = iprop_solid_1_cond
          rmu = zero
          rlm = zero
          rlm2mu = zero
          con = mat_prop(mater,icond,1)
          return
        case default
          call error ('getdiff ', 'ERROR: index can not be set!',0)
        end select
c
        if (matflg(2,mater) == 0) then ! Shear Law: Constant Viscosity
          rmu = mat_prop(mater,ivisc,1)
        else
          call error ('getdiff ', 'ERROR: Constant Viscosity is supported ONLY!', 0)
        endif
c
        if (matflg(3,mater) == 0) then ! Bulk Viscosity Law: Constant Bulk Viscosity
          rlm = -pt66 * rmu
        else
          call error ('getdiff ', 'ERROR: Constant Bulk Viscosity is supported ONLY!', 0)
        endif
c
        rlm2mu = rlm + two * rmu
c
        con = mat_prop(mater,icond,1)
c
c        if (iLES .gt. 0 .or. iRANS.eq.0) then
c          call error ('getdiff ', 'ERROR: Turbulence Viscosity is NOT supported!', 0)
c        endif
c
      end subroutine getdiff
c
      subroutine getDiffSclr(rmu, rlm, rlm2mu, con, npro, mater)
c
        use number_def_m
        use matdat_def_m
        use sclrs_m
c
        implicit none
c
        integer, intent(in) :: npro, mater
        real*8, dimension(npro), intent(out) :: rmu,rlm,rlm2mu,con
c
c
        rmu = scdiff(isclr)
c
        rlm = -pt66 * rmu
        rlm2mu = rlm + two * rmu
c
      end subroutine getdiffSclr
