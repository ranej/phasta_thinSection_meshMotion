c
c  rho                   : density
c  p                     : pressure
c  T                     : temperature
c  ei                    : internal energy
c  h                     : enthalpy
c  alphaP                : expansivity
c  betaT                 : isothermal compressibility
c  cp                    : specific heat at constant pressure
c  rk                    : kinetic energy
c
      module e3_param_m
c
c----------------------------------------
c    allocate/deallocate common parameters
c    used for DG interface only
c----------------------------------------
        use propar_m
        implicit none
        abstract interface
          subroutine getthm2
          end subroutine getthm2
        end interface
        integer :: nqpt
        integer :: mater0, mater1, mater
        real*8, dimension(:), pointer :: rho, pres, T, cp, rk
        real*8, dimension(:), pointer :: ei, p, h, cv
     &,                                  alphaP, betaT, gamb, c
     &,                                  vap_frac
        procedure(getthm2), pointer :: getthm6_ptr, getthm7_ptr
        procedure(e3_malloc), pointer :: e3_malloc_ptr
        procedure(e3_mfree), pointer :: e3_mfree_ptr
        contains 
        subroutine e3_malloc
          allocate(rho(npro))
          allocate(pres(npro))
          allocate(T(npro))
          allocate(cp(npro))
          allocate(rk(npro))
          allocate(ei(npro))
          allocate(p(npro))
          allocate(h(npro))
          allocate(cv(npro))
          allocate(alphap(npro))
          allocate(betaT(npro))
          allocate(gamb(npro))
          allocate(c(npro))
          allocate(vap_frac(npro))
        end subroutine e3_malloc
        subroutine e3_mfree
          deallocate(rho,pres,T,cp,rk)
          deallocate(ei,p,h,cv,alphap,betaT,gamb,c)
          deallocate(vap_frac)
        end subroutine e3_mfree
      end module e3_param_m
c
      module e3Sclr_param_m
c
        use e3_param_m
c
        implicit none
c
        real*8, dimension(:), pointer :: Sclr
c
        contains
c
        subroutine e3Sclr_malloc
          call e3_malloc
          allocate(Sclr(npro))
        end subroutine e3Sclr_malloc
c
        subroutine e3Sclr_mfree
          call e3_mfree
          deallocate(Sclr)
        end subroutine e3Sclr_mfree
c
      end module e3Sclr_param_m
