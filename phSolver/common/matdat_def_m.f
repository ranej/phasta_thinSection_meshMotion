      module matdat_def_m
c
        use global_const_m
c
        implicit none
c
        integer :: nummat
     &,            matflg(MAXMAT,MAXTS)
     &,            mat_tag(MAXMAT,MAXTS)
     &,            mat_eos(MAXMAT,MAXTS)
        integer, parameter :: ieos_ideal_gas = 1
     &,                       ieos_ideal_gas_2 = 2
     &,                       ieos_ideal_gas_mixture = 3
     &,                       ieos_liquid_1  = 4
     &,                       ieos_solid_1   = 10
c
        integer, parameter :: iprop_ideal_gas_mw     = 1  ! molecular weight
     &,                       iprop_ideal_gas_gamma  = 2
     &,                       iprop_gas_visc         = 3
     &,                       iprop_gas_cond         = 4
     &,                       iprop_liquid_1_rho_ref = 1
     &,                       iprop_liquid_1_p_ref   = 2
     &,                       iprop_liquid_1_T_ref   = 3
     &,                       iprop_liquid_1_cv      = 4
     &,                       iprop_liquid_1_alphaP  = 5
     &,                       iprop_liquid_1_betaT   = 6
     &,                       iprop_liquid_1_visc    = 7
     &,                       iprop_liquid_1_cond    = 8
c
c......Adding the solid capability    
     &,                       iprop_solid_1_rho_ref  = 1
     &,                       iprop_solid_1_p_ref    = 2
     &,                       iprop_solid_1_T_ref    = 3
     &,                       iprop_solid_1_cv       = 4
     &,                       iprop_solid_1_alphaP   = 5
     &,                       iprop_solid_1_bulkMod  = 6
     &,                       iprop_solid_1_shearMod = 7
     &,                       iprop_solid_1_cond     = 8
C      &,                       iprop_solid_1_shearMod = 7
C      &,                       iprop_solid_1_bulkMod  = 8

c
        integer :: surface_tension_flag
        integer :: datelas_volume_YM
c
        real*8  :: datmat(MAXMAT,MAXPROP,MAXTS)
        real*8  :: mat_prop(MAXMAT,MAXPROP,MAXTS)
        real*8  :: surface_tension_coeff
        real*8  :: datelas(1,2)
c
c NOTE: This should only be an integer array BUT
c       since vtk reader does not support int yet,
c       so, we cast the integer numbers into this 
c       real array!!!
        real*8, pointer :: mattype_interior(:)
c
        logical :: mexist
c
        common /matdat/ matflg, mat_tag, mat_eos, nummat, mexist
     &,   mat_prop, datmat, datelas
     &,   surface_tension_coeff, surface_tension_flag
     &,   datelas_volume_YM
c
c
      end module matdat_def_m
