      module weighted_normal_data_m
c-------------------------------------------------------------------------------
c  testing the weighted normal for the flow solver
c-------------------------------------------------------------------------------
        implicit none
c... define the weighted normal for 0 and 1 phases
        real*8, dimension(:,:), allocatable :: w_normal_global
        real*8, dimension(:), allocatable :: length_temp
        integer, dimension(:), pointer :: calc_factor_temp
        real*8, dimension(:,:,:), pointer :: w_normal_l0, w_normal_l1 
c                        
      end module weighted_normal_data_m
