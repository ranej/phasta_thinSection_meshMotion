      module block_m
c
c----------------------------------------
c    aims to query the data of each
c    interior block
c----------------------------------------
        use blkdat_m
        use propar_m
        use elmpar_m
        use intpt_m
        use shpdat_m
        use eqn_state_m
        use e3_param_m
        use matdat_def_m
        use pointer_data
        implicit none
c
        interface 
          subroutine set_block(iblk)
            integer, intent(in) :: iblk
          end subroutine
        end interface 
c
        integer, dimension(:,:), pointer :: ien
        procedure(set_block), pointer :: set_block_ptr
c
      contains
c
      subroutine set_interior_block(iblk)
        integer, intent(in) :: iblk
        npro   = lcblk(1,iblk+1) - lcblk(1,iblk)
        lcsyst = lcblk(3,iblk)
        mattyp = lcblk(7,iblk)
        nshl   = lcblk(10,iblk)
        ngauss = nint(lcsyst)
        call set_block_pointers(mattyp)
        ien => mien(iblk)%p
      end subroutine
c
      subroutine set_block_pointers(mattyp)
        use e3_solid_m
        integer, intent(in) :: mattyp
        e3_malloc_ptr => e3_malloc
        e3_mfree_ptr => e3_mfree
        select case (mat_eos(mattyp,1))
        case (ieos_ideal_gas,ieos_ideal_gas_2)
          getthm6_ptr => getthm6_ideal_gas
          getthm7_ptr => getthm7_ideal_gas
        case (ieos_liquid_1)
          getthm6_ptr => getthm6_liquid_1
          getthm7_ptr => getthm7_liquid_1
        case (ieos_solid_1)
          e3_malloc_ptr => e3_malloc_solid
          e3_mfree_ptr => e3_mfree_solid
          getthm6_ptr => getthm6_solid_1
          getthm7_ptr => getthm7_solid_1
        case default
          call error ('getthm  ', 'wrong material', mattyp)
        end select
      end subroutine set_block_pointers
c
      end module block_m
