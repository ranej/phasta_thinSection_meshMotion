        subroutine e3if_setparam
     &  (
     &    nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_,
     &    npro_,ndof_,nsd_,nflow_,ipord_,nqpt_
     &  )
c
          use e3if_param_m
c
          implicit none
c
          integer, intent(in) :: nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_
          integer, intent(in) :: npro_,ndof_,nsd_,nflow_,ipord_,nqpt_
c
          nshg  = nshg_
          nshl0 = nshl0_
          nshl1 = nshl1_
          nenl0 = nenl0_
          nenl1 = nenl1_
          lcsyst0 = lcsyst0_
          lcsyst1 = lcsyst1_
          npro  = npro_
          ndof  = ndof_
          nflow = nflow_
          ipord = ipord_
          nqpt  = nqpt_
c
        end subroutine e3if_setparam
c
        subroutine e3if_setparam2
     &  (
     &    egmassif00,egmassif01,egmassif10,egmassif11,time_
     &  )
c
          use dgifinp_m
          use e3if_param_m
          use matdat_def_m
c
          implicit none
c
          real*8, dimension(:,:,:), allocatable, target, intent(in) :: egmassif00,egmassif01,egmassif10,egmassif11
          real*8, intent(in) :: time_
c
          time  = time_
c
          if (allocated(egmassif00)) then
            egmass00 => egmassif00
            egmass01 => egmassif01
            egmass10 => egmassif10
            egmass11 => egmassif11
          else
            nullify(egmass00)
            nullify(egmass01)
            nullify(egmass10)
            nullify(egmass11)
          endif
c
          s = dgif_s
          e = dgif_e
	  emu = dgif_emu
	  ek = dgif_ek
          length_h = dgif_h
c
          
c
          lhs_dg = lhs
c
        end subroutine e3if_setparam2
