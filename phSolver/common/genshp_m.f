      module genshp_m
c
c----------------------------------------
c    aims to generate the shape function
c    specially for DG interface elements
c----------------------------------------
        use global_const_m
        use number_def_m
        use elmpar_m
        use shpdat_m
        use genpar_m
        use blkdat_m
        use intpt_m
c
        implicit none
c
        contains

        subroutine genshpif
     &(            shpif
     &,            shgif
     &)
c
          real*8, dimension(MAXTOP,      maxsh, MAXQPT), intent(out) :: shpif   ! shape function
          real*8, dimension(MAXTOP, nsd, maxsh, MAXQPT), intent(out) :: shgif   ! shape function gradient 
c
          integer :: i, iblk, id, itpif
c
c... loop over the interface blocks and 
c    generate shape functions and gradients
c
          do iblk = 1, nelblif
c
            nshl0   = lcblkif(iblkif_nshl0,iblk)
            nshl1   = lcblkif(iblkif_nshl1,iblk)
            id      = lcblkif(iblkif_topology,iblk)
c
            do i = 1, nintif(id)
              call shpTet(ipord, Qptif(itp_tet,1:3,i), shpif(itp_tet,:,i), shgif(itp_tet,:,:,i))
              call shp6w (ipord, Qptif(itp_wedge,1:3,i), shpif(itp_wedge,:,i),shgif(itp_wedge,:,:,i))
            enddo
c
            shgif(itp_tet,:,:,1:nintif(id)) = 
     &        shgif(itp_tet,:,:,1:nintif(id)) / two 
c
          enddo
c
        end subroutine genshpif
c
      end module genshp_m
