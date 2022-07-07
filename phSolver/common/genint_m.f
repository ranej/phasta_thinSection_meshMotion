      module genint_m
c
c----------------------------------------
c    aims to generate the quadrature point
c    specially for DG interface elements
c----------------------------------------
        integer, dimension(4), parameter :: nint_tet (4) = (/1,4,16,29/)
        integer, dimension(5), parameter :: nint_quad(5) = (/1,4,9,16,25/)
c
        integer, dimension(4), parameter :: nint_tri (4) = (/1,3,6,12/)
c
        contains
c
        subroutine genint_if
c
          use number_def_m
          use intpt_m
          use intdat_m
c
          implicit none
c
          integer :: itp,ierr
          real*8, allocatable :: tmpQptif(:,:), tmpQwtif(:)
c
c ... Tet
c
          nintif = nint_tri(intg(3,1))
c
          allocate(tmpQptif(4,nintif(nint_tri(intg(3,1)))))
          allocate(tmpQwtif(nintif(nint_tri(intg(3,1)))))
c
          itp = itp_tet
c
          call symtri(nintif(itp),tmpQptif,tmpQwtif,ierr) 
c
          Qptif (itp,1:4,1:nintif(itp)) = tmpQptif(1:4,1:nintif(itp))
          Qwtif (itp,1:nintif(itp))     = tmpQwtif(1:nintif(itp))
c
c.... adjust quadrature weights to be consistent with the
c     design of tau. 
c
          Qwtif(itp,:)  = two*Qwtif(itp,:)
c
c ... Wedge
c
          itp = itp_wedge
c
          call symtri(nintif(itp),tmpQptif,tmpQwtif,ierr) 
c
          Qptif(itp,1:2,2:nintif(3)) = tmpQptif(1:2,1:nintif(3)-1)
          Qptif(itp,1:2,1) = tmpQptif(1:2,nintif(3))
c     
c     wedges want the third entry to be zeta=-1 (not t=1-r-s)
c     4th entry not used
c     
          Qptif(itp,3:4,1:nintif(3)) = -1
c$$$       Qptb(3,1:4,1:nintb(3)) = tmpQptb(1:4,1:nintb(3))
          Qwtif(itp,1:nintif(3)) = tmpQwtif(1:nintif(3))
c
          deallocate(tmpQptif,tmpQwtif)
c
        end subroutine genint_if
c
      end module genint_m
