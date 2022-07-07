      module mattype_m

       use workfc_m
c
        implicit none
c
        integer, pointer :: neltp_mattype(:)
        integer, target, allocatable :: mattype(:), mattypeif(:,:)
        integer, pointer :: ientmp(:,:), ienif0tmp(:,:), ienif1tmp(:,:), ibcbtmp(:,:)
        real*8,  pointer :: bcbtmp(:,:)
c
      contains
c
      subroutine count_elem_mattype(mattype,neltp,mat_tag,nummat)
c
        integer, intent(in) :: nummat,neltp
        integer, intent(in) :: mattype(neltp)
        integer, intent(in) :: mat_tag(nummat)
c
        logical :: mat_added
        integer :: iel,imat, total_count
c
c count the number of elements with same material type
c
          neltp_mattype = 0
          total_count = 0
c
          do iel = 1, neltp
            do imat = 1,nummat
              if (mattype(iel) == mat_tag(imat)) then
                neltp_mattype(imat) = neltp_mattype(imat) + 1
                total_count = total_count + 1
              endif
            enddo
          enddo
c
          if (total_count /= neltp) then
            write(*,'(a,i4,a,4i6)') '[',myrank,'] ERROR: count_elem_mattypex. total_count, neltp, neltp_mattype(1:2):',
     &       total_count,neltp,neltp_mattype
            write(*,*) '[',myrank,'] mat_tag:',mat_tag
            write(*,*) '[',myrank,'] mattype:',mattype
            call error ('genblk  ', 'Could not count material type correctly', total_count)
          endif
c
      end subroutine count_elem_mattype
c
      subroutine get_npro_mattype(npro,mattype,ientp,isave,mat_tag,neltp,nshl,ibksz)
c
        integer, intent(out) :: npro
        integer, intent(inout) :: isave
        integer, intent(in)  :: mat_tag,neltp,nshl,ibksz
        integer, intent(in)  :: mattype(neltp),ientp(neltp,nshl)
c
        integer :: ipro
c
        npro = 0
        do 
          if (mattype(isave) == mat_tag) then
            npro = npro + 1
            ientmp(npro,1:nshl) = ientp(isave,1:nshl)
          endif
          isave = isave + 1
          if (npro == ibksz .or. isave>neltp ) exit
        enddo
c
      end subroutine get_npro_mattype
c
      end module mattype_m
