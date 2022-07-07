      module asadj_m
c
c----------------------------------------
c    aims to modify the sparse matrix data
c    structure for interface
c----------------------------------------
       implicit none
c
      contains
c
        subroutine Asadj_if(row_fill_list, adjcnt, ienif0, ienif1,
     &                      npro,nshl0,nshl1,nshg,nnz)
c
          integer, intent(in) :: npro,nshl0,nshl1,nshg,nnz
          integer, dimension(nshg, 15*nnz), intent(inout) :: row_fill_list
          integer, dimension(nshg),         intent(inout) :: adjcnt
          integer, dimension(npro,nshl0), intent(in) :: ienif0
          integer, dimension(npro,nshl1), intent(in) :: ienif1
c
          integer :: i,j
          integer :: ndlist(nshl0+nshl1)
c
          do i = 1, npro
c
c... add nodes to the list
c
            do j = 1, nshl0
              ndlist(j) = ienif0(i, j)
            enddo
            do j = 1, nshl1
              ndlist(j+nshl0) = ienif1(i, j)
            enddo
      if (any(ndlist(:).gt.nshg))then
        write(*,*) 'WARNING: nno > nsgh!!!'
      endif
c
            do j = 1, nshl0+nshl1
              call set_this_row
     &            (j,nshl0+nshl1,nshg,15*nnz,ndlist,adjcnt,row_fill_list(ndlist(j),:))
            enddo
c
          enddo
c
        end subroutine Asadj_if
c
        subroutine set_this_row(j,n,nshg,maxlen,ndlist,adjcnt,this_row_fill_list)
c
          integer, dimension(maxlen), intent(inout) :: this_row_fill_list
          integer, dimension(nshg), intent(inout) :: adjcnt
          integer, intent(in) :: j, n, ndlist(n), nshg, maxlen
c
          integer :: jnd, jlngth, knd, k, ibroke, l
c
          jnd    = ndlist(j)
          jlngth = adjcnt(jnd)
c
          do k = 1, n
            knd = ndlist(k)
            ibroke = 0
            do l = 1, jlngth
              if (this_row_fill_list(l).eq.knd) then
                ibroke = 1
                exit
              endif
            enddo
c
            if (ibroke .eq. 0) then
              jlngth = jlngth + 1
              if (jlngth .gt. maxlen) then
                 write(*,*) 'increase overflow factor in genadj'
                 stop
              endif
              this_row_fill_list(jlngth)=knd ! add unique entry to list
            endif
c
            adjcnt(jnd) = jlngth
c
          enddo
c
        end subroutine set_this_row
c
      end module asadj_m
