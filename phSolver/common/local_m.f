      module local_m
c
c----------------------------------------
c    aims to perform a vector gather/scatter
c    operation for DG interface
c----------------------------------------
        use mtimer2_m
        implicit none
c
      contains
c
        subroutine local
     & (
     &    global, rlocal, ientmp, n, code,
     &    nshg,   nshl,  npro,   ipord
     & )

          real*8, dimension(nshg,n),          intent(inout) :: global
          real*8,  dimension(:,:,:), pointer, intent(inout) :: rlocal
          integer, dimension(:,:), pointer,   intent(in)    :: ientmp
          integer,                            intent(in)    :: n
          character*8,                        intent(in)    :: code
          integer,                            intent(in)    :: nshg, nshl, npro, ipord
c
          integer :: i,j,nel
          integer, dimension(npro,nshl) :: ien
c
c... correct the sign for the higher orders
c
          if (ipord > 2) then
            ien = abs(ientmp)
          else
            ien = ientmp
          endif
c
c.... ------------------------->  'assembling '  <----------------------
c
          if (code .eq. 'scatter ') then
c
c.... scatter the data (possible collisions)
c
            do j = 1, n
              do i = 1, nshl
                do nel = 1,npro
                  global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                                 + rlocal(nel,i,j)
                enddo
              enddo
            enddo
c
c.... transfer and flop counts
c
            sbytes = sbytes + n*nshl*npro
            flops = flops + n*nshl*npro
c
            return
          endif
c
c.... --------------------------->  error  <---------------------------
c
          call error ('local   ', code, 0)
c
        end subroutine local
c
        subroutine localx
     &  (
     &    global, rlocal, ien,  nsd, code,
     &    nshg,   nenl,   npro
     &  )
          real*8, dimension(nshg,nsd),        intent(in)    :: global
          real*8,  dimension(:,:,:), pointer, intent(inout) :: rlocal
          integer, dimension(:,:), pointer,   intent(in)    :: ien
          integer,                            intent(in)    :: nsd
          character*8,                        intent(in)    :: code
          integer,                            intent(in)    :: nshg, nenl, npro
c
          integer :: i,j
c
c... ----------------------------> 'localization' <-----------------------
c
          if (code .eq. 'gather  ') then
c
c... gather the data
c
            do j = 1, nsd
              do i = 1, nenl
                rlocal(:,i,j) = global(ien(:,i),j)
              enddo
            enddo
c
c... transfer count
c
            gbytes = gbytes + nsd*nenl*npro
c
c... return
c
            return
          endif
c
c... ---------------------------> 'error' <-------------------------
c
          call error ('local   ', code, 0)
c
        end subroutine localx
c
        subroutine localy
     & (
     &    global, rlocal, ientmp, n, code,
     &    nshg,   nshl,  npro,   ipord
     & )

          real*8, dimension(nshg,n),          intent(in)    :: global
          real*8,  dimension(:,:,:), pointer, intent(inout) :: rlocal
          integer, dimension(:,:), pointer,   intent(in)    :: ientmp
          integer,                            intent(in)    :: n
          character*8,                        intent(in)    :: code
          integer,                            intent(in)    :: nshg, nshl, npro, ipord
c
          integer :: ishp,ivar
          integer, dimension(npro,nshl) :: ien
c
c... correct the sign for the higher orders
c
          if (ipord > 2) then
            ien = abs(ientmp)
          else
            ien = ientmp
          endif
c
c... -----------------------> localization <-------------------
c
          if (code .eq. 'gather  ') then
c
CAD      rlocal = yl={P, u, v, w, T, scalar1, ...}
CAD	 global = y = {u, v, w, P, T, scalar1, ...}
c
CAD      Put Pressure in the first slot of yl
CAD      Put u,v,w in the slots 2,3,4 of yl 
c
            do ishp = 1,nshl
              rlocal(:,ishp,1)   = global(ien(:,ishp),4)
              rlocal(:,ishp,2:4) = global(ien(:,ishp),1:3)
            enddo
c
CAD      Fill in the remaining slots with T, and additional scalars
c
            if (n > 4) then
              do ivar = 5,n
                do ishp = 1,nshl
                  rlocal(:,ishp,ivar) = global(ien(:,ishp),ivar)
                enddo
              enddo
            endif
c
c.... transfer count
c
          gbytes = gbytes + n*nshl*npro
c
c.... return
c
c          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
           write(*,*) 'do not use localy here'
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
           write(*,*) 'do not use localy here'
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end subroutine localy
c
      end module local_m
