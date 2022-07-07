      subroutine genbkif (ibksz)
c
c----------------------------------------
c    aims to read data from geombc and
c    generate the interface blocks
c----------------------------------------
        use pointer_data
        use phio
        use iso_c_binding
        use mattype_m
c
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk
c
        integer, intent(in) :: ibksz
c
        integer :: iel, iblk, itpblk,
     &             n,neltp,ioprdl,nnface,
     &             iientpsiz0,iientpsiz1
        integer, target :: intfromfile(50)
        integer, target, allocatable :: ientp(:,:)
        integer, dimension(ibksz) :: mater
        integer :: descriptor, descriptorG, GPID, color
        integer ::  numparts, writeLock
        integer :: ierr_io, numprocs
        integer, target :: itpblktot,ierr
        integer, parameter :: ione = 1, itwo = 2, iseven = 7, inine = 9
        character*255 :: fname1, fname2, tempchar
        character(len=30) :: dataInt
        integer :: tmpnenl0,tmpnenl1,tmpnshl0,tmpnshl1,tmplcsyst0,tmplcsyst1,iptr
        dataInt = c_char_'integer'//c_null_char
c
c
        iel = 1              ! element index
        itpblk = nelblif     ! number of topological blocks
        nelblif = 0         ! set a counter for the interface element blocks (size of ibksz)
c
        ! Get the total number of different interface topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
C
        itpblktot = 0
        call phio_readheader(fhandle,
     &   c_char_'number of interface tpblocks' // char(0),
     &   c_loc(itpblktot), ione, dataInt, iotype) 
c
        if (itpblktot == -1) then 
          ! The field 'total number of different interface tpblocks' was not found in the geombc file.
          ! Scan all the geombc file for the 'connectivity interface' fields to get this information.
          iblk=0
          neltp=0
          do while(neltp .ne. -1) 

            ! intfromfile is reinitialized to -1 every time.
            ! If connectivity interior@xxx is not found, then 
            ! readheader will return intfromfile unchanged

            intfromfile(:)=-1
            iblk = iblk+1
            if(input_mode.ge.1) then
              write (fname2,"('connectivity interface',i1)") iblk
            else
              write (fname2,"('connectivity interface linear tetrahedron tetrahedron')") 
!              write (fname2,"('connectivity interface?')") 
            endif

            !write(*,*) 'rank, fname2',myrank, trim(adjustl(fname2))
            call phio_readheader(fhandle, fname2 // char(0),
     &       c_loc(intfromfile), iseven, dataInt, iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of interface topologies: ',itpblktot
        endif

        mattyp0 = 0
        mattyp1 = 1
        ndofl  = ndof
        nsymdl = nsymdf      ! ????
c
c        iblk_loop: do iblk = 1, maxtopif
        iblk_loop: do iblk = 1, itpblktot
c
           writeLock=0;
           if(input_mode.ge.1) then
             write (fname2,"('connectivity interface',i1)") iblk
           else
c             select case (iblk)
c             case (itpif_tet_tet)
c               write (fname2,"('connectivity interface linear tetrahedron tetrahedron')") 
c             case (itpif_wedge_tet)
c               write (fname2,"('connectivity interface linear wedge tetrahedron')") 
c             case (itpif_tet_wedge)
c               write (fname2,"('connectivity interface linear tetrahedron wedge')") 
c             case (itpif_wedge_wedge)
c               write (fname2,"('connectivity interface linear wedge wedge')") 
c             case default
c               cycle iblk_loop
c             end select
      write (fname2,"('connectivity interface?')") 
           endif
c
           ! Synchronization for performance monitoring, as some parts do not include some topologies
c           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           intfromfile(:)=-1
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), inine, dataInt, iotype)
c
           neltp  = intfromfile(1)       ! number of pair of elements in this block
           tmpnenl0  = intfromfile(2)       ! number of nodes in element 0
           tmpnenl1  = intfromfile(3)       ! number of nodes in element 1
           ipordl = intfromfile(4)       ! polynomial order
           tmpnshl0  = intfromfile(5)
           tmpnshl1  = intfromfile(6)
           nnface = intfromfile(7)       ! number of nodes on the interface
           tmplcsyst0= intfromfile(8)       ! element type 0
           tmplcsyst1= intfromfile(9)       ! element type 1
c      write(*,*) 'BKIF: iblk, intfromfile: ',iblk,intfromfile(1:9)
c
           if (neltp<0) then
              writeLock=1;
              cycle iblk_loop
           endif
c
c
c... reads all the connectivity data in one array
c
           iientpsiz = neltp*(tmpnshl0+tmpnshl1)
           allocate (ientp(neltp,(tmpnshl0+tmpnshl1)))
           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(ientp), iientpsiz, dataInt, iotype)
c
          if (nummat <= 1) then
            write(*,'(a,i4,a)') '[',myrank,'] ERROR (genbkif): Number of materials is <= 1'
            call error ('genbkif  ', 'Number of Materials', nummat)
          endif
c
           if(input_mode.ge.1) then
             write(fname2,"('material type interface',i1)") iblk
           else
c             select case (iblk)
c             case (itpif_tet_tet)
c               write(fname2,"('material type interface linear tetrahedron tetrahedron')")
c             case (itpif_wedge_tet)
c               write(fname2,"('material type interface linear wedge tetrahedron')")
c             case (itpif_tet_wedge)
c               write(fname2,"('material type interface linear tetrahedron wedge')")
c             case (itpif_wedge_wedge)
c               write(fname2,"('material type interface linear wedge wedge')")
c             end select
      write(fname2,"('material type interface?')")
           endif
c
c           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           intfromfile(:)=-1
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), itwo, dataInt, iotype)
           allocate(mattypeif(intfromfile(1),intfromfile(2)))
           imattypesiz = intfromfile(1)*intfromfile(2)
           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(mattypeif), imattypesiz, dataInt, iotype)
c
c ... material tag first always takes mattype0 and ienif0
c
          mattype0 = 1
          mattype1 = 2
c
          if     (mattypeif(1,1) == mat_tag(1,1)) then
            nshl0 = tmpnshl0
            nshl1 = tmpnshl1
            lcsyst0 = tmplcsyst0
            lcsyst1 = tmplcsyst1
            nenl0 = tmpnenl0
            nenl1 = tmpnenl1
          elseif (mattypeif(1,2) == mat_tag(1,1)) then
            nshl0 = tmpnshl1
            nshl1 = tmpnshl0
            lcsyst0 = tmplcsyst1
            lcsyst1 = tmplcsyst0
            nenl0 = tmpnenl1
            nenl1 = tmpnenl0
          else
            write(*,*) 'ERROR: in genbkif. material type is wrong!'
            call error ('genbkif  ', '', 0)
          endif
c
          allocate(ienif0tmp(neltp,nshl0))
          allocate(ienif1tmp(neltp,nshl1))
c
          if(writeLock==0) then
c
c... make blocks of elements
c
          iptr = 1
c
          blocks_loop:   do n = 1, neltp, ibksz
c
            npro = 0
c
            do
              npro = npro + 1
              if     (mattypeif(iptr,1) == mat_tag(1,1)) then
                ienif0tmp(npro,1:nshl0) = ientp(iptr,1:nshl0)
                ienif1tmp(npro,1:nshl1) = ientp(iptr,nshl0+1:nshl0+nshl1)
              elseif (mattypeif(iptr,2) == mat_tag(1,1)) then
                ienif0tmp(npro,1:nshl0) = ientp(iptr,nshl1+1:nshl1+nshl0)
                ienif1tmp(npro,1:nshl1) = ientp(iptr,1:nshl1)
              endif
              iptr = iptr + 1
              if (npro == ibksz .or. iptr > neltp) exit
            enddo
c
c      write(*,11) iel, lcsyst0, lcsyst1,nshl0,nshl1, mattype0, mattype1
11    format(' iel, lcsyst0, lcsyst1, nshl0, nshl1, mattype0, mattype1: ',7i4)
c
            nelblif = nelblif + 1
            lcblkif(1,nelblif) = iel
c            lcblkif(2,nelblif) = iopen   ! ??? see genblk.f
            lcblkif(3,nelblif) = lcsyst0  ! local coordinate system
            lcblkif(4,nelblif) = lcsyst1
            lcblkif(5,nelblif) = ipordl
            lcblkif(6,nelblif) = nenl0
            lcblkif(7,nelblif) = nenl1
            lcblkif(8,nelblif) = nnface
            lcblkif(9,nelblif) = mattype0
            lcblkif(10,nelblif) = mattype1
            lcblkif(11,nelblif) = ndof
            lcblkif(12,nelblif) = nsymdl ! ????
            lcblkif(iblkif_nshl0,nelblif) = nshl0
            lcblkif(iblkif_nshl1,nelblif) = nshl1
            lcblkif(iblkif_topology,nelblif) = iftp_map(lcsyst0,lcsyst1)
c      write(*,*) 'IF TOPOLOGY: ',iftp_map(lcsyst0,lcsyst1)
c
c... allocate memory for stack arrays
c
            allocate(mienif0(nelblif)%p(npro,nshl0))
            allocate(mienif1(nelblif)%p(npro,nshl1))
c
c... save the interface elements block
c
            call gensavif(mienif0(nelblif)%p,mienif1(nelblif)%p)
c
            iel = iel + npro
c
          enddo blocks_loop
          endif
c
          deallocate(ientp)
          deallocate(mattypeif)
          deallocate(ienif0tmp,ienif1tmp)
c
        enddo iblk_loop
c
        lcblkif(1,nelblif+1) = iel
c
        return
c
      end subroutine genbkif
