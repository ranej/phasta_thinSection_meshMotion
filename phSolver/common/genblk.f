        subroutine genblk (IBKSZ)
c
c----------------------------------------------------------------------
c
c  This routine reads the interior elements and generates the
c  appropriate blocks.
c
c Zdenek Johan, Fall 1991.
c----------------------------------------------------------------------
c
        use pointer_data
        use phio
        use iso_c_binding
        use mattype_m
        use matdat_def_m
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk

        integer, target, allocatable :: ientp(:,:)
        integer, pointer :: ieMaptmp(:)
        integer mater(ibksz)
        integer, target :: intfromfile(50) ! integers read from headers
        character*255 :: fname1
        integer :: descriptor, descriptorG, GPID, color
        integer ::  numparts, writeLock
        integer :: ierr_io, numprocs
        integer, target :: itpblktot,ierr,iseven
        integer :: numel_ct
        character*255 fname2
        character(len=30) :: dataInt
        dataInt = c_char_'integer'//c_null_char
        numparts = numpe !This is the common settings. Beware if you try to compute several parts per process
        ione=1
        itwo=2
        iseven=7
        ieleven=11
        iel=1
        itpblk=nelblk

        ! Get the total number of different interior topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
C
C        itpblktot=1  ! hardwired to montopology for now
        call phio_readheader(fhandle,
     &   c_char_'number of interior tpblocks' // char(0),
     &   c_loc(itpblktot), ione, dataInt, iotype) 
C
        if (itpblktot == -1) then 
          ! The field 'total number of different interior tpblocks' was not found in the geombc file.
          ! Scan all the geombc file for the 'connectivity interior' fields to get this information.
          iblk=0
          neltp=0
          do while(neltp .ne. -1) 

            ! intfromfile is reinitialized to -1 every time.
            ! If connectivity interior@xxx is not found, then 
            ! readheader will return intfromfile unchanged

            intfromfile(:)=-1
            iblk = iblk+1
            if(input_mode.ge.1) then
              write (fname2,"('connectivity interior',i1)") iblk
            else
              write (fname2,"('connectivity interior?')") 
            endif

            !write(*,*) 'rank, fname2',myrank, trim(adjustl(fname2))
            call phio_readheader(fhandle, fname2 // char(0),
     &       c_loc(intfromfile), iseven, dataInt, iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of interior topologies: ',itpblktot
        endif

        nelblk=0
        mattyp = 0
        ndofl = ndof
        nsymdl = nsymdf
        numel_ct = 0
        allocate(mattype_interior(numel))

        iblk_loop: do iblk = 1, itpblktot
c        iblk_loop: do iblk = 1, maxtop
c
           writeLock=0;
c
           if(input_mode.ge.1) then
             write (fname2,"('connectivity interior',i1)") iblk
           else
c             select case (iblk)
c             case (itp_tet)
c               write (fname2,"('connectivity interior linear tetrahedron')") 
c            case (itp_wedge)
c              write (fname2,"('connectivity interior linear wedge')") 
c            case default
c               cycle iblk_loop
c             end select
      write (fname2,"('connectivity interior?')") 
           endif

           ! Synchronization for performance monitoring, as some parts do not include some topologies
C           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           intfromfile(:)=-1
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), iseven, dataInt, iotype)
           neltp  =intfromfile(1)
           nenl   =intfromfile(2)
           ipordl =intfromfile(3)
           nshl   =intfromfile(4)
           ijunk  =intfromfile(5)
           ijunk  =intfromfile(6)
           lcsyst =intfromfile(7)
c      write(*,*) ' BLK: iblk, intfromfile:',iblk,intfromfile(1:7) 
C
           if (neltp<0) then
              writeLock=1;
              cycle iblk_loop
           endif
C
           allocate (ientp(neltp,nshl))
           allocate (ientmp(ibksz,nshl))
           allocate (ieMaptmp(ibksz))
           allocate (neltp_mattype(nummat))
           iientpsiz=neltp*nshl

           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(ientp), iientpsiz, dataInt, iotype)
c
           if(input_mode.ge.1) then
             write(fname2,"('material type interior',i1)") iblk
           else
c             select case (iblk)
c             case (itp_tet)
c               write (fname2,"('material type interior linear tetrahedron')") 
c             case (itp_wedge)
c               write (fname2,"('material type interior linear wedge')") 
c             end select
      write (fname2,"('material type interior?')") 
           endif
c
C           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           intfromfile(:)=-1
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), 1, dataInt, iotype)
c
           allocate (mattype(intfromfile(1)))
c
           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(mattype), intfromfile(1), dataInt, iotype)

           if(writeLock==0) then
c
c.... mattype_interior is the global material type for each element
c     it is used for visualization
c
           mattype_interior(numel_ct+1:numel_ct+intfromfile(1)) = mattype(:)

c ... count elements with the same mattype
c
             call count_elem_mattype(mattype,neltp,mat_tag(1:nummat,1),nummat)
c
             material_loop: do imattype = 1, nummat
c
               iptr = 1
c
               blocks_loop: do
c
                 npro = 0
c
                 do 
                   if (mattype(iptr) == mat_tag(imattype,1)) then
                     npro = npro + 1
                     ientmp(npro,1:nshl) = ientp(iptr,1:nshl)
                     ieMaptmp(npro) = numel_ct + iptr
                   endif
                   iptr = iptr + 1
                   if (npro == ibksz .or. iptr > neltp) exit
                 enddo
c
                 if (npro == 0) cycle material_loop
c
                 nelblk = nelblk + 1
                 lcblk(1,nelblk)  = iel
                 lcblk(3,nelblk)  = lcsyst
                 lcblk(4,nelblk)  = ipordl
                 lcblk(5,nelblk)  = nenl
                 lcblk(6,nelblk)  = nfacel
                 lcblk(7,nelblk)  = imattype
                 lcblk(8,nelblk)  = ndofl
                 lcblk(9,nelblk)  = nsymdl 
                 lcblk(10,nelblk) = nshl ! # of shape functions per elt
c
c.... allocate memory for stack arrays
c
                 allocate (mien(nelblk)%p(npro,nshl))
                 allocate (mxmudmi(nelblk)%p(npro,maxsh))
                 if(usingpetsc.eq.0) then
                    allocate (mienG(nelblk)%p(1,1))
                 else
                    allocate (mienG(nelblk)%p(npro,nshl))
                 endif
                 ! note mienG will be passed to gensav but nothing filled if not 
                 ! using PETSc so this is safe
c
c.... attach map array to pointer
                 allocate (mieMap(nelblk)%p(npro))
                 mieMap(nelblk)%p(1:npro) = ieMaptmp(1:npro)
c
c.... save the element block
c
c                n1=n
c                n2=n+npro-1
c                mater=1   ! all one material for now
c                call gensav (ientp(n1:n2,1:nshl),
c     &                       mater,           mien(nelblk)%p,
c     &                       mienG(nelblk)%p,
c     &                       mmat(nelblk)%p)
                 call gensav (mien(nelblk)%p,mienG(nelblk)%p)
                 iel=iel+npro
                 if (iptr > neltp) exit blocks_loop
               enddo blocks_loop
             enddo material_loop
             numel_ct = numel_ct + intfromfile(1)
           endif ! writeLock
           deallocate(ientp,ientmp)
           deallocate(ieMaptmp)
           deallocate(mattype,neltp_mattype)
        enddo iblk_loop
c
        if (numel_ct .ne. numel)
     &    call error ('genblk global material type  ','numel_ct  ',numel_ct)
c
        lcblk(1,nelblk+1) = iel
        return
1000    format(a80,//,
     &  ' N o d a l   C o n n e c t i v i t y',//,
     &  '   Elem  ',/,
     &  '  Number  ',7x,27('Node',i2,:,2x))
1100    format(2x,i5,6x,27i8)
        end
