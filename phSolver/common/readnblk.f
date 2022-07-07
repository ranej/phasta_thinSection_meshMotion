c  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
c
c    module readarrays ("Red Arrays") -- contains the arrays that
c     are read in from binary files but not immediately blocked 
c     through pointers.
c
c    subroutine readnblk ("Reed and Block") -- allocates space for
c     and reads data to be contained in module readarrays.  Reads
c     all remaining data and blocks them with pointers.
c
      module m2gfields
        integer, allocatable :: m2gClsfcn(:,:)
        real*8, allocatable :: m2gParCoord(:,:)
      end module

      module rigidBodyReadData
        integer              :: rbUseReadData
        integer, allocatable :: rbIDs(:)
        integer, allocatable :: rbMTs(:)
        integer, allocatable :: rbFlags(:)
        real*8,  allocatable :: rbParamRead(:,:)
      end module

      module interfaceflag
        integer, allocatable :: ifFlag(:)
      end module

      module BLparameters
        real*8, allocatable  :: BLflt(:)
        real*8, allocatable  :: BLgr(:)
        integer, allocatable :: BLtnv(:)
        integer, allocatable :: BLlist(:)
        integer*8, allocatable :: BLflag(:)
      end module

      module TSparameters
        integer, allocatable :: TStnv(:)
        integer, allocatable :: TSlist(:)
        integer*8, allocatable :: TSflag(:)
      end module

      module readarrays

      use m2gfields
      use interfaceflag
      use rigidBodyReadData
      use BLparameters
      use TSparameters

      real*8, allocatable :: point2x(:,:)
      real*8, allocatable :: qold(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: acold(:,:)
      real*8, allocatable :: xn(:,:)
      real*8, allocatable :: xdotold(:,:)
      real*8, allocatable :: umesh(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable :: BCinp(:,:)

      integer, allocatable :: point2ilwork(:)
!      integer, allocatable :: fncorp(:)
      integer, allocatable :: twodncorp(:,:)
      integer, allocatable :: nBC(:)
      integer, allocatable :: point2iper(:)
      integer, target, allocatable :: point2ifath(:)
      integer, target, allocatable :: point2nsons(:)

      end module

      subroutine readnblk
      use iso_c_binding
      use readarrays
      use fncorpmod
      use phio
      use phstr
      use syncio
      use posixio
      use streamio
      use solid_m
      include "common.h"

      real*8, target, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, target, allocatable :: uread(:,:), xnread(:,:)
      real*8, target, allocatable :: xdotoldread(:,:), umeshread(:,:)
      real*8, target, allocatable :: BCinpread(:,:)
      real*8, target, allocatable :: tmpBLDbl(:)
      real*8, target, allocatable :: tmpTSDbl(:)
      real*8, target, allocatable :: tmpm2gParCoord(:,:)
      real*8 :: iotime
      integer, target, allocatable :: iperread(:), iBCtmpread(:)
      integer, target, allocatable :: ilworkread(:), nBCread(:)
      integer, target, allocatable :: fncorpread(:)
      integer, target, allocatable :: tmpBLInt(:), tmpBLlist(:)
      integer, target, allocatable :: tmpTSInt(:), tmpTSlist(:)
      integer, target, allocatable :: tmpm2gClsfcn(:,:)
      integer, target, allocatable :: tmpifFlag(:)
      integer, target, allocatable :: tmprbIDs(:), tmprbMTs(:), tmprbFlags(:)
      real*8, target, allocatable  :: tmprbParamRead(:,:)
      integer fncorpsize
      character*10 cname2, cname2nd
      character*8 mach2
      character*30 fmt1
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      integer igeomBAK, ierr
      integer, target :: intfromfile(50) ! integers read from headers
      integer :: descriptor, descriptorG, GPID, color, nfields
      integer ::  numparts, nppf
      integer :: ierr_io, numprocs, itmp, itmp2
      integer :: ignored
      integer :: fileFmt
      integer :: numm2g, ixsiz
      integer :: listcounter, ioffset, ngc, itnv, basevID
      integer :: numel_ct
      character*255 fname2, temp2
      character*64 temp1
      type(c_ptr) :: handle
      character(len=1024) :: dataInt, dataDbl
      integer, target, allocatable :: itemp(:)
      dataInt = c_char_'integer'//c_null_char
      dataDbl = c_char_'double'//c_null_char
c
c.... determine the step number to start with
c
      open(unit=72,file='numstart.dat',status='old')
      read(72,*) irstart
      close(72)
      lstep=irstart ! in case restart files have no fields

      fname1='geombc.dat'
      fname1= trim(fname1)  // cname2(myrank+1)

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
      write (fmt1,"('(''restart.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(myrank+1)
c
c.... open input files
c.... input the geometry parameters
c
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

      itwo=2
      ione=1
      ieleven=11
      itmp = int(log10(float(myrank+1)))+1

      iotime = TMRC()
      if( input_mode .eq. -1 ) then
        call streamio_setup_read(fhandle, geomRestartStream)
      else if( input_mode .eq. 0 ) then
        call posixio_setup(fhandle, c_char_'r')
      else if( input_mode .ge. 1 ) then
        call syncio_setup_read(nsynciofiles, fhandle)
      end if
      call phio_constructName(fhandle, 
     &        c_char_'geombc' // char(0), fname1)
      call phio_openfile(fname1, fhandle);

      call phio_readheader(fhandle,c_char_'number of nodes' // char(0),
     & c_loc(numnp),ione, dataInt, iotype)

      call phio_readheader(fhandle,c_char_'number of modes' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior elements' // char(0),
     &  c_loc(numel),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of boundary elements' // char(0),
     &  c_loc(numelb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interface elements' // char(0),
     &  c_loc(numelif),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'maximum number of element nodes' // char(0),
     &  c_loc(nen),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior tpblocks' // char(0),
     &  c_loc(nelblk),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of boundary tpblocks' // char(0),
     & c_loc(nelblb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of interface tpblocks' // char(0),
     & c_loc(nelblif),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of nodes with Dirichlet BCs' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of shape functions' // char(0),
     & c_loc(ntopsh),ione, dataInt, iotype)
c
c.... calculate the maximum number of boundary element nodes
c     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   ','nen     ',nen)
      endif
c
c.... setup some useful constants
c
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
      if(matflg(1,1).lt.0) then
         nflow = nsd + 1
      else
         nflow = nsd + 2
      endif 
c
      ndof   = nsd + 2
      nsclr=impl(1)/100
      ndof=ndof+nsclr           ! number of sclr transport equations to solve
c      
      if (iALE .eq. 2) then     ! Mesh-elastic is ON
         nelas   = nsd              ! FOR mesh-elastic 
         ndofBC  = ndof + I3nsd     ! dimension of BC array
     &          + nelas + I3nsd    ! add nelas for mesh-elastic solve
         ndofBC2 = 3+2+4+7+8    ! (assuming 4 scalars to be ON) and 8 is for ec11 ec12 ec13 em1 ec21 ec22 ec23 em2 
      else
         nelas   = 0
         ndofBC  = ndof + I3nsd     ! dimension of BC array
         ndofBC2 = ndof + 7 
      endif
c
      ndiBCB = 2                ! dimension of iBCB array
      ndBCB  = ndof + 1         ! dimension of BCB array
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
c
c.... ----------------------> Communication tasks <--------------------
c
      if(numpe > 1) then
         call phio_readheader(fhandle,
     &    c_char_'size of ilwork array' // char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         call phio_readheader(fhandle,
     &    c_char_'ilwork' //char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         allocate( point2ilwork(nlwork) )
         allocate( ilworkread(nlwork) )
         call phio_readdatablock(fhandle, c_char_'ilwork' // char(0),
     &      c_loc(ilworkread), nlwork, dataInt , iotype)

         point2ilwork = ilworkread
         call ctypes (point2ilwork)

       if((usingPETSc.eq.1).or.(svLSFlag.eq.1)) then
         fncorpsize = nshg
         allocate(fncorp(fncorpsize))
         call gen_ncorp(fncorp, ilworkread, nlwork, fncorpsize)
!
! the  following code finds the global range of the owned nodes
!
         maxowned=0
         minowned=maxval(fncorp)
         do i = 1,nshg      
            if(fncorp(i).gt.0) then  ! don't consider remote copies
               maxowned=max(maxowned,fncorp(i))
               minowned=min(minowned,fncorp(i))
            endif
         enddo
!
!  end of global range code
!
         call commuInt(fncorp, point2ilwork, 1, 'out')
         ncorpsize = fncorpsize 
       endif
       if(svLSFlag.eq.1) then
         allocate(ltg(ncorpsize))
         ltg(:)=fncorp(:)
       endif
      else
           allocate(ltg(nshg))
           do i =1,nshg
             ltg(i)=i
           enddo
           nlwork=1
           allocate( point2ilwork(1))
           nshg0 = nshg
      endif


      itwo=2

      call phio_readheader(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(intfromfile),itwo, dataDbl, iotype)
      numnp=intfromfile(1)
      allocate( point2x(numnp,nsd) )
      allocate( xread(numnp,nsd) )
      ixsiz=numnp*nsd
      call phio_readdatablock(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(xread),ixsiz, dataDbl, iotype)
      point2x = xread
c
c.... read mesh to geom fields
c
        intfromfile=0
        call phio_readheader(fhandle,
     &   c_char_'m2g classification' // char(0),
     &   c_loc(intfromfile),itwo, dataInt, iotype)
        numm2g = intfromfile(1)
        ixsiz=numnp*3 ! dim and tag and discrete
        if (numm2g > 0) then
          if (numm2g .ne. numnp) then
            write(*,*) "size of m2g field is not consistent with numnp"
            call error ('readnblk  ', 'numnp   ', numnp)
          endif
          allocate( tmpm2gClsfcn(numnp,3) )
          allocate( m2gClsfcn(numnp,3) )
          mesh2geom = 1
          call phio_readdatablock(fhandle,
     &     c_char_'m2g classification' // char(0),
     &     c_loc(tmpm2gClsfcn),ixsiz, dataInt, iotype)
          m2gClsfcn = tmpm2gClsfcn
          deallocate( tmpm2gClsfcn )
        else
          allocate( m2gClsfcn(1,1) )
          mesh2geom = 0
          m2gClsfcn = 0
        endif
c...
        intfromfile=0;
        call phio_readheader(fhandle,
     &   c_char_'m2g parametric coordinate' // char(0),
     &   c_loc(intfromfile),itwo, dataDbl, iotype)
        numm2g = intfromfile(1)
        ixsiz=numnp*2 ! par1 and par2
        if (numm2g > 0) then
          if (numm2g .ne. numnp) then
            write(*,*) "size of m2g field is not consistent with numnp"
            call error ('readnblk  ', 'numnp   ', numnp)
          endif
          allocate( tmpm2gParCoord(numnp,2) )
          allocate( m2gParCoord(numnp,2) )
          mesh2geom = 1 ! may not needed
          call phio_readdatablock(fhandle,
     &     c_char_'m2g parametric coordinate' // char(0),
     &     c_loc(tmpm2gParCoord),ixsiz, dataDbl, iotype)
          m2gParCoord = tmpm2gParCoord
          deallocate( tmpm2gParCoord )
        else
          allocate( m2gParCoord(1,1) )
          mesh2geom = 0 ! may not needed
          m2gParCoord = 0
        endif
c
c.... end read mesh to geom fields
c
c
c.... read in and block out the connectivity
c
      call genblk (IBKSIZ)
c
c.... read the boundary condition mapping array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      allocate( nBC(nshg) )

      allocate( nBCread(nshg) )

      call phio_readdatablock(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nBCread), nshg, dataInt, iotype)

      nBC=nBCread
c
c.... read the temporary iBC array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      if ( numpbc > 0 ) then
        allocate( iBCtmp(numpbc) )
        allocate( iBCtmpread(numpbc) )
      else
        allocate( iBCtmp(1) )
        allocate( iBCtmpread(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(iBCtmpread), numpbc, dataInt, iotype)

      if ( numpbc > 0 ) then
         iBCtmp=iBCtmpread
      else  ! sometimes a partition has no BC's
         deallocate( iBCtmpread)
         iBCtmp=0
      endif
c
c.... read boundary condition data
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(intfromfile),ione, dataDbl, iotype)

      if ( numpbc > 0 ) then
c         if(intfromfile(1).ne.(ndof+7)*numpbc) then
c           warning='WARNING more data in BCinp than needed: keeping 1st'
c           write(*,*) warning, ndof+7, numpbc,intfromfile(1),(ndof+7)*numpbc
c         endif
         allocate( BCinp(numpbc,ndofBC2) )
         nsecondrank=intfromfile(1)/numpbc
         allocate( BCinpread(numpbc,nsecondrank) )
         iBCinpsiz=intfromfile(1)
      else
         allocate( BCinp(1,ndofBC2) )
         allocate( BCinpread(0,0) ) !dummy
         iBCinpsiz=intfromfile(1)
      endif

      call phio_readdatablock(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(BCinpread), iBCinpsiz, dataDbl, iotype)

      if ( numpbc > 0 ) then
         BCinp(:,1:ndofBC2)=BCinpread(:,1:ndofBC2)
      else  ! sometimes a partition has no BC's
         deallocate(BCinpread)
         BCinp=0
      endif
c
c.... read global interface flag
c
      ione=1
      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'DG interface flag' // char(0),
     & c_loc(intfromfile), ione, dataInt, iotype)

      allocate( ifFlag(nshg) )
      if ( intfromfile(1) .gt. 0 ) then
        if ( intfromfile(1) .ne. nshg ) then
          call error ('readnblk  ', 'size of interface flag ', intfromfile(1))
        endif
        allocate( tmpifFlag(nshg) )
      else
        allocate( tmpifFlag(1) )
      endif

      call phio_readdatablock(fhandle,
     & c_char_'DG interface flag' // char(0),
     & c_loc(tmpifFlag), intfromfile(1), dataInt, iotype)

      if ( intfromfile(1) .gt. 0 ) then
         ifFlag = tmpifFlag
         deallocate( tmpifFlag )
      else  ! sometimes a partition has no interface
         ifFlag = zero
         deallocate( tmpifFlag )
      endif
c
c--------------------- read rigid body tag ------------------
c
c.... read IDs and model tags
        ione=1
        intfromfile=0
        if (numrbs > 0) then
          allocate( tmprbIDs(numrbs) )
          allocate( tmprbMTs(numrbs) )
          allocate( rbIDs(numrbs) )
          allocate( rbMTs(numrbs) )
          rbIDs = -1
          rbMTs = -1
c
          call phio_readheader(fhandle,
     &     c_char_'rigid body IDs' // char(0),
     &     c_loc(intfromfile),ione, dataInt, iotype)
          if(intfromfile(1) .ne. numrbs) then
            call error ('readnblk  ', 'num of rbs not equal input'
     &                  , intfromfile(1))
          endif
          call phio_readdatablock(fhandle,
     &     c_char_'rigid body IDs' // char(0),
     &     c_loc(tmprbIDs),numrbs, dataInt, iotype)
c
          call phio_readheader(fhandle,
     &     c_char_'rigid body MTs' // char(0),
     &     c_loc(intfromfile),ione, dataInt, iotype)
          if(intfromfile(1) .ne. numrbs) then
            call error ('readnblk  ', 'num of rbs not equal input'
     &                  , intfromfile(1))
          endif
          call phio_readdatablock(fhandle,
     &     c_char_'rigid body MTs' // char(0),
     &     c_loc(tmprbMTs),numrbs, dataInt, iotype)
c
          do i = 1, numrbs
            do j = 1, numrbs
              if(tmprbIDs(i) .eq. rbsTags(j)) then
                rbIDs(j) = tmprbIDs(i)
                rbMTs(j) = tmprbMTs(i)
              endif
            enddo
          enddo
          deallocate( tmprbIDs )
          deallocate( tmprbMTs )
        else
          allocate( rbIDs(1) )
          allocate( rbMTs(1) )
          rbIDs = -1
          rbMTs = -1
          numrbs = 0 ! make sure it is not negative
        endif
c
c.... read tag for each vertex
      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'rigid body tag' // char(0),
     & c_loc(intfromfile), ione, dataInt, iotype)

      if ( intfromfile(1) .gt. 0 ) then
        if ( intfromfile(1) .ne. numnp ) then
          call error ('readnblk  ', 'size of rigid body tag ', intfromfile(1))
        endif
        allocate( tmprbFlags(numnp) )
        allocate( rbFlags(numnp) )
        rbFlags = 0
        call phio_readdatablock(fhandle,
     &   c_char_'rigid body tag' // char(0),
     &   c_loc(tmprbFlags), intfromfile(1), dataInt, iotype)
        do i = 1, numnp
          do j = 1, numrbs
            if(tmprbFlags(i) .eq. rbsTags(j)) rbFlags(i) = j
          enddo
        enddo
        deallocate( tmprbFlags )
      else
        allocate( rbFlags(1) )
        rbFlags = 0
      endif
c
c--------------------- end read rigid body tag --------------
c
c--------------------- read the layered mesh parameters ------------------
c
c.... first layer thickness
      ione=1
      numgc = 0
      call phio_readheader(fhandle,
     & c_char_'first layer thickness' // char(0),
     & c_loc(numgc),ione, dataDbl, iotype)
      if ( numgc > 0 ) then
        allocate( tmpBLDbl(numgc) )
        allocate( BLflt(numgc) )
      else
        allocate( tmpBLDbl(1) )
        allocate( BLflt(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'first layer thickness' // char(0),
     & c_loc(tmpBLDbl), numgc, dataDbl, iotype)

      if ( numgc > 0 ) then
         BLflt = tmpBLDbl
         deallocate( tmpBLDbl )
      else  ! sometimes a partition has no BL
         BLflt = 0
         deallocate( tmpBLDbl )
      endif
c
c.... growth ratio
      ione=1
      call phio_readheader(fhandle,
     & c_char_'growth ratio' // char(0),
     & c_loc(numgc),ione, dataDbl, iotype)
      if ( numgc > 0 ) then
        allocate( tmpBLDbl(numgc) )
        allocate( BLgr(numgc) )
      else
        allocate( tmpBLDbl(1) )
        allocate( BLgr(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'growth ratio' // char(0),
     & c_loc(tmpBLDbl), numgc, dataDbl, iotype)

      if ( numgc > 0 ) then
         BLgr = tmpBLDbl
         deallocate( tmpBLDbl )
      else  ! sometimes a partition has no BL
         BLgr = 0
         deallocate( tmpBLDbl )
      endif
c
c.... total number of layers
      ione=1
      call phio_readheader(fhandle,
     & c_char_'number of vertices on growth curve' // char(0),
     & c_loc(numgc),ione, dataInt, iotype)
      if ( numgc > 0 ) then
        allocate( tmpBLInt(numgc) )
        allocate( BLtnv(numgc) )
      else
        allocate( tmpBLInt(1) )
        allocate( BLtnv(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'number of vertices on growth curve' // char(0),
     & c_loc(tmpBLInt), numgc, dataInt, iotype)

      if ( numgc > 0 ) then
         BLtnv = tmpBLInt
         deallocate( tmpBLInt )
      else  ! sometimes a partition has no BL
         BLtnv = 0
         deallocate( tmpBLInt )
      endif
c
c.... growth curve connectivity
      ione=1
      numgcnp = 0
      call phio_readheader(fhandle,
     & c_char_'list of vertices on growth curve' // char(0),
     & c_loc(numgcnp),ione, dataInt, iotype)
      if ( numgcnp > 0 ) then
        allocate( tmpBLlist(numgcnp) )
        allocate( BLlist(numgcnp) )
      else
        allocate( tmpBLlist(1) )
        allocate( BLlist(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'list of vertices on growth curve' // char(0),
     & c_loc(tmpBLlist), numgcnp, dataInt, iotype)

      if ( numgcnp > 0 ) then
         BLlist = tmpBLlist
         deallocate( tmpBLlist )
      else  ! sometimes a partition has no BL
         BLlist = 0
         deallocate( tmpBLlist )
      endif
c
c.... growth curve base vertex flag
      allocate( BLflag(nshg) )
      BLflag = zero
      if (numgc > 0) then
        listcounter = 0
        ioffset = 1 ! the ID starts from 1 in phasta
        do ngc = 1, numgc
          itnv = BLtnv(ngc) ! number of vertices on this growth curve
          basevID = BLlist(listcounter + 1) + ioffset
          BLflag(basevID) = 1
          listcounter = listcounter + itnv ! update counter
        enddo
      endif
c
      if ( numpe > 1 ) then
        call commuInt(BLflag, point2ilwork, 1, 'out')
      endif
c
c------------------ end read the layered mesh parameters ------------------
c
cDebug
c--------------------- read the Thin section layered mesh parameters ------------------
c
c.... total number of layers
      ione=1
      numts=0
      call phio_readheader(fhandle,
     & c_char_'number of vertices on thin section stack' // char(0),
     & c_loc(numts),ione, dataInt, iotype)
      if ( numts > 0 ) then
        allocate( tmpTSInt(numts) )
        allocate( TStnv(numts) )
      else
        allocate( tmpTSInt(1) )
        allocate( TStnv(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'number of vertices on thin section stack' // char(0),
     & c_loc(tmpTSInt), numts, dataInt, iotype)

      if ( numts > 0 ) then
         TStnv = tmpTSInt
         deallocate( tmpTSInt )
      else  ! sometimes a partition has no TS
         TStnv = 0
         deallocate( tmpTSInt )
      endif
c
c.... thin section stack connectivity
      ione=1
      numtsnp = 0
      call phio_readheader(fhandle,
     & c_char_'list of vertices on thin section stack' // char(0),
     & c_loc(numtsnp),ione, dataInt, iotype)
      if ( numtsnp > 0 ) then
        allocate( tmpTSlist(numtsnp) )
        allocate( TSlist(numtsnp) )
      else
        allocate( tmpTSlist(1) )
        allocate( TSlist(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'list of vertices on thin section stack' // char(0),
     & c_loc(tmpTSlist), numtsnp, dataInt, iotype)

      if ( numtsnp > 0 ) then
         TSlist = tmpTSlist
         deallocate( tmpTSlist )
      else  ! sometimes a partition has no TS
         TSlist = 0
         deallocate( tmpTSlist )
      endif
c
c.... growth curve base vertex flag
      allocate( TSflag(nshg) )
      TSflag = zero
      if (numts > 0) then
        listcounter = 0
        ioffset = 1 ! the ID starts from 1 in phasta
        do nTS = 1, numts
          itnv = TStnv(nTS) ! number of vertices on this growth curve
          basevID = TSlist(listcounter + 1) + ioffset
          TSflag(basevID) = 1
          listcounter = listcounter + itnv ! update counter
        enddo
      endif
c
      if ( numpe > 1 ) then
        call commuInt(TSflag, point2ilwork, 1, 'out')
      endif
c
c------------------ end read the thin section layered mesh parameters ------------------
c
cDebug
c.... read periodic boundary conditions
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(nshg), ione, dataInt, iotype)
      allocate( point2iper(nshg) )
      allocate( iperread(nshg) )
      call phio_readdatablock(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(iperread), nshg, dataInt, iotype)
      point2iper=iperread
c
c.... generate the boundary element blocks
c
      call genbkb (ibksiz)
c
      call genbkif (ibksiz)
c
c  Read in the nsons and ifath arrays if needed
c
c  There is a fundamental shift in the meaning of ifath based on whether
c  there exist homogenous directions in the flow.  
c
c  HOMOGENOUS DIRECTIONS EXIST:  Here nfath is the number of inhomogenous
c  points in the TOTAL mesh.  That is to say that each partition keeps a 
c  link to  ALL inhomogenous points.  This link is furthermore not to the
c  sms numbering but to the original structured grid numbering.  These 
c  inhomogenous points are thought of as fathers, with their sons being all
c  the points in the homogenous directions that have this father's 
c  inhomogeneity.  The array ifath takes as an arguement the sms numbering
c  and returns as a result the father.
c
c  In this case nsons is the number of sons that each father has and ifath
c  is an array which tells the 
c
c  NO HOMOGENOUS DIRECTIONS.  In this case the mesh would grow to rapidly
c  if we followed the above strategy since every partition would index its
c  points to the ENTIRE mesh.  Furthermore, there would never be a need
c  to average to a node off processor since there is no spatial averaging.
c  Therefore, to properly account for this case we must recognize it and
c  inerrupt certain actions (i.e. assembly of the average across partitions).
c  This case is easily identified by noting that maxval(nsons) =1 (i.e. no
c  father has any sons).  Reiterating to be clear, in this case ifath does
c  not point to a global numbering but instead just points to itself.
c
      nfath=1  ! some architectures choke on a zero or undeclared
                 ! dimension variable.  This sets it to a safe, small value.
      if(((iLES .lt. 20) .and. (iLES.gt.0))
     &                   .or. (itwmod.gt.0)  ) then ! don't forget same
                                                    ! conditional in proces.f
                                                    ! needed for alloc
         ione=1
         if(nohomog.gt.0) then
            call phio_readheader(fhandle,
     &       c_char_'number of father-nodes' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            call phio_readheader(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            allocate (point2nsons(nfath))

            call phio_readdatablock(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(point2nsons),nfath, dataInt, iotype)

            call phio_readheader(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(nshg), ione, dataInt, iotype);

            allocate (point2ifath(nshg))

            call phio_readdatablock(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(point2ifath), nshg, dataInt, iotype)
     
            nsonmax=maxval(point2nsons)
         else  ! this is the case where there is no homogeneity
               ! therefore ever node is a father (too itself).  sonfath
               ! (a routine in NSpre) will set this up but this gives
               ! you an option to avoid that.
            nfath=nshg
            allocate (point2nsons(nfath))
            point2nsons=1
            allocate (point2ifath(nshg))
            do i=1,nshg
               point2ifath(i)=i
            enddo
            nsonmax=1
         endif
      else
         allocate (point2nsons(1))
         allocate (point2ifath(1))
      endif

      call phio_closefile(fhandle);
      iotime = TMRC() - iotime
      if (myrank.eq.master) then
        write(*,*) 'time to read geombc (seconds)', iotime
      endif
c
c-------------------------------------------------------------------------
c------------------------- Read restart files ----------------------------
c-------------------------------------------------------------------------
c
      iotime = TMRC()
      if( input_mode .eq. -1 ) then
        call streamio_setup_read(fhandle, geomRestartStream)
      else if( input_mode .eq. 0 ) then
        call posixio_setup(fhandle, c_char_'r')
      else if( input_mode .ge. 1 ) then
        call syncio_setup_read(nsynciofiles, fhandle)
      end if
      call phio_constructName(fhandle,
     &        c_char_'restart' // char(0), fnamer)
      call phstr_appendInt(fnamer, irstart)
      call phstr_appendStr(fnamer, c_char_'.'//c_null_char)
      call phio_openfile(fnamer, fhandle);

      ithree=3

      itmp = int(log10(float(myrank+1)))+1

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'solution' // char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
c
c.... read the values of primitive variables into q
c
      allocate( qold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         if(ndof2.ne.ndof) then

         endif
        if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( qread(nshg,ndof2) )
         iqsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'solution' // char(0),
     &    c_loc(qread),iqsiz, dataDbl,iotype)
         qold(:,1:ndof)=qread(:,1:ndof)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='Solution is set to zero (with p and T to one)'
            else
               warning='Solution is set to zero'
            endif
            write(*,*) warning
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'time derivative of solution' // char(0),
     & c_loc(intfromfile), ithree, dataInt, iotype)
      allocate( acold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)

         if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( acread(nshg,ndof2) )
         acread=zero
         iacsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'time derivative of solution' // char(0),
     &    c_loc(acread), iacsiz, dataDbl,iotype)
         acold(:,1:ndof)=acread(:,1:ndof)
         deallocate(acread)
      else
         if (myrank.eq.master) then
            warning='Time derivative of solution is set to zero (SAFE)'
            write(*,*) warning
         endif
         acold=zero
      endif
c
c... read solid part
c
      solid_p%is_active = any(mat_eos(:,1) .eq. ieos_solid_1)
c
      if (solid_p%is_active) 
     &  call read_restart_solid
c
c read in ALE stuff
c read in coordinate at n time step      
      intfromfile=0
      call phio_readheader(fhandle, 
     & c_char_'motion_coords' //char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
      allocate( xn(numnp,nsd) )
      if(intfromfile(1).ne.0) then 
         numnp2=intfromfile(1)
         nsd2=intfromfile(2)
         lstep=intfromfile(3)
         
         if (numnp2 .ne. numnp) 
     &        call error ('restar  ', 'numnp   ', numnp)
c     
         allocate( xnread(numnp,nsd2) )
         xnread=zero

         iacsiz=numnp*nsd2
         call phio_readdatablock(fhandle,
     &    c_char_'motion_coords' // char(0),
     &    c_loc(xnread), iacsiz, dataDbl,iotype)
         xn(:,1:nsd)=xnread(:,1:nsd)
         deallocate(xnread)
      else
         if (myrank.eq.master) then
           warning='Mesh coordinates set to original coordinates (SAFE)'
            write(*,*) warning
         endif
         xn=point2x
      endif
c read in xdotold
c      fname1='xdot?'
c      intfromfile=0
c      call phio_readheader(fhandle, 
c     & c_char_'xdot' //char(0), 
c     & c_loc(intfromfile), ithree, dataInt, iotype)
c      call readheader(irstin,fname1,intfromfile,
c     &     ithree,'integer', iotype)
c      allocate( xdotold(numnp,nsd) )
c      if(intfromfile(1).ne.0) then 
c         numnp2=intfromfile(1)
c         nsd2=intfromfile(2)
c         lstep=intfromfile(3)
c         
c         if (numnp2 .ne. numnp) 
c     &        call error ('restar  ', 'numnp   ', numnp)
c     
c         allocate( xdotoldread(numnp,nsd2) )
c         xdotoldread=zero
c
c         iacsiz=numnp*nsd2
c
c         call phio_readdatablock(fhandle,
c     &    c_char_'xdot' // char(0),
c     &    c_loc(xdotoldread), iacsiz, dataDbl,iotype)
c         call readdatablock(irstin,fname1,xdotoldread,iacsiz,
c     &                   'double',iotype)
c         xdotold(:,1:nsd)=xdotoldread(:,1:nsd)
c         deallocate(xdotoldread)
c      else
c         if (myrank.eq.master) then
c            warning='Time derivative of mesh disp is set to zero (SAFE)'
c            write(*,*) warning
c         endif
c         xdotold=zero
c      endif
c read in umesh
      intfromfile=0
       call phio_readheader(fhandle, 
     & c_char_'mesh_vel' //char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
       allocate( umesh(numnp,nsd) )
       if(intfromfile(1).ne.0) then 
          numnp2=intfromfile(1)
          nsd2=intfromfile(2)
          lstep=intfromfile(3)
          
          if (numnp2 .ne. numnp) 
     &        call error ('restar  ', 'numnp   ', numnp)
      
          allocate( umeshread(numnp,nsd2) )
          umeshread=zero
 
          iacsiz=numnp*nsd2
          call phio_readdatablock(fhandle,
     &    c_char_'mesh_vel' // char(0),
     &    c_loc(umeshread), iacsiz, dataDbl,iotype)
          umesh(:,1:nsd)=umeshread(:,1:nsd)
          deallocate(umeshread)
       else
          if (myrank.eq.master) then
             warning='Mesh velocity is set to zero (SAFE)'
             write(*,*) warning
          endif
         umesh=zero
       endif
c
c.... end read ALE stuff
c
c
c.... read in rigid body data
c
       rbUseReadData = 0
       if (numrbs .gt. 0) then
         intfromfile=0
         call phio_readheader(fhandle,
     &   c_char_'rbParams' //char(0),
     &   c_loc(intfromfile), itwo, dataInt, iotype)
c
         allocate( rbParamRead(numrbs, rbParamSize) )
         rbParamRead = zero
c
         if(intfromfile(1).ne.0) then
           if (intfromfile(1) .ne. numrbs)
     &       call error ('restar  ', 'rbParams size 1', numrbs)
           if (intfromfile(2) .ne. rbParamSize)
     &       call error ('restar  ', 'rbParams size 2', rbParamSize)
c
           rbUseReadData = 1 ! turn on use read data flag
           allocate( tmprbParamRead(numrbs, rbParamSize) )
           tmprbParamRead = zero
c
           iacsiz=numrbs * rbParamSize
           call phio_readdatablock(fhandle,
     &     c_char_'rbParams' // char(0),
     &     c_loc(tmprbParamRead), iacsiz, dataDbl, iotype)
c
           rbParamRead(:,:) = tmprbParamRead(:,:)
c
           deallocate( tmprbParamRead )
c
         else
           if (myrank.eq.master) then
             warning='Rigid body data is set to zero (SAFE)'
             write(*,*) warning
           endif
         endif
c
       endif ! end if numrbs greater than 0
c
c.... end read rigid body data
c
cc
cc.... read the header and check it against the run data
cc
      if (ideformwall.eq.1) then

          intfromfile=0
          call phio_readheader(fhandle,
     &     c_char_'displacement' // char(0),
     &     c_loc(intfromfile), ithree, dataInt, iotype)

         nshg2=intfromfile(1)
         ndisp=intfromfile(2)
         lstep=intfromfile(3)
         if(ndisp.ne.nsd) then
            warning='WARNING ndisp not equal nsd'
            write(*,*) warning , ndisp
         endif
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
c
c.... read the values of primitive variables into uold
c
         allocate( uold(nshg,nsd) )
         allocate( uread(nshg,nsd) )
         
         iusiz=nshg*nsd

         call phio_readdatablock(fhandle,
     &    c_char_'displacement' // char(0),
     &    c_loc(uread), iusiz, dataDbl, iotype)

         uold(:,1:nsd)=uread(:,1:nsd)
       else
         allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
       endif
c
c.... close c-binary files
c
      call phio_closefile(fhandle)
      iotime = TMRC() - iotime
      if (myrank.eq.master) then
        write(*,*) 'time to read restart (seconds)', iotime
      endif

      deallocate(xread)
      if ( numpbc > 0 )  then
         deallocate(bcinpread)
         deallocate(ibctmpread)
      endif
      deallocate(iperread)
      if(numpe.gt.1)
     &     deallocate(ilworkread)
      deallocate(nbcread)

      return
 994  call error ('input   ','opening ', igeomBAK)
 995  call error ('input   ','opening ', igeomBAK)
 997  call error ('input   ','end file', igeomBAK)
 998  call error ('input   ','end file', igeomBAK)
      end
c
c No longer called but kept around in case....
c
      subroutine genpzero(iBC)

      use pointer_data
      include "common.h"
      integer iBC(nshg)
c
c....  check to see if any of the nodes have a dirichlet pressure
c
      pzero=1
      if (any(btest(iBC,2))) pzero=0  
      do iblk = 1, nelblb
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         do i=1, npro
            iBCB1=miBCB(iblk)%p(i,1)
c     
c.... check to see if any of the nodes have a Neumann pressure 
c     but not periodic (note that 
c     
            if(btest(iBCB1,1)) pzero=0
         enddo
c     
c.... share results with other processors
c     
         pzl=pzero
         if (numpe .gt. 1)
     &        call MPI_ALLREDUCE (pzl, pzero, 1,
     &        MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
      enddo
      return
      end
