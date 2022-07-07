        subroutine ElmGMRElas (x,       disp,    shp,     shgl,
     &                         iBC,     BC,      shpb,    shglb,
     &                         shpif,   elasres, elasBDiag,
     &                         iper,    ilwork,  elaslhsK,
     &                         col,     row,     meshq)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual
c vector, and the preconditioning matrix, for mesh-elastic solve
c
c----------------------------------------------------------------------
c
         use pointer_data
         use timedataC
         use readarrays ! read BLflt, BLgr, BLtnv, BLlist
         use m2gfields ! read m2g fields
         use interfaceflag ! read ifFlag
c
        include "common.h"
        include "mpif.h"
c
        integer col(nshg+1), row(nnz*nshg)
        real*8  elaslhsK(nelas*nelas,nnz_tot),
     &          meshq(numel),   meshV(numel)
        real*8, allocatable :: meshqbk(:),  meshVbk(:)
c
        real*8  gcnormal(nshg, nsd)
        real*8  tsnormal(nshg, nsd)
c
        dimension x(numnp,nsd),        disp(numnp,nsd),
     &            xtmp(numnp,nsd),     iBC(nshg),
     &            BC(nshg,ndofBC),
     &            elasres(nshg,nelas),
     &            elasBDiag(nshg,nelas,nelas),
     &            iper(nshg)
c
        real*8  flag(numnp),   dispCommu(numnp,4)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)
c
        real*8, dimension(maxtop,    maxsh,maxqpt) :: shpif
c
        dimension ilwork(nlwork)
c
        integer errorcount(2)
c
        integer listcounter, ngc, itnv, nv, vID, vID2,
     &          ioffset, basevID, secondvID, nts, endvID
        real*8  iflt, igr, igrexp, tmp
        real*8  dist 
        real*8  inormal(nsd)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: Estiff(:,:,:)
c
c.... ------------------->   layer base elements   <-------------------
c
      gcnormal = zero
      if (numgc .gt. 0) then
c
c.... calculate the normal of each growth curve based on new boundary positions
c
        xtmp = x + disp
c
c.... loop over the boundary elements
c
        boundary_blocks: do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,             shpb(lcsyst,1:nshl,:),
     &                         mienb(iblk)%p,    miBCB(iblk)%p,
     &                         gcnormal)
c
        enddo boundary_blocks ! end loop the boundary elements
c
c.... loop over the interface elements
c
        interface_blocks: do iblk = 1, nelblif
c
          iel     = lcblkif(1, iblk)
          npro    = lcblkif(1,iblk+1) - iel
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          ipord   = lcblkif(5, iblk)    ! polynomial order
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          nenbl   = lcblkif(8, iblk)    ! number of vertices on boundary
          mater0  = lcblkif(9, iblk)
          mater1  = lcblkif(10,iblk)
          nshl0   = lcblkif(iblkif_nshl0,iblk)
          nshl1   = lcblkif(iblkif_nshl1,iblk)
          itpid   = lcblkif(iblkif_topology,iblk)
          ngaussif = nintif(itpid)
c
c.... hack; DG interface doesn't support 2nd and higher order
c
          if(ipord .eq. 1) then
            nshlb   = nenbl;
          else if(ipord .gt. 1) then
            write(*,*) "need to implement for higher order"
            call error('elmgmrelas  ','higher order', ipord)
          endif
c
c.... the 0 side
c
          lcsyst = lcsyst0
          nenl = nenl0
          nshl = nshl0
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif0(iblk)%p, miBCB(iblk)%p,
     &                         gcnormal)
c
c.... end of the 0 side
c
c
c.... the 1 side
c
          lcsyst = lcsyst1
          nenl = nenl1
          nshl = nshl1
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif1(iblk)%p, miBCB(iblk)%p,
     &                         gcnormal)
c
c.... end of the 1 side
c
        enddo interface_blocks ! end loop the interface elements
c
c.... end calculation of growth curve normal
c
      endif ! end if numgc greater than 0
c
c.... communication
c
      if (numpe > 1) then
        call commu (gcnormal  , ilwork, nsd  , 'in ')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call commu (gcnormal  , ilwork, nsd  , 'out')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
c
c.... ---------------->   Re-position layered mesh   <-----------------
c
c.... if no growth curve, skip whole re-positioning method
c
      flag = zero
      if (numgc .gt. 0) then
c
c.... loop over growth curves
c
        listcounter = 0
        ioffset = 1 ! the ID starts from 1 in phasta
        do ngc = 1, numgc
          itnv = BLtnv(ngc) ! number of vertices on this growth curve
          basevID = BLlist(listcounter + 1) + ioffset
c
c.... precaution
c
          if (itnv .lt. 2) then
            listcounter = listcounter + itnv
            cycle ! not loop over vertices
          endif
c.... no BC on base, skip this growth curve
          if (ibits(iBC(basevID),14,3) .eq. 0) then
            cycle ! not loop over vertices
          endif
c
c.... prepare other paramteres
c
          iflt = BLflt(ngc) ! first layer thickness of this growth curve
          igr  = BLgr(ngc)  ! growth ratio of this growth curve
c
c.... if the first and second vertex on boundary, use current normal
c     this assumes that the growth curve is exposed to a flat boundary
c
          if (mesh2geom .eq. 1) then
            secondvID = BLlist(listcounter + 2) + ioffset
c
            if ( (m2gClsfcn(basevID,1).ne.3) .and.
     &           (m2gClsfcn(secondvID,1).ne.3) ) then
              gcnormal(basevID,:) = x(secondvID,:) - x(basevID,:)
            endif
          endif
c
c... if base face option is set to be interface-only,
c    use current normal if base is not on inteface
          if (gcBaseOpt .eq. 1) then
            if (ifFlag(basevID) .ne. 1) then
              secondvID = BLlist(listcounter + 2) + ioffset
              gcnormal(basevID,:) = x(secondvID,:) - x(basevID,:)
            endif
          endif
c
          tmp  = sqrt( gcnormal(basevID,1) * gcnormal(basevID,1)
     &               + gcnormal(basevID,2) * gcnormal(basevID,2)
     &               + gcnormal(basevID,3) * gcnormal(basevID,3) )
          inormal(:) = gcnormal(basevID,:) / tmp
c
c.... loop over vertices on this growth curve
c
          igrexp = 1.0
          do nv = 2, itnv
            vID = BLlist(listcounter + nv) + ioffset
            vID2= BLlist(listcounter + nv - 1) + ioffset ! the previous one
            xtmp(vID,:) = xtmp(vID2,:) + iflt * inormal(:) * igrexp  ! igr**(nv-2)
            disp(vID,:) = xtmp(vID,:) - x(vID,:)
            igrexp = igrexp * igr

c... set flag
            call iBCelas_to_dbl(iBC(basevID), flag(vID))
          enddo ! over this growth curve
c
          listcounter = listcounter + itnv ! update counter
c
c.... end loop growth curves
c
        enddo
c
c.... end if growth curve
c
      endif ! end if numgc greater than 0
c
c.... -------------------->     communication     <--------------------
c
c.... combine disp + flag into dispCommu
c
      do i = 1, numnp
        dispCommu(i,1:3) = disp(i,:)
        dispCommu(i,4)   = flag(i)
      enddo
c
c.... turn on commu flag
c
      layerCommuFlag = 1
c
c.... commu
c
      if (numpe > 1) then
        call commu (dispCommu  , ilwork, 4  , 'in ')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call commu (dispCommu  , ilwork, 4  , 'out')
      endif
c
c.... turn off commu flag
c
      layerCommuFlag = 0
c
c.... separate dispCommu into disp and iBC; update BC
c
      do i = 1, numnp
        disp(i,:) = dispCommu(i,1:3)
        flag(i)   = dispCommu(i,4)
        call dbl_to_iBCelas(flag(i), iBC(i))
        call setBLbc(disp(i,:), iBC(i), BC(i, ndof+2:ndof+5))
      enddo
c
c.... end re-position layered mesh
c


c... Debug
c.... ------------------->   Thin section base elements   <-------------------
c
      tsnormal = zero
      if (numts .gt. 0) then
c
c.... calculate the normal of each thin section stack based on new boundary positions
c
        xtmp = x + disp
c
c.... loop over the boundary elements
c
        boundary_blocks_thin: do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble non-unit normal
c
          call calc_ts_normal (xtmp,             shpb(lcsyst,1:nshl,:),
     &                         mienb(iblk)%p,    miBCB(iblk)%p,
     &                         tsnormal)
c
        enddo boundary_blocks_thin ! end loop the boundary elements
c
c.... loop over the interface elements
c
        interface_blocks_thin: do iblk = 1, nelblif
c
          iel     = lcblkif(1, iblk)
          npro    = lcblkif(1,iblk+1) - iel
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          ipord   = lcblkif(5, iblk)    ! polynomial order
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          nenbl   = lcblkif(8, iblk)    ! number of vertices on boundary
          mater0  = lcblkif(9, iblk)
          mater1  = lcblkif(10,iblk)
          nshl0   = lcblkif(iblkif_nshl0,iblk)
          nshl1   = lcblkif(iblkif_nshl1,iblk)
          itpid   = lcblkif(iblkif_topology,iblk)
          ngaussif = nintif(itpid)
c
c.... hack; DG interface doesn't support 2nd and higher order
c
          if(ipord .eq. 1) then
            nshlb   = nenbl;
          else if(ipord .gt. 1) then
            write(*,*) "need to implement for higher order"
            call error('elmgmrelas  ','higher order', ipord)
          endif
c
c.... the 0 side
c
          lcsyst = lcsyst0
          nenl = nenl0
          nshl = nshl0
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble non-unit normal
c
          call calc_ts_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif0(iblk)%p, miBCB(iblk)%p,
     &                         tsnormal)
c
c.... end of the 0 side
c
c
c.... the 1 side
c
          lcsyst = lcsyst1
          nenl = nenl1
          nshl = nshl1
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble non-unit normal
c
          call calc_ts_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif1(iblk)%p, miBCB(iblk)%p,
     &                         tsnormal)
c
c.... end of the 1 side
c
        enddo interface_blocks_thin ! end loop the interface elements
c
c.... end calculation of growth curve normal
c
      endif ! end if numts greater than 0
c
c.... communication
c
      if (numpe > 1) then
        call commu (tsnormal  , ilwork, nsd  , 'in ')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call commu (tsnormal  , ilwork, nsd  , 'out')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
c
c.... ---------------->   Re-position thin section stack mesh   <-----------------
c
c.... if no thin section stack, skip whole re-positioning method
c
      flag = zero
      if (numts .gt. 0) then
c        print*, 'thin section stack:', numts
c
c.... loop over growth curves
c
        listcounter = 0
        ioffset = 1 ! the ID starts from 1 in phasta
        do nts = 1, numts
          itnv = TStnv(nts) ! number of vertices on this growth curve
          basevID = TSlist(listcounter + 1) + ioffset
          endvID = TSlist(listcounter + itnv) + ioffset  ! last vertex on growth curve
c          print*, 'itnv', itnv
c          print*,'basevID', basevID, x(basevID,:)
c          print*, 'endvID', endvID, x(endvID,:) 
          dist = sqrt( (x(basevID,1) - x(endvID,1))**2
     &                 + (x(basevID,2) - x(endvID,2))**2
     &                 + (x(basevID,3) - x(endvID,3))**2 )
          dist = dist/(itnv-1)        ! dist between each layer (const) 
c
c.... precaution
c
          if (itnv .lt. 2) then
            listcounter = listcounter + itnv
            cycle ! not loop over vertices
          endif
c.... no BC on base, skip this growth curve
          if (ibits(iBC(basevID),14,3) .eq. 0) then
            cycle ! not loop over vertices
          endif
c
c.... prepare other paramteres
c
c          iflt = BLflt(ngc) ! first layer thickness of this growth curve
c          igr  = BLgr(ngc)  ! growth ratio of this growth curve
c
c.... if the first and second vertex on boundary, use current normal
c     this assumes that the growth curve is exposed to a flat boundary
c
          if (mesh2geom .eq. 1) then
            secondvID = TSlist(listcounter + 2) + ioffset
c
            if ( (m2gClsfcn(basevID,1).ne.3) .and.
     &           (m2gClsfcn(secondvID,1).ne.3) ) then
              tsnormal(basevID,:) = x(secondvID,:) - x(basevID,:)
            endif
          endif
c
c... if base face option is set to be interface-only,
c    use current normal if base is not on inteface
c          if (gcBaseOpt .eq. 1) then
c            if (ifFlag(basevID) .ne. 1) then
c              secondvID = BLlist(listcounter + 2) + ioffset
c              gcnormal(basevID,:) = x(secondvID,:) - x(basevID,:)
c            endif
c          endif
c
          tmp  = sqrt( tsnormal(basevID,1) * tsnormal(basevID,1)
     &               + tsnormal(basevID,2) * tsnormal(basevID,2)
     &               + tsnormal(basevID,3) * tsnormal(basevID,3) )
          inormal(:) = tsnormal(basevID,:) / tmp
c
c.... loop over vertices on this growth curve
c
          igrexp = 1.0
          do nv = 2, itnv
            vID = TSlist(listcounter + nv) + ioffset
            vID2= TSlist(listcounter + nv - 1) + ioffset ! the previous one
c            dist = sqrt( (x(vID,1) - x(vID2,1))**2
c     &                 + (x(vID,2) - x(vID2,2))**2
c     &                 + (x(vID,3) - x(vID2,3))**2 ) 
c            print*, 'vID1:', vID, x(vID,:)
c            print*, 'vID2:', vID2, x(vID2,:), dist
            xtmp(vID,:) = xtmp(vID2,:) + dist * inormal(:);
c            xtmp(vID,:) = xtmp(vID2,:) + iflt * inormal(:) * igrexp  ! igr**(nv-2)
c            print*, 'disp before:', disp(vID,:)
            disp(vID,:) = xtmp(vID,:) - x(vID,:)
c            print*, 'xtmp:', xtmp(vID,:)
c            print*, 'x:', x(vID,:)
c            print*, 'inormal:', inormal
c            print*, 'disp after:', disp(vID,:)
c... set flag
            call iBCelas_to_dbl(iBC(basevID), flag(vID))
          enddo ! over this growth curve
c
          listcounter = listcounter + itnv ! update counter
c
c.... end loop growth curves
c
        enddo
c
c.... end if growth curve
c
      endif ! end if numTS greater than 0
c
c.... -------------------->     communication     <--------------------
c
c.... combine disp + flag into dispCommu
c
      do i = 1, numnp
        dispCommu(i,1:3) = disp(i,:)
        dispCommu(i,4)   = flag(i)
      enddo
c
c.... turn on commu flag
c
      layerCommuFlag = 1
c
c.... commu
c
      if (numpe > 1) then
        call commu (dispCommu  , ilwork, 4  , 'in ')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call commu (dispCommu  , ilwork, 4  , 'out')
      endif
c
c.... turn off commu flag
c
      layerCommuFlag = 0
c
c.... separate dispCommu into disp and iBC; update BC
c
      do i = 1, numnp
        disp(i,:) = dispCommu(i,1:3)
        flag(i)   = dispCommu(i,4)
        call dbl_to_iBCelas(flag(i), iBC(i))
        call setTSbc(disp(i,:), iBC(i), BC(i, ndof+2:ndof+5))
      enddo
c
c.... end re-position thin section layered mesh




c.... -------------------->   interior elements   <--------------------
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        elasres = zero
        errorcount = zero
        if (lhs. eq. 1)    elaslhsK  = zero
        if (iprec .ne. 0)  elasBDiag = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk            ! used in timeseries
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
          ndofelas = nshl * nelas
c
          allocate (Estiff(npro,ndofelas,ndofelas))
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          allocate (meshqbk(npro))
          allocate (meshVbk(npro))
c
          Estiff = zero
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
c.... Shape measure. Calculate the shape quality
c
          call shpMeasure(x, mien(iblk)%p, tmpshp, tmpshgl,
     &                    meshqbk,  meshVbk,  errorcount )
c
c.... map local element to global
          do i = 1, npro
            meshq(mieMap(iblk)%p(i)) = meshqbk(i)
          enddo
c
          call AsIGMRElas (x,             disp,
     &                     tmpshp,        tmpshgl,
     &                     mien(iblk)%p,  elasres,
     &                     elasBDiag,     Estiff,
     &                     meshqbk,       meshVbk   )
c
c.... satisfy the BCs on the implicit LHS
c
          call bc3LHSElas (iBC, BC(:, ndof+2:ndof+5),
     &                     mien(iblk)%p, Estiff)
c
c.... Fill-up the global sparse LHS mass matrix
c
          call fillsparseElas( mien(iblk)%p, Estiff,
     &                         elaslhsK, row, col)
c
          deallocate ( Estiff )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
          deallocate ( meshqbk )
          deallocate ( meshVbk )
c
c.... end of interior element loop
c
       enddo
c
      if (errorcount(1) .lt. 0 .or. errorcount(2) .lt. 0) then
        write(*,*) errorcount(1), " elements Meshq larger than one; ",
     &             errorcount(2), " elements Meshq smaller than zero."
      endif
c
      if(iabc==1) then               ! are there any axisym BCs
          call rotabc(elasres(1,2), iBC,  'in ')
      endif
c
c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
          call commu (elasres  , ilwork, nelas  , 'in ')
c
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c
          if (iprec .ne. 0) then
             call commu (elasBDiag,    ilwork,
     &                   nelas*nelas,  'in ')
          endif
      endif
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResElas (iBC,     BC(:, ndof+2:ndof+5),
     &                 elasres, iper,    ilwork)
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
      if (iprec .ne. 0) then
         call bc3BDgElas (iBC,       BC(:, ndof+2:ndof+5),
     &                    elasBDiag, iper,    ilwork)
      endif
c
c.... return
c
      return
      end

