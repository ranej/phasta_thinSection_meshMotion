        subroutine updateSnapSurfBC ( x,       disp_snap,
     &                                iBC,     BC)
c
         use m2gfields ! read m2g fields
         use pointer_data
         use core_snap
c
        include "common.h"
c
        real*8    x(numnp,nsd),    disp_snap(numnp,nsd)
        dimension iBC(nshg),       BC(nshg,4)
        integer   face_snap(numnp)
        integer   answer
c
        if ( snapSurfFlag .eq. 1 ) then ! double check
c
        face_snap = zero
        if (snapNumFaceTags .gt. 0) then ! user specify face tags
          do i = 1, nshg
            if (ibits(iBC(i),14,3) .eq. 7) then
              cycle
            endif
            if (m2gClsfcn(i,1) .eq. 3 .or. m2gClsfcn(i,1) .eq. 0) then ! region or vertex
              cycle
            else if (m2gClsfcn(i,1) .eq. 2) then ! face
              do j = 1, snapNumFaceTags
                if (m2gClsfcn(i,2) .eq. snapFaceTags(j)) then
                  face_snap(i) = 1
                  exit
                endif
              enddo
            else if (m2gClsfcn(i,1) .eq. 1) then ! edge
              answer = zero
              do j = 1, snapNumFaceTags
                call core_is_in_closure(m2gClsfcn(i,1), m2gClsfcn(i,2),
     &                                 2,              snapFaceTags(j),
     &                                 answer)
                if (answer .ne. 0) then
                  face_snap(i) = 1
                  exit
                endif
              enddo
            endif ! end check dimension
          enddo ! end loop nshg
        else ! auto-detect snapping face
          do i = 1, nshg
            if(m2gClsfcn(i,1) .lt. 3 .and. m2gClsfcn(i,1) .gt. 0
     &         .and. ibits(iBC(i),14,3) .ne. 7) then
              face_snap(i) = 1
            endif
          enddo
        endif
c
        call resetSnapBC(x, disp_snap, iBC, BC, face_snap)
c
        endif ! end check snap flag
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine resetSnapBC ( x,       disp_snap,
     &                           iBC,     BC,        face_snap)
c
        use iso_c_binding
        use core_snap
c
        include "common.h"
c
        real*8    x(numnp,nsd),  disp_snap(numnp,nsd)
        dimension iBC(nshg),     BC(nshg, 4)
        integer   face_snap(numnp)
        real*8    mag,                rad
        integer   i,    j
        real*8    x_tmp_1(numnp), x_tmp_2(numnp), x_tmp_3(numnp)
        real*8    x_crt_1(numnp), x_crt_2(numnp), x_crt_3(numnp)
c
        x_tmp_1 = x(:,1) + disp_snap(:,1)
        x_tmp_2 = x(:,2) + disp_snap(:,2)
        x_tmp_3 = x(:,3) + disp_snap(:,3)
        if(elasFDC .gt. 0 .and. numrbs .eq. 0) then
          write(*,*), 'Calling core_get_pos_on_surf_discrete'
          call core_get_pos_on_surf_discrete (elasFDC, x_tmp_1, x_tmp_2,
     &             x_tmp_3, numnp, face_snap, x_crt_1, x_crt_2, x_crt_3)
        else 
          call core_get_pos_on_surf (x_tmp_1, x_tmp_2, x_tmp_3, numnp,
     &                             face_snap, x_crt_1, x_crt_2, x_crt_3)
        endif
c       
        do i = 1, nshg
c... if face_snap = 1
          if (face_snap(i) .ne. 0) then
c... update iBC and BC
            iBC(i) = ibset(iBC(i), 14)
            iBC(i) = ibset(iBC(i), 15)
            iBC(i) = ibset(iBC(i), 16)
            BC(i,1)= x_crt_1(i) - x(i,1)
            BC(i,2)= x_crt_2(i) - x(i,2)
            BC(i,3)= x_crt_3(i) - x(i,3)
          endif ! if face_snap is on
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine updateTDComp1BC ( x,       shpb,
     &                               iBC,     BC)
c
         use pointer_data
c
        include "common.h"
c
        real*8    x(numnp,nsd)
        dimension iBC(nshg),       BC(nshg,4),
     &            shpb(MAXTOP,maxsh,MAXQPT)
        integer   surfID_TDComp1(nshg)
        real*8    normal(nshg, nsd)
c
        if ( timeDepComp1Flag .eq. 1) then ! double check
c
        write(*,*) "Modify comp1_elas BC on surfID",timeDepComp1ID
c
        surfID_TDComp1 = zero
        normal = zero
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
          call calc_TDComp1(x,  shpb, mienb(iblk)%p,  miBCB(iblk)%p,
     &                      surfID_TDComp1, normal)
c
        enddo boundary_blocks ! end loop the boundary elements
c
        call resetTDComp1BC(x, surfID_TDComp1, normal, iBC, BC)
c
        endif ! end check TDComp1 flag
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine calc_TDComp1 (x,    shpb,     ienb,    iBCB,
     &                           surfID_TDComp1, normal)
c
        include "common.h"
c
        integer   surfID_TDComp1(nshg)
        dimension x(numnp,nsd),       shpb(nshl,ngaussb),
     &            iBCB(npro,ndiBCB),  ienb(npro,nshl)
        real*8    normal(nshg, nsd)
        integer   calc_factor(npro)
c
c... collect surf ID = timeDepComp1ID
c
        calc_factor(:) = 0
        do inode = 1, npro
          if (iBCB(inode,2) .eq. timeDepComp1ID) then
            calc_factor(inode) = 1
            do i = 1, nenbl ! only loop over vtx on the boundary
              surfID_TDComp1(ienb(inode,i)) = iBCB(inode,2)
            enddo
          endif
        enddo
c
c.... assemble the normal vector
c
        call calc_normal(x, shpb, calc_factor, ienb, iBCB, normal)
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine resetTDComp1BC (x, surfID_TDComp1, normal, iBC, BC)
c
        include "common.h" ! access timeDepComp1Mag
c
        real*8    x(numnp,nsd)
        integer   surfID_TDComp1(nshg)
        real*8    normal(nshg, nsd),    tmp
        dimension iBC(nshg),     BC(nshg, 4)
        integer   maxDir(1)
c
        do i = 1, nshg
c... if surf ID is timeDepComp1ID
          if (surfID_TDComp1(i) .eq. timeDepComp1ID) then
c... if BC code is comp1 (iBC = 1 or 2 or 4)
            if ((ibits(iBC(i),14,3) .eq. 1) .or.
     &          (ibits(iBC(i),14,3) .eq. 2) .or.
     &          (ibits(iBC(i),14,3) .eq. 4)) then
c... calculate normal and normalize it
              tmp = normal(i,1)*normal(i,1)
     &            + normal(i,2)*normal(i,2)
     &            + normal(i,3)*normal(i,3)
              normal(i,:) = normal(i,:) / sqrt(tmp)
              maxDir = maxloc(abs(normal(i,:)))
              select case (maxDir(1))
c... if x of normal is the max
              case (1)
              iBC(i) = ibset(iBC(i), 14)
              iBC(i) = ibclr(iBC(i), 15)
              iBC(i) = ibclr(iBC(i), 16)
              BC(i,1)= timeDepComp1Mag(i) / normal(i,1)
              BC(i,2)=        normal(i,2) / normal(i,1)
              BC(i,3)=        normal(i,3) / normal(i,1)
c... if y of normal is the max
              case (2)
              iBC(i) = ibclr(iBC(i), 14)
              iBC(i) = ibset(iBC(i), 15)
              iBC(i) = ibclr(iBC(i), 16)
              BC(i,1)= timeDepComp1Mag(i) / normal(i,2)
              BC(i,2)=        normal(i,1) / normal(i,2)
              BC(i,3)=        normal(i,3) / normal(i,2)
c... if z of normal is the max
              case (3)
              iBC(i) = ibclr(iBC(i), 14)
              iBC(i) = ibclr(iBC(i), 15)
              iBC(i) = ibset(iBC(i), 16)
              BC(i,1)= timeDepComp1Mag(i) / normal(i,3)
              BC(i,2)=        normal(i,1) / normal(i,3)
              BC(i,3)=        normal(i,2) / normal(i,3)
              end select
            endif ! if BC code is comp1
          endif ! if equal to timeDepComp1ID
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c

