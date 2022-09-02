       subroutine timeDependBCElas(x, iBC, BC, BC_flow, umeshold)
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      include "common.h"
      include "mpif.h"

      real*8    x(numnp,nsd)
      real*8    umeshold(numnp,nsd)
      dimension iBC(nshg),        BC(nshg,4), BC_flow(nshg,4)
c

      if (numrbs .gt. 0)  then
        call rigidBodyBCElas(x, iBC, BC, BC_flow)
      endif
c
      if (elasFDC .gt. 0 .or. elasSICC .gt. 0) then
        call prescribedBCElas(x, iBC, BC, BC_flow, umeshold)
      endif
c
      return
      end
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      subroutine rigidBodyBCElas(x,  iBC,  BC,  BC_flow)
c
        use rigidBodyReadData
        use rigidBodyForce
c
        include "common.h"
        include "mpif.h"
c
        real*8    x(numnp,nsd)
c        real*8    umeshold(numnp,nsd)
        dimension iBC(nshg), BC(nshg,4), BC_flow(nshg,4)
        real*8    t
        real*8, dimension(nsd) :: ra, rc, r_disp
c
        do j = 1,numrbs
          if (rbsMM(j).eq.3 ) then
            call calc_rbMotionPrescribed
          else
            call calc_rbMotion
          endif
        enddo
           
c
c.... loop over mesh vertices
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (rbFlags(i) .gt. 0) ) then
c.... calculate the displacement due to rotation
            t = rbAng(rbFlags(i))/180.0*pi
            ra = rb_prop(rbFlags(i),5:7)
            rc = rb_prop(rbFlags(i),8:10)
            r_disp(1) =
     &               (cos(t)+ra(1)*ra(1)*(1-cos(t))-1)*(x(i,1)-rc(1))
     &         + (ra(1)*ra(2)*(1-cos(t))-ra(3)*sin(t))*(x(i,2)-rc(2))
     &         + (ra(1)*ra(3)*(1-cos(t))+ra(2)*sin(t))*(x(i,3)-rc(3))
            r_disp(2) =
     &           (ra(2)*ra(1)*(1-cos(t))+ra(3)*sin(t))*(x(i,1)-rc(1))
     &             + (cos(t)+ra(2)*ra(2)*(1-cos(t))-1)*(x(i,2)-rc(2))
     &         + (ra(2)*ra(3)*(1-cos(t))-ra(1)*sin(t))*(x(i,3)-rc(3))
            r_disp(3) =
     &           (ra(3)*ra(1)*(1-cos(t))-ra(2)*sin(t))*(x(i,1)-rc(1))
     &         + (ra(3)*ra(2)*(1-cos(t))+ra(1)*sin(t))*(x(i,2)-rc(2))
     &             + (cos(t)+ra(3)*ra(3)*(1-cos(t))-1)*(x(i,3)-rc(3))
c.... set corresponding mesh elas BC
            BC(i,1:3) = rbDisp(rbFlags(i),1:3) + r_disp(1:3)
c.... update flow BC
            BC_flow(i,1:3) = rbVel(rbFlags(i),1:3) + r_disp(1:3)/Delt(1)
          endif
        enddo ! end loop numnp
c
      return
      end
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
       subroutine prescribedBCElas(x, iBC, BC, BC_flow, umeshold)
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      use m2gfields ! read m2g fields
      use core_snap
      use interfaceflag
      use core_rigid_body
      use rotatingBandForce
c
      include "common.h"
      include "mpif.h"
c
c
      real*8    x(numnp,nsd)
      real*8    disp(numnp,nsd+1)
      real*8    umeshold(numnp,nsd)
      dimension iBC(nshg),        BC(nshg,4), BC_flow(nshg,4)
      integer   casenumber
      real*8    acc
      real*8    totalForce(3),    objMass
      dimension Forin(3),         Forout(3)
c
c.... dynamic origin, x translation, rotation frequence
      real*8    dyn_org,   xtsl,  rotf
      real*8    dyn_lnt,   shrk,  shrkfactor
      integer   answer
      real*8    cent(numrbs, 3)
c

c.... precribed BC for rifling case
      real*8    ang, x_center(nshg,3), rotx(nshg,3)
      integer   rot_band_ftag, ftag1, ftag2, ftag3, ftag4
      integer   etag1, etag2, etag3, etag4
      integer   maxDir(1)
      real*8    pnormal(nshg,3), inormal(nshg,3)
      real*8    m2gD(nshg,3)
      real*8    cos_theta, sin_theta
      integer   counter
      real*8    Forcetmp(numnp,3)

      casenumber = 0
      if (elasFDC .gt. 0) then
        if (myrank .eq. master) then
          write(*,*) "use Force-driven case:", elasFDC
        endif
        casenumber = elasFDC
      endif

      if (elasSICC .gt. 0) then
        if (myrank .eq. master) then
          write(*,*) "use set-in-code case:", elasSICC
        endif
        casenumber = elasSICC
      endif
c
c.... Update BC value based on geom and iBC
c
      if ( casenumber .eq. 1 ) then
        acc = 60000.0
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)-0.8) .lt. 0.79) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 0.004225)) then
c.... 0.004225 = 0.065^2
            BC(i,1)   = acc * lstep * Delt(1) * Delt(1)
            BC(i,2)   = zero
            BC(i,3)   = zero
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 1
c
c.... Update BC value based on total force on the object
c
      if ( casenumber .eq. 2 ) then
c.... debugging {
c        do j = 1,numrbs
c          call core_get_centroid(rbsTags(j), cent(j,:))
c        enddo
c.... debugging }
        totalForce(:) = zero
        objMass = 15.0 ! kg
c.... bcast force to all processors
        if (numpe > 1) then
          do j = 1, nsrfCM
            Forin  = (/ Force(1,j), Force(2,j), Force(3,j) /)
            call MPI_BCAST (Forin(1), 3, MPI_DOUBLE_PRECISION,
     &                      master,      MPI_COMM_WORLD,ierr)
            Force(1:3,j) = Forin(1:3)
          enddo
        endif
c.... collect total force
        do j = 1,nsrfCM
          totalForce(:) = totalForce(:) + Force(:,j)
        enddo ! end collect total force

        acc = totalForce(1) / objMass
        if (myrank .eq. master) then
          write(*,*) "projectile acceleration is:", acc
        endif

        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)-0.8) .lt. 0.79) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 0.004225)) then
c.... 0.004225 = 0.065^2
            BC(i,1)   = (umeshold(i,1) + acc * Delt(1)) * Delt(1)
            BC(i,2)   = zero
            BC(i,3)   = zero

c.... >>> hard-coding update flow velocity boundry condition on mortar surfaces
            BC_flow(i,1:3) = BC(i,1:3) / Delt(1)
c.... <<< hard-coding

          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 2
c
c.... Update BC value based on geom and iBC
c.... rotate + trans + shrink
c
      if ( casenumber .eq. 3 ) then
        do i = 1,numnp
c          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
          if ( (ibits(iBC(i),3,3) .eq. 7) .and.
     &         (abs(x(i,1)) .lt. 0.2) .and.
     &         (abs(x(i,2)) .lt. 0.2) .and.
     &         (abs(x(i,3)) .lt. 0.2) ) then
c
            disp(i,1) = -0.005 * (x(i,1)) !+ 0.02
     &                + (x(i,1)) * (cos(pi/600) - 1.0)
     &                - (x(i,2)) *  sin(pi/600)
            disp(i,2) = !-0.005 * (x(i,2)) !+ 0.01
     &                + (x(i,2)) * (cos(pi/600) - 1.0)
     &                + (x(i,1)) *  sin(pi/600)
            disp(i,3) = -0.005 * x(i,3)
c
            BC(i,:)   = disp(i,:)
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 3
c
c.... Update BC value based on geom and iBC
c.... twist one end of a cylinder
c
      if ( casenumber .eq. 4 ) then
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 0.501) .and.
     &         (x(i,1) .gt. 0.499) ) then
c
            disp(i,1) = 0.0
            disp(i,2) = x(i,2) * (cos(pi/100) - 1.0)
     &                - x(i,3) *  sin(pi/100)
            disp(i,3) = x(i,3) * (cos(pi/100) - 1.0)
     &                + x(i,2) *  sin(pi/100)
c
            BC(i,:)   = disp(i,:)
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 4
c
c.... test case 5
c.... original mov+rot+shrk
c
      if ( casenumber .eq. 5 ) then
        if (lstep .eq. 0) then
        dyn_org    = 2.625e-3
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 6.65e-3) .and. (x(i,1) .gt. -1.4e-3) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 2.89e-6) ) then
            if( x(i,2) .ge. 0.6e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! top tail
                disp(i,1) = -0.05 * (x(i,1)-0.5e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.05 * x(i,3)
              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! top head
                disp(i,1) = -0.05 * (x(i,1)-4.75e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.05 * x(i,3)
              else ! top middle
                disp(i,1) = -0.05 * (x(i,1)-2.625e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.05 * x(i,3)
              endif ! end top--switch head middle tail
            else if( x(i,2) .le. -0.6e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! bottom tail
                disp(i,1) = -0.05 * (x(i,1)-0.5e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.05 * x(i,3)
              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! bottom head
                disp(i,1) = -0.05 * (x(i,1)-4.75e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.05 * x(i,3)
              else ! bottom middle
                disp(i,1) = -0.05 * (x(i,1)-2.625e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.05 * x(i,3)
              endif ! end bottom--switch head middle tail
            else ! middle
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! middle tail
                disp(i,1) = -0.05 * (x(i,1)-0.5e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
                disp(i,3) = -0.05 * x(i,3)
              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! middle head
                disp(i,1) = -0.05 * (x(i,1)-4.75e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
                disp(i,3) = -0.05 * x(i,3)
              else ! middle middle
                disp(i,1) = -0.05 * (x(i,1)-2.625e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
                disp(i,3) = -0.05 * x(i,3)
              endif ! end middle--switch head middle tail
            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
        endif ! first time step
      endif ! end if case 5
c
c.... end test case 5
c
c
c.... test case 6
c.... original mov+rot+shrk
c
      if ( casenumber .eq. 6 ) then
        xtsl     = 2.0000000000000e-4
        dyn_org  = 2.6250000000000e-3 + DBLE(lstep) * xtsl
        if (myrank .eq. master) then
          write(*,*) "current lstep:", lstep, "dyn_org:", dyn_org
        endif
        rotf     = 90.0
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 6.65e-3) .and. (x(i,1) .gt. -1.4e-3) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 2.89e-6) ) then
            if( x(i,2) .ge. 0.6e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! top tail
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) + xtsl
     &                    + (x(i,1)-dyn_org+2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org+2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! top head
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) + xtsl
     &                    + (x(i,1)-dyn_org-2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org-2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              else ! top middle
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) + xtsl
     &                    + (x(i,1)-dyn_org) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -0.6e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! bottom tail
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) + xtsl
     &                    - (x(i,1)-dyn_org+2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org+2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! bottom head
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) + xtsl
     &                    - (x(i,1)-dyn_org-2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org-2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              else ! bottom middle
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) + xtsl
     &                    - (x(i,1)-dyn_org) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! middle tail
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) + xtsl
c     &                    + (x(i,1)-0.5e-3) * (cos(pi/rotf) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/rotf) - 1.0)
c     &                    + (x(i,1)-0.5e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! middle head
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) + xtsl
c     &                    + (x(i,1)-4.75e-3) * (cos(pi/rotf) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/rotf) - 1.0)
c     &                    + (x(i,1)-4.75e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              else ! middle middle
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) + xtsl
c     &                    + (x(i,1)-2.625e-3) * (cos(pi/rotf) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/rotf) - 1.0)
c     &                    + (x(i,1)-2.625e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * x(i,3)

              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 6
c
c.... end test case 6
c
c
c.... test case 7
c
      if ( casenumber .eq. 7 ) then
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &        (x(i,1) .lt. 0.029) .and.
     &        (x(i,1) .gt.-0.009) ) then
c
            disp(i,1) = 0.001
            disp(i,2) = 0.0
            disp(i,3) = 0.0
c
            BC(i,:)   = disp(i,:)
          endif
        enddo ! end loop numnp
      endif ! end if case 7
c
c.... end test case 7
c
c
c.... test case 8
c.... mov+rot+shrk non-uniformly
c
      if ( casenumber .eq. 8 ) then
        xtsl     = 2.0000000000000e-4
        dyn_org  = 2.6250000000000e-3 + DBLE(lstep) * xtsl
        dyn_lnt  = 2.0000000000000e-3 * 0.9**lstep
        if (myrank .eq. master) then
          write(*,*) "current lstep:", lstep, "dyn_org:", dyn_org
        endif
        rotf     = 90.0
        shrkfactor = 2.0
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 6.65e-3) .and. (x(i,1) .gt. -1.4e-3) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 2.89e-6) ) then
            if( x(i,2) .ge. 0.6e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! top tail
                shrk = 0.5 - abs(x(i,1)-dyn_org+2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) + xtsl
     &                    + (x(i,1)-dyn_org+2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org+2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! top head
                shrk = 0.5 - abs(x(i,1)-dyn_org-2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) + xtsl
     &                    + (x(i,1)-dyn_org-2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org-2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else ! top middle
                shrk = 0.5 - abs(x(i,1)-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) + xtsl
     &                    + (x(i,1)-dyn_org) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org) *  sin(pi/rotf)
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -0.6e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! bottom tail
                shrk = 0.5 - abs(x(i,1)-dyn_org+2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) + xtsl
     &                    - (x(i,1)-dyn_org+2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org+2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! bottom head
                shrk = 0.5 - abs(x(i,1)-dyn_org-2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) + xtsl
     &                    - (x(i,1)-dyn_org-2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org-2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else ! bottom middle
                shrk = 0.5 - abs(x(i,1)-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) + xtsl
     &                    - (x(i,1)-dyn_org) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org) *  sin(pi/rotf)
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! middle tail
                shrk = 0.5 - abs(x(i,1)-dyn_org+2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) + xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)) ! + 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! middle head
                shrk = 0.5 - abs(x(i,1)-dyn_org-2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) + xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)) ! + 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else ! middle middle
                shrk = 0.5 - abs(x(i,1)-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) + xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)) ! + 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 8
c
c.... end test case 8
c
c.... test case 9
c.... Update mesh2geom to identify vertex
c.... twist one end of a cylinder
c
      if ( casenumber .eq. 9 ) then
        do i = 1,numnp
c
          call core_is_in_closure(m2gClsfcn(i,1), m2gClsfcn(i,2),
     &                            3,              108,
     &                            answer)
c
          if (answer .ne. 0) then
c            disp(i,1) = 0.0015
c            disp(i,2) = 0.0
c            disp(i,3) = 0.0
            disp(i,1) = 0.0015
            disp(i,2) = x(i,2) * (cos(pi/100) - 1.0)
     &                - x(i,3) *  sin(pi/100)
            disp(i,3) = x(i,3) * (cos(pi/100) - 1.0)
     &                + x(i,2) *  sin(pi/100)
c
            BC(i,:)   = disp(i,:)
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 9
c
c.... test case 10
c.... drive interface toward -x direction in stefan case
c
      if ( casenumber .eq. 10 ) then
        do i = 1,numnp
          if ( ifFlag(i) .eq. 1 ) then ! interface node
            disp(i,1) = -1.0 * Delt(1)
            disp(i,2) =  0.0 * Delt(1)
            disp(i,3) =  0.0 * Delt(1)
            BC(i,:) = disp(i,:)
          endif
        enddo ! end loop numnp
      endif  ! end case 10
c
c.... test case 11
c.... prescribe 6-grain displacement
c     mimic real simulation deformation
c     non-uniformly shrink; no rotation
c
      if ( casenumber .eq. 11 ) then
        xtsl    = 5.3333333333333e-3 * 0.98**lstep
        dyn_org = (5.3333333333333e-3 - 5.3333333333333e-3*0.98**lstep)/
     &            (1.0-0.98)
        dyn_lnt = 2.0000000000000e1 * 0.98**lstep
        shrkfactor = 2.0
c
        if (myrank .eq. master) then
          write(*,*) "current lstep:", lstep, "xtsl:", xtsl
        endif
c
        do i = 1,numnp
          if ( ifFlag(i) .eq. 1 ) then ! interface node
            if( x(i,2) .gt. 0.0 ) then ! top
              if ( x(i,1) .le. -1.25 ) then ! top tail
                shrk = 0.5 - abs(x(i,1)+2.5+dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)+2.5+dyn_org) - xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)-0.75) !+ 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else if ( x(i,1) .ge. 1.25 ) then ! top head
                shrk = 0.5 - abs(x(i,1)-2.5-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)-2.5-dyn_org) + xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)-0.75) !+ 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else ! top middle
                shrk = 0.5 - abs(x(i,1))/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1))
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)-0.75) !+ 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              endif ! end top--switch head middle tail
            else ! bottom
              if ( x(i,1) .le. -1.25 ) then ! bottom tail
                shrk = 0.5 - abs(x(i,1)+2.5+dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)+2.5+dyn_org) - xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)+0.75) !- 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else if ( x(i,1) .ge. 1.25 ) then ! bottom head
                shrk = 0.5 - abs(x(i,1)-2.5-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)-2.5-dyn_org) + xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)+0.75) !- 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else ! bottom middle
                shrk = 0.5 - abs(x(i,1))/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1))
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)+0.75) !- 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              endif ! end bottom--switch head middle tail
            endif ! end if switch top bottom
            BC(i,1:3)   = disp(i,1:3)
          endif
        enddo ! end loop numnp
      endif  ! end case 11
c
c.... test case 12
c.... shrink non-uniformlyic
      if ( casenumber .eq. 12 ) then
        xtsl     = 2.0000000000000e-4
        dyn_org  = 2.6250000000000e-3 ! + DBLE(lstep) * xtsl
        dyn_lnt  = 2.0000000000000e-3 * 0.9**lstep
        if (myrank .eq. master) then
          write(*,*) "current lstep:", lstep, "dyn_org:", dyn_org
        endif
        rotf     = 90.0
        shrkfactor = 2.0
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 6.65e-3) .and. (x(i,1) .gt. -1.4e-3) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 2.89e-6) ) then
            if( x(i,2) .ge. 0.6e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! top tail
                shrk = 0.5 - abs(x(i,1)-dyn_org+2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! top head
                shrk = 0.5 - abs(x(i,1)-dyn_org-2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else ! top middle
                shrk = 0.5 - abs(x(i,1)-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -0.6e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! bottom tail
                shrk = 0.5 - abs(x(i,1)-dyn_org+2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! bottom head
                shrk = 0.5 - abs(x(i,1)-dyn_org-2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else ! bottom middle
                shrk = 0.5 - abs(x(i,1)-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! middle tail
                shrk = 0.5 - abs(x(i,1)-dyn_org+2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org+2.125e-3) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)) ! + 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! middle head
                shrk = 0.5 - abs(x(i,1)-dyn_org-2.125e-3)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org-2.125e-3) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)) ! + 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              else ! middle middle
                shrk = 0.5 - abs(x(i,1)-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.1 * (x(i,1)-dyn_org) !+ xtsl
                disp(i,2) = -0.1 * (1.0 + shrk) * (x(i,2)) ! + 0.01
                disp(i,3) = -0.1 * (1.0 + shrk) * x(i,3)

              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 12
c
c.... end test case 12
c
c.... test case 13
c.... prescribe 6-grain displacement
c     non-uniformly shrink; no rotation
c
      if ( casenumber .eq. 13 ) then
        xtsl    = 5.3333333333333e-3 * 0.98**lstep
        dyn_org = (5.3333333333333e-3 - 5.3333333333333e-3*0.98**lstep)/
     &            (1.0-0.98)
        dyn_lnt = 2.0000000000000e1 * 0.85**lstep
        shrkfactor = 2.0
c
        if (myrank .eq. master) then
          write(*,*) "current lstep:", lstep, "xtsl:", xtsl
        endif
c
        do i = 1,numnp
          if ( ifFlag(i) .eq. 1 ) then ! interface node
            if( x(i,2) .gt. 0.0 ) then ! top
              if ( x(i,1) .le. -1.25 ) then ! top tail
                shrk = 0.5 - abs(x(i,1)+2.5+dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)+2.5+dyn_org) - xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)-0.75) !+ 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else if ( x(i,1) .ge. 1.25 ) then ! top head
                shrk = 0.5 - abs(x(i,1)-2.5-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)-2.5-dyn_org) + xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)-0.75) !+ 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else ! top middle
                shrk = 0.5 - abs(x(i,1))/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1))
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)-0.75) !+ 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              endif ! end top--switch head middle tail
            else ! bottom
              if ( x(i,1) .le. -1.25 ) then ! bottom tail
                shrk = 0.5 - abs(x(i,1)+2.5+dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)+2.5+dyn_org) - xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)+0.75) !- 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else if ( x(i,1) .ge. 1.25 ) then ! bottom head
                shrk = 0.5 - abs(x(i,1)-2.5-dyn_org)/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1)-2.5-dyn_org) + xtsl
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)+0.75) !- 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              else ! bottom middle
                shrk = 0.5 - abs(x(i,1))/dyn_lnt
                shrk = shrkfactor * shrk
c
                disp(i,1) = -0.02 * (x(i,1))
                disp(i,2) = -0.02 * (1.0 + shrk) * (x(i,2)+0.75) !- 0.01
                disp(i,3) = -0.02 * (1.0 + shrk) * x(i,3)
              endif ! end bottom--switch head middle tail
            endif ! end if switch top bottom
            BC(i,1:3)   = disp(i,1:3)
          endif
        enddo ! end loop numnp
      endif  ! end case 13
c
c.... test case 14
c

      if ( casenumber .eq. 14 ) then
        do i = 1,numnp
          if  (ibits(iBC(i),14,3) .eq. 7) then
c
            disp(i,1) = 0.0
            disp(i,2) = -0.01/5
            disp(i,3) = 0.0
c
            BC(i,:)   = disp(i,:)
          endif
        enddo ! end loop numnp
      endif ! end if case 14
c
c.... end test case 14

c
c.... test case 15
c
      if ( casenumber .eq. 15 ) then
        do i = 1,numnp
          if  ((ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,3) .gt. 7.5)) then

            disp(i,1) = 0
            disp(i,2) = 0
            disp(i,3) = 0.1
c
            BC(i,:)   = disp(i,:)
          endif
        enddo ! end loop numnp
      endif ! end if case 14
c
c.... end test case 15

c
c.... test case 16
c


      if ( casenumber .eq. 16 ) then

        rot_band_ftag = 511
        ftag1 = 78
        ftag2 = 82
        etag1 = 324
        etag2 = 417
        etag3 = 331
        etag4 = 402
        ftag3 = 521
        ftag4 = 502
        
                
c For 1 degree rotation
        ang = 1
        ang = ang*pi/180

        cos_theta = 9.9984769515639127e-01
        sin_theta = 1.7452406437283512e-02

c        cos_theta = 9.9939082701909576e-01
c        sin_theta = 3.4899496702500969e-02


        do i = 1, numnp

          x_center(i,1) = x(i,1) 
          x_center(i,2) = 0.0077
          x_center(i,3) = -4.9403542382592772e-02

          if (m2gClsfcn(i,1) .ne. 3) then
            call core_is_in_closure(m2gClsfcn(i,1), m2gClsfcn(i,2),
     &                            2,              rot_band_ftag,
     &                            answer)
c
            if (answer .ne. 0) then
              rotx(i,1) = x(i,1)
              rotx(i,2) = x_center(i,2) + (x(i,2) - x_center(i,2))
     &          *cos_theta - (x(i,3) - x_center(i,3))*sin_theta
              rotx(i,3) = x_center(i,3) + (x(i,2) - x_center(i,2))
     &          *sin_theta + (x(i,3) - x_center(i,3))*cos_theta
              iBC(i) = ibset(iBC(i), 14)
              iBC(i) = ibset(iBC(i), 15)
              iBC(i) = ibset(iBC(i), 16)
              BC(i,1) = rotx(i,1) - x(i,1)
              BC(i,2) = rotx(i,2) - x(i,2)
              BC(i,3) = rotx(i,3) - x(i,3)

            endif
          endif




c Prescribe BC on two curved faces and 4 edges
          if( (m2gClsfcn(i,1) .eq. 2 .and. m2gClsfcn(i,2). eq. ftag1) .or. 
     &      (m2gClsfcn(i,1) .eq. 2 .and. m2gClsfcn(i,2). eq. ftag2) .or.
     &      (m2gClsfcn(i,1) .eq. 1 .and. m2gClsfcn(i,2). eq. etag1) .or.
     &      (m2gClsfcn(i,1) .eq. 1 .and. m2gClsfcn(i,2). eq. etag2) .or.
     &      (m2gClsfcn(i,1) .eq. 1 .and. m2gClsfcn(i,2). eq. etag3) .or.
     &      (m2gClsfcn(i,1) .eq. 1 .and. m2gClsfcn(i,2). eq. etag4) )
     &      then
            
            pnormal(i,:) = x(i,:) - x_center(i,:)
            tmp  = sqrt( pnormal(i,1) * pnormal(i,1)
     &                 + pnormal(i,2) * pnormal(i,2)
     &                 + pnormal(i,3) * pnormal(i,3) )
            inormal(i,:) = pnormal(i,:) / (tmp) 
            
            maxDir  = maxloc(abs(inormal(i,:)))
            select case(maxDir(1))
c... if x of normal is the max
            case (1)
            iBC(i) = ibset(iBC(i), 14)
            iBC(i) = ibclr(iBC(i), 15)
            iBC(i) = ibclr(iBC(i), 16)
            BC(i,1)= 0d0 / pnormal(i,1)
            BC(i,2)= pnormal(i,2) / pnormal(i,1)
            BC(i,3)= pnormal(i,3) / pnormal(i,1)
c... if y of normal is the max
            case (2)
            if (m2gClsfcn(i,1) .eq. 2) then
              iBC(i) = ibclr(iBC(i), 14)
              iBC(i) = ibset(iBC(i), 15)
              iBC(i) = ibclr(iBC(i), 16)
              BC(i,1)= 0d0 / pnormal(i,2)
              BC(i,2)= pnormal(i,1) / pnormal(i,2)
              BC(i,3)= pnormal(i,3) / pnormal(i,2)
            else
              iBC(i) = ibset(iBC(i), 14)
              iBC(i) = ibset(iBC(i), 15)
              iBC(i) = ibclr(iBC(i), 16)
              BC(i,1)= 0d0
              BC(i,2)= 0d0
              BC(i,3)= 0d0 / pnormal(i,2)
              BC(i,4)= pnormal(i,3) / pnormal(i,2)              
            endif    
c... if z of normal is the max
            case (3)
            if (m2gClsfcn(i,1) .eq. 2) then
              iBC(i) = ibclr(iBC(i), 14)
              iBC(i) = ibclr(iBC(i), 15)
              iBC(i) = ibset(iBC(i), 16)
              BC(i,1)= 0d0 / pnormal(i,3)
              BC(i,2)= pnormal(i,1) / pnormal(i,3)
              BC(i,3)= pnormal(i,2) / pnormal(i,3)
            else
              iBC(i) = ibset(iBC(i), 14)
              iBC(i) = ibclr(iBC(i), 15)
              iBC(i) = ibset(iBC(i), 16)
              BC(i,1)= 0d0
              BC(i,2)= 0d0
              BC(i,3)= 0d0 / pnormal(i,2)
              BC(i,4)= pnormal(i,2) / pnormal(i,3)              
            endif    
            end select
          endif
        
        enddo
      endif

c
c.... end test case 16
c
      if ( casenumber .eq. 17 ) then


        rotBandForcetmp(:,:) = zero
        
        do i = 1, numnp
          counter = 0
c.. counter is used for averaging forces on vertex shared by multiple
c.. rotating bands/faces

          if (m2gClsfcn(i,1) .ne. 3) then
            do j = 1,numRotBands
              do k = 1,numRotBandFaceTags
                call core_is_in_closure(m2gClsfcn(i,1), m2gClsfcn(i,2),
     &                            2,              rotBandTag(j,k),
     &                            answer)          
                if (answer. ne. 0) then
                  Forcetmp(i,:) = Forcetmp(i,:) + rotBandForce(j,:)
                  counter=counter+1                                   
                endif
              enddo
            enddo
          endif
c... Average Forces on shared vertices
          if (counter .gt. 0)
            Forcetmp(i,:) = Forcetmp(i,:)/counter
          endif
         enddo        

c... end loop over vertices

      end if
c
c.... end test case 17
c

      return
      end
c

