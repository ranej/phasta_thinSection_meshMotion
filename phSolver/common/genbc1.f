        subroutine genBC1 (BCtmp,  iBC,  BC)
c
c----------------------------------------------------------------------
c
c This subroutine adjusts the boundary conditions to accommodate for 
c the velocity constraints in the non-axes directions. It copies the
c reduced constraint parameters in BC.
c
c input:
c  BCtmp (nshg,6+5*I3nsd) : input BC parameters (density, temperature,
c                             pressure, (nsd-1)(nsd+1) velocity params,
c                             upto 4 scalar params)
c  iBC   (nshg)           : boundary condition code
c
c output:
c  BC    (nshg,ndofBC)    : the constraint eq's parameters
c
c
c Farzin Shakib, Winter 1986.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      include "common.h"
        
c
        dimension BCtmp(nshg,ndofBC2),    iBC(nshg),
     &            BC(nshg,ndofBC),        tmpbc(4)
c
        dimension tmp(nshg)
c
        integer elas3, elas4, elas5, elas6
c
c.... scalars
c
        do isclr=1,nsclr
           where (btest(iBC,5+isclr)) BC(:,6+isclr) = BCtmp(:,12+isclr)
        enddo
c
c.... set up the thermodynamic properties
c
        where (btest(iBC,0)) BC(:,1) = BCtmp(:,1) ! density
        where (btest(iBC,1)) BC(:,2) = BCtmp(:,2) ! temperature
        where (btest(iBC,2)) BC(:,1) = BCtmp(:,3) ! pressure
c
c.... if the velocity in the x1-direction is specified
c
        where (ibits(iBC,3,3) .eq. 1)
          tmp     = BCtmp(:,4)**2 + BCtmp(:,5)**2 + BCtmp(:,6)**2
          BC(:,3) = tmp * BCtmp(:,7) / BCtmp(:,4)
          BC(:,4) =       BCtmp(:,5) / BCtmp(:,4)
          BC(:,5) =       BCtmp(:,6) / BCtmp(:,4)
        endwhere
c
c.... if the velocity in the x2-direction is specified
c
        where (ibits(iBC,3,3) .eq. 2)
          tmp     = BCtmp(:,4)**2 + BCtmp(:,5)**2 + BCtmp(:,6)**2
          BC(:,3) = tmp * BCtmp(:,7) / BCtmp(:,5)
          BC(:,4) =       BCtmp(:,4) / BCtmp(:,5)
          BC(:,5) =       BCtmp(:,6) / BCtmp(:,5)
        endwhere
c
c.... if the two velocities are specified (x1 & x2-direction)
c
c
c  Protect against user flipping the order of x1 and x2 in 
c  the vector 1 and vector 2.  Without this it will blow up.
c
        do i=1,nshg
          if(ibits(iBC(i),3,3) .eq. 3 .and. 
     &       (BCtmp(i,4).eq.0 .or. BCtmp(i,9).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,4:7)
              BCtmp(i,4:7)=BCtmp(i,8:11)
              BCtmp(i,8:11)=tmpbc(1:4)
          endif
        enddo
        where (ibits(iBC,3,3) .eq. 3)
          tmp         = sqrt (BCtmp(:, 4)**2 + BCtmp(:, 5)**2
     &                                       + BCtmp(:, 6)**2)
          BCtmp(:, 4) = BCtmp(:, 4) / tmp
          BCtmp(:, 5) = BCtmp(:, 5) / tmp
          BCtmp(:, 6) = BCtmp(:, 6) / tmp
          BCtmp(:, 7) = BCtmp(:, 7) * tmp
c
          tmp         = sqrt (BCtmp(:, 8)**2 + BCtmp(:, 9)**2
     &                                       + BCtmp(:,10)**2)
          BCtmp(:, 8) = BCtmp(:, 8) / tmp
          BCtmp(:, 9) = BCtmp(:, 9) / tmp
          BCtmp(:,10) = BCtmp(:,10) / tmp
          BCtmp(:,11) = BCtmp(:,11) * tmp
c
          BCtmp(:, 4) = BCtmp(:, 9) * BCtmp(:, 4)
     &                - BCtmp(:, 5) * BCtmp(:, 8)
          BCtmp(:, 6) = BCtmp(:, 9) * BCtmp(:, 6)
     &                - BCtmp(:, 5) * BCtmp(:,10)
          BCtmp(:, 7) = BCtmp(:, 9) * BCtmp(:, 7)
     &                - BCtmp(:, 5) * BCtmp(:,11)
          BC(:,3)     = BCtmp(:, 7) / BCtmp(:, 4)
          BC(:,4)     = BCtmp(:, 6) / BCtmp(:, 4)
c
          BCtmp(:, 9) = BCtmp(:, 4) * BCtmp(:, 9) 
          BCtmp(:,10) = BCtmp(:, 4) * BCtmp(:,10)
     &                - BCtmp(:, 8) * BCtmp(:, 6)
          BCtmp(:,11) = BCtmp(:, 4) * BCtmp(:,11)
     &                - BCtmp(:, 8) * BCtmp(:, 7)
          BC(:,5)     = BCtmp(:,11) / BCtmp(:, 9)
          BC(:,6)     = BCtmp(:,10) / BCtmp(:, 9)
        endwhere
c
c.... if the velocity in the x3-direction is specified
c
        if (nsd .eq. 3) then
        where (ibits(iBC,3,3) .eq. 4)
          tmp     = BCtmp(:,4)**2 + BCtmp(:,5)**2 + BCtmp(:,6)**2
          BC(:,3) = tmp * BCtmp(:,7) / BCtmp(:,6)
          BC(:,4) =       BCtmp(:,4) / BCtmp(:,6)
          BC(:,5) =       BCtmp(:,5) / BCtmp(:,6)
        endwhere
        endif
c
c.... if two velocities are specified (x1 & x3-direction)
c
        if (nsd .eq. 3) then
c
c  Protect against user flipping the order of x1 and x3 in 
c  the vector 1 and vector 2.  Without this it will blow up.
c
        do i=1,nshg
          if(ibits(iBC(i),3,3) .eq. 5 .and.
     &       (BCtmp(i,4).eq.0 .or. BCtmp(i,10).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,4:7)
              BCtmp(i,4:7)=BCtmp(i,8:11)
              BCtmp(i,8:11)=tmpbc(1:4)
           endif
        enddo
        where (ibits(iBC,3,3) .eq. 5)
          tmp         = sqrt (BCtmp(:, 4)**2 + BCtmp(:, 5)**2
     &                                       + BCtmp(:, 6)**2)
          BCtmp(:, 4) = BCtmp(:, 4) / tmp
          BCtmp(:, 5) = BCtmp(:, 5) / tmp
          BCtmp(:, 6) = BCtmp(:, 6) / tmp
          BCtmp(:, 7) = BCtmp(:, 7) * tmp
c
          tmp         = sqrt (BCtmp(:, 8)**2 + BCtmp(:, 9)**2
     &                                       + BCtmp(:,10)**2)
          BCtmp(:, 8) = BCtmp(:, 8) / tmp
          BCtmp(:, 9) = BCtmp(:, 9) / tmp
          BCtmp(:,10) = BCtmp(:,10) / tmp
          BCtmp(:,11) = BCtmp(:,11) * tmp
c
          BCtmp(:, 4) = BCtmp(:,10) * BCtmp(:, 4)
     &                - BCtmp(:, 6) * BCtmp(:, 8)
          BCtmp(:, 5) = BCtmp(:,10) * BCtmp(:, 5)
     &                - BCtmp(:, 6) * BCtmp(:, 9)
          BCtmp(:, 7) = BCtmp(:,10) * BCtmp(:, 7)
     &                - BCtmp(:, 6) * BCtmp(:,11)
          BC(:,3)     = BCtmp(:, 7) / BCtmp(:, 4)
          BC(:,4)     = BCtmp(:, 5) / BCtmp(:, 4)
c
          BCtmp(:, 9) = BCtmp(:, 4) * BCtmp(:, 9)
     &                - BCtmp(:, 8) * BCtmp(:, 5)
          BCtmp(:,10) = BCtmp(:, 4) * BCtmp(:,10)
          BCtmp(:,11) = BCtmp(:, 4) * BCtmp(:,11)
     &                - BCtmp(:, 8) * BCtmp(:, 7)
          BC(:,5)     = BCtmp(:,11) / BCtmp(:,10)
          BC(:,6)     = BCtmp(:, 9) / BCtmp(:,10)
        endwhere
        endif
c
c.... if two velocities are specified (x2 & x3-direction)
c
        if (nsd .eq. 3) then
c
c  Protect against user flipping the order of x2 and x3 in 
c  the vector 1 and vector 2.  Without this it will blow up.
c
        do i=1,nshg
          if(ibits(iBC(i),3,3) .eq. 6 .and. (
     &       BCtmp(i,5).eq.0 .or. BCtmp(i,10).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,4:7)
              BCtmp(i,4:7)=BCtmp(i,8:11)
              BCtmp(i,8:11)=tmpbc(1:4)
           endif
        enddo
        where (ibits(iBC,3,3) .eq. 6)
          tmp         = sqrt (BCtmp(:, 4)**2 + BCtmp(:, 5)**2
     &                                       + BCtmp(:, 6)**2)
          BCtmp(:, 4) = BCtmp(:, 4) / tmp
          BCtmp(:, 5) = BCtmp(:, 5) / tmp
          BCtmp(:, 6) = BCtmp(:, 6) / tmp
          BCtmp(:, 7) = BCtmp(:, 7) * tmp
c
          tmp         = sqrt (BCtmp(:, 8)**2 + BCtmp(:, 9)**2
     &                                       + BCtmp(:,10)**2)
          BCtmp(:, 8) = BCtmp(:, 8) / tmp
          BCtmp(:, 9) = BCtmp(:, 9) / tmp
          BCtmp(:,10) = BCtmp(:,10) / tmp
          BCtmp(:,11) = BCtmp(:,11) * tmp
c
          BCtmp(:, 4) = BCtmp(:,10) * BCtmp(:, 4)
     &                - BCtmp(:, 6) * BCtmp(:, 8)
          BCtmp(:, 5) = BCtmp(:,10) * BCtmp(:, 5)
     &                - BCtmp(:, 6) * BCtmp(:, 9)
          BCtmp(:, 7) = BCtmp(:,10) * BCtmp(:, 7)
     &                - BCtmp(:, 6) * BCtmp(:,11)
          BC(:,3)     = BCtmp(:, 7) / BCtmp(:, 5)
          BC(:,4)     = BCtmp(:, 4) / BCtmp(:, 5)
c
          BCtmp(:, 8) = BCtmp(:, 5) * BCtmp(:, 8)
     &                - BCtmp(:, 9) * BCtmp(:, 4) 
          BCtmp(:,10) = BCtmp(:, 5) * BCtmp(:,10)
          BCtmp(:,11) = BCtmp(:, 5) * BCtmp(:,11)
     &                - BCtmp(:, 9) * BCtmp(:, 7)
          BC(:,5)     = BCtmp(:,11) / BCtmp(:,10)
          BC(:,6)     = BCtmp(:, 8) / BCtmp(:,10)
        endwhere
        endif
c
c.... if all velocities are specified
c

        if (nsd .eq. 3) then
        where (ibits(iBC,3,3) .eq. 7)
          BC(:,3) = BCtmp(:,7) * BCtmp(:,4)
          BC(:,4) = BCtmp(:,7) * BCtmp(:,5)
          BC(:,5) = BCtmp(:,7) * BCtmp(:,6)
        endwhere
        endif

c      write(*,*) ' In genbc1: ibc, bc'
c      do i = 1,nshg
c        if (ibits(iBC(i),3,3) .eq. 7) then
c          write(*,'(a,x,i4,2x,12f6.2)')'comp3:', iBC(i),BCtmp(i,:)
c          write(*,'(a,x,i4,2x,6f10.2)')'comp3:', iBC(i),BC(i,:)
c        endif
c      enddo
c
c.... end
c

c.... if we have the mesh elastic solver
      if (iALE .eq. 2) then 
c
c----------------Re-organize Mesh Elastic Velocities-----------------
c.... convert mesh velocity to mesh displacement
c
        elas3 = ndof + 2
        elas4 = ndof + 3
        elas5 = ndof + 4
        elas6 = ndof + 5
c
c.... if the velocity in the x1-direction is specified
c
        where (ibits(iBC,14,3) .eq. 1)
          tmp     = BCtmp(:,17)**2 + BCtmp(:,18)**2 + BCtmp(:,19)**2
          BC(:,elas3) = tmp * BCtmp(:,20) / BCtmp(:,17) * Delt(1)
          BC(:,elas4) =       BCtmp(:,18) / BCtmp(:,17)
          BC(:,elas5) =       BCtmp(:,19) / BCtmp(:,17)
        endwhere
c
c.... if the velocity in the x2-direction is specified
c
        where (ibits(iBC,14,3) .eq. 2)
          tmp     = BCtmp(:,17)**2 + BCtmp(:,18)**2 + BCtmp(:,19)**2
          BC(:,elas3) = tmp * BCtmp(:,20) / BCtmp(:,18) * Delt(1)
          BC(:,elas4) =       BCtmp(:,17) / BCtmp(:,18)
          BC(:,elas5) =       BCtmp(:,19) / BCtmp(:,18)
        endwhere
c
c.... if the two velocities are specified (x1 & x2-direction)
c
c
c  Protect against user flipping the order of x1 and x2 in 
c  the vector 1 and vector 2.  Without this it will blow up.
c
        do i=1,nshg
          if(ibits(iBC(i),14,3) .eq. 3 .and. 
     &       (BCtmp(i,17).eq.0 .or. BCtmp(i,22).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,17:20)
              BCtmp(i,17:20)=BCtmp(i,21:24)
              BCtmp(i,21:24)=tmpbc(1:4)
          endif
        enddo
        where (ibits(iBC,14,3) .eq. 3)
          tmp         = sqrt (BCtmp(:, 17)**2 + BCtmp(:, 18)**2
     &                                       + BCtmp(:, 19)**2)
          BCtmp(:, 17) = BCtmp(:, 17) / tmp
          BCtmp(:, 18) = BCtmp(:, 18) / tmp
          BCtmp(:, 19) = BCtmp(:, 19) / tmp
          BCtmp(:, 20) = BCtmp(:, 20) * tmp
c
          tmp         = sqrt (BCtmp(:, 21)**2 + BCtmp(:, 22)**2
     &                                       + BCtmp(:,23)**2)
          BCtmp(:, 21) = BCtmp(:, 21) / tmp
          BCtmp(:, 22) = BCtmp(:, 22) / tmp
          BCtmp(:, 23) = BCtmp(:, 23) / tmp
          BCtmp(:, 24) = BCtmp(:, 24) * tmp
c
          BCtmp(:, 17) = BCtmp(:, 22) * BCtmp(:, 17)
     &                 - BCtmp(:, 18) * BCtmp(:, 21)
          BCtmp(:, 19) = BCtmp(:, 22) * BCtmp(:, 19)
     &                 - BCtmp(:, 18) * BCtmp(:, 23)
          BCtmp(:, 20) = BCtmp(:, 22) * BCtmp(:, 20)
     &                 - BCtmp(:, 18) * BCtmp(:, 23)
          BC(:,elas3)  = BCtmp(:, 20) / BCtmp(:, 17) * Delt(1)
          BC(:,elas4)  = BCtmp(:, 19) / BCtmp(:, 17)
c
          BCtmp(:, 22) = BCtmp(:, 17) * BCtmp(:, 22) 
          BCtmp(:, 23) = BCtmp(:, 17) * BCtmp(:, 23)
     &                 - BCtmp(:, 21) * BCtmp(:, 19)
          BCtmp(:, 24) = BCtmp(:, 17) * BCtmp(:, 24)
     &                 - BCtmp(:, 21) * BCtmp(:, 20)
          BC(:,elas5)  = BCtmp(:, 24) / BCtmp(:, 22) * Delt(1)
          BC(:,elas6)  = BCtmp(:, 23) / BCtmp(:, 22)
        endwhere
c
c.... if the velocity in the x3-direction is specified
c
        if (nsd .eq. 3) then
        where (ibits(iBC,3,3) .eq. 4)
          tmp     = BCtmp(:,17)**2 + BCtmp(:,18)**2 + BCtmp(:,19)**2
          BC(:,elas3) = tmp * BCtmp(:,20) / BCtmp(:,19) * Delt(1)
          BC(:,elas4) =       BCtmp(:,17) / BCtmp(:,19)
          BC(:,elas5) =       BCtmp(:,18) / BCtmp(:,19)
        endwhere
        endif
c
c.... if two velocities are specified (x1 & x3-direction)
c
        if (nsd .eq. 3) then
c
c  Protect against user flipping the order of x1 and x3 in 
c  the vector 1 and vector 2.  Without this it will blow up.
c
        do i=1,nshg
          if(ibits(iBC(i),14,3) .eq. 5 .and.
     &       (BCtmp(i,17).eq.0 .or. BCtmp(i,23).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,17:20)
              BCtmp(i,17:20)=BCtmp(i,21:24)
              BCtmp(i,21:24)=tmpbc(1:4)
           endif
        enddo
        where (ibits(iBC,14,3) .eq. 5)
          tmp         = sqrt (BCtmp(:, 17)**2 + BCtmp(:, 18)**2
     &                                       + BCtmp(:, 19)**2)
          BCtmp(:, 17) = BCtmp(:, 17) / tmp
          BCtmp(:, 18) = BCtmp(:, 18) / tmp
          BCtmp(:, 19) = BCtmp(:, 19) / tmp
          BCtmp(:, 20) = BCtmp(:, 20) * tmp
c
          tmp         = sqrt (BCtmp(:, 21)**2 + BCtmp(:, 22)**2
     &                                       + BCtmp(:, 23)**2)
          BCtmp(:, 21) = BCtmp(:, 21) / tmp
          BCtmp(:, 22) = BCtmp(:, 22) / tmp
          BCtmp(:, 23) = BCtmp(:, 23) / tmp
          BCtmp(:, 24) = BCtmp(:, 24) * tmp
c
          BCtmp(:, 17) = BCtmp(:, 23) * BCtmp(:, 17)
     &                 - BCtmp(:, 19) * BCtmp(:, 21)
          BCtmp(:, 18) = BCtmp(:, 23) * BCtmp(:, 18)
     &                 - BCtmp(:, 19) * BCtmp(:, 22)
          BCtmp(:, 20) = BCtmp(:, 23) * BCtmp(:, 20)
     &                 - BCtmp(:, 19) * BCtmp(:, 24)
          BC(:,elas3)  = BCtmp(:, 20) / BCtmp(:, 17) * Delt(1)
          BC(:,elas4)  = BCtmp(:, 18) / BCtmp(:, 17)
c
          BCtmp(:, 22) = BCtmp(:, 17) * BCtmp(:, 22)
     &                 - BCtmp(:, 21) * BCtmp(:, 18)
          BCtmp(:, 23) = BCtmp(:, 17) * BCtmp(:, 23)
          BCtmp(:, 24) = BCtmp(:, 17) * BCtmp(:, 24)
     &                 - BCtmp(:, 21) * BCtmp(:, 20)
          BC(:,elas5)  = BCtmp(:, 24) / BCtmp(:, 23) * Delt(1)
          BC(:,elas6)  = BCtmp(:, 22) / BCtmp(:, 23)
        endwhere
        endif
c
c.... if two velocities are specified (x2 & x3-direction)
c
        if (nsd .eq. 3) then
c
c  Protect against user flipping the order of x2 and x3 in 
c  the vector 1 and vector 2.  Without this it will blow up.
c
        do i=1,nshg
          if(ibits(iBC(i),14,3) .eq. 6 .and. (
     &       BCtmp(i,18).eq.0 .or. BCtmp(i,23).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,17:20)
              BCtmp(i,17:20)=BCtmp(i,21:24)
              BCtmp(i,21:24)=tmpbc(1:4)
           endif
        enddo
        where (ibits(iBC,14,3) .eq. 6)
          tmp         = sqrt (BCtmp(:, 17)**2 + BCtmp(:, 18)**2
     &                                       + BCtmp(:, 19)**2)
          BCtmp(:, 17) = BCtmp(:, 17) / tmp
          BCtmp(:, 18) = BCtmp(:, 18) / tmp
          BCtmp(:, 19) = BCtmp(:, 19) / tmp
          BCtmp(:, 20) = BCtmp(:, 20) * tmp
c
          tmp         = sqrt (BCtmp(:, 21)**2 + BCtmp(:, 22)**2
     &                                       + BCtmp(:,23)**2)
          BCtmp(:, 21) = BCtmp(:, 21) / tmp
          BCtmp(:, 22) = BCtmp(:, 22) / tmp
          BCtmp(:, 23) = BCtmp(:, 23) / tmp
          BCtmp(:, 24) = BCtmp(:, 24) * tmp
c
          BCtmp(:, 17) = BCtmp(:, 23) * BCtmp(:, 17)
     &                 - BCtmp(:, 19) * BCtmp(:, 21)
          BCtmp(:, 18) = BCtmp(:, 23) * BCtmp(:, 18)
     &                 - BCtmp(:, 19) * BCtmp(:, 22)
          BCtmp(:, 20) = BCtmp(:, 23) * BCtmp(:, 20)
     &                 - BCtmp(:, 19) * BCtmp(:, 23)
          BC(:,elas3)  = BCtmp(:, 20) / BCtmp(:, 18) * Delt(1)
          BC(:,elas4)  = BCtmp(:, 17) / BCtmp(:, 18)
c
          BCtmp(:, 21) = BCtmp(:, 18) * BCtmp(:, 21)
     &                 - BCtmp(:, 22) * BCtmp(:, 17) 
          BCtmp(:, 23) = BCtmp(:, 18) * BCtmp(:, 23)
          BCtmp(:, 24) = BCtmp(:, 18) * BCtmp(:, 24)
     &                 - BCtmp(:, 22) * BCtmp(:, 20)
          BC(:,elas5)  = BCtmp(:, 24) / BCtmp(:, 23) * Delt(1)
          BC(:,elas6)  = BCtmp(:, 21) / BCtmp(:, 23)
        endwhere
        endif
c
c.... if all velocities are specified
c

        if (nsd .eq. 3) then
        where (ibits(iBC,14,3) .eq. 7)
          BC(:,elas3) = BCtmp(:,20) * BCtmp(:,17) * Delt(1)
          BC(:,elas4) = BCtmp(:,20) * BCtmp(:,18) * Delt(1)
          BC(:,elas5) = BCtmp(:,20) * BCtmp(:,19) * Delt(1)
        endwhere
        endif
      endif ! end if iALE = 2
c----------------End Mesh Elastic Velocity---------------
c
c.... end
        return
        end
