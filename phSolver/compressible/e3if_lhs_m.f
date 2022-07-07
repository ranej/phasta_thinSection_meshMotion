      module e3if_lhs_m
c
c----------------------------------------
c    detail subroutines to computer the
c    LHS for interface elements
c----------------------------------------
        use e3if_param_m
        use e3if_func_m
c
        implicit none
c
      contains
c
c      subroutine set_lhs_matrices
c
c        integer :: isd,jsd,n,p,q,r
c
c        do q = 1,nflow
c          do p = 1,nflow
c
c            A0Na_0(:,p,q) = zero
c            A0Na_1(:,p,q) = zero
c
c            do n = 1,nshl0
c              A0Na_0(:,p,q) = A0Na_0(:,p,q) + A0_0(:,p,q)*shp0(:,n)
c            enddo
c
c            do n = 1,nshl1
c              A0Na_1(:,p,q) = A0Na_1(:,p,q) + A0_1(:,p,q)*shp1(:,n)
c            enddo
c
c            do isd = 1,nsd
c
c              AiNa0  (:,isd,p,q) = zero
c              KijNaj0(:,isd,p,q) = zero
c
c              do n = 1,nshl0
c                AiNa0  (:,isd,p,q) = AiNa0(:,isd,p,q) + Ai0(:,isd,p,q)*shp0(:,n)
c                do jsd = 1,nsd
c                  KijNaj0(:,isd,p,q) = KijNaj0(:,isd,p,q) + Kij0(:,isd,jsd,p,q)*shg0(:,n,jsd)
c                enddo
c              enddo
c
c              AiNa1  (:,isd,p,q) = zero
c              KijNaj1(:,isd,p,q) = zero
c
c              do n = 1,nshl1
c                AiNa1  (:,isd,p,q) = AiNa1(:,isd,p,q) + Ai1(:,isd,p,q)*shp1(:,n)
c                do jsd = 1,nsd
c                  KijNaj1(:,isd,p,q) = KijNaj1(:,isd,p,q) + Kij1(:,isd,jsd,p,q)*shg1(:,n,jsd)
c                enddo
c              enddo
c
c            enddo
c          enddo
c        enddo
c
c        do q = 1,nflow
c          do p = 1,nflow
c            do isd = 1,nsd
c
c              KijNajC0(:,isd,p,q) = zero
c              KijNajC1(:,isd,p,q) = zero
c
c              do r = 1,nflow
c                KijNajC0(:,isd,p,q) = KijNajC0(:,isd,p,q) + KijNaj0(:,isd,p,r)*cmtrx(:,r,q)
c                KijNajC1(:,isd,p,q) = KijNajC1(:,isd,p,q) + KijNaj1(:,isd,p,r)*cmtrx(:,r,q)
c              enddo
c
c            enddo
c          enddo
c        enddo
c
c        return
c
c      end subroutine set_lhs_matrices
c
      subroutine calc_egmass_ (
     &                        egmass
     &,                       shp0, shp1, shg0, shg1
     &,                       Ai1
     &,                       Kij0, Kij1
     &,                       ni0, ni1, WdetJ0
     &,                       nshl0, nshl1
     &                        )
        real*8, dimension(:,:,:),   intent(inout) :: egmass
        real*8, dimension(:,:),     intent(in)    :: shp0, shp1
        real*8, dimension(:,:,:),   intent(in)    :: shg0, shg1
        real*8, dimension(:,:,:,:), intent(in)    :: Ai1
        real*8, dimension(:,:,:,:,:), intent(in) :: Kij0, Kij1
        real*8, dimension(:,:),     intent(in)    :: ni0, ni1
        real*8, dimension(:),       intent(in)    :: WdetJ0
        integer, intent(in)  :: nshl0,nshl1
c
        integer :: i,j,i0,j0,il,jl,iflow,jflow,kflow,isd,jsd
        real*8 :: this_mu(npro,nflow,nflow)
        real*8, dimension(npro) :: Ai1Na1ni0,Kij1Naj1ni0,Kij0Nbj0CNa1ni1,CNb0ni0CNa1ni1
c
c---------> ATTENTION <--------
c  This still needs to be confirmed...
c
      this_mu = zero
      this_mu(:,2,2) = mu(:,2) ! mu
      this_mu(:,3,3) = mu(:,3) ! "
      this_mu(:,4,4) = mu(:,4) ! "
      this_mu(:,5,5) = mu(:,5) ! kappa
c
c---------> END ATTANETION <-----
c
        do i = 1,nshl0
          i0 = nflow*(i-1)
c
          do j = 1,nshl0
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                Ai1Na1ni0       = zero
                Kij1Naj1ni0     = zero
                Kij0Nbj0CNa1ni1 = zero
                CNb0ni0CNa1ni1  = zero
c
                do isd = 1,nsd
                  Ai1Na1ni0 = Ai1Na1ni0 + Ai1(:,isd,iflow,jflow)*shp1(:,j)*ni0(:,isd)
                  do jsd = 1,nsd
                    Kij1Naj1ni0 = Kij1Naj1ni0 + Kij1(:,isd,jsd,iflow,jflow)*shg1(:,j,jsd)*ni0(:,isd)
                    do kflow = 1,nflow
                      Kij0Nbj0CNa1ni1 = Kij0Nbj0CNa1ni1 + Kij0(:,isd,jsd,iflow,kflow) 
     &                                                  * shg0(:,i,jsd)
     &                                                  * cmtrx(:,kflow,jflow)
     &                                                  * shp1(:,j)
     &                                                  * ni1(:,isd)
                    enddo
                  enddo
!E! NOTE: Have to check this (penalty term...)
                  do kflow = 1,nflow
                      CNb0ni0CNa1ni1 = CNb0ni0CNa1ni1 + cmtrx(:,iflow,kflow) * shp0(:,i) * ni0(:,isd) 
     &                                                * cmtrx(:,kflow,jflow) * shp1(:,j) * ni1(:,isd)
                  enddo
                enddo
c
                egmass(:,il,jl) = egmass(:,il,jl)
     &          - ( 
     &              pt50 * shp0(:,i) * (Ai1Na1ni0 - Kij1Naj1ni0)
     &          +   pt50 * s         * Kij0Nbj0CNa1ni1
     &          +   e * this_mu(:,iflow,jflow) /length_h * CNb0ni0CNa1ni1
     &            ) * WdetJ0
c
              enddo
            enddo
c
          enddo
        enddo
c
      end subroutine calc_egmass_
c
c
      subroutine calc_egmass_fix_sign (
     &                        egmass
     &,                       shp0, shp1, shg0, shg1
     &,                       Ai1
     &,                       Kij0, Kij1
     &,                       ni0, ni1, WdetJ0
     &,                       nshl0, nshl1
     &,                       code )
        real*8, dimension(:,:,:),   intent(inout) :: egmass
        real*8, dimension(:,:),     intent(in)    :: shp0, shp1
        real*8, dimension(:,:,:),   intent(in)    :: shg0, shg1
        real*8, dimension(:,:,:,:), intent(in)    :: Ai1
        real*8, dimension(:,:,:,:,:), intent(in) :: Kij0, Kij1
        real*8, dimension(:,:),     intent(in)    :: ni0, ni1
        real*8, dimension(:),       intent(in)    :: WdetJ0
        integer, intent(in)  :: nshl0,nshl1
c
        integer :: i,j,i0,j0,il,jl,iflow,jflow,kflow,isd,jsd
        real*8 :: this_mu(npro,nflow,nflow)
        real*8, dimension(npro) :: Ai1Na1ni0,Kij1Naj1ni0,Kij0Nbj0CNa1ni1,CNb0ni0CNa1ni1
        real*8 :: factor
        character*4 code
c------------------------------------------------------------------------------
c  calculation the linearization of RHS w.r.t primitive variables
c  control variable:
c     code  = 'same' deal with linearization for same phase, 
c                    (\partial G_b^alpha/ \partial Y_a^alpha)
c           = 'diff' deal with linearization for different phase,
c                    (\partial G_b^alpha/ \partial Y_a^beta)
c------------------------------------------------------------------------------
c---------> ATTENTION <--------
c  This still needs to be confirmed...
c
      this_mu = zero
      this_mu(:,2,2) = mu(:,2) ! mu
      this_mu(:,3,3) = mu(:,3) ! "
      this_mu(:,4,4) = mu(:,4) ! "
      this_mu(:,5,5) = mu(:,5) ! kappa
c
c---------> END ATTANETION <-----
c---determine the sign of term associated with advetive flux
      if ( code == 'same') then
        factor = -one
      elseif (code == 'diff') then
        factor = one
      else
        call error ('cal_egmass', 'wrong code', 0)
      endif
c
        do i = 1,nshl0
          i0 = nflow*(i-1)
c
          do j = 1,nshl1
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                Ai1Na1ni0       = zero
                Kij1Naj1ni0     = zero
                Kij0Nbj0CNa1ni1 = zero
                CNb0ni0CNa1ni1  = zero
c
                do isd = 1,nsd
                  Ai1Na1ni0 = Ai1Na1ni0 + Ai1(:,isd,iflow,jflow)*shp1(:,j)*ni0(:,isd)
                  do jsd = 1,nsd
                    Kij1Naj1ni0 = Kij1Naj1ni0 + Kij1(:,isd,jsd,iflow,jflow)*shg1(:,j,jsd)*ni0(:,isd)
                    do kflow = 1,nflow
                      Kij0Nbj0CNa1ni1 = Kij0Nbj0CNa1ni1 + Kij0(:,isd,jsd,iflow,kflow) 
     &                                                  * shg0(:,i,jsd)
     &                                                  * cmtrx(:,kflow,jflow)
     &                                                  * shp1(:,j)
     &                                                  * ni0(:,isd)
                    enddo
                  enddo
                enddo
!E! NOTE: Have to check this (penalty term...)
                do kflow = 1,nflow
                    CNb0ni0CNa1ni1 = CNb0ni0CNa1ni1 + cmtrx(:,iflow,kflow) * shp0(:,i)
     &                                              * mu(:,kflow)
     &                                              * cmtrx(:,kflow,jflow) * shp1(:,j)
                enddo
c
                egmass(:,il,jl) = egmass(:,il,jl)
     &          + ( 
     &              pt50 * shp0(:,i) * ( factor*Ai1Na1ni0 - Kij1Naj1ni0)
     &          -   factor * pt50 * s * Kij0Nbj0CNa1ni1
     &          -   factor * e /length_h * CNb0ni0CNa1ni1
     &            ) * WdetJ0
c
              enddo
            enddo
c
          enddo
        enddo
c
      end subroutine calc_egmass_fix_sign                    
c
c      subroutine calc_egmass( egmass00, egmass01,
c     &                         A0_0, A0_1, Ai0, Ai1,
c     &                         Kij0, Kij1,
c     &                         AiNa0_, AiNa1_, KijNaj0_, KijNaj1_,
c     &                         KijNajC0, KijNajC1,
c     &                         shp0, n0, n1, WdetJ0,
c     &                         prop,
c     &                         nshl0, nshl1)
c        
c        real*8, dimension(:,:,:), intent(inout) :: egmass00, egmass01
c        real*8, dimension(:,:,:,:), intent(in)  :: AiNa0_, AiNa1_, KijNaj0_, KijNaj1_, KijNajC0, KijNajC1
c        real*8, dimension(:,:,:), intent(in):: A0_0,A0_1
c        real*8, dimension(:,:,:,:), intent(in) :: Ai0, Ai1
c        real*8, dimension(:,:,:,:,:), pointer, intent(in) :: Kij0, Kij1
c        real*8, dimension(:,:), intent(in) :: shp0,n0, n1
c        real*8, dimension(:), intent(in) :: WdetJ0
c        type(prop_t), dimension(:), pointer, intent(in) :: prop
c        integer, intent(in) :: nshl0,nshl1
c
c        integer :: i,j,p,q,inode0,inode1,isd,jsd
c        integer :: i0,j0,il,jl,iflow,jflow
c        real*8 :: this_mu(npro,nflow,nflow)
c        real*8 :: AiNa0(npro), KijNaj0(npro), KijNbjCn0(npro)
c        real*8 :: AiNa1(npro), KijNaj1(npro), KijNbjCn1(npro)
c
c---------> ATTENTION <--------
c  This still needs to be confirmed...
c
c      this_mu = zero
c      this_mu(:,2,2) = mu(:,2) ! mu
c      this_mu(:,3,3) = mu(:,3) ! "
c      this_mu(:,4,4) = mu(:,4) ! "
c      this_mu(:,5,5) = mu(:,5) ! kappa
c
c---------> END ATTANETION <-----
c
c        do i = 1,nshl0
c          i0 = nflow*(i-1)
c
c          do j = 1,nshl0
c            j0 = nflow*(j-1)
c
c            do iflow = 1,nflow
c              do jflow = 1,nflow
c
c                il = i0 + iflow
c                jl = j0 + jflow
c
c                AiNa0 = zero
c                KijNaj0 = zero
c                KijNbjCn0 = zero
c
c                do isd = 1,nsd
c                  AiNa0 = AiNa0 + pt50 * Ai0(:,isd,iflow,jflow) * n0(:,isd) * shp0(:,j)
c                  KijNbjCn0 = KijNbjCn0 + pt50 * s * KijNajC0(:,isd,iflow,jflow) * n0(:,isd)
c                  do jsd = 1,nsd
c                    KijNaj0 = KijNaj0 + pt50 * Kij0(:,isd,jsd,iflow,jflow) * shg0(:,j,jsd) * n0(:,isd)
c                  enddo
c                enddo
c
c                egmass00(:,il,jl) = egmass00(:,il,jl)
c     &          - ( 
c     &            + AiNa0 
c     &            - KijNaj0
c     &            + KijNbjCn0 * shp0(:,j)
cCc     &            + pt50 * A0_0(:,iflow,jflow)*alpha_LF*shp0(:,j)
cc     &            + pt50 * (A0_0(:,iflow,jflow)+A0_1(:,iflow,jflow))*alpha_LF*shp0(:,j)
c     &            + e*this_mu(:,iflow,jflow)/length_h * ctc(:,iflow,jflow) * shp0(:,j)
c    &            ) * shp0(:,i) * WdetJ0
c             enddo
c            enddo
c
c          enddo
c
c          do j = 1,nshl1
c            j0 = nflow*(j-1)
c
c            do iflow = 1,nflow
c              do jflow = 1,nflow
c
c                il = i0 + iflow
c                jl = j0 + jflow
c
c                AiNa1 = zero
c                KijNaj1 = zero
c                KijNbjCn0 = zero
c
c                do isd = 1,nsd
c                  AiNa1 = AiNa1 + pt50 * Ai1(:,isd,iflow,jflow) * n0(:,isd) * shp1(:,j)
c                  KijNbjCn0 = KijNbjCn0 + pt50 * s * KijNajC0(:,isd,iflow,jflow) * n0(:,isd)
c                  do jsd = 1,nsd
c                    KijNaj1 = KijNaj1 + pt50 * Kij1(:,isd,jsd,iflow,jflow) * shg1(:,j,jsd) * n0(:,isd)
c                  enddo
c                enddo
c
c                egmass01(:,il,jl) = egmass01(:,il,jl)
c     &          - ( 
c     &            + AiNa1
c     &            - KijNaj1
c     &            - KijNbjCn0 * shp1(:,j)
CCc     &            - pt50 * A0_1(:,iflow,jflow)*alpha_LF*shp1(:,j)
cc     &            - pt50 * (A0_0(:,iflow,jflow)+A0_1(:,iflow,jflow))*alpha_LF*shp1(:,j)
c     &            - e*this_mu(:,iflow,jflow)/length_h * ctc(:,iflow,jflow) * shp1(:,j)
c     &            ) * shp0(:,i) * WdetJ0
c              enddo
c            enddo
c
c          enddo
c
c        enddo
c
c      end subroutine calc_egmass
c
      end module e3if_lhs_m
