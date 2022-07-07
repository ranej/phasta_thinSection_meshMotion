      module e3if_diff_m
c
c----------------------------------------
c    aims to calculate the diffusive coefficient
c    matrix for interface elements
c----------------------------------------
        use e3if_param_m
c
        implicit none
c
      contains
c
        subroutine calc_stiff (prop,var,mater)
c
          type(prop_t), dimension(:), pointer, intent(inout) :: prop
          type(var_t),  dimension(:), pointer, intent(in)    :: var
          integer, intent(in) :: mater
c
          integer :: i,j
          real*8, dimension(npro) :: rmu, rlm, rlm2mu, con
          real*8, dimension(:), pointer :: u1, u2, u3
c
          call getdiff(rmu, rlm, rlm2mu, con, npro, mater)
c
c          rmu    => prop%rmu
c          rlm    => prop%rlm
c          rlm2mu => prop%rlm2mu
c          con    => prop%con
c
          u1     => var%y(2)
          u2     => var%y(3)
          u3     => var%y(4)
c
          do i = 1,nsd*nflow
            do j = 1,nsd*nflow
              prop%stiff(i,j) = zero
            enddo
          enddo
c
c.... K11
c
          prop%stiff(2, 2) = rlm2mu
          prop%stiff(3, 3) = rmu
          prop%stiff(4, 4) = rmu
          prop%stiff(5, 2) = rlm2mu * u1
          prop%stiff(5, 3) = rmu    * u2
          prop%stiff(5, 4) = rmu    * u3
          prop%stiff(5, 5) = con
c
c     
c.... K12
c     
          prop%stiff(2, 8) = rlm
          prop%stiff(3, 7) = rmu
          prop%stiff(5, 7) = rmu    * u2
          prop%stiff(5, 8) = rlm    * u1
c     
c.... K13
c     
          prop%stiff(2,14) = rlm
          prop%stiff(4,12) = rmu
          prop%stiff(5,12) = rmu    * u3
          prop%stiff(5,14) = rlm    * u1
           
c     
c.... K21
c     
          prop%stiff(7, 3) = rmu
          prop%stiff(8, 2) = rlm
          prop%stiff(10, 2) = rlm    * u2
          prop%stiff(10, 3) = rmu    * u1
           
c     
c.... K22
c     
          prop%stiff(7, 7) = rmu
          prop%stiff(8, 8) = rlm2mu
          prop%stiff(9, 9) = rmu
          prop%stiff(10, 7) = rmu    * u1
          prop%stiff(10, 8) = rlm2mu * u2
          prop%stiff(10, 9) = rmu    * u3
          prop%stiff(10,10) = con
c     
c.... K23
c     
          prop%stiff(8,14) = rlm
          prop%stiff(9,13) = rmu
          prop%stiff(10,13) = rmu    * u3
          prop%stiff(10,14) = rlm    * u2
c     
c.... K31
c     
          prop%stiff(12, 4) = rmu
          prop%stiff(14, 2) = rlm
          prop%stiff(15, 2) = rlm    * u3
          prop%stiff(15, 4) = rmu    * u1
c     
c.... K32
c     
          prop%stiff(13, 9) = rmu
          prop%stiff(14, 8) = rlm
          prop%stiff(15, 8) = rlm    * u3
          prop%stiff(15, 9) = rmu    * u2
c     
c.... K33
c     
          prop%stiff(12,12) = rmu
          prop%stiff(13,13) = rmu
          prop%stiff(14,14) = rlm2mu
          prop%stiff(15,12) = rmu    * u1
          prop%stiff(15,13) = rmu    * u2
          prop%stiff(15,14) = rlm2mu * u3
          prop%stiff(15,15) = con
c
        end subroutine calc_stiff
c
      end module e3if_diff_m
