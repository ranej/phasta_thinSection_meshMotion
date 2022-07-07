      subroutine e3ivar_solid(g1yi_,g2yi_,g3yi_,almBi_,alfBi_,gamBi_)
        use e3_solid_func_m
        implicit none
        real*8, dimension(npro,nflow), target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, intent(in) :: almBi_, alfBi_, gamBi_
c
        almBi = almBi_
        alfBi = alfBi_
        gamBi = gamBi_  
c
        dudx => g1yi_(:,2:4)
        dudy => g2yi_(:,2:4)
        dudz => g3yi_(:,2:4)
c
        call calc_solid
c
      end subroutine e3ivar_solid
c
      subroutine e3bvar_solid(g1yi_,g2yi_,g3yi_,almBi_,alfBi_,gamBi_)
        use e3_solid_func_m
        implicit none
        real*8, dimension(npro,nflow), target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, intent(in) :: almBi_, alfBi_, gamBi_
c
        almBi = almBi_
        alfBi = alfBi_
        gamBi = gamBi_  
c
        dudx => g1yi_(:,2:4)
        dudy => g2yi_(:,2:4)
        dudz => g3yi_(:,2:4)
c
        call calc_solid_bdy
c
      end subroutine e3bvar_solid
