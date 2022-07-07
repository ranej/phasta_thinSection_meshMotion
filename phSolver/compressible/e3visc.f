        subroutine e3visc (g1yi,    g2yi,    g3yi,
     &                     dxidx,
     &                     rmu,     rlm,     rlm2mu,
     &                     u1,      u2,      u3,      
     &                     ri,      rmi,     stiff,
     &                     con,     rlsli,   compK)
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of viscous and heat fluxes
c to both RHS and LHS.
c
c input:
c  g1yi   (npro,nflow)         : grad-y in direction 1
c  g2yi   (npro,nflow)         : grad-y in direction 2
c  g3yi   (npro,nflow)         : grad-y in direction 3
c  u1     (npro)              : x1-velocity component
c  u2     (npro)              : x2-velocity component
c  u3     (npro)              : x3-velocity component
c  rmu    (npro)              : viscosity
c  rlm    (npro)              : Lambda
c  rlm2mu (npro)              : Lambda + 2 Mu
c  con    (npro)              : Conductivity
c  cp     (npro)              : specific heat at constant pressure
c  rlsli  (npro,6)            : Resolved Reynolds stresses
c output:
c  ri     (npro,nflow*(nsd+1)) : partial residual
c  rmi    (npro,nflow*(nsd+1)) : partial modified residual
c  stiff  (npro,nsd*nflow,nsd*nflow) : stiffness matrix
c  compK (npro,10)             : K_ij in (eq.134) A new ... III 
c
c
c Zdenek Johan, Summer 1990. (Modified from e2visc.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c
      use e3_param_m
      use e3_solid_m
      use solid_m
c
      include "common.h"
c
c     passed arrays
c
      dimension g1yi(npro,nflow),           g2yi(npro,nflow),
     &     g3yi(npro,nflow),   
     &     Sclr(npro),                dxidx(npro,nsd,nsd),
     &     u1(npro),                  u2(npro),
     &     u3(npro),
     &     ri(npro,nflow*(nsd+1)),     rmi(npro,nflow*(nsd+1))
c
c
      dimension rmu(npro),                 rlm(npro),
     &     rlm2mu(npro),                   con(npro),
     &     stiff(npro,3*nflow,3*nflow),
     &     rlsli(npro,6),                  compK(npro,10),
     &     f1(npro), f2(npro), f3(npro), f4(npro), 
     &     f5(npro), f6(npro)
c     &,     rk(npro) 
c.....Used for solid calculation
      real*8, dimension(npro) :: d_temp1,d_temp2,d_temp3
      real*8, dimension(npro,6) :: bq_af
c.....
      ttim(23) = ttim(23) - secs(0.0)
c
c... dynamic model is now being accounted for in getdiff.f
c
c
c.... --------------------->  Diffusivity Matrix  <-------------------
c
      if (lhs .eq. 1) then

c......if this is a solid blk, change the K_ij matrix
        if (mat_eos(mater,1).eq.ieos_solid_1)then
c
        d_temp1(:) = 2.0/3.0*d(:,1) - 1.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
        d_temp2(:) = -1.0/3.0*d(:,1) + 2.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
        d_temp3(:) = -1.0/3.0*d(:,1) - 1.0/3.0*d(:,2) + 2.0/3.0*d(:,3)
c.... K11
c
         stiff(:, 2, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     $                    ( (4.0/3.0)*d(:,1) - (5.0/3.0)* (almBi)**2 * d_temp1 )  
         stiff(:, 2, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     $                    (-2.0/3.0)*d(:,6)
         stiff(:, 2, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    (-2.0/3.0)*d(:,5)
         stiff(:, 3, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    (d(:,6)- (5.0/3.0)* (almBi)**2 * d(:,6))
         stiff(:, 3, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,1) )
         stiff(:, 4, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    (d(:,5)- (5.0/3.0)* (almBi)**2 * d(:,5))
         stiff(:, 4, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,1))
         stiff(:, 5, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,1) * u1 + d(:,6) * u2 + d(:,5) * u3 - 
     &                    (5.0/3.0)* (almBi)**2 *( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
         stiff(:, 5, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( (-2.0/3.0)*d(:,6) * u1 + d(:,1) * u2 )
         stiff(:, 5, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( (-2.0/3.0)*d(:,5) * u1 + d(:,1) * u3 )
         stiff(:, 5, 5) = con ! notice the + or -
c     
c.... K12
c     
         stiff(:, 2, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,6) )  
         stiff(:, 2, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,2) - (5.0/3.0)* (almBi)**2 * d_temp1 )
         stiff(:, 2, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,4) )
         stiff(:, 3, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,2))
         stiff(:, 3, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,6) - (5.0/3.0)* (almBi)**2 * d(:,6) )
         stiff(:, 4, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,4) )
         stiff(:, 4, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( 0.0 - (5.0/3.0)* (almBi)**2 * d(:,5) )    
         stiff(:, 4, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,6) )
          stiff(:, 5, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,6) * u1 + d(:,2) * u2 + d(:,4) * u3 )
          stiff(:, 5, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,2) * u1 + d(:,6) * u2 - (5.0/3.0)* (almBi)**2 *
     &                    ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
         stiff(:, 5, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,4) * u1 + d(:,6) * u3 )
c
c.... K13
c     
         stiff(:, 2,12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0 /3.0)*d(:,5) )  
         stiff(:, 2,13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0 /3.0)*d(:,4) )
         stiff(:, 2,14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0 /3.0)*d(:,3) - (5.0/3.0)* (almBi)**2 * d_temp1)
         stiff(:, 3,12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,4))
         stiff(:, 3,13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,5) )
         stiff(:, 3,14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( - (5.0/3.0)* (almBi)**2 * d(:,6) )   
         stiff(:, 4,12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,3) )  
         stiff(:, 4,14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,5)- (5.0/3.0)* (almBi)**2 * d(:,5) )
         stiff(:, 5,12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0 /3.0)*d(:,5) * u1 + d(:,4) * u2 + d(:,3) * u3 )
         stiff(:, 5,13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0 /3.0)*d(:,4) * u1 + d(:,5) * u2  )
         stiff(:, 5,14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( (-2.0 /3.0)*d(:,3) * u1 + d(:,5) * u3 - (5.0/3.0)* (almBi)**2 *
     &                    ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3))
c     
c.... K21
c     
         stiff(:, 7, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    (d(:,6)- (5.0/3.0)* (almBi)**2 * d(:,6))
         stiff(:, 7, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*d(:,1)
         stiff(:, 8, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ((-2.0/3.0)*d(:,1) - (5.0/3.0)* (almBi)**2 * d_temp2 )  
        stiff(:, 8, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     (4.0/3.0)*d(:,6)
         stiff(:, 8, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     (-2.0/3.0)*d(:,5)
         stiff(:, 9, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (- (5.0/3.0)* (almBi)**2 * d(:,4))
         stiff(:, 9, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    d(:,5)
         stiff(:, 9, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    d(:,6)
         stiff(:, 10, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,6) * u1 + (-2.0/3.0)*d(:,1) * u2 - (5.0/3.0)* (almBi)**2 *
     &                    ( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3) )
         stiff(:, 10, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,1) * u1 + (4.0/3.0)*d(:,6) * u2 +d(:,5) * u3)
         stiff(:, 10, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ((-2.0/3.0)*d(:,5) * u2 +d(:,6) * u3)
c     
c.... K22
     
         stiff(:, 7, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,2) )
         stiff(:, 7, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,6) -(5.0/3.0)* (almBi)**2 * d(:,6) )
         stiff(:, 8, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( -(2.0/3.0)*d(:,6) )
         stiff(:, 8, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,2) -(5.0/3.0)* (almBi)**2 * d_temp2 )
         stiff(:, 8, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,4) )   
         stiff(:, 9, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,4) -(5.0/3.0)* (almBi)**2 * d(:,4))
         stiff(:, 9, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,2) )
         stiff(:, 10, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,2) * u1 + (-2.0/3.0)*d(:,6) * u2 )
         stiff(:, 10, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,6) * u1 + (4.0/3.0)*d(:,2) * u2 +d(:,4) * u3 -
     &                    (5.0/3.0)* (almBi)**2 *( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3))
         stiff(:, 10, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ((-2.0/3.0)*d(:,4) * u2 +d(:,2) * u3)
         stiff(:, 10, 10) = con ! notice the + or -
c     
c.... K23
c     
         stiff(:, 7, 12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,4) )
         stiff(:, 7, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,5) )
         stiff(:, 7, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( 0.0 -(5.0/3.0)* (almBi)**2 * d(:,6) )
         stiff(:, 8, 12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( -(2.0/3.0)*d(:,5) )
         stiff(:, 8, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,4) )
         stiff(:, 8, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,3) -(5.0/3.0)* (almBi)**2 * d_temp2 )   
         stiff(:, 9, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,3) )
         stiff(:, 9, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,4) -(5.0/3.0)* (almBi)**2 * d(:,4))
         stiff(:, 10, 12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,4) * u1 + (-2.0/3.0)*d(:,5) * u2 )
         stiff(:, 10, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,5) * u1 + (4.0/3.0)*d(:,4) * u2 +d(:,3) * u3 )
         stiff(:, 10, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ((-2.0/3.0)*d(:,3) * u2 +d(:,4) * u3 -
     &                    (5.0/3.0)* (almBi)**2 *( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3))
c     
c.... K31
c
         stiff(:, 12, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,5) - (5.0/3.0)* (almBi)**2 *d(:,5) )
         stiff(:, 12, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,1) )
         stiff(:, 13, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( 0.0 - (5.0/3.0)* (almBi)**2 *d(:,4))
         stiff(:, 13, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,5) )
         stiff(:, 13, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,6) )
         stiff(:, 14, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (-2.0/3.0)*d(:,1) - (5.0/3.0)* (almBi)**2 *d_temp3 )   
         stiff(:, 14, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( (-2.0/3.0)*d(:,6) )
         stiff(:, 14, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( (4.0/3.0)*d(:,5) )
         stiff(:, 15, 2) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,5) * u1 + (-2.0/3.0)*d(:,1) * u3 -
     &                    (5.0/3.0)* (almBi)**2 *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3) )
         stiff(:, 15, 3) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,5) * u2 +(-2.0/3.0)*d(:,6) * u3 )     
         stiff(:, 15, 4) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,1) * u1 + (4.0/3.0)*d(:,5) * u3 +d(:,6) * u2 )
c
c.... K32
c     
         stiff(:, 12, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,4) )
         stiff(:, 12, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( 0.0 - (5.0/3.0)* (almBi)**2 *d(:,5) )
         stiff(:, 12, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,6) )
         stiff(:, 13, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,4) - (5.0/3.0)* (almBi)**2 *d(:,4) )
         stiff(:, 13, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,2) )
         stiff(:, 14, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( -(2.0/3.0)*d(:,6) )
         stiff(:, 14, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( -(2.0/3.0)*d(:,2) - (5.0/3.0)* (almBi)**2 *d_temp3)
         stiff(:, 14, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,4) )   
         stiff(:, 15, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,4) * u1 + (-2.0/3.0)*d(:,6) * u3 )
         stiff(:, 15, 8) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,4) * u2 + (-2.0/3.0)*d(:,2) * u3 -
     &                    (5.0/3.0)* (almBi)**2 *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3))
         stiff(:, 15, 9) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,6) * u1 + d(:,2) * u2 + (4.0/3.0)*d(:,4) * u3)
c     
c.... K33
c     
         stiff(:, 12, 12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,3) )
         stiff(:, 12, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( d(:,5)-(5.0/3.0)* (almBi)**2 * d(:,5)  )
         stiff(:, 13, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,3) )
         stiff(:, 13, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                     ( d(:,4) -(5.0/3.0)* (almBi)**2 * d(:,4))
         stiff(:, 14, 12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( -(2.0/3.0)*d(:,5) )
         stiff(:, 14, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( -(2.0/3.0)*d(:,4) )
         stiff(:, 14, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*
     &                    ( (4.0/3.0)*d(:,3) -(5.0/3.0)* (almBi)**2 * d_temp3 )   
         stiff(:, 15, 12) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,3) * u1 + (-2.0/3.0)*d(:,5) * u3 )
         stiff(:, 15, 13) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    ( d(:,3) * u2 + (-2.0/3.0)*d(:,4) * u3 )
         stiff(:, 15, 14) = ShearMod * (almBi * det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)* 
     &                    (d(:,5) * u1 + d(:,4) * u2 + (4.0/3.0)*d(:,3) * u3 -
     &                    (5.0/3.0)* (almBi)**2 *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3))
         stiff(:, 15, 15) = con ! notice the + or -
c         stiff(:,:,:) = -stiff(:,:,:)
c
c
c......end of modification for solid
        else
c......if the blocks are liquid and gas
c.... K11
c
         stiff(:, 2, 2) = rlm2mu
         stiff(:, 3, 3) = rmu
         stiff(:, 4, 4) = rmu
         stiff(:, 5, 2) = rlm2mu * u1
         stiff(:, 5, 3) = rmu    * u2
         stiff(:, 5, 4) = rmu    * u3
         stiff(:, 5, 5) = con
c     
c.... K12
c     
         stiff(:, 2, 8) = rlm
         stiff(:, 3, 7) = rmu
         stiff(:, 5, 7) = rmu    * u2
         stiff(:, 5, 8) = rlm    * u1
c     
c.... K13
c     
         stiff(:, 2,14) = rlm
         stiff(:, 4,12) = rmu
         stiff(:, 5,12) = rmu    * u3
         stiff(:, 5,14) = rlm    * u1
           
c     
c.... K21
c     
         stiff(:, 7, 3) = rmu
         stiff(:, 8, 2) = rlm
         stiff(:,10, 2) = rlm    * u2
         stiff(:,10, 3) = rmu    * u1
           
c     
c.... K22
c     
         stiff(:, 7, 7) = rmu
         stiff(:, 8, 8) = rlm2mu
         stiff(:, 9, 9) = rmu
         stiff(:,10, 7) = rmu    * u1
         stiff(:,10, 8) = rlm2mu * u2
         stiff(:,10, 9) = rmu    * u3
         stiff(:,10,10) = con
c     
c.... K23
c     
         stiff(:, 8,14) = rlm
         stiff(:, 9,13) = rmu
         stiff(:,10,13) = rmu    * u3
         stiff(:,10,14) = rlm    * u2
c     
c.... K31
c     
         stiff(:,12, 4) = rmu
         stiff(:,14, 2) = rlm
         stiff(:,15, 2) = rlm    * u3
         stiff(:,15, 4) = rmu    * u1
c     
c.... K32
c     
         stiff(:,13, 9) = rmu
         stiff(:,14, 8) = rlm
         stiff(:,15, 8) = rlm    * u3
         stiff(:,15, 9) = rmu    * u2
c     
c.... K33
c     
         stiff(:,12,12) = rmu
         stiff(:,13,13) = rmu
         stiff(:,14,14) = rlm2mu
         stiff(:,15,12) = rmu    * u1
         stiff(:,15,13) = rmu    * u2
         stiff(:,15,14) = rlm2mu * u3
         stiff(:,15,15) = con
c
        endif 
c....end for multilple materials
      endif

      if (itau .ge. 10) then     ! non-diag tau, buld compK
      
c     
c.... calculate element factors
c     
         f1 = dxidx(:,1,1) * dxidx(:,1,1) +
     &           dxidx(:,2,1) * dxidx(:,2,1) +
     &           dxidx(:,3,1) * dxidx(:,3,1)
         f2 = dxidx(:,1,1) * dxidx(:,1,2) +
     &           dxidx(:,2,1) * dxidx(:,2,2) +
     &           dxidx(:,3,1) * dxidx(:,3,2)
         f3 = dxidx(:,1,2) * dxidx(:,1,2) +
     &           dxidx(:,2,2) * dxidx(:,2,2) +
     &           dxidx(:,3,2) * dxidx(:,3,2)
         f4 = dxidx(:,1,1) * dxidx(:,1,3) +
     &           dxidx(:,2,1) * dxidx(:,2,3) +
     &           dxidx(:,3,1) * dxidx(:,3,3)
         f5 = dxidx(:,1,2) * dxidx(:,1,3) +
     &           dxidx(:,2,2) * dxidx(:,2,3) +
     &           dxidx(:,3,2) * dxidx(:,3,3)
         f6 = dxidx(:,1,3) * dxidx(:,1,3) +
     &           dxidx(:,2,3) * dxidx(:,2,3) +
     &           dxidx(:,3,3) * dxidx(:,3,3)
c     
c.... correction for tetrahedra (invariance w.r.t triangular coord)
c     
         if (lcsyst .eq. 1) then
            f1 = ( f1 + (dxidx(:,1,1) +
     &           dxidx(:,2,1) + dxidx(:,3,1))**2 ) * pt39
            f2 = ( f2 + (dxidx(:,1,1) +
     &              dxidx(:,2,1) + dxidx(:,3,1)) *
     &              (dxidx(:,1,2) +
     &              dxidx(:,2,2) + dxidx(:,3,2)) ) * pt39
            f3 = ( f3 + (dxidx(:,1,2) +
     &              dxidx(:,2,2) + dxidx(:,3,2))**2 ) * pt39
            f4 = ( f4 + (dxidx(:,1,1) +
     &              dxidx(:,2,1) + dxidx(:,3,1)) *
     &              (dxidx(:,1,3) +
     &              dxidx(:,2,3) + dxidx(:,3,3)) ) * pt39
            f5 = ( f5 + (dxidx(:,1,2) +
     &              dxidx(:,2,2) + dxidx(:,3,2)) *
     &              (dxidx(:,1,3) +
     &              dxidx(:,2,3) + dxidx(:,3,3)) ) * pt39
            f6 = ( f6 + (dxidx(:,1,3) +
     &              dxidx(:,2,3) + dxidx(:,3,3))**2 ) * pt39
c     
       !      flops = flops + 36*npro
         endif
c     
c.... correction for wedges
c     
         if ((lcsyst .eq. 3) .or. (lcsyst .eq. 4)) then
            f1 = ( dxidx(:,1,1) * dxidx(:,1,1) +
     &           dxidx(:,2,1) * dxidx(:,2,1) +
     &           ( dxidx(:,1,1) + dxidx(:,2,1) )**2 ) * pt57
     &           + dxidx(:,3,1) * dxidx(:,3,1)
            f2 = ( dxidx(:,1,1) * dxidx(:,1,2) +
     &           dxidx(:,2,1) * dxidx(:,2,2) +
     &           ( dxidx(:,1,1) + dxidx(:,2,1) ) *
     &           ( dxidx(:,1,2) + dxidx(:,2,2) ) ) * pt57
     &           + dxidx(:,3,1) * dxidx(:,3,2)
            f3 = ( dxidx(:,1,2) * dxidx(:,1,2) +
     &           dxidx(:,2,2) * dxidx(:,2,2) +
     &           ( dxidx(:,1,2) + dxidx(:,2,2) )**2 ) * pt57
     &           + dxidx(:,3,2) * dxidx(:,3,2)
            f4 = ( dxidx(:,1,1) * dxidx(:,1,3) +
     &           dxidx(:,2,1) * dxidx(:,2,3) +
     &           ( dxidx(:,1,1) + dxidx(:,2,1) ) *
     &           ( dxidx(:,1,3) + dxidx(:,2,3) ) ) * pt57
     &           + dxidx(:,3,1) * dxidx(:,3,3)
            f5 = ( dxidx(:,1,2) * dxidx(:,1,3) +
     &           dxidx(:,2,2) * dxidx(:,2,3) +
     &           ( dxidx(:,1,2) + dxidx(:,2,2) ) *
     &           ( dxidx(:,1,3) + dxidx(:,2,3) ) ) * pt57
     &           + dxidx(:,3,2) * dxidx(:,3,3)
            f6 = ( dxidx(:,1,3) * dxidx(:,1,3) +
     &           dxidx(:,2,3) * dxidx(:,2,3) +
     &           ( dxidx(:,1,3) + dxidx(:,2,3) )**2 ) * pt57
     &           + dxidx(:,3,3) * dxidx(:,3,3)
c     
       !      flops = flops + 51*npro
         endif
c     
c.... calculate compact K matrix in local parent coordinates
c.... equation 134 in "A new ... III" only w/ K^tilde_jj. Might need
c.... complete Kij.
      
         compK(:, 1) = f1 * T * rlm2mu + f3 * T * rmu
     &        + f6 * T * rmu
c     
         compK(:, 2) = f2 * T * (rlm + rmu)
         compK(:, 3) = f1 * T * rmu + f3 * T * rlm2mu
     &        + f6 * T * rmu
c     
         compK(:, 4) = f4 * T * (rlm + rmu)
         compK(:, 5) = f5 * T * (rlm + rmu)
         compK(:, 6) = f1 * T * rmu + f3 * T * rmu
     &        + f6 * T * rlm2mu
c     
         compK(:, 7) = f1 * T * rlm2mu  * u1 + f2 * T * (rlm + rmu) * u2
     &        + f3 * T * rmu * u1 + f4 * T * (rlm + rmu) * u3
     &        + f6 * T * rmu * u1
         compK(:, 8) = f1 * T * rmu * u2 + f2 * T * (rlm + rmu) * u1
     &        + f3 * T * rlm2mu  * u2 + f5 * T * (rlm + rmu) * u3
     &        + f6 * T * rmu * u2
         compK(:, 9) = f1 * T * rmu * u3 + f3 * T * rmu * u3
     &        + f4 * T * (rlm + rmu) * u1 + f5 * T * (rlm + rmu) * u2
     &        + f6 * T * rlm2mu  * u3

         rk=pt5*(u1**2+u2**2+u3**2)
         
         compK(:,10) = f1 * T * (con    * T + two * rmu * rk + (rlm +
     &        rmu) * u1**2) + f2 * T * (rlm + rmu) * two * u1 * u2 
     &        + f3 * T * (con    * T + two * rmu * rk + (rlm + rmu) *
     &        u2**2) + f4 * T * (rlm + rmu) * two * u1 * u3 
     &        + f5 * T * (rlm + rmu) * two * u2 * u3 + f6 * T * (con
     &        * T + two * rmu * rk + (rlm + rmu) * u3**2)  
c     
c.... flop count
c
    !      flops = flops + 86*npro
c     
c.... end of GLS
c     
      
      endif
c     
c.... --------------------------->  RHS  <-----------------------------
c     
c.... compute diffusive fluxes and add them to ri and rmi
c
c....if this is the solid block, modify the momentum and energy equation
        if (mat_eos(mater,1).eq.ieos_solid_1)then
c
         bq_af(:,:)= b_af(iblk_solid)%p(:,intp_s,:)
c.... diffusive flux in x1-direction

c         rmi(:,1) = zero ! already initialized
         rmi(:,2) = (det_baf)**(-5.0/6.0) * ShearMod * (1.0/3.0) * 
     &              (2 * bq_af(:,1) - bq_af(:,2) - bq_af(:,3))
         rmi(:,3) =  (det_baf)**(-5.0/6.0) * ShearMod * 
     &              bq_af(:,6)
         rmi(:,4) =  (det_baf)**(-5.0/6.0) * ShearMod *
     &              bq_af(:,5)
         rmi(:,5) =  (det_baf)**(-5.0/6.0) * ShearMod * 
     &              ( (1.0/3.0) * ( 2 * bq_af(:,1) - bq_af(:,2) - bq_af(:,3) ) *u1 + 
     &              bq_af(:,6) * u2 + bq_af(:,5) * u3 )
     &               + con      * g1yi(:,5)
c		 
        ri (:,2:5) = ri (:,2:5) + rmi(:,2:5)
c       rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x2-direction
c
c       rmi(:, 6) = zero
        rmi(:, 7) =  (det_baf)**(-5.0/6.0) * ShearMod * 
     &               bq_af(:,6)
        rmi(:, 8) =  (det_baf)**(-5.0/6.0) * ShearMod *
     &               (1.0/3.0) *(-bq_af(:,1) + 2.0* bq_af(:,2) - bq_af(:,3) )
        rmi(:, 9) =  (det_baf)**(-5.0/6.0) * ShearMod * bq_af(:,4)
        rmi(:,10) =  (det_baf)**(-5.0/6.0) * ShearMod * 
     &               (bq_af(:,6) * u1 + (1.0/3.0) * ( -bq_af(:,1)+
     &               2.0* bq_af(:,2) - bq_af(:,3) ) *u2 +
     &               bq_af(:,4) * u3 )
     &               + con      * g2yi(:,5)

c
      ri (:,7:10) = ri (:,7:10) + rmi(:,7:10)
c     rmi(:,7:10) = rmi(:,7:10) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x3-direction
c
c       rmi(:,11) = zero
        rmi(:,12) = (det_baf)**(-5.0/6.0) * ShearMod * 
     &               bq_af(:,5)
        rmi(:,13) = (det_baf)**(-5.0/6.0) * ShearMod * 
     &               bq_af(:,4)
        rmi(:,14) =  (det_baf)**(-5.0/6.0) * ShearMod *
     &               (1.0/3.0) *(-bq_af(:,1) - bq_af(:,2) + 2.0* bq_af(:,3) )
        rmi(:,15) =  (det_baf)**(-5.0/6.0) * ShearMod * 
     &               (bq_af(:,5) * u1 + (1.0/3.0) * ( -bq_af(:,1) - bq_af(:,2) +
     &               2.0* bq_af(:,3) ) *u3 +
     &               bq_af(:,4) * u2 )
     &               + con      * g3yi(:,5)

c
       ri (:,12:15) = ri (:,12:15) + rmi(:,12:15)
c!      flops = flops + 74*npro
      else
c .... other material(liquid or gas)
c.... diffusive flux in x1-direction
c
c       rmi(:,1) = zero ! already initialized
        rmi(:,2) =  rlm2mu      * g1yi(:,2) 
     &               +      rlm * g2yi(:,3) 
     &               +      rlm * g3yi(:,4)
     &               -      rlsli(:,1)
        rmi(:,3) =  rmu         * g1yi(:,3) 
     &               +      rmu * g2yi(:,2) 
     &               -      rlsli(:,4)
        rmi(:,4) =  rmu         * g1yi(:,4)
     &               +      rmu * g3yi(:,2)
     &               -      rlsli(:,5)
        rmi(:,5) =  rlm2mu * u1 * g1yi(:,2) + rmu * u2 * g1yi(:,3)
     &                                   +    rmu * u3 * g1yi(:,4)
     &               + rmu * u2 * g2yi(:,2) + rlm * u1 * g2yi(:,3)
     &               + rmu * u3 * g3yi(:,2) + rlm * u1 * g3yi(:,4)
     &               + con      * g1yi(:,5)

c
      ri (:,2:5) = ri (:,2:5) + rmi(:,2:5)
c     rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x2-direction
c
c       rmi(:, 6) = zero
        rmi(:, 7) =       rmu * g1yi(:,3) 
     &             +      rmu * g2yi(:,2)
     &             -      rlsli(:,4)
        rmi(:, 8) =       rlm * g1yi(:,2)
     &             +   rlm2mu * g2yi(:,3)
     &             +      rlm * g3yi(:,4)
     &             -      rlsli(:,2)
        rmi(:, 9) =       rmu * g2yi(:,4)
     &             +      rmu * g3yi(:,3)
     &             -      rlsli(:,6)
        rmi(:,10) =  rlm * u2 * g1yi(:,2) +    rmu * u1 * g1yi(:,3)
     &             + rmu * u1 * g2yi(:,2) + rlm2mu * u2 * g2yi(:,3)
     &             + rmu * u3 * g2yi(:,4)
     &             + rmu * u3 * g3yi(:,3) +    rlm * u2 * g3yi(:,4)
     &             +    con   * g2yi(:,5)
c
      ri (:,7:10) = ri (:,7:10) + rmi(:,7:10)
c     rmi(:,7:10) = rmi(:,7:10) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x3-direction
c
c       rmi(:,11) = zero
        rmi(:,12) =       rmu * g1yi(:,4)
     &             +      rmu * g3yi(:,2)
     &             -      rlsli(:,5)
        rmi(:,13) =       rmu * g2yi(:,4)
     &             +      rmu * g3yi(:,3)
     &             -      rlsli(:,6)
        rmi(:,14) =       rlm * g1yi(:,2)
     &             +      rlm * g2yi(:,3)
     &             +   rlm2mu * g3yi(:,4)
     &             -   rlsli(:,3)
        rmi(:,15) =     rlm * u3 * g1yi(:,2) + rmu * u1 * g1yi(:,4)
     &             +    rlm * u3 * g2yi(:,3) + rmu * u2 * g2yi(:,4)
     &             +    rmu * u1 * g3yi(:,2) + rmu * u2 * g3yi(:,3)
     &             + rlm2mu * u3 * g3yi(:,4)
     &             +    con      * g3yi(:,5) 
c
      ri (:,12:15) = ri (:,12:15) + rmi(:,12:15)
c
      endif
c
c!      flops = flops + 74*npro
c
c  stiff for preconditioner has been eliminated
c  all preconditioner stuff is in e3bdg.f
c

      ttim(23) = ttim(23) + secs(0.0)
      
c
c.... return
c
        return
        end
c     
c     
c     
      subroutine e3viscSclr (g1yti,    g2yti,    g3yti,
     &     rmu,      Sclr,     rho,
     &     rti,      rmti,     stifft)
c     
c----------------------------------------------------------------------
c     
c     This routine calculates the contribution of viscous and heat fluxes
c     to both RHS and LHS.
c     
c     input:
c     g1yti  (npro)              : grad-y in direction 1
c     g2yti  (npro)              : grad-y in direction 2
c     g3yti  (npro)              : grad-y in direction 3
c     rmu    (npro)              : viscosity
c     Sclr   (npro)              : scalar variable
c     
c     output:
c     rti     (npro,nsd+1)       : partial residual
c     rmti    (npro,nsd+1)       : partial modified residual
c     stifft  (npro,nsd,nsd)     : stiffness matrix
c     
c     
c     
c     Zdenek Johan, Summer 1990. (Modified from e2visc.f)
c     Zdenek Johan, Winter 1991. (Fortran 90)
c     Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c     
      use turbSA                ! for saSigma
      include "common.h"
c     
c     passed arrays
c     
      dimension g1yti(npro),           g2yti(npro),
     &     g3yti(npro),           rmu(npro),
     &     Sclr(npro),            rho(npro),
     &     rti(npro,nsd+1),       rmti(npro,nsd+1),
     &     stifft(npro,3,3)

      ttim(23) = ttim(23) - tmr()

      if ((ilset.ne.zero) .and. (isclr.lt.3)) return 
c     
c.... --------------------------->  RHS  <-----------------------------
c     
c.... --------------------->  Diffusivity Matrix  <-------------------
c     
c      if (lhs .eq. 1) then

         stifft = zero
c     
c.... K11
c     
         stifft(:,1,1)=rmu
         if (iRANS .lt. 0) then
            stifft(:,1,1)=saSigmaInv*rho*((rmu/rho)+max(zero,Sclr))
         endif
c     
c.... K22
c     
         stifft(:,2,2)=stifft(:,1,1)
c     
c.... K33
c     
         stifft(:,3,3)=stifft(:,1,1)

c      endif
c     
c.... --------------------------->  RHS  <-----------------------------
c     
c.... compute diffusive fluxes and add them to ri and rmi
c     
c.... diffusive fluxes
c     
      rmti(:,1) = stifft(:,1,1) * g1yti(:) 
      rmti(:,2) = stifft(:,2,2) * g2yti(:) 
      rmti(:,3) = stifft(:,3,3) * g3yti(:) 

      rti (:,:) = rti (:,:) + rmti(:,:)
c     rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c     
c!      flops = flops + 74*npro
c     
      ttim(23) = ttim(23) + tmr()

c     
c.... return
c     
      return
      end
