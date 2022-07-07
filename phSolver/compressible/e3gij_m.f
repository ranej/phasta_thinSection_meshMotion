      module e3gij_m
c-----------------------------------------------------------------------
c get the metric tensor g_{ij}=xi_{k,i} xi_{k,j},namely gijd,
c and its inverse g^{ij} = (g_{ij})^(-1), namely giju
c-----------------------------------------------------------------------
        implicit none
c
c
      contains
        subroutine e3gijd( dxidx, npro, nsd, lcsyst, gijd )
c-----------------------------------------------------------------------
c get the metric tensor g_{ij}=xi_{k,i} xi_{k,j}. 
c could be improved by use modules to take care of npro, nsd,lcsyst 
c-----------------------------------------------------------------------        
        
          implicit none
c   
        
          real*8  dxidx(npro,nsd,nsd),  gijd(npro,6),
     &            tmp1(npro),           tmp2(npro),
     &            tmp3(npro)
          integer :: npro,nsd,lcsyst    
          real*8 :: c1, c2 
c  form metric tensor g_{ij}=xi_{k,i} xi_{k,j}.  It is a symmetric
c  tensor so we only form 6 components and use symmetric matrix numbering.
c  (d for down since giju=[gijd]^{-1})
c  (Note FARZIN and others use numbering of 1,2,3 being diagonal 456 off)
          if (lcsyst .ge. 2) then
  
            gijd(:,1) = dxidx(:,1,1) * dxidx(:,1,1)
     &                + dxidx(:,2,1) * dxidx(:,2,1)
     &                + dxidx(:,3,1) * dxidx(:,3,1)
c
            gijd(:,2) = dxidx(:,1,1) * dxidx(:,1,2)
     &                + dxidx(:,2,1) * dxidx(:,2,2)
     &                + dxidx(:,3,1) * dxidx(:,3,2)
c
            gijd(:,3) = dxidx(:,1,2) * dxidx(:,1,2)
     &                + dxidx(:,2,2) * dxidx(:,2,2)
     &                + dxidx(:,3,2) * dxidx(:,3,2)
c
            gijd(:,4) = dxidx(:,1,1) * dxidx(:,1,3)
     &                + dxidx(:,2,1) * dxidx(:,2,3)
     &                + dxidx(:,3,1) * dxidx(:,3,3)
c
            gijd(:,5) = dxidx(:,1,2) * dxidx(:,1,3)
     &                + dxidx(:,2,2) * dxidx(:,2,3)
     &                + dxidx(:,3,2) * dxidx(:,3,3)
c
            gijd(:,6) = dxidx(:,1,3) * dxidx(:,1,3)
     &               + dxidx(:,2,3) * dxidx(:,2,3)
     &               + dxidx(:,3,3) * dxidx(:,3,3)
c
          else   if (lcsyst .eq. 1) then   
c
c  There is an invariance problem with tets 
c  It is fixed by the following modifications to gijd 
c

            c1 = 1.259921049894873D+00
            c2 = 6.299605249474365D-01
c
            tmp1(:) = c1 * dxidx(:,1,1) + c2 *(dxidx(:,2,1)+dxidx(:,3,1))
            tmp2(:) = c1 * dxidx(:,2,1) + c2 *(dxidx(:,1,1)+dxidx(:,3,1))
            tmp3(:) = c1 * dxidx(:,3,1) + c2 *(dxidx(:,1,1)+dxidx(:,2,1))
            gijd(:,1) = dxidx(:,1,1) * tmp1
     1                + dxidx(:,2,1) * tmp2
     2                + dxidx(:,3,1) * tmp3
c
            tmp1(:) = c1 * dxidx(:,1,2) + c2 *(dxidx(:,2,2)+dxidx(:,3,2))
            tmp2(:) = c1 * dxidx(:,2,2) + c2 *(dxidx(:,1,2)+dxidx(:,3,2))
            tmp3(:) = c1 * dxidx(:,3,2) + c2 *(dxidx(:,1,2)+dxidx(:,2,2))
            gijd(:,2) = dxidx(:,1,1) * tmp1
     1                + dxidx(:,2,1) * tmp2
     2                + dxidx(:,3,1) * tmp3
c
            gijd(:,3) = dxidx(:,1,2) * tmp1
     1                + dxidx(:,2,2) * tmp2
     2                + dxidx(:,3,2) * tmp3
c
            tmp1(:) = c1 * dxidx(:,1,3) + c2 *(dxidx(:,2,3)+dxidx(:,3,3))
            tmp2(:) = c1 * dxidx(:,2,3) + c2 *(dxidx(:,1,3)+dxidx(:,3,3))
            tmp3(:) = c1 * dxidx(:,3,3) + c2 *(dxidx(:,1,3)+dxidx(:,2,3))
            gijd(:,4) = dxidx(:,1,1) * tmp1
     1                + dxidx(:,2,1) * tmp2
     2                + dxidx(:,3,1) * tmp3
c
            gijd(:,5) = dxidx(:,1,2) * tmp1
     1                + dxidx(:,2,2) * tmp2
     2                + dxidx(:,3,2) * tmp3
c
            gijd(:,6) = dxidx(:,1,3) * tmp1
     1                + dxidx(:,2,3) * tmp2
     2                + dxidx(:,3,3) * tmp3
c
          else  
c  This is just the hex copied from above.  I have
c  to find my notes on invariance factors for wedges
c         

            gijd(:,1) = dxidx(:,1,1) * dxidx(:,1,1)
     &                + dxidx(:,2,1) * dxidx(:,2,1)
     &                + dxidx(:,3,1) * dxidx(:,3,1)
c
            gijd(:,2) = dxidx(:,1,1) * dxidx(:,1,2)
     &                + dxidx(:,2,1) * dxidx(:,2,2)
     &                + dxidx(:,3,1) * dxidx(:,3,2)
c
            gijd(:,3) = dxidx(:,1,2) * dxidx(:,1,2)
     &                + dxidx(:,2,2) * dxidx(:,2,2)
     &                + dxidx(:,3,2) * dxidx(:,3,2)
c
            gijd(:,4) = dxidx(:,1,1) * dxidx(:,1,3)
     &                + dxidx(:,2,1) * dxidx(:,2,3)
     &                + dxidx(:,3,1) * dxidx(:,3,3)
c
            gijd(:,5) = dxidx(:,1,2) * dxidx(:,1,3)
     &                + dxidx(:,2,2) * dxidx(:,2,3)
     &                + dxidx(:,3,2) * dxidx(:,3,3)
c
            gijd(:,6) = dxidx(:,1,3) * dxidx(:,1,3)
     &                + dxidx(:,2,3) * dxidx(:,2,3)
     &                + dxidx(:,3,3) * dxidx(:,3,3)
          endif
c
          return
        end subroutine e3gijd
c
        subroutine e3giju (giju, dxidx, npro, nsd, lcsyst)
c-----------------------------------------------------------------------
c get the inverse of metric tensor g^{ij}= (g_{ij})^(-1)
c could be improved by use modules to take care of npro, nsd,lcsyst   
c-----------------------------------------------------------------------        
          use number_def_m
          implicit none
c
          real*8, dimension(npro,nsd,nsd),intent(in) :: dxidx
          integer,intent(in) :: npro,nsd,lcsyst
          real*8, dimension(npro,6),intent(out) :: giju
c
          real*8, dimension(npro,6) :: gijdu
          real*8, dimension(npro,6) :: gijd
          real*8, dimension(npro) ::  detgijI
c
c.... for the notation of different numbering
c     
           call e3gijd( dxidx, npro, nsd, lcsyst, gijd )

           gijdu(:,1)=gijd(:,1)
           gijdu(:,2)=gijd(:,3)
           gijdu(:,3)=gijd(:,6)
           gijdu(:,4)=gijd(:,2)
           gijdu(:,5)=gijd(:,4)
           gijdu(:,6)=gijd(:,5)
c     
c     
           detgijI = one/(
     &          gijdu(:,1) * gijdu(:,2) * gijdu(:,3)
     &          - gijdu(:,1) * gijdu(:,6) * gijdu(:,6)
     &          - gijdu(:,4) * gijdu(:,4) * gijdu(:,3)
     &          + gijdu(:,4) * gijdu(:,5) * gijdu(:,6) * two
     &          - gijdu(:,5) * gijdu(:,5) * gijdu(:,2) 
     &          )
           giju(:,1) = detgijI * (gijdu(:,2)*gijdu(:,3) 
     &               - gijdu(:,6)**2)
           giju(:,2) = detgijI * (gijdu(:,1)*gijdu(:,3) 
     &               - gijdu(:,5)**2)
           giju(:,3) = detgijI * (gijdu(:,1)*gijdu(:,2)
     &               - gijdu(:,4)**2)
           giju(:,4) = detgijI * (gijdu(:,5)*gijdu(:,6)
     &               - gijdu(:,4)*gijdu(:,3) )
           giju(:,5) = detgijI * (gijdu(:,4)*gijdu(:,6)
     &               - gijdu(:,5)*gijdu(:,2) )
           giju(:,6) = detgijI * (gijdu(:,4)*gijdu(:,5)
     &               - gijdu(:,1)*gijdu(:,6) )                             
        
        
        end subroutine e3giju
c
      end module e3gij_m
