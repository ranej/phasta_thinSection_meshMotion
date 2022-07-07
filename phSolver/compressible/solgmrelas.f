      subroutine SolGMRElas (x,       disp,    iBC,     BC,
     &                       col,     row,     meshq,
     &                       HBrg,    eBrg,    yBrg,
     &                       Rcos,    Rsin,    iper,    ilwork,
     &                       shp,     shgl,    shpb,    shglb,
     &                       shpif,   Dy)
c
c----------------------------------------------------------------------
c
c  This is the preconditioned GMRES driver routine for the mesh-elastic solve
c
c input:
c  x      (numnp,nsd)            : node coordinates
c  disp   (numnp,nsd)            : mesh displacement
c  iBC    (nshg)                 : BC codes
c  BC     (nshg,ndofBC)          : BC constraint parameters
c  shp(b) (nen,maxsh,melCat)     : element shape functions (boundary)
c  shgl(b)(nsd,nen,maxsh,melCat) : local gradients of shape functions
c  
c output:
c  Dy     (nshg,nelas)           : mesh displacement solution
c
c----------------------------------------------------------------------
c
      use pointer_data
c        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c
      integer col(nshg+1), row(nnz*nshg)
      real*8  elaslhsK(nelas*nelas,nnz_tot), 
     &        meshq(numel)
c
      real*8  beta,  summed
c
      dimension x(numnp,nsd),             disp(numnp,nsd),
     &          iBC(nshg),                BC(nshg,ndofBC),
     &          elasres(nshg,nelas), 
     &          elasBDiag(nshg,nelas,nelas),
     &          HBrg(Kspace+1,Kspace),    eBrg(Kspace+1),
     &          yBrg(Kspace),             Rcos(Kspace),
     &          Rsin(Kspace),             ilwork(nlwork),
     &          iper(nshg)
c
      dimension Dy(nshg,nelas),  
     &          temp(nshg,nelas),
     &          uBrg(nshg,nelas,Kspace+1)
c
      dimension shp(MAXTOP,maxsh,MAXQPT),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT)
c
      real*8, dimension(maxtop,    maxsh,maxqpt) :: shpif
c
c.... *******************>> Element Data Formation <<******************
c
c.... form the LHS matrices, the residual vector, and the block
c     diagonal preconditioner
c
      call ElmGMRElas (x,       disp,    shp,     shgl,
     &                 iBC,     BC,      shpb,    shglb,
     &                 shpif,   elasres, elasBDiag,
     &                 iper,    ilwork,  elaslhsK,
     &                 col,     row,     meshq)
c
c.... **********************>>    EBE - GMRES    <<********************
c
CC      call timer ('Solver  ')
c
c.... ------------------------> Initialization <-----------------------
c
c
c.... LU decompose the block diagonals
c
      if (iprec .ne. 0) then
         call i3LU (elasBDiag, elasres, nelas,  'LU_Fact ')
         if (numpe > 1) then
            call commu (elasBDiag  , ilwork, nelas*nelas  , 'out')
         endif
      endif
c
c.... block diagonal precondition residual
c
      call i3LU (elasBDiag, elasres, nelas, 'forward ')
c
c.... initialize Dy
c
      Dy = zero
c
c.... Pre-precondition the LHS mass matrix and set up the sparse 
c     preconditioners
c

      if(lhs.eq.1) call Spsi3pre (elasBDiag, elaslhsK, nelas, col, row)
c     
c.... copy res in uBrg(1)
c     
      uBrg(:,:,1) = elasres
c     
c.... calculate norm of residual
c
      temp  = elasres**2

      call sumgat (temp, nelas, summed)
      unorm = sqrt(summed)
c
c.... check if GMRES iterations are required
c
      iKs    = 0
      lGMRES = 0
c
c.... if we are down to machine precision, dont bother solving
c
      if (unorm .lt. 100.*epsM**2) goto 3000 
c
c.... set up tolerance of the Hessenbergs problem
c
      epsnrm = etolelas * unorm
c
c.... ------------------------>  GMRES Loop  <-------------------------
c
c.... loop through GMRES cycles
c
      do 2000 mGMRES = 1, nGMRES
         lGMRES = mGMRES - 1
c
         if (lGMRES .gt. 0) then
c
c.... if GMRES restarts are necessary, calculate  R - A x
c
c
c.... right precondition Dy
c
            temp = Dy           
c
c.... perform the A x product
c
            call SparseAp (iper, ilwork, iBC,      nelas,
     &                     col,  row,    elaslhsK, temp)
c
c.... subtract A x from residual and calculate the norm
c           
            temp = elasres - temp
            uBrg(:,:,1) = temp
c
c.... calculate the norm
c
            temp  = temp**2
            call sumgat (temp, nelas, summed)
            unorm = sqrt(summed)
c     
c.... flop count
c     
            flops = flops + nelas*nshg+nshg
c     
         endif
c
c.... set up RHS of the Hessenbergs problem
c
         call clear (eBrg, Kspace+1)
         eBrg(1) = unorm
c
c.... normalize the first Krylov vector
c
         uBrg(:,:,1) = uBrg(:,:,1) / unorm
c
c.... loop through GMRES iterations
c
         do 1000 iK = 1, Kspace
            iKs = iK

            uBrg(:,:,iKs+1) = uBrg(:,:,iKs)
c
c.... Au product  ( u_{i+1} <-- EGmass u_{i+1} )
c
            call SparseAp (iper,      ilwork, iBC,
     &                     nelas,     col,    row, 
     &                     elaslhsK,  uBrg(:,:,iKs+1) )
c
c.... orthogonalize and get the norm
c
            do jK = 1, iKs+1  
c
               if (jK .eq. 1) then
c
                  temp = uBrg(:,:,iKs+1) * uBrg(:,:,1) ! {u_{i+1}*u_1} vector 
                  call sumgat (temp, nelas, beta) ! sum vector=(u_{i+1},u_1)
c
               else
c
c project off jK-1 vector
c
                  uBrg(:,:,iKs+1)=uBrg(:,:,iKs+1)-beta * uBrg(:,:,jK-1)
c
                  temp = uBrg(:,:,iKs+1) * uBrg(:,:,jK) !{u_{i+1}*u_j} vector
                  call sumgat (temp, nelas, beta) ! sum vector=(u_{i+1},u_j)
c
               endif
c
               HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix
c
            enddo
c
            flops = flops + (3*iKs+1)*nelas*numnp+(iKs+1)*numnp
c
c  the last inner product was with what was left of the vector (after
c  projecting off all of the previous vectors
c
            unorm           = sqrt(beta)
            HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band
c
c.... normalize the Krylov vector
c
            uBrg(:,:,iKs+1) = uBrg(:,:,iKs+1) / unorm ! normalize the next Krylov
c vector
c
c.... construct and reduce the Hessenberg Matrix
c  since there is only one subdiagonal we can use a Givens rotation to 
c  rotate off each subdiagonal AS IT IS FORMED.   We do this because it
c  allows us to check progress of solution and quit when satisfied.  Note
c  that all future K vects will put a subdiagonal in the next column so
c  there is no penalty to work ahead as  the rotation for the next vector
c  will be unaffected by this rotation.
        
c     
c     H Y = E ========>   R_i H Y = R_i E
c     
            do jK = 1, iKs-1
               tmp            =  Rcos(jK) * HBrg(jK,  iKs) +
     &                           Rsin(jK) * HBrg(jK+1,iKs)
               HBrg(jK+1,iKs) = -Rsin(jK) * HBrg(jK,  iKs) +
     &                           Rcos(jK) * HBrg(jK+1,iKs)
               HBrg(jK,  iKs) =  tmp
            enddo
c     
            tmp            = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
            Rcos(iKs)      = HBrg(iKs,  iKs) / tmp
            Rsin(iKs)      = HBrg(iKs+1,iKs) / tmp
            HBrg(iKs,  iKs)= tmp
            HBrg(iKs+1,iKs)= zero
c     
c.... rotate eBrg    R_i E
c     
            tmp        = Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
            eBrg(iKs+1)=-Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
            eBrg(iKs)  = tmp
c     
c.... check for convergence
c     
            ntotGMelas = ntotGMelas + 1
            echeck=abs(eBrg(iKs+1))
            if (echeck .le. epsnrm.and. iKs .ge. minIters) exit
c     
c.... end of GMRES iteration loop
c     
 1000    continue
c
c.... ------------------------->   Solution   <------------------------
c
c.... if converged or end of Krylov space
c
c.... solve for yBrg
c
         do jK = iKs, 1, -1
            yBrg(jK) = eBrg(jK) / HBrg(jK,jK)
            do lK = 1, jK-1
               eBrg(lK) = eBrg(lK) - yBrg(jK) * HBrg(lK,jK)
            enddo
         enddo
c     
c.... update Dy
c
         do jK = 1, iKs
            Dy = Dy + yBrg(jK) * uBrg(:,:,jK)
         enddo
c     
c.... flop count
c
         flops = flops + (3*iKs+1)*nelas*nshg
c
c.... check for convergence
c     
        echeck=abs(eBrg(iKs+1))
        if (echeck .le. epsnrm) exit
        if(myrank.eq.master) write(*,*)'solver tolerance %satisfaction',
     &  (one-echeck*etolelas/epsnrm)/(one-etolelas)*100

c     
c.... end of mGMRES loop
c
 2000 continue
c     
c.... ------------------------>   Converged   <------------------------
c     
c.... if converged
c     
 3000 continue

c
c     
c.... back block-diagonal precondition the results 
c
      call i3LU (elasBDiag, Dy, nelas, 'backward')
c     
c
c.... output the statistics
c
      call rstatElas (elasres, ilwork) 
c    
c.... stop the timer
c     
 3002 continue                  ! no solve just res.
CC      call timer ('Back    ')
c     
c.... end
c     
      return
      end
