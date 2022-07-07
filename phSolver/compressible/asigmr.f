        subroutine AsIGMR (y,       ac,      x,          xmudmi,   
     &                     shp,     shgl,    ien,     
     &                     mater,   res,     rmes,    
     &                     BDiag,   qres,    EGmass,   
     &                     rerr,    umesh)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use rlssave     ! Use the resolved Leonard stresses at the nodes.
        use timedataC    ! time series
        use specialBC    ! get ytarget to localize and send down
        include "common.h"
c
        dimension y(nshg,ndofl),            ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,MAXQPT),  
     &            shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),  
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),         BDiag(nshg,nflow,nflow),
     &            qres(nshg,idflx)
      integer, intent(in) :: mater

c
        dimension ycl(npro,nshl,ndofl),     acl(npro,nshl,ndof),
     &            xl(npro,nenl,nsd),        ytargetl(npro,nshl,nflow),
     &            rl(npro,nshl,nflow),      rml(npro,nshl,nflow),
     &            BDiagl(npro,nshl,nflow,nflow),
     &            ql(npro,nshl,idflx)
c        
        dimension  xmudmi(npro,ngauss)
        dimension sgn(npro,nshl),  EGmass(npro,nedof,nedof)
c
        dimension umesh(numnp, nsd),  uml(npro,nshl,nsd)
c
        dimension rlsl(npro,nshl,6) 
        real*8 rerrl(npro,nshl,6), rerr(nshg,10)
c
#define debug 0
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
c
c.... gather the variables
c
        call localy(y,      ycl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call local (qres,   ql,     ien,    idflx,  'gather  ')
        call local (umesh,  uml,    ien,    nsd,    'gather  ')

        if(matflg(5,1).ge.4 )
     &   call localy (ytarget,   ytargetl,  ien,   nflow,  'gather  ')


        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero
        BDiagl = zero

        if(ierrcalc.eq.1) rerrl = zero
        ttim(31) = ttim(31) - secs(0.0)

        call e3  (ycl,     ycl,     acl,     shp,
     &            shgl,    xl,      rl,      rml,   xmudmi,
     &            BDiagl,  ql,      sgn,     rlsl,  EGmass,
     &            rerrl,   ytargetl, uml)

        ttim(31) = ttim(31) + secs(0.0)
c
#if debug==1
      i = 97
      do iel = 1,npro
        do n = 1,nshl
          if (ien(iel,n) .ne. i) cycle
          write(*,20) 'ASIGMR res before: ',ien(iel,n),res(ien(iel,n),:)
        enddo
      enddo
      do iel=1,npro
        do n = 1,nshl
          if (ien(iel,n) .ne. i) cycle
          write(*,10) 'ASIGMR rl:  ',iel,n,ien(iel,n),rl(iel,n,:)
        enddo
      enddo
#endif
c
c.... assemble the residual and modified residual
c
        call local (res,    rl,     ien,    nflow,  'scatter ')
c
#if debug ==1
      do iel = 1,npro
        do n = 1,nshl
          if (ien(iel,n) .ne. i) cycle
          write(*,20) 'ASIGMR res after: ',ien(iel,n),res(ien(iel,n),:)
        enddo
      enddo
#endif
10    format(a12,3i6,5e24.16)
20    format(a12,1i6,5e24.16)
c
        if ( ierrcalc .eq. 1 ) then
           call local (rerr, rerrl,  ien, 6, 'scatter ')
        endif
c
c.... extract and assemble the Block-Diagonal (see note in elmgmr, line 280)
c
        if (iprec .ne. 0) then 
           do i = 1, nshl
              do j = 1, nflow
                 i0 = (i - 1) * nflow + j
                 do k = 1, nflow
                    j0 = (i - 1) * nflow + k
                    BDiagl(:,i,j,k) = EGmass(:,i0,j0)
                 enddo
              enddo
           enddo
           call local (BDiag,  BDiagl, ien, nflow*nflow, 'scatter ')
        endif
        
c
c... call timeseries
c

        if (exts) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(ycl,xl,ien,sgn)
           endif
        endif
        
c
c.... end
c
        return
        end
c
c
c
        subroutine AsIGMRSclr (y,       ac,
     &                         x,       elDwl,
     &                         shp,     shgl,      ien,     
     &                         mater,   rest,      rmest,    
     &                         qrest,   EGmasst,   Diag,
     &                         umesh)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use turbSA
        include "common.h"
c
        dimension y(nshg,ndof),             
     &            ac(nshg,ndof),
     &            x(numnp,nsd),              
     &            shp(nshl,MAXQPT),        shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),
     &            mater(npro),            rest(nshg),
     &            rmest(nshg),            Diag(nshg),
     &            qrest(nshg)

c
        dimension ycl(npro,nshl,ndof),       
     &            acl(npro,nshl,ndof),    dwl(npro,nenl),
     &            xl(npro,nenl,nsd),      Diagl(npro,nshl),
     &            rtl(npro,nshl),         rmtl(npro,nshl),
     &            qtl(npro,nshl),         sgn(npro,nshl)
c        
        dimension EGmasst(npro,nshape, nshape)
        real*8    elDwl(npro)
c
        real*8, dimension(numnp,nsd), intent(in) :: umesh
c
        real*8, dimension(npro,nshl,nsd) :: uml
#define debugsclr 0
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        call getsgn(ien,sgn)
c
c
c.... gather the variables
c
        call localy (y,       ycl,      ien,    ndof,  'gather  ')
        call localy (ac,      acl,     ien,    ndof,  'gather  ')
        call localx (x,       xl,      ien,    nsd,   'gather  ')
        call local (umesh,  uml,    ien,    nsd,    'gather  ')
c       call local (qrest,   qtl,     ien,    1,     'gather  ')
        if (iRANS .lt. 0) then
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rtl     = zero
        Diagl   = zero
        
        ttim(31) = ttim(31) - tmr()

        call e3Sclr (ycl,     acl,  
     &               dwl,    elDwl,   shp,
     &               sgn,    shgl,    xl,
     &               rtl,    rmtl,
     &               qtl,    EGmasst,
     &               uml )

        ttim(31) = ttim(31) + tmr()
c
#if debugsclr == 1
      i = 172
      do iel = 1,npro
        do n = 1,nshl
          if (ien(iel,n) .ne. i) cycle
          write(*,20) 'ASIGMRSclr rest before: ',ien(iel,n),rest(ien(iel,n))
        enddo
      enddo
      do iel = 1,npro
        do n = 1,nshl
          if (ien(iel,n) .ne. i) cycle
          write(*,10) 'ASIGMRSclr rtl:  ',iel,n,ien(iel,n),rtl(iel,n)
        enddo
      enddo
#endif
c
c.... assemble the residual and modified residual
c
        call local (rest,    rtl,     ien,    1,  'scatter ')
c
#if debugsclr == 1
      i = 172
      do iel = 1,npro
        do n = 1,nshl
          if (ien(iel,n) .ne. i) cycle
          write(*,20) 'ASIGMRSclr rest after: ',ien(iel,n),rest(ien(iel,n))
        enddo
      enddo
#endif
c
10    format(a18,3i6,e24.16)
20    format(a18,1i6,e24.16)
c
c.... extract and assemble the Diagonal
c
        if (iprec .ne. 0) then 
           do i=1,nshl
              Diagl(:,i)=EGmassT(:,i,i)
           enddo
           call local(Diag, Diagl, ien, 1, 'scatter ')
        endif
c
c.... end
c
        return
        end


