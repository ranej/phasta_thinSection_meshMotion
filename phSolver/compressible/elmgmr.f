        subroutine ElmGMRe (y,         ac,        x,         
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       rmes,      BDiag,
     &                     iper,      ilwork,    EGmass, rerr)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
        use timedataC
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),               
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),      BDiag(nshg,nflow,nflow),
     &            iper(nshg),           EGmass(numel,nedof,nedof)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg, idflx),     rmass(nshg)
c
        dimension ilwork(nlwork)

        real*8 Bdiagvec(nshg,nflow), rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

	ttim(80) = ttim(80) - secs(0.0)
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        ires   = 1
c
        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        rmass = zero
        
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIq (y,                x,                       
     &               tmpshp,              
     &               tmpshgl,
     &               mien(iblk)%p,     mxmudmi(iblk)%p,
     &               qres,                   
     &               rmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
       
c
c.... take care of periodic boundary conditions
c

       call qpbc( rmass, qres, iBC,  iper, ilwork )       
c
      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        res    = zero
        rmes   = zero ! to avoid trap_uninitialized
        if (lhs. eq. 1)   EGmass = zero
        if (iprec .ne. 0) BDiag = zero
        flxID = zero

c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk          ! used in timeseries
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,
     &                 tmpshp,
     &                 tmpshgl,             mien(iblk)%p,
     &                 mmat(iblk)%p,        res,
     &                 rmes,                BDiag,
     &                 qres,                EGmass(iel:inum,:,:),
     &                 rerr)
c
c.... satisfy the BC's on the implicit LHS
c     
          call bc3LHS (iBC,                  BC,  mien(iblk)%p, 
     &                 EGmass(iel:inum,:,:)  ) 

          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (y,                       x,
     &                 tmpshpb,                 tmpshglb, 
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     rmes)

          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
        enddo
c
      ttim(80) = ttim(80) + secs(0.0)
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
c      if(iabc==1)               !are there any axisym bc's
c     &     call rotabc(res(1,2), iBC, BC, nflow,  'in ')
      if(iabc==1) then               !are there any axisym bc's
          call rotabc(res(1,2), iBC,  'in ')
c          Bdiagvec(:,1)=BDiag(:,1,1)
c          Bdiagvec(:,2)=BDiag(:,2,2)
c          Bdiagvec(:,3)=BDiag(:,3,3)
c          Bdiagvec(:,4)=BDiag(:,4,4)
c          Bdiagvec(:,5)=BDiag(:,5,5)
c          call rotabc(Bdiagvec(1,2), iBC, BC, 2,  'in ')
c          BDiag(:,:,:)=zero
c          BDiag(:,1,1)=Bdiagvec(:,1)
c          BDiag(:,2,2)=Bdiagvec(:,2)
c          BDiag(:,3,3)=Bdiagvec(:,3)
c          BDiag(:,4,4)=Bdiagvec(:,4)
c          BDiag(:,5,5)=Bdiagvec(:,5)
       endif

c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (res  , ilwork, nflow  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        if(iprec .ne. 0) call commu (BDiag, ilwork, nflow*nflow, 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (y,  iBC,  BC,  res,  iper, ilwork)
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
      if (iprec .ne. 0) then
         call bc3BDg (y,  iBC,  BC,  BDiag, iper, ilwork)
      endif
c
c.... return
c
      call timer ('Back    ')
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc       SPARSE   
c_______________________________________________________________

        subroutine ElmGMRs (y,         ac,        x,         
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     shpif,     shgif,
     &                     res,       rmes,      BDiag,
     &                     iper,      ilwork,    lhsK,  
     &                     col,       row,       rerr,     umesh)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
        use e3_param_m
        use e3if_param_m
        use e3if_geom_m
        use e3if_func_m
        use solid_data_m
        use e3_solid_func_m
        use timedataC
        use if_velocity_m
        use if_global_m
        use eqn_state_m
        use e3_solid_m
        use probe_m
        use ifbc_def_m
        use weighted_normal_data_m !for weighted normal
        use interfaceflag ! for weighted normal
        use dgifinp_m, only: i_w_normal     
c
        include "common.h"
        include "mpif.h"
#define debug 0
c
#if debug==1 
      integer imax(5),idbg
#endif
c
        interface 
          subroutine e3if_setparam
     &    (
     &      nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_,
     &      npro_,ndof_,nsd_,nflow_,ipord_,nqpt_
     &    )
            use e3if_param_m
            implicit none
            integer, intent(in) :: nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_
            integer, intent(in) :: npro_,ndof_,nsd_,nflow_,ipord_,nqpt_
          end subroutine e3if_setparam
          subroutine e3if_setparam2
     &    (
     &     egmassif00,egmassif01,egmassif10,egmassif11,
     &     time_
     &    )
            use dgifinp_m
            use e3if_param_m
            implicit none
            real*8, dimension(:,:,:), allocatable, target, intent(in) :: egmassif00,egmassif01,egmassif10,egmassif11
            real*8, intent(in) :: time_
          end subroutine e3if_setparam2
          subroutine asidgif_geom
     &    (
     &     x,shpif0,shpif1,shgif0,shgif1,
     &     qwtif, qwtif0, qwtif1,
     &     ienif0, ienif1
     & )
            use hierarchic_m
            use local_m
            use e3if_geom_m
            use if_global_m
            implicit none
            real*8, intent(in) :: x(nshg,nsd)
            real*8, dimension(nshl0,nqpt),intent(in)   :: shpif0
            real*8, dimension(nshl1,nqpt),intent(in)   :: shpif1
            real*8, dimension(nsd,nshl0,nqpt), intent(in) :: shgif0
            real*8, dimension(nsd,nshl1,nqpt), intent(in) :: shgif1
            real*8, intent(in) :: qwtif(nqpt), qwtif0(nqpt), qwtif1(nqpt)
            integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
          end subroutine asidgif_geom
          subroutine asidgif
     &    (
     &     res,
     &     y,        x,       umesh,
     &     shpif0,   shpif1,  shgif0,  shgif1,
     &     qwtif, qwtif0,   qwtif1,
     &     ienif0,   ienif1,
     &     BDiag)
          use hierarchic_m
          use local_m
          use e3if_m
          use e3if_geom_m
          use if_global_m
          use conpar_m
          use weighted_normal_data_m, only:w_normal_l0, w_normal_l1, w_normal_global !for weighted normal on the interface
          use genpar_m, only: iprec
          implicit none
          real*8, dimension(nshg,nflow), intent(inout) :: res
          real*8, dimension(nshg,ndof),  intent(in)    :: y
          real*8, dimension(nshg,nsd),   intent(in)    :: x
          real*8, dimension(nshg, nsd), intent(inout) :: umesh
          real*8, dimension(nshl0,nqpt),intent(in)   :: shpif0
          real*8, dimension(nshl1,nqpt),intent(in)   :: shpif1
          real*8, dimension(nsd,nshl0,nqpt), intent(in)  :: shgif0
          real*8, dimension(nsd,nshl1,nqpt), intent(in)  :: shgif1
          real*8, dimension(nqpt), intent(in) :: qwtif, qwtif0, qwtif1
          real*8, dimension(nshg,nflow,nflow), intent(inout) :: BDiag
          integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
          end subroutine asidgif
          subroutine fillsparse_if
     &    ( lhsk,
     &      ienif0,ienif1,
     &      col,row,
     &      egmass,
     &      npro,
     &      nedof0, nedof1,
     &      nflow,nshg,nnz,nnz_tot)
            implicit none
            real*8, intent(inout) :: lhsK(nflow*nflow,nnz_tot)
            integer, dimension(:,:), pointer, intent(in) :: ienif0,ienif1
            integer, intent(in) :: col(nshg+1), row(nnz*nshg)
            real*8, intent(in) :: egmass(npro,nedof0,nedof1)
            integer, intent(in) :: nflow,nshg,nnz,nnz_tot,npro,nedof0,nedof1
          end subroutine fillsparse_if
        end interface
c
        integer col(nshg+1), row(nnz*nshg)
        real*8 lhsK(nflow*nflow,nnz_tot)
        
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),      BDiag(nshg,nflow,nflow),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        real*8, dimension(maxtop,    maxsh,maxqpt) :: shpif
        real*8, dimension(maxtop,nsd,maxsh,maxqpt) :: shgif
c
        dimension qres(nshg, idflx),     rmass(nshg)
c
        dimension ilwork(nlwork)
c  
        dimension umesh(numnp, nsd)
c
        real*8 Bdiagvec(nshg,nflow), rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
        real*8, allocatable :: EGmass(:,:,:)
c
        real*8, dimension(:,:,:), allocatable :: egmassif00,egmassif01,egmassif10,egmassif11
        real*8 :: length
c
        integer :: nedof0,nedof1
        integer, pointer, dimension(:,:) :: ienif0,ienif1
c
        ttim(80) = ttim(80) - secs(0.0)
c
c.... set up the timer
c
        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        ires   = 1
c
        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        rmass = zero
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nshl   = lcblk(10,iblk)
          mater  = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
          e3_malloc_ptr => e3_malloc
          e3_mfree_ptr => e3_mfree
c
          select case (mat_eos(mater,1))
          case (ieos_ideal_gas,ieos_ideal_gas_2)
            getthm6_ptr => getthm6_ideal_gas
            getthm7_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthm6_ptr => getthm6_ideal_gas_mixture
            getthm7_ptr => getthm7_ideal_gas_mixture
          case (ieos_liquid_1)
            getthm6_ptr => getthm6_liquid_1
            getthm7_ptr => getthm7_liquid_1
          case (ieos_solid_1)
            getthm6_ptr => getthm6_solid_1
            getthm7_ptr => getthm7_solid_1
            iblk_solid = iblk 
            e3_malloc_ptr => e3_malloc_solid
            e3_mfree_ptr => e3_mfree_solid
          case default
            call error ('getthm  ', 'wrong material', mater)
          end select
c
          if (associated(e3_malloc_ptr)) call e3_malloc_ptr
c
          call AsIq (y,                x,                       
     &               tmpshp,              
     &               tmpshgl,
     &               mien(iblk)%p,     mxmudmi(iblk)%p,
     &               qres,                   
     &               rmass)
c
          if (associated(e3_mfree_ptr)) call e3_mfree_ptr
c
          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
       
c
c.... take care of periodic boundary conditions
c

       call qpbc( rmass, qres, iBC, iper, ilwork )       
c
      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        res    = zero
        rmes   = zero ! to avoid trap_uninitialized
        if (lhs. eq. 1) lhsK = zero
        if (iprec .ne. 0) BDiag = zero
        flxID = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk          ! used in timeseries
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mater  = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c
          if(lhs.eq.1) then
             allocate (EGmass(npro,nedof,nedof))
             EGmass = zero
          else
             allocate (EGmass(1,1,1))
          endif

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
c
          e3_malloc_ptr => e3_malloc
          e3_mfree_ptr => e3_mfree
c
          select case (mat_eos(mater,1))
          case (ieos_ideal_gas,ieos_ideal_gas_2)
            getthm6_ptr => getthm6_ideal_gas
            getthm7_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthm6_ptr => getthm6_ideal_gas_mixture
            getthm7_ptr => getthm7_ideal_gas_mixture
          case (ieos_liquid_1)
            getthm6_ptr => getthm6_liquid_1
            getthm7_ptr => getthm7_liquid_1
          case (ieos_solid_1)
            getthm6_ptr => getthm6_solid_1
            getthm7_ptr => getthm7_solid_1
            iblk_solid = iblk 
            e3_malloc_ptr => e3_malloc_solid
            e3_mfree_ptr => e3_mfree_solid
          case default
            call error ('getthm  ', 'wrong material', mater)
          end select
c
          if (associated(e3_malloc_ptr)) call e3_malloc_ptr
c
          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,
     &                 tmpshp,
     &                 tmpshgl,             mien(iblk)%p,
     &                 mater,               res,
     &                 rmes,                BDiag,
     &                 qres,                EGmass,
     &                 rerr,                umesh)
c
          if(lhs.eq.1) then
c
c.... satisfy the BCs on the implicit LHS
c
             call bc3LHS (iBC,                  BC,  mien(iblk)%p, 
     &                    EGmass  ) 
c
c.... Fill-up the global sparse LHS mass matrix
c
             call fillsparseC( mien(iblk)%p, EGmass,
     1                        lhsK, row, col)
          endif
c
          if (associated(e3_mfree_ptr)) call e3_mfree_ptr
c
          deallocate ( EGmass )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mater  = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          if(lhs.eq.1 .and. iLHScond >= 1) then
             allocate (EGmass(npro,nshl,nshl))
             EGmass = zero
          else
             allocate (EGmass(1,1,1))
          endif
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)
c
          e3_malloc_ptr => e3_malloc
          e3_mfree_ptr => e3_mfree
c
          select case (mat_eos(mater ,1))
          case (ieos_ideal_gas,ieos_ideal_gas_2)
            getthm6_ptr => getthm6_ideal_gas
            getthm7_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthm6_ptr => getthm6_ideal_gas_mixture
            getthm7_ptr => getthm7_ideal_gas_mixture
          case (ieos_liquid_1)
            getthm6_ptr => getthm6_liquid_1
            getthm7_ptr => getthm7_liquid_1
          case (ieos_solid_1)
            getthm6_ptr => getthm6_solid_1
            getthm7_ptr => getthm7_solid_1
            iblk_solid = iblk 
            e3_malloc_ptr => e3_malloc_solid
            e3_mfree_ptr => e3_mfree_solid
          case default
            call error ('getthm  ', 'wrong material', mater )
          end select
c
          if (associated(e3_malloc_ptr)) call e3_malloc_ptr
c
          call AsBMFG (y,                       x,
     &                 tmpshpb,                 tmpshglb, 
     &                 mienb(iblk)%p,           mater ,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     rmes, 
     &                 EGmass,                  umesh)
c
          if(lhs == 1 .and. iLHScond > 0) then
            call fillSparseC_BC(mienb(iblk)%p, EGmass, 
     &                   lhsk, row, col)
          endif
c
          if (associated(e3_mfree_ptr)) call e3_mfree_ptr
c
          deallocate (EGmass)
          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
        enddo   !end of boundary element loop
c
#if debug==1
      imax = maxloc(res,1)
      idbg = imax(1)
      write(*,'(a22,i6,5e24.16)') 'ELMGMR: after boundary',idbg,res(idbg,:)
#endif
c
      ttim(80) = ttim(80) + secs(0.0)
c
c.... -------------------->   interface elements   <--------------------
c
c... loop over the interface element blocks
c    The first loop is for interface outward normal vector calculations on the interface
c    The second loop is for residual calculations...
c
        if (surface_tension_flag .eq. 1)  then
          allocate(if_kappa(nshg,nsd+1))
          if_kappa = zero
        else
          nullify(if_kappa)
        endif
c
        if_normal = zero
c
c... allocation and initialization of the weighted normal
        if (i_w_normal .eq. 1) then
           allocate(w_normal_global(nshg,nsd))
           allocate(length_temp(nshg))
c
           w_normal_global = zero
           length_temp = zero
        endif
c          
        if_blocks1: do iblk = 1, nelblif
c
          iel     = lcblkif(1, iblk)
          npro    = lcblkif(1,iblk+1) - iel
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          ipord   = lcblkif(5, iblk)    ! polynomial order
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          nenbl_if= lcblkif(8, iblk)    ! number of vertices on the interface
          mater0  = lcblkif(9, iblk)
          mater1  = lcblkif(10,iblk)
          nshl0   = lcblkif(iblkif_nshl0,iblk)
          nshl1   = lcblkif(iblkif_nshl1,iblk)
          itpid   = lcblkif(iblkif_topology,iblk)
          ngaussif = nintif(itpid)
c
          call e3if_setparam
     &    (
     &     nshg,nshl0,nshl1,nenl0,nenl1,lcsyst0,lcsyst1,
     &     npro,ndof,nsd,nflow,ipord,ngaussif
     &    )
     &   
c
          call e3if_geom_malloc
c
          ienif0 => mienif0(iblk)%p
          ienif1 => mienif1(iblk)%p
c
          call asidgif_geom
     &   (
     &    x,
     &    shpif(lcsyst0,1:nshl0,:), 
     &    shpif(lcsyst1,1:nshl1,:), 
     &    shgif(lcsyst0,1:nsd,1:nshl0,:),
     &    shgif(lcsyst1,1:nsd,1:nshl1,:),
     &    qwtif(itpid,:), qwtif(lcsyst0,:), qwtif(lcsyst1,:),
     &    ienif0, ienif1
     & )
c
c... preparation of getting weighted normal starts, major part of codes are copied
c    from the way to calculate gcnorml in elmgmrelas
c
          if (i_w_normal .eq. 1) then         
c... allocation and initialization of the factor
            allocate(calc_factor_temp(npro))
c
            calc_factor_temp(:) = 1 
c... copy from elmgmrelas, defining some global parameters:
c.... hack; DG interface doesn't support 2nd and higher order
c
            if(ipord .eq. 1) then
              nshlb   = nenbl_if ! only work with linear element
            else if(ipord .gt. 1) then
              write(*,*) "need to implement for higher order"
              call error('elmgmrelas  ','higher order', ipord)
            endif
c
c.... the 0 side
c
            lcsyst = lcsyst0 ! passing to the global variable, could be improved
            nenl = nenl0 ! same as above
            nshl = nshl0 ! same as above
            nenbl = nenbl_if ! same as above
c
            if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl_if ! may not be necessary
            ngaussb = ngaussif ! passing to the global variable, could be improved
c... calculate and assemble non-unit normal for side 0
            call calc_normal(x, shpif(lcsyst0,1:nshl0,:), 
     &                       calc_factor_temp, mienif0(iblk)%p, 
     &                       miBCB(iblk)%p, w_normal_global)
c.... end of the 0 side
c
c.... the 1 side
c
            lcsyst = lcsyst1 ! passing to the global variable, could be improved
            nenl = nenl1 ! same as above
            nshl = nshl1 ! same as above
            nenbl = nenbl_if ! same as above
c
            if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl_if ! may not be necessary
            ngaussb = ngaussif ! passing to the global variable, could be improved
c
c.... compute and assemble non-unit normal for side 1
c
            call calc_normal(x, shpif(lcsyst1,1:nshl1,:),
     &                       calc_factor_temp, mienif1(iblk)%p,
     &                       miBCB(iblk)%p, w_normal_global)
c
c.... end of the 1 side              
c... ends of getting weighted normal
          endif
c
          call e3if_geom_mfree
c
c... deallocation for weighted normal
          if (i_w_normal .eq. 1) then
            deallocate(calc_factor_temp)
          endif 
c
        enddo if_blocks1
c
        if (i_w_normal .eq. 1) then        
c...communication of the weighted normal
          if (numpe > 1) then
            call commu (w_normal_global, ilwork, nsd  , 'in ')
            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
            call commu (w_normal_global, ilwork, nsd  , 'out')
            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          endif
c...normalize the weighted normal
          do inode = 1, nshg
            if ( ifFlag(inode) .eq. 1 ) then
              length_temp(inode) = sqrt( w_normal_global(inode,1)
     &                           * w_normal_global(inode,1)
     &                           + w_normal_global(inode,2)
     &                           * w_normal_global(inode,2)
     &                           + w_normal_global(inode,3)
     &                           * w_normal_global(inode,3) )
              do isd = 1, nsd
                w_normal_global(inode,isd) = w_normal_global(inode,isd) 
     &                                     / length_temp(inode)
              enddo
            endif
          enddo
c... changing from inward normal to outward normal
          w_normal_global = - w_normal_global                      
        endif        
c
        if (numpe > 1) then
          call commu (if_normal(:,1:3), ilwork, nsd, 'in ')
          if (associated(if_kappa)) then
            call commu (if_kappa(:,1:nsd), ilwork, nsd, 'in ')
            call commu (if_kappa(:,nsd+1), ilwork, 1, 'in ')
          endif
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c
        do inode = 1,nshg
          length = sqrt(dot_product(if_normal(inode,1:nsd),if_normal(inode,1:nsd)))
          if (length > zero) then
            if_normal(inode,1:nsd) = if_normal(inode,1:nsd) / length
          endif
          if (associated(if_kappa)) then
            if (if_kappa(inode,nsd+1) > zero) 
C     &        if_kappa(inode,1:nsd) = if_kappa(inode,1:nsd)/if_kappa(inode,nsd+1)
     &        if_kappa(inode,1:nsd) = pt50*if_kappa(inode,1:nsd)/if_kappa(inode,nsd+1)
          endif
        enddo
c
        if (numpe > 1 .and. associated(if_kappa)) then
          call commu (if_kappa(:,1:nsd), ilwork, nsd, 'out')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c
c        call calc_kappa_error(x,lcblkif(1,:),nelblif,nsd,nshg)
c
        sum_vi_area = zero
        ifbc = zero    
c
        if_blocks: do iblk = 1, nelblif
c
c... set up the parameters
c
          iblkts  = iblk                ! used in time series
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          iel     = lcblkif(1, iblk)
          lelCat  = lcblkif(2, iblk)    ! ??? NOT USED?
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          iorder  = lcblkif(5, iblk)    ! polynomial order
          nshl0   = lcblkif(iblkif_nshl0,iblk)
          nshl1   = lcblkif(iblkif_nshl1,iblk)
          itpid   = lcblkif(iblkif_topology,iblk)
          mater0  = lcblkif(9, iblk)
          mater1  = lcblkif(10,iblk)
          ndof    = lcblkif(11,iblk)
          nsymdl  = lcblkif(12,iblk)    ! ???
          npro    = lcblkif(1,iblk+1) - iel
          inum    = iel + npro - 1
          ngaussif = nintif(itpid)
c
c... remember that the 0 side goes to the first material in the solver.inp
c...           and the 1 side goes to the second material in the solver.inp
c
            ienif0 => mienif0(iblk)%p
            ienif1 => mienif1(iblk)%p
c
c... set equations of state
c
          get_vap_frac0 => e3if_empty
c
          select case (mat_eos(mater0,1))
          case (ieos_ideal_gas)
            getthmif0_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthmif0_ptr => getthm7_ideal_gas_mixture
            get_vap_frac0 => get_vapor_fraction0
          case (ieos_liquid_1)
            getthmif0_ptr => getthm7_liquid_1
          case default
            call error ('getthm  ', 'WE DO NOT SUPPORT THIS MATERIAL (3)', mater1)
          end select
c
c          get_vap_frac1 => e3if_empty
c
          select case (mat_eos(mater1,1))
          case (ieos_ideal_gas)
            getthmif1_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
c            getthmif1_ptr => getthm7_ideal_gas_mixture
c            get_vap_frac1 => get_vapor_fraction1
            write(*,*) 'ERROR: Please put mixture on the 0 side (first in the input)'
            call error ('getthm  ', 'WE DO NOT SUPPORT THIS MATERIAL (3)', mater1)
          case (ieos_liquid_1)
            getthmif1_ptr => getthm7_liquid_1
          case default
            call error ('getthm  ', 'WE DO NOT SUPPORT THIS MATERIAL (3)', mater1)
          end select
c
c... compute and assemble the residual and tangent matrix
c
          if (lhs .eq. 1) then
            nedof0 = nflow*nshl0
            nedof1 = nflow*nshl1
            allocate (egmassif00(npro,nedof0,nedof0))
            allocate (egmassif01(npro,nedof0,nedof1))
            allocate (egmassif10(npro,nedof1,nedof0))
            allocate (egmassif11(npro,nedof1,nedof1))
            egmassif00 = zero
            egmassif01 = zero
            egmassif10 = zero
            egmassif11 = zero
          else
            allocate (egmassif00(1,1,1))
            allocate (egmassif01(1,1,1))
            allocate (egmassif10(1,1,1))
            allocate (egmassif11(1,1,1))
          endif

          call e3if_setparam
     &    (
     &     nshg,nshl0,nshl1,nenl0,nenl1,lcsyst0,lcsyst1,
     &     npro,ndof,nsd,nflow,ipord,ngaussif
     &    )
     &   
c
          call e3if_setparam2
     &    (
     &     egmassif00,egmassif01,egmassif10,egmassif11,
     &     real(lstep+1,8)*delt(1)
     &    )
c
          call e3if_geom_malloc
          call e3if_malloc
c... allocation of w_normal_l0 and w_normal_l1 for weighted normal
          if (i_w_normal .eq. 1) then
            allocate(w_normal_l0(npro,nshl0,nsd))
            allocate(w_normal_l1(npro,nshl1,nsd))
          endif          
c
c      do i=1,nshg
c        print*, i,res(i,:)
c      enddo
          call asidgif
     &    (
     &      res,
     &      y, x, umesh,
     &      shpif(lcsyst0,1:nshl0,:), 
     &      shpif(lcsyst1,1:nshl1,:), 
     &      shgif(lcsyst0,1:nsd,1:nshl0,:),
     &      shgif(lcsyst1,1:nsd,1:nshl1,:),
     &      qwtif(itpid,:), qwtif(lcsyst0,:), qwtif(lcsyst1,:),
     &      ienif0, ienif1,
     &      BDiag)
c
          if (lhs .eq. 1) then
c
c.... Fill-up the global sparse LHS mass matrix
c
            call fillsparse_if(lhsk,ienif0,ienif0,col,row,egmassif00,npro,nedof0,nedof0,nflow,nshg,nnz,nnz_tot)
            call fillsparse_if(lhsk,ienif1,ienif0,col,row,egmassif10,npro,nedof1,nedof0,nflow,nshg,nnz,nnz_tot)
            call fillsparse_if(lhsk,ienif0,ienif1,col,row,egmassif01,npro,nedof0,nedof1,nflow,nshg,nnz,nnz_tot)
            call fillsparse_if(lhsk,ienif1,ienif1,col,row,egmassif11,npro,nedof1,nedof1,nflow,nshg,nnz,nnz_tot)
c
          endif
c
          call e3if_geom_mfree
          call e3if_mfree
c
          deallocate (egmassif00)
          deallocate (egmassif01)
          deallocate (egmassif10)
          deallocate (egmassif11)
c... deallocation for weighted normal
          if (i_w_normal .eq. 1) then
            deallocate(w_normal_l0, w_normal_l1)
          endif          
c
        enddo if_blocks
c
#if debug==1
      write(*,'(a22,i6,5e24.16)') 'ELMGMR: after interface',idbg,res(idbg,:)
#endif
c
        if (associated(if_kappa)) then
          deallocate (if_kappa)
          nullify(if_kappa)
        endif
c... deallocation of the arrays used in weighted normal
        if (i_w_normal .eq. 1) then
          deallocate(w_normal_global)
          deallocate(length_temp)
        endif        
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
      if(iabc==1) then               !are there any axisym bc's
          call rotabc(res(1,2), iBC,  'in ')
c          Bdiagvec(:,1)=BDiag(:,1,1)
c          Bdiagvec(:,2)=BDiag(:,2,2)
c          Bdiagvec(:,3)=BDiag(:,3,3)
c          Bdiagvec(:,4)=BDiag(:,4,4)
c          Bdiagvec(:,5)=BDiag(:,5,5)
c          call rotabc(Bdiagvec(1,2), iBC,  'in ')
c          BDiag(:,:,:)=zero
c          BDiag(:,1,1)=Bdiagvec(:,1)
c          BDiag(:,2,2)=Bdiagvec(:,2)
c          BDiag(:,3,3)=Bdiagvec(:,3)
c          BDiag(:,4,4)=Bdiagvec(:,4)
c          BDiag(:,5,5)=Bdiagvec(:,5)
       endif

c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
        call commu_rbForce

        call commu (res  , ilwork, nflow  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        if(iprec .ne. 0) call commu (BDiag, ilwork, nflow*nflow, 'in ')
      endif
c
c------> BEGIN DEBUG <---------
c
c      write(*,998) '[',myrank,'] in elmgmr AFTER commu.'
c      do i = 1,nshg
c        if (abs(x(i,1)-x1)<tol .and. abs(x(i,2)-x2)<tol .and. abs(x(i,3)-x3)<tol) then
c        if (any(abs(res(i,:))>1.e-6)) then
c          write(*,997) myrank,i,x(i,:),y(i,:)
c          write(*,999) myrank,i,x(i,:),res(i,:)
c        end if
c      enddo
c
c      do irank = 0,numpe-1
c        call MPI_Barrier (MPI_COMM_WORLD,ierr)
c        if (irank == myrank) then
c          numtask = ilwork(1)
c          itkbeg = 1
c          m = 0
c!          write(*,990) myrank,numtask
c          do itask = 1, numtask
c            m = m + 1
c            iother = ilwork (itkbeg + 3)
c            numseg = ilwork (itkbeg + 4)
c!            write(*,991) myrank,ilwork(itkbeg+1:itkbeg+5)
c!      if (myrank == 0 .and. iother == 1 .or.
c!     &    myrank == 1 .and. iother == 0) then
c        do is = 1,numseg
c          isgbeg = ilwork(itkbeg + 3 + 2*is)
c          lenseg = ilwork(itkbeg + 4 + 2*is)
c          isgend = isgbeg + lenseg - 1
c          do isg = isgbeg,isgend
c!            write(*,801) myrank,is,isg,lenseg,x(isg,:)!,res(isg,:)
c          enddo
c        enddo
c!      endif
c            itkbeg = itkbeg + 4 + 2*numseg
c          enddo
c        endif
c      enddo
c
801   format('[',i2,'] is,isgbeg,lenseg,isg,x:',3i4,x,3f8.3,x,5e24.16)
990   format('[',i2,'] numtask:',i3)
991   format('[',i2,'] itag, iacc, iother, numseg, isgbeg:',5i6)
997   format('[',i2,'] i,x,y:  ',i6,3f7.3,5e24.16)
998   format(a,i2,a)
999   format('[',i2,'] i,x,res:',i6,3f7.3,5e24.16)
c
c--------> END DEBUG <--------------
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (y,  iBC,  BC,  res,  iper, ilwork)
c
c------> BEGIN DEBUG <-----------
c
!      write(*,998) '[',myrank,'] in elmgmr AFTER bc3res.'
!      do i = 1,nshg
!        if (x(i,1) < 0.1001 .and. x(i,1) > 0.0999 
!     & .and. x(i,2) < 0.0901 .and. x(i,2) > 0.0899
!     & .and. x(i,3) < 0.0101 .and. x(i,3) > 0.0099
!     &)
!     &   write(*,999) '[',myrank,'] :',i,x(i,:),res(i,:)
!      enddo
c------> END DEBUG <-----------
c
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
c  This code fragment would extract the "on processor diagonal
c      blocks". lhsK alread has the BC's applied to it (using BC3lhs), 
c      though it was on an ebe basis. For now, we forgo this and still 
c      form BDiag before BC3lhs, leaving the need to still apply BC's
c      below.  Keep in mind that if we used the code fragment below we
c      would still need to make BDiag periodic since BC3lhs did not do
c      that part.
c
      if (iprec .ne. 0) then
c$$$         do iaa=1,nshg
c$$$            k = sparseloc( row(col(iaa)), col(iaa+1)-colm(iaa), iaa )
c$$$     &       + colm(iaa)-1
c$$$            do idof=1,nflow
c$$$               do jdof=1,nflow
c$$$                  idx=idof+(jdof-1)*nflow
c$$$                  BDiag(iaa,idof,jdof)=lhsK(idx,k)
c$$$               enddo
c$$$            enddo
c$$$         enddo
         call bc3BDg (y,  iBC,  BC,  BDiag, iper, ilwork)
      endif
c
c.... return
c
      call timer ('Back    ')
      return
      end
c
c

c
        subroutine ElmGMRSclr(y,      ac,
     &                        x,      elDw,         
     &                        shp,    shgl,   iBC,
     &                        BC,     shpb,   shglb,
     &                        rest,   rmest,  Diag,
     &                        iper,   ilwork, EGmasst,
     &                        umesh )
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
        use eqn_state_m
        use e3Sclr_param_m
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),              
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            rest(nshg),           Diag(nshg),
     &            rmest(nshg),          BDiag(nshg,nflow,nflow),
     &            iper(nshg),           EGmasst(numel,nshape,nshape)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qrest(nshg),          rmasst(nshg)
c
        dimension ilwork(nlwork)
        dimension elDw(numel)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
        real*8, allocatable :: elDwl(:)
        real*8, dimension(numnp,nsd), intent(in) :: umesh
c
	ttim(80) = ttim(80) - tmr()
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        intrul = intg  (1,itseq)
        intind = intpt (intrul)
c
        ires   = 1

c       if (idiff>=1) then  ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
c        qrest = zero
c        rmasst = zero
c        
c        do iblk = 1, nelblk
c
c.... set up the parameters
c
c          iel    = lcblk(1,iblk)
c          lelCat = lcblk(2,iblk)
c          lcsyst = lcblk(3,iblk)
c          iorder = lcblk(4,iblk)
c          nenl   = lcblk(5,iblk)   ! no. of vertices per element
c          mattyp = lcblk(7,iblk)
c          ndofl  = lcblk(8,iblk)
c          nsymdl = lcblk(9,iblk)
c          npro   = lcblk(1,iblk+1) - iel 
c
c          nintg  = numQpt (nsd,intrul,lcsyst)
c          
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass
c
c          call AsIq (y,                x,                       
c     &               shp(1,intind,lelCat), 
c     &               shgl(1,intind,lelCat), 
c     &               mien(iblk)%p,    mxmudmi(iblk)%p,  
c     &               qres,           rmass )
c       
c       enddo
c       
c
c.... compute qi for node A, i.e., qres <-- qres/rmass
c
c       if (numpe > 1) then
c          call commu (qres  , ilwork, (ndof-1)*nsd  , 'in ')
c          call commu (rmass , ilwork,  1            , 'in ')
c       endif
c
c.... take care of periodic boundary conditions
c
c       call qpbc( rmass, qres, iBC, iper )       
c       
c       rmass = one/rmass
c       
c       do i=1, (nflow-1)*nsd
c          qres(:,i) = rmass*qres(:,i)
c       enddo
c
c       if(numpe > 1) then
c          call commu (qres, ilwork, (nflow-1)*nsd, 'out')    
c       endif
c
c      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        rest    = zero
        rmest   = zero ! to avoid trap_uninitialized
        if (lhs .eq. 1)   EGmasst = zero
        if (iprec. ne. 0)   Diag  = zero 
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mater  = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          allocate (elDwl(npro))
c
          select case (mat_eos(mater,1))
          case (ieos_ideal_gas,ieos_ideal_gas_2)
            getthm7_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthm7_ptr => getthm7_ideal_gas_mixture
          case (ieos_liquid_1)
            getthm7_ptr => getthm7_liquid_1
          case default
            call error ('getthm  ', 'wrong material', mater)
          end select
c
          call e3Sclr_malloc
c
          call AsIGMRSclr(y,                   
     &                    ac,
     &                    x,               elDwl,                   
     &                    tmpshp,          tmpshgl,
     &                    mien(iblk)%p,
     &                    mmat(iblk)%p,    rest,
     &                    rmest,               
     &                    qrest,           EGmasst(iel:inum,:,:),
     &                    Diag ,
     &                    umesh)
c
 
c.... satisfy the BC's on the implicit LHS
c     
          call bc3LHSSclr (iBC, mien(iblk)%p, EGmasst(iel:inum,:,:) )
c
          elDw(iel:inum)=elDwl(1:npro)
          deallocate ( elDwl )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
          call e3Sclr_mfree
c
c.... end of interior element loop
c
       enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... set up parameters
c
        intrul = intg   (2,itseq)
        intind = intptb (intrul)
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mater  = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)
c
          select case (mat_eos(mater ,1))
          case (ieos_ideal_gas,ieos_ideal_gas_2)
            getthm7_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthm7_ptr => getthm7_ideal_gas_mixture
          case (ieos_liquid_1)
            getthm7_ptr => getthm7_liquid_1
          case default
            call error ('getthm  ', 'wrong material', mater )
          end select
c
          call e3Sclr_malloc
c
          call AsBMFGSclr (y,                  x,
     &                     tmpshpb,
     &                     tmpshglb, 
     &                     mienb(iblk)%p,      mmatb(iblk)%p,
     &                     miBCB(iblk)%p,      mBCB(iblk)%p,
     &                     rest,               rmest)
c
          deallocate ( tmpshpb )
          deallocate ( tmpshglb )
c
          call e3Sclr_mfree

c.... end of boundary element loop
c
        enddo


      ttim(80) = ttim(80) + tmr()
c
c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (rest  , ilwork, 1  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        if(iprec .ne. 0) call commu (Diag, ilwork, 1, 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (y,  iBC,  BC,  rest,  iper, ilwork)
c
c.... satisfy the BCs on the preconditioner
c
      call bc3BDgSclr (iBC, Diag, iper, ilwork)
c
c.... return
c
      call timer ('Back    ')
      return
      end

