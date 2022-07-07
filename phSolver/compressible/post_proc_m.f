        module post_param_m
          integer :: post_proc_loop
          real*8, allocatable :: meshCFL(:)
          real*8, allocatable :: VMS_error(:,:)
          real*8, allocatable :: meshCFLblk(:)
          real*8, allocatable :: VMS_errorblk(:,:)
        end module
c
        subroutine malloc_post_param
c
          use post_param_m
          use conpar_m
c
          allocate( meshCFL(numel) )
          allocate( VMS_error(numel, nflow) )
          post_proc_loop = 0
          meshCFL = zero
          VMS_error = zero
        end subroutine malloc_post_param
c
        subroutine release_post_param
c
          use post_param_m
c
          if (allocated(meshCFL))
     &      deallocate( meshCFL )
          if (allocated(VMS_error))
     &      deallocate( VMS_error )
c
        end subroutine release_post_param
c
        subroutine ElmPost (y,         ac,        x,
     &                      shp,       shgl,      iBC,
     &                      BC,        shpb,      shglb,
     &                      shpif,     shgif,
     &                      res,       iper,      ilwork,
     &                      rerr,      umesh)
c
c----------------------------------------------------------------------
c
c This routine loop over element blocks to collect necessary
c post-processing variables.
c
c Main part is copied from ElmGMRs.
c
c Fan Yang. March 2019.
c
c----------------------------------------------------------------------
c
        use pointer_data
        use e3_param_m
        use timedataC
        use eqn_state_m
        use e3_solid_m
        use probe_m
        use post_param_m
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),
     &            iBC(nshg),
     &            BC(nshg,ndofBC),
     &            res(nshg,nflow),
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
        real*8  rerr(nshg,10)
c
        real*8  tempMomtError
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
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
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
c
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
c.... initialize the arrays
c
        res    = zero
        meshCFL = zero
        VMS_error = zero
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
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          allocate (meshCFLblk(npro))
          allocate (VMS_errorblk(npro,nflow))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
          meshCFLblk = zero
          VMS_errorblk = zero
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
          call AsIPost (y,                   ac,
     &                  x,                   mxmudmi(iblk)%p,
     &                  tmpshp,
     &                  tmpshgl,             mien(iblk)%p,
     &                  mater,               res,
     &                  qres,
     &                  rerr,                umesh)

c.... map local meshCFL to global
          do i = 1, npro
            meshCFL(mieMap(iblk)%p(i)) = meshCFLblk(i)
          enddo
c
          if (errorEstimation .ge. 1) then
            do i = 1, npro
c.... map local VMS_error to global
              do j = 1, nflow
                VMS_error(mieMap(iblk)%p(i),j) = VMS_errorblk(i,j)
              enddo
c.... record the max error
              if (VMS_errorblk(i,1) .gt. errorMaxMass)
     &            errorMaxMass = VMS_errorblk(i,1)
              tempMomtError = sqrt(
     &                           VMS_errorblk(i,2)*VMS_errorblk(i,2)+
     &                           VMS_errorblk(i,3)*VMS_errorblk(i,3)+
     &                           VMS_errorblk(i,4)*VMS_errorblk(i,4) )
              if (tempMomtError .gt. errorMaxMomt)
     &            errorMaxMomt = tempMomtError
              if (VMS_errorblk(i,5) .gt. errorMaxEngy)
     &            errorMaxEngy = VMS_errorblk(i,5)
            enddo
          endif
c
          if (associated(e3_mfree_ptr)) call e3_mfree_ptr
c
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
          deallocate ( meshCFLblk )
          deallocate ( VMS_errorblk )
c
c.... end of interior element loop
c
       enddo
c
      return
      end
c
c
        subroutine AsIPost (y,       ac,      x,       xmudmi,
     &                      shp,     shgl,    ien,
     &                      mater,   res,     qres,
     &                      rerr,    umesh)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the necessary variables
c for post-processing
c
c Main part is copied from AsIGMRs.
c
c----------------------------------------------------------------------
c
        use rlssave     ! Use the resolved Leonard stresses at the nodes.
        use timedataC    ! time series
        use specialBC    ! get ytarget to localize and send down
        use post_param_m
        include "common.h"
c
        dimension y(nshg,ndofl),            ac(nshg,ndofl),
     &            x(numnp,nsd),
     &            shp(nshl,MAXQPT),
     &            shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),
     &            res(nshg,nflow),
     &            qres(nshg,idflx)
      integer, intent(in) :: mater
c
        dimension ycl(npro,nshl,ndofl),     acl(npro,nshl,ndof),
     &            xl(npro,nenl,nsd),        ytargetl(npro,nshl,nflow),
     &            rl(npro,nshl,nflow),      rml(npro,nshl,nflow),
     &            BDiagl(npro,nshl,nflow,nflow),
     &            ql(npro,nshl,idflx)
c
        dimension xmudmi(npro,ngauss)
        dimension sgn(npro,nshl),  EGmass(npro,nedof,nedof)
c
        dimension umesh(numnp, nsd),  uml(npro,nshl,nsd)
c
        dimension rlsl(npro,nshl,6)
        real*8 rerrl(npro,nshl,6), rerr(nshg,10)
c
c
        EGmass = zero
c
c.... create the matrix of mode signs for the hierarchic basis
c     functions.
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
c
        if(matflg(5,1).ge.4 )
     &   call localy (ytarget,   ytargetl,  ien,   nflow,  'gather  ')
c
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
c
c.... turn on post-processing flag
        post_proc_loop = 1
c
        call e3  (ycl,     ycl,     acl,     shp,
     &            shgl,    xl,      rl,      rml,   xmudmi,
     &            BDiagl,  ql,      sgn,     rlsl,  EGmass,
     &            rerrl,   ytargetl, uml)
c
c.... turn off post-processing flag
        post_proc_loop = 0
c
c.... assemble the residual and modified residual
c
        call local (res,    rl,     ien,    nflow,  'scatter ')
c
        if ( ierrcalc .eq. 1 ) then
           call local (rerr, rerrl,  ien, 6, 'scatter ')
        endif
c
c.... end
c
        return
        end
c
