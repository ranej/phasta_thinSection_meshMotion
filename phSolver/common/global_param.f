c
c.... NOTE:  IF YOU CHANGE ANY OF THESE, YOU HAVE TO CHANGE THEM IN
c                  common_c.h ALSO
c
      module global_const_m
        use iso_c_binding
        implicit none
        integer, parameter     :: MAXBLK = 50000
        integer, parameter     :: MAXSH = 32, NSD = 3 , NSDSQ = 9
c
c  The five types of region topology are  1= Tet, 2=Hex, 3= Wedge (tri-start),
c                                         4= Wedge (quad-first) 5=pyramid
c
c  The two types of face topology are  1= tri, 2=quad
c
        integer, parameter     :: MAXTOP = 6, MAXSURF=1000
c
        integer, parameter     :: itp_tet          = 1
     &,                           itp_hex          = 2
     &,                           itp_wedge        = 3
     &,                           itp_wedge_tri    = 3
     &,                           itp_wedge_quad   = 4
     &,                           itp_pyramid      = 5
     &,                           itp_pyramid_quad = 5
     &,                           itp_pyramid_tri  = 6
     &,                           itp_tri     = 1
     &,                           itp_quad    = 2
c
c...  The twelve different topological interface region are:
c
        integer, parameter     :: MAXTOPIF = 4
c
c  sharing a tri face:
c
c  1= tet-tet 
c  2= tet-wedge
c  3= wedge-tet
c  4= wedge-wedge
c
        integer, parameter     :: itpif_tet_tet     = 1
     &,                           itpif_tet_wedge   = 2
     &,                           itpif_wedge_tet   = 3
     &,                           itpif_wedge_wedge = 4
c
        integer, dimension(MAXTOP,MAXTOP), parameter :: iftp_map = 
     &                reshape( (/itpif_tet_tet,itpif_tet_wedge,-1,-1,-1,-1,
     &                  -1,-1,-1,-1,-1,-1,
     &                  itpif_wedge_tet,-1,itpif_wedge_wedge,-1,-1,-1,
     &                  -1,-1,-1,-1,-1,-1,
     &                  -1,-1,-1,-1,-1,-1,
     &                  -1,-1,-1,-1,-1,-1/) , (/MAXTOP,MAXTOP/) )
c
        integer, parameter :: MAXQPT = 125
        integer, parameter :: MAXTS = 100
c
        integer, parameter :: MAXMAT  = 6, MAXPROP = 10
c
        real*8, parameter :: Ru = 8.314d0    ! Universal gas constant [J/mol.K]
c
      end module global_const_m
c
c----------------------------------------------------------------------
c
c.... common /conpar/   : input constants
c
c numnp         : number of nodal points
c numel         : number of elements
c numelb        : number of boundary elements
c numpbc        : number of nodes having a boundary condition
c nen           : maximum number of element nodes
c nfaces        : maximum number of element faces
c nsd           : number of space dimensions
c numflx        : number of flux boundary nodes
c ndof          : number of degrees of freedom per node
c iALE          : ALE formulation flag
c iSOLID        : Solid formulation flag
c icoord        : coordinate system flag
c navier        : Navier-Stokes calculation flag
c irs           : restart option 
c iexec         : execute flag
c necho         : input echo parameter
c ichem         : equilibrium chemistry flag (for outchem.step dump)
c iRK           : Runge-Kutta flag
c nshg          : global number of shape functions (degrees of freedom,
c                 or equations). Computed from the specified p-order,
c                 the number of edges, and the number of faces (in the
c                 entire mesh)
c
c----------------------------------------------------------------------
c
      module conpar_m
        use iso_c_binding
        implicit none
        integer(c_int), target :: numnp, numel,  numelb, numelif,
     &                  numpbc,   nen,    nfaces,
     &                  numflx,   ndof,   iALE, iSOLID,
     &                  icoord,   navier,
     &                  irs,      iexec,  necho,  ichem,  iRK,    nedof,
     &                  ndofelas, nshg,   nnz,    istop,  nflow,  nelas,
     &                  nnz_tot,  idtn,
     &                  ncorpsize, iownnodes, usingpetsc,
     &                  elasModel, elasFDC, elasSICC,  mesh2geom
        common /conpar/ numnp, numel,  numelb, numelif,
     &                  numpbc,   nen,    nfaces,
     &                  numflx,   ndof,   iALE, iSOLID,
     &                  icoord,   navier,
     &                  irs,      iexec,  necho,  ichem,  iRK,    nedof,
     &                  ndofelas, nshg,   nnz,    istop,  nflow,  nelas,
     &                  nnz_tot,  idtn,
     &                  ncorpsize, iownnodes, usingpetsc,
     &                  elasModel, elasFDC, elasSICC, mesh2geom
      end module conpar_m
c
c----------------------------------------------------------------------
c
c.... common /laymesh/   : layered mesh
c
c numgc           : number of growth curves
c numgcnp         : total number of nodal points of all growth curves
c gcBaseOpt       : input base face option for repositioning method
c layerCommuFlag  : the flag used for layer mesh master assembly in commu
c blfactor        : multiplied with stiffness of wedge element in elas solver
c
      module laymesh_m
        use iso_c_binding
        implicit none
        real(c_double), target :: blfactor
        integer(c_int), target :: numgc,   numgcnp,   gcBaseOpt
        integer :: layerCommuFlag = 0
        common /laymesh/ blfactor,  numgc,   numgcnp,   gcBaseOpt,
     &                   layerCommuFlag
      end module laymesh_m
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c.... common /laymeshts/   : thin section mesh
c
c numts           : number of thin section stacks
c numtsnp         : total number of nodal points of all thin section stacks
c TSBaseOpt       : input base face option for repositioning method
c layerTSCommuFlag  : the flag used for layer mesh master assembly in commu
c TSfactor        : multiplied with stiffness of wedge element in elas solver
c
      module laymeshts_m
        use iso_c_binding
        implicit none
        real(c_double), target :: tsfactor
        integer(c_int), target :: numts,   numtsnp,   tsBaseOpt
        integer :: layertsCommuFlag = 0
        common /laymeshts/ tsfactor,  numts,   numtsnp,   tsBaseOpt,
     &                   layertsCommuFlag
      end module laymeshts_m
c
c----------------------------------------------------------------------
c
c.... common /snapmesh/   : snap mesh nodes back the surface
c
c snapSurfFlag    : flag for snap back method
c snapNumFaceTags : number of face tags; default is zero
c snapFaceTags    : face tags for the target surface
c
      module snapmesh_m
        use iso_c_binding
        use global_const_m
        implicit none
        integer(c_int), target, dimension(MAXTS) :: snapFaceTags
        integer(c_int), target :: snapSurfFlag,     snapNumFaceTags,
     &                            timeDepComp1Flag, timeDepComp1ID
        real*8, allocatable    :: timeDepComp1Mag(:)
        common /snapmesh/ snapSurfFlag,  snapNumFaceTags,  snapFaceTags,
     &                    timeDepComp1Flag, timeDepComp1ID
      end module snapmesh_m
c
c----------------------------------------------------------------------
c
c.... common /meshquality/   : mesh quality and auto adaptation trigger
c
c autoTrigger  : flag used to turn on/off auto mesh adaptation trigger option
c volMeshqTol  : threshold for volume mesh quality
c faceMeshqTol : threshold for mesh quality of triangle face
c                in boundary layered mesh
c triggerNow   : flag used to terminate solver at this time step
c
      module meshquality_m
        use iso_c_binding
        implicit none
        real(c_double), target :: volMeshqTol,  faceMeshqTol,
     &                            errorTolMass, errorTolMomt,
     &                            errorTolEngy, errorMaxMass,
     &                            errorMaxMomt, errorMaxEngy
        integer(c_int), target :: autoTrigger,     triggerNow, 
     &                            errorEstimation, errorTriggerEqn
        common /meshquality/  volMeshqTol,  faceMeshqTol,
     &                        autoTrigger,  triggerNow,
     &                        errorTolMass, errorTolMomt,
     &                        errorTolEngy, errorMaxMass,
     &                        errorMaxMomt, errorMaxEngy,
     &                        errorEstimation,
     &                        errorTriggerEqn
      end module meshquality_m
c
c----------------------------------------------------------------------
c
c.... common /rigidbody/   : rigid body parameters, properties, constraints
c
c numrbs         : number of rigid bodies
c rbsTags        : tags of rigid bodies
c rbsMM          : rigid body motion mode for each rigid body
c rb_prop        : properties and constraints of rigid bodies
c rb_commuMotion : communicate motion option flag
c
      module rigidbody_m
        use iso_c_binding
        use global_const_m
        implicit none
        real(c_double), target :: rb_prop(MAXTS, MAXTS)
        integer(c_int), target :: rbsTags(MAXTS)
        integer(c_int), target :: rbsMM(MAXTS)
        integer(c_int)         :: numrbs
        integer(c_int)         :: rb_commuMotion
        integer                :: rbParamSize = 14
        integer(c_int), target :: rbsPcase(MAXTS)
        common /rigidbody/     rb_prop,   numrbs,   rbsTags,
     &                         rbsMM,     rb_commuMotion, rbsPcase
      end module rigidbody_m
c
c----------------------------------------------------------------------
c
c.... common /rotatingband/   : rotating band parameters
c
c numRotBands         : number of rotating bands
c numRotBandFaceTags  : number of faces for each rotating band
c rotBandTag          : tags of faces for each rotating band
c
      module rotatingband_m
        use iso_c_binding
        use global_const_m
        implicit none
        integer(c_int)         :: numRotBands
        integer(c_int)         :: numRotBandFaceTags
        integer(c_int), target :: rotBandTag(MAXTS, MAXTS)
        integer(c_int), target :: rotBandFO(MAXTS)
        common /rotatingband/    numRotBands,   numRotBandFaceTags,  
     &                         rotBandTag, rotBandFO
      end module rotatingband_m
c
c----------------------------------------------------------------------
c.... common /timdat/   : time data
c
c time          : current run time
c CFLfld        : CFL number for fluid flow
c CFLsld        : CFL number for structural heating
c Dtgl          : inverse of global time step
c Dtmax         : maximum delta-time
c alpha         : trapezoidal rule parameter
c etol          : epsilon tolerance for GMRES
c lstep         : current time step
c ifunc         : func. eval. counter (=niter*(lstep-lstep0) + iter)
c itseq         : sequence number
c istep         : step number (reseted at the beginning of the run)
c iter          : current iteration number
c nitr          : number of multi-corrector iterations of flow solve
c                 for the current stagger
c
      module timdat_m
        use iso_c_binding
        implicit none
        integer(c_int) :: lstep, ifunc, itseq, istep, iter, nitr, iCFLworst, lskeep
        real(c_double) :: time, CFLfld, CFLsld, Dtgl, Dtmax, alpha, etol,
     &            almi, alfi, gami, almBi, alfBi, gamBi, flmpl, flmpr, dtol(2)
        common /timdat/ time,   CFLfld, CFLsld, Dtgl,   Dtmax,  alpha,
     &                  etol,   lstep,  ifunc,  itseq,  istep,  iter,
     &                  nitr,   almi,   alfi,   gami,   
     &                  almBi,  alfBi,  gamBi,
     &                  flmpl,  flmpr,
     &                  dtol, iCFLworst, lskeep
      end module timdat_m
c
c----------------------------------------------------------------------
c
c.... common /elmpar/   : element parameters
c
c lelCat        : element category (P1, Q1, P2, Q2, etc.)
c lcsyst        : element coordinate system
c iorder        : element order (=k for Pk and Qk)
c nenb          : number of element nodes per boundary sides
c maxsh         : total number integration points
c maxshb        : total number integration points of boundary elements
c nelblk        : number of element blocks
c nelblb        : number of boundary element blocks
c nelblif       : number of interface element blocks
c ndofl         : number of degrees of freedom (for current block)
c nsymdl        : number of d.o.f for symm. storage (for current block)
c nenl          : number of element nodes (for current block)
c nfacel        : number of element faces (for current block)
c nenbl         : number of boundary element nodes
c intind        : integration data index
c nintg         : number of integration points
c mattyp        : material type ( = 0 for fluid; = 1 for solid )
c
      module elmpar_m
        use iso_c_binding
        use global_const_m
        implicit none
        integer, target ::  nelblk, nelblb, nelblif
        integer :: lelCat, lcsyst, iorder, nenb,   
     &                  ndofl,  nsymdl, nenl,   nfacel,
     &                  nenl0,  nenl1,  lcsyst0, lcsyst1,
     &                  nenbl,  intind, mattyp,
     &                  mattyp0, mattyp1
     &,                 itpid                              ! element topology id
        common /elmpar/ lelCat, lcsyst, iorder, nenb,   
     &                  nelblk, nelblb, nelblif,
     &                  ndofl,  nsymdl, nenl,   nfacel,
     &                  nenbl,  intind, mattyp
      end module elmpar_m
c
c----------------------------------------------------------------------
c
c.... common /blkdat/   : blocking data
c
c lcblk  (10,MAXBLK+1) : blocking data for the interior elements
c lcblkb (10,MAXBLK+1) : blocking data for the boundary elements
c lcblkif (14,MAXBLK+1) : blocking data for the interface elements
c
      module blkdat_m
        use global_const_m
        implicit none
        integer, parameter :: iblk_mattype  = 7
     &,                       iblkif_nshl0    = 13
     &,                       iblkif_nshl1    = 14
     &,                       iblkif_topology = 15
        integer :: lcblk  (10,MAXBLK+1),
     &             lcblkb (10,MAXBLK+1),
     &             lcblkif(15,MAXBLK+1)
      end module blkdat_m
c
c----------------------------------------------------------------------
c
c.... common /intdat/   : integration data
c
c intg  (2,MAXTS) : integration parameters
c intpt (3)       : integration pointers
c intptb(3)       : integration pointers of boundary elements
c
      module intdat_m
        use iso_c_binding
        use global_const_m
        implicit none
        integer(c_int) :: intg(3,MAXTS)
        integer :: intpt(3), intptb(3)
        common /intdat/ intg,intpt,intptb
      end module intdat_m
c
c----------------------------------------------------------------------
      module intpt_m
        use global_const_m
        implicit none
c
c.... hierarchic basis functions
c
        real*8 :: Qpt (MAXTOP ,4,MAXQPT), Qwt (MAXTOP ,MAXQPT), 
     &            Qptb(MAXTOP,4,MAXQPT),  Qwtb(MAXTOP,MAXQPT), 
     &            Qptif (MAXTOP,4,MAXQPT), Qwtif (MAXTOP,MAXQPT),
     &            Qptif0(MAXTOP,4,MAXQPT), Qwtif0(MAXTOP,MAXQPT),
     &            Qptif1(MAXTOP,4,MAXQPT), Qwtif1(MAXTOP,MAXQPT)
        integer ::    nint(MAXTOP),           nintb(MAXTOP),
     &                nintif(MAXTOPIF),
     &                ngauss,                 ngaussb,        ngaussif,
     &                intp,
     &                maxnint
      end module intpt_m
c
      module shpdat_m
        use iso_c_binding
        implicit none
	      integer, target ::  ntopsh, nfath
        integer :: nshape, nshapeb, nshapeif, maxshb,
     &             nshl, nshlb, nshl0, nshl1,      ! nshl0,1 for interface
     &             nsonmax
        common /shpdat/ nshape, nshapeb, nshapeif, maxshb,
     &                  nshl, nshlb, nshl0, nshl1,      ! nshl0,1 for interface
     &                  nfath,  ntopsh,  nsonmax
      end module shpdat_m  
c
c----------------------------------------------------------------------
c
c.... common /genpar/   : control parameters
c
c E3nsd         : NSD .eq. 3 flag; 0. for 2D, 1. for 3D
c I3nsd         : NSD .eq. 3 flag; 0  for 2D, 1  for 3D
c nsymdf        : number of d.o.f.s in symm. storage (= ndof*(ndof+1)/2)
c ndofBC        : dimension size of the boundary condition array BC
c ndiBCB        : dimension size of the boundary condition array iBCB
c ndBCB         : dimension size of the boundary condition array BCB
c Jactyp        : Jacobian type flag
c jump          : jump term computation flag
c ires          : residual type computation flag
c iprec         : block-diagonal preconditioner flag
c iprev         : ypl array allocation flag
c ibound        : boundary element flag
c idiff         : diffusive flux vector flag
c                 ( = 0 not used; = 1 global reconstruction )
c itau          : type of tau to be used
c iLHScond      : add contributiosn from the heat flux BC to the LHS 
c                 tangency matrix. 
c ndofBC2       : dimension size of the boundary condition array BC before constraint
c
      module genpar_m
        use iso_c_binding
        implicit  none
        integer :: I3nsd, nsymdf, ndofBC, ndiBCB, ndBCB, Jactyp, jump, ires, iprec,
     &    iprev, ibound, idiff, lhs, itau, ipord, ipred, lstres, iepstm, ibksiz, 
     &    iabc, isurf, idflx, 
     &    EntropyPressure, irampViscOutlet, istretchOutlet, iremoveStabTimeTerm, 
     &    iLHScond, ndofBC2
        real*8 :: E3nsd, dtsfct, taucfct, Bo
        common /genpar/ E3nsd,  I3nsd,  nsymdf, ndofBC, ndiBCB, ndBCB,
     &                  Jactyp, jump,   ires,   iprec,  iprev,  ibound,
     &                  idiff,  lhs,    itau,   ipord,  ipred,  lstres,
     &                  iepstm, dtsfct, taucfct, ibksiz, iabc, isurf,
     &                  idflx,  Bo,     EntropyPressure, irampViscOutlet,
     &                  istretchOutlet, iremoveStabTimeTerm, iLHScond,
     &                  ndofBC2
      end module genpar_m
c
c----------------------------------------------------------------------
c
c.... common /inpdat/   : time sequence input data
c
c epstol (MAXTS)  : tolerance for GMRES solvers
c etolelas        : tolerance for Mesh Elas solvers
c Delt   (MAXTS)  : global time step
c CFLfl  (MAXTS)  : CFL number for fluid flow
c CFLsl  (MAXTS)  : CFL number for structural heating
c nstep  (MAXTS)  : number of time steps
c niter  (MAXTS)  : number of iterations per time step
c impl   (MAXTS)  : solver flag
c iturb  (MAXTS)  : turbulence model flag
c rhoinf (MAXTS)  : time integration spectral radius paramter
c                             (0=Gears       1= trapezoidal rule)
c LHSupd (MAXTS)  : LHS/preconditioner update
c loctim (MAXTS)  : local time stepping flag
c
      module inpdat_m
        use iso_c_binding
        use global_const_m
        implicit none
        integer :: leslib, svLSFlag, svLSType
        integer, dimension(MAXTS) :: nstep, niter, impl, loctim
        integer, dimension(6) :: LHSupd
        real*8, dimension(6) :: epstol
        real*8  :: etolelas
        real*8, dimension(MAXTS) :: Delt, CFLfl, CFLsl, rhoinf, rhoinfS,
     &                              rhoinf_B,    rhoinf_rb
        real*8, dimension(MAXTS,2) :: deltol
        common /inpdat/ epstol,  etolelas, Delt,    CFLfl,
     &                  CFLsl,   nstep,    niter,
     &                  impl,    rhoinf,   rhoinfS,
     &                  rhoinf_B,          rhoinf_rb,
     &                  LHSupd,  loctim,  deltol, 
     &                  leslib,     svLSFlag,   svLSType
      end module inpdat_m
c
c----------------------------------------------------------------------
c
c....common /propar/    : processor related information
c
c npro          : number of virtual processors for the current block
c
      module propar_m
        use iso_c_binding
        implicit none
        integer(c_int) :: npro
        common /propar/ npro
      end module propar_m
c
      module mmatpar_m
        use iso_c_binding
        implicit none
        real(c_double) :: Rgas, gamma, gamma1, mw
        real(c_double) :: pr,     Planck, Stefan, Nh,     Rh,
     &                    s0, const,  xN2,    xO2,
     &                    yN2,    yO2,    Msh(5), cpsh(5),s0sh(5),h0sh(5),
     &                    Rs(5),  cps(5), cvs(5), h0s(5), Trot(5),sigs(5),
     &                    Tvib(5),g0s(5), dofs(5),ithm
        common /mmatpar/ Rgas, gamma, gamma1, mw, pr
      end module mmatpar_m
c
      module mtimer1_m
        use iso_c_binding
        character(len=8) :: ccode(13)
        common /mtimer1/ ccode
      end module mtimer1_m
c
      module mtimer2_m
        use iso_c_binding
        implicit none
        integer(c_int) :: iclock, icd, icode, icode2, icode3
        real(c_double) :: flops,gbytes,sbytes
        common /mtimer2/ flops,gbytes,sbytes,iclock,icd,icode,icode2,icode3
      end module mtimer2_m
c
      module timer3_m
        use iso_c_binding
        implicit none
        integer(c_int) :: nacess(11)
        real(c_double) :: cpu(11), cpu0(11)
        common /timer3/ cpu,cpu0,nacess
      end module timer3_m
c
      module dgifinp_m
        use iso_c_binding
        implicit none
        integer, parameter :: no_ramp = 1, linear_ramp = 2
        integer, parameter :: no_vi = 1, 
     &                        const_vi = 2, 
     &                        vieilles_burning=3, 
     &                        clausius_clapeyron=4,
     &                        cavitation=5
        integer(c_int) :: phase_change_model, vi_ramping, i_w_normal,
     &                    i_if_dc
        real(c_double) :: ramp_time, vi_mag, dgif_alpha, dgif_beta, dgif_s, dgif_e, dgif_h,
     &                    dgif_emu, dgif_ek, if_e_dc, if_reaction_heat  
        real(c_double) :: burn_rate_exp, burn_rate_coeff, burn_rate_pref
        real(c_double) :: hfg_liquid, mw_liquid, T_boil_liquid
        common /dgifinp/ phase_change_model,vi_ramping,i_w_normal,i_if_dc,
     &                   ramp_time,vi_mag,dgif_s,dgif_e,dgif_emu,dgif_ek,
     &                   dgif_h, if_e_dc, if_reaction_heat,
     &                   hfg_liquid, mw_liquid, T_boil_liquid,
     &                   burn_rate_exp, burn_rate_coeff, burn_rate_pref
      end module dgifinp_m
c
c----------------------------------------------------------------------
c
c.... common /outpar/   : output parameters
c
c ro            : density     rescaling factor for output
c vel           : velocity    rescaling factor for output
c temper        : temperature rescaling factor for output
c press         : pressure    rescaling factor for output
c entrop        : entropy     rescaling factor for output
c ntout         : number of steps between consecutive printouts
c ioform        : output I/O format
c
      module outpar_m
        use iso_c_binding
        implicit none
        integer(c_int) :: ntout,ioform,iowflux,iofieldv,ioybar,nstepsincycle,nphasesincycle,
     &    ncycles_startphaseavg,ivort,icomputevort,nsynciofiles,nsynciofieldswriterestart,
     &    iv_rankpercore,iv_corepernode,input_mode,output_mode,conservation_probe,
     &    write_residual, imeshCFL
        real(c_double) :: ro,vel,temper,press,entrop
        character(len=80) :: iotype
        common /outpar/ ro,     vel,    temper, press,  entrop, ntout,
     &                  ioform, iowflux, iofieldv, iotype, ioybar,
     &                  nstepsincycle, nphasesincycle, 
     &                  ncycles_startphaseavg, ivort, icomputevort,
     &                  nsynciofiles, nsynciofieldswriterestart, 
     &                  iv_rankpercore, iv_corepernode, 
     &                  input_mode, output_mode, conservation_probe,
     &                  write_residual, imeshCFL
      end module outpar_m
c
      module workfc_m
        use iso_c_binding
        implicit none
        integer(c_int) :: master, numpe, myrank
	      common /workfc/ master, numpe, myrank
      end module workfc_m
c
c----------------------------------------------------------------------
c
c.... common /io    /   : io channels
c
c iin           : input  (main parameters)          [INPUT.DAT]
c igeom         : input  (problem geometry)         [GEOM.DAT]
c ipar          : in/out (spectral mapping)         [PARTITION.DAT]
c ibndc         : input  (problem boundary cond.)   [BC.DAT]
c imat          : input  (element material types)   [MATERIAL.DAT]
c iecho         : output (echo of input)            [ECHO.DAT]
c iout          : output (result output)            [OUTPUT.lstep]
c ichmou        : output (chemistry output)         [OUTCHM.lstep]
c irstin        : input  (input restart)            [RESTAR.INP]
c irstou        : output (output restart)           [RESTAR.OUT]
c ihist         : output (history output)           [HISTOR.DAT]
c iflux         : output (boundary flux)            [FLUX.lstep]
c ierror        : output (error messages)           [ERROR.DAT]
c itable        : input  (equilibrium chemistry)    [TABLE.DAT]
c iforce        : output (aerodynamic forces)       [FORCES.DAT]
c
      module mio_m
        use iso_c_binding
        implicit none
        integer(c_int) :: iin,    igeom,  ipar,   ibndc,  imat,   iecho,
     &                  iout,   ichmou, irstin, irstou, ihist,  iflux,
     &                  ierror, itable, iforce, igraph, itime, iconserv
        common /mio    / iin,    igeom,  ipar,   ibndc,  imat,   iecho,
     &                  iout,   ichmou, irstin, irstou, ihist,  iflux,
     &                  ierror, itable, iforce, igraph, itime, iconserv
      end module mio_m
c
c----------------------------------------------------------------------
c
c.... common /ioname/   : io file names
c
c fin           : input.dat
c fgeom         : geom.dat
c fpar          : partition.dat
c fbndc         : bc.dat
c fmat          : material.dat
c fecho         : echo.dat
c frstin        : restar.inp
c frstou        : restar.out
c fhist         : histor.dat
c ferror        : error.dat
c ftable        : table.dat
c fforce        : forces.dat
c
      module mioname_m
        use iso_c_binding
        implicit none
        character*80    fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &                  frstin, frstou, fhist,  ferror, ftable, fforce,
     &                  fgraph, ftime, fconserv
        common /mioname/ fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &                  frstin, frstou, fhist,  ferror, ftable, fforce,
     &                  fgraph, ftime, fconserv
      end module mioname_m
c
      module sclrs_m
        use iso_c_binding
        implicit none
        real(c_double) :: scdiff(5), tdecay 
        integer(c_int) :: nsclr,isclr,nsolt
     &,                   nosource,consrv_sclr_conv_vel
        common /sclrs/ scdiff,tdecay,nsclr,isclr,nsolt,nosource,
     &            consrv_sclr_conv_vel
      end module sclrs_m
