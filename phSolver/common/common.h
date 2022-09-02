c----------------------------------------------------------------------
c
c This file contains the common blocks and the data declaration needed
c for the routines.
c
c Input variables that have been previously declared in common_c.h have to be
c re-declared here, in a consistant block. 
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use number_def_m
      use matdat_def_m
      use dgifinp_m
      use global_const_m
      use conpar_m
      use laymesh_m
      use laymeshts_m
      use snapmesh_m
      use meshquality_m
      use rigidbody_m
      use rotatingband_m
      use timedependbcflow_m
      use timdat_m
      use elmpar_m
      use blkdat_m
      use intpt_m
      use intdat_m
      use shpdat_m
      use genpar_m
      use inpdat_m
      use propar_m
      use mmatpar_m
      use mtimer1_m
      use mtimer2_m
      use timer3_m
      use outpar_m
      use workfc_m
      use mio_m
      use mioname_m
      use sclrs_m

	IMPLICIT REAL*8 (a-h,o-z)
c
c the common block nomodule holds all the things which have been removed
c from different modules
     
        integer seqsize, stepseq
        integer spongecontinuity, spongemomentum1, spongemomentum2
        integer spongeenergy, spongemomentum3
        integer*8 nshgt,minowned, maxowned
	common /amgvarr/strong_eps,ramg_eps,ramg_relax,ramg_trunc,
     &              ramg_chebyratio
	common /amgvari/irun_amg,irun_amg_prec,
     &                  iamg_verb,
     &                  iamg_neg_sten,iamg_nlevel,
     &                  iamg_c_solver,
     &                  iamg_init,
     &        iamg_setup_frez,
     &        iamg_interp,maxnev,maxncv,iamg_smoother,mlsdeg,
     &        iamg_reduce

        common /nomodule/ bcttimescale,ValueListResist(0:MAXSURF),
     &  rhovw,thicknessvw, evw, rnuvw, rshearconstantvw, betai,
     &  icardio, itvn, ipvsq, numResistSrfs, nsrflistResist(0:MAXSURF),
     &  numImpSrfs, nsrflistImp(0:MAXSURF),impfile,
     &  numRCRSrfs, nsrflistRCR(0:MAXSURF),ircrfile,
     &  ideformwall, iwallmassfactor, iwallstiffactor, iviscflux 
        common /sequence/ seqsize, stepseq(100)
	common /fronts/ maxfront, nlwork
	common /newdim/ nshgt, minowned,maxowned, numper, nshg0
	common /timer4/ birth, death, comtim
        common /extrat/ ttim(100)
        common /spongevar/ zoutSponge, radSponge, zinSponge,
     &            grthOSponge,grthISponge,betamax,
     &            spongecontinuity, spongemomentum1, spongemomentum2,
     &            spongeenergy, spongemomentum3
        common /turbvar/ eles,ylimit(3,9), rampmdot(2,3),
     &                   rmutarget, pzero,  wtavei, 
     &                   dtavei, dke,  fwr1, flump, DES_SA_hmin,
     &                   ierrcalc, ihessian, itwmod, ngaussf,idim,
     &                   nlist, nintf(MAXTOP)
        common /turbvari/iRANS, iLES, idistcalc, isubmod, ifproj,
     &                   i2filt, modlstats, idis, nohomog,
     &                   ierrsmooth, iramp
        common /mpistats/iISend, iISendScal, iIRecv, iIRecvScal, 
     &                   iWaitAll,iWaitAllScal, iAllR, iAllRScal,
     &                   impistat, impistat2, rmpitmr,
     &                   rISend, rISendScal, rIRecv, rIRecvScal, 
     &                   rWaitAll, rWaitAllScal, rAllR, rAllRScal, 
     &                   rCommu, rCommuScal

        common /memstats/rheap,rheapavail,rstack,rstackavail,rshared,
     &                   rpersist,rguard,rmmap,rmemstats

        common /spebcvr/ irscale, intpres, plandist,
     &            thetag, ds, tolerence, radcyl, rbltin, rvscal

c 
c.... common blocks
c
c nsrflist is a binary switch that tells us if a given srfID should be
c included in the consistent flux calculation.  It starts from zero
c since we need to be able to handle/ignore surfaces with no srfID attached
c
c flxID(numfluxes,nIDs+1)
c numfluxes = area, mass, fx, fy, fz, heat, scalar_flux_{1,2,3,4}
c nIDs currently set to MAXSURF, each surface has its own
c
        common /aerfrc/ flxID(10,0:MAXSURF), flxIDsclr(4,MAXSURF),
     &                  Force(3,MAXSURF),
     &                  HFlux(MAXSURF),    nsrfCM,
     &                  nsrflist(0:MAXSURF),
     &                  isrfIM,  irankfilesforce(0:MAXSURF)
c
        common /astore/ a(100000)
c
        common /mbndnod/ mnodeb(9,8,3)
c
	      integer, target ::  nlwork
c
c...........................................................................
        common /ctrlvari/ iI2Binlet, isetOutPres, isetInitial
        
        common /Ductvari/  BlowingVelDuct,
     &                    BlowingIniMdotDuct,
     &                    BlowingFnlMdotDuct,
     &                    suctionVbottom, 
     &                    suctionVside_lower,
     &                    suctionVside_upper, 
     &                    suctionVtop, 
     &                    blowerVelocity, 
     &                    blowerTemperature, 
     &                    blowerEV,
     &                    isetOutletID, 
     &                    isetInitial_Duct,
     &                    isetInlet_Duct, 
     &                    isetSuctionID_Duct,
     &                    isetBlowerID_Duct,  
     &                    iDuctgeometryType, 
     &                    iStraightPrint,
     &                    isetEV_IC_BC,
     &                    isetEVramp,
     &                    isetBlowing_Duct,
     &                    ifixBlowingVel_Duct, 
     &                    nBlowingStepsDuct
        real*8 inletVelX
        common /ctrlvar/  inletVelX,   outPres1, 
     &                    xvel_ini,    yvel_ini,    zvel_ini,
     &                    temp_ini,    pres_ini,    evis_ini

        common /Ductvar/   evis_IC_BC,
     &                    EVrampXmin, 
     &                    EVrampXmax,
     &                    EVrampMin,
     &                    EVrampMax
c...........................................................................

c
        common /alevar/ raleF,   raleA,   raleX,   raleY,   raleLx,
     &                  raleLy, raleRx,   raleRy,  ialeD,   ialeT
c
        common /levlset/ epsilon_ls, epsilon_lsd, dtlset, iLSet, 
     &                   ivconstraint, iExpLSSclr1, iExpLSSclr2

c 
        common /datpnt/ mshp,   mshgl,  mwght,  mshpb,  mshglb, mwghtb,
     &                  mmut,   mrhot,  mxst
c
        common /melmcat/ mcsyst, melCat, nenCat(8,3),    nfaCat(8,3)
c
c
        common /mintpar/ indQpt(3,3,4),  numQpt(3,3,4),
     &                  intmax
c /*         common /andres/ fwr1,ngaussf,idim,nlist */

        common /itrpar/ eGMRES, lGMRES, lGMRESs, iKs, iKss,    ntotGM,
     &                  ntotGMs, ntotGMelas
c
        common /point / mbeg,   mend,   mprec
c
        common /precis/ epsM,   iabres
c
        common /resdat/ resfrt, resfrts
c
        common /solpar/ imap,   ivart,  iDC,    iPcond, Kspace, nGMRES,
     &                  iconvflow, iconvsclr, idcsclr(2)
c
        common /msympar/ indsym(5,5)
c
        common /timpar/ LCtime, ntseq
c
        common /incomp/ numeqns(100), minIters, maxIters, 
     &                  iprjFlag,     nPrjs,    ipresPrjFlag, nPresPrjs,
     &                  prestol,      statsflow(6), statssclr(6),
     &                  iverbose
c
        character*80    title,  ititle
        common /title / title,  ititle
c
        character*8     machin
        parameter     ( machin = 'RS/6000 ' )
        parameter     ( machfl = 4 )
 
c
c----------------------------------------------------------------------
c
c.... element pointers
c
c mmat   (MAXBLK)  : pointer to interior element material number
c mmatb  (MAXBLK)  : pointer to boundary element material number
c mien   (MAXBLK)  : pointer to ien array
c mienb  (MAXBLK)  : pointer to ienb array
c miBCB  (MAXBLK)  : pointer to iBCB array
c mDt    (MAXBLK)  : pointer to Dt array
c mDC    (MAXBLK)  : pointer to DC array
c mBCB   (MAXBLK)  : pointer to BCB array
c mstiff (MAXBLK)  : pointer to stiff array
c
c----------------------------------------------------------------------
c
c.... common /aerfrc/   : aerodynamic forces
c
c Force(3)      : components of the aerodynamic forces
c HFlux         : total heat flux
c
c----------------------------------------------------------------------
c
c.... common /astore/   : the dynamic memory allocation area
c
c a(...)        : the blank array used for front-end data storage
c
c----------------------------------------------------------------------
c
c.... common /bndnod/   : boundary nodes of boundary elements
c
c mnodeb (9,8,3) : boundary nodes of each element category and dimension
c
c.... common /datpnt/   : front-end data pointers
c
c mshp          : pointer to shape-functions 
c mshgl         : pointer to local-grad-shape-functions
c mwght         : pointer to quadrature weights
c mshpb         : pointer to shape-functions of boundary elements
c mshglb        : pointer to local-grad-shape-functions of bound. elem.
c mwghtb        : pointer to quadrature weights of bound. elements
c mmut          : pointer to table mu  = mu  (p,T)
c mrhot         : pointer to table rho = rho (p,T)
c mxst          : pointer to table xs  = xs  (p,T)
c
c----------------------------------------------------------------------
c
c.... common /elmcat/   : element category information
c
c mcsyst        : maximum number of element coordinate system
c melCat        : maximum number of element categories
c nenCat (8,3)  : number of nodes for each category and dimension
c nfaCat (8,3)  : number of faces for each category and dimension
c
c----------------------------------------------------------------------
c
c.... common /shpdat/   : hierarchic shape function quadrature data
c
c Qpt  (3,MAXQPT)  : interior element quadrature points (xi,eta,zeta)
c Qwt  (MAXQPT)    : interior element quad. weights
c Qptb (2,MAXQPT)  : boundary element quad. pnts.
c Qwtb (MAXQPT)    : boundary element quad. weights
c nshape           : number of interior element shape functions
c nshapeb          :   "    "  boundary  "        "       "
c nshapeif         :   "    "  interface "        "       "
c ngauss           : number of interior element integration points
c ngaussb          :   "    "  boundary  "        "           "
c ngaussif         :   "    "  interface "        "           "
c----------------------------------------------------------------------
c
c.... common /intpar/   : integration parameters
c
c Qpt   (4,*)   : xi, eta, zeta, weight of quadrature points
c indQpt(3,3,4) : index to quadrature points for a given rule
c numQpt(3,3,4) : number of quadrature points for a given rule
c intmax        : number of allowable spatial integ. points per nsd
c
c----------------------------------------------------------------------
c
c.... common /itrpar/   : Preconditioned GMRES parameters
c
c eGMRES        : finite difference interval
c lGMRES        : number of GMRES cycles
c iKs           : current Krylov vector
c ntotGM        : total number of GMRES iterations
c
c----------------------------------------------------------------------
c
c.... common /itrpnt/   : Preconditioned GMRES array pointers
c
c mHBrg         : pointer to Hessenberg matrix
c meBrg         : pointer to Hessenberg's RHS matrix
c myBrg         : pointer to minimize solution matrix
c mRcos         : pointer to Rotation Cosine of QR algorithm
c mRsin         : pointer to Rotation Sine   of QR algorithm
c
c----------------------------------------------------------------------
c
c.... common /matpar/   : material constants
c
c pr            : Prandtl number
c Planck        : Planck's constant
c Stefan        : Stefan's constant (for radiation)
c Nh            : Avogadro's number
c Rh            : universal gas constant
c Rgas          : specific gas constant
c gamma         : specific heat ratio
c gamma1        : gamma - 1
c s0            : reference specific entropy
c const         : special constant
c xN2           : mole fraction of diatomic nitrogen
c xO2           : mole fraction of diatomic oxygen
c yN2           : mole fraction of diatomic nitrogen
c yO2           : mole fraction of diatomic oxygen
c Msh  (5)      : molar mass of species
c cpsh (5)      : molar heat at constant pressure of species
c s0sh (5)      : molar reference entropy of species
c h0sh (5)      : molar heat of formation of species
c Rs   (5)      : specific gas constant of species
c cps  (5)      : specific heat at constant pressure of species
c cvs  (5)      : specific heat at constant volume of species
c h0s  (5)      : specific heat of formation of species
c Trot (5)      : characteristic rotational temperature of species
c sigs (5)      : symmetry factor of species
c Tvib (5)      : characteristic vibrational temperature of species
c g0s  (5)      : ground degeneracy of electronic energy
c dofs (5)      : degrees of freedom of species
c ithm          : thermodynamic property flag
c
c----------------------------------------------------------------------
c
c.... common /matdat/   : material data
c
c datmat (3,5,2) : material data
c matflg (5,100)   : material type flag
c nummat           : number of materials
c mexist           : flag indicating the presence of MATERIAL.DAT
c
c----------------------------------------------------------------------
c
c.... common /point /   : dynamic storage pointer management data
c
c mbeg          : pointer to the beginning of the free storage
c mend          : pointer to the end of the storage
c mprec         : precision of the floating point data
c
c----------------------------------------------------------------------
c
c.... common /precis/   : finite difference interval data
c
c epsM          : square root of machine precision
c iabres        : absolute value residual flag
c
c----------------------------------------------------------------------
c
c....common /resdat/    : residual statistics data
c
c resfrt        : first residual of convergence
c
c----------------------------------------------------------------------
c
c.... common /solpar/   : solution parameters
c
c imap          : permutation mapping flag
c ivart         : variational formulation type
c iDC           : DC type
c iPcond        : type of preconditioner
c Kspace        : dimension of Krylov space
c nGMRES        : maximum number of GMRES iterations
c
c----------------------------------------------------------------------
c
c.... common /sympar/   : symmetric storage parameters
c
c indsym (5,5)  : mapping from 2D storage to symmetric one
c
c----------------------------------------------------------------------
c
c.... common /timpar/   : time integration parameters
c
c LCtime        : local time stepping flag
c ntseq         : number of time sequences
c
c----------------------------------------------------------------------
c
c.... common /timer1/   : timer parameters
c.... common /timer2/   : timer parameters
c.... common /timer3/   : timer parameters
c
c ccode(13)     : timing entities codes
c flops         : flop counter
c gbytes        : byte counter for gather operation
c sbytes        : byte counter for scatter operation
c iclock        : wall-clock time (in milliseconds)
c icd           : number of timing entities
c icode         : current timer code
c icode2        : last timer code
c icode3        : next-to-last timer code
c cpu(11)       : cpu time of each entity
c cpu0(11)      : initial cpu time of each entity
c nacess(11)    : number of times each entity is accessed
c
c----------------------------------------------------------------------
c
c.... common /title /   : problem title
c
c title         : problem title
c ititle        : problem title (with form feed)
c
c----------------------------------------------------------------------
c
c.... common /avging / : nfath
c 
c nfath         : total number of global fathers over which certain
c                 quantities will be averaged
c 
c----------------------------------------------------------------------
c
c.... parameters        : machine data
c
c machin        : machine type
c                  (set parameter)
c machfl        : single precision floating point lenght in bytes
c                  (set parameter)
c
c----------------------------------------------------------------------
c
c.... parameters        : useful constants
c
c zero          : 0.0
c pt125         : 0.125
c pt25          : 0.25
c pt33          : 0.33 (1/3)
c pt39          : 2^(-4/3)
c pt5           : 0.5
c pt57          : 1/sqrt(3)
c pt66          : 0.66 (2/3)
c pt75          : 0.75
c one           : 1.0
c sqt2          : sqrt(2)
c onept5        : 1.5
c two           : 2.0
c three         : 3.0
c four          : 4.0
c five          : 5.0
c pi            : the magical number :-)
c 
c---------------------------------------------------------------------- 
c 
c Zdenek Johan, Winter 1991.
c 
c----------------------------------------------------------------------
