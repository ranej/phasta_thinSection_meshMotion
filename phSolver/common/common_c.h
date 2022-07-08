/* Routine contains the structures for reading the user input through
 input_fform.cpp. The default values for all these variables are defined in
 input.config.

 Input variables that have been previously declared in common.h have to be
 re-declared here, in a consistant structure.*/

#include <FCMangle.h>

#define workfc FortranCInterface_GLOBAL_(workfc,WORKFC)
#define fronts FortranCInterface_GLOBAL_(fronts,FRONTS)
#define newdim FortranCInterface_GLOBAL_(newdim,NEWDIM)
#define timer4 FortranCInterface_GLOBAL_(timer4,TIMER4)
#define extrat FortranCInterface_GLOBAL_(extrat,EXTRAT)
#define spongevar FortranCInterface_GLOBAL_(spongevar,SPONGEVAR)
#define turbvar FortranCInterface_GLOBAL_(turbvar,TURBVAR)
#define turbvari FortranCInterface_GLOBAL_(turbvari,TURBVARI)
#define mpistats FortranCInterface_GLOBAL_(mpistats,MPISTATS)
#define memstats FortranCInterface_GLOBAL_(memstats,MEMSTATS)
#define spebcvr FortranCInterface_GLOBAL_(spebcvr,SPEBCVR)
#define aerfrc FortranCInterface_GLOBAL_(aerfrc,AERFRC)
#define laymesh FortranCInterface_GLOBAL_(laymesh,LAYMESH)
#define laymeshts FortranCInterface_GLOBAL_(laymeshts,LAYMESHTS)
#define snapmesh FortranCInterface_GLOBAL_(snapmesh,SNAPMESH)
#define meshquality FortranCInterface_GLOBAL_(meshquality,MESHQUALITY)
#define rigidbody FortranCInterface_GLOBAL_(rigidbody,RIGIDBODY)
#define rotatingband FortranCInterface_GLOBAL_(rotatingband,ROTATINGBAND)
#define astore FortranCInterface_GLOBAL_(astore,ASTORE)
#define conpar FortranCInterface_GLOBAL_(conpar,CONPAR)
#define ctrlvari FortranCInterface_GLOBAL_(ctrlvari,CTRLVARI)
#define ductvari FortranCInterface_GLOBAL_(ductvari,DUCTVARI)
#define ctrlvar FortranCInterface_GLOBAL_(ctrlvar,CTRLVAR)
#define ductvar FortranCInterface_GLOBAL_(ductvar,DUCTVAR)
#define alevar FortranCInterface_GLOBAL_(alevar, ALEVAR)
#define shpdat FortranCInterface_GLOBAL_(shpdat,SHPDAT)
#define datpnt FortranCInterface_GLOBAL_(datpnt,DATPNT)
#define elmpar FortranCInterface_GLOBAL_(elmpar,ELMPAR)
#define genpar FortranCInterface_GLOBAL_(genpar,GENPAR)
#define inpdat FortranCInterface_GLOBAL_(inpdat,INPDAT)
#define intdat FortranCInterface_GLOBAL_(intdat,INTDAT)
#define mio FortranCInterface_GLOBAL_(mio,MIO)
#define mioname FortranCInterface_GLOBAL_(mioname,MIONAME)
#define itrpar FortranCInterface_GLOBAL_(itrpar,ITRPAR)
#define itrpnt FortranCInterface_GLOBAL_(itrpnt,ITRPNT)
#define matdat FortranCInterface_GLOBAL_(matdat,MATDAT)
#define dgifinp FortranCInterface_GLOBAL_(dgifinp,DGIFINP)
#define mmatpar FortranCInterface_GLOBAL_(mmatpar,MMATPAR)
#define outpar FortranCInterface_GLOBAL_(outpar,OUTPAR)
#define point FortranCInterface_GLOBAL_(point,POINT)
#define precis FortranCInterface_GLOBAL_(precis,PRECIS)
#define propar FortranCInterface_GLOBAL_(propar,PROPAR)
#define resdat FortranCInterface_GLOBAL_(resdat,RESDAT)
#define solpar FortranCInterface_GLOBAL_(solpar,SOLPAR)
#define timdat FortranCInterface_GLOBAL_(timdat,TIMDAT)
#define timpar FortranCInterface_GLOBAL_(timpar,TIMPAR)
#define incomp FortranCInterface_GLOBAL_(incomp,INCOMP)
#define mtimer1 FortranCInterface_GLOBAL_(mtimer1,MTIMER1)
#define mtimer2 FortranCInterface_GLOBAL_(mtimer2,MTIMER2)
#define timer3 FortranCInterface_GLOBAL_(timer3,TIMER3)
#define title FortranCInterface_GLOBAL_(title,TITLE)
#define sclrs FortranCInterface_GLOBAL_(sclrs,SCLRS)
#define levlset FortranCInterface_GLOBAL_(levlset,LEVLSET)
#define nomodule FortranCInterface_GLOBAL_(nomodule,NOMODULE)
#define sequence FortranCInterface_GLOBAL_(sequence,SEQUENCE)
#define amgvarr FortranCInterface_GLOBAL_(amgvarr,AMGVARR)
#define amgvari FortranCInterface_GLOBAL_(amgvari,AMGVARI)

#define MAXBLK   50000
#define MAXSURF  1000  
#define MAXTS   100
#define MAXTOP   6
#define MAXTOPIF 12
#define MAXMAT  6
#define MAXPROP 10
#define MAXQPT   125
#define MAXSH    32
#define NSD      3
#define NSDSQ    9
#define machin   'RS/6000'
#define machfl   4
#define zero   0.0000000000000000000000000000000d0
#define pt125   0.1250000000000000000000000000000d0
#define pt25   0.2500000000000000000000000000000d0
#define pt33   0.3333333333333333333333333333333d0
#define pt39   0.3968502629920498686879264098181d0
#define pt5   0.5000000000000000000000000000000d0
#define pt57   0.5773502691896257645091487805020d0
#define pt66   0.6666666666666666666666666666667d0
#define pt75   0.7500000000000000000000000000000d0
#define one   1.0000000000000000000000000000000d0
#define sqt2   1.4142135623730950488016887242097d0
#define onept5   1.5000000000000000000000000000000d0
#define two   2.0000000000000000000000000000000d0
#define three   3.0000000000000000000000000000000d0
#define four   4.0000000000000000000000000000000d0
#define five   5.0000000000000000000000000000000d0
#define pi   3.1415926535897932384626433832795d0
#define inv1024sq 9.5367431640625e-7

#define ieos_ideal_gas 1
#define ieos_ideal_gas_2 2
#define ieos_ideal_gas_mixture 3
#define ieos_liquid_1  4
#define ieos_solid_1 10

#define idg_no_vi 1
#define idg_const_vi 2
#define idg_vieilles_burning 3
#define idg_clausius_clapeyron 4
#define idg_cavitation 5

#define idg_no_ramp 1
#define idg_linear_ramp 2

#ifdef __cplusplus
extern "C" {
#endif
  extern struct { 
    int master;
    int numpe;
    int myrank;
  } workfc ;

  extern struct { 
    int maxfront;
    int nlwork;
  } fronts ;

  extern struct { 
    long long int nshgt;
    long long int minowned;
    long long int maxowned;
    int numper;
    //int nshgt;
    int nshg0;
  } newdim ;

  extern struct { 
    double birth;
    double death;
    double comtim;
  } timer4 ;

  extern struct { 
    double ttim[100];
  } extrat ;

  extern struct {
    double zoutsponge, radsponge, zinsponge, grthosponge, grthisponge;
    double betamax;
    int spongecontinuity, spongemomentum1, spongemomentum2;
    int spongeenergy, spongemomentum3;
  } spongevar ;

  extern struct {
    double eles;
    double ylimit[9][3]; /* 9 = 5 + 4 = puvwT + 4Scalars */
    double rampmdot[3][2];
    double rmutarget;
    double pzero;
    double wtavei;
    double dtavei;
    double dke;
    double fwr1;
    double flump;
    double DES_SA_hmin;
    int ierrcalc;
    int ihessian;
    int itwmod;
    int ngaussf;
    int idim;
    int nlist;
    int nintf[MAXTOP];
  } turbvar ;

  extern struct {
    int irans, iles, idistcalc, isubmod;
    int ifproj;
    int i2filt;
    int modlstats;
    int idis;
    int nohomog;
    int ierrsmooth;
    int iramp;

/*      int itwmod; */
/*      double rtavei; */
/*      int ierrcalc; */
  } turbvari ;

  extern struct { 
    int iISend;
    int iISendScal;
    int iIRecv;
    int iIRecvScal;
    int iWaitAll;
    int iWaitAllScal;
    int iAllR;
    int iAllRScal;
    int impistat;
    int impistat2;
    double rmpitmr;
    double rISend;
    double rISendScal;
    double rIRecv;
    double rIRecvScal;
    double rWaitAll;
    double rWaitAllScal;
    double rAllR;
    double rAllRScal;
    double rCommu;
    double rCommuScal;
  } mpistats ;

  extern struct { 
    double rheap;
    double rheapavail;
    double rstack;
    double rstackavail;
    double rshared;
    double rpersist;
    double rguard;
    double rmmap;
  } memstats ;


  extern struct { 
    int irscale;
    int intpres;
    double plandist;
    double thetag;
    double ds;
    double tolerence;
    double radcyl;
    double rbltin;
    double rvscal;
  } spebcvr ;

  extern struct {
    double scdiff[5];
    double tdecay;
    int nsclr, isclr,nsolt, nosource;
    int consrv_sclr_conv_vel;
  } sclrs;

  extern struct {
    double flxID[MAXSURF+1][10] ;
    double flxIDsclr[MAXSURF][4];
    double Force[MAXSURF][3];
    double HFlux[MAXSURF];
    int nsrfCM;
    int nsrflist[MAXSURF+1];
    int isrfIM;
    int irankfilesforce[MAXSURF+1];
  } aerfrc ;

  extern struct {
    double blfactor;
    int numgc;
    int numgcnp;
    int gcBaseOpt;
  } laymesh ;

  extern struct {
    double tsfactor;
    int numts;
    int numtsnp;
    int tsBaseOpt;
  } laymeshts ;

  extern struct {
    int snapSurfFlag;
    int snapNumFaceTags;
    int snapFaceTags[MAXTS];
    int timeDepComp1Flag;
    int timeDepComp1ID;
  } snapmesh ;

  extern struct {
    double volMeshqTol;
    double faceMeshqTol;
    int autoTrigger;
    int triggerNow;
    double errorTolMass;
    double errorTolMomt;
    double errorTolEngy;
    double errorMaxMass;
    double errorMaxMomt;
    double errorMaxEngy;
    int errorEstimation;
    int errorTriggerEqn;
  } meshquality ;

  extern struct {
    double rb_prop[MAXTS][MAXTS];
    int numrbs;
    int rbsTags[MAXTS];
    int rbsMM[MAXTS];
    int rb_commuMotion;
    int rbsPcase[MAXTS];
  } rigidbody ;


  extern struct {
    int numRotBands;
    int numRotBandFaceTags;
    int rotBandTag[MAXTS][MAXTS];
  } rotatingband ;
  

  extern struct { 
    double a[100000];
  } astore ;

  extern struct { 
    int numnp;
    int numel;
    int numelb;
    int numelif;
    int numpbc;
    int nen;
    int nfaces;
    int numflx;
    int ndof;
    int iALE;
    int iSOLID;
    int icoord;
    int navier;
    int irs;
    int iexec;
    int necho;
    int ichem;
    int iRK;
    int nedof;
    int ndofelas;
    int nshg;
    int nnz;
    int istop;
    int nflow;
    int nelas;
    int nnz_tot;
    int idtn;
    int ncorpsize;
    int iownnodes;
    int usingpetsc;
    int elasModel;
    int elasFDC;
    int mesh2geom;
  } conpar ;
 
/*chen Sep 25 2009  Flow Control Parameters*/
  extern struct{
    int iI2Binlet;
    int isetOutPres;
    int isetInitial;
  } ctrlvari;

	extern struct{
		double BlowingVelDuct; 
		double BlowingIniMdotDuct;
		double BlowingFnlMdotDuct;
		double suctionVbottom;
		double suctionVside_lower;
		double suctionVside_upper;
		double suctionVtop;
		double blowerVelocity;
		double blowerTemperature;
		double blowerEV;
		int isetOutletID;
		int isetInitial_Duct;
		int isetInlet_Duct;
		int isetSuctionID_Duct;
		int isetBlowerID_Duct;
		int iDuctgeometryType;
		int iStraigtPrint;
		int isetEV_IC_BC;
		int isetEVramp;
		int isetBlowing_Duct;
		int ifixBlowingVel_Duct;  
		int nBlowingStepsDuct;
	}ductvari;

  extern struct{
    double inletVelX;
    double outPres1; 
    double xvel_ini;
    double yvel_ini;
    double zvel_ini;
    double temp_ini;
    double pres_ini;
    double evis_ini;
  } ctrlvar;

	extern struct{
		double evis_IC_BC;
		double EVrampXmin;
		double EVrampXmax;
		double EVrampMin;
		double EVrampMax;
	} ductvar;

  extern struct {
    double raleF;
    double raleA;
    double raleX;
    double raleY;
    double raleLx;
    double raleLy;
    double raleRx;
    double raleRy;
    int ialeD;
    int ialeT;
  } alevar;
 
  extern struct { 
    double epsilon_ls;
    double epsilon_lsd;
    double dtlset;
    int iLSet;
    int ivconstraint;
    int iExpLSSclr1;
    int iExpLSSclr2;
  } levlset;

  extern struct { 
    int nshape;
    int nshapeb;
    int nshapeif;
    int maxshb;
    int nshl;
    int nshlb;
    int nshl0;
    int nshl1;
    int nfath;
    int ntopsh;
    int nsonmax;
  } shpdat ;

  extern struct { 
    int mshp;
    int mshgl;
    int mwght;
    int mshpb;
    int mshglb;
    int mwghtb;
    int mmut;
    int mrhot;
    int mxst;
  } datpnt ;

  extern struct { 
    int lelCat;
    int lcsyst;
    int iorder;
    int nenb;
    int nelblk;
    int nelblb;
    int nelblif;
    int ndofl;
    int nsymdl;
    int nenl;
    int nfacel;
    int nenbl;
    int intind;
    int mattyp;
  } elmpar ;

  extern struct { 
    double E3nsd;
    int I3nsd;
    int nsymdf;
    int ndofBC;
    int ndiBCB;
    int ndBCB;
    int Jactyp;
    int jump;
    int ires;
    int iprec;
    int iprev;
    int ibound;
    int idiff;
    int lhs;
    int itau;
    int ipord;
    int ipred;
    int lstres;
    int iepstm;
    double dtsfct;
    double taucfct;
    int ibksiz;
    int iabc;
    int isurf;
    int idflx;
    double Bo;
    int EntropyPressure;
    int irampViscOutlet;
    int istretchOutlet;
    int iremoveStabTimeTerm;
    int iLHScond;
    int ndofBC2;
  } genpar ;

  extern struct { 
    double epstol[6];  /* 1+ max number of scalars  (beginning of the
                          end of time sequences) */
    double etolelas;
    double Delt[MAXTS];
    double CFLfl[MAXTS];
    double CFLsl[MAXTS];
    int nstep[MAXTS];
    int niter[MAXTS];
    int impl[MAXTS];
    double rhoinf[MAXTS];
    double rhoinfS[MAXTS];
    double rhoinf_B[MAXTS];//for solid
    double rhoinf_rb[MAXTS];//for rigid body motion
    int LHSupd[6];
    int loctim[MAXTS];
    double deltol[2][MAXTS];
    int leslib;
    int svLSFlag;
    int svLSType;
  } inpdat ;

  extern struct { 
    int iin;
    int igeom;
    int ipar;
    int ibndc;
    int imat;
    int iecho;
    int iout;
    int ichmou;
    int irstin;
    int irstou;
    int ihist;
    int iflux;
    int ierror;
    int itable;
    int iforce;
    int igraph;
    int itime;
    int iconserv;
  } mio ;

  extern struct { 
    double fin;
    double fgeom;
    double fpar;
    double fbndc;
    double fmat;
    double fecho;
    double frstin;
    double frstou;
    double fhist;
    double ferror;
    double ftable;
    double fforce;
    double fgraph;
    double ftime;
    double fconserv;
  } mioname ;

  extern struct { 
    double eGMRES;
    int lGMRES;
    int lGMRESs;
    int iKs;
    int iKss;
    int ntotGM;
    int ntotGMs;
    int ntotGMelas;
  } itrpar ;

  extern struct { 
    int matflg[MAXTS][MAXMAT];
    int mat_tag[MAXTS][MAXMAT];
    int mat_eos[MAXTS][MAXMAT];
    int nummat;
    int mexist;
    double mat_prop[MAXTS][MAXPROP][MAXMAT];
    double datmat[MAXTS][MAXPROP][MAXMAT];
    double datelas[2][1];
    double surface_tension_coeff;
    int surface_tension_flag;
    int datelas_volume_YM;
  } matdat ;

  extern struct {
    int phase_change_model;
    int vi_ramping;
    int i_w_normal;
    int i_if_dc;
    double ramp_time;
    double vi_mag;
    double s;
    double e;
    double emu;
    double ek;
    double h;
    double if_e_dc;
    double if_reaction_heat;
    double hfg_liquid, mw_liquid, T_boil_liquid;
    double burn_rate_exp, burn_rate_coeff, burn_rate_pref;
  } dgifinp;

  extern struct { 
    double Rgas, gamma, gamma1, mw;
    double pr;
    /* double pr, Planck, Stephan, Nh, Rh, Rgas;*/
    /*, const, xN2, xO2;*/
    /*double yN2,    yO2,    Msh[5], cpsh[5],s0sh[5],h0sh[5];*/
    /*double Rs[5],  cps[5], cvs[5], h0s[5], Trot[5],sigs[5];*/
    /*double Tvib[5],g0s[5], dofs[5],ithm;*/
  } mmatpar ;

  extern struct { 
    double ro;
    double vel;
    double temper;
    double press;
    double entrop;
    int ntout;
    int ioform;
    int iowflux;
    int iofieldv;
    char iotype[80];
    int ioybar;
    int nstepsincycle;
    int nphasesincycle;
    int ncycles_startphaseavg;
    int ivort;
    int icomputevort;
    int nsynciofiles;
    int nsynciofieldswriterestart;
    int iv_rankpercore;
    int iv_corepernode; 
    int input_mode; //FIXME -1:streams, 0:posix, >0:syncio
    int output_mode; //FIXME -1:streams, 0:posix, >0:syncio
    int conservation_probe;
    int write_residual;
    int imeshCFL;
    /*  int iostats; */
/*      int ipresref; */
  } outpar ;

  extern struct { 
    int mbeg;
    int mend;
    int mprec;
  } point ;

  extern struct { 
    double epsM;
    int iabres;
  } precis ;

  extern struct { 
    int npro;
  } propar ;

  extern struct { 
    double resfrt;
    double resfrts;
  } resdat ;

  extern struct { 
    int imap;
    int ivart;
    int iDC;
    int iPcond;
    int Kspace;
    int nGMRES;
    int iconvflow;
    int iconvsclr;
    int idcsclr[2];
  } solpar ;

  extern struct { 
    double time;
    double CFLfld;
    double CFLsld;
    double Dtgl;
    double Dtmax;
    double alpha;
    double etol;
    int lstep;  // read from numstart.dat and incremented every time step
    int ifunc;
    int itseq;
    int istep; //  how many steps (starting from 0 each run)
    int iter;
    int nitr;
    double almi;
    double alfi;
    double gami;
    double almBi; //for solid
    double alfBi; //for solid
    double gamBi; //for solid
    double flmpl;
    double flmpr;
    double dtol[2];
    int iCFLworst;
    int lskeep;
  } timdat ;

  extern struct { 
    int LCtime;
    int ntseq;
  } timpar ;

  extern struct { 
    int numeqns[100];
    int minIters;
    int maxIters;
    int iprjFlag;
    int nPrjs;
    int ipresPrjFlag;
    int nPresPrjs;
    double prestol;
    double statsflow[6];
    double statssclr[6];
    int iverbose;
  } incomp ;

  extern struct { 
    double ccode[13];
  } mtimer1 ;

  extern struct { 
    double flops;
    double gbytes;
    double sbytes;
    int iclock;
    int icd;
    int icode;
    int icode2;
    int icode3;
  } mtimer2 ;

  extern struct { 
    double cpu[11];
    double cpu0[11];
    int nacess[11];
  } timer3 ;

  extern struct { 
    double title;
    int ititle;
  } title ;

  extern struct {
    int intg[MAXTS][3];
  }intdat;

  extern struct {
    double bcttimescale;    
    double ValueListResist[MAXSURF+1];
    double rhovw;
    double thicknessvw;
    double evw;
    double rnuvw;
    double rshearconstantvw;
    double betai;
    int icardio;
    int itvn;
    int ipvsq;
    int numResistSrfs;
    int nsrflistResist[MAXSURF+1];
    int numImpSrfs;
    int nsrflistImp[MAXSURF+1];
    int impfile;
    int numRCRSrfs;
    int nsrflistRCR[MAXSURF+1];
    int ircrfile;
    int ideformwall;  
    int iwallmassfactor;
    int iwallstiffactor;
    int iviscflux;   
 } nomodule;

  extern struct {
    int seqsize;
    int stepseq[100];
  } sequence;

  extern struct {
    double strong_eps;      /* strong criterion Stuben factor    */
    double ramg_eps;        /* AMG convergence eps               */
    double ramg_relax;       /* relaxation factor Gauss-Seidel/Jac*/
    double ramg_trunc;      /* truncation select */
    double ramg_chebyratio; /* Eigen ratio for chebyshev smoothing */
 } amgvarr ;
  
  extern struct {
    int irun_amg;           /* Employ AMG feature solfar.f      */
    int irun_amg_prec;      /* Run AMG as preconditioner to CG */
    int iamg_verb;          /* amg verbosity flag                */
    int iamg_neg_sten;      /* neg only stencil or neg and pos   */
    int iamg_nlevel;        /* number of levels 2-V etc.         */
    int iamg_c_solver;     /* solve fine level iter. method     */
    int iamg_init;           /* setup flag */
    int iamg_setup_frez;    /* how many solfars to re setup amg */
    int iamg_interp;        /* interpolation select */
    int maxnev;             /* total eigenvectors used for ggb*/
    int maxncv;             /* total iterative vectors for ggb*/
    int iamg_smoother;      /* Smoother type */
    int mlsdeg;             /* Polynomial Smoothing (MLS) degree */
    int iamg_reduce;        /* Run a reduced case */
 } amgvari ;

#ifdef __cplusplus
}
#endif
