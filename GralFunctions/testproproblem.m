aw1=471;aw2=200;af1=471;af2=471;nw=3;


************************************************************************
*            : Modeling Story: Structure    *
*            -----------------------------------------                 *
************************************************************************
************************************************************************
*DATE: 04 06 2019									  
*BY: CAROLINA FRANCO 								  
************************************************************************
*  TEST 1
*
*           NUMERICAL MODEL OF A multiframe STORY
*             PERIODIC BOUNDARY CONDITIONS
*
* In this numerical model, the geometry of one story                 *
* structure will be built up by using periodic boundary conditions.    *
* That means that some degrees of freedoms are going to be constrained * 
* to the same value. Transverse displacement at the top will be allowed* 
* As well as rotation and vertical displacement but guaranteeing the   *
* same for both ends.*
************************************************************************
*DEFINYING THE DIMENSION IN 2D SPACE AND THE TYPE OF GEOMETRY
*______________________________________________________________________*
*DEFINITION OF ONLY ONE STORY
************************************************************************
OPTION ECHO 0 ; 
*FUNCTION TO CREATE THE GRID/ GRID
DEBPROC GRID  inip*POINT  lf*FLOTTANT  lw*FLOTTANT  nw*ENTIER;
        TABLN        = TABLE;
		TABUN        = TABLE;
             REPE BOUCL nw;
			 TABLN. &BOUCL= table;
                  vectr1 = POIN ((&BOUCL-1.)*lf) 0.;
                  TABLN. &BOUCL . 1 = inip PLUS vectr1; 
             FIN BOUCL;
			 REPE BOUCUP nw;
			 TABUN. &BOUCUP= table;
                  vectr2 = ((&BOUCUP-1.)*lf) (lw);
                  TABUN. &BOUCUP . 1 = inip PLUS vectr2; 
             FIN BOUCUP;
FINPROC TABLN TABUN;

*FUNCTION TO DEFINE THE NODES AND ELEMENTS OF THE STORY WALLS;
DEBPROC WALLS TABLN*TABLE TABUN*TABLE ;			
OPTI COUL TURQ;		 
		nx = DIME TABLN;
        ny = DIME (TABLN. 1);
        TABW = TABLE;
        REPE BOUCX nx ;
             TABW. &BOUCX = TABLE;
             REPE BOUCY ny;
                  POILN = TABLN. &BOUCX . &BOUCY;
                  POIUN = TABUN. &BOUCX . &BOUCY;
                  ELE  = DROI POILN POIUN;
                  TABW. &BOUCX . &BOUCY = ELE;
                  SI ((&BOUCX*&BOUCY) 'EGA' 1);
                      MAILW = ELE;
                  SINON;
                      MAILW = MAILW ET ELE;
                  FINSI;
             FIN BOUCY;
        FIN BOUCX;
FINPROC TABW MAILW;

*FUNCTION TO DEFINE THE NODES AND ELEMENTS OF THE STORY FLOORS;
DEBPROC FLOORS TABUN*TABLE;
OPTI COUL ROSE;
		npx = (DIME TABUN) - 1;
        TABFX = TABLE;
        REPE BOUCX npx ;
             TABFX. &BOUCX = TABLE;
                  POINL = TABUN. &BOUCX . 1;
                  POINR = TABUN. (&BOUCX+1) . 1;
                  FLOORX = DROI POINL POINR;
                  SI (&BOUCX 'EGA' 1);
                      MAILF = FLOORX;
                  SINON;
                      MAILF = MAILF ET FLOORX;
                  FINSI;
                  TABFX. &BOUCX . 1 = FLOORX;
         FIN BOUCX; 
FINPROC TABFX MAILF;

*FUNCTION TO DEFINE THE NODES AND ELEMENTS OF THE ALL THE STORIES;
DEBPROC STRUCT inip*POINT lf*FLOTTANT lw*FLOTTANT nw*ENTIER   N*ENTIER; 
        TABSTRUC     = TABLE;
        TABSTRUC.WAL = TABLE;
        TABSTRUC.FLO = TABLE;
		REPE BOUCN N;
			 UPi = POIN 0. ((&BOUCN-1)*lw);
			 inipi= inip PLUS UPi;
			 TABLNi1 TABUNi1 = GRID inipi lf lw nw;
             WALLPi  MAILWi = WALLS TABLNi1 TABUNi1;
             FLOORPi MAILFi = FLOORS TABUNi1;
             TABSTRUC.WAL.&BOUCN = wallpi;
             TABSTRUC.FLO.&BOUCN = floorpi;
             SI (&BOUCN 'EGA' 1);
                MESHSTRU = MAILWi ET MAILFi;
             SINON;
                MESHSTRU = MESHSTRU ET MAILWi ET MAILFi;              
             FINSI; 
        FIN BOUCN;  
FINPROC TABSTRUC MESHSTRU;

*...................................................................*
*          COMPUTATION FILE                                         *
*...................................................................*

OPTI 'DIME' 2 'ELEM' 'SEG2' MODE TRID ;
OPTI EPSI LINEAIRE;
OPTI TRAC X;

*LENGTHS AND MESH DENSITY DEFINITION;
lw=3000.;     								
lf=3000.;									
N=1;									
dens1 = 0.8;
dens dens1;			
							

inip= POIN 0. 0.;	


TABSTRUC MESHSTRU = STRUCT inip lf lw nw N; 

*TRAC (MESHSTRU ) 'TITR' '[0] STRCUTURAL MODEL';
* Elimination OF SUPERPOSED POINTS
* ---------------------------------
elim MESHSTRU 1.E-4;

*ELEMENTS RECOVERY
*________________________________
*IF THERE ARE DIFFERENT PROPERTIES BETWEEN INNER AND OUTER ELEMENTS



BOOL1= (nw 'EGA' 2);
BOOL2= (nw 'EGA' 3);
SI BOOL1;

MAILW1= (TABSTRUC.WAL . 1 . 1 . 1);
MAILF1= (TABSTRUC.FLO . 1 . 1 . 1);
MAILF2= MAILF1;
MAILW2= MESHSTRU ELEM TURQ;
MAILW2= DIFF MAILW2 MAILW1;
SINON;
    SI BOOL2;
    MAILW1= (TABSTRUC.WAL . 1 . 1 . 1) ET (TABSTRUC.WAL . 1 . (nw) . 1);
    MAILF1= (TABSTRUC.FLO . 1 . 1 . 1);
    MAILW2= MESHSTRU ELEM TURQ;
    MAILW2= DIFF MAILW2 MAILW1;
    MAILF2=  MESHSTRU ELEM ROSE;
    MAILF2=  DIFF MAILF2 MAILF1;

   SINON;
  MAILW1= (TABSTRUC.WAL . 1 . 1 . 1) ET (TABSTRUC.WAL . 1 . (nw) . 1);
  MAILF1= (TABSTRUC.FLO . 1 . 1 . 1) ET (TABSTRUC.FLO . 1 . (nw-1) . 1);
  MAILW2= MESHSTRU ELEM TURQ;
   MAILW2= DIFF MAILW2 MAILW1;
   MAILF2=  MESHSTRU ELEM ROSE;
   MAILF2=  DIFF MAILF2 MAILF1;
   FINSI;
FINSI;


*_______________________________

MESHSTRU=MAILW1 ET MAILW2 ET MAILF1 ET MAILF2;
ELIM MESHSTRU 1.E-4;


*BOUNDARY CONDITION POINTS RECOVERY 	(Points at Y=0 and Y=lw)

*Y=0.; Pbase1 = 0. Y; Pbase2 = 1. Y; ni=lw*N; Pn1 = 0. ni; Pn2 = 1. ni ;*
*MAILBCb = MESHSTRU POIN 'DROIT' Pbase1 Pbase2 0.000000001;*
*MAILBCn = MESHSTRU POIN 'DROIT' Pn1 Pn2 0.000000001;*
*________________________________*
FGROUP1i        = TABLE;
FGROUP1f        = TABLE;
REPE BOUCX nw ;
FGROUP1i. &BOUCX= (TABSTRUC.WAL . 1 . &BOUCX . 1) POIN 'INITIAL';
FGROUP1f. &BOUCX= (TABSTRUC.WAL . 1 . &BOUCX . 1) POIN 'FINAL';
SI (&BOUCX 'EGA' 1);
                      MAILBCb = FGROUP1i. &BOUCX;
                      MAILBCn = FGROUP1f. &BOUCX;
                  SINON;
                      MAILBCb = MAILBCb ET FGROUP1i. &BOUCX;
                      MAILBCn = MAILBCn ET FGROUP1f. &BOUCX;
                  FINSI;
FIN BOUCX;
*________________________________*

ALLPOIN= MAILBCb ET MAILBCn;

*TRAC (MESHSTRU ET MAILBCb) 'TITR' '[1] STRCUTURAL MODEL';
*LIST MESHSTRU;




*...................................................................*
*          Structural Behavior Definition                           *
*...................................................................*
* Properties
* -------------------------------
h=1000.;
Arw1   = aw1*h;
Arw2   = aw2*h;
Arf1= af1*h;
Arf2= af2*h;
IZw1  = h*(aw1**3)/12;
IZw2  = h*(aw2**3)/12;
IZf1  = h*(af1**3)/12;
IZf2  = h*(af2**3)/12;
IYw1  = aw1*(h**3)/12;
IYw2  = aw2*(h**3)/12;
IYf1  = af1*(h**3)/12;
IYf2  = af2*(h**3)/12;
tor1=Iyw1+Izw1;
tor2=Iyw2+Izw2;
tor3=Iyf1+Izf1;
tor4=Iyf2+Izf2;
E  = 30000;
Nu1  = 0.2;
Rho1 = 2.3E-9;


* MODEL TYPE DEFINITION
* --------------------
*MODEL
OPTI DIME 2 ELEM SEG2 MODE TRID;

MOD1 = MODEL MAILW1 MECANIQUE ELASTIQUE ISOTROPE POUT ;
MOD2 = MODEL MAILW2 MECANIQUE ELASTIQUE ISOTROPE POUT ;
MOD3 = MODEL MAILF1 MECANIQUE ELASTIQUE ISOTROPE POUT ;
MOD4 = MODEL MAILF2 MECANIQUE ELASTIQUE ISOTROPE POUT ;

SI BOOL1;
MODSTRU= mod1 et mod2 et mod3;
SINON;
MODSTRU= mod1 et mod2 et mod3 et mod4;
FINSI;




* MATERIAL TYPE DEFINITION
* --------------------
*CONCRETE
MAT1 = MATE MOD1 'YOUN' E 'NU' Nu1 'RHO ' Rho1;
MAT2 = MATE MOD2 'YOUN' E 'NU' Nu1 'RHO ' Rho1;
MAT3 = MATE MOD3 'YOUN' E 'NU' Nu1 'RHO ' Rho1;
MAT4 = MATE MOD4 'YOUN' E 'NU' Nu1 'RHO ' Rho1;


SI BOOL1;
MATSTRUC= MAT1 et MAT2 et MAT3;
SINON;
MATSTRUC= MAT1 et MAT2 et MAT3 et MAT4;
FINSI;

* GEOMETRY DEFINITION
* --------------------
*WALLS AND FLOORS
CASTRU1 = CARA MOD1 'SECT' Arw1 'INRY' Iyw1 'INRZ' Izw1 'TORS' tor1 ;
CASTRU2 = CARA MOD2 'SECT' Arw2 'INRY' Iyw2 'INRZ' Izw2 'TORS' tor2;
CASTRU3 = CARA MOD3 'SECT' Arf1  'INRY' Iyf1 'INRZ' Izf1 'TORS' tor3;
CASTRU4 = CARA MOD4 'SECT' Arf2  'INRY' Iyf2 'INRZ' Izf2 'TORS' tor4;


SI BOOL1;
CASTRUC= CASTRU1 ET CASTRU2 ET CASTRU3;
SINON;
CASTRUC= CASTRU1 ET CASTRU2 ET CASTRU3 et CASTRU4;
FINSI;


PROSTRUC = MATSTRUC ET CASTRUC;

*OPTI 'DONNE' 5;

*.......................................................................*
*                  RESTRAINTS DEFINITION                                *
*.......................................................................*

* Periodic Conditions (Condition B);
* -------------------------------------------------
PLANecl=BLOQ MESHSTRU 'UZ' 'RX' 'RY';
*TOP;
NTDISP= BLOQ MAILBCn 'UX';
*BOTTOM
BDISP = BLOQ MAILBCb 'UX';

REPE BOUCD nw;
POICONS=(FGROUP1i. &BOUCD) ET (FGROUP1f. &BOUCD) ;
        BLOR= RELA 'ENSE' 'RZ' POICONS;
        SI (&BOUCD 'EGA' 1);
                BLORY = BLOR;
             SINON;
                BLORY=BLORY ET BLOR;              
        FINSI; 
FIN BOUCD;
*.........................UNIT TRANSVERSE DISPLACEMENT................*
*Unit Transverse Displacement HERE
TDISPL = lw;
CLTDISPL= DEPI NTDISP TDISPL;
VDISP = BLOQ ALLPOIN 'UY';
CONSTU= PLANecl ET NTDISP ET BDISP ET VDISP ET BLORY;


*.......................................................................*
*                  MECHANICAL COMPUTATIONS                              *
*.......................................................................*

* Mass and Stiffness Matrices
* --------------------------------------------------------------
RIGSTRU = RIGI MODSTRU PROSTRUC;
MASSTRU = MASS MODSTRU PROSTRUC;

* Global Stiffness Matrix
* -------------------------------------------------------------
RIGTOT  = RIGSTRU ET CONSTU; 

DEU= RESO RIGTOT CLTDISPL;

*LIST DE1;

REA1= REAC DEU RIGTOT;
FX= EXCO 'FX' REA1;

OPTI 'SORT' 'SHEARFORCE_b.inp';
SORT 'AVS' FX;


DEF_U= DEFO MESHSTRU DEU 'BLEU';

DEF_INI=DEFO MESHSTRU DEU 0. 'GRIS';
*TRAC (DEF_INI ET DEF_U);

FIN;




