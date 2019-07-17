/** \file fondam.h

This file contains all the information needed to work in a given
system of units.  Contains fundamental constants used thorough the
code.
*/

/* R_MU : ideal gas constant divided by mean molecular weight (R/\mu) */

/* Note: the mean molecular weight assumed below has the fiducial
   value of 2.4 g/mol for protoplanetary disks. */

//Define constant
#define        PI  3.141592653589793
#define        NA  6.022e23    //Avogadro constant

//Scale free
#define      G_SF  1.0
#define  MSTAR_SF  1.0
#define     R0_SF  1.0
#define   R_MU_SF  1.0
#define   CSH2_SF  1.0e-19/(5.2*1.49597871e11)/(5.2*1.49597871e11) //CSH2_SF = (CSH2_MKS / R0_MKS*R0_MKS) * R0_SF*R0_SF    
#define   MUH2_SF  0.002/1.9891e30   //MUH2_SF = (MUH2_MKS / (MSTAR_MKS/mol)) * (MSTAR_SF/mol) = (MUH2_MKS / MSTAR_MKS) * MSTAR_SF 
#define    MU0_SF  1.0
#define MEARTH_SF  5.973332e24/1.9891e30      //MEARTH_SF = (MEARTH_MKS / MSTAR_MKS) * MSTAR_SF     
#define KE_SF      1.0e-3*(1.9891e30/(5.2*1.49597871e11)/(5.2*1.49597871e11))

#define      G_MKS  6.674e-11
#define  MSTAR_MKS  1.9891e30
#define     R0_MKS  (5.2*1.49597871e11)
#define   R_MU_MKS  3460.0
#define   CSH2_MKS  1.0e-19 //cross section in m²
#define   MUH2_MKS  0.002  //molecular mass of H2 0.002 kg/mol    
#define    MU0_MKS  1.25663706143591e-6   //B in Tesla
#define MEARTH_MKS  5.973332e24
#define KE_MKS      1.0e-3    //opacity 0.001 m2/kg

#define      G_CGS  6.674e-8
#define  MSTAR_CGS  1.9891e33
#define     R0_CGS  (5.2*1.49597871e13)
#define   R_MU_CGS  36149835.0
#define   CSH2_CGS  1.0e-15 //cross section in cm²
#define   MUH2_CGS  2.0   //molecular mass of H2 2 g/mol
#define    MU0_CGS  12.5663706143591   //B in Gauss
#define MEARTH_CGS  5.973332e27
#define KE_CGS      0.01      //opacity 0.01 cm2/g  

#if !(defined(MKS) || defined (CGS))
#define		G	(G_SF)
#define     MSTAR       (MSTAR_SF)
#define        R0       (R0_SF)
#define      R_MU       (R_MU_SF)
#define      CSH2       (CSH2_SF)
#define      MUH2       (MUH2_SF)
#define       MU0       (MU0_SF)
#define    MEARTH       (MEARTH_SF)
#define    KE           (KE_SF)
#endif

//International System
#ifdef MKS
#define		G	(G_MKS)
#define     MSTAR       (MSTAR_MKS)
#define        R0       (R0_MKS)
#define      R_MU       (R_MU_MKS)
#define      CSH2       (CSH2_MKS)
#define      MUH2       (MUH2_MKS)
#define       MU0       (MU0_MKS)
#define    MEARTH       (MEARTH_MKS)
#define    KE           (KE_MKS)
#endif

//cgs
#ifdef CGS
#define		G	(G_CGS)
#define     MSTAR       (MSTAR_CGS)
#define        R0       (R0_CGS)
#define      R_MU       (R_MU_CGS)
#define      CSH2       (CSH2_CGS)
#define      MUH2       (MUH2_CGS)
#define       MU0       (MU0_CGS)
#define    MEARTH       (MEARTH_CGS)
#define    KE           (KE_CGS)
#endif


// Stefan's constant
#define STEFANK (5.6705e-5*pow(R_MU/R_MU_CGS,4.0)*pow(G/G_CGS,-2.5)*pow(MSTAR/MSTAR_CGS,-1.5)*pow(R0/R0_CGS,-0.5))

// Speed of light
#define C0      (2.99792458e10*sqrt(G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS))

// Cosmological microwave background's temperature
#define TCMB    (2.73*(G*MSTAR/R0/R_MU)/(G_CGS*MSTAR_CGS/R0_CGS/R_MU_CGS))

//#define TCMB    (30.0*(G*MSTAR/R0/R_MU)/(G_CGS*MSTAR_CGS/R0_CGS/R_MU_CGS))

#define THRESHOLD_STELLAR_MASS 0.05*MSTAR //Our arbitrary threshold to consider an object as stellar.
