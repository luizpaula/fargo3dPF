//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Computeaccretionsolidsequilibrium_cpu(real dt, int kp) {

//<USER_DEFINED>
  INPUT(Sigmam);
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(Sigmam);
//<\USER_DEFINED>

//<EXTERNAL>
  real* sigmam = Sigmam->field_cpu;
  real* cs = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* x = Sys->x_cpu;
  real* y = Sys->y_cpu;
  real* z = Sys->z_cpu;
  real* m = Sys->mass_cpu;
  real* mcore = Sys->masscore_cpu;
  real* sigmarms = Sys->sigmarms_cpu;
  real* eccrms = Sys->eccrms_cpu;
  real* incrms = Sys->incrms_cpu;
  real* rcore = Sys->rcore_cpu;
  real* rcap = Sys->rcap_cpu;
  real* accsolids = Sys->accsolids_cpu;
  real* masscrit = Sys->masscrit_cpu;
  real* massenv  = Sys->massenv_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ; 
//<\EXTERNAL>

//<INTERNAL>
  real delta;
  real rpl;
  real RHpl;
  real dMdA;
  int kpl;
  int i;
  int j;
  int k;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real ASPECTRATIO(1);
// real B(1);
// real GAMMA(1);
// real RICE(1);
// real RHOROCKY(1);
// real RHOICE(1);
// real RM(1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++) {
#endif
//<#>
           rpl =  sqrt(x[kp]*x[kp] + y[kp]*y[kp] + z[kp]*z[kp]);
           RHpl = rpl*pow(m[kp]/(3.0*MSTAR),(1.0/3.0));

#ifdef SPHERICAL
           delta =  ymed(j)*sin(zmed(k)) - rpl;
#endif
#ifdef CYLINDRICAL
           delta =  ymed(j) - rpl;
#endif
#ifdef CARTESIAN
           delta =  sqrt(XC*XC+YC*YC) - rpl;
#endif

           if (delta < 0) delta = (-1.0)*delta;

           if (delta < (B/2.0)*RHpl)
              {
              dMdA = accsolids[kp]*dt/(2.0*PI*rpl*B*RHpl);
              sigmam[l] = sigmam[l] - dMdA;
              if (sigmam[l] < 0.0) sigmam[l] = 0.0;
              }
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>

//<LAST_BLOCK>

#ifdef __GPU
  real* x = Sys->x_cpu;
  real* y = Sys->y_cpu;
  real* z = Sys->z_cpu;
  real* m = Sys->mass_cpu;
  real* mcore = Sys->masscore_cpu;
  real* sigmarms = Sys->sigmarms_cpu;
  real* eccrms = Sys->eccrms_cpu;
  real* incrms = Sys->incrms_cpu;
  real* rcore = Sys->rcore_cpu;
  real* rcap = Sys->rcap_cpu;
  real* accsolids = Sys->accsolids_cpu;
  real* masscrit = Sys->masscrit_cpu;
  real* massenv  = Sys->massenv_cpu;
#endif

  real year;
  real rcy;
  real dy;
  real temp = 0.0;
  int index;
  real r;
  real RH;
  real Porbital;
  real ecc_red, inc_red;
  real Beta;
  real IFBETA;
  real IGBETA;
  real Phigh;
  real Pmed;
  real Plow;
  real Pcoll;
  real Mesc;
  real fcap = 0.1;
  real rhog;
  real dp;
  real vrel;
  real vk;
  real sigmacrit;
  real L;
  real vslocal;
  real densg;
  real Pgas;
  real RBondi;
  real Tgas;
  real W;
  real sigmat;
  real gamma;
  real c1;
  real c2;

  real mpl;

  #ifdef MKS
     year = 31536000.0;
  #endif
  #ifdef CGS
     year = 31536000.0;
  #endif
  #if !(defined(MKS) || defined (CGS))
     year = 31536000.0/(sqrt(R0_MKS*R0_MKS*R0_MKS/G_MKS/MSTAR_MKS));
  #endif

  r = sqrt(x[kp]*x[kp] + y[kp]*y[kp] + z[kp]*z[kp]);
  RH = r*pow(m[kp]/(3.0*MSTAR),(1.0/3.0));

  #ifdef Y
     dy = Ymin(1) - Ymin(0);
     #ifdef CARTESIAN/*y -> y*/
        index = (int)((y[kp]-Ymin(0))/dy);
     #endif
     #ifdef CYLINDRICAL/*y -> r*/
        rcy = sqrt(x[kp]*x[kp] + y[kp]*y[kp]);
        index = (int)((rcy-Ymin(0))/dy);
     #endif
     #ifdef SPHERICAL/*y -> r*/
        index = (int)((r-Ymin(0))/dy);
     #endif
  #endif

  #ifndef __GPU
  //for MPI situation - subgrids
  if (index < 0) index = 0;
  if (index > Ny) index = 0;
  #endif

  if (index != 0) 
     {
     sigmarms[kp]  = reduction_full_SUM(Sigmam,index,index+1,NGHZ,Nz+NGHZ)/NX;
     vslocal = reduction_full_SUM(Energy,index,index+1,NGHZ,Nz+NGHZ)/NX;
     densg  = reduction_full_SUM(Density,index,index+1,0,1)/NX;
     }
  else 
     {
     sigmarms[kp] = 0.0;
     vslocal = 0.0;
     densg = 0.0;
     }

  #ifndef __GPU
    MPI_Allreduce(&sigmarms[kp],&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    sigmarms[kp] = temp;
    MPI_Allreduce(&vslocal,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    vslocal = temp;
    MPI_Allreduce(&densg,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   
    densg = temp;
  #endif

  /*VOLUMETRIC DENSITY OF GAS IN CELL CONTAINING THE PLANET -- ARMITAGE - PLANET FORMATION - EQ. 2.9*/ 
  #ifdef Z
     rhog = densg;
     densg = (rhog/(1.0/sqrt(2.0*PI)))*(ASPECTRATIO*r);
  #else 
     rhog = (1.0/sqrt(2.0*PI))*(densg/(ASPECTRATIO*r));
  #endif 

  #ifdef ISOTHERMAL
     vslocal = vslocal;
     Pgas = rhog*vslocal*vslocal;
     gamma = 5.0/3.0;
  #endif
  #ifdef ADIABATIC
     vslocal = sqrt(GAMMA*(GAMMA - 1.0)*vslocal/dens[l]); 
     Pgas = (GAMMA-1.0)*vslocal;
     gamma = GAMMA;
  #endif
  #ifdef POLYTROPIC
     vslocal = sqrt(vslocal*GAMMA*pow(dens[l],GAMMA-1.0));
     Pgas = vslocal*pow(vslocal,GAMMA);
     gamma = GAMMA;
  #endif

  /*RADIUS OF CORE*/ 
  rcore[kp] = pow(3.0*mcore[kp]/(4.0*PI*RHOM),1.0/3.0); 

  /*RADIUS OF CAPTURE*/
  if (r < RICE) dp = RHOROCKY;  //after line ice
  else dp = RHOICE; //before line ice
  vk = sqrt(G*MSTAR/r);
  //vrel = vk*sqrt(5.0*eccrms[kp]*eccrms[kp]/8.0 + 1.0*incrms[kp]*incrms[kp]/2.0);
  //sigmacrit = (dp*RM*vrel*vrel)/(3.0*G*mcore[kp]*rhog);
  sigmacrit = ((6.0 + eccrms[kp]*eccrms[kp]*(RH/r)*(RH/r))*RM*dp)/(9.0*RH)/rhog;
  L = (G*mcore[kp]/rcore[kp])*accsolids[kp];
  RBondi = (G*m[kp]/(gamma*vslocal*vslocal));
  Tgas = (vslocal*vslocal)/R_MU;
  W = (3.0*KE*L*Pgas)/(64.0*PI*STEFANK*G*mcore[kp]*pow(Tgas,4.0));
  sigmat = (1.0/(5.0*W));

  if (sigmacrit < sigmat)
     {
     c1 = 1.0 + (2.0*W*(sigmacrit - 1.0) + log(sigmacrit))/(gamma);
     rcap[kp] = RBondi/c1;
     }
  else
     {
     c1 = 1.0 + 2.0*W*(sigmat-1.0) + log(sigmat);
     c2 = 1.0/c1 + (4.0/gamma)*pow(4.0*W,1.0/3.0)*(pow(sigmacrit,1.0/3.0) - pow(sigmat,1.0/3.0));
     rcap[kp] = RBondi/c2;
     }

  if (rcap[kp] < rcore[kp]) rcap[kp] = rcore[kp];
  if (rcap[kp] > RH) rcap[kp] = RH;

   /*PERIOD ORBITAL*/ 
   Porbital = (2.0*PI)*pow((r*r*r)/(G*(MSTAR+mcore[kp])),0.5);

   /*WE CALCULATE THE MASS OF THE PLANETESIMAl**/ 
   mpl = dp*(4.0/3.0)*PI*pow(RM,3.0);

   /*WE CALCULATE THE ECCENTRICITY AND INCLINATION OF THE PLANETESIMAl (see Fortier e.t al. 2012 - Eq. 56 and 57*/ 
   eccrms[kp] = 1.7*(pow(mpl,1.0/15.0)*pow(mcore[kp],1.0/3.0)*pow(RHOM,2.0/15.0))*pow(B,-1.0/5.0)*pow(CD,-1.0/5.0)*pow(rhog,-1.0/5.0)*pow(MSTAR,-1.0/3.0)*pow(r,-1.0/5.0);
   incrms[kp] = eccrms[kp]/2.0;

   /*ECCENTRICITY AND INCLINATION REDUCED (see Fortier et. al., 2012)*/  
   ecc_red = r*eccrms[kp]/RH;
   inc_red = r*incrms[kp]/RH;

   /*VALOR OF BETA (see Fortier et. al., 2012)*/ 
   Beta = inc_red/ecc_red;

   /*IF e IG (Fortier et. al, 2013 - Eq 26 and Eq 27)*/ 
   IFBETA = (1.0 + 0.95925*Beta + 0.77251*Beta*Beta)/(Beta*(0.13142 + 0.12295*Beta));
   IGBETA = (1.0 + 0.39996*Beta)/(Beta*(0.0369 + 0.04833*Beta + 0.00687*Beta*Beta));

   /*VALOR FOR PROBABILITIES (Fortier et. al, 2013 - Eq 23, 24 and Eq 27)*/ 
   Phigh = (pow((rcap[kp] + RM),2.0)/(2.0*PI*RH*RH))*(IFBETA + (6.0*RH*IGBETA)/((rcap[kp] + RM)*ecc_red*ecc_red));
   Pmed = (pow((rcap[kp] + RM),2.0)/(4.0*PI*RH*RH*inc_red))*(17.3 + (232.0*RH)/(rcap[kp] + RM));
   Plow = 11.3*pow(((rcap[kp] + RM)/RH),0.5); 

   /*VALOR FOR PCOLL (Fortier et. al, 2013 - Eq 28)*/ 
   Pcoll = pow(Phigh,-2.0) + pow(Plow,-2.0);//PCOLL (Fortier et. al, 2013 - Eq 28)
   Pcoll = pow(Pcoll,-0.5);
   Pcoll = MIN(Pmed,Pcoll);

   if (accsolids[kp] < 0.0) accsolids[kp] = 0.0; 
   else accsolids[kp] = (2.0*PI*sigmarms[kp]*RH*RH/Porbital)*Pcoll;

   Mesc = (2.0*MSTAR*rcore[kp])/(sqrt(fcap)*r);
   if (mcore[kp] > Mesc) accsolids[kp] = 0.0; //scattering

   mcore[kp] = mcore[kp] + accsolids[kp]*dt;
   m[kp] = m[kp] + accsolids[kp]*dt;

   if (massenv[kp] == 0.0) masscrit[kp] = 10.0*pow(accsolids[kp]/((MEARTH/year)*10.0e-6),0.25)*MEARTH;
   else masscrit[kp] = masscrit[kp];
//<\LAST_BLOCK>
}
