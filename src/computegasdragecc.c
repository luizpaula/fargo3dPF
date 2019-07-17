//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Computegasdragecc_cpu() {
  
//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  INPUT(Rhogas);
  INPUT(Reymol);
  INPUT(Ecc);
  INPUT(Inc);
  OUTPUT(Dedt);
//<\USER_DEFINED>

//<EXTERNAL>
  real* cs  = Energy->field_cpu;
  real* dens  = Density->field_cpu;
  real* rhogas  = Rhogas->field_cpu;
  real* ecc  = Ecc->field_cpu;
  real* inc  = Inc->field_cpu;
  real* reymol  = Reymol->field_cpu;
  real* dedt  = Dedt->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>
  
//<INTERNAL>
  real vs;
  real numberH2;
  real lambda;
  real p;
  real q;
  real hg;
  real vk;
  real nfactor;
  real dist;
  real tdrag;
  real c0;
  real c1;
  real c2;
  real c3;
  real c4;
  float c5;
  real epsilon = 1.211;
  int i;
  int j;
  int k;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real SIGMASLOPE(1);
// real FLARINGINDEX(1);
// real ASPECTRATIO(1);
// real RHOM(1);
// real RM(1);
// real CD(1);
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

#ifdef ISOTHERMAL
       vs = cs[l];
#endif
#ifdef ADIABATIC
       vs = sqrt(gamma*(gamma-1.0)*cs[l]/dens[l]); 
#endif
#ifdef POLYTROPIC
       vs = sqrt(cs[l]*gamma*pow(dens[l],gamma-1.0));
#endif	

       numberH2 = (NA/MUH2)*rhogas[l];
       lambda = (1.0/CSH2)/numberH2;

       /*WE CALCULATE nfactor (see Takeuchi 2002 - Eq.17)*/
#ifdef z
       p = SIGMASLOPE;
#else
       p = (SIGMASLOPE - FLARINGINDEX);
#endif
       q = 2.0*FLARINGINDEX - 3.0;

#ifdef SPHERICAL
       vk = sqrt(G*MSTAR/ymed(j));
       dist = ymed(j)*sin(zmed(k));
       hg = (ASPECTRATIO*dist)*pow((dist/R0),FLARINGINDEX);
       nfactor = -pow(hg/dist,2.0)*(p + q + FLARINGINDEX*(ymed(j)*cos(zmed(k))*(ymed(j)*cos(zmed(k))/hg); 
#endif
#ifdef CYLINDRICAL
       vk = sqrt(G*MSTAR/sqrt(ymed(j)*ymed(j)+zmed(k)*zmed(k)));
       dist = ymed(j);
       hg = (ASPECTRATIO*dist)*pow((dist/R0),FLARINGINDEX);
       nfactor = -pow(hg/dist,2.0)*(p + q + FLARINGINDEX*zmed(k)*zmed(k)/hg);  
#endif
#ifdef CARTESIAN
       vk = sqrt(G*MSTAR/sqrt(XC*XC+YC*YC+ZC*ZC));
       dist = sqrt(XC*XC+YC*YC);
       hg = (ASPECTRATIO*dist)*pow((dist/R0),FLARINGINDEX);
       nfactor = -pow(hg/dist,2.0)*(p + q + FLARINGINDEX*zmed(k)*zmed(k)/hg);
#endif

       if((reymol[l] >= 20.0) && (RM >= lambda))
         {
         tdrag = (8.0*RHOM*RM)/(3.0*CD*vk*rhogas[l]);
         c0 = -(1.0/tdrag)*ecc[l];
         c1 = (9.0/4.0)*nfactor*nfactor;
         c2 = (9.0/(4.0*PI))*epsilon*epsilon*ecc[l]*ecc[l];
         c3 = (1.0/PI)*inc[l]*inc[l];
         c4 = c1 + c2 + c3;
         c5 = c0*sqrt(c4);
         dedt[l] = c5;
         }

       if((reymol[l] < 20.0) && (RM >= lambda))
         {
         dedt[l] = -(3.0/4.0)*(lambda/(RHOM*RM*RM));
         dedt[l] *= vs*rhogas[l];
         dedt[l] *= ecc[l];
         }

       if(RM < lambda)
         {
         dedt[l] = -(1.0/2.0)*(1.0/(RHOM*RM))*vs;
         dedt[l] *= rhogas[l]*ecc[l];
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
}
