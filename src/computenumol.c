//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Computenumol_cpu() {
  
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Rhogas);
  INPUT(Energy);
  OUTPUT(Numol);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* rhogas = Rhogas->field_cpu;
  real* cs = Energy->field_cpu;
  real* numol  = Numol->field_cpu;
  real gamma = GAMMA;
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
  int i;
  int j;
  int k;
//<\INTERNAL>

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
       vs = (1.0/3.0)*cs[l];
#endif
#ifdef ADIABATIC
       vs = (1.0/3.0)*sqrt(gamma*(gamma-1.0)*cs[l]/dens[l]); 
#endif
#ifdef POLYTROPIC
       vs = (1.0/3.0)*sqrt(cs[l]*gamma*pow(dens[l],gamma-1.0));
#endif	

       numberH2 = (NA/MUH2)*rhogas[l];
       lambda = (1.0/CSH2)/numberH2;
       numol[l] = lambda*vs;
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
