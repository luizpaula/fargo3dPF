//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Computereymol_cpu() {
  
//<USER_DEFINED>
  INPUT(Ecc);
  INPUT(Inc);
  INPUT(Numol);
  OUTPUT(Reymol);
//<\USER_DEFINED>

//<EXTERNAL>
  real* reymol  = Reymol->field_cpu;
  real* numol  = Numol->field_cpu;
  real* ecc = Ecc->field_cpu;
  real* inc = Inc->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>
  
//<INTERNAL>
  real vk;
  float vrel;
  real c1;
  real c2;
  real c3;
  int i;
  int j;
  int k;
  real r;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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

#ifdef SPHERICAL
       vk = sqrt(G*MSTAR/ymed(j));
#endif
#ifdef CYLINDRICAL
       vk = sqrt(G*MSTAR/sqrt(ymed(j)*ymed(j)+zmed(k)*zmed(k)));
#endif
#ifdef CARTESIAN
       vk = sqrt(G*MSTAR/sqrt(XC*XC+YC*YC+ZC*ZC));
#endif

      c1 = ecc[l];
      c1 = c1*ecc[l];
      c1 = (5.0/8.0)*c1;

      c2 = inc[l];
      c2 = c2*inc[l];
      c2 = (1.0/2.0)*c2;
    
      c3 = c1 + c2;

      vrel = vk*sqrt(c3);
      reymol[l] = (RM*vrel)/numol[l];
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
