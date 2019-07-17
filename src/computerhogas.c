//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Computerhogas_cpu() {
  
//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Rhogas);
//<\USER_DEFINED>

//<EXTERNAL>
  real* rhogas  = Rhogas->field_cpu;
  real* dens = Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>
  
//<INTERNAL>
  int i;
  int j;
  int k;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real ASPECTRATIO(1);
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

#ifdef Z
        rhogas[l] = dens[l];
#else
#ifdef SPHERICAL
	rhogas[l] = (1.0/sqrt(2.0*PI))*dens[l]/(ASPECTRATIO*ymed(j)*sin(zmed(k)));
#endif
#ifdef CYLINDRICAL
	rhogas[l] = (1.0/sqrt(2.0*PI))*dens[l]/(ASPECTRATIO*ymed(j));
#endif
#ifdef CARTESIAN
	rhogas[l] = (1.0/sqrt(2.0*PI))*dens[l]/(ASPECTRATIO*sqrt(XC*XC+YC*YC));
#endif
#endif

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
