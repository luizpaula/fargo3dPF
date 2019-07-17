//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Updateinc_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Didt);
  INPUT(Inc);
  OUTPUT(Inc);
  OUTPUT(Incquad);
//<\USER_DEFINED>

//<EXTERNAL>
  real* inc = Inc->field_cpu;
  real* didt = Didt->field_cpu;
  real* incquad = Incquad->field_cpu;
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
  int ll;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
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
      inc[l] = inc[l] + didt[l]*dt;
      incquad[l] = inc[l]*inc[l];
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
  int NbPlanets = Sys->nb;
  int kp;
  real* x = Sys->x_cpu;
  real* y = Sys->y_cpu;
  real* z = Sys->z_cpu;
  real* mass = Sys->mass_cpu;
  real* incrms = Sys->incrms_cpu;

  real r;
  real rpl;
  real dy;
  real RH;
  real dist;
  real temp = 0.0;
  int index;

  for (kp = 0; kp < NbPlanets; kp++)
     {
     r = sqrt(x[kp]*x[kp] + y[kp]*y[kp] + z[kp]*z[kp]);
     RH = r*pow(mass[kp]/(3.0*MSTAR),(1.0/3.0));

     #ifdef Y
        dy = Ymin(1) - Ymin(0);
        #ifdef CARTESIAN/*y -> y*/
           index = (int)((y[kp]-Ymin(0))/dy);
        #endif
        #ifdef CYLINDRICAL/*y -> r*/
           rpl = sqrt(x[kp]*x[kp] + y[kp]*y[kp]);
           index = (int)((rpl-Ymin(0))/dy);
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
        incrms[kp] = reduction_full_SUM(Incquad,index,index+1,NGHZ,Nz+NGHZ);
        incrms[kp] = sqrt(incrms[kp]/NX);
        }
     else 
        {
        incrms[kp] = 0.0;
        }

     #ifndef __GPU
       MPI_Allreduce(&incrms[kp],&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   
       incrms[kp] = temp;
     #endif
    }
//<\LAST_BLOCK>
}
