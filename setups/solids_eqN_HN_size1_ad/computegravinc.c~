//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Computegravinc_cpu() {
  
//<USER_DEFINED>
  INPUT(Ecc);
  INPUT(Inc);
  INPUT(Didt);
  OUTPUT(Didt);
//<\USER_DEFINED>

//<EXTERNAL>
  real* ecc  = Ecc->field_cpu;
  real* inc  = Inc->field_cpu;
  real* didt  = Didt->field_cpu;
  int Nbpl = Sys->nb;
  real* xpl = Sys->x_cpu;
  real* ypl = Sys->y_cpu;
  real* zpl = Sys->z_cpu;
  real* mpl = Sys->mass_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>
  
//<INTERNAL>
  int kp;
  real r;
  real Porbital;
  real RH;
  real ecc_red;
  real inc_red;
  real Beta;
  real IQVS;
  real Lambda;
  real QVS;
  real delta;
  real f_delta;
  float c0;
  int i;
  int j;
  int k;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real B(1);
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

       for (kp = 0; kp < Nbpl; kp++)
           {
           r =  sqrt(xpl[kp]*xpl[kp] + ypl[kp]*ypl[kp] + zpl[kp]*zpl[kp]);
           Porbital = (2.0*PI)*sqrt((r*r*r)/(G*(MSTAR+mpl[kp])));
           RH = r*pow(mpl[kp]/(3.0*MSTAR),(1.0/3.0));

           ecc_red = r*ecc[l]/RH; 
           inc_red = r*inc[l]/RH;

           Beta = inc_red/ecc_red;

           IQVS = (0.71946 - Beta)/(0.21239 + 0.49764*Beta + 0.14369*Beta*Beta);

           Lambda = inc_red*(ecc_red*ecc_red + inc_red*inc_red)/12.0;

           QVS = (((4.0*inc_red*inc_red)+(0.2*inc_red*ecc_red*ecc_red*ecc_red))/(10.0*Lambda*Lambda*ecc_red));
           QVS = QVS*log(1.0 +  10.0*Lambda*Lambda*ecc_red);
           QVS = QVS + (72.0*IQVS/(PI*ecc_red*inc_red))*log(1.0 + Lambda*Lambda);

#ifdef SPHERICAL
           delta =  ymed(j)*sin(zmed(k)) - r;
#endif
#ifdef CYLINDRICAL
           delta =  ymed(j) - r;
#endif
#ifdef CARTESIAN
           delta =  sqrt(XC*XC+YC*YC) - r;
#endif

           if (delta < 0) delta = (-1.0)*delta;

           f_delta = pow(delta/(5.0*RH),5.0);
           f_delta = pow(1.0 + f_delta,-1.0);

           c0 = (1.0/2.0)*f_delta*(mpl[kp]/(3.0*B*MSTAR*Porbital))*QVS;

           didt[l] += c0/inc[l];
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
