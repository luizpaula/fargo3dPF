//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#define LEFT  0
#define RIGHT 1
#define DOWN  2
#define UP    3
//<\INCLUDES>

void boundary_ymin_cpu () {

//<USER_DEFINED>
INPUT(Density);
INPUT(Vx);
INPUT(Vy);
INPUT(Sigmam);
INPUT(Ecc);
INPUT(Inc);
OUTPUT(Density);
OUTPUT(Vx);
OUTPUT(Vy);
OUTPUT(Sigmam);
OUTPUT(Ecc);
OUTPUT(Inc);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int kact;
  int kgh;
  int lgh;
  int lghs;
  int lact;
  int lacts;
  int lacts_null;
//<\INTERNAL>

//<EXTERNAL>
  real* density = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* sigmam = Sigmam->field_cpu;
  real* ecc = Ecc->field_cpu;
  real* inc = Inc->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = NGHY;
  int size_z = Nz+2*NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real dx = Dx;
  real alpha = ALPHA;
  real aspectratio = ASPECTRATIO;
  real flaringindex = FLARINGINDEX;
  real sigmaslope = SIGMASLOPE;
  real omegaframe = OMEGAFRAME;
//<\EXTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	lgh = l;
	lghs = l;
	lact = i + (2*nghy-j-1)*pitch + k*stride;
	lacts = i + (2*nghy-j)*pitch + k*stride;
	lacts_null = i + nghy*pitch + k*stride;
	jgh = j;
	jact = (2*nghy-j-1);

	density[lgh] = density[lact];
	vx[lgh] = (sqrt(6.674e-11*1.9891e30/ymed(jgh)/ymed(jgh)/ymed(jgh)))*ymed(jgh)*sqrt(1.0+pow(aspectratio,2.0)*pow(ymed(jgh)/(5.2*1.49597871e11),2.0*flaringindex)*(2.0*flaringindex - 1.0 - sigmaslope)) - omegaframe*ymed(jgh);
	vy[lghs] = -3.0*alpha*ymed(jgh)*(sqrt(6.674e-11*1.9891e30/ymed(jgh)/ymed(jgh)/ymed(jgh)))*pow(aspectratio,2.0)*(flaringindex - sigmaslope + 1.0)*pow(ymed(jgh)/(5.2*1.49597871e11),2.0*flaringindex);
	vy[lacts_null] = -3.0*alpha*ymed(jgh)*(sqrt(6.674e-11*1.9891e30/ymed(jgh)/ymed(jgh)/ymed(jgh)))*pow(aspectratio,2.0)*(flaringindex - sigmaslope + 1.0)*pow(ymed(jgh)/(5.2*1.49597871e11),2.0*flaringindex);
	sigmam[lgh] = sigmam[lact];
	ecc[lgh] = ecc[lact];
	inc[lgh] = inc[lact];
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
