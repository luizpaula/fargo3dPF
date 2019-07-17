//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Accretiongas_cpu(real dt, int kp){
   
//<USER_DEFINED>
  INPUT(Vx);
  INPUT(Vy);
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(DeltaM);
  OUTPUT(Density);
//<\USER_DEFINED>

//<EXTERNAL>
  real* cs = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* Vt = Vx->field_cpu;
  real* Vr = Vy->field_cpu;
  real* massacc = DeltaM->field_cpu;
  real* x  = Sys->x_cpu;
  real* y  = Sys->y_cpu;
  real* z  = Sys->z_cpu;
  real* mt  = Sys->mass_cpu;
  real* vxplanet = Sys->vx_cpu;
  real* vyplanet = Sys->vy_cpu;
  boolean* Feel  = Sys->FeelDisk;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX-1; 
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
//<\EXTERNAL>

//<INTERNAL>
  real vslocal;
  real RBondi;
  real Raij;
  float Trij;
  real facc = 1;
  real facc2 = 1;
  float planetdistance;
  real rroche;
  real Hij;
  real smoothing;
  real dist;
  real vxcel;
  real vycel;
  real vtcel;
  real vrcel;
  real rij;
  real Haij;
  real c1;
  real c2;
  real c3;
  real c4;
  real c5;
  real c6;
  float c7;
  real deltaVrel;
  real c8;
  int i;
  int j;
  int k;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real Szj(Ny+2*NGHY);
// real Szk(Nz+2*NGHZ);
// real GAMMA(1);
// real OMEGAFRAME(1);
// boolean AccretionRunaway(1);
// real ratio_KH_max(1);
// real ASPECTRATIO(1);
// real ROCHESMOOTHING(1);
// real FLARINGINDEX(1);
// real THICKNESSSMOOTHING(1);
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
           planetdistance = sqrt(x[kp]*x[kp]+y[kp]*y[kp]+z[kp]*z[kp]);
           rroche = planetdistance*pow(mt[kp]/(3.0*MSTAR),(1.0/3.0));
           Hij = ASPECTRATIO*pow(planetdistance/R0,FLARINGINDEX)*planetdistance;

	   if (ROCHESMOOTHING != 0)
	     smoothing = rroche*ROCHESMOOTHING;
	   else
	     smoothing = Hij*THICKNESSSMOOTHING;

	   smoothing*=smoothing;

           #ifdef ISOTHERMAL
              vslocal = cs[l];
           #endif
           #ifdef ADIABATIC
              vslocal = sqrt(GAMMA*(GAMMA - 1.0)*cs[l]/dens[l]); 
           #endif
           #ifdef POLYTROPIC
              vslocal = sqrt(cs[l]*GAMMA*pow(dens[l],GAMMA-1.0));
           #endif

            #ifdef CARTESIAN  
               vxcel = Vt[l];
               vycel = Vr[l];       
            #endif

            #ifdef CYLINDRICAL 
               vtcel = 0.5*(Vt[l] + Vt[lxp]) + ymed(j)*OMEGAFRAME;
               vrcel = 0.5*(Vr[l] + Vr[lyp]);
               vxcel = (vrcel*XC - vtcel*YC)/ymed(j);
               vycel = (vrcel*YC + vtcel*XC)/ymed(j);
            #endif

            #ifdef SPHERICAL
               //working
            #endif

            c1 = vxcel - vxplanet[kp];
            c2 = vycel - vyplanet[kp];
            c3 = pow(c1,2.0);
            c4 = pow(c2,2.0);
            deltaVrel = sqrt(c3 + c4);
            c5 = pow(vslocal,2.0);
            c6 = pow(deltaVrel,2.0);
            c7 = c6+c5;
            RBondi = (2.0*G/c7)*mt[kp];

            Raij = MIN(RBondi,rroche);
 
	    dist =  ((XC-x[kp])*(XC-x[kp])+
		    (YC-y[kp])*(YC-y[kp])+
		    (ZC-z[kp])*(ZC-z[kp]));
            dist = sqrt(dist);

            if (dist < Raij)
               {
               #ifndef Z 
                  rij = sqrt(dist*dist + smoothing);
               #else
                  rij = sqrt(dist*dist);
               #endif
               c8 = rij*rij*rij;
               Trij = 2.0*PI*sqrt(c8/(G*mt[kp]));
               facc = (1/Trij)*dt;
               Haij = sqrt(Raij*Raij - dist*dist);

               #ifndef Z
                  if (Haij < Hij) facc2 = Haij/Hij;
                  else facc2 = 1.0;
               #else
                  facc2 = 1.0;
               #endif

               massacc[l] = facc*facc2*SurfZ(j,k)*dens[l];

               if (AccretionRunaway == YES)
                  {
                  if (dens[l] > 1.0e-10) 
                     {
                     dens[l] = dens[l] - (massacc[l]/SurfZ(j,k));
                     }
                  else
                     {
                     dens[l] = dens[l];
                     }
                  }
               else
                  {
                  if (dens[l] > 1.0e-10)
                     {
                     dens[l] = dens[l] - ratio_KH_max*(massacc[l]/SurfZ(j,k));
                     }
                  else  
                     {
                     dens[l] = dens[l];
                     }
                  }
               }
            else massacc[l] = 0.0;
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
  real* mass  = Sys->mass_cpu;
  real* vxp = Sys->vx_cpu;
  real* vyp = Sys->vy_cpu;
  real* accgas = Sys->accgas_cpu;
  real* massenv  = Sys->massenv_cpu;

  real* x1  = Sys->x_cpu;
  real* y1  = Sys->y_cpu;
  real* z1  = Sys->z_cpu;

  real deltaM;
  real temp = 0;
  real b = 8.0;
  real c = 2.5;
  real year;
  real tauKH;
  real dMdtKH;
  real dMdtmax;
  real dPxPlanet = 0.0;
  real dPyPlanet = 0.0;
  real PxPlanet = 0.0;
  real PyPlanet = 0.0;

  #ifdef MKS
     year = 31536000.0;
  #endif
  #ifdef CGS
     year = 31536000.0;
  #endif
  #if !(defined(MKS) || defined (CGS))
     year = 31536000.0/(sqrt(R0_MKS*R0_MKS*R0_MKS/G_MKS/MSTAR_MKS));
  #endif

  deltaM  = reduction_full_SUM(DeltaM,NGHY,Ny+NGHY,NGHZ,Nz+NGHZ);

  #ifndef __GPU
    MPI_Allreduce(&deltaM,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   
    deltaM = temp;
  #endif

  tauKH = pow(10.0,b)*pow(mass[kp]/MEARTH,-c);
  tauKH = tauKH*year;
  dMdtKH = mass[kp]/tauKH;
  dMdtmax = deltaM/dt;

  if (dMdtKH <= dMdtmax)
     {
      accgas[kp] = dMdtKH;
      ratio_KH_max = dMdtKH/dMdtmax;
      AccretionRunaway = NO;
      }
  else
      {
      accgas[kp] = deltaM/dt; 
      ratio_KH_max = 1.0;
      AccretionRunaway = YES;
      }

  massenv[kp] = massenv[kp] + accgas[kp]*dt;
  mass[kp] = mass[kp] + accgas[kp]*dt;
//<\LAST_BLOCK>
}
