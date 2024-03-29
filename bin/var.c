#define __LOCAL
#include "../src/fargo3d.h"
#undef __LOCAL

void InitVariables() {
  init_var("RHOICE", (char*)&RHOICE, REAL, NO, "2.37e8");
  init_var("B", (char*)&B, REAL, NO, "10.0");
  init_var("SIGMASLOPE", (char*)&SIGMASLOPE, REAL, NO, "0.5");
  init_var("FDG", (char*)&FDG, REAL, NO, "0.03");
  init_var("THICKNESSSMOOTHING", (char*)&THICKNESSSMOOTHING, REAL, NO, "0.6");
  init_var("FRI", (char*)&FRI, REAL, NO, "0.25");
  init_var("YMAX", (char*)&YMAX, REAL, NO, "2.0");
  init_var("YMIN", (char*)&YMIN, REAL, NO, "0.1");
  init_var("DT", (char*)&DT, REAL, NO, "0.314159265359");
  init_var("ROCHESMOOTHING", (char*)&ROCHESMOOTHING, REAL, NO, "0.0");
  init_var("XMIN", (char*)&XMIN, REAL, NO, "-3.14159265358979323844");
  init_var("CD", (char*)&CD, REAL, NO, "1.0");
  init_var("XMAX", (char*)&XMAX, REAL, NO, "3.14159265358979323844");
  init_var("RHOM", (char*)&RHOM, REAL, NO, "7.57e8");
  init_var("RM", (char*)&RM, REAL, NO, "1.28549e-10");
  init_var("ECCENTRICITY", (char*)&ECCENTRICITY, REAL, NO, "0.0");
  init_var("ALPHA", (char*)&ALPHA, REAL, NO, "3.00e-3");
  init_var("ASPECTRATIO", (char*)&ASPECTRATIO, REAL, NO, "0.05");
  init_var("RICE", (char*)&RICE, REAL, NO, "0.5769");
  init_var("SIGMA0", (char*)&SIGMA0, REAL, NO, "2.66e-3");
  init_var("RHOROCKY", (char*)&RHOROCKY, REAL, NO, "7.57e8");
  init_var("NX", (char*)&NX, INT, NO, "583");
  init_var("NY", (char*)&NY, INT, NO, "251");
  init_var("NINTERM", (char*)&NINTERM, INT, NO, "30");
  init_var("NTOT", (char*)&NTOT, INT, NO, "3000000");
  init_var("PLANETCONFIG", (char*)&PLANETCONFIG, STRING, NO, "planets/earth.cfg");
  init_var("SETUP", (char*)&SETUP, STRING, NO, "solids_eqY_size1	");
  init_var("FRAME", (char*)&FRAME, STRING, NO, "C");
  init_var("OUTPUTDIR", (char*)&OUTPUTDIR, STRING, NO, "@outputs/solids_eqY_size1");
  init_var("INDIRECTTERM", (char*)&INDIRECTTERM, BOOL, NO, "1");
  init_var("EXCLUDEHILL", (char*)&EXCLUDEHILL, BOOL, NO, "1");
  init_var("ZMAX", (char*)&ZMAX, REAL, NO, "1.0");
  init_var("SEMIMAJORAXIS", (char*)&SEMIMAJORAXIS, REAL, NO, "0.0");
  init_var("ZMIN", (char*)&ZMIN, REAL, NO, "0.0");
  init_var("ETA", (char*)&ETA, REAL, NO, "0.0");
  init_var("GAMMA", (char*)&GAMMA, REAL, NO, "1.66666667");
  init_var("NOISE", (char*)&NOISE, REAL, NO, "0.0");
  init_var("MASSTAPER", (char*)&MASSTAPER, REAL, NO, "0.0");
  init_var("PLANETMASS", (char*)&PLANETMASS, REAL, NO, "0.0");
  init_var("VMAX", (char*)&VMAX, REAL, NO, "1.0");
  init_var("INNERRESONANCE", (char*)&INNERRESONANCE, REAL, NO, "0.5");
  init_var("RELEASEDATE", (char*)&RELEASEDATE, REAL, NO, "0.0");
  init_var("FLARINGINDEX", (char*)&FLARINGINDEX, REAL, NO, "0.0");
  init_var("VMIN", (char*)&VMIN, REAL, NO, "0.0");
  init_var("CS", (char*)&CS, REAL, NO, "1.0");
  init_var("OORTA", (char*)&OORTA, REAL, NO, "-0.75");
  init_var("INCLINATION", (char*)&INCLINATION, REAL, NO, "0.0");
  init_var("CFL", (char*)&CFL, REAL, NO, "0.44");
  init_var("ORBITALRADIUS", (char*)&ORBITALRADIUS, REAL, NO, "0.0");
  init_var("OMEGAFRAME", (char*)&OMEGAFRAME, REAL, NO, "0.0");
  init_var("NU", (char*)&NU, REAL, NO, "0.0");
  init_var("RELEASERADIUS", (char*)&RELEASERADIUS, REAL, NO, "0.0");
  init_var("OUTERRESONANCE", (char*)&OUTERRESONANCE, REAL, NO, "1.5");
  init_var("VERTICALDAMPING", (char*)&VERTICALDAMPING, REAL, NO, "0.0");
  init_var("NZ", (char*)&NZ, INT, NO, "1");
  init_var("NSNAP", (char*)&NSNAP, INT, NO, "0");
  init_var("REALTYPE", (char*)&REALTYPE, STRING, NO, "Standard  #float or double");
  init_var("PLOTLINE", (char*)&PLOTLINE, STRING, NO, "field[:,:,0]");
  init_var("SPACING", (char*)&SPACING, STRING, NO, "default");
  init_var("COORDINATES", (char*)&COORDINATES, STRING, NO, "standard");
  init_var("FIELD", (char*)&FIELD, STRING, NO, "gasdens");
  init_var("CMAP", (char*)&CMAP, STRING, NO, "cubehelix");
  init_var("FUNCARCHFILE", (char*)&FUNCARCHFILE, STRING, NO, "std/func_arch.cfg");
  init_var("ASPECT", (char*)&ASPECT, STRING, NO, "auto");
  init_var("WRITEENERGY", (char*)&WRITEENERGY, BOOL, NO, "0");
  init_var("WRITEENERGYRAD", (char*)&WRITEENERGYRAD, BOOL, NO, "0");
  init_var("LOG", (char*)&LOG, BOOL, NO, "0");
  init_var("WRITEECC", (char*)&WRITEECC, BOOL, NO, "0");
  init_var("WRITEDIVERGENCE", (char*)&WRITEDIVERGENCE, BOOL, NO, "0");
  init_var("DISK", (char*)&DISK, BOOL, NO, "0");
  init_var("WRITETAU", (char*)&WRITETAU, BOOL, NO, "0");
  init_var("WRITEVX", (char*)&WRITEVX, BOOL, NO, "0");
  init_var("WRITEVY", (char*)&WRITEVY, BOOL, NO, "0");
  init_var("WRITEVZ", (char*)&WRITEVZ, BOOL, NO, "0");
  init_var("WRITEDENSITY", (char*)&WRITEDENSITY, BOOL, NO, "0");
  init_var("PLANETHEATING", (char*)&PLANETHEATING, BOOL, NO, "0");
  init_var("VTK", (char*)&VTK, BOOL, NO, "0");
  init_var("WRITEINC", (char*)&WRITEINC, BOOL, NO, "0");
  init_var("WRITESIGMAM", (char*)&WRITESIGMAM, BOOL, NO, "0");
  init_var("PERIODICY", (char*)&PERIODICY, BOOL, NO, "0");
  init_var("PERIODICZ", (char*)&PERIODICZ, BOOL, NO, "0");
  init_var("WRITEBX", (char*)&WRITEBX, BOOL, NO, "0");
  init_var("WRITEBY", (char*)&WRITEBY, BOOL, NO, "0");
  init_var("WRITEBZ", (char*)&WRITEBZ, BOOL, NO, "0");
  init_var("COLORBAR", (char*)&COLORBAR, BOOL, NO, "1");
  init_var("AUTOCOLOR", (char*)&AUTOCOLOR, BOOL, NO, "1");
}