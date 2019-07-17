
#include "fargo3d.h"

void rescale () {
RHOICE *= MSTAR/(R0*R0*R0);
DT *= sqrt(R0*R0*R0/G/MSTAR);
RHOM *= MSTAR/(R0*R0*R0);
RM *= R0;
RICE *= R0;
SIGMA0 *= MSTAR/(R0*R0);
RHOROCKY *= MSTAR/(R0*R0*R0);
}
