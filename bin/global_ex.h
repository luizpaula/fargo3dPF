 /* This file was created automatically during compilation
from global.h. Do not edit. See python script
"global.py" for details. */ 

extern int CPU_Rank;
extern int CPU_Number;
extern boolean CPU_Master;
extern boolean VxIsResidual;
extern boolean LogGrid;
extern boolean GuidingCenter;
extern boolean Corotating;
extern boolean Restart;
extern boolean Restart_Full;
extern boolean Stockholm;
extern boolean Merge;
extern boolean Merge_All;
extern boolean MonitorIntegral;
extern boolean TimeInfo;
extern boolean EverythingOnCPU;
extern boolean ForwardOneStep;
extern boolean Vtk2dat;
extern boolean Dat2vtk;
extern boolean PostRestart;
extern boolean OnlyInit;
extern boolean EarlyOutputRename;
extern boolean RedefineOptions;
extern boolean DeviceFileSpecified;
extern boolean StretchOldOutput;
extern boolean ThereArePlanets;
extern boolean ThereIsACentralBinary;
extern boolean AccretionRunaway;
extern real ratio_KH_max;
extern real PhysicalTimeInitial;
extern real PhysicalTime;
extern real XAxisRotationAngle;
extern char NewOutputdir[1024];
extern char DefaultOut[1024];
extern char DeviceFile[1024];
extern int TimeStep;
extern int FullArrayComms;
extern int ContourComms;
extern int DeviceManualSelection;
extern int ArrayNb;
extern int StretchNumber;
extern int BinaryStar1;
extern int BinaryStar2;
extern real InnerBorder;
extern real OuterBorder;
extern TimeProcess t_speedup_cpu;
extern TimeProcess t_speedup_gpu;
extern int t_speedup_count;
extern real time_speedup_cpu;
extern real time_speedup_gpu;
extern PlanetarySystem *Sys;
extern Point DiskOnPrimaryAcceleration;
extern Point IndirectTerm;
extern real step_time;
extern real localforce[12];
extern real globalforce[12];
extern Grid Gridd;
extern Field *Vx;
extern Field *Vy;
extern Field *Vz;
extern Field *Vx_temp;
extern Field *Vy_temp;
extern Field *Vz_temp;
extern Field *Slope;
extern Field *Mpx;
extern Field *Mpy;
extern Field *Mpz;
extern Field *Mmx;
extern Field *Mmy;
extern Field *Mmz;
extern Field *Pot;
extern Field *DivRho;
extern Field *DensStar;
extern Field *Qs;
extern Field *Density;
extern Field *Energy;
extern Field *Pressure;
extern Field *Sigmam;
extern Field *Ecc;
extern Field *Inc;
extern Field *Rhogas;
extern Field *Numol;
extern Field *Reymol;
extern Field *Dedt;
extern Field *Didt;
extern Field *DeltaM;
extern Field *Eccquad;
extern Field *Incquad;
extern Field *QL;
extern Field *QR;
extern Field *LapPPA;
extern Field2D *VxMed;
extern Field2D *Vxhy;
extern Field2D *Vxhyr;
extern Field2D *Vxhz;
extern Field2D *Vxhzr;
extern Field2D *Reduction2D;
extern Field2D *Eta_profile_xi;
extern Field2D *Eta_profile_xizi;
extern Field2D *Eta_profile_zi;
extern FieldInt2D *Nxhy;
extern FieldInt2D *Nxhz;
extern FieldInt2D *Nshift;
extern Field *Bx;
extern Field *By;
extern Field *Bz;
extern Field *B1_star;
extern Field *B2_star;
extern Field *V1_star;
extern Field *V2_star;
extern Field *Slope_b1;
extern Field *Slope_v1;
extern Field *Slope_b2;
extern Field *Slope_v2;
extern Field *Emfx;
extern Field *Emfy;
extern Field *Emfz;
extern Field *Divergence;
extern Field2D *Density0;
extern Field2D *Vx0;
extern Field2D *Vy0;
extern Field2D *Vz0;
extern Field2D *Energy0;
extern int Ncpu_x;
extern int Ncpu_y;
extern Buffer Bfd;
extern Buffer Bfu;
extern Buffer Bfl;
extern Buffer Bfr;
extern Buffer Bfcdl;
extern Buffer Bfcdr;
extern Buffer Bfcul;
extern Buffer Bfcur;
extern real Dx;
extern real *Xmin;
extern real *Ymin;
extern real *Zmin;
extern real *Xmed;
extern real *Ymed;
extern real *Zmed;
extern real *InvDiffXmed;
extern real *InvDiffYmed;
extern real *InvDiffZmed;
extern real *Sxj;
extern real *Sxk;
extern real *Syj;
extern real *Syk;
extern real *Szj;
extern real *Szk;
extern real *InvVj;
extern real shift_buffer[MAX1D];
extern real *Dx_d;
extern real *Xmin_d;
extern real *Ymin_d;
extern real *Zmin_d;
extern real *Sxj_d;
extern real *Sxk_d;
extern real *Syj_d;
extern real *Syk_d;
extern real *Szj_d;
extern real *Szk_d;
extern real *InvVj_d;
extern real shift_buffer_d[MAX1D];
extern int Nx;
extern int Ny;
extern int Nz;
extern int J;
extern int K;
extern int Y0;
extern int Z0;
extern int Stride_cpu;
extern int Stride_gpu;
extern int Pitch_cpu;
extern int Pitch_gpu;
extern int Pitch_Int_gpu;
extern int Pitch2D;
extern int Stride;
extern int ix;
extern int iy;
extern int ycells;
extern int zcells;
extern int y0cell;
extern int z0cell;
extern Field *ListOfGrids;
extern int Bounl;
extern int Bounr;
extern int Bounu;
extern int Bound;
extern real Xplanet;
extern real Yplanet;
extern real Zplanet;
extern real VXplanet;
extern real VYplanet;
extern real VZplanet;
extern real MplanetVirtual;
extern real SigmarmsVirtual;
extern real EccrmsVirtual;
extern real IncrmsVirtual;
extern real McoreVirtual;
extern real RcoreVirtual;
extern real RcapVirtual;
extern real AccretionsolidsVirtual;
extern real MenvVirtual;
extern real McritVirtual;
extern real AccretiongasVirtual;
extern MPI_Status fargostat;
extern real OMEGAFRAME0;
extern int Fscan;
extern long VtkPosition;
extern void (*Denshom)(int);
extern void (*Computerhogas)();
extern void (*Computenumol)();
extern void (*Computegasdragecc)();
extern void (*Computegravecc)();
extern void (*Computegasdraginc)();
extern void (*Computegravinc)();
extern void (*Computereymol)();
extern void (*Updateecc)(real);
extern void (*Updateinc)(real);
extern void (*Computeaccretionsolids)(real,int);
extern void (*Computeaccretionsolidsequilibrium)(real,int);
extern void (*Accretiongas)(real,int);
extern void (*ComputePressureFieldIso)();
extern void (*ComputePressureFieldAd)();
extern void (*ComputePressureFieldPoly)();
extern void (*SubStep1_x)(real);
extern void (*SubStep1_y)(real);
extern void (*SubStep1_z)(real);
extern void (*SubStep2_a)(real);
extern void (*SubStep2_b)(real);
extern void (*SubStep3)(real);
extern void (*DivideByRho)(Field*);
extern void (*VanLeerX_a)(Field*);
extern void (*VanLeerX_b)(real,Field*,Field*,Field*);
extern void (*VanLeerY_a)(Field*);
extern void (*VanLeerY_b)(real,Field*,Field*);
extern void (*VanLeerZ_a)(Field*);
extern void (*VanLeerZ_b)(real,Field*,Field*);
extern void (*momenta_x)();
extern void (*momenta_y)();
extern void (*momenta_z)();
extern void (*reduction_SUM)(Field*,int,int,int,int);
extern void (*reduction_MIN)(Field*,int,int,int,int);
extern void (*UpdateX)(real,Field*,Field*,Field*);
extern void (*UpdateY)(real,Field*,Field*);
extern void (*UpdateZ)(real,Field*,Field*);
extern void (*UpdateDensityX)(real,Field*,Field*);
extern void (*UpdateDensityY)(real,Field*);
extern void (*UpdateDensityZ)(real,Field*);
extern void (*NewVelocity_x)();
extern void (*NewVelocity_y)();
extern void (*NewVelocity_z)();
extern void (*AdvectSHIFT)(Field*,FieldInt2D*);
extern void (*ComputeResidual)(real);
extern void (*ChangeFrame)(int,Field*,Field2D*);
extern void (*Potential)();
extern void (*CorrectVtheta)(real);
extern void (*cfl)();
extern void (*_ComputeForce)(real,real,real,real,real);
extern void (*copy_velocities)(int);
extern void (*VanLeerX_PPA_a)(Field*);
extern void (*VanLeerX_PPA_b)(Field*);
extern void (*VanLeerX_PPA_steep)(Field*);
extern void (*VanLeerX_PPA_c)(Field*);
extern void (*VanLeerX_PPA_d)(real,Field*,Field*,Field*);
extern void (*VanLeerX_PPA_d_2d)(real,Field*,Field*,Field2D*);
extern void (*mon_dens)();
extern void (*mon_momx)();
extern void (*mon_momy)();
extern void (*mon_momz)();
extern void (*mon_torq)();
extern void (*mon_reynolds)();
extern void (*mon_maxwell)();
extern void (*mon_bxflux)();
extern void (*comm)();
extern void (*ComputeSlopes)(int,int,int,Field*,Field*);
extern void (*_ComputeStar)(real,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*);
extern void (*_ComputeEmf)(real,int,int,int,int,int,int,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*);
extern void (*_UpdateMagneticField)(real,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*);
extern void (*_LorentzForce)(real,int,int,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*,Field*,Field*);
extern void (*_Resist)(int,int,int,int,int,int,int,int,int,Field*,Field*,Field*,Field2D*);
extern void (*EMF_Upstream_Integrate)(real);
extern void (*StockholmBoundary)(real);
extern void (*visctensor_cart)();
extern void (*addviscosity_cart)(real);
extern void (*visctensor_cyl)();
extern void (*addviscosity_cyl)(real);
extern void (*visctensor_sph)();
extern void (*addviscosity_sph)(real);
extern void (*boundary_ymin)();
extern void (*boundary_zmin)();
extern void (*boundary_ymax)();
extern void (*boundary_zmax)();
extern void (*Fill_GhostsX)();
extern void (*CheckMuteY)();
extern void (*CheckMuteZ)();
extern void (*SetupHook1)();
extern void (*__WriteField)();
extern void (*__Restart)(Field*,int);
