/* This file stores the global variables */
// Avoid redefine
#ifndef VAR_H
#define VAR_H
// Pointers:
extern double **Points;                          // N x 3 array;
extern int **Surfaces;                           // N x 3 array; integer point to the vertices.
extern double **SurfaceNorms;                    // N x 3 array;
extern double *SurfaceAreas;                     // N x 1 array;
extern double *Specularity;                      // N x 1 array;
extern int *InSurfaces;                          // k x 1 array; integer point to the surfaces.
extern int *OutSurfaces;                         // k x 1 array; integer point to the surfaces.
extern double *Displacement;                     // N x 1 array; for calculating distance.
extern double *dist2surf;                        // N x 1 array: distance to surface.
extern int *SurfaceMarker;                       // N x 1 array; mark the surface types.
extern double *CumuAreasL, *CumuAreasR;          // N x 1 array; Surface Cumulative Areas for choosing the starting surface. L for left; R for right.

// Final results
extern double effectiveTC;                       // Effective Thermal Conductivity.    Need to intialize once.

// Rand Seed
extern long RandSeed[6];              // Random number generator seed;

// Calculation Parameters:
extern int NumPoints;                            // Number of Points;
extern int NumSurfaces;                          // Number of Surfaces;
extern int NumPeriods;                           // Simulation period number;
extern int NumInSurfaces;                        // Number of Surfaces in inflow.
extern int NumOutSurfaces;                       // Number of Surfaces in outflow.
extern long long NumParticles;                         // Number of Particles in simulation.re-intialize for different bulkMFP
extern int NumTrans;                             // Number of transmitted Particles.  re-intialize for different bulkMFP
extern int NumReflect;                           // Number of reflected Particles.    re-intialize for different bulkMFP
extern int SBinPerCell;                          // NUmber of Spatial bins per unit cell.
extern double TotalLength;                      // Input the sample length.Thus the inflow and outflow surfaces are parallel.
extern double UnitCellLength;                     // Input the unit cell length.
extern double Spec;                              // Input the specularity of the outer surfaces;
extern double Eeff;

// Structure;
typedef struct particle{
	double x;                             // location @ x;
	double y;                             // location @ y;
	double z;                             // relative location @ z;
	double realz;                         // the realistic z location.
	double tx;                            // direction @ x;
	double ty;                            // direction @ y;
	double tz;                            // direction @ z;
	int Surface;                          // collision surface.
	int NumMode;                          // the mode of the particle.
	int Pnum;                             // Currently in which unit cell. [0,NumPeriod);
	int sign;                             // the sign is determined by the sign of (T - Tref).
}particle;

// Random State
typedef struct{
	double x10, x11, x12;
	double x20, x21, x22;
}state_t;

// Phonon Property
extern int NumBins;                             // Number of frequency bins.
extern double *freq;                            // Unit for frequency: rad/s: frequency at the center of the bin.
extern double *DOS;                             // Unit: none
extern double *Vg;                              // Unit: m/s
extern double *dfreq;                           // Unit: rad/s
extern double *RTime;                           // Unit: s;
extern double *RRTime;                           // Unit: s;
extern int *Polar;                             // Phonon Branch: 1 for LA; 2 for TA;
extern double *de_dT;                           // for cumulative function calculation

// Thermal Property
extern double *Qz;                             // Nt * Nz Matrix: to record the contribution of each particle.
extern double *T;                              // Nt * Nz Matrix: to record the contribution of each particle.
extern double Tleft;                            // Left boundary temperature.
extern double Tright;                           // Right boundary temperature.
extern double Tinit;                            // Initialized Temperature in the center region, which is expected to be the average of Tleft and Tright, the same as Tref.
extern double Tref;                             // Referenced Temperature for the deviational algorithm.
extern double C;                                // heat capacity

// Cumulative functions
extern double *cumu_b;                            // for the bulk region.
extern double *cumu_v;                            // for heat flux at the boundary.
extern double *cumu_col;                          // collision term

// deviational energy
extern double enrgInit;                           // deviational energy for intialization in the bulk region.
extern double enrgLeft;                           // deviational energy for the Left Contact.
extern double enrgRight;                          // deviational energy for the Right Contact.
extern double enrgTot;                            // Total deviational energy.

// Geometrical Parameters. Firstly apply to a simple thin film case.
extern double Volume;                            // Total Volume of the system.
extern double *Vs;                               // Volume for the spatial cell. Only consider one direction discretization.
extern int Ztot;                                 // Total Number of Spatial Cells.
extern double *Zmins, *Zmaxs,*Zcenters;           // Help with finding the spatial cell. The range of the Bins along z direction is [Zmin,Zmax];
extern double AreaLeft;                          // Left Surface Area.
extern double AreaRight;                         // Right Surface Area.
extern double Xlen;                              // Total Length at X direction.
extern double Ylen;                              // Total Length at Y direction.
extern double Zlen;                              // Total Length at Z direction.

// Fixed Constant;
extern const double INF;                       // Represent the infinity.
extern const double PI;                        // Constant PI.
extern const double kB;                        // Boltzmann constant
extern const double hbar;                      // reduced Planck constant.
extern const double uconvert;            // convert length scale from nm to m.


// timing
extern double elapsedtime;                // for timing.
#endif
