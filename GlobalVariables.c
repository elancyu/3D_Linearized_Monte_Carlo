/* This file stores the global variables */

// Pointers:
double **Points;                          // N x 3 array;
int **Surfaces;                           // N x 3 array; integer point to the vertices.
double **SurfaceNorms;                    // N x 3 array;
double *SurfaceAreas;                     // N x 1 array;
double *Specularity;                      // N x 1 array;
int *InSurfaces;                          // k x 1 array; integer point to the surfaces.
int *OutSurfaces;                         // k x 1 array; integer point to the surfaces.
double *Displacement;                     // N x 1 array; for calculating distance.
double *dist2surf;                        // N x 1 array: distance to surface.
int *SurfaceMarker;                       // N x 1 array; mark the surface types.
double *CumuAreasL, *CumuAreasR;          // N x 1 array; Surface Cumulative Areas for choosing the starting surface. L for Left; R for Right.

// Final results
double effectiveTC;                       // Effective Thermal Conductivity.    Need to intialize once.

// RandSeed
long RandSeed[6];              // Random Seed

// Calculation Parameters:
int NumPoints;                            // Number of Points;
int NumSurfaces;                          // Number of Surfaces;
int NumPeriods;                           // Simulation period number;
int NumInSurfaces;                        // Number of Surfaces in inflow.
int NumOutSurfaces;                       // Number of Surfaces in outflow.
long long NumParticles;                         // Number of Particles in simulation.re-intialize for different bulkMFP
int NumTrans;                             // Number of transmitted Particles.  re-intialize for different bulkMFP
int NumReflect;                           // Number of reflected Particles.    re-intialize for different bulkMFP
int SBinPerCell;                         // Number of Spatial Bins Per Unit Cell.
double TotalLength;                      // Thus the inflow and outflow surfaces are parallel. Total Len = Unit Len * NumPeriods
double UnitCellLength;                   // Input the unit cell length.
double Spec=0;                            // Input specularity.
double Eeff;                              // Effective Energy for the simulation phonon bundle.

// Phonon Property
int NumBins;                             // Number of frequency bins.
double *freq;                            // Unit for frequency: rad/s: frequency at the center of the bin.
double *DOS;                             // Unit: none
double *Vg;                              // Unit: m/s
double *dfreq;                           // Unit: rad/s
double *RTime;                           // Unit: s;
double *RRTime;                           // Unit: s;
int *Polar;                              // The branch; 1 for LA; 2 for TA.
double *de_dT;                           // for cumulative function calculation.

// Thermal Property
double *Qz;                             // Nt * Nz Matrix: to record the contribution of each particle.
double *T;                              // Nt * Nz Matrix: to record the contribution of each particle.
double Tleft;                            // Left boundary temperature.
double Tright;                           // Right boundary temperature.
double Tinit;                            // Initialized Temperature in the center region, which is expected to be the average of Tleft and Tright, the same as Tref.
double Tref;                             // Referenced Temperature for the deviational algorithm.
double C;                                // heat capacity

// Cumulative functions
double *cumu_b;                            // for the bulk region.
double *cumu_v;                            // for heat flux at the boundary.
double *cumu_col;                          // collision term

// deviational energy
double enrgInit;                           // deviational energy for intialization in the bulk region.
double enrgLeft;                           // deviational energy for the Left Contact.
double enrgRight;                          // deviational energy for the Right Contact.
double enrgTot;                            // Total deviational energy.

// Geometrical Parameters. Firstly apply to a simple thin film case.
double Volume;                            // Total Volume of the system.
double *Vs;                               // Volume for the spatial cell. Only consider one direction discretization.
int Ztot;                                 // Total Number of Spatial Cells.
double *Zmins, *Zmaxs, *Zcenters;         // Help with finding the spatial cell. The range of the Bins along z direction is [Zmin,Zmax];
double AreaLeft;                          // Left Surface Area.
double AreaRight;                         // Right Surface Area.
double Xlen;                              // Total Length at X direction.
double Ylen;                              // Total Length at Y direction.
double Zlen;                              // Total Length at Z direction.

// Fixed Constant;
double const INF=1e20;                         // Represent the infinity.
double const PI=3.141592653;                   // Constant PI from UIUC: truncated at 9th digit.
double const kB = 1.38e-23;                    // Boltzmann Constant, unit: J/K
double const hbar = 1.05457180014e-34;         // reduced Planck constant, unit: SI unit.
double const uconvert = 1e-9;            // convert nm to m: multiply this constant.

double elapsedtime;
