/* This file lists all the functions */

#include"GlobalVariables.h"

// read in the geometry and build the list.
// input: the input file name/address.
// output: Check file to ensure correctness of inputs.
int GeoBuilder(char *fname);

// read in the geometry and build the list.
// input: the input file name/address.
// output: Check file to ensure correctness of inputs.
void SpatialCellBuilder(char *fname);

// allocate memory for Statistics.
// input: None
// output: None.
void StatBuilder(void);

// read in the material property.
// input: the input file name/address.
// output: Check file to ensure correctness of inputs.
void PhononBuilder(char *fname);

// free all the allocated memory.
// input: argc.
// output: none.
void FreeMem(int argc);

// calculate the thermal conductivity.
// input:The bulk MFP spectra.
// output: effective thermal conductivity to file.
void PerformMC();

// Choose the Mode Number for the particle.
// input: 1.type;
// output: reflected or transmitted?
int ChooseMode(int type);

// Locate which unit cell the location is in.
// input: the location
// output: Cell Num
int LocateCellNum(double z);

// build cumulative function for choosing starting surface
// input: none; prerequisite: geobuilder;
// output: none; Area CF is computed
void buildCumuAreas(void);

// calculate mean heat flux and the effective thermal conductivity.
// input: none;
// output: none;
void Logger(void);

// RNG Initializer
// input: randseed
// output none
void InitRand(unsigned long long RandSeed);

// Random number Generator
// input: none
// output: a uniformly distributed random number in [0,1)
double RandR();

// tell if the point is inside an triangle in 2D
// input: Surface Number, point position.
// output: true or false: 1 or 0;
int IsInsideTriangle2D(double px, double py, int SurfaceNum);

// Calculate the area of triangle.
// input: three vertices
// output: area
double TriArea(double p1[3], double p2[3], double p3[3]);
