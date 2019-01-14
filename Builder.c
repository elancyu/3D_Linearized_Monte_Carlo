/* This C will read in geometry then store it */
// input: geometry.in
// output: none.

/* description for input file: the surface is consist of points
and points are in triplet of real number. First block is points
Second block is surface. and Final block is periodicity.     */



#include"GlobalVariables.h"
#include"Functions.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

int GeoBuilder(char* fname)
{
	FILE *fpgeo;                // geometry input file.
	FILE *fpchecker;            // geometry checker.
	int i, j, surfnum;          // loop integer variable.
	int pn, sn;                 // point number, surface number
	double a, b, c;             // for Norm and Area Calculation.
	int p1, p2, p3;             // vertices for triangle.
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double temp;
	// open the input file with the file name
	if ((fpgeo = fopen(fname, "a+")) == NULL)
	{
		printf("cannot open geometry file\n");
		return 0;
	}
	// Read random seed
	fscanf(fpgeo,"%d",&RandSeed);
	InitRand(RandSeed);
	// read in points order(x,y,z);allocmem;readin points;Points is a 2D array.
	fscanf(fpgeo,"%d",&NumPoints);
	
	Points = (double**)malloc(NumPoints*sizeof(double*));
	for (i =  0; i < NumPoints; i++)
		Points[i] = (double*)malloc(3*sizeof(double));
	for (i = 0; i < NumPoints; i++)
	{
		fscanf(fpgeo,"%d",&pn);
		for (j = 0; j < 3; j++)
			fscanf(fpgeo, "%lf",&Points[pn][j]);
	}
	
	// read in surfaces order(p1,p2,p3,spec)
	fscanf(fpgeo,"%d",&NumSurfaces);
	Surfaces = (int**)malloc(NumSurfaces*sizeof(int*));
	SurfaceMarker = (int*)malloc(NumSurfaces*sizeof(int));
	SurfaceNorms = (double**)malloc(NumSurfaces*sizeof(double*));
	SurfaceAreas = (double*)malloc(NumSurfaces*sizeof(double));
	Displacement = (double*)malloc(NumSurfaces*sizeof(double));
	dist2surf = (double*)malloc(NumSurfaces*sizeof(double));
	Specularity = (double*)malloc(NumSurfaces*sizeof(double));
	for (i = 0; i < NumSurfaces; i++)
	{
		Surfaces[i] = (int*)malloc(3*sizeof(int));
		SurfaceNorms[i] = (double*)malloc(3*sizeof(double));
		SurfaceAreas[i] = 1.0;
	}
	for (i = 0; i < NumSurfaces; i++)
	{
		fscanf(fpgeo,"%d",&sn);
		for (j = 0; j < 3; j++)
		{
			fscanf(fpgeo,"%d",&Surfaces[sn][j]);
			SurfaceNorms[sn][j] = 0.0;           // set to a safe value for the norms.
		}
		fscanf(fpgeo,"%d",&SurfaceMarker[sn]);   // SurfaceMarker: mark the surface type.
		if (SurfaceMarker[sn]==0)
			Specularity[sn] = Spec;              // input outer surface specularity.
		else
			Specularity[sn] = 1.0;               // For periodic or in/out surfaces.
		// SurfaceMarker: 0: normal surfaces; 1: inflow surfaces; 2:outflow surfaces.
	}
	
	// Force define the inflow and outflow surfaces.
	fscanf(fpgeo,"%d",&NumInSurfaces);
	InSurfaces = (int*)malloc(NumInSurfaces*sizeof(int));
	for (i = 0; i < NumInSurfaces; i++)
		fscanf(fpgeo, "%d", &InSurfaces[i]);
	fscanf(fpgeo, "%d",&NumOutSurfaces);
	OutSurfaces = (int*)malloc(NumOutSurfaces*sizeof(int));
	for (i = 0; i < NumOutSurfaces; i++)
		fscanf(fpgeo, "%d", &OutSurfaces[i]);
	// read in period num
	fscanf(fpgeo,"%d", &NumPeriods);
	// read in unitcell length (in nm)
	fscanf(fpgeo, "%lf", &UnitCellLength);
	// Initialize the total length
	TotalLength = NumPeriods*UnitCellLength;
	// read in the total number of particle in simulation.
	fscanf(fpgeo, "%d", &NumParticles);
	// Readin the temperatures
	fscanf(fpgeo, "%lf", &Tleft);
	fscanf(fpgeo, "%lf", &Tright);
	fscanf(fpgeo, "%lf", &Tinit);
	fscanf(fpgeo, "%lf", &Tref);
	
	// First Time calculate the Surfaces Normals,Surface Area and Surface_Displacement.
	for (i = 0; i < NumSurfaces; i++)
	{
		// S = 1/2 * cross(a,b). And Norm = cross(a,b)/module.
		// get points and data
		p1 = Surfaces[i][0];
		p2 = Surfaces[i][1];
		p3 = Surfaces[i][2];
		x1 = Points[p1][0];
		y1 = Points[p1][1];
		z1 = Points[p1][2];
		x2 = Points[p2][0];
		y2 = Points[p2][1];
		z2 = Points[p2][2];
		x3 = Points[p3][0];
		y3 = Points[p3][1];
		z3 = Points[p3][2];
		// calculate the cross
		a = (y1-y2)*(z1-z3)-(z1-z2)*(y1-y3);
		b = (z1-z2)*(x1-x3)-(x1-x2)*(z1-z3);
		c = (x1-x2)*(y1-y3)-(y1-y2)*(x1-x3);
		temp = sqrt(a*a+b*b+c*c);
		
		// Nomalized Normal.
		SurfaceNorms[i][0] = a/temp;
		SurfaceNorms[i][1] = b/temp;
		SurfaceNorms[i][2] = c/temp;
		SurfaceAreas[i] = temp/2;

		// Displacement
		Displacement[i] = -(x1*SurfaceNorms[i][0]+y1*SurfaceNorms[i][1]+z1*SurfaceNorms[i][2]);
	}
	// close FILE pointer;
	fclose(fpgeo);
	// Calculate the left and right surface areas.
	AreaLeft = 0; AreaRight = 0;
	for (i = 0; i < NumInSurfaces; i++)
	{
		surfnum = InSurfaces[i];
		AreaLeft += SurfaceAreas[surfnum];
	}
	for (i = 0; i < NumOutSurfaces; i++)
	{
		surfnum = OutSurfaces[i];
		AreaRight += SurfaceAreas[surfnum];
	}
	// Output the data for checker;
	if ((fpchecker=fopen("./CheckGeo.out","a+"))==NULL)
		return -1;                 // output some error info
	fprintf(fpchecker,"-----------Output for Comparison----------\n");
	for (i = 0; i < 6; i++)
		fprintf(fpchecker,"RandSeed:%d ",RandSeed[i]);
	fprintf(fpchecker,"\nNumPoints:%d  Unit: m\n",NumPoints);
	fprintf(fpchecker,"Rank  | P1 | P2 | P3 | Specularity | Mark\n");
	for (i = 0; i<NumPoints; i++)
		fprintf(fpchecker, "%d %.6e %.6e %.6e\n",i, Points[i][0],Points[i][1],Points[i][2]);
	fprintf(fpchecker,"\n");
	fprintf(fpchecker,"NumSurfaces:%d\n",NumSurfaces);
	for (i = 0; i < NumSurfaces; i++)
		fprintf(fpchecker,"%d  %d %d %d  %.3lf  %d\n",i, Surfaces[i][0],Surfaces[i][1],Surfaces[i][2],Specularity[i],SurfaceMarker[i]);
	
	fprintf(fpchecker,"NumLeftSurfaces:%d\n",NumInSurfaces);
	for (i = 0 ; i < NumInSurfaces; i++)
		fprintf(fpchecker, "%d\n", InSurfaces[i]);
	fprintf(fpchecker,"NumRightSurfaces:%d\n",NumOutSurfaces);
	for (i = 0; i < NumOutSurfaces; i++)
		fprintf(fpchecker, "%d\n",OutSurfaces[i]);
	fprintf(fpchecker,"Number of Periods:%d\n",NumPeriods);
	fprintf(fpchecker,"Unit Cell Length:%.3e m\n", UnitCellLength);
	fprintf(fpchecker,"Sample Length:%.3e m\n",TotalLength);
	fprintf(fpchecker,"Particle Number:%d\n", NumParticles);
	fprintf(fpchecker,"Left Temperature:%.3lf K\n", Tleft);
	fprintf(fpchecker,"Right Temperature:%.3lf K\n", Tright);
	fprintf(fpchecker,"Initial Temperature:%.3lf K\n", Tinit);
	fprintf(fpchecker,"Referenced Temperature:%.3lf K\n", Tref);
	fprintf(fpchecker,"-------------------END--------------------\n");
	fclose(fpchecker);
	// Build the Area Cumulative Function for inflow surfaces.
	buildCumuAreas();

	printf("Building Geometry Succeed!\n");
	// Some other operations? Output the codes to generate the 3D plot?
	/* This need only be done once. End of function*/
}

void PhononBuilder(char *fname)
{
	FILE *fp, *fpchecker, *fpcumu;
	int i;
	double x;           // intermediate variable
	fp = fopen(fname,"a+");
	// 1st Number: Number of Bins
	fscanf(fp,"%d",&NumBins);
	freq = (double*)malloc(NumBins*sizeof(double));
	dfreq = (double*)malloc(NumBins*sizeof(double));
	DOS = (double*)malloc(NumBins*sizeof(double));
	Vg = (double*)malloc(NumBins*sizeof(double));
	RTime = (double*)malloc(NumBins*sizeof(double));
	de_dT = (double*)malloc(NumBins * sizeof(double));
	Polar = (int*)malloc(NumBins*sizeof(int));
	cumu_b = (double*)malloc(NumBins*sizeof(double));
	cumu_v = (double*)malloc(NumBins*sizeof(double));
	cumu_col = (double*)malloc(NumBins*sizeof(double));
	for (i = 0; i < NumBins; i++)
	{
		fscanf(fp, "%lf", &freq[i]);
		fscanf(fp, "%lf", &DOS[i]);
		fscanf(fp, "%lf", &Vg[i]);
		fscanf(fp, "%lf", &dfreq[i]);
		fscanf(fp, "%lf", &RTime[i]);      // temporarily use the Relaxation time at 300K. May turn to the formula instead.
		fscanf(fp, "%d", &Polar[i]);       // Inidcator of Polarization.
		cumu_b[i] = 0.0;
		cumu_v[i] = 0.0;
		cumu_col[i] = 0.0;
	}
	fclose(fp);
	// Build the cumulative function for the particle generation.
	for (i = 0; i < NumBins; i++)
	{
		x = hbar*freq[i]/(kB*Tref);
		de_dT[i] = x*x*kB*exp(x)/((exp(x)-1)*(exp(x) - 1));
		// for cumulative function calculation.
	}
	// cumulative functions calculation.
	i = 0;
	cumu_b[0] = DOS[i]*dfreq[i]*de_dT[i];
	cumu_v[0] = DOS[i]*dfreq[i]*de_dT[i]*Vg[i];
	cumu_col[0] = DOS[i]*dfreq[i]*de_dT[i]/RTime[i];
	
	// iterations
	for (i = 1; i < NumBins; i++)
	{
		cumu_b[i] = cumu_b[i-1] + DOS[i]*dfreq[i]*de_dT[i];
		cumu_v[i] = cumu_v[i-1] + DOS[i]*dfreq[i]*de_dT[i]*Vg[i];
		cumu_col[i] = cumu_col[i-1] + DOS[i]*dfreq[i]*de_dT[i]/RTime[i];
	}
	C = cumu_b[NumBins-1];                            // heat capacity.

	// Write Out Cumu functions.
	fpcumu = fopen("./CumuFunctions.out","a+");
	fprintf(fpcumu,"  CumuBulk  |   CumuVel   | CumuCol  | de_dT\n");
	for (i = 0; i < NumBins; i++)
		fprintf(fpcumu,"%e  %e  %e  %e\n",cumu_b[i], cumu_v[i], cumu_col[i], de_dT[i]);
	fclose(fpcumu);
	// Energy calculation for the Left, Right and Bulk regions. Remains to be corrected.
	// enrgInit = Volume*C*fabs(Tref-Tinit);
	enrgLeft = AreaLeft*cumu_v[NumBins-1]*fabs(Tleft - Tref)/4;
	enrgRight = AreaRight*cumu_v[NumBins-1]*fabs(Tright - Tref)/4;
	// Total deviational energy.
	enrgTot = enrgLeft + enrgRight; //enrgInit + enrgLeft + enrgRight;
	// effective deviational energy.
	Eeff = enrgTot/NumParticles;
	// Normalization of the cumulative functions
	for (i = 0; i < NumBins; i++)
	{
		cumu_b[i] = cumu_b[i]/cumu_b[NumBins-1];
		cumu_v[i] = cumu_v[i]/cumu_v[NumBins-1];
		cumu_col[i] = cumu_col[i]/cumu_col[NumBins-1];
	}
	// Cumulative functions finished.

	// Output the info for check
	fpchecker = fopen("./CheckerPhonon.out","a+");
	fprintf(fpchecker,"---------------------Material_Property_for_Check-----------------------\n");
	fprintf(fpchecker,"NumBins:%d\n",NumBins);
	fprintf(fpchecker,"freq[rad/s] | DOS | Vg[m/s] | dfreq | Rtime[s] | Branch\n");
	for (i = 0; i < NumBins; i++)
		fprintf(fpchecker,"%e  %e  %e  %e  %e  %d\n", freq[i], DOS[i], Vg[i], dfreq[i], RTime[i], Polar[i]);
	fclose(fpchecker);
}
// build the cumulative function for Inflow Surfaces Area.
// Done before call transmission function, After buiding geometry
void buildCumuAreas()
{
	int i, surfnum;
	double sum;
	CumuAreasL = (double*)malloc(NumInSurfaces*sizeof(double));
	CumuAreasR = (double*)malloc(NumOutSurfaces*sizeof(double));

	sum = 0;
	for (i = 0; i < NumInSurfaces; i++)
	{
		surfnum = InSurfaces[i];
		sum = sum + SurfaceAreas[surfnum];
		CumuAreasL[i] = sum;
	}
	// Normalize
	for (i = 0; i < NumInSurfaces; i++)
		CumuAreasL[i] = CumuAreasL[i]/sum;
	// Right Surfaces
	sum = 0;
	for (i = 0; i < NumOutSurfaces; i++)
	{
		surfnum = OutSurfaces[i];
		sum = sum + SurfaceAreas[surfnum];
		CumuAreasR[i] = sum;
	}
	// Normalize
	for (i = 0; i < NumOutSurfaces; i++)
		CumuAreasR[i] = CumuAreasR[i]/sum;
	// done!
}

// Build the Spatial Cell Volume && range
void SpatialCellBuilder(char *fname)
{
	int i;
	FILE *fp, *fpchecker;
	fp = fopen(fname,"r");
	Volume = 0;                                    // Total Volume of the system
	// FORM: Z_min, Zmax, Volume.  Unit: SI, m.
	fscanf(fp,"%d",&Ztot);
	Vs = (double*)malloc(Ztot*sizeof(double));
	Zmins = (double*)malloc(Ztot*sizeof(double));
	Zmaxs = (double*)malloc(Ztot*sizeof(double));
	Zcenters = (double*)malloc(Ztot*sizeof(double));
	for (i = 0; i < Ztot; i++)
	{
		fscanf(fp,"%lf", &Zmins[i]);
		fscanf(fp,"%lf",&Zmaxs[i]);
		fscanf(fp,"%lf",&Vs[i]);
		Zcenters[i] = 0.5*(Zmins[i]+Zmaxs[i]);
		Volume += Vs[i];                          // sum up for the total volume of the system.
	}
	fclose(fp);
	// No need to discretize along X or Y directions.
	// Output the readin info.
	fpchecker = fopen("CheckSpatialCells.out","a+");
	fprintf(fpchecker,"-----------------Spatial_Cells_for_Check------------------\n");
	fprintf(fpchecker, "Num.of.Cells in Z direction: %d\n",Ztot);
	for (i = 0; i < Ztot; i++)
		fprintf(fpchecker,"%e(m)  %e(m)  %e(m^3)\n",Zmins[i], Zmaxs[i], Vs[i]);
	fprintf(fpchecker,"Total Volume:%e m^3\n", Volume);
	printf("Spatial Cells Build Succeed!\n");
	fclose(fpchecker);
}

void StatBuilder()
{
	int i, j;
	Qz = (double*)malloc(Ztot*sizeof(double*));
	T = (double*)malloc(Ztot*sizeof(double*));
	for (i = 0; i < Ztot; i++)
	{
		T[i] = 0.0;
		Qz[i] = 0.0;
	}
	printf("Statistic Build Succeed\n");
}
