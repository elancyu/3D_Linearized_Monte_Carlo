/* Free Memory */
// End of the program: cleanup the memory space
// Input:
// Output None

#include<stdio.h>
#include<stdlib.h>
#include"Functions.h"
#include"GlobalVariables.h"

void FreeMem(int Indicator)
{
	//inidicator = 0: doesn't calculate TC; indicator = 1: Calculate TC directly.
	int i;
	for (i = 0; i < NumPoints; i++)
		free(Points[i]);
	free(Points);
	for (i = 0; i < NumSurfaces; i++)
	{
		free(SurfaceNorms[i]);
		free(Surfaces[i]);
	}
	free(SurfaceNorms);
	free(Surfaces);
	free(SurfaceAreas);
	free(Specularity);
	free(InSurfaces);
	free(OutSurfaces);
	free(Displacement);
	free(SurfaceMarker);
	free(CumuAreasL);
	free(CumuAreasR);
	free(cumu_b);
	free(cumu_v);
	free(cumu_col);
	free(freq);
	free(DOS);
	free(RTime);
	free(dfreq);
	free(de_dT);
	free(Polar);
	free(Zmins);
	free(Zmaxs);
	free(Zcenters);
	free(T);
	free(Qz);
}
