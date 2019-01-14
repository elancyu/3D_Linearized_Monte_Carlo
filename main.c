#include<stdio.h>
#include<stdlib.h>
#include"GlobalVariables.h"
#include"Functions.h"

int main()//int argc, char* argv[])
{
	// depends on the parameters input to decide what to do.
	// input: ./RayTracing indicator , geofile,  phonon info.
	char *pGeometry;
	char pMaterial[100];
	char *pSpatial;
	int np = 3;  // testing value.
	pGeometry = "geo.in";// argv[1];
	pSpatial = "SpatialCells.in";
	GeoBuilder(pGeometry);
	SpatialCellBuilder(pSpatial);
	StatBuilder();
	sprintf(pMaterial,"Si%d.in",(int)Tref);
	PhononBuilder(pMaterial);
	PerformMC(np);
	Logger();
	printf("Calculation Done!\n");
	return 1;
}
