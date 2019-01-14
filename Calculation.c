/* This block of code read in bulk thermal conductivity MFP spectra and
launch the Ray Tracing to calculate the effective MFP for nanostructure
, and complish the integration procedure to get the thermal conductivity*/

// input file: bulktcspectra.in
// output file: Effective BS rate plot.
#include<stdio.h>
#include<stdlib.h>
#include"GlobalVariables.h"
#include"Functions.h"
#include<string.h>
#include<math.h>
#include<time.h>

void PerformMC()
{
	double start, finish;
	// needs a local copy of surface norm and displacement!!!
	int i, j, Nz, out;
	double L1, L2, Ltot, dist2scat, dist2next, deltat;               // for steady sampling.
	double x1, y1, z1, rez, tmp;
	double Ri, R1,R2, theta, phi;
	particle ptc;                    // the particle for ray tracing MC.
	long fracnum, ti;
	// for raytracing
	double a,b,c,d,e, cosangle;
	double ux,uy,uz,vx,vy,vz,wx,wy,wz,u,v;
	double position[3];
	double dist, dotprod, thetaRot, phiRot;
	int hitface, p1, p2, p3, nCold, nCnew;
	double s1, s2, s3, s4, c1, c2, c3, c4;
	// for interpolation.
	double ca, cb, cc;
	int surfnum, L, R, M;
	FILE* fp;
	// Local initialization of statistics
	start = clock();
	fracnum = NumParticles * 5 / 100;
	// Perform Ray Tracing like MC calculation.
	for (ti = 0; ti < NumParticles; ti++)
	{
		out = 0;                          // still inside the simulation domain.
		Ri = RandR();
		// only the master thread output the progress
		if (ti%fracnum==0)
			printf("%d%% completed\n",ti/fracnum*5);
		if (Ri < enrgLeft/enrgTot)   // Left Surfaces.
		{
			// Generate the particle in the left contact.
			// Firstly Choose the surface and then determine the initial position.
			Ri = RandR();
			L = -1; R = NumInSurfaces - 1;
			while(R - L > 1)
			{
				M = floor((L+R)/2);
				if (Ri < CumuAreasL[M])
					R = M;
				else
					L = M;
			}
			surfnum = InSurfaces[R];                    // the surface number of the in surface.
			// Now determine the location of the particle.
			R1 = RandR();
			R2 = RandR();
			ptc.Surface = surfnum;
			p1 = Surfaces[surfnum][0];                 // 1st point.
			p2 = Surfaces[surfnum][1];                 // 2nd point.
			p3 = Surfaces[surfnum][2];                 // 3rd point.
			// coefficients
			ca = 1 - sqrt(R1);
			cb = sqrt(R1)*(1-R2);
			cc = R2*sqrt(R1);
			// Initialize the Positions.
			ptc.x = ca*Points[p1][0] + cb*Points[p2][0] + cc*Points[p3][0];
			ptc.y = ca*Points[p1][1] + cb*Points[p2][1] + cc*Points[p3][1];
			ptc.z = ca*Points[p1][2] + cb*Points[p2][2] + cc*Points[p3][2];

			// Convert to global coordinate.
			ptc.Pnum = 0;
			ptc.realz = ptc.Pnum*UnitCellLength + ptc.z;

			// Generate the traveling direction of the particle
			// theta = acos(1-RandR());        // pure geometric, no inclusion of velocity.
			theta = asin(sqrt(RandR()));     // inclusion of velocity.
			phi = 2*PI*RandR();                       // PI need to be defined.
			ptc.tx = sin(theta)*cos(phi);
			ptc.ty = sin(theta)*sin(phi);
			ptc.tz = cos(theta);                      // update traveling direction.
			// Choose the mode.
			ptc.NumMode = ChooseMode(2);
			// Choose the initial time.
			ptc.sign = fabs(Tleft - Tref)/(Tleft - Tref);
		}
		else                                          // Right Surface.
		{
			// Firstly Choose the surface and then determine the initial position.
			Ri = RandR();
			L = -1; R = NumInSurfaces - 1;
			while(R - L > 1)
			{
				M = floor((L+R)/2);
				if (Ri < CumuAreasR[M])
					R = M;
				else
					L = M;
			}
			surfnum = OutSurfaces[R];                    // the surface number of the in surface.
			// Now determine the location of the particle.
			R1 = RandR();
			R2 = RandR();
			ptc.Surface = surfnum;
			p1 = Surfaces[surfnum][0];                 // 1st point.
			p2 = Surfaces[surfnum][1];                 // 2nd point.
			p3 = Surfaces[surfnum][2];                 // 3rd point.

			// coefficients
			ca = 1 - sqrt(R1);
			cb = sqrt(R1)*(1-R2);
			cc = R2*sqrt(R1);
			// Initialize the Positions.
			ptc.x = ca*Points[p1][0] + cb*Points[p2][0] + cc*Points[p3][0];
			ptc.y = ca*Points[p1][1] + cb*Points[p2][1] + cc*Points[p3][1];
			ptc.z = ca*Points[p1][2] + cb*Points[p2][2] + cc*Points[p3][2];
			// Convert to the global coordinate.
			ptc.Pnum = NumPeriods - 1;
			ptc.realz = ptc.Pnum*UnitCellLength + ptc.z;
			// Generate the travelling direction of the particle
			// theta = acos(1-RandR());              // pure geometric, no inclusion of velocity
			theta = asin(sqrt(RandR()));            // inclusion of velocity
			phi = 2*PI*RandR();                       // PI need to be defined.
			ptc.tx = sin(theta)*cos(phi);
			ptc.ty = sin(theta)*sin(phi);
			ptc.tz = -cos(theta);                      // update traveling direction.
			ptc.NumMode = ChooseMode(2);
			// Choose the initial time.
			ptc.sign = fabs(Tright - Tref)/(Tright - Tref);
		}

		// Mark the time location.
		// firstly we calculate the next location.
		// after that, we change the mode.
		while (!out)
		{
			dist2scat = -log(RandR())*RTime[ptc.NumMode]*Vg[ptc.NumMode];        // distance to next intrinsic scattering.

			// This block of code: find out the surface that it may scatter.
			for (i = 0; i < NumSurfaces; i++)
			{
				// calculate the angle cosine value
				a = SurfaceNorms[i][0];
				b = SurfaceNorms[i][1];
				c = SurfaceNorms[i][2];
				cosangle = a*ptc.tx+b*ptc.ty+c*ptc.tz;
				// check if the surface normal is the same direction as the particle.
				if (cosangle>0)
				{
					SurfaceNorms[i][0] = -a;
					SurfaceNorms[i][1] = -b;
					SurfaceNorms[i][2] = -c;
					Displacement[i] = -Displacement[i];
					cosangle = SurfaceNorms[i][0]*ptc.tx+SurfaceNorms[i][1]*ptc.ty+SurfaceNorms[i][2]*ptc.tz;
				}
				// exclude the previous surface and parallel surface
				if (i!=ptc.Surface && cosangle!=0)
				{
					// Calculate the distance to surface along the path.
					a = SurfaceNorms[i][0];
					b = SurfaceNorms[i][1];
					c = SurfaceNorms[i][2];
					// traveling distance along traveling direction.
					dist2surf[i] = -(a*ptc.x+b*ptc.y+c*ptc.z+Displacement[i])/cosangle;
					// dist2surf < 0 means it will not hit on the surface.
					if (dist2surf[i]<0)
						dist2surf[i] = INF;
					else
					{
						// Find the intercept point location.
						position[0] = ptc.x + dist2surf[i]*ptc.tx;
						position[1] = ptc.y + dist2surf[i]*ptc.ty;
						position[2] = ptc.z + dist2surf[i]*ptc.tz;
						// check if the intercepting point is outside the triangle.
						p1 = Surfaces[i][0];
						p2 = Surfaces[i][1];
						p3 = Surfaces[i][2];
						// U
						ux = Points[p2][0] - Points[p1][0];
						uy = Points[p2][1] - Points[p1][1];
						uz = Points[p2][2] - Points[p1][2];
						// V
						vx = Points[p3][0] - Points[p1][0];
						vy = Points[p3][1] - Points[p1][1];
						vz = Points[p3][2] - Points[p1][2];
						// W
						wx = position[0] - Points[p1][0];
						wy = position[1] - Points[p1][1];
						wz = position[2] - Points[p1][2];
						// a
						a = ux*ux + uy*uy + uz*uz;
						// b
						b = ux*vx + uy*vy + uz*vz;
						// c
						c = vx*vx + vy*vy + vz*vz;
						// d
						d = wx*ux + wy*uy + wz*uz;
						// e
						e = wx*vx + wy*vy + wz*vz;
						// u & v: coefficients
						u = (c*d-b*e)/(a*c-b*b);
						v = (a*e-b*d)/(a*c-b*b);
						if (u>=0 && v>=0 && u+v<=1)
							;
						else
							dist2surf[i] = INF;
					}
				}
				else
					dist2surf[i] = INF;
			}   // end of for loop: find the collided surface.
			// find the nearest surface
			dist2next = dist2surf[0];
			hitface = 0;
			for (i = 1; i < NumSurfaces; i++)
			{
				if (dist2surf[i]<dist2next)
				{
					dist2next = dist2surf[i];
					hitface = i;
				}
			}

			// Update the new location of the particle.
			if (dist2scat < dist2next)                      // intrinsic scattering: distance to next intrinsic scattering smaller than boundary scattering.
				dist = dist2scat;
			else
				dist = dist2next;

			x1 = ptc.x + dist*ptc.tx;
			y1 = ptc.y + dist*ptc.ty;
				// here the z value should be the real z value in global coordinate.
			z1 = ptc.z + dist*ptc.tz;
			rez = ptc.realz + z1 - ptc.z;

			// Stage 1: find out the num cell of the previous location and the new location
			nCold = LocateCellNum(ptc.realz);
			nCnew = LocateCellNum(rez);
			// Case 1: both are within the same cell
			if (nCold == nCnew)
			{
				Ltot = rez - ptc.realz;
				deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
				T[nCold] += deltat*Eeff*ptc.sign/(C*Vs[nCold]);
				Qz[nCold] += Ltot*Eeff*ptc.sign/Vs[nCold];
			}
			// Case 2: nCnew > nCold: ptc.tz > 0;
			else if (nCnew > nCold)
			{
				// in the cell nCold
				i = nCold;
				Ltot = Zmaxs[i] - ptc.realz;
				deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
				T[i] += deltat*Eeff*ptc.sign/(C*Vs[i]);
				Qz[i] += Ltot*Eeff*ptc.sign/Vs[i];
				// in the central cells
				for (i = nCold + 1; i <nCnew; i++)
				{
					Ltot = Zmaxs[i] - Zmins[i];
					deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
					T[i] += deltat*Eeff*ptc.sign/(C*Vs[i]);
					Qz[i] += Ltot*Eeff*ptc.sign/Vs[i];
				}
				// in the cell nCnew
				i = nCnew;
				Ltot = rez - Zmins[i];
				deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
				T[i] += deltat*Eeff*ptc.sign/(C*Vs[i]);
				Qz[i] += Ltot*Eeff*ptc.sign/Vs[i];
			}
			// case 3: nCnew < nCold: ptc.tz < 0;
			else
			{
				// in the cell nCnew
				i = nCnew;
				Ltot = rez - Zmaxs[i];
				deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
				T[i] += deltat*Eeff*ptc.sign/(C*Vs[i]);
				Qz[i] += Ltot*Eeff*ptc.sign/Vs[i];
				// in the central cells
				for (i = nCnew + 1; i < nCold; i++)
				{
					Ltot = Zmins[i] - Zmaxs[i];
					deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
					T[i] += deltat*Eeff*ptc.sign/(C*Vs[i]);
					Qz[i] += Ltot*Eeff*ptc.sign/Vs[i];
				}
				// in the cell nCold
				i = nCold;
				Ltot = Zmins[i] - ptc.realz;
				deltat = Ltot/(Vg[ptc.NumMode]*ptc.tz);
				T[i] += deltat*Eeff*ptc.sign/(C*Vs[i]);
				Qz[i] += Ltot*Eeff*ptc.sign/Vs[i];
			}
			// avoid the loop over all the spatial cells.

			// Update the position of the particle
			ptc.x = x1;
			ptc.y = y1;
			ptc.z = z1;
			ptc.realz = rez;
			// intrinsic scattering
			if (dist2scat < dist2next)                      // intrinsic scattering.
			{
				// randomly set the new traveling direction
				theta = acos(1-2*RandR());              // Math functions like asin and sqrt need to be defined.
				phi = 2*PI*RandR();                       // PI need to be defined.
				ptc.tx = sin(theta)*cos(phi);
				ptc.ty = sin(theta)*sin(phi);
				ptc.tz = cos(theta);                      // update traveling direction.
				ptc.Surface = -1;                         // set the last hit surface as -1;
				ptc.NumMode = ChooseMode(3);               // choose the mode after collision.
			}
			// collide with surfaces.
			else if (SurfaceMarker[hitface]==1)            // SurfaceMarker=1: In-flow Surfaces.
			{
				if (ptc.Pnum == 0)                       // already in the left most unit cell and exit from left boundary
					out = -1;
				else
				{
					ptc.Pnum = ptc.Pnum - 1;           // update the current unit cell location
					ptc.z = UnitCellLength;            // Move it from the N-th Unit left boundary to the (N-1)th Unit right boundary; no need to change realz.
					for (j = 0; j < NumOutSurfaces; j++)  // Find out the current surfaces.
					{
						if (IsInsideTriangle2D(ptc.x, ptc.y, OutSurfaces[j]))
						{
							ptc.Surface = OutSurfaces[j];
							break;
						}
					}
				}
			}
			// type II
			else if (SurfaceMarker[hitface]==2)            // SurfaceMarker=2: Out-flow Surfaces.
			{
				if (ptc.Pnum == NumPeriods - 1)         // alreay in the right most unit cell and exit from the right boundary. no need to change realz.
					out = 1;
				else
				{
					ptc.Pnum = ptc.Pnum + 1;           // Update the current unit cell
					ptc.z = 0;                         // Move it from the right boundary to the left boundary
					for (j = 0; j < NumInSurfaces; j++)
					{
						if (IsInsideTriangle2D(ptc.x, ptc.y, InSurfaces[j]))   // find the corresponding surfaces
						{
							ptc.Surface = InSurfaces[j];
							break;
						}
					}
				}
		 }
				// 	Regular surfaces.
			else
			{
				ptc.Surface = hitface;
				// specular
				if (RandR()<=Specularity[hitface])
				{
					// change direction.
					dotprod = -2*(ptc.tx*SurfaceNorms[hitface][0]+ptc.ty*SurfaceNorms[hitface][1]+ptc.tz*SurfaceNorms[hitface][2]);
					ptc.tx += dotprod*SurfaceNorms[hitface][0];
					ptc.ty += dotprod*SurfaceNorms[hitface][1];
					ptc.tz += dotprod*SurfaceNorms[hitface][2];
					// Nothing else.
				}
				// diffuse
				else
				{
					// Firstly generate the direction and the multiply rotation matrix.
					// original angles.
					theta = asin(sqrt(RandR()));                  //
					phi = 2*PI*RandR();                           // PI is constant.
					// X,Y,Z partial directions in original Coordinate.
					// Rotation Angles
					phiRot = atan2(SurfaceNorms[hitface][1],SurfaceNorms[hitface][0]);    // thetaRot = atan2(ny,nx); atan2 to be defined.
					thetaRot = acos(SurfaceNorms[hitface][2]);                                // phiRot = asin(nz);
					// the cosine and sine value for these four angles
					s1 = sin(theta);
					c1 = cos(theta);
					s2 = sin(phi);
					c2 = cos(phi);
					s3 = sin(thetaRot);
					c3 = cos(thetaRot);
					s4 = sin(phiRot);
					c4 = cos(phiRot);
					// Directly have the direction after reflection.
					ptc.tx = (s1*c2*c3*c4-s1*s2*s4+c1*s3*c4);
					ptc.ty = (s1*c2*c3*s4+s1*s2*c4+c1*s3*s4);
					ptc.tz = (-s1*c2*s3+c1*c3);
					// hit the surface, not need to update the position.
				}
			}
		}
	}// for loop
	finish = clock();
	elapsedtime = (finish - start)/ CLOCKS_PER_SEC;
}


// Calculate triangle areas
double TriArea(double p1[3], double p2[3], double p3[3])
{
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double a, b, c, area;
	x1 = p1[0]; y1 = p1[1]; z1 = p1[2];
	x2 = p2[0]; y2 = p2[1]; z2 = p2[2];
	x3 = p3[0]; y3 = p3[1]; z3 = p3[2];
	a = (y1-y2)*(z1-z3)-(z1-z2)*(y1-y3);
	b = (z1-z2)*(x1-x3)-(x1-x2)*(z1-z3);
	c = (x1-x2)*(y1-y3)-(y1-y2)*(x1-x3);
	area = sqrt(a*a+b*b+c*c);
	return (area/2);
}

// ChooseMode: May need a double check.
int ChooseMode(int type)
{
	int L,R, M;            // for binary search.
	double Rm;          // the random number for choosing the mode
	Rm = RandR();
	// type I:Bulk Region; type II:Left Boundary & Right Boundary
	if (type==1)
	{
		L = -1;
		R = NumBins - 1;
		while (R - L > 1)
		{
			M = floor((L+R)/2);
			if (Rm < cumu_b[M])
				R = M;
			else
				L = M;
		}
		return R;
	}
	if (type==2)
	{
		L = -1;
		R = NumBins - 1;
		while (R - L > 1)
		{
			M = floor((L+R)/2);
			if (Rm < cumu_v[M])
				R = M;
			else
				L = M;
		}
		return R;
	}
	if (type==3)
	{
		L = -1;
		R = NumBins - 1;
		while (R - L > 1)
		{
			M = floor((L + R) / 2);
			if (Rm < cumu_col[M])
				R = M;
			else
				L = M;
		}
		return R;
	}
}


// This function is used to output the results. Namely the heat flux, Temperature and effective thermal conductivity.
void Logger()
{
	int i, j;
	double kappa;                           // effective thermal conductivity.
	double TLeft, TRight, deltaT;           // Extrapolated boundary temperature.
	double Qmean = 0.0;                           // average heat flux.
	FILE *fp;

	// Calculate the effective thermal conductivity and plot the average heat flux.

	fp = fopen("./Steady.dat","a+");
	fprintf(fp,"Variables = \"Position\",\"Temperature,K\",\"Heat flux\"\n");
	fprintf(fp,"Zone T = \"Last Moment\", I = %d, DataPacking = Point\n",Ztot);
	for (i = 0; i < Ztot; i++)
	{
		Qmean = Qmean + Qz[i];
		fprintf(fp,"%e  %e  %e\n",Zcenters[i], T[i]+Tref, Qz[i]);
	}
	fclose(fp);

	TLeft = 1.5*T[0] - 0.5*T[1];
	TRight = 1.5*T[Ztot-1] - 0.5*T[Ztot-2];
	deltaT = TLeft - TRight;
	Qmean = Qmean/Ztot;
	kappa = Qmean*TotalLength/deltaT;
	fp = fopen("./EffectiveTC.txt","a+");
	fprintf(fp,"Effective TC:%.6f\n",kappa);
	fprintf(fp,"elasped time:%.3lf s\n",elapsedtime);
	fclose(fp);
}

void CalRRTime()
{
	// The formulae and coefficients from A.J.Minnich PhD thesis.
	// Reciprocal Relaxation Time is combined with Matthiessen's Rule.
	int i, polar;
	// LA branch;
	double AL, AT,AI, Theta, w;
	AL = 2e-19;
	AT = 1.2e-19;
	AI = 3e-45;
	Theta = 80;
	// TA branch;
	for (i = 0; i < NumBins; i++)
	{
		polar = Polar[i];
		w = freq[i];                           // angular frequency.
		RRTime[i] = AI*w*w*w*w;                // Impurity.
		if (polar==1)                          // LA branch
			RRTime[i] = RRTime[i] + AL*w*w*pow(Tref,1.49)*exp(-Theta/Tref);
		if (polar==2)                          // TA branch
			RRTime[i] = RRTime[i] + AT*w*w*pow(Tref,1.65)*exp(-Theta/Tref);
	}
}

// IsInsideTriangle2D
// Input: position of the point, the surface number.
// Output: true:1; false:0
int IsInsideTriangle2D(double px, double py, int SurfaceNum)
{
	int p1, p2, p3;               // three vertices of the triangle.
	int n = SurfaceNum;           // surface number.
	double a1,a2,b1,b2,c1,c2;     // for telling the inside or not
	double t1, t2, t3;
	p1 = Surfaces[n][0];
	p2 = Surfaces[n][1];
	p3 = Surfaces[n][2];
	a1 = px - Points[p1][0];
	a2 = py - Points[p1][1];
	b1 = px - Points[p2][0];
	b2 = py - Points[p2][1];
	c1 = px - Points[p3][0];
	c2 = py - Points[p3][1];
	t1 = (a1*b2-a2*b1);
	t2 = (b1*c2-b2*c1);
	t3 = (c1*a2-c2*a1);
	if (t1*t2>0 && t1*t3>0)
		return 1;
	else
		return 0;
}

// use binary search to locate the unit cell number of the location.
int LocateCellNum(double z)
{
	int CellNum;
	int L, R, M;
	// use Zmaxs to locate
	L = -1;
	R = Ztot-1;
	while (R - L > 1)
	{
		M = (L+R)/2;
		if (z < Zmaxs[M])
			R = M;
		else
			L = M;
	}
	return R;
}
