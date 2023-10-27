/* ********************************************************************
 * Copyright (c) 2019 Biological Modeling and Visualization Laboratory,
 * University of Calgary. All rights reserved.
 * ********************************************************************/

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <lpfgall.h>
#include "bspline.h"
#include "markers.h"

#define DRAW_STEP 1

#define SHOW_ACTIVE_RING            1
#define SHOW_PRIMORDIA_CENTERS      0
#define SHOW_RHO_SCALE              0.3
#define SHOW_PRIMORDIA_BOUNDARIES   1
#define BP                          64 // number of boundary points
#define RING_WIDTH                  0.3
#define BOUNDARY_WIDTH              0.3


/* constants */
const double RadToDeg = 180.0 / M_PI;
const double TwoPi = 2. * M_PI;

const double dphi = 0.001; // angular step for scanning intervals, in radians
#define ELLIPTIC_RING 0 // enable elliptic active ring
#define ELLIPTIC_RING_BP int(TwoPi/dphi)
#define ELLIPTIC_RING_FUNC(dr,ar,phi) 0.002*std::min(dr-ar,ar)*std::sin(2.*phi-(160.0-20.0)/RadToDeg)

/* constants for repelling forces */
const double k1 = 0.02266; // - k1 * (distance - thresholdDistance)

const double DAGE = 1.0;

/* Global variables */
const double contourStep0 = 0.5*0.001;
double contourStep = contourStep0; // interpolation step through contours per dProgress step


const double dProgress = 1.;
double Progress;	 // progress of development
				     // needed to interpolate between initial and final profile
double normProgress; // normalized progress (from 0 to 1)

double L;	 	 	// length of the generating curve

int n;		// number of primordia formed so far;

double thresholdDistance;
const double initialThresholdDistance = 3.5*M_PI ; // Multiplier in threshold for insertion
#define USER_FUNC_TH_DIST(t)		(func(thresholdDistanceFunc,t))

/* Interpolate between the initial and final generating curve ("profile") */

const int NUM_CNTRS = 7; // number of contours
const double CONTOUR_SCALE = 100.0;
const int CLOSED_BSPLINE = 1;

const int CONTOUR_WEIGHTS_SPLINE_DEGREE = 2;
BSplineBasisFunc contourWeights(NUM_CNTRS, CONTOUR_WEIGHTS_SPLINE_DEGREE, CLOSED_BSPLINE);

const int MARKER_WEIGHTS_SPLINE_DEGREE = 2;
BSplineBasisFunc markerWeights(NUM_CNTRS, MARKER_WEIGHTS_SPLINE_DEGREE, CLOSED_BSPLINE);

void resetWeights (double alpha) {
	contourWeights.computeFuncs(alpha);
	markerWeights.computeFuncs(alpha);
}

BSplineBasisFunc primordiaMarkerWeights(0,MARKER_WEIGHTS_SPLINE_DEGREE,1);

vector<double> activeRingMarkers;
vector<vector<double>> primordiaMarkers;
const int NUMBER_OF_BRACT_MARKERS = 2;  

void ComputeMarkers (void) {

    activeRingMarkers.clear();

    int maxLandmarks = myCurveMarkers(NUM_CNTRS)-1; // -1 for active ring marker

    primordiaMarkers = vector<vector<double>> (maxLandmarks+1,vector<double>(NUM_CNTRS,0.0)); // + 1 for s = 1

    // make an array of active ring markers and another of the primordia markers
   	for (int cntr = 1; cntr <= NUM_CNTRS; cntr++) {
        int numMarkers = myCurveMarkers(cntr);
        activeRingMarkers.push_back(myCurveMarker(cntr,numMarkers));

        primordiaMarkers[0][cntr-1] = 1.0;
        for (int i = 1; i < numMarkers; i++)
            primordiaMarkers[i][cntr-1] = myCurveMarker(cntr,i);
    } 

    // propagate landmarks to earlier contours
    for (int i = 0; i < primordiaMarkers.size(); i++) {
        // find the first non-zero entry
        int zeroIndex = primordiaMarkers[i].size() - 1;
        while(primordiaMarkers[i][zeroIndex-1] > 0. && zeroIndex-1 >= 0)
               --zeroIndex;

        // if there are zero entries
        if (zeroIndex-1 > 0) {
            // find the ratio of expansion between the new landmark and existing one
            double ratio = primordiaMarkers[i][zeroIndex] / primordiaMarkers[i-1][zeroIndex];
            // propagate the ratio to the earlier contours
            for (int k = 0; k < zeroIndex; k++) 
                    primordiaMarkers[i][k] = ratio * primordiaMarkers[i-1][k];
        }
    }

    // add landmarks at s = 0 for each contour
    primordiaMarkers.push_back(vector<double>(NUM_CNTRS,0.0));

    // create the B-spline weights for the primordia markers
    primordiaMarkerWeights = BSplineBasisFunc(primordiaMarkers.size(), MARKER_WEIGHTS_SPLINE_DEGREE, 1);
}



V3d profileXYZ(double ns) // C(t,s) in description
{
    extern BSplineBasisFunc contourWeights;
	V3d res(0,0,0);
	for (int i = 1; i <= NUM_CNTRS; i++) {
		res += contourWeights.valueAtFunc(i-1) * ((V3d) curveXYZ(i,ns)); // curveXYZ(k,s) is C_k(s) in description
	}
	
	return CONTOUR_SCALE * res;
}

V3d profileNormal(double ns) // not in description, for visualization only
{
    extern BSplineBasisFunc contourWeights;
	V3d res(0,0,0);
	for (int i = 1; i <= NUM_CNTRS; i++) {
		res += contourWeights.valueAtFunc(i-1) * ((V3d) curveNormal(i,ns)); 
	}
	
	return -res;
}

double hatNS(double u) // s(t,u) in description: weighted sum over all s_k(u)
{
    extern BSplineBasisFunc contourWeights;
    extern BSplineBasisFunc primordiaMarkerWeights;
    extern vector<vector<double>> primordiaMarkers;
    
    primordiaMarkerWeights.computeFuncs(u);

    double s = 0.;
    for (int i = 1; i <= NUM_CNTRS; i++) {
        double newNS = 0.;
        for (int j = 0; j < primordiaMarkers.size(); j++) {
            newNS += primordiaMarkerWeights.valueAtFunc(j) *
                    primordiaMarkers[primordiaMarkers.size()-1-j][i-1];
        }
        s += contourWeights.valueAtFunc(i-1) * newNS;
	}
    return s;
}


V3d hatProfileXYZ(double u) // \hat{C}(t,u) in description
{
    double ns = hatNS(u);
    return profileXYZ(ns);
}

double findU(double sgen) { // find u where s^{gen}(t) = s(t,u)

    // start from base of receptacle u=1, and find closest 'u' where s = hatNS(u)
    double u = 1.;
    while (sgen < hatNS(u) && u > 0.)
        u -= 1e-3;
    return u;
}

// interpolation function for active ring
double active_ring_newNS()  // s^{gen}(t) in description
{
    extern BSplineBasisFunc markerWeights;
    extern vector<double> activeRingMarkers;
    double new_ns = 0.;
    for (int i = 0; i < NUM_CNTRS; i++) {
        new_ns += markerWeights.valueAtFunc(i) * activeRingMarkers[i];
    }
	return new_ns;
}

double marker_newNS(int markerIndex) // trajectory of landmark. not in description
{
    extern BSplineBasisFunc markerWeights;
    extern vector<vector<double>> primordiaMarkers;

    double new_ns = 0.;
    for (int i = 0; i < NUM_CNTRS; i++) {
        new_ns += markerWeights.valueAtFunc(i) * primordiaMarkers[markerIndex][i];
    }
	return new_ns;
}



// Compute the length of the current (interpolated) generating curve
double CurveLength()
{
	const double nds = 0.0005; // normalized arc length used for integrating the generating curve length
	double AccumulatedLength = 0.;
	double ns = 0.;
	for (ns=0; ns<1.0-nds; ns+=nds)
  		AccumulatedLength += (profileXYZ(ns+nds) - profileXYZ(ns)).Length();
	AccumulatedLength += (profileXYZ(1.) - profileXYZ(ns)).Length();
	return AccumulatedLength;
}


struct primordium{
	double rho;	// primordium radius
	double ns;	// normalized curvilinear coordinate
    double inserted_u; // material point coordinate at time of insertion
	double phi;	// angular coordinate
	double h;	// vertical cylindrical coordinate
	double r;	// radial cylindrical cordinate
	double inserted_r;
	V3d pos;	// Cartesian position
	double age;
	int neighbour_left;
	int neighbour_right;
	double inserted_dist;
    bool stable;
};

primordium P[2000];		// array of primordia centers formed so far
primordium B[2000][BP];	// array of primordia represented by their boundaries.

primordium tp;	// trial primordium
#if ELLIPTIC_RING==1
primordium tpB[ELLIPTIC_RING_BP]; // active ring boundary
#endif

// Find new cylindrical and Cartesian coordinates, given 
// curvilinear coordinate p.s and angular coordinate p.phi.
void SyncCoordinates(primordium *p)
{
    V3d pt = profileXYZ(p->ns);

	p->r = pt.x;
	p->h = pt.y;

	p->pos.x = p->r*cos(p->phi); // convert to cartesian coordinates
	p->pos.z = p->r*sin(p->phi);
	p->pos.y = p->h;
	
	return;
}

double Distance (const primordium& P1, const primordium& P2, double C_arclen, double C_slen)
{
    // compute arc length between primordia using the smallest angle between them
	double angle = std::fabs(P1.phi - P2.phi);
    angle = std::min(angle,TwoPi-angle);
    double arcLen = (P1.r+P2.r)*0.5 * angle;

    extern double L;
    double dS = L * (P1.ns - P2.ns);
    double dist = std::sqrt(C_arclen*arcLen*arcLen + C_slen*dS*dS);

	return dist;
}


// Check if test primordium tp collides with previous primordia 

int HofmeisterTest(primordium *tp)
{
	extern int n;
	double min_dist1 = 1e11;
	double min_dist2 = 1e11+1.0;
	double min_dist3 = 1e11+2.0;
	int min_i1 = n+1;
	int min_i2 = n+2;
	int min_i3 = n+3;

	double dist = 0;
	for (int i = 0; i<n; i++) {	// iterate over primordia

		dist = Distance(*tp,P[i],1.,1.);
		if (dist < thresholdDistance)
			return 1;

		// check if this is the closest primordium
		if (dist < min_dist1) {
			min_dist3 = min_dist2;
			min_i3 = min_i2;
			min_dist2 = min_dist1;
			min_i2 = min_i1;
			min_dist1 = dist;
			min_i1 = i;
		} // else check if this is the second closest
		else if (dist < min_dist2) {
			min_dist3 = min_dist2;
			min_i3 = min_i2;
			min_dist2 = dist;
			min_i2 = i;
		}
		else if (dist < min_dist3) {
			min_dist3 = dist;
			min_i3 = i;
		}
		
	}

    // if there is one primordium or none, the min indicies should all be zero
    // this simplifies assignment of neighbours below
    if (n < 2) {
        min_i1 = min_i2 = min_i3 = 0;
    } 
    else {
        // if the three distances are almost the same,
        // the incipient primordium is tangential to three existing primordia
		if (std::fabs(min_dist2 - min_dist3) < 1e-5 &&
			std::fabs(min_dist1 - min_dist2) < 1e-5 && 
			std::fabs(min_dist1 - min_dist3) < 1e-5) {

			Printf ("Three tangential primordia: %d %d %d\n", min_i1+1,min_i2+1,min_i3+1);

            // pick the two youngest
			if (P[min_i1].age > P[min_i2].age && P[min_i1].age > P[min_i3].age) { // min_i1 is oldest
				min_i1 = min_i3;
			}
			else if (P[min_i2].age > P[min_i1].age && P[min_i2].age > P[min_i3].age) { // min_i2 is oldest
				min_i2 = min_i3;
			}
			//else { // min_i3 is oldest
			//}
		}
	}

    // save index of left and right neighbours (anti-clockwise and clockwise, respectively)
    // based on the angle between them
    double angle1 = P[min_i1].phi - tp->phi;
    angle1 -= TwoPi * std::floor(angle1/TwoPi); // wrap the angle from 0 to TwoPi

    // if the angle >= pi, min1 is CCW neighbour
    if (angle1 >= M_PI) {

        tp->neighbour_left = min_i1;
        P[min_i1].neighbour_right = n;

        tp->neighbour_right = min_i2;
        P[min_i2].neighbour_left = n;

    } else { // else min2 is CCW neighbour

        tp->neighbour_left = min_i2;
        P[min_i2].neighbour_right = n;

        tp->neighbour_right = min_i1;
        P[min_i1].neighbour_left = n;
    }

	return 0;	// if all primordia were far enough, collision has not occurred
}


int bractPrimordia; // number of bract primordia
std::vector<V2d> bpTable;
FILE *primsInitialFile = NULL;
FILE *primsStableFile = NULL;
FILE *primsFinalFile = NULL;
char primsInitialFilename[32] = "bellis_initial.prims";
char primsStableFilename[32] = "bellis_stable.prims";
char primsFinalFilename[32] = "bellis_final.prims";

module Draw();  // draws everything
module DrawReceptacle();
module DrawRing(double);

derivation length: 0;

Axiom: Draw() ;

Start:
{
	// read knot values if any
	FILE *posKnotFile = NULL;
	if ((posKnotFile = fopen("knotPositions.dat", "r")) != NULL) {
		int knot;
		double x,y,z;
		std::vector<double> knots;
		while (!feof(posKnotFile)) {
			if (fscanf(posKnotFile, "%d: %lg %lg %lg", &knot, &x, &y, &z) == 4) {
				knots.push_back(x);
			}
		}
		fclose(posKnotFile);
		if (!contourWeights.SetKnotSeq(knots)) {
            Printf ("Did not set contourWeights knots\n");
			contourWeights.ResetKnotSeq(CLOSED_BSPLINE);
		}
		if (!markerWeights.SetKnotSeq(knots)) {
            Printf ("Did not set markerWeights knots\n");
			markerWeights.ResetKnotSeq(CLOSED_BSPLINE);
		}
	}

    ComputeMarkers();

	Progress = 0.0;
	normProgress = 0.0;
	resetWeights(normProgress);

	thresholdDistance = initialThresholdDistance * USER_FUNC_TH_DIST(0.);

	L = CurveLength() ;
	
	tp.ns = active_ring_newNS();
	tp.rho = 0.;	// trial primordium
	tp.age = 0.;
	tp.neighbour_left = -1; 
	tp.neighbour_right = -1; 
	tp.inserted_dist = thresholdDistance;
    tp.stable = false;

	SyncCoordinates (&tp);

#if ELLIPTIC_RING==1
    double tpBphi = M_PI;
    for (int j=0; j<ELLIPTIC_RING_BP; j++) {
        tpB[j].ns = tp.ns;
        tpB[j].phi = tpBphi;
        tpBphi += TwoPi/double(ELLIPTIC_RING_BP);
		if (tpBphi > TwoPi)
			tpBphi -= TwoPi;
        SyncCoordinates(&tpB[j]);
    }
#endif

	n = 0; // there is no initial primordium
	bractPrimordia = 0;

	// create table for boundary points
	for (int j=0; j<BP; j++) {
		V2d pt;
		pt.x = sin(double(j)*TwoPi/double(BP));
		pt.y = cos(double(j)*TwoPi/double(BP));
		bpTable.push_back(pt);
	}

    primsInitialFile = fopen(primsInitialFilename, "w");
    fprintf(primsInitialFile, "%.5d\n", 0);
    primsStableFile = fopen(primsStableFilename, "w");
    fprintf(primsStableFile, "%.5d\n", 0);
}

StartEach: {

	// Simulate growth of the meristem
	// once most primordia are placed (active ring position is near 0)
	// reduce contourStep to smooth growth of flower head
	if (tp.ns < 0.01) {
		if (contourStep > contourStep0*0.1)
			contourStep *= 0.993;
	}

	normProgress += dProgress*contourStep;
	if (normProgress > 1.)
		normProgress = 1.;
	resetWeights(normProgress);

	// Calculate new generative curve length
	L = CurveLength() ;

	// Advance position of primordia
	for (int i=0; i<n; i++) {
		// advance primordia age
		P[i].age += DAGE * dProgress;

		// advance primordia centers
        P[i].ns = hatNS(P[i].inserted_u);

		// primordia shrink until last bract is produced:
		// save the thresholdDistance until the last bract is produced
        if (n <= 34) // changed to n <= 34 from (tp.ns > marker_newNS(NUMBER_OF_BRACT_MARKERS+1))
			P[i].inserted_dist = thresholdDistance;

		// recompute radius
		P[i].rho = P[i].r * std::sin(std::min(P[i].inserted_dist*0.5/P[i].r,M_PI));

#if SHOW_PRIMORDIA_BOUNDARIES==1
		// advance primordia boundaries
		if (n <= 34) { // changed to n <= 34 from (tp.ns > marker_newNS(NUMBER_OF_BRACT_MARKERS+1)) {
			// recompute the boundary points, delaying growth of boundary
			P[i].inserted_r = P[i].r;
			for (int j=0; j<BP; j++) {
				B[i][j].ns = P[i].ns + 0.95*P[i].rho * bpTable[j].x / L;
				B[i][j].phi = P[i].phi + 0.95*P[i].rho * bpTable[j].y / P[i].r;
				SyncCoordinates(&B[i][j]);
            }
		} else {
			// delay growth of boundary upto some age (as above)
			if (P[i].age < 10.) { // changed from 5 to 10
				P[i].inserted_r = P[i].r;
                for (int j=0; j<BP; j++) {
                    B[i][j].ns = P[i].ns + 0.95*P[i].rho * bpTable[j].x / L;
                    B[i][j].phi = P[i].phi + 0.95*P[i].rho * bpTable[j].y / P[i].r;
                    SyncCoordinates(&B[i][j]);
                }
		  	} else {
				for (int j=0; j<BP; j++) {
                    //check if 'u' has to be computed
					if (B[i][j].inserted_u == -1.) {
                        B[i][j].inserted_u = findU(P[i].ns + 0.95*P[i].rho * bpTable[j].x / L);
					} 
					B[i][j].phi = P[i].phi + 0.95*P[i].rho* bpTable[j].y / P[i].inserted_r;
                    B[i][j].ns = hatNS(B[i][j].inserted_u);
                    SyncCoordinates(&B[i][j]);
				}
		  	}
		}
#endif
	}

	// Apply lateral displacement (right is clockwise, left is counter clockwise)
	for (int i = 0; i < n; i++) {

		int left_i = P[i].neighbour_left;
		int right_i = P[i].neighbour_right;

        // only incipient primordium can move laterally
        // identified by being younger than both of its neighbours
		// extra check: if primordium has moved far from active ring,
		// it should not move to its older neighbour anymore.
        if (n > 1 && P[i].ns < 1.5*tp.ns && P[i].age <= P[left_i].age && P[i].age <= P[right_i].age) {

            double dPhi = 0.;

            if (P[left_i].age > P[right_i].age) { // CCW is older

                double arcLenL = Distance(P[i],P[left_i],1.,0.);

                double radiusL = 0.5*(P[i].r+P[left_i].r);

                double d = Distance(P[i],P[left_i],1.,1.);
                double a = arcLenL / d;

                // -k1*(-1) so lateral displacement is in CCW direction
                dPhi = a * k1 * (arcLenL - a * thresholdDistance) / radiusL;

            } else { // CW is older

                double arcLenR = Distance(P[i],P[right_i],1.,0.);

                double radiusR = 0.5*(P[i].r+P[right_i].r);

                double d = Distance(P[i],P[right_i],1.,1.);
                double a = arcLenR / d;

                dPhi = - a * k1 * (arcLenR - a * thresholdDistance) / radiusR;
            }

            // error checking: the lateral displacement should not be too large
            if (std::abs(dPhi) > 0.1) {
                Printf ("Warning: change in phi is > 0.1. P[%d]: %g\n", i+1, std::abs(dPhi));
                PauseSimulation();
            }

            // update the angle
            double newPhi = P[i].phi - dPhi * dProgress;
            newPhi -= TwoPi * std::floor(newPhi/TwoPi); // wrap newPhi from 0 to TwoPi
        
            P[i].phi = newPhi;
        } else {
            if (!P[i].stable) {
			    //Printf("Stable primordium %d %g\n", i+1, P[i].phi*RadToDeg);
        		SyncCoordinates(&P[i]); // Sync here because of outputting the radius
                fprintf(primsStableFile, "%d %d %g %g %g %d\n", i, StepNo(), P[i].ns, P[i].phi, (P[i].pos - B[i][0].pos).Length(), P[i].ns > marker_newNS(NUMBER_OF_BRACT_MARKERS+1));
                P[i].stable = true;
            }
        }

		SyncCoordinates(&P[i]);
	}	

	// Insert new primordium if they do not collide with previous ones
    if (tp.ns > 0.01) { // stop inserting primordia when the active ring has filled the receptacle

        double age = 0.;
#if ELLIPTIC_RING==0
        for (double phi = M_PI; phi < M_PI+TwoPi; phi += dphi) { 
            tp.phi = phi;
#else
        for (int j=0; j<ELLIPTIC_RING_BP; j++) {
            tp.phi = tpB[j].phi;
            tp.ns = tpB[j].ns;
#endif
            if (tp.phi > TwoPi)
                tp.phi -= TwoPi;

            if (HofmeisterTest(&tp) == 0) {
                Printf("Producing primordium %d %g\n", n+1, tp.phi*RadToDeg);

                // first five landmarks are for the bracts
                // so anything before the first ray floret marker is a bract 
                if (tp.ns > marker_newNS(NUMBER_OF_BRACT_MARKERS+1)) {
                    ++bractPrimordia;
                }
                Printf ("Bract primordia: %d\n", bractPrimordia);

                tp.age = age;
                tp.inserted_dist = thresholdDistance;
                tp.stable = false;
                tp.rho = tp.r * std::sin(std::min(thresholdDistance*0.5/tp.r,M_PI));
                tp.inserted_u = findU(tp.ns);
                SyncCoordinates (&tp);

                P[n] = tp;

#if SHOW_PRIMORDIA_BOUNDARIES==1
                for (int j=0; j<BP; j++) {
                    B[n][j].inserted_u = -1;
                    B[n][j].ns = tp.ns + 0.95*tp.rho * bpTable[j].x / L;
                    B[n][j].phi = tp.phi + 0.95*tp.rho * bpTable[j].y / tp.r;
                    SyncCoordinates(&B[n][j]);
                }
#endif
                fprintf(primsInitialFile, "%d %d %g %g %g %d\n", n, StepNo(), P[n].ns, P[n].phi, (P[n].pos - B[n][0].pos).Length(), P[n].ns > marker_newNS(NUMBER_OF_BRACT_MARKERS+1));

                ++n;
            }
        }
    }

	// Advance position of the active ring
	tp.ns = active_ring_newNS();
	SyncCoordinates(&tp);
#if ELLIPTIC_RING==1
    for (int j=0; j<ELLIPTIC_RING_BP; j++) {
        tpB[j].ns = tp.ns + ELLIPTIC_RING_FUNC(P[0].r,tp.r,tpB[j].phi);
        SyncCoordinates(&tpB[j]);
    }
#endif

	// stop simulation when the active ring has filled entire receptacle
	if (tp.ns < 0.005 && normProgress >= 1.) {
        // write number of primordia to begin of file
        fseek(primsInitialFile, 0, SEEK_SET);
        fprintf(primsInitialFile, "%.5d\n", n);
        fclose(primsInitialFile);

        fseek(primsStableFile, 0, SEEK_SET);
        fprintf(primsStableFile, "%.5d\n", n);
        fclose(primsStableFile);

        primsFinalFile = fopen(primsFinalFilename, "w");
        fprintf(primsFinalFile, "%.5d\n", n);
	    for (int i = 0; i < n; i++) {
            fprintf(primsFinalFile, "%d %d %g %g %g %d\n", i, StepNo(), P[i].ns, P[i].phi, (P[i].pos - B[i][0].pos).Length(), P[i].ns > marker_newNS(NUMBER_OF_BRACT_MARKERS+1));
        }
        fclose(primsFinalFile);
 
		Printf ("Simulation done. Pausing.\n");
		PauseSimulation();
	}

	// reduce the threshold distance for insertion of new primordia
	thresholdDistance = initialThresholdDistance * USER_FUNC_TH_DIST(normProgress);

	// advance time
	Progress += dProgress;

	if (StepNo() % DRAW_STEP*int(1./dProgress) == 0) {
		DisplayFrame();
	}

}

interpretation:
maximum depth: 2;

Draw():
{
    nproduce SetColor(15) DrawReceptacle;

#if SHOW_ACTIVE_RING==1
	// Draw active ring
#if ELLIPTIC_RING==0
	nproduce SB SetColor(254) MoveTo(0.,tp.pos.y,0.) SetCoordinateSystem(1.) DrawRing(tp.r) EB;
#else

	nproduce SB	CurrentContour(0) SetWidth (RING_WIDTH)  SetColor(254) ; 
		// draw the boundary as a generalized cylinder
		V3f h = (tpB[1].pos - tpB[0].pos).Normalize();
		V3f start_h = h; // save the initial heading
		nproduce SetHead3f(h) MoveTo3d(tpB[0].pos) StartGC;
		for (int j=0; j<ELLIPTIC_RING_BP-1; j++) {
			h = (tpB[j+1].pos - tpB[j].pos);
			float l = h.Length();
			h /= l;
			nproduce SetHead3f(h) F(l);	// iterate over remaining points
		}
		h = (tpB[0].pos - tpB[ELLIPTIC_RING_BP-1].pos);
		float l = h.Length();
		h /= l;
		nproduce SetHead3f(h) F(l);
		nproduce SetHead3f(start_h) F(0); // close the gen. cylinder given starting heading
		nproduce EndGC;
	nproduce EB;
#endif
#endif

	int colour;
#if SHOW_PRIMORDIA_BOUNDARIES==1
	nproduce SB;
	for (int i=0; i<n; i++) { 	// select color according to primordium number
	nproduce CurrentContour(0) SetWidth (BOUNDARY_WIDTH) ; 
		colour = 1;
		//if (n < 21)
		//	colour = 5; // was set to 6 to highlight the first 21 primorida
		//else if (i < bractPrimordia)
        if (i < 34)  
			colour = 5;
		nproduce SetColor(colour) ; 

		// draw the boundary as a generalized cylinder
		V3f h = (B[i][1].pos - B[i][0].pos).Normalize();
		V3f start_h = h; // save the initial heading
		nproduce SetHead3f(h) MoveTo3d(B[i][0].pos) StartGC;	// move to lat bounarhy point
		for (int j=0; j<BP-1; j++) {
			h = (B[i][j+1].pos - B[i][j].pos);
			float l = h.Length();
			h /= l;
			nproduce SetHead3f(h) F(l);	// iterate over remaining points
		}
		h = (B[i][0].pos - B[i][BP-1].pos);
		float l = h.Length();
		h /= l;
		nproduce SetHead3f(h) F(l);
		nproduce SetHead3f(start_h) F(0); // close the gen. cylinder given starting heading
		nproduce EndGC;
	}
	nproduce EB;
#endif

	// Draw primordia
#if SHOW_PRIMORDIA_CENTERS==1
    nproduce SetColor(1);
    for (int i=0; i<n; i++) { 	// select color according to primordium number
	    nproduce MoveTo3d(P[i].pos) Sphere(P[i].rho * SHOW_RHO_SCALE);
	}
#endif

	produce;	
}

DrawRing(s): {
	const double NUM_SEGS = 36.;
	const double theta = 360. / NUM_SEGS;
	const double step = 2.*M_PI*s / NUM_SEGS;
	nproduce SetWidth(RING_WIDTH) MoveTo(s,0.,0.) Up(90.) Right(0.5*theta) StartGC Left(theta);
	for (int i = 0; i < int(NUM_SEGS-1) ; i++) {
		nproduce F(step) Left(theta) ;
	}
	produce F(step) EndGC;
}

DrawReceptacle() : {

  // here we draw the surface of revolution using polygons, SP ... EP  

  const int TRANSVERSE_SECS = 80;
  const int EXCLUDING_APICAL_SEC = int(double(TRANSVERSE_SECS)*0.02);
  const int LONGITUDINAL_SECS = 40;

  nproduce SB;

  for (int j = 0; j < LONGITUDINAL_SECS; ++j) {
    double phi1 = double(j) / double(LONGITUDINAL_SECS) * 2. * M_PI;
    double phi2 = double(j+1) / double(LONGITUDINAL_SECS) * 2. * M_PI;

    V3d p1, n1, t1;
    V3d p2, n2, t2;
    V3d p3, n3, t3;
    V3d p4, n4, t4;

    // for first transverse section, compute the bottom two points
    double s1 = 1.0;
   
    // first pt, normal and tangent at this value of 's' on the curve
    V3d pc1 = profileXYZ(s1);
    p1.Set(pc1.x*cos(phi1), pc1.y, -pc1.x*sin(phi1));
    V3d nc1 = profileNormal(s1);
    n1.Set(nc1.x*cos(phi1), nc1.y, -nc1.x*sin(phi1));
    n1 = n1.Normalize();
    t1.Set(-nc1.y*cos(phi1), nc1.x, nc1.y*sin(phi1));
    t1 = t1.Normalize();
  
    // next point, same s but next phi
    V3d pc4 = pc1;//profileXYZ(s1);
    p4.Set(pc4.x*cos(phi2), pc4.y, -pc4.x*sin(phi2));
    V3d nc4 = nc1;//profileNormal(s1);
    n4.Set(nc4.x*cos(phi2), nc4.y, -nc4.x*sin(phi2));
    n4 = n4.Normalize();
    t4.Set(-nc4.y*cos(phi2), nc4.x, nc4.y*sin(phi2));
    t4 = t4.Normalize();
      
    for (int i = TRANSVERSE_SECS-1; i > EXCLUDING_APICAL_SEC; i--) {
      // get the next two points of the next transverse section 
      double s2 = double(i) / double(TRANSVERSE_SECS);

      V3d pc2 = profileXYZ(s2);
      p2.Set(pc2.x*cos(phi1), pc2.y, -pc2.x*sin(phi1));
      V3d nc2 = profileNormal(s2);
      n2.Set(nc2.x*cos(phi1), nc2.y, -nc2.x*sin(phi1));
      n2 = n2.Normalize();
      t2.Set(-nc2.y*cos(phi1), nc2.x, nc2.y*sin(phi1));
      t2 = t2.Normalize();
  
      V3d pc3 = pc2;//profileXYZ(s2);
      p3.Set(pc3.x*cos(phi2), pc3.y, -pc3.x*sin(phi2));
      V3d nc3 = nc2;//profileNormal(s2,contourWeights);
      n3.Set(nc3.x*cos(phi2), nc3.y, -nc3.x*sin(phi2));
      n3 = n3.Normalize();
      t3.Set(-nc3.y*cos(phi2), nc3.x, nc3.y*sin(phi2));
      t3 = t3.Normalize();
  
      // draw a quad (normal vector set using turtle 'up' vector, but make sure to set 'heading' to tangent)
      nproduce	SP
  		MoveTo(p1.x,p1.y,p1.z) SetHead(t1.x,t1.y,t1.z,n1.x,n1.y,n1.z) PP
  		MoveTo(p2.x,p2.y,p2.z) SetHead(t2.x,t2.y,t2.z,n2.x,n2.y,n2.z) PP
  		MoveTo(p3.x,p3.y,p3.z) SetHead(t3.x,t3.y,t3.z,n3.x,n3.y,n3.z) PP
  		MoveTo(p4.x,p4.y,p4.z) SetHead(t4.x,t4.y,t4.z,n4.x,n4.y,n4.z) PP 
    		EP;

      // copy the top points of this quad to the bottom points of the next quad 
      p1 = p2;
      n1 = n2;
      t1 = t2;
      p4 = p3;
      n4 = n3;
      t4 = t3;
    }

    // draw a triangle as the last element (get apex point at tip [s==1, phi==0])
    pc1 = profileXYZ(0.);
    p1.Set(pc1.x, pc1.y, 0.);
    nc1 = profileNormal(0.);
    n1.Set(nc1.x, nc1.y, 0.);
    n1 = n1.Normalize();
    t1.Set(-nc1.y, nc1.x, 0.);
    t1 = t1.Normalize();

    nproduce	SP
  		MoveTo(p1.x,p1.y,p1.z) SetHead(t1.x,t1.y,t1.z,n1.x,n1.y,n1.z) PP
  		MoveTo(p3.x,p3.y,p3.z) SetHead(t3.x,t3.y,t3.z,n3.x,n3.y,n3.z) PP
  		MoveTo(p2.x,p2.y,p2.z) SetHead(t2.x,t2.y,t2.z,n2.x,n2.y,n2.z) PP
		EP;

  }
  nproduce EB;

  produce;
}

