/* ********************************************************************
 * Copyright (c) 2019 Biological Modeling and Visualization Laboratory,
 * University of Calgary. All rights reserved.
 * ********************************************************************/

#include <cmath>
#include <lpfgall.h>
#include "bspline.h"
#include "markers.h"

#define OUTPUT_CONTOURS 0       // set to 1 to output contours. should be 0 when playing with parameters...

#define CURVE_SCALE 	3.0		// scaling factor for the graphically-defined apex profile
#define SPHERE_RADIUS 	0.007	
#define LINE_WIDTH      0.007 

#define dns 			0.001 	// normalized arc-length distance betwen consecutive points
#if OUTPUT_CONTOURS==1
#define DALPHA			0.01*0.32 	// interpolation step
#else
#define DALPHA          0.01
#endif


const int NUM_CNTRS = 7; // number of contours

/* Interpolate between the contour curves */
// Specifically, given:
// - graphically defined profile (contour) curves,
//   both parametrized by normalized arc-lenth coordinate ns;
// - a weighting of the curves 
// - a value of coordinate ns in [0,1]
// return position of the interpolated point.

V3d hatProfileXYZ(double ns)
{
    extern BSplineBasisFunc contourWeights;
    extern BSplineBasisFunc primordiaMarkerWeights;
    extern vector<vector<double>> primordiaMarkers;
    
    primordiaMarkerWeights.computeFuncs(ns);

    V3d v(0,0,0);
    for (int i = 1; i <= NUM_CNTRS; i++) {
        double newNS = 0.;
        for (int j = 0; j < primordiaMarkers.size(); j++) {
            newNS += primordiaMarkerWeights.valueAtFunc(j) * primordiaMarkers[j][i-1];
        }

		v += contourWeights.valueAtFunc(i-1) * ((V3d) curveXYZ(i,newNS));
	}
    return CURVE_SCALE * v;
}

V3d pure_profileXYZ(double ns) // C(t,s) in description
{
    extern BSplineBasisFunc contourWeights;
	V3d res(0,0,0);
	for (int i = 1; i <= NUM_CNTRS; i++) {
		res += contourWeights.valueAtFunc(i-1) * ((V3d) curveXYZ(i,ns)); // curveXYZ(k,s) is C_k(s) in description
	}
	
	return CURVE_SCALE * res;
}


double newNS(int markerIndex)
{
    extern BSplineBasisFunc markerWeights;
	extern vector<vector<double>> primordiaMarkers;
	double new_ns = 0.;	
	for (int i = 0; i < NUM_CNTRS; i++) {
		new_ns += markerWeights.valueAtFunc(i) * primordiaMarkers[markerIndex][i];
	}	
	return new_ns;
}

V3d profileXYZ(int markerIndex)
{
	double new_ns = newNS(markerIndex);
	return pure_profileXYZ(new_ns);
}


// interpolation functions for active ring
double active_ring_newNS() // s^{gen}(t) in description
{
    extern BSplineBasisFunc markerWeights;
    extern vector<double> activeRingMarkers;
    double new_ns = 0.;
    for (int i = 0; i < activeRingMarkers.size(); i++) {
        new_ns += markerWeights.valueAtFunc(i) * activeRingMarkers[i];
    }
	return new_ns;
}

V3d active_ring_profileXYZ() // G(t) in description
{
	double new_ns = active_ring_newNS();	
	return pure_profileXYZ(new_ns);
}


module Recompute();	
module DrawAll(double);
module DrawData();
module DrawBackground();
module OutputContours();

derivation length: 2;

Axiom:	Recompute();

const int CONTOUR_WEIGHTS_SPLINE_DEGREE = 2;
BSplineBasisFunc contourWeights(NUM_CNTRS, CONTOUR_WEIGHTS_SPLINE_DEGREE, 1);

const int MARKER_WEIGHTS_SPLINE_DEGREE = 2;
BSplineBasisFunc markerWeights(NUM_CNTRS, MARKER_WEIGHTS_SPLINE_DEGREE, 1);

void resetWeights (double alpha) {
	contourWeights.computeFuncs(alpha);
	markerWeights.computeFuncs(alpha);
}

BSplineBasisFunc primordiaMarkerWeights(0,MARKER_WEIGHTS_SPLINE_DEGREE,1);

vector<double> activeRingMarkers;
vector<vector<double>> primordiaMarkers;

production:

Recompute() :
{
	Printf ("Recompute...");

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

	Printf ("Done.\n");

#if OUTPUT_CONTOURS
    nproduce OutputContours();
#endif

	produce DrawBackground() DrawData() DrawAll(DALPHA);
}

DrawAll(dalpha) :
{
	nproduce SetWidth(LINE_WIDTH*2.) SetColor(1);

	// draw interpolated contours
	double alpha = 0.;
	while (alpha < 1.) {
		resetWeights(alpha);
        nproduce MoveTo3d(hatProfileXYZ(0.));
		for (double ns= dns; ns<1.; ns+= dns) {
			nproduce LineTo3d(hatProfileXYZ(ns));
        }
        nproduce LineTo3d(hatProfileXYZ(1.));
		alpha += 30.0 * dalpha * func(output_dalpha,alpha);
	}

	// draw lines connecting markers across contours
	nproduce SetWidth(SPHERE_RADIUS*0.8 * 2.);
	int colour = 2;
	for (int i = 0; i < primordiaMarkers.size(); ++i) {

        int numMarkers = primordiaMarkers[i].size();

        // draw starting point
        resetWeights(0.);
        V3d v = profileXYZ(i) + V3d(0.,0.,0.01);
        nproduce SetColor(colour) MoveTo3d(v) ;

		// line to next marker (increment alpha by dalpha)
        // change colour when the line intersects the active ring line
		double nsAR = myCurveMarker(1,myCurveMarkers(1));
		double alpha = 0.;
		while (alpha < 1.) {
			resetWeights(alpha);
			// change colour once marker passes the active ring
			double new_nsAR = active_ring_newNS();	
			if (newNS(i) >= new_nsAR)
				nproduce SetColor(colour);
			else
				nproduce SetColor(colour+5);

			v = profileXYZ(i) + V3d(0.,0.,0.01);
			nproduce LineTo3d(v) ;
			alpha += 1.0 * dalpha * func(output_dalpha,alpha);
		}

        // draw ending point
        resetWeights(1.);
        v = profileXYZ(i) + V3d(0.,0.,0.01);
        nproduce LineTo3d(v) ;//Sphere(SPHERE_RADIUS*0.8);
	}

	// draw the interpolated active ring marker
	nproduce SetWidth(LINE_WIDTH*1.3 * 2.);

    resetWeights(0.);
    V3d v = active_ring_profileXYZ() + V3d(0.,0.,0.015);
	nproduce SetColor(3) MoveTo3d(v) ;//Sphere(SPHERE_RADIUS*1.3);

	for (double alpha = 0.; alpha < 1.; alpha += dalpha) {
		resetWeights(alpha);
        v = active_ring_profileXYZ() + V3d(0.,0.,0.015);
		nproduce LineTo3d(v) ;//Sphere(SPHERE_RADIUS*1.3);
	}

	resetWeights(1.);
    v = active_ring_profileXYZ() + V3d(0.,0.,0.015);
	nproduce LineTo3d(v) ;//Sphere(SPHERE_RADIUS*1.3);

 	produce;
}

DrawData() :
{
	nproduce SetWidth(LINE_WIDTH*1.1*2.) ;

	for (int cntr = 1; cntr <= NUM_CNTRS; ++cntr) {
	
		nproduce SB SetColor(4) MoveTo3d(CURVE_SCALE * (V3d) curveXYZ(cntr,0.)) Right(90) StartGC;
		for (double ns=dns; ns<1.; ns+= dns)
			nproduce LineTo3d(CURVE_SCALE * (V3d) curveXYZ(cntr,ns) );
		nproduce LineTo3d(CURVE_SCALE * (V3d) curveXYZ(cntr,1.)) EndGC EB;

		// draw marker for the bracts on curve
		nproduce SetColor(2);
		int numMarkers = myCurveMarkers(cntr);
		for (int i = 0; i < numMarkers-1; i++) {
			double ns = myCurveMarker(cntr,i+1);
			V3d v = CURVE_SCALE * (V3d) curveXYZ(cntr,ns) + V3d(0.,0.,0.01);
			nproduce MoveTo3d(v) ;//Sphere(SPHERE_RADIUS*1.4 * 1.5);
		}

		// draw the active ring marker on curve
		double ns = myCurveMarker(cntr,numMarkers);
		V3d v = CURVE_SCALE * (V3d) curveXYZ(cntr,ns);
		nproduce SetColor(6) MoveTo3d(v) ;
	}

	produce;
}

DrawBackground():
{
	float w = 0.86*CURVE_SCALE; 
	float h = 0.86*CURVE_SCALE; 
	nproduce SetColor(12); 
    nproduce CurrentContour(flat_contour); 
    nproduce CurrentTexture(0); 
    nproduce TextureVCoeff(1.0/h) SetWidth(w);
    nproduce SB MoveTo(0.,0.,0.) StartGC() F(h) EndGC() EB; 
    produce CurrentTexture(-1) CurrentContour(0);
}

OutputContours():
{
	FILE *output_file = NULL;

	if ((output_file = fopen("bellis.cntrs", "w")) == NULL) {
		Printf ("Error: cannot open bellis.cntrs for writing contour animation file\n");
		produce;
	}

    Printf ("Save to file...");

	// output the number of contours that will be written to the file
	double alpha = 0.;
	int num_contours = 0;
	while (alpha < 1.) {
		alpha += DALPHA * func(output_dalpha,alpha);
		++num_contours;
	}
	fprintf (output_file, "%d\n", int(num_contours)+1);


	// set the maximum and minimum number of steps to take between markers
	// that are used in a linear function of the arc length between markers
	const double MAX_STEPS = 6.;
	const double MIN_STEPS = 2.;
	const double len_du = 0.0001; // to compute arc lengths
	double arclen[128], steps[128], steps_du[128];

	// compute arc lengths between markers on final contour
	resetWeights(1.);
	double max_arc_len = -1.;
	double min_arc_len = 1e11;
    int num_markers = primordiaMarkers.size();
	arclen[0] = 0.f;
	for (int i = num_markers-1; i > 0; --i) {
		arclen[i] = 0.;
		for (double u = primordiaMarkers[i][NUM_CNTRS-1]; u < primordiaMarkers[i-1][NUM_CNTRS-1]; u += len_du) {
			arclen[i] +=
			(pure_profileXYZ(u+len_du) - pure_profileXYZ(u)).Length();
		}
		if (arclen[i] > max_arc_len)
			max_arc_len = arclen[i];
		if (arclen[i] < min_arc_len)
			min_arc_len = arclen[i];
	}
	// estimate number of steps to take between markers on the final contour
	int num_points = 0;
	for (int i = num_markers-1; i > 0; --i) {
		steps[i] = (MAX_STEPS-MIN_STEPS)/(max_arc_len-min_arc_len)
			* (arclen[i]-min_arc_len) + MIN_STEPS;
		num_points += steps[i];
	}
	
	fprintf (output_file, "%d\n", num_points+1);
	// output which of the points are landmarks --- 1 is landmark, 0 is not
	for (int i = num_markers-1; i > 0; --i) {
		fprintf (output_file, "1\n");
		for (int j = 0; j < int(steps[i])-1; ++j) {
			fprintf (output_file, "0\n");
		}
	}
	fprintf (output_file, "1\n");

	// output pts along the contour over "time" (changing alpha)
	alpha = 0.;
	while (alpha < 1.) {

		resetWeights(alpha);
		alpha += DALPHA * func(output_dalpha,alpha);

		// compute arc lengths between the markers on the current contour
		arclen[0] = 0.f;
		for (int i = num_markers-1; i > 0; --i) {
			arclen[i] = 0.;
			//for (double u = primordiaMarkers[i][NUM_CNTRS-1]; u < primordiaMarkers[i-1][NUM_CNTRS-1]; u += len_du) {
			for (double u = newNS(i); u < newNS(i-1); u += len_du) {
				arclen[i] += (pure_profileXYZ(u+len_du) - pure_profileXYZ(u)).Length();
			}
			steps_du[i] = arclen[i] / std::floor(steps[i]);
		}

		// output the marker and the points in between
		for (int i = 0; i < num_markers; ++i) {
			//double u = primordiaMarkers[num_markers-1-i][NUM_CNTRS-1];
            double u = newNS(num_markers-1-i);
            V3d pt = pure_profileXYZ(u);
			nproduce MoveTo(pt.x,pt.y,0.0) SetColor(5) Sphere(0.0125);
			fprintf (output_file, "%g %g %g\n", pt.x, pt.y, u);

			// step along curve adding new points
			for (int step = 0; step < int(steps[num_markers-1-i]-1); ++step) {
				// accumulate arc length to find 'u' parameter of next point
				double a = 0.f;
				while (a < steps_du[num_markers-1-i] && u < 1.) {
					a += (pure_profileXYZ(u+len_du) - pure_profileXYZ(u)).Length();
					u += len_du;
				}
				V3d pt = pure_profileXYZ(u);
				nproduce MoveTo(pt.x,pt.y,0.0) SetColor(5) Sphere(0.01);
				fprintf (output_file, "%g %g %g\n", pt.x, pt.y, u);
			}
		}
	}

	// output the final curve at alpha = 1
	alpha = 1.f;
	resetWeights(alpha);

	// compute arc lengths between markers on final curve at alpha = 1
	arclen[0] = 0.f;
	for (int i = num_markers-1; i > 0; --i) {
		arclen[i] = 0.;
		for (double u = primordiaMarkers[i][NUM_CNTRS-1]; u < primordiaMarkers[i-1][NUM_CNTRS-1]; u += len_du) {
			arclen[i] += (pure_profileXYZ(u+len_du) - pure_profileXYZ(u)).Length();
		}
		steps_du[i] = arclen[i] / std::floor(steps[i]);
	}

	// output the pt on the final curve at alpha = 1
	for (int i = 0; i < num_markers; ++i) {
		//double u = primordiaMarkers[num_markers-1-i][NUM_CNTRS-1];
        double u = newNS(num_markers-1-i);
		V3d pt = pure_profileXYZ(u);
		nproduce MoveTo(pt.x,pt.y,0.0) SetColor(5) Sphere(0.02);
		fprintf (output_file, "%g %g %g\n", pt.x, pt.y, u);

		// step along curve adding new points
		for (int step = 0; step < int(steps[num_markers-1-i]-1); ++step) {
			// accumulate arc length to find 'u' parameter of next point
			double a = 0.f;
			while (a < steps_du[num_markers-1-i] && u < 1.) {
				a += (pure_profileXYZ(u+len_du) - pure_profileXYZ(u)).Length();
				u += len_du;
			}
			V3d pt = pure_profileXYZ(u);
			nproduce MoveTo(pt.x,pt.y,0.0) SetColor(5) Sphere(0.01);
			fprintf (output_file, "%g %g %g\n", pt.x, pt.y, u);
		}

	}

	// output the interpolated active ring marker
	double ns = myCurveMarker(1,myCurveMarkers(1));
	alpha = 0.;
	while (alpha < 1.) {
	    resetWeights(alpha);
		alpha += DALPHA * func(output_dalpha,alpha);
		V3d pt = active_ring_profileXYZ();
		nproduce MoveTo(pt.x,pt.y,0.0) SetColor(4) Sphere(0.02);
		//fprintf (output_file, "%g: ", active_ring_newNS(ns, markerWeights));
		fprintf (output_file, "%g %g\n", pt.x, pt.y);
	}
	alpha = 1.f;
	resetWeights(alpha);
	V3d pt = active_ring_profileXYZ();
	nproduce MoveTo(pt.x,pt.y,0.0) SetColor(4) Sphere(0.01);
	//fprintf (output_file, "%g: ", active_ring_newNS(ns, markerWeights));
	fprintf (output_file, "%g %g\n", pt.x, pt.y);


	fclose(output_file);
    Printf ("Done\n");
	produce;
}

