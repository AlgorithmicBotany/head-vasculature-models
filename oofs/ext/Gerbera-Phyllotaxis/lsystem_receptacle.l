/* ********************************************************************
 * Copyright (c) 2019 Biological Modeling and Visualization Laboratory,
 * University of Calgary. All rights reserved.
 * ********************************************************************/

#include <cmath>
#include <lpfgall.h>
#include "bspline.h"
#include "markers.h"

const double DT = 0.01;
double Progress;		// progress of development (from 0 to 1)
double L;

const double x_scale = 5.0; // scale of x-axis in graph of B-spline basis functions
const double SPHERE_RADIUS = 0.01;

const double CURVE_SCALE = 3.0;

// don't forget to remove the knots file, if NUM_CNTRS changes
const int NUM_CNTRS = 11; // number of contours


V3d profileXYZ(double ns) // C(t,s) in description
{
    extern BSplineBasisFunc contourWeights;
	V3d res(0,0,0);
	for (int i = 1; i <= NUM_CNTRS; i++) {
		res += contourWeights.valueAtFunc(i-1) * ((V3d) curveXYZ(i,ns)); // curveXYZ(k,s) is C_k(s) in description
	}
	
	return CURVE_SCALE * res;
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

// interpolation functions for active ring
double active_ring_newNS()  // s^{gen}(t) in description
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
	return profileXYZ(new_ns);
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

V3d marker_profileXYZ(int markerIndex) // trajectory of landmark. not in description
{
	double new_ns = marker_newNS(markerIndex);
	return profileXYZ(new_ns);
}

double CurveLength()
{
	const double nds = 0.0005; // normalized arc length used for integrating the generating curve length
	double AccumulatedLength = 0.;
	double ns = 0.;
	for (ns=0; ns<1.0-nds; ns+=nds)
  		AccumulatedLength +=
			 (hatProfileXYZ(ns+nds) - hatProfileXYZ(ns)).Length();
	AccumulatedLength += (hatProfileXYZ(1.) - hatProfileXYZ(ns)).Length();
	return AccumulatedLength;
}

vector<double> nsSamples;
void resetNSSamples (void) {

    extern double L;

    double minLen = L*0.005;

    nsSamples.clear();

	const double nds = 0.0001; 
	double ns = 1.; // ns=1 is at bottom of receptacle
    while (ns > 0.01) {
	  double AccumulatedLength = 0.;
      while (AccumulatedLength < minLen && ns > 0.) {
  		  AccumulatedLength += (profileXYZ(ns-nds) - profileXYZ(ns)).Length();
          ns -= nds;
      }
      nsSamples.push_back(ns);
    }

    return;
}

#define Wait 1
#define SelectOrAdd 2
#define Move 3
#define TransferInfo 4
#define WriteInfo 5
#define Init 6
int Phase;	// simulation phase
MouseStatus ms;
int AlreadySelected = 0;
int PlayForward = 0;
FILE *posKnotFile = NULL;

// the degree of contourWeights and markerWeights has to be the same because they use the same knot values
const int CONTOUR_WEIGHTS_SPLINE_DEGREE = 2;
BSplineBasisFunc contourWeights(NUM_CNTRS, CONTOUR_WEIGHTS_SPLINE_DEGREE, 1);
// if MARKER_WEIGHTS_SPLINE_DEGREE != CONTOUR_WEIGHTS_SPLINE_DEGREE, markerWeights uses uniform knot values
const int MARKER_WEIGHTS_SPLINE_DEGREE = 2;
BSplineBasisFunc markerWeights(NUM_CNTRS, MARKER_WEIGHTS_SPLINE_DEGREE, 1);

void resetWeights (double alpha) {
	contourWeights.computeFuncs(alpha);
	markerWeights.computeFuncs(alpha);
}

std::vector<double> knots;

void updateWeights (void) {
// called to generate samples of the contour weights (basis functions)
	contourWeights.SetKnotSeq(knots);
	markerWeights.SetKnotSeq(knots);

	contourWeights.generateSamples();
}


BSplineBasisFunc primordiaMarkerWeights(0,MARKER_WEIGHTS_SPLINE_DEGREE,1);

vector<double> activeRingMarkers;
vector<vector<double>> primordiaMarkers;

const int NUMBER_OF_BRACT_MARKERS = 5;
const int NUMBER_OF_TRANS_MARKERS = 5;

struct primordium{
	double ns;	// normalized curvilinear coordinate
	double phi;	// angular coordinate
	double h;	// vertical cylindrical coordinate
	double r;	// radial cylindrical cordinate
	V3d pos;	// Cartesian position

	void SyncCoordinates(void) {
		V3d pt = profileXYZ(ns);
		r = pt.x;
		h = pt.y;
		pos.x = r*cos(phi); // convert to cartesian coordinates
		pos.z = r*sin(phi);
		pos.y = h;
	}
};
primordium tp;	// trial primordium

Start: {
	static bool firstTime = true;
	// scale the first contour by 19% (estimated by overlapping image sent by Teng)
	if (firstTime) {
		curveScale(1, 0.19, 0.19, 0.19);
		firstTime = false;
	}

	Progress = 0.0;
    resetWeights(0.);

	UseView(View1);
	UseView(View2);
	UseView(View3);

	Phase = Init;
	AlreadySelected = 0;
	PlayForward = 1;
}

StartEach: 
{
	ms = GetMouseStatus();
	switch (Phase) {
	case Init:
		break;
	case Wait: 
		if (ms.lbPushed) {
			Phase = SelectOrAdd;
			AlreadySelected = 0;
		}
		break;
	case SelectOrAdd:
		Phase = Move;
		break;
	case Move:
		if (ms.lbReleased) {
			Phase = TransferInfo;
			posKnotFile = fopen("knotPositions.dat", "w");
		}
		break;
	case TransferInfo:
		Phase = WriteInfo;
		break;
	case WriteInfo:
		Phase = Wait;
		fclose(posKnotFile);
		break;
	}
	UseGroup(Phase);
}	

EndEach: {
	if (PlayForward && Progress < 1.0) {
		resetWeights(Progress);
		Progress = std::min (Progress + DT, 1.);

		L = CurveLength();

        resetNSSamples();

		// Advance position of the active ring
		tp.ns = active_ring_newNS();
		tp.SyncCoordinates();

	} else if (PlayForward && Progress >= 1.0) {
		Progress = 1.; // Was set to 0, which restarts at the beginning
		resetWeights(Progress);

		L = CurveLength();

        resetNSSamples();

		// Advance position of the active ring
		tp.ns = active_ring_newNS();
		tp.SyncCoordinates();
	}
}

module Recompute();
module Draw();
module DrawReceptacle();
module PosKnot(V3d, int, int);
module ProgressBar(double,int);
module PlayButton();

derivation length: 1;

axiom: Recompute() Draw()
	;


production:

group Init:

Recompute() :
{
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


	// set position of active ring
	resetWeights(Progress);
	L = CurveLength();
	tp.ns = active_ring_newNS();
	tp.SyncCoordinates();

	produce;
}


Draw():
{
	if ((posKnotFile = fopen("knotPositions.dat", "r")) == NULL) {
		contourWeights.CopyKnotSeq(knots);
		for (int knot = 0; knot < knots.size(); knot++) {
			nproduce PosKnot(V3d(knots[knot]*x_scale,-0.1,0),0,knot);
  		}
	} else {
		int knot;
		double x,y,z;
		knots.clear();
		while (!feof(posKnotFile)) {
			if (fscanf(posKnotFile, "%d: %lg %lg %lg", &knot, &x, &y, &z) == 4) {
				knots.push_back(x);
				nproduce PosKnot(V3d(x*x_scale,y,z),0,knot);
			}
		}
		fclose(posKnotFile);
		if (!contourWeights.SetKnotSeq(knots)) {
			contourWeights.CopyKnotSeq(knots);
			for (int knot = 0; knot < knots.size(); knot++) {
				nproduce PosKnot(V3d(knots[knot]*x_scale,-0.1,0),0,knot);
  			}
		}
	}
	updateWeights();

	nproduce ProgressBar(0.,0) PlayButton;

	Phase = Wait;
	produce Draw();
}

group Wait:

ProgressBar(u,selected) :
{
	if (PlayForward)
		u = Progress * x_scale;
	produce ProgressBar(u,selected);
}


group SelectOrAdd:

MouseIns() PosKnot(pos, selected, s) :
{
	if (AlreadySelected)
		produce PosKnot(pos,0,s);

	AlreadySelected = 1;
	produce PosKnot(pos, 1, s) ;
}

MouseIns() ProgressBar(u,selected) :
{
	if (AlreadySelected)
		produce ProgressBar(u,0);

	AlreadySelected = 1;
	PlayForward = 0;
	produce ProgressBar(u,1);
}

MouseIns() PlayButton :
{
	if (AlreadySelected)
		produce PlayButton;

	AlreadySelected = 1;
	if (PlayForward)
		PlayForward = 0;
	else
		PlayForward = 1;
	produce PlayButton;
}

group Move:

PosKnot(pos, selected, s) :
{
	if (selected) {

		double dx = pos.x - ms.atMiddle.x;

		knots[s] -= dx*0.1;
		// clamp in [0,1]
		if (knots[s] < 0.0)
		  knots[s] = 0.0;
		if (knots[s] > 1.0)
		  knots[s] = 1.0;

		// clamp to prev or next knot
		if (s+1 < knots.size() && knots[s] > knots[s+1])
		  knots[s] = knots[s+1];
		if (s-1 > 0 && knots[s] < knots[s-1])
		  knots[s] = knots[s-1];

		// update position
		V3d pt = V3d(knots[s]*x_scale,-0.1,0);

		updateWeights();

		produce PosKnot(pt, 1, s); 
        }
}

ProgressBar(u,selected) :
{
	if (selected) {

		double dx = u - ms.atMiddle.x;

		u -= dx*0.5;
		if (u < 0.0)
			u = 0.0;
		if (u > 1.0*x_scale)
			u = 1.0*x_scale;

		Progress = u/x_scale;

		resetWeights(Progress);
		L = CurveLength();
		tp.ns = active_ring_newNS();
		tp.SyncCoordinates();

		produce ProgressBar(u, 1);

	}
}

group TransferInfo:

PosKnot(pos, selected, s) :
{
	fprintf(posKnotFile, "%d: %g %g %g\n", s, pos.x/x_scale, pos.y, pos.z);
	if (selected) {
		produce PosKnot(pos, 0, s) ;
	}
}

ProgressBar(u, selected) :
{
	if (selected) {
		produce ProgressBar(u,0);
	}
}

group 0:

interpretation:
maximum depth: 2;

vgroup View3:

module DrawBand(double,double);

Draw():
{
	resetWeights(Progress);

	// Draw receptacle
	nproduce SB SetColor(4) DrawReceptacle() EB;

	// Draw active ring
	nproduce SetColor(16) SetCoordinateSystem(1.005) DrawBand(tp.ns,0.02);

	produce;
}

DrawBand(ns,w): {

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
    double s1 = ns;
   
    // first pt, normal and tangent at this value of 's' on the curve
    V3d pc1 = profileXYZ(s1);
    p1.Set(pc1.x*cos(phi1), pc1.y, -pc1.x*sin(phi1));
    V3d nc1 = profileNormal(s1);
    n1.Set(nc1.x*cos(phi1), nc1.y, -nc1.x*sin(phi1));
    n1 = -n1.Normalize();
    t1.Set(-nc1.y*cos(phi1), nc1.x, nc1.y*sin(phi1));
    t1 = -t1.Normalize();
  
    // next point, same s but next phi
    V3d pc4 = pc1;
    p4.Set(pc4.x*cos(phi2), pc4.y, -pc4.x*sin(phi2));
    V3d nc4 = nc1;
    n4.Set(nc4.x*cos(phi2), nc4.y, -nc4.x*sin(phi2));
    n4 = -n4.Normalize();
    t4.Set(-nc4.y*cos(phi2), nc4.x, nc4.y*sin(phi2));
    t4 = -t4.Normalize();
      
    const double ds = 0.001;
    double AccumulatedLength = (profileXYZ(s1+ds) - profileXYZ(s1)).Length();
    double s = s1 + ds; 
    while (AccumulatedLength < w && s < 1.) {
      // get the next two points of the next transverse section 
      double s2 = s;

      V3d pc2 = profileXYZ(s2);
      p2.Set(pc2.x*cos(phi1), pc2.y, -pc2.x*sin(phi1));
      V3d nc2 = profileNormal(s2);
      n2.Set(nc2.x*cos(phi1), nc2.y, -nc2.x*sin(phi1));
      n2 = -n2.Normalize();
      t2.Set(-nc2.y*cos(phi1), nc2.x, nc2.y*sin(phi1));
      t2 = -t2.Normalize();
  
      V3d pc3 = pc2;
      p3.Set(pc3.x*cos(phi2), pc3.y, -pc3.x*sin(phi2));
      V3d nc3 = nc2;
      n3.Set(nc3.x*cos(phi2), nc3.y, -nc3.x*sin(phi2));
      n3 = -n3.Normalize();
      t3.Set(-nc3.y*cos(phi2), nc3.x, nc3.y*sin(phi2));
      t3 = -t3.Normalize();

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

      // move the arclength coordinate
      AccumulatedLength += (profileXYZ(s+ds) - profileXYZ(s)).Length();
      s += ds;

    }
  }
  nproduce EB;

  produce;
}

DrawReceptacle() : {

  // the receptacle can be drawn 

  // here we draw the surface of revolution using polygons, SP ... EP --- no calls to StartGC / EndGC 

  const int LONGITUDINAL_SECS = 80;

  nproduce SB;

  for (int j = 0; j < LONGITUDINAL_SECS; ++j) {
    double phi1 = double(j) / double(LONGITUDINAL_SECS) * 2. * M_PI;
    double phi2 = double(j+1) / double(LONGITUDINAL_SECS) * 2. * M_PI;

    V3d p1, n1, t1;
    V3d p2, n2, t2;
    V3d p3, n3, t3;
    V3d p4, n4, t4;

    // for first transverse section, compute the bottom two points
    double s1 = 1.;
   
    // first pt, normal and tangent at this value of 's' on the curve
    V3d pc1 = profileXYZ(s1);
    p1.Set(pc1.x*cos(phi1), pc1.y, -pc1.x*sin(phi1));
    V3d nc1 = profileNormal(s1);
    n1.Set(nc1.x*cos(phi1), nc1.y, -nc1.x*sin(phi1));
    n1 = n1.Normalize();
    t1.Set(-nc1.y*cos(phi1), nc1.x, nc1.y*sin(phi1));
    t1 = t1.Normalize();
  
    // next point, same s but next phi
    V3d pc4 = pc1;
    p4.Set(pc4.x*cos(phi2), pc4.y, -pc4.x*sin(phi2));
    V3d nc4 = nc1;
    n4.Set(nc4.x*cos(phi2), nc4.y, -nc4.x*sin(phi2));
    n4 = n4.Normalize();
    t4.Set(-nc4.y*cos(phi2), nc4.x, nc4.y*sin(phi2));
    t4 = t4.Normalize();


    for (int i = 0; i < nsSamples.size() ; ++i) {

      double s2 = nsSamples[i];

      V3d pc2 = profileXYZ(s2);
      p2.Set(pc2.x*cos(phi1), pc2.y, -pc2.x*sin(phi1));
      V3d nc2 = profileNormal(s2);
      n2.Set(nc2.x*cos(phi1), nc2.y, -nc2.x*sin(phi1));
      n2 = n2.Normalize();
      t2.Set(-nc2.y*cos(phi1), nc2.x, nc2.y*sin(phi1));
      t2 = t2.Normalize();
  
      V3d pc3 = pc2;
      p3.Set(pc3.x*cos(phi2), pc3.y, -pc3.x*sin(phi2));
      V3d nc3 = nc2;
      n3.Set(nc3.x*cos(phi2), nc3.y, -nc3.x*sin(phi2));
      n3 = n3.Normalize();
      t3.Set(-nc3.y*cos(phi2), nc3.x, nc3.y*sin(phi2));
      t3 = t3.Normalize();

      int colour = 5; // disc
      if (s2 >= marker_newNS(1))
        colour = 7; // stem
      if (s2 < marker_newNS(1) && s2 >= marker_newNS(NUMBER_OF_BRACT_MARKERS))
        colour = 6; // bract
      if (s2 < marker_newNS(NUMBER_OF_BRACT_MARKERS) && s2 >= marker_newNS(NUMBER_OF_BRACT_MARKERS+NUMBER_OF_TRANS_MARKERS))
        colour = 9; // trans
      nproduce SetColor(colour);

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

    nproduce SP
 		    MoveTo(p1.x,p1.y,p1.z) SetHead(t1.x,t1.y,t1.z,n1.x,n1.y,n1.z) PP
  		    MoveTo(p3.x,p3.y,p3.z) SetHead(t3.x,t3.y,t3.z,n3.x,n3.y,n3.z) PP
  		    MoveTo(p2.x,p2.y,p2.z) SetHead(t2.x,t2.y,t2.z,n2.x,n2.y,n2.z) PP
		    EP;

  }
  nproduce EB;

  produce;
}

vgroup View2:

Draw():
{
  const double LINE_SCALE = 1.2;

  nproduce MoveTo(0,0,0) SetColor(4) SetWidth(0.015*LINE_SCALE) ;

  // W axis
  nproduce SB /*Right(90) f(-0.0) Left(90)*/ Down(90) f(-0.9) Up(90) Sphere(0.025*0.5) F(1.0) 
            F(0.05) SetWidth(0.025*3.*LINE_SCALE) StartGC SetWidth(0.025*0.1) F(0.1) EndGC EB;

  // U axis
  nproduce SB Down(90) f(-0.9) Up(90) Right(90) F(x_scale)
              F(0.05) SetWidth(0.025*3.*LINE_SCALE) StartGC SetWidth(0.025*0.1) F(0.1) EndGC EB;

  nproduce SB SetColor(255) MoveTo(x_scale*0.45,-0.25,0.) Label("time") EB;
  nproduce SB SetColor(255) MoveTo(-.35,0.5,0.) Label("B(t)") EB;

  // compute the basis functions
  const int Usteps = 20 * NUM_CNTRS;

  // show the basis functions
  nproduce SetWidth(0.02*LINE_SCALE) SB;

  for (int i = 0; i < NUM_CNTRS; i++) {
    
    nproduce SetColor(i+32) MoveTo3d(V3d(0.0,contourWeights.sampleAt(0,i),double(i)/double(NUM_CNTRS)));
    nproduce TurnAround() StartGC;
    for (int ui = 0; ui <= Usteps; ui++) {
      double u = ((double)ui)/((double)Usteps);

      if (u > knots[CONTOUR_WEIGHTS_SPLINE_DEGREE] &&
          u <= knots[knots.size()-1-CONTOUR_WEIGHTS_SPLINE_DEGREE])
        nproduce LineTo3d(V3d(u*x_scale,contourWeights.sampleAt(ui,i),double(i)/double(2.*NUM_CNTRS)));
      else
        nproduce MoveTo3d(V3d(u*x_scale,contourWeights.sampleAt(ui,i),double(i)/double(2.*NUM_CNTRS)));
    }
    nproduce EndGC;
  }

  produce EB;
}

PosKnot(pos, selected, scount) : {
  produce MoveTo3d(V3d(knots[scount]*x_scale,-0.06,0)) SetColor(1) Sphere(0.025);
}

ProgressBar(u, selected) : {
  if (PlayForward)
    u = Progress * x_scale;
  produce SB SetColor(5) SetWidth(0.03) MoveTo(u,0,0.9) StartGC F(0.95) EndGC Sphere(0.05) EB;
}

PlayButton() : {
produce;
  const double size = 0.35;
  if (PlayForward)
    produce SB SetColor(16) MoveTo(x_scale*0.5+size*0.5,-0.5,0) Sphere(size*0.5)
               SetColor(255) MoveTo(x_scale*0.5+size,-0.5,0) Label("STOP") EB;
  else
    produce SB SetColor(17) MoveTo(x_scale*0.5+size*0.5,-0.5,0) Sphere(size*0.5)
               SetColor(255) MoveTo(x_scale*0.5+size,-0.5,0) Label("PLAY") EB;
}

vgroup View1:

module DrawBackground();
module DrawData();

const double dns = 0.005; 	// normalized arc-length distance betwen consecutive points
const double LINE_WIDTH	= 0.015;
const double DALPHA = 0.02; 	// interpolation step

Draw():
{
	nproduce SetWidth(LINE_WIDTH*0.5) SetColor(2) ;
	double alpha = 0.;
	while (alpha < 1.) {
		resetWeights(alpha);
        nproduce MoveTo3d(profileXYZ(0.));
		for (double ns= dns; ns<1.; ns+= dns) {
			nproduce LineTo3d(profileXYZ(ns));
        }
        nproduce LineTo3d(profileXYZ(1.));
		alpha += DALPHA;
	}

	// draw approximating contour
	resetWeights(Progress);
	nproduce SetWidth(LINE_WIDTH*2.) SetColor(5) ;
	nproduce SB MoveTo3d(profileXYZ(0.)) Right(90) StartGC ;
	for (double s = dns; s < 1.; s += dns) {
	  	nproduce LineTo3d(profileXYZ(s));
	}
	nproduce LineTo3d(profileXYZ(1.)) EndGC EB;

	// draw lines connecting markers across contours
	nproduce SetWidth(LINE_WIDTH*0.5);
	int colour = 10;
	for (int i = 0; i < primordiaMarkers.size(); ++i) {

        int numMarkers = primordiaMarkers[i].size();

        // draw starting point
        resetWeights(0.);
        V3d v = marker_profileXYZ(i);
        nproduce SetColor(colour) MoveTo3d(v);

		// line to next marker (increment alpha by dalpha)
        // change colour when the line intersects the active ring line
		double nsAR = myCurveMarker(1,myCurveMarkers(1));
		double alpha = 0.;
		while (alpha < 1.) {
			resetWeights(alpha);
			// change colour once marker passes the active ring
			double new_nsAR = active_ring_newNS();	
			if (marker_newNS(i) >= new_nsAR)
				nproduce SetColor(colour);
			else
				nproduce SetColor(colour+1);

			v = marker_profileXYZ(i);
			nproduce LineTo3d(v);
			alpha += DALPHA;
		}

        // draw ending point
        resetWeights(1.);
        v = marker_profileXYZ(i);
        nproduce LineTo3d(v);
	}

	// draw the interpolated active ring marker
	nproduce SetWidth(LINE_WIDTH*1.25);

    resetWeights(0.);
    V3d v = active_ring_profileXYZ();
	nproduce SetColor(16) MoveTo3d(v);

	for (double alpha = 0.; alpha < 1.; alpha += DALPHA) {
		resetWeights(alpha);
        v = active_ring_profileXYZ();
		nproduce LineTo3d(v);
	}

	resetWeights(1.);
    v = active_ring_profileXYZ();
	nproduce LineTo3d(v);

	produce DrawBackground() DrawData() ;
}

DrawData() :
{
	nproduce SetWidth(LINE_WIDTH*1.1) ;

	for (int cntr = 1; cntr <= NUM_CNTRS; ++cntr) {
	
		nproduce SetColor(cntr+32-1) MoveTo3d(CURVE_SCALE * (V3d) curveXYZ(cntr,0.));
		for (double ns=dns; ns<1.; ns+= dns)
			nproduce LineTo3d(CURVE_SCALE * (V3d) curveXYZ(cntr,ns) );
		nproduce LineTo3d(CURVE_SCALE * (V3d) curveXYZ(cntr,1.));

		// draw marker for the bracts on curve
		int numMarkers = myCurveMarkers(cntr);
		for (int i = 0; i < numMarkers-1; i++) {
			double ns = myCurveMarker(cntr,i+1);
			V3d v = CURVE_SCALE * (V3d) curveXYZ(cntr,ns);
			nproduce MoveTo3d(v) Sphere(SPHERE_RADIUS*1.4);
		}

		// draw the active ring marker on curve
		double ns = myCurveMarker(cntr,numMarkers);
		V3d v = CURVE_SCALE * (V3d) curveXYZ(cntr,ns);
		nproduce SetColor(3) MoveTo3d(v);// Sphere(SPHERE_RADIUS*1.5);
	}

	produce;
}

DrawBackground():
{
	const double SCALE = 0.0034;
	produce SB SetColor(12) MoveTo(-1293.*SCALE*0.5,0.,-0.01) CurrentTexture(0) Surface(0,SCALE) CurrentTexture(-1) EB;
}

