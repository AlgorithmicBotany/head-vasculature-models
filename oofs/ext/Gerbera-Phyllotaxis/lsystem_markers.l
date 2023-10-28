/* ********************************************************************
 * Copyright (c) 2019 Biological Modeling and Visualization Laboratory,
 * University of Calgary. All rights reserved.
 * ********************************************************************/

#include <cmath>
#include <cstdio>
#include <lpfgall.h>

const int NUM_CNTRS = 11; // number of contours

Start: {


  FILE *outfile = fopen("markers.h", "w");

  fprintf (outfile, "#ifndef __MARKERS_H\n");
  fprintf (outfile, "#define __MARKERS_H\n"); 

  fprintf (outfile, "double __dataMarkers[%d][%d] = \n", NUM_CNTRS, 40);
  fprintf (outfile, "{\n");

  for (int cntr = 1; cntr <= NUM_CNTRS; cntr++) {

    int numMarkers = curveMarkers(cntr);

    fprintf (outfile, "{\n%d,\n", numMarkers);

    for (int i = 1; i <= numMarkers; i++) {
      fprintf (outfile, "%g,\n", curveMarker(cntr,i));
    }

    fprintf (outfile, "}");
    if (cntr < NUM_CNTRS)
      fprintf (outfile, ",");
    fprintf (outfile, "\n");

  }

  fprintf (outfile, "};\n");



  fprintf (outfile, "double myCurveMarkers (int cntr) {\n"
                    "  return __dataMarkers[cntr-1][0];\n"
                    "}\n");

  fprintf (outfile, "double myCurveMarker (int cntr, int i) {\n"
                    "  return __dataMarkers[cntr-1][i];\n"
                    "}\n");

  fprintf (outfile, "#endif\n");

  Printf ("\nDone generating markers.h\n");
  exit(0);
}


derivation length: 0;

Axiom: ;




