#ifndef __MARKERS_H
#define __MARKERS_H
double __dataMarkers[7][40] = 
{
{
2,
1,
0.914578,
},
{
2,
0.7615,
0.670094,
},
{
3,
0.630337,
0.56343,
0.468665,
},
{
6,
0.673937,
0.571516,
0.508164,
0.453069,
0.392259,
0.081909,
},
{
7,
0.718675,
0.600347,
0.529835,
0.474592,
0.416086,
0.371114,
0,
},
{
17,
0.829094,
0.724934,
0.653097,
0.598898,
0.542105,
0.497668,
0.449712,
0.401992,
0.357358,
0.311463,
0.265128,
0.220096,
0.185779,
0.150545,
0.118022,
0.091229,
0,
},
{
18,
0.883139,
0.800179,
0.719477,
0.648791,
0.582521,
0.517491,
0.469626,
0.424017,
0.37905,
0.332558,
0.289275,
0.247404,
0.204589,
0.164206,
0.1241,
0.087006,
0.046487,
0,
}
};
double myCurveMarkers (int cntr) {
  return __dataMarkers[cntr-1][0];
}
double myCurveMarker (int cntr, int i) {
  return __dataMarkers[cntr-1][i];
}
#endif
