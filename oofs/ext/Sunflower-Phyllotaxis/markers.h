#ifndef __MARKERS_H
#define __MARKERS_H
double __dataMarkers[9][60] = 
{
{
2,
0.836909,
0.567023,
},
{
2,
0.835239,
0.494324,
},
{
3,
0.846615,
0.633988,
0.468451,
},
{
5,
0.859463,
0.632225,
0.535581,
0.451522,
0.398903,
},
{
9,
0.855477,
0.670623,
0.567035,
0.502535,
0.464758,
0.428724,
0.397661,
0.367582,
0.245054,
},
{
12,
0.880082,
0.712965,
0.614735,
0.543259,
0.488164,
0.449589,
0.418034,
0.385316,
0.346239,
0.308837,
0.281746,
0.178223,
},
{
16,
0.896269,
0.696706,
0.594729,
0.536679,
0.494978,
0.464551,
0.436296,
0.412363,
0.384145,
0.356314,
0.330376,
0.308447,
0.284154,
0.26366,
0.243672,
0.07746,
},
{
18,
0.896896,
0.711629,
0.601219,
0.544785,
0.499007,
0.465108,
0.436669,
0.416782,
0.394813,
0.375233,
0.352944,
0.336111,
0.317165,
0.298173,
0.277293,
0.256107,
0.235751,
0,
},
{
35,
0.813514,
0.622846,
0.552563,
0.513998,
0.488685,
0.469351,
0.448156,
0.415965,
0.390532,
0.363907,
0.33938,
0.318315,
0.29949,
0.277658,
0.259573,
0.243878,
0.227012,
0.211113,
0.194409,
0.179966,
0.164354,
0.149779,
0.137111,
0.125838,
0.114164,
0.103478,
0.091338,
0.079387,
0.06722,
0.055304,
0.043207,
0.032404,
0.021583,
0.010603,
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
