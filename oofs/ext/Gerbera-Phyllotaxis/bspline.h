/* ********************************************************************
 * Copyright (c) 2019 Biological Modeling and Visualization Laboratory,
 * University of Calgary. All rights reserved.
 * ********************************************************************/

#include <vector>

using std::vector;

class BSpline
{
  public:
    BSpline();
    ~BSpline();
    void SetCtrlPoints (const vector<V3d> &ctrl_pts, int degree, int type);
    void AddCtrlPoint (const V3d &new_pt, int index);
    void SetKnotSeq (const vector<double> &new_knots);
    void CopyKnotSeq (vector<double> &dest_knots);
    V3d C (double u);
    V3d dCdu (double u);
    V3d deBoor (double u, int degree, const vector<V3d> &ctrl_pts, const vector<double> &knot_seq);

  private:
    vector<V3d> ctrl_pts;
    vector<V3d> dctrl_pts;
    vector<double> knot_seq;
    vector<double> dknot_seq;
    int degree;
    int type;
};

BSpline::BSpline() {}
BSpline::~BSpline() {}

void BSpline::SetCtrlPoints (const vector<V3d> &new_ctrl_pts, int new_degree, int new_type)
// called to set the control points or change them 
{
  type = new_type;
  
  if (new_degree > 0)
    degree = new_degree;
  else
    degree = 1;

  ctrl_pts.erase(ctrl_pts.begin(),ctrl_pts.end());
  ctrl_pts = new_ctrl_pts;

  // generate the knot vector
  knot_seq.erase(knot_seq.begin(),knot_seq.end());
  int knot_size = degree + ctrl_pts.size() + 1;
  for (int i = 0; i < knot_size; i++)
    if (type == 0) { // open curve (or uniform curve)
      knot_seq.push_back(double(i)/(double)(knot_size-1)); // open curve
    }
    else { // clamped (or open uniform curve)
      if (i < degree+1)
        knot_seq.push_back(0.0);
      else if (i >= knot_size-degree-1)
        knot_seq.push_back(1.0);
      else
        knot_seq.push_back((double) (i-degree) / (double) (knot_size-1-2*degree));
    }

  // generate the ctrl pts for first derivative
  dctrl_pts.erase(dctrl_pts.begin(),dctrl_pts.end());
  for (int i = 0; i < ctrl_pts.size()-1; i++) {
    dctrl_pts.push_back(degree / (knot_seq[i+degree+1] - knot_seq[i+1]) * (ctrl_pts[i+1] - ctrl_pts[i]));
  }
  // the first derivative's knot sequence has the first and last knots removed
  dknot_seq.erase(dknot_seq.begin(),dknot_seq.end());
  dknot_seq = knot_seq;
  dknot_seq.pop_back(); // remove last element
  dknot_seq.erase(dknot_seq.begin()); // remove first element

  return;
}

void BSpline::AddCtrlPoint (const V3d &new_pt, int not_used_index)
{
  int index = 0;
  while (index < ctrl_pts.size() && new_pt.x > ctrl_pts[index].x) {
	  ++index;
  }
  ctrl_pts.insert(ctrl_pts.begin()+index,new_pt);

  // generate the knot vector
  knot_seq.erase(knot_seq.begin(),knot_seq.end());
  int knot_size = degree + ctrl_pts.size() + 1;
  for (int i = 0; i < knot_size; i++)
    if (type == 0) { // open curve (or uniform curve)
      knot_seq.push_back(double(i)/(double)(knot_size-1)); // open curve
    }
    else { // clamped (or open uniform curve)
      if (i < degree+1)
        knot_seq.push_back(0.0);
      else if (i >= knot_size-degree-1)
        knot_seq.push_back(1.0);
      else
        knot_seq.push_back((double) (i-degree) / (double) (knot_size-1-2*degree));
    }

  // generate the ctrl pts for first derivative
  dctrl_pts.erase(dctrl_pts.begin(),dctrl_pts.end());
  for (int i = 0; i < ctrl_pts.size()-1; i++) {
    dctrl_pts.push_back(degree / (knot_seq[i+degree+1] - knot_seq[i+1]) * (ctrl_pts[i+1] - ctrl_pts[i]));
  }
  // the first derivative's knot sequence has the first and last knots removed
  dknot_seq.erase(dknot_seq.begin(),dknot_seq.end());
  dknot_seq = knot_seq;
  dknot_seq.pop_back(); // remove last element
  dknot_seq.erase(dknot_seq.begin()); // remove first element


}

void BSpline::SetKnotSeq (const vector<double> &new_knots)
{
  if (knot_seq.size() != new_knots.size()) {
    Printf ("Error: size of new knot sequence does not match current one: %d != %d", knot_seq.size(), new_knots.size());
    return;
  }

  knot_seq = new_knots;

  for (int i = 0; i < knot_seq.size(); i++) {
    if (knot_seq[i] < 0.0)
      knot_seq[i] = 0.0;
    if (knot_seq[i] > 1.0)
      knot_seq[i] = 1.0;
  }

  // regenerate the ctrl pts for first derivative, given the new knot sequence
  dctrl_pts.erase(dctrl_pts.begin(),dctrl_pts.end());
  for (int i = 0; i < ctrl_pts.size()-1; i++) {
    dctrl_pts.push_back(degree / (knot_seq[i+degree+1] - knot_seq[i+1]) * (ctrl_pts[i+1] - ctrl_pts[i]));
  }
  // regenerate the first derivative's knot sequence
  dknot_seq.erase(dknot_seq.begin(),dknot_seq.end());
  dknot_seq = knot_seq;
  dknot_seq.pop_back(); // remove last element
  dknot_seq.erase(dknot_seq.begin()); // remove first element

  return;
}

void BSpline::CopyKnotSeq (vector<double> &dest_knots)
{
  dest_knots.erase(dest_knots.begin(),dest_knots.end());

  dest_knots = knot_seq;

  return;
}

V3d BSpline::deBoor (double u, int degree, const vector<V3d> &ctrl_pts, const vector<double> &knot_seq)
{
  double omega;
  int h, s, knot;

  if (knot_seq.empty())
    return V3d(0,0,0);

  // find u in the knot sequence, u in [u_k,u_{k+1})
  for (knot = 0; knot < knot_seq.size()-1; knot++) {
    if ((u >= knot_seq[knot]) && (u < knot_seq[knot+1]))
      break;
  }

  // clamp the knot (an open curve is defined on [u_p,u_{m-p}] = [u_{degree},u{num_ctrl_pts - degree}]
  if (knot < degree){
    knot = degree;
    u = knot_seq[knot];
  }
  if (knot >= knot_seq.size()-1-degree){
    knot = knot_seq.size()-1-degree;
    u = knot_seq[knot];
  }

  // if u != u_k, insert u "degree" times
  h = degree;
  s = 0;
  // else if u_k has multiplicity s, insert u "degree-s" times
  if (fabs(u - knot_seq[knot]) < 1e-5) {
    s = 1;
    for (int i = knot; i > 0; --i) {
      if (knot_seq[i] == knot_seq[i-1])
        ++s;
      else
        break;
    }
    if (s > degree)
      s = degree;
    h = degree - s;
  }

  // copy the affected control points
  vector<V3d> pts(degree+1-s);
  for (int i = s; i <= degree; i++) {
      pts[degree-i] = ctrl_pts[knot-i];
  }

  // insert "knot" h times
  for (int r = 1; r <= h; r++) {
    for (int i = knot - s; i >= knot - degree + r; --i) {
      omega = (u - knot_seq[i]) / (knot_seq[i+degree-r+1] - knot_seq[i]);
      pts[i-(knot-degree)] = (1.0-omega)*pts[i-1-(knot-degree)] + omega*pts[i-(knot-degree)];
    }
  }

  return pts[degree-s];
}

V3d BSpline::C (double u)
// return the point on the b-spline curve, C(u)
{
  return deBoor(u, degree, ctrl_pts, knot_seq);
}

V3d BSpline::dCdu (double u)
// return the point on the first derivative of the b-spline curve, dC(u)/du
{
  return deBoor(u, degree-1, dctrl_pts, dknot_seq);
}

/* -------------------------------------------------------------------------------------- */


class BSplineBasisFunc
{
  public:
    BSplineBasisFunc(int num_funcs, int new_degree, int open_or_closed);
    ~BSplineBasisFunc();

    bool SetKnotSeq (const vector<double> &new_knots);
    void CopyKnotSeq (vector<double> &dest_knots);
    void ResetKnotSeq (int open_or_closed);

    void evaluateAt (vector<double> &basisWeights, double u) const;
    void computeFuncs (double u);
    double valueAtFunc (int index) const;
    
    void generateSamples (int num_samples);
    void generateSamples (void);
    double sampleAt (int sampleIndex, int funcIndex) const;

  private:
    int degree;
    vector<double> basis_funcs;
    vector<double> knot_seq;
    vector<vector<double>> samples;
};

BSplineBasisFunc::BSplineBasisFunc(int num_funcs, int new_degree, int open_or_closed) {

  if (new_degree > 0)
    degree = new_degree;
  else
    degree = 1;

  if (num_funcs <= degree)
    num_funcs = degree+1;

  // generate the knot vector
  knot_seq.clear();//erase(knot_seq.begin(),knot_seq.end());
  int knot_size = degree + num_funcs + 1;
  for (int i = 0; i < knot_size; i++) {
    if (open_or_closed == 0) { // open curve (or uniform curve)
      knot_seq.push_back(double(i)/(double)(knot_size-1)); // open curve
    }
    else { // clamped (or open uniform curve)
      if (i < degree+1)
        knot_seq.push_back(0.0);
      else if (i >= knot_size-degree-1)
        knot_seq.push_back(1.0);
      else
        knot_seq.push_back((double) (i-degree) / (double) (knot_size-1-2*degree));
    }
  }

  // set basis funcs at u = 0
  basis_funcs.resize(num_funcs);
  computeFuncs(0.);

  // set samples
  samples = vector<vector<double>>(20*num_funcs+1, vector<double>(num_funcs));
}

BSplineBasisFunc::~BSplineBasisFunc() {}

bool BSplineBasisFunc::SetKnotSeq (const vector<double> &new_knots)
{
  if (knot_seq.size() != new_knots.size()) {
    Printf ("Error: size of new knot sequence does not match current one: %d != %d\n", knot_seq.size(), new_knots.size());
    return false;
  }

  knot_seq = new_knots;

  for (int i = 0; i < knot_seq.size(); i++) {
    if (knot_seq[i] < 0.0)
      knot_seq[i] = 0.0;
    if (knot_seq[i] > 1.0)
      knot_seq[i] = 1.0;
  }

  return true;
}

void BSplineBasisFunc::CopyKnotSeq (vector<double> &dest_knots)
{
  dest_knots.erase(dest_knots.begin(),dest_knots.end());

  dest_knots = knot_seq;

  return;
}

void BSplineBasisFunc::ResetKnotSeq (int open_or_closed)
{
  // generate the knot vector
  knot_seq.clear();
  int knot_size = degree + basis_funcs.size() + 1;
  for (int i = 0; i < knot_size; i++) {
    if (open_or_closed == 0) { // open curve (or uniform curve)
      knot_seq.push_back(double(i)/(double)(knot_size-1)); // open curve
    }
    else { // clamped (or open uniform curve)
      if (i < degree+1)
        knot_seq.push_back(0.0);
      else if (i >= knot_size-degree-1)
        knot_seq.push_back(1.0);
      else
        knot_seq.push_back((double) (i-degree) / (double) (knot_size-1-2*degree));
    }
  }
}


void BSplineBasisFunc::evaluateAt(vector<double> &basisWeights, double u) const {

  int knot;
  // find u in the knot sequence, u in [u_k,u_{k+1})
  for (knot = 0; knot < knot_seq.size()-1; knot++) {
    if ((u >= knot_seq[knot]) && (u < knot_seq[knot+1]))
      break;
  }

  // clamp the knot 
  if (knot < degree){
    knot = degree;
    u = knot_seq[knot];
  }
  if (knot >= knot_seq.size()-1-degree){
    u = knot_seq[knot];
    knot = knot_seq.size()-1-degree - 1; // NUM_PTS - 1
  }

  std::fill(basisWeights.begin(),basisWeights.end(),0);
  basisWeights[knot] = 1.0;

  for (int d = 1; d <= degree; d++) {

    if (knot_seq[knot+1] - knot_seq[knot-d+1] > 0.0)
      basisWeights[knot-d] = (knot_seq[knot+1]-u) / (knot_seq[knot+1] - knot_seq[knot-d+1])
                             * basisWeights[knot-d+1];
    else
      basisWeights[knot-d] = 0.0;
    
    for (int i = knot-d+1; i <= knot-1; i++) { // compute internal terms 
      double first = 0.0;
      if (knot_seq[i+d] - knot_seq[i] > 0.0)
        first = (u - knot_seq[i]) / (knot_seq[i+d] - knot_seq[i]);

      double second = 0.0;
      if (knot_seq[i+1+d] - knot_seq[i+1] > 0.0)
        second = (knot_seq[i+1+d]-u) / (knot_seq[i+1+d] - knot_seq[i+1]);

      basisWeights[i] = first * basisWeights[i] + second * basisWeights[i+1];

    }

    if (knot_seq[knot+d] - knot_seq[knot] > 0.0)
      basisWeights[knot] = (u - knot_seq[knot]) / (knot_seq[knot+d] - knot_seq[knot])
                      * basisWeights[knot];
    else
      basisWeights[knot] = 0.0;
  }

  return;
}

void BSplineBasisFunc::computeFuncs(double u) {

  evaluateAt(basis_funcs, u);

  return;
}

double BSplineBasisFunc::valueAtFunc(int index) const {
  if (index < 0 || index >= basis_funcs.size())
    return -1.; // -1 signifies an error
  return basis_funcs[index];
}

void BSplineBasisFunc::generateSamples(int num_samples) {

  // resize if needed
  if (samples.size() != num_samples) {
    samples = vector<vector<double>>(num_samples, vector<double>(basis_funcs.size()));
  }

  generateSamples();
  
  return;
}

void BSplineBasisFunc::generateSamples(void) {

  // compute samples
  for (int ui = 0; ui < samples.size(); ui++) {
    double u = ((double)ui)/((double)samples.size()-1);

    evaluateAt(samples[ui], u); // set weight at 'u' for each basis function
  }

  return;
}

double BSplineBasisFunc::sampleAt(int sampleIndex, int funcIndex) const {
  if (funcIndex < 0 || funcIndex >= basis_funcs.size())
    return -1.; // -1 signifies an error
  if (sampleIndex < 0 || sampleIndex >= samples.size())
    return -1.;
  return samples[sampleIndex][funcIndex];
}
