//
// The "fitEllipseTGraphRMM" macro uses the "ROOT::Math::Minimizer"
// interface for fitting an ellipse to a set of data points from a TGraph
//
// To try this macro, in a ROOT prompt, do:
// .L fitEllipseTGraphRMM.cxx // or ".L fitEllipseTGraphRMM.cxx++"
// fitEllipseTGraphRMM(TestGraphRMM());
// for (Int_t i=0; i<10; i++) { fitEllipseTGraphRMM(); gSystem->Sleep(333); }
//
// Last update: Thu Jul 31 18:00:00 UTC 2014
//
// Changes:
// 2014.07.31 - (initial version)
//

#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <cmath>
#include <iostream>

//
// ellipse_fcn calculates the "normalized distance" from the ellipse's center
//
// x, y = point for which one wants to calculate the "normalized distance"
// x0 = ellipse's "x" center
// y0 = ellipse's "y" center
// a = ellipse's "semimajor" axis along "x" (> 0)
// b = ellipse's "semiminor" axis along "y" (> 0)
// theta = ellipse's axes rotation angle (-45 ... 135 degrees)
//
Double_t ellipse_fcn(Double_t x, Double_t y,
                     Double_t x0, Double_t y0,
                     Double_t a, Double_t b,
                     Double_t theta) // (in degrees)
{
  Double_t v = 9999.9;
  if ((a == 0.0) || (b == 0.0)) return v; // just a precaution
  // shift the center
  x -= x0;
  y -= y0;
  // un-rotate the axes
  theta *= TMath::Pi() / 180.0; // degrees -> radians
  v = x;
  x = x * std::cos(theta) + y * std::sin(theta);
  y = y * std::cos(theta) - v * std::sin(theta);
  // "scale" axes
  x /= a;
  y /= b;
  // calculate the "normalized distance"
  v = x * x + y * y;
  v = std::sqrt(v);
  return v;
}

//
// x[0] = "x"
// x[1] = "y"
// params[0] = ellipse's "x" center ("x0")
// params[1] = ellipse's "y" center ("y0")
// params[2] = ellipse's "semimajor" axis along "x" ("a" > 0)
// params[3] = ellipse's "semiminor" axis along "y" ("b" > 0)
// params[4] = ellipse's axes rotation angle ("theta" = -45 ... 135 degrees)
//
Double_t ellipse_fcn(const Double_t *x, const Double_t *params)
{
  return ellipse_fcn(x[0], x[1], // "x", "y"
                     params[0], params[1], // "x0", "y0"
                     params[2], params[3], // "a", "b"
                     params[4]); // "theta" (in degrees)
}

//
// the TGraph to be fitted (used by the ellipse_TGraph_chi2 function below)
//
TGraph *ellipse_TGraph = ((TGraph *)0);
//
// x[0] = ellipse's "x" center ("x0")
// x[1] = ellipse's "y" center ("y0")
// x[2] = ellipse's "semimajor" axis along "x" ("a" > 0)
// x[3] = ellipse's "semiminor" axis along "y" ("b" > 0)
// x[4] = ellipse's axes rotation angle ("theta" = -45 ... 135 degrees)
//
Double_t ellipse_TGraph_chi2(const Double_t *x)
{
  Double_t v = 0.0;
  if (!ellipse_TGraph) return v; // just a precaution
  for (Int_t i = 0; i < ellipse_TGraph->GetN(); i++) {
    Double_t r = ellipse_fcn((ellipse_TGraph->GetX())[i], // "x"
                             (ellipse_TGraph->GetY())[i], // "y"
                             x[0], x[1], // "x0", "y0"
                             x[2], x[3], // "a", "b"
                             x[4]); // "theta" (in degrees)
    r -= 1.0; // ellipse's "radius" in "normalized coordinates" is always 1
    v += r * r;
  }
  return v;
}

//
// http://root.cern/drupal/content/numerical-minimization#multidim_minim
// http://root.cern/root/html534/tutorials/fit/NumericalMinimization.C.html
//
ROOT::Math::Minimizer *ellipse_TGraph_minimize(TGraph *g)
{
  if (!g) return 0; // just a precaution
  if (g->GetN() < 6) return 0; // just a precaution

  // set the TGraph to be fitted (used by the ellipse_TGraph_chi2 function)
  ellipse_TGraph = g;

  // create minimizer giving a name and (optionally) a name of the algorithm
#if 0 /* 0 or 1 */
  ROOT::Math::Minimizer* m =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
#elif 0 /* 0 or 1 */
  ROOT::Math::Minimizer* m =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
#elif 1 /* 0 or 1 */
  ROOT::Math::Minimizer* m =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
#else /* 0 or 1 */
  ROOT::Math::Minimizer* m =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
#endif /* 0 or 1 */

  // set tolerance, etc. ...
  m->SetMaxFunctionCalls(1000000); // for Minuit and Minuit2
  m->SetMaxIterations(100000); // for GSL
  m->SetTolerance(0.001); // edm
#if 1 /* 0 or 1 */
  m->SetPrintLevel(1);
#endif /* 0 or 1 */

  // create function wrapper for minimizer (a IMultiGenFunction type)
  ROOT::Math::Functor f(&ellipse_TGraph_chi2, 5);

  m->SetFunction(f);

  m->Clear(); // just a precaution

  // estimate all initial values (note: good initial values
  // are CRUCIAL for the minimizing procedure to succeed)
  Double_t xmin, xmax, ymin, ymax;
  Double_t x0, y0, a, b, theta;
  ellipse_TGraph->ComputeRange(xmin, ymin, xmax, ymax);
  x0 = (xmax + xmin) / 2.0;
  y0 = (ymax + ymin) / 2.0;
  a = (ellipse_TGraph->GetX())[0] - x0;
  b = (ellipse_TGraph->GetY())[0] - y0;
  theta = ((std::abs(b) > 9999.9 * std::abs(a)) ? 9999.9 : (b / a));
  a = a * a + b * b;
  b = a;
  for (Int_t i = 1; i < ellipse_TGraph->GetN(); i++) {
    Double_t dx = (ellipse_TGraph->GetX())[i] - x0;
    Double_t dy = (ellipse_TGraph->GetY())[i] - y0;
    Double_t d = dx * dx + dy * dy;
    // try to keep "a" > "b"
    if (a < d) {
      a = d;
      theta = ((std::abs(dy) > 9999.9 * std::abs(dx)) ? 9999.9 : (dy / dx));
    }
    if (b > d) b = d;
  }
  a = std::sqrt(a); if (!(a > 0)) a = 0.001;
  b = std::sqrt(b); if (!(b > 0)) b = 0.001;
  theta = std::atan(theta) * 180.0 / TMath::Pi();
  if (theta < -45.0) theta += 180.0; // "theta" = -45 ... 135 degrees

  // set the variables to be minimized
  m->SetVariable(0, "x0", x0, (xmax - xmin) / 100.0);
  m->SetVariable(1, "y0", y0, (ymax - ymin) / 100.0);
  m->SetVariable(2, "a", a, a / 100.0);
  m->SetVariable(3, "b", b, b / 100.0);
  m->SetVariable(4, "theta", theta, 1);

#if 1 /* 0 or 1 */
  // set the variables' limits
  m->SetVariableLimits(0, xmin, xmax);
  m->SetVariableLimits(1, ymin, ymax);
  if (theta < 45.0) {
    if (a < ((xmax - xmin) / 2.0)) a = (xmax - xmin) / 2.0;
    if (b < ((ymax - ymin) / 2.0)) b = (ymax - ymin) / 2.0;
  } else {
    if (a < ((ymax - ymin) / 2.0)) a = (ymax - ymin) / 2.0;
    if (b < ((xmax - xmin) / 2.0)) b = (xmax - xmin) / 2.0;
  }
  m->SetVariableLimits(2, 0, a * 3.0);
  m->SetVariableLimits(3, 0, b * 3.0);
  // m->SetVariableLimits(4, theta - 30.0, theta + 30.0); // theta -+ 30
  m->SetVariableLimits(4, theta - 45.0, theta + 45.0); // theta -+ 45
#endif /* 0 or 1 */

  // do the minimization
  m->Minimize();

#if 0 /* 0 or 1 */
  const Double_t *xm = m->X();
  std::cout << "Minimum ( "
            << xm[0] << " , " << xm[1] << " , " // "x0", "y0"
            << xm[2] << " , " << xm[3] << " , " // "a", "b"
            << xm[4] << " ) = " // "theta" (in degrees)
            << m->MinValue() // it's equal to ellipse_TGraph_chi2(xm)
            << std::endl;
#endif /* 0 or 1 */

  return m;
}

//
// creates a test TGraph with an ellipse
//
TGraph *TestGraphRMM(Bool_t randomize = kFALSE) {
  Int_t i;

  // define the test ellipse
  Double_t x0 = 4; // ellipse's "x" center
  Double_t y0 = 3; // ellipse's "y" center
  Double_t a = 2; // ellipse's "semimajor" axis along "x" (> 0)
  Double_t b = 1; // ellipse's "semiminor" axis along "y" (> 0)
  Double_t theta = 100; // ellipse's axes rotation angle (-45 ... 135 degrees)

  // gRandom->SetSeed(0);
  if (randomize) {
    x0 = 10.0 - 20.0 * gRandom->Rndm();
    y0 = 10.0 - 20.0 * gRandom->Rndm();
    a = 0.5 + 4.5 * gRandom->Rndm();
    b = 0.5 + 4.5 * gRandom->Rndm();
    theta = 180.0 - 360.0 * gRandom->Rndm();
  }

  const Int_t n = 100; // number of points
  Double_t x[n], y[n];
  Double_t dt = TMath::TwoPi() / Double_t(n);
  Double_t tmp;
  theta *= TMath::PiOver2() / 90.0; // degrees -> radians
  for (i = 0; i < n; i++) {
    x[i] = a * (std::cos(dt * Double_t(i)) + 0.1 * gRandom->Rndm() - 0.05);
    y[i] = b * (std::sin(dt * Double_t(i)) + 0.1 * gRandom->Rndm() - 0.05);
    // rotate the axes
    tmp = x[i];
    x[i] = x[i] * std::cos(theta) - y[i] * std::sin(theta);
    y[i] = y[i] * std::cos(theta) + tmp * std::sin(theta);
    // shift the center
    x[i] += x0;
    y[i] += y0;
  }

  // create the test TGraph
  TGraph *g = ((TGraph *)(gROOT->FindObject("g")));
  if (g) delete g;
  g = new TGraph(n, x, y);
  g->SetNameTitle("g", "test ellipse");

  return g;
}

//
// "ROOT Script" entry point (the same name as the "filename's base")
//
void fitEllipseTGraphRMM(TGraph *g = ((TGraph *)0))
{
  if (!g) g = TestGraphRMM(kTRUE); // create a "random" ellipse

#if 0 /* 0 or 1 */
  // create the "ellipse" TF2 (just for fun)
  TF2 *ellipse = ((TF2 *)(gROOT->GetListOfFunctions()->FindObject("ellipse")));
  if (ellipse) delete ellipse;
  ellipse = new TF2("ellipse", ellipse_fcn, -1, 1, -1, 1, 5);
  ellipse->SetMaximum(2.0); // just for nice graphics
  ellipse->SetParNames("x0", "y0", "a", "b", "theta");
  ellipse->SetParameters(0.4, 0.3, 0.2, 0.1, 10);
#endif /* 0 or 1 */

  // fit the TGraph
  ROOT::Math::Minimizer *m = ellipse_TGraph_minimize(g);

#if 1 /* 0 or 1 */
  // draw everything
  TCanvas *c = ((TCanvas *)(gROOT->GetListOfCanvases()->FindObject("c")));
  if (c) { c->Clear(); } else { c = new TCanvas("c", "c"); }
  c->SetGrid(1, 1);
  g->Draw("A*");
  if ( m && (!(m->Status())) ) {
    const Double_t *xm = m->X();
    TEllipse *e = new TEllipse(xm[0], xm[1], // "x0", "y0"
                               xm[2], xm[3], // "a", "b"
                               0, 360,
                               xm[4]); // "theta" (in degrees)
    e->SetFillStyle(0); // hollow
    e->Draw();
  }
  c->Modified(); c->Update(); // make sure it's really drawn
#endif /* 0 or 1 */

  delete m; // "cleanup"

  return;
}

// end of file fitEllipseTGraphRMM.cxx by Silesius Anonymus
