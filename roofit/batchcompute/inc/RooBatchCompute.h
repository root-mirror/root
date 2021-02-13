#ifndef ROOFIT_BATCHCOMPUTE_ROOBATCHCOMPUTE_H
#define ROOFIT_BATCHCOMPUTE_ROOBATCHCOMPUTE_H

#include "RooSpan.h"
#include "RooVDTHeaders.h"
#include "RunContext.h"
#include "BracketAdapter.h"
#include "DllImport.h" //for R__EXTERN, needed for windows

class RooAbsReal;
class RooListProxy;

namespace RooBatchCompute {
  /**
   * \brief The interface which should be implemented to provide optimised evaluateSpan() functionality for RooFit PDF classes.
   *
   * This interface contains the signatures of the compute functions of every PDF that has an optimised implementation available.
   * These are the functions that perform the actual computations in batches.
   * \see dispatch, RooBatchComputeClass, RF_ARCH
   */ 
  class RooBatchComputeInterface {
  public:
    virtual ~RooBatchComputeInterface() = default;
    virtual RooSpan<double> computeArgusBG(const RooAbsReal*, RunContext&, RooSpan<const double> m, RooSpan<const double> m0, RooSpan<const double> c, RooSpan<const double> p) = 0;
    virtual void computeBernstein(size_t batchSize, double * __restrict output, const double * __restrict const xData, double xmin, double xmax, std::vector<double> coef) = 0;
    virtual RooSpan<double> computeBifurGauss(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> sigmaL, RooSpan<const double> sigmaR) = 0;
    virtual RooSpan<double> computeBukin(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> Xp, RooSpan<const double> sigp, RooSpan<const double> xi, RooSpan<const double> rho1, RooSpan<const double> rho2) = 0;
    virtual RooSpan<double> computeBreitWigner(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> width) = 0;
    virtual RooSpan<double> computeCBShape(const RooAbsReal*, RunContext&, RooSpan<const double> m, RooSpan<const double> m0, RooSpan<const double> sigma, RooSpan<const double> alpha, RooSpan<const double> n) = 0;
    virtual void computeChebychev(size_t batchSize, double * __restrict output, const double * __restrict const xData, double xmin, double xmax, std::vector<double> coef) = 0;
    virtual RooSpan<double> computeChiSquare(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> ndof) = 0;
    virtual RooSpan<double> computeDstD0BG(const RooAbsReal*, RunContext&, RooSpan<const double> dm, RooSpan<const double> dm0, RooSpan<const double> C, RooSpan<const double> A, RooSpan<const double> B) = 0;
    virtual RooSpan<double> computeExponential(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> c) = 0;
    virtual RooSpan<double> computeGamma(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> gamma, RooSpan<const double> beta, RooSpan<const double> mu) = 0;
    virtual RooSpan<double> computeGaussian(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> sigma) = 0;
    virtual RooSpan<double> computeJohnson(const RooAbsReal*, RunContext&, RooSpan<const double> mass, RooSpan<const double> mu, RooSpan<const double> lambda, RooSpan<const double> gamma, RooSpan<const double> delta, double massThreshold) = 0;
    virtual RooSpan<double> computeLandau(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> sigma) = 0;
    virtual RooSpan<double> computeLognormal(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> m0, RooSpan<const double> k) = 0;
    virtual RooSpan<double> computeNovosibirsk(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> peak, RooSpan<const double> width, RooSpan<const double> tail) = 0;
    virtual RooSpan<double> computePoisson(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> mean, bool protectNegative, bool noRounding) = 0;
    virtual void computePolynomial(size_t batchSize, double * __restrict output, const double * __restrict const xData, int lowestOrder, std::vector<BracketAdapterWithMask> &coef) = 0;
    virtual RooSpan<double> computeVoigtian(const RooAbsReal*, RunContext&, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> width, RooSpan<const double> sigma) = 0;
  };

  R__EXTERN RooBatchComputeInterface* dispatch;
}

#endif
