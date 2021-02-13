// RooBatchCompute library created September 2020 by Emmanouil Michalainas
#include "RooBatchCompute.h"
#include "RooMath.h"

#include <complex>

namespace RooBatchCompute {
  /**
   * \brief Contains the part of the code of the RooBatchCompute Library that needs to be compiled for every different cpu architecture.
   *
   * RF_ARCH is a macro that is defined by cmake. The macro gets a different name for each copy of the library, namely
   * GENERIC, SSE4, AVX, AVX2, AVX512, CUDA. This ensures that name clashes are avoided.
   * \see RooBatchComputeInterface, RooBatchComputeClass, dispatch
   */
  namespace RF_ARCH {

    struct ArgusBGComputer {
      template<class Tm, class Tm0, class Tc, class Tp>
      void run(size_t batchSize, double * __restrict output, Tm M, Tm0 M0, Tc C, Tp P ) const
      {
        for (size_t i=0; i<batchSize; i++) {
          const double t = M[i]/M0[i];
          const double u = 1 - t*t;
          output[i] = C[i]*u + P[i]*fast_log(u);
        }
        for (size_t i=0; i<batchSize; i++) {
          if (M[i] >= M0[i]) output[i] = 0.0;
          else output[i] = M[i]*fast_exp(output[i]);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void startComputationBernstein(size_t batchSize, double * __restrict output, const double * __restrict const xData, double xmin, double xmax, std::vector<double> coef)
    {
      constexpr size_t block = 128;
      const int nCoef = coef.size();
      const int degree = nCoef-1;
      double X[block], _1_X[block], powX[block], pow_1_X[block];
      double *Binomial = new double[nCoef+5];
      //Binomial stores values c(degree,i) for i in [0..degree]

      Binomial[0] = 1.0;
      for (int i=1; i<=degree; i++) {
        Binomial[i] = Binomial[i-1]*(degree-i+1)/i;
      }

      for (size_t i=0; i<batchSize; i+=block) {
        const size_t stop = (i+block > batchSize) ? batchSize-i : block;

        //initialization
        for (size_t j=0; j<stop; j++) {
          powX[j] = pow_1_X[j] = 1.0;
          X[j] = (xData[i+j]-xmin) / (xmax-xmin);
          _1_X[j] = 1-X[j];
          output[i+j] = 0.0;
        }

        //raising 1-x to the power of degree
        for (int k=2; k<=degree; k+=2)
          for (size_t j=0; j<stop; j++)
            pow_1_X[j] *= _1_X[j]*_1_X[j];

        if (degree%2 == 1)
          for (size_t j=0; j<stop; j++)
            pow_1_X[j] *= _1_X[j];

        //inverting 1-x ---> 1/(1-x)
        for (size_t j=0; j<stop; j++)
          _1_X[j] = 1/_1_X[j];

        for (int k=0; k<nCoef; k++) {
          for (size_t j=0; j<stop; j++) {
            output[i+j] += coef[k]*Binomial[k]*powX[j]*pow_1_X[j];

            //calculating next power for x and 1-x
            powX[j] *= X[j];
            pow_1_X[j] *= _1_X[j];
          }
        }
      }
      delete[] Binomial;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BifurGaussComputer {
      template<class Tx, class Tm, class Tsl, class Tsr>
      void run(size_t batchSize, double * __restrict output, Tx X, Tm M, Tsl SL, Tsr SR) const
      {
        for (size_t i=0; i<batchSize; i++) {
          const double arg = X[i]-M[i];
          output[i] = arg / ((arg < 0.0)*SL[i] + (arg >= 0.0)*SR[i]);
        }

        for (size_t i=0; i<batchSize; i++) {
          if (X[i]-M[i]>1e-30 || X[i]-M[i]<-1e-30) {
            output[i] = fast_exp(-0.5*output[i]*output[i]);
          }
          else {
            output[i] = 1.0;
          }
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BukinComputer {
      template<class Tx, class TXp, class TSigp, class Txi, class Trho1, class Trho2>
      void run(size_t batchSize, double * __restrict output, Tx X, TXp XP, TSigp SP, Txi XI, Trho1 R1, Trho2 R2) const
      {
        const double r3 = log(2.0);
        const double r6 = exp(-6.0);
        const double r7 = 2*sqrt(2*log(2.0));

        for (size_t i=0; i<batchSize; i++) {
          const double r1 = XI[i]*fast_isqrt(XI[i]*XI[i]+1);
          const double r4 = 1/fast_isqrt(XI[i]*XI[i]+1);
          const double hp = 1 / (SP[i]*r7);
          const double x1 = XP[i] + 0.5*SP[i]*r7*(r1-1);
          const double x2 = XP[i] + 0.5*SP[i]*r7*(r1+1);

          double r5 = 1.0;
          if (XI[i]>r6 || XI[i]<-r6) r5 = XI[i]/fast_log(r4+XI[i]);

          double factor=1, y=X[i]-x1, Yp=XP[i]-x1, yi=r4-XI[i], rho=R1[i];
          if (X[i]>=x2) {
            factor = -1;
            y = X[i]-x2;
            Yp = XP[i]-x2;
            yi = r4+XI[i];
            rho = R2[i];
          }

          output[i] = rho*y*y/Yp/Yp -r3 + factor*4*r3*y*hp*r5*r4/yi/yi;
          if (X[i]>=x1 && X[i]<x2) {
            output[i] = fast_log(1 + 4*XI[i]*r4*(X[i]-XP[i])*hp) / fast_log(1 +2*XI[i]*( XI[i]-r4 ));
            output[i] *= -output[i]*r3;
          }
          if (X[i]>=x1 && X[i]<x2 && XI[i]<r6 && XI[i]>-r6) {
            output[i] = -4*r3*(X[i]-XP[i])*(X[i]-XP[i])*hp*hp;
          }
        }
        for (size_t i=0; i<batchSize; i++) {
          output[i] = fast_exp(output[i]);
        }
      }
    };

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BreitWignerComputer {
      template<class Tx, class Tmean, class Twidth>
      void run(size_t batchSize, double * __restrict output, Tx X, Tmean M, Twidth W) const
      {
        for (size_t i=0; i<batchSize; i++) {
          const double arg = X[i]-M[i];
          output[i] = 1 / (arg*arg + 0.25*W[i]*W[i]);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct CBShapeComputer {
      template<class Tm, class Tm0, class Tsigma, class Talpha, class Tn>
      void run(	size_t batchSize, double * __restrict output, Tm M, Tm0 M0, Tsigma S, Talpha A, Tn N) const
      {
        for (size_t i=0; i<batchSize; i++) {
          const double t = (M[i]-M0[i]) / S[i];
          if ( (A[i]>0 && t>=-A[i]) || (A[i]<0 && -t>=A[i]) ) {
            output[i] = -0.5*t*t;
          } else {
            output[i] = N[i] / (N[i] -A[i]*A[i] -A[i]*t);
            output[i] = fast_log(output[i]);
            output[i] *= N[i];
            output[i] -= 0.5*A[i]*A[i];
          }
        }

        for (size_t i=0; i<batchSize; i++) {
          output[i] = fast_exp(output[i]);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void startComputationChebychev(size_t batchSize, double * __restrict output, const double * __restrict const xData, double xmin, double xmax, std::vector<double> coef)
  {
    constexpr size_t block = 128;
    const size_t nCoef = coef.size();
    double prev[block][2], X[block];

    for (size_t i=0; i<batchSize; i+=block) {
      size_t stop = (i+block >= batchSize) ? batchSize-i : block;

      // set a0-->prev[j][0] and a1-->prev[j][1]
      // and x tranfsformed to range[-1..1]-->X[j]
      for (size_t j=0; j<stop; j++) {
        prev[j][0] = output[i+j] = 1.0;
        prev[j][1] = X[j] = (xData[i+j] -0.5*(xmax + xmin)) / (0.5*(xmax - xmin));
      }

      for (size_t k=0; k<nCoef; k++) {
        for (size_t j=0; j<stop; j++) {
          output[i+j] += prev[j][1]*coef[k];

          //compute next order
          const double next = 2*X[j]*prev[j][1] -prev[j][0];
          prev[j][0] = prev[j][1];
          prev[j][1] = next;
        }
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct ChiSquareComputer {
    template<class T_x, class T_ndof>
    void run(size_t batchSize, double * __restrict output, T_x X, T_ndof N) const
    {
      if ( N.isBatch() ) {
        for (size_t i=0; i<batchSize; i++) {
          if (X[i] > 0) {
            output[i] = 1/std::tgamma(N[i]/2.0);
          }
        }
      }
      else {
        // N is just a scalar so bracket adapter ignores index.
        const double gamma = 1/std::tgamma(N[2019]/2.0);
        for (size_t i=0; i<batchSize; i++) {
          output[i] = gamma;
        }
      }

      constexpr double ln2 = 0.693147180559945309417232121458;
      const double lnx0 = std::log(X[0]);
      for (size_t i=0; i<batchSize; i++) {
        double lnx;
        if ( X.isBatch() ) lnx = fast_log(X[i]);
        else lnx = lnx0;

        double arg = (N[i]-2)*lnx -X[i] -N[i]*ln2;
        output[i] *= fast_exp(0.5*arg);
      }
    }
  };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct DstD0BGComputer {
      template<class Tdm, class Tdm0, class TC, class TA, class TB>
      void run(size_t batchSize, double * __restrict output, Tdm DM, Tdm0 DM0, TC C, TA A, TB B) const
      {
        for (size_t i=0; i<batchSize; i++) {
          const double ratio = DM[i] / DM0[i];
          const double arg1 = (DM0[i]-DM[i]) / C[i];
          const double arg2 = A[i]*fast_log(ratio);
          output[i] = (1 -fast_exp(arg1)) * fast_exp(arg2) +B[i]*(ratio-1);
        }

        for (size_t i=0; i<batchSize; i++) {
          if (output[i]<0) output[i] = 0;
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct ExponentialComputer {
      template<class Tx, class Tc>
      void run(size_t n, double* __restrict output, Tx x, Tc c) const
      {
        for (size_t i = 0; i < n; ++i) {
          output[i] = fast_exp(x[i]*c[i]);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct GammaComputer {
      template<class Tx, class Tgamma, class Tbeta, class Tmu>
      void run (size_t batchSize, double * __restrict output, Tx X, Tgamma G, Tbeta B, Tmu M) const
      {
        constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
        for (size_t i=0; i<batchSize; i++) {
          if (X[i]<M[i] || G[i] <= 0.0 || B[i] <= 0.0) {
            output[i] = NaN;
          }
          if (X[i] == M[i]) {
            output[i] = ((G[i]==1.0) ? 1. : 0.)/B[i];
          }
          else {
            output[i] = 0.0;
          }
        }

        if (G.isBatch()) {
          for (size_t i=0; i<batchSize; i++) {
            if (output[i] == 0.0) {
            output[i] = -std::lgamma(G[i]);
            }
          }
        }
        else {
        double gamma = -std::lgamma(G[2019]);
        for (size_t i=0; i<batchSize; i++) {
          if (output[i] == 0.0) {
            output[i] = gamma;
            }
          }
        }

        for (size_t i=0; i<batchSize; i++) {
          if (X[i] != M[i]) {
            const double invBeta = 1/B[i];
            double arg = (X[i]-M[i])*invBeta;
            output[i] -= arg;
            arg = fast_log(arg);
            output[i] += arg*(G[i]-1);
            output[i] = fast_exp(output[i]);
            output[i] *= invBeta;
          }
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///Actual computations for the batch evaluation of the Gaussian.
  ///May vectorise over x, mean, sigma, depending on the types of the inputs.
  ///\note The output and input spans are assumed to be non-overlapping. If they
  ///overlap, results will likely be garbage.
    struct GaussianComputer {
      template<class Tx, class TMean, class TSig>
      void run(size_t n, double* __restrict output, Tx x, TMean mean, TSig sigma) const
      {
        for (std::size_t i=0; i<n; ++i) {
          const double arg = x[i]-mean[i];
          const double halfBySigmaSq = -0.5 / (sigma[i]*sigma[i]);
          output[i] = fast_exp(arg*arg*halfBySigmaSq);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct JohnsonComputer {
      ///Actual computations for the batch evaluation of the Johnson.
      ///May vectorise over observables depending on types of inputs.
      ///\note The output and input spans are assumed to be non-overlapping. If they
      ///overlap, results will likely be garbage.
      template<class TMass, class TMu, class TLambda, class TGamma, class TDelta>
      void run(size_t n, double* __restrict output, TMass mass, TMu mu, TLambda lambda, TGamma gamma, TDelta delta) const
      {
        const double sqrt_twoPi = sqrt(TMath::TwoPi());

        for (size_t i=0; i<n; ++i) {
          const double arg = (mass[i]-mu[i]) / lambda[i];
      #ifdef R__HAS_VDT
          const double asinh_arg = fast_log(arg + 1/fast_isqrt(arg*arg+1));
      #else
          const double asinh_arg = asinh(arg);
      #endif
          const double expo = gamma[i] + delta[i]*asinh_arg;
          const double result = delta[i]*fast_exp(-0.5*expo*expo)*fast_isqrt(1. +arg*arg) / (sqrt_twoPi*lambda[i]);

          const double passThrough = mass[i] >= massThreshold;
          output[i] = result*passThrough;
        }
      }
      const double massThreshold;
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct LandauComputer {
      /* Actual computation of Landau(x,mean,sigma) in a vectorization-friendly way
       * Code copied from function landau_pdf (math/mathcore/src/PdfFuncMathCore.cxx)
       * and rewritten to take advantage for the most popular case
       * which is -1 < (x-mean)/sigma < 1. The rest cases are handled in scalar way
       */
      template<class Tx, class Tmean, class Tsigma>
      void run(size_t batchSize, double* __restrict output, Tx X, Tmean M, Tsigma S) const
      {
        const double p1[5] = {0.4259894875,-0.1249762550, 0.03984243700, -0.006298287635,   0.001511162253};
        const double q1[5] = {1.0         ,-0.3388260629, 0.09594393323, -0.01608042283,    0.003778942063};

        const double p2[5] = {0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411,   0.0001283617211};
        const double q2[5] = {1.0         , 0.7428795082, 0.3153932961,   0.06694219548,    0.008790609714};

        const double p3[5] = {0.1788544503, 0.09359161662,0.006325387654, 0.00006611667319,-0.000002031049101};
        const double q3[5] = {1.0         , 0.6097809921, 0.2560616665,   0.04746722384,    0.006957301675};

        const double p4[5] = {0.9874054407, 118.6723273,  849.2794360,   -743.7792444,      427.0262186};
        const double q4[5] = {1.0         , 106.8615961,  337.6496214,    2016.712389,      1597.063511};

        const double p5[5] = {1.003675074,  167.5702434,  4789.711289,    21217.86767,     -22324.94910};
        const double q5[5] = {1.0         , 156.9424537,  3745.310488,    9834.698876,      66924.28357};

        const double p6[5] = {1.000827619,  664.9143136,  62972.92665,    475554.6998,     -5743609.109};
        const double q6[5] = {1.0         , 651.4101098,  56974.73333,    165917.4725,     -2815759.939};

        const double a1[3] = {0.04166666667,-0.01996527778, 0.02709538966};
        const double a2[2] = {-1.845568670,-4.284640743};

        const double NaN = std::nan("");
        const size_t block=256;
        double v[block];

        for (size_t i=0; i<batchSize; i+=block) { //CHECK_VECTORISE
          const size_t stop = (i+block < batchSize) ? block : batchSize-i ;

          for (size_t j=0; j<stop; j++) { //CHECK_VECTORISE
            v[j] = (X[i+j]-M[i+j]) / S[i+j];
            output[i+j] = (p2[0]+(p2[1]+(p2[2]+(p2[3]+p2[4]*v[j])*v[j])*v[j])*v[j]) /
                          (q2[0]+(q2[1]+(q2[2]+(q2[3]+q2[4]*v[j])*v[j])*v[j])*v[j]);
          }

          for (size_t j=0; j<stop; j++) { //CHECK_VECTORISE
            const bool mask = S[i+j] > 0;
            /*  comparison with NaN will give result false, so the next
             *  loop won't affect output, for cases where sigma <=0
             */
            if (!mask) v[j] = NaN;
            output[i+j] *= mask;
          }

          double u, ue, us;
          for (size_t j=0; j<stop; j++) { //CHECK_VECTORISE
            // if branch written in way to quickly process the most popular case -1 <= v[j] < 1
            if (v[j] >= 1) {
              if (v[j] < 5) {
                output[i+j] = (p3[0]+(p3[1]+(p3[2]+(p3[3]+p3[4]*v[j])*v[j])*v[j])*v[j]) /
                         (q3[0]+(q3[1]+(q3[2]+(q3[3]+q3[4]*v[j])*v[j])*v[j])*v[j]);
              } else if (v[j] < 12) {
                  u   = 1/v[j];
                  output[i+j] = u*u*(p4[0]+(p4[1]+(p4[2]+(p4[3]+p4[4]*u)*u)*u)*u) /
                          (q4[0]+(q4[1]+(q4[2]+(q4[3]+q4[4]*u)*u)*u)*u);
              } else if (v[j] < 50) {
                  u   = 1/v[j];
                  output[i+j] = u*u*(p5[0]+(p5[1]+(p5[2]+(p5[3]+p5[4]*u)*u)*u)*u) /
                           (q5[0]+(q5[1]+(q5[2]+(q5[3]+q5[4]*u)*u)*u)*u);
              } else if (v[j] < 300) {
                  u   = 1/v[j];
                  output[i+j] = u*u*(p6[0]+(p6[1]+(p6[2]+(p6[3]+p6[4]*u)*u)*u)*u) /
                           (q6[0]+(q6[1]+(q6[2]+(q6[3]+q6[4]*u)*u)*u)*u);
              } else {
                  u   = 1 / (v[j] -v[j]*std::log(v[j])/(v[j]+1) );
                  output[i+j] = u*u*(1 +(a2[0] +a2[1]*u)*u );
              }
            } else if (v[j] < -1) {
                if (v[j] >= -5.5) {
                  u   = std::exp(-v[j]-1);
                  output[i+j] = std::exp(-u)*std::sqrt(u)*
                    (p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*v[j])*v[j])*v[j])*v[j])/
                    (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*v[j])*v[j])*v[j])*v[j]);
                } else  {
                    u   = std::exp(v[j]+1.0);
                    if (u < 1e-10) output[i+j] = 0.0;
                    else {
                      ue  = std::exp(-1/u);
                      us  = std::sqrt(u);
                      output[i+j] = 0.3989422803*(ue/us)*(1+(a1[0]+(a1[1]+a1[2]*u)*u)*u);
                    }
                }
              }
          }
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct LognormalComputer {
      template<class Tx, class Tm0, class Tk>
      void run(size_t batchSize, double* __restrict output, Tx X, Tm0 M0, Tk K) const
      {
        const double rootOf2pi = 2.506628274631000502415765284811;
        for (size_t i=0; i<batchSize; i++) {
          double lnxOverM0 = fast_log(X[i]/M0[i]);
          double lnk = fast_log(K[i]);
          if (lnk<0) lnk = -lnk;
          double arg = lnxOverM0/lnk;
          arg *= -0.5*arg;
          output[i] = fast_exp(arg) / (X[i]*lnk*rootOf2pi);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct NovosibirskComputer {
      /* TMath::ASinH(x) needs to be replaced with ln( x + sqrt(x^2+1))
       * argasinh -> the argument of TMath::ASinH()
       * argln -> the argument of the logarithm that replaces AsinH
       * asinh -> the value that the function evaluates to
       *
       * ln is the logarithm that was solely present in the initial
       * formula, that is before the asinh replacement
       */
      template<class Tx, class Twidth, class Tpeak, class Ttail>
      void run(size_t batchSize, double * __restrict output, Tx X, Tpeak P, Twidth W, Ttail T) const
      {
        constexpr double xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )
        for (size_t i=0; i<batchSize; i++) {
          double argasinh = 0.5*xi*T[i];
          double argln = argasinh + 1/fast_isqrt(argasinh*argasinh +1);
          double asinh = fast_log(argln);

          double argln2 = 1 -(X[i]-P[i])*T[i]/W[i];
          double ln    = fast_log(argln2);
          output[i] = ln/asinh;
          output[i] *= -0.125*xi*xi*output[i];
          output[i] -= 2.0/xi/xi*asinh*asinh;
        }

        //faster if you exponentiate in a seperate loop (dark magic!)
        for (size_t i=0; i<batchSize; i++) {
          output[i] = fast_exp(output[i]);
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct PoissonComputer {
      template<class Tx, class TMean>
      void run(const size_t n, double* __restrict output, Tx x, TMean mean) const
      {
        for (size_t i = 0; i < n; ++i) { //CHECK_VECTORISE
          const double x_i = noRounding ? x[i] : floor(x[i]);
          // The std::lgamma yields different values than in the scalar implementation.
          // Need to check which one is more accurate.
          // output[i] = std::lgamma(x_i + 1.);
          output[i] = TMath::LnGamma(x_i + 1.);
        }

        for (size_t i = 0; i < n; ++i) {
          const double x_i = noRounding ? x[i] : floor(x[i]);
          const double logMean = fast_log(mean[i]);
          const double logPoisson = x_i * logMean - mean[i] - output[i];
          output[i] = fast_exp(logPoisson);

          // Cosmetics
          if (x_i < 0.)
            output[i] = 0.;
          else if (x_i == 0.) {
            output[i] = 1./fast_exp(mean[i]);
          }
          if (protectNegative && mean[i] < 0.)
            output[i] = 1.E-3;
        }
      }
      const bool protectNegative;
      const bool noRounding;
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void startComputationPolynomial(size_t batchSize, double* __restrict output, const double* __restrict const X, int lowestOrder, std::vector<BracketAdapterWithMask>& coefList )
    {
      const int nCoef = coefList.size();
      if (nCoef==0 && lowestOrder==0) {
        for (size_t i=0; i<batchSize; i++) {
          output[i] = 0.0;
        }
      }
      else if (nCoef==0 && lowestOrder>0) {
        for (size_t i=0; i<batchSize; i++) {
          output[i] = 1.0;
        }
      } else {
        for (size_t i=0; i<batchSize; i++) {
          output[i] = coefList[nCoef-1][i];
        }
      }
      if (nCoef == 0) return;

      /* Indexes are in range 0..nCoef-1 but coefList[nCoef-1]
       * has already been processed. In order to traverse the list,
       * with step of 2 we have to start at index nCoef-3 and use
       * coefList[k+1] and coefList[k]
       */
      for (int k=nCoef-3; k>=0; k-=2) {
        for (size_t i=0; i<batchSize; i++) {
          double coef1 = coefList[k+1][i];
          double coef2 = coefList[ k ][i];
          output[i] = X[i]*(output[i]*X[i] + coef1) + coef2;
        }
      }
      // If nCoef is odd, then the coefList[0] didn't get processed
      if (nCoef%2 == 0) {
        for (size_t i=0; i<batchSize; i++) {
          output[i] = output[i]*X[i] + coefList[0][i];
        }
      }
      //Increase the order of the polynomial, first by myltiplying with X[i]^2
      if (lowestOrder == 0) return;
      for (int k=2; k<=lowestOrder; k+=2) {
        for (size_t i=0; i<batchSize; i++) {
          output[i] *= X[i]*X[i];
        }
      }
      const bool isOdd = lowestOrder%2;
      for (size_t i=0; i<batchSize; i++) {
        if (isOdd) output[i] *= X[i];
        output[i] += 1.0;
      }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct VoigtianComputer {
      template<class Tx, class Tmean, class Twidth, class Tsigma>
      void run(size_t batchSize, double * __restrict output, Tx X, Tmean M, Twidth W, Tsigma S) const
      {
        const double invSqrt2 = 0.707106781186547524400844362105;
        for (size_t i=0; i<batchSize; i++) {
          const double arg = (X[i]-M[i])*(X[i]-M[i]);
          if (S[i]==0.0 && W[i]==0.0) {
            output[i] = 1.0;
          } else if (S[i]==0.0) {
            output[i] = 1/(arg+0.25*W[i]*W[i]);
          } else if (W[i]==0.0) {
            output[i] = fast_exp(-0.5*arg/(S[i]*S[i]));
          } else {
            output[i] = invSqrt2/S[i];
          }
        }

        for (size_t i=0; i<batchSize; i++) {
          if (S[i]!=0.0 && W[i]!=0.0) {
            if (output[i] < 0) output[i] = -output[i];
            const double factor = W[i]>0.0 ? 0.5 : -0.5;
            std::complex<Double_t> z( output[i]*(X[i]-M[i]) , factor*output[i]*W[i] );
            output[i] *= RooMath::faddeeva(z).real();
          }
        }
      }
    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class RooBatchComputeClass : public RooBatchComputeInterface {
      private:

        struct AnalysisInfo {
          size_t batchSize=SIZE_MAX;
          bool canDoHighPerf=true;
        };
        /// Small helping function that determines the sizes of the batches and executes either
        /// * A high-performance and optimized compute function instance, if the observable is a span and all parameters are const scalars
        /// * A less optimized one-fits-all compute function instance that covers every other (rare) scenario
        AnalysisInfo analyseInputSpans(std::vector<RooSpan<const double>> parameters)
        {
          AnalysisInfo ret;
          if (parameters[0].size()<=1) ret.canDoHighPerf=false;
          else ret.batchSize = std::min(ret.batchSize, parameters[0].size());
          for (size_t i=1; i<parameters.size(); i++)
            if (parameters[i].size()>1)
            {
              ret.canDoHighPerf=false;
              ret.batchSize = std::min(ret.batchSize, parameters[i].size());
            }
          return ret;
        }

        /// Templated function that works for every PDF: does the necessary preprocessing and launches
        /// the correct overload of the actual computing function.
        template <class Computer_t, typename Arg_t, typename... Args_t>
        RooSpan<double> startComputation(const RooAbsReal* caller, RunContext& evalData, Computer_t computer, Arg_t first, Args_t... rest)
        {
          AnalysisInfo info = analyseInputSpans({first, rest...});
          RooSpan<double> output = evalData.makeBatch(caller, info.batchSize);

          if (info.canDoHighPerf) computer.run(info.batchSize, output.data(), first, BracketAdapter<double>(rest[0])...);
          else                    computer.run(info.batchSize, output.data(), BracketAdapterWithMask(first), BracketAdapterWithMask(rest)...);

          return output;
        }

      public:
        RooBatchComputeClass() {
          // Set the dispatch pointer to this instance of the library upon loading
          RooBatchCompute::dispatch = this;
        }
        RooSpan<double> computeArgusBG(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> m, RooSpan<const double> m0, RooSpan<const double> c, RooSpan<const double> p)  override {
          return startComputation(caller, evalData, ArgusBGComputer{}, m, m0, c, p);
        }
        void computeBernstein(size_t batchSize, double * __restrict output, const double * __restrict const xData, double xmin, double xmax, std::vector<double> coef)  override {
          startComputationBernstein(batchSize, output, xData, xmin, xmax, coef);
        }
        RooSpan<double> computeBifurGauss(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> sigmaL, RooSpan<const double> sigmaR)  override {
          return startComputation(caller, evalData, BifurGaussComputer{}, x, mean, sigmaL, sigmaR);
        }
        RooSpan<double> computeBukin(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> Xp, RooSpan<const double> sigp, RooSpan<const double> xi, RooSpan<const double> rho1, RooSpan<const double> rho2)  override {
          return startComputation(caller, evalData, BukinComputer{}, x, Xp, sigp, xi, rho1, rho2);
        }
        RooSpan<double> computeBreitWigner(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> width)  override {
          return startComputation(caller, evalData, BreitWignerComputer{}, x, mean, width);
        }
        RooSpan<double> computeCBShape(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> m, RooSpan<const double> m0, RooSpan<const double> sigma, RooSpan<const double> alpha, RooSpan<const double> n)  override {
          return startComputation(caller, evalData, CBShapeComputer{}, m, m0, sigma, alpha, n);
        }
        void computeChebychev(size_t batchSize, double * __restrict output, const double * __restrict const xData, double xmin, double xmax, std::vector<double> coef)  override {
          startComputationChebychev(batchSize, output, xData, xmin, xmax, coef);
        }
        RooSpan<double> computeChiSquare(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> ndof)  override {
          return startComputation(caller, evalData, ChiSquareComputer{}, x, ndof);
        }
        RooSpan<double> computeDstD0BG(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> dm, RooSpan<const double> dm0, RooSpan<const double> C, RooSpan<const double> A, RooSpan<const double> B)  override {
          return startComputation(caller, evalData, DstD0BGComputer{}, dm, dm0, C, A, B);
        }
        RooSpan<double> computeExponential(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> c)  override {
          return startComputation(caller, evalData, ExponentialComputer{}, x, c);
        }
        RooSpan<double> computeGamma(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> gamma, RooSpan<const double> beta, RooSpan<const double> mu)  override {
          return startComputation(caller, evalData, GammaComputer{}, x, gamma, beta, mu);
        }
        RooSpan<double> computeGaussian(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> sigma)  override {
          return startComputation(caller, evalData, GaussianComputer{}, x, mean, sigma);
        }
        RooSpan<double> computeJohnson(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> mass, RooSpan<const double> mu, RooSpan<const double> lambda, RooSpan<const double> gamma, RooSpan<const double> delta, double massThreshold)  override {
          return startComputation(caller, evalData, JohnsonComputer{massThreshold}, mass, mu, lambda, gamma, delta);
        }
        RooSpan<double> computeLandau(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> sigma)  override {
          return startComputation(caller, evalData, LandauComputer{}, x, mean, sigma);
        }
        RooSpan<double> computeLognormal(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> m0, RooSpan<const double> k)  override {
          return startComputation(caller, evalData, LognormalComputer{}, x, m0, k);
        }
        RooSpan<double> computeNovosibirsk(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> peak, RooSpan<const double> width, RooSpan<const double> tail)  override {
          return startComputation(caller, evalData, NovosibirskComputer{}, x, peak, width, tail);
        }
        RooSpan<double> computePoisson(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> mean, bool protectNegative, bool noRounding)  override {
          return startComputation(caller, evalData, PoissonComputer{protectNegative, noRounding}, x, mean);
        }
        void computePolynomial(size_t batchSize, double* __restrict output, const double* __restrict const X, int lowestOrder, std::vector<BracketAdapterWithMask>& coefList)  override {
          startComputationPolynomial(batchSize, output, X, lowestOrder, coefList);
        }
        RooSpan<double> computeVoigtian(const RooAbsReal* caller, RunContext& evalData, RooSpan<const double> x, RooSpan<const double> mean, RooSpan<const double> width, RooSpan<const double> sigma)  override {
          return startComputation(caller, evalData, VoigtianComputer{}, x, mean, width, sigma);
        }
    }; // End class RooBatchComputeClass

    /// Static object to trigger the constructor which overwrites the dispatch pointer.
    static RooBatchComputeClass computeObj;

  } //End namespace RF_ARCH
} //End namespace RooBatchCompute
