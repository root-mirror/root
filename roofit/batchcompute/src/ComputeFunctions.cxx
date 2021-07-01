// Auto-vectorizable compute functions written by Emmanouil Michalainas, summer 2019

#include "RooBatchCompute.h"
#include "Batches.h"
#include "RooVDTHeaders.h"

#include <complex>

#ifdef __CUDACC__
  #define BEGIN blockDim.x*blockIdx.x + threadIdx.x
  #define STEP  blockDim.x*gridDim.x
  // for CUDA, we include RooMath.cxx and ComputeFunctions.cxx in ComputeFunctions.cu
#else
  #define BEGIN  0
  #define STEP   1
  #include "RooMath.h"
#endif // #ifdef __CUDACC__

namespace RooBatchCompute {
namespace RF_ARCH {

__global__ void computeAddPdf(Batches batches)
{
  const int nPdfs = batches.getNExtraArgs();
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = batches.extraArg(0)*batches[0][i];
  for (int pdf=1; pdf<nPdfs; pdf++)
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] += batches.extraArg(pdf)*batches[pdf][i];
}

__global__ void computeArgusBG(Batches batches)
{
  Batch m=batches[0], m0=batches[1], c=batches[2], p=batches[3], normVal=batches[4];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double t = m[i]/m0[i];
    const double u = 1 - t*t;
    batches._output[i] = c[i]*u + p[i]*fast_log(u);
  }
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    if (m[i] >= m0[i]) batches._output[i] = 0.0;
    else batches._output[i] = m[i]*fast_exp(batches._output[i])/normVal[i];
  }
}

__global__ void computeBernstein(Batches batches)
{
  const int nCoef = batches.getNExtraArgs()-2;
  const int degree = nCoef-1;
  const double xmin = batches.extraArg(nCoef);
  const double xmax = batches.extraArg(nCoef+1);
  Batch xData=batches[0], normVal=batches[1];
  
  double Binomial[maxExtraArgs] = {1.0}; //Binomial stores values c(degree,i) for i in [0..degree]
  for (int i=1; i<=degree; i++)
    Binomial[i] = Binomial[i-1]*(degree-i+1)/i;

  if (STEP==1)
  {
    double X[bufferSize], _1_X[bufferSize], powX[bufferSize], pow_1_X[bufferSize];
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    {
      powX[i] = pow_1_X[i] = 1.0;
      X[i] = (xData[i]-xmin) / (xmax-xmin);
      _1_X[i] = 1-X[i];
      batches._output[i] = 0.0;
    }

    //raising 1-x to the power of degree
    for (int k=2; k<=degree; k+=2)
      for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
        pow_1_X[i] *= _1_X[i]*_1_X[i];

    if (degree%2 == 1)
      for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
        pow_1_X[i] *= _1_X[i];

    //inverting 1-x ---> 1/(1-x)
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      _1_X[i] = 1/_1_X[i];

    for (int k=0; k<nCoef; k++)
      for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      {
        batches._output[i] += batches.extraArg(k)*Binomial[k]*powX[i]*pow_1_X[i];

        //calculating next power for x and 1-x
        powX[i] *= X[i];
        pow_1_X[i] *= _1_X[i];
      }
  }
  else
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    {
      batches._output[i] = 0.0;
      const double X = (xData[i]-xmin) / (xmax-xmin);
      double powX=1.0, pow_1_X=1.0;
      for (int k=1; k<=degree; k++) pow_1_X *= 1-X;
      const double _1_X = 1/(1-X);
      for (int k=0; k<nCoef; k++)
      {
        batches._output[i] += batches.extraArg(k)*Binomial[k]*powX*pow_1_X;
        powX *= X;
        pow_1_X *= _1_X;
      }
    }

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

__global__ void computeBifurGauss(Batches batches)
{
  Batch X=batches[0], M=batches[1], SL=batches[2], SR=batches[3], normVal=batches[4];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    double arg = X[i]-M[i];
    if (arg<0) arg /= SL[i];
    else arg /= SR[i];
    batches._output[i] = fast_exp(-0.5*arg*arg)/normVal[i];
  }
}

__global__ void computeBreitWigner(Batches batches)
{
  Batch X=batches[0], M=batches[1], W=batches[2], normVal=batches[3];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double arg = X[i]-M[i];
    batches._output[i] = 1 / (arg*arg + 0.25*W[i]*W[i]) / normVal[i];
  }
}
  
__global__ void computeBukin(Batches batches)
{
  Batch X=batches[0], XP=batches[1], SP=batches[2], XI=batches[3], R1=batches[4], R2=batches[5], normVal=batches[6];
  const double r3 = log(2.0);
  const double r6 = exp(-6.0);
  const double r7 = 2*sqrt(2*log(2.0));

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double r1 = XI[i]*fast_isqrt(XI[i]*XI[i]+1);
    const double r4 = 1/fast_isqrt(XI[i]*XI[i]+1);
    const double hp = 1 / (SP[i]*r7);
    const double x1 = XP[i] + 0.5*SP[i]*r7*(r1-1);
    const double x2 = XP[i] + 0.5*SP[i]*r7*(r1+1);

    double r5 = 1.0;
    if (XI[i]>r6 || XI[i]<-r6) r5 = XI[i]/fast_log(r4+XI[i]);

    double factor=1, y=X[i]-x1, Yp=XP[i]-x1, yi=r4-XI[i], rho=R1[i];
    if (X[i]>=x2) 
    {
      factor = -1;
      y = X[i]-x2;
      Yp = XP[i]-x2;
      yi = r4+XI[i];
      rho = R2[i];
    }

    batches._output[i] = rho*y*y/Yp/Yp -r3 + factor*4*r3*y*hp*r5*r4/yi/yi;
    if (X[i]>=x1 && X[i]<x2) 
    {
      batches._output[i] = fast_log(1 + 4*XI[i]*r4*(X[i]-XP[i])*hp) / fast_log(1 +2*XI[i]*( XI[i]-r4 ));
      batches._output[i] *= -batches._output[i]*r3;
    }
    if (X[i]>=x1 && X[i]<x2 && XI[i]<r6 && XI[i]>-r6) 
      batches._output[i] = -4*r3*(X[i]-XP[i])*(X[i]-XP[i])*hp*hp;
    
  }
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = fast_exp(batches._output[i]) / normVal[i];
}

__global__ void computeCBShape(Batches batches)
{
  Batch M=batches[0], M0=batches[1], S=batches[2], A=batches[3], N=batches[4], normVal=batches[5];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double t = (M[i]-M0[i]) / S[i];
    if ( (A[i]>0 && t>=-A[i]) || (A[i]<0 && -t>=A[i]) )
      batches._output[i] = -0.5*t*t;
    else 
    {
      batches._output[i] = N[i] / (N[i] -A[i]*A[i] -A[i]*t);
      batches._output[i] = fast_log(batches._output[i]);
      batches._output[i] *= N[i];
      batches._output[i] -= 0.5*A[i]*A[i];
    }
  }
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = fast_exp(batches._output[i])/normVal[i];
}

__global__ void computeChebychev(Batches batches)
{
  Batch xData=batches[0], normVal=batches[1];
  const int nCoef = batches.getNExtraArgs()-2;
  const double xmin=batches.extraArg(nCoef);
  const double xmax=batches.extraArg(nCoef+1);
  
  if (STEP==1)
  {
    double prev[bufferSize][2], X[bufferSize];

    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP) 
    {
      // set a0-->prev[i][0] and a1-->prev[i][1]
      // and x tranfsformed to range[-1..1]-->X[i]
      prev[i][0] = batches._output[i] = 1.0;
      prev[i][1] = X[i] = 2*( xData[i]-0.5*(xmax+xmin) ) / (xmax-xmin);
    }
    for (int k=0; k<nCoef; k++)
      for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP) 
      {
        batches._output[i] += prev[i][1]*batches.extraArg(k);

        //compute next order
        const double next = 2*X[i]*prev[i][1] -prev[i][0];
        prev[i][0] = prev[i][1];
        prev[i][1] = next;
      }
  }
  else
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    {
      double prev0=1.0, prev1 = 2*( xData[i]-0.5*(xmax+xmin) ) / (xmax-xmin), X=prev1;
      batches._output[i]=1.0;
      for (int k=0; k<nCoef; k++)
      {
        batches._output[i] += prev1*batches.extraArg(k);
        
        //compute next order
        const double next = 2*X*prev1 -prev0;
        prev0 = prev1;
        prev1 = next;
      }
    }
    
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

__global__ void computeChiSquare(Batches batches)
{
  Batch X=batches[0], normVal=batches[1];
  const double ndof = batches.extraArg(0);
  const double gamma = 1/std::tgamma(ndof/2.0);
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i]=gamma/normVal[i];


  constexpr double ln2 = 0.693147180559945309417232121458;
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    double arg = (ndof-2)*fast_log(X[i]) -X[i] -ndof*ln2;
    batches._output[i] *= fast_exp(0.5*arg);
  }
}

__global__ void computeDstD0BG(Batches batches)
{
  Batch DM=batches[0], DM0=batches[1], C=batches[2], A=batches[3], B=batches[4], normVal=batches[5];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double ratio = DM[i] / DM0[i];
    const double arg1 = (DM0[i]-DM[i]) / C[i];
    const double arg2 = A[i]*fast_log(ratio);
    batches._output[i] = (1 -fast_exp(arg1)) * fast_exp(arg2) +B[i]*(ratio-1);
  }

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    if (batches._output[i]<0) batches._output[i] = 0;
    else batches._output[i] /= normVal[i];
}

__global__ void computeExponential(Batches batches)
{
  Batch x=batches[0], c=batches[1], normVal=batches[2];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = fast_exp(x[i]*c[i])/normVal[i];
}

__global__ void computeGamma (Batches batches)
{
  Batch X=batches[0], G=batches[1], B=batches[2], M=batches[3], normVal=batches[4];
  double gamma = -std::lgamma(G[0]);
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    if (X[i] == M[i]) batches._output[i] = (G[i]==1.0)/B[i];
    else if (G.isItVector()) batches._output[i] = -std::lgamma(G[i]);
    else batches._output[i]=gamma;

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    if (X[i] != M[i]) 
    {
      const double invBeta = 1/B[i];
      double arg = (X[i]-M[i])*invBeta;
      batches._output[i] -= arg;
      arg = fast_log(arg);
      batches._output[i] += arg*(G[i]-1);
      batches._output[i] = fast_exp(batches._output[i]);
      batches._output[i] *= invBeta;
    }
    
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Actual computations for the batch evaluation of the Gaussian.
///May vectorise over x, mean, sigma, depending on the types of the inputs.
///\note The output and input spans are assumed to be non-overlapping. If they
///overlap, results will likely be garbage.
__global__ void computeGaussian(Batches batches) 
{
  auto x=batches[0], mean=batches[1], sigma=batches[2], normVal=batches[3];
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double arg = x[i]-mean[i];
    const double halfBySigmaSq = -0.5 / (sigma[i]*sigma[i]);
    batches._output[i] = fast_exp(arg*arg*halfBySigmaSq)/normVal[i];
  }
}

__global__ void computeNegativeLogarithms(Batches batches)
{
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = -fast_log(batches[0][i]);
  // Multiply by weights if they exist
  if (batches.extraArg(0))
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] *= batches[1][i];

}

__global__ void computeJohnson(Batches batches)
{
  Batch mass=batches[0], mu=batches[1], lambda=batches[2], gamma=batches[3], delta=batches[4], normVal=batches[5];
  const double sqrt_twoPi = sqrt(2*M_PI);
  const double massThreshold=batches.extraArg(0);

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double arg = (mass[i]-mu[i]) / lambda[i];
#ifdef R__HAS_VDT
    const double asinh_arg = fast_log(arg + 1/fast_isqrt(arg*arg+1));
#else
    const double asinh_arg = asinh(arg);
#endif
    const double expo = gamma[i] + delta[i]*asinh_arg;
    const double result = delta[i]*fast_exp(-0.5*expo*expo)*fast_isqrt(1. +arg*arg) / (sqrt_twoPi*lambda[i]);

    const double passThrough = mass[i]>=massThreshold;
    batches._output[i] = result*passThrough/normVal[i];
  }
}

/* Actual computation of Landau(x,mean,sigma) in a vectorization-friendly way
 * Code copied from function landau_pdf (math/mathcore/src/PdfFuncMathCore.cxx)
 * and rewritten to enable vectorization.
 */
__global__ void computeLandau(Batches batches)
{
  auto case0 = [](double x)
  {
    const double a1[3] = {0.04166666667,-0.01996527778, 0.02709538966};
    const double u = fast_exp(x+1.0);
    return 0.3989422803*fast_exp( -1/u -0.5*(x+1) )*(1+(a1[0]+(a1[1]+a1[2]*u)*u)*u);
  };
  auto case1 = [](double x)
  {
    constexpr double p1[5] = {0.4259894875,-0.1249762550, 0.03984243700, -0.006298287635,   0.001511162253};
    constexpr double q1[5] = {1.0         ,-0.3388260629, 0.09594393323, -0.01608042283,    0.003778942063};
    const double u = fast_exp(-x-1);
    return fast_exp( -u -0.5*(x+1) ) *
      (p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*x)*x)*x)*x) / (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*x)*x)*x)*x);
  };
  auto case2 = [](double x)
  {
    constexpr double p2[5] = {0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411,   0.0001283617211};
    constexpr double q2[5] = {1.0         , 0.7428795082, 0.3153932961,   0.06694219548,    0.008790609714};
    return (p2[0]+(p2[1]+(p2[2]+(p2[3]+p2[4]*x)*x)*x)*x) / (q2[0]+(q2[1]+(q2[2]+(q2[3]+q2[4]*x)*x)*x)*x);
  };
  auto case3 = [](double x)
  {
    constexpr double p3[5] = {0.1788544503, 0.09359161662,0.006325387654, 0.00006611667319,-0.000002031049101};
    constexpr double q3[5] = {1.0         , 0.6097809921, 0.2560616665,   0.04746722384,    0.006957301675};
    return (p3[0]+(p3[1]+(p3[2]+(p3[3]+p3[4]*x)*x)*x)*x) / (q3[0]+(q3[1]+(q3[2]+(q3[3]+q3[4]*x)*x)*x)*x);
  };
  auto case4 = [](double x)
  {
    constexpr double p4[5] = {0.9874054407, 118.6723273,  849.2794360,   -743.7792444,      427.0262186};
    constexpr double q4[5] = {1.0         , 106.8615961,  337.6496214,    2016.712389,      1597.063511};
    const double u = 1/x;
    return u*u*(p4[0]+(p4[1]+(p4[2]+(p4[3]+p4[4]*u)*u)*u)*u) / (q4[0]+(q4[1]+(q4[2]+(q4[3]+q4[4]*u)*u)*u)*u);
  };
  auto case5 = [](double x)
  {
    constexpr double p5[5] = {1.003675074,  167.5702434,  4789.711289,    21217.86767,     -22324.94910};
    constexpr double q5[5] = {1.0         , 156.9424537,  3745.310488,    9834.698876,      66924.28357};
    const double u = 1/x;
    return u*u*(p5[0]+(p5[1]+(p5[2]+(p5[3]+p5[4]*u)*u)*u)*u) / (q5[0]+(q5[1]+(q5[2]+(q5[3]+q5[4]*u)*u)*u)*u);
  };
  auto case6 = [](double x)
  {
    constexpr double p6[5] = {1.000827619,  664.9143136,  62972.92665,    475554.6998,     -5743609.109};
    constexpr double q6[5] = {1.0         , 651.4101098,  56974.73333,    165917.4725,     -2815759.939};
    const double u = 1/x;
    return u*u*(p6[0]+(p6[1]+(p6[2]+(p6[3]+p6[4]*u)*u)*u)*u) / (q6[0]+(q6[1]+(q6[2]+(q6[3]+q6[4]*u)*u)*u)*u);
  };
  auto case7 = [](double x)
  {
    const double a2[2] = {-1.845568670,-4.284640743};
    const double u = 1 / ( x-x*fast_log(x)/(x+1) );
    return u*u*(1+(a2[0]+a2[1]*u)*u);
  };
  
  Batch X=batches[0], M=batches[1], S=batches[2], normVal=batches[3];
  
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = (X[i]-M[i]) / S[i];

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    if (S[i] <= 0.0) batches._output[i] = 0;
    else if (batches._output[i]<-5.5) batches._output[i] = case0(batches._output[i]);
    else if (batches._output[i]<-1.0) batches._output[i] = case1(batches._output[i]);
    else if (batches._output[i]< 1.0) batches._output[i] = case2(batches._output[i]);
    else if (batches._output[i]< 5.0) batches._output[i] = case3(batches._output[i]);
    else if (batches._output[i]<12.0) batches._output[i] = case4(batches._output[i]);
    else if (batches._output[i]<50.0) batches._output[i] = case5(batches._output[i]);
    else if (batches._output[i]<300.) batches._output[i] = case6(batches._output[i]);
    else                             batches._output[i] = case7(batches._output[i]);

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

__global__ void computeLognormal(Batches batches)
{
  Batch X=batches[0], M0=batches[1], K=batches[2], normVal=batches[3];
  const double rootOf2pi = 2.506628274631000502415765284811;
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    double lnxOverM0 = fast_log(X[i]/M0[i]);
    double lnk = fast_log(K[i]);
    if (lnk<0) lnk = -lnk;
    double arg = lnxOverM0/lnk;
    arg *= -0.5*arg;
    batches._output[i] = fast_exp(arg) / (X[i]*lnk*rootOf2pi*normVal[i]);
  }
}

/* TMath::ASinH(x) needs to be replaced with ln( x + sqrt(x^2+1))
 * argasinh -> the argument of TMath::ASinH()
 * argln -> the argument of the logarithm that replaces AsinH
 * asinh -> the value that the function evaluates to
 *
 * ln is the logarithm that was solely present in the initial
 * formula, that is before the asinh replacement
 */
__global__ void computeNovosibirsk(Batches batches)
{
  Batch X=batches[0], P=batches[1], W=batches[2], T=batches[3], normVal=batches[4];
  constexpr double xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    double argasinh = 0.5*xi*T[i];
    double argln = argasinh + 1/fast_isqrt(argasinh*argasinh +1);
    double asinh = fast_log(argln);

    double argln2 = 1 -(X[i]-P[i])*T[i]/W[i];
    double ln    = fast_log(argln2);
    batches._output[i] = ln/asinh;
    batches._output[i] *= -0.125*xi*xi*batches._output[i];
    batches._output[i] -= 2.0/xi/xi*asinh*asinh;
  }

  //faster if you exponentiate in a seperate loop (dark magic!)
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = fast_exp(batches._output[i])/normVal[i];
}

__global__ void computePoisson(Batches batches)
{
  Batch x=batches[0], mean=batches[1], normVal=batches[2];
  bool protectNegative = batches.extraArg(0);
  bool noRounding      = batches.extraArg(1);
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double x_i = noRounding ? x[i] : floor(x[i]);
    batches._output[i] = std::lgamma(x_i + 1.);
  }

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double x_i = noRounding ? x[i] : floor(x[i]);
    const double logMean = fast_log(mean[i]);
    const double logPoisson = x_i * logMean - mean[i] - batches._output[i];
    batches._output[i] = fast_exp(logPoisson);

    // Cosmetics
    if (x_i < 0)       batches._output[i] = 0;
    else if (x_i == 0) batches._output[i] = 1/fast_exp(mean[i]);

    if (protectNegative && mean[i] < 0) batches._output[i] = 1.E-3;
  }
  
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

__global__ void computePolynomial(Batches batches)
{
  Batch X=batches[0], normVal=batches[1];
  const int nCoef = batches.getNExtraArgs()-1;
  const int lowestOrder = batches.extraArg(nCoef);
  if (nCoef==0)
  {
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] = (lowestOrder>0.0) / normVal[i];
    return;
  }
  else
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] = batches.extraArg(nCoef-1);

  /* Indexes are in range 0..nCoef-1 but coefList[nCoef-1]
   * has already been processed. In order to traverse the list,
   * with step of 2 we have to start at index nCoef-3 and use
   * coefList[k+1] and coefList[k]
   */
  for (int k=nCoef-3; k>=0; k-=2)
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] = X[i]*(batches._output[i]*X[i] + batches.extraArg(k+1)) + batches.extraArg(k);
      
  // If nCoef is even, then the coefList[0] didn't get processed
  if (nCoef%2 == 0)
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] = batches._output[i]*X[i] + batches.extraArg(0);

  //Increase the order of the polynomial, first by myltiplying with X[i]^2
  if (lowestOrder == 0) goto finale;
  for (int k=2; k<=lowestOrder; k+=2)
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] *= X[i]*X[i];

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    if (lowestOrder%2 == 1) batches._output[i] *= X[i];
    batches._output[i] += 1.0;
  }
  
  finale:
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

__global__ void computeProdPdf(Batches batches)
{
  const int nPdfs = batches.extraArg(0);
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] = 1/batches[nPdfs][i]; // normalization
  for (int pdf=0; pdf<nPdfs; pdf++)
    for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
      batches._output[i] *= batches[pdf][i];
}

__global__ void computeVoigtian(Batches batches)
{
  Batch X=batches[0], M=batches[1], W=batches[2], S=batches[3], normVal=batches[4];
  const double invSqrt2 = 0.707106781186547524400844362105;
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
  {
    const double arg = (X[i]-M[i])*(X[i]-M[i]);
    if (S[i]==0.0 && W[i]==0.0) batches._output[i] = 1.0;
    else if (S[i]==0.0) batches._output[i] = 1 / (arg +0.25*W[i]*W[i]);
    else if (W[i]==0.0) batches._output[i] = fast_exp( -0.5*arg/(S[i]*S[i]) );
    else                batches._output[i] = invSqrt2/S[i];
  }

  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    if (S[i]!=0.0 && W[i]!=0.0) 
    {
      if (batches._output[i] < 0) batches._output[i] = -batches._output[i];
      const double factor = W[i]>0.0 ? 0.5 : -0.5;
      std::complex<double> z( batches._output[i]*(X[i]-M[i]) , factor*batches._output[i]*W[i] );
      batches._output[i] *= RooMath::faddeeva(z).real();
    }
    
  for (size_t i=BEGIN; i<batches.getNEvents(); i+=STEP)
    batches._output[i] /= normVal[i];
}

std::vector<void(*)(Batches)> getFunctions()
{
  return {computeAddPdf, computeArgusBG, computeBernstein, computeBifurGauss, computeBreitWigner,
          computeBukin, computeCBShape, computeChebychev, computeChiSquare, computeDstD0BG, 
          computeExponential, computeGaussian, computeGamma, computeJohnson, computeLandau,
          computeLognormal, computeNegativeLogarithms, computeNovosibirsk, computePoisson, 
          computePolynomial, computeProdPdf, computeVoigtian};
}
} // End namespace RF_ARCH
} //End namespace RooBatchCompute
