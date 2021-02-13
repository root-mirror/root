// @(#)root/mathmore:$Id$
// Authors: L. Moneta, M. Slawinska 10/2007

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2007 ROOT Foundation,  CERN/PH-SFT                   *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for class Integrator
//
//
#ifndef ROOT_Math_Integrator
#define ROOT_Math_Integrator

#include "Math/AllIntegrationTypes.h"

#include "Math/IntegratorOptions.h"

#include "Math/IFunction.h"

#include "Math/VirtualIntegrator.h"

#include <memory>
#include <vector>
#include <string>


/**
@defgroup NumAlgo Numerical Algorithms

Numerical Algorithm classes from the \ref MathCore and \ref MathMore libraries.

@ingroup MathCore
@ingroup MathMore

*/


/**

@defgroup Integration Numerical Integration

Classes for numerical integration of functions.
These classes provide algorithms for integration of one-dimensional functions, with several adaptive and non-adaptive methods
and for integration of multi-dimensional function using an adaptive method or MonteCarlo Integration (GSLMCIntegrator).
The basic classes ROOT::Math::IntegratorOneDim provides a common interface for the one-dimensional methods while the class
ROOT::Math::IntegratorMultiDim provides the interface for the multi-dimensional ones.
The methods can be configured (e.g  setting the default method with its defult parameters) using the ROOT::Math::IntegratorOneDimOptions and
ROOT::Math::IntegratorMultiDimOptions classes.

@ingroup  NumAlgo

*/



namespace ROOT {
namespace Math {




//____________________________________________________________________________________________
/**

User Class for performing numerical integration of a function in one dimension.
It uses the plug-in manager to load advanced numerical integration algorithms from GSL, which reimplements the
algorithms used in the QUADPACK, a numerical integration package written in Fortran.

Various types of adaptive and non-adaptive integration are supported. These include
integration over infinite and semi-infinite ranges and singular integrals.

The integration type is selected using the Integration::type enumeration
in the class constructor.
The default type is adaptive integration with singularity
(ADAPTIVESINGULAR or QAGS in the QUADPACK convention) applying a Gauss-Kronrod 21-point integration rule.
In the case of ADAPTIVE type, the integration rule can also be specified via the
Integration::GKRule. The default rule is 31 points.

In the case of integration over infinite and semi-infinite ranges, the type used is always
ADAPTIVESINGULAR applying a transformation from the original interval into (0,1).

The ADAPTIVESINGULAR type is the most sophicticated type. When performances are
important, it is then recommened to use the NONADAPTIVE type in case of smooth functions or
 ADAPTIVE with a lower Gauss-Kronrod rule.

For detailed description on GSL integration algorithms see the
<A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_16.html#SEC248">GSL Manual</A>.


@ingroup Integration

*/


class IntegratorOneDim {

public:

   typedef IntegrationOneDim::Type Type;   // for the enumerations defining the types

    // constructors


    /**
        Constructor of one dimensional Integrator, default type is adaptive

       @param type   integration type (adaptive, non-adaptive, etc..)
       @param absTol desired absolute Error
       @param relTol desired relative Error
       @param size maximum number of sub-intervals
       @param rule  Gauss-Kronrod integration rule (only for GSL kADAPTIVE type)

       Possible type values are : kGAUSS (simple Gauss method), kADAPTIVE (from GSL), kADAPTIVESINGULAR (from GSL), kNONADAPTIVE (from GSL)
       Possible rule values are kGAUS15 (rule = 1), kGAUS21( rule = 2), kGAUS31(rule =3), kGAUS41 (rule=4), kGAUS51 (rule =5), kGAUS61(rule =6)
       lower rules are indicated for singular functions while higher for smooth functions to get better accuracies

       NOTE: When the default values are passed, the values used are taken from the default defined in ROOT::Math::IntegratorOneDimOptions
    */
    explicit
    IntegratorOneDim(IntegrationOneDim::Type type = IntegrationOneDim::kDEFAULT, double absTol = -1, double relTol = -1, unsigned int size = 0, unsigned int rule = 0) :
       fIntegrator(0), fFunc(0)
   {
      fIntegrator = CreateIntegrator(type, absTol, relTol, size, rule);
   }

   /**
       Constructor of one dimensional Integrator passing a function interface

       @param f      integration function (1D interface). It is copied inside
       @param type   integration type (adaptive, non-adaptive, etc..)
       @param absTol desired absolute tolerance. The algorithm will stop when either the absolute OR the relative tolerance are satisfied.
       @param relTol desired relative tolerance
       @param size maximum number of sub-intervals
       @param rule Gauss-Kronrod integration rule (only for GSL ADAPTIVE type)

       NOTE: When no values are passed, the values used are taken from the default defined in ROOT::Math::IntegratorOneDimOptions
    */
   explicit
   IntegratorOneDim(const IGenFunction &f, IntegrationOneDim::Type type = IntegrationOneDim::kDEFAULT, double absTol = -1, double relTol = -1, unsigned int size = 0, int rule = 0) :
      fIntegrator(0), fFunc(0)
   {
      fIntegrator = CreateIntegrator(type, absTol, relTol, size, rule);
      SetFunction(f,true);
   }

   /**
        Template Constructor of one dimensional Integrator passing a generic function object

       @param f      integration function (any C++ callable object implementing operator()(double x)
       @param type   integration type (adaptive, non-adaptive, etc..)
       @param absTol desired absolute tolerance. The algorithm will stop when either the absolute OR the relative tolerance are satisfied.
       @param relTol desired relative tolerance
       @param size maximum number of sub-intervals
       @param rule Gauss-Kronrod integration rule (only for GSL ADAPTIVE type)

       NOTE: When no values are passed, the values used are taken from the default defined in ROOT::Math::IntegratorOneDimOptions

    */

   template<class Function>
   explicit
   IntegratorOneDim(Function & f, IntegrationOneDim::Type type = IntegrationOneDim::kDEFAULT, double absTol = -1, double relTol = -1, unsigned int size = 0, int rule = 0) :
      fIntegrator(0), fFunc(0)
   {
      fIntegrator = CreateIntegrator(type, absTol, relTol, size, rule);
      SetFunction(f);
   }

   /// destructor (will delete contained pointers)
   virtual ~IntegratorOneDim() {
      if (fIntegrator) delete fIntegrator;
      if (fFunc) delete fFunc;
   }

   // disable copy constructur and assignment operator

private:
   IntegratorOneDim(const IntegratorOneDim &) : fIntegrator(0), fFunc(0) {}
   IntegratorOneDim & operator=(const IntegratorOneDim &) { return *this; }

public:


   // template methods for generic functors

   /**
      method to set the a generic integration function
      @param f integration function. The function type must implement the assigment operator, <em>  double  operator() (  double  x ) </em>

   */


   template<class Function>
   inline void SetFunction(Function & f);

   /**
       set one dimensional function for 1D integration
    */
   void SetFunction  (const IGenFunction &f, bool copy = false) {
      if (!fIntegrator) return;
      if (copy) {
         if (fFunc) delete fFunc;
         fFunc = f.Clone();
         fIntegrator->SetFunction(*fFunc);
         return;
      }
      fIntegrator->SetFunction(f);
   }


   /**
      Set integration function from a multi-dim function type.
      Can be used in case of having 1D function implementing the generic interface
       @param f      integration function
       @param icoord index of coordinate on which the integration is performed
       @param x  array of the passed variables values. In case of dim=1 a 0 can be passed
   */
   void SetFunction(const IMultiGenFunction &f, unsigned int icoord , const double * x );

    // integration methods using a function

    /**
       evaluate the Integral of a function f over the defined interval (a,b)
       @param f integration function. The function type must be a C++ callable object implementing operator()(double x)
       @param a lower value of the integration interval
       @param b upper value of the integration interval
    */
   template<class Function>
   double Integral(Function & f, double a, double b);


    /**
       evaluate the Integral of a function f over the defined interval (a,b)
       @param f integration function. The function type must implement the mathlib::IGenFunction interface
       @param a lower value of the integration interval
       @param b upper value of the integration interval
    */
   double Integral(const IGenFunction & f, double a, double b) {
     SetFunction(f,false);
      return Integral(a,b);
   }


    /**
      evaluate the Integral of a function f over the infinite interval (-inf,+inf)
       @param f integration function. The function type must be a C++ callable object implementing operator()(double x)
    */
//    template<class Function>
//    double Integral(const Function & f);

   /**
      evaluate the Integral of a function f over the infinite interval (-inf,+inf)
      @param f integration function. The function type must implement the mathlib::IGenFunction interface
   */
   double Integral(const IGenFunction & f) {
      SetFunction(f,false);
      return Integral();
   }


    /**
      evaluate the Integral of a function f over the semi-infinite interval (a,+inf)
      @param f integration function. The function type must be a C++ callable object implementing operator()(double x)
      @param a lower value of the integration interval
    */
//    template<class Function>
//    double IntegralUp(Function & f, double a);

   /**
      evaluate the Integral of a function f over the semi-infinite interval (a,+inf)
      @param f integration function. The function type must implement the mathlib::IGenFunction interface
      @param a lower value of the integration interval

   */
   double IntegralUp(const IGenFunction & f, double a ) {
     SetFunction(f,false);
      return IntegralUp(a);
   }

    /**
      evaluate the Integral of a function f over the over the semi-infinite interval (-inf,b)
      @param f integration function. The function type must be a C++ callable object implementing operator()(double x)
      @param b upper value of the integration interval
    */
//    template<class Function>
//    double IntegralLow(Function & f, double b);

   /**
      evaluate the Integral of a function f over the over the semi-infinite interval (-inf,b)
      @param f integration function. The function type must implement the mathlib::IGenFunction interface
      @param b upper value of the integration interval
   */
   double IntegralLow(const IGenFunction & f, double b ) {
     SetFunction(f,false);
      return IntegralLow(b);
   }

   /**
      evaluate the Integral of a function f with known singular points over the defined Integral (a,b)
      @param f integration function. The function type must be a C++ callable object implementing operator()(double x)
      @param pts vector containing both the function singular points and the lower/upper edges of the interval. The vector must have as first element the lower edge of the integration Integral ( \a a) and last element the upper value.

   */
   template<class Function>
   double Integral(Function & f, const std::vector<double> & pts );

   /**
      evaluate the Integral of a function f with known singular points over the defined Integral (a,b)
      @param f integration function. The function type must implement the mathlib::IGenFunction interface
      @param pts vector containing both the function singular points and the lower/upper edges of the interval. The vector must have as first element the lower edge of the integration Integral ( \a a) and last element the upper value.

   */
   double Integral(const IGenFunction & f, const std::vector<double> & pts ) {
     SetFunction(f,false);
      return Integral(pts);
   }

   /**
      evaluate the Cauchy principal value of the integral of  a function f over the defined interval (a,b) with a singularity at c
      @param f integration function. The function type must be a C++ callable object implementing operator()(double x)
      @param a lower value of the integration interval
      @param b upper value of the integration interval
      @param c position of singularity

   */
   template<class Function>
   double IntegralCauchy(Function & f, double a, double b, double c);

   /**
      evaluate the Cauchy principal value of the integral of  a function f over the defined interval (a,b) with a singularity at c
       @param f integration function. The function type must implement the mathlib::IGenFunction interface
       @param a lower value of the integration interval
       @param b upper value of the integration interval
       @param c position of singularity

   */
   double IntegralCauchy(const IGenFunction & f, double a, double b, double c) {
     SetFunction(f,false);
      return IntegralCauchy(a,b,c);
   }



   // integration method using cached function

   /**
      evaluate the Integral over the defined interval (a,b) using the function previously set with Integrator::SetFunction method
      @param a lower value of the integration interval
      @param b upper value of the integration interval
   */

   double Integral(double a, double b) {
      return fIntegrator == 0 ? 0 : fIntegrator->Integral(a,b);
   }


   /**
      evaluate the Integral over the infinite interval (-inf,+inf) using the function previously set with Integrator::SetFunction method.
   */

   double Integral( ) {
      return fIntegrator == 0 ? 0 : fIntegrator->Integral();
   }

   /**
      evaluate the Integral of a function f over the semi-infinite interval (a,+inf) using the function previously set with Integrator::SetFunction method.
      @param a lower value of the integration interval
   */
   double IntegralUp(double a ) {
      return fIntegrator == 0 ? 0 : fIntegrator->IntegralUp(a);
   }

   /**
      evaluate the Integral of a function f over the over the semi-infinite interval (-inf,b) using the function previously set with Integrator::SetFunction method.
      @param b upper value of the integration interval
   */
   double IntegralLow( double b ) {
      return fIntegrator == 0 ? 0 : fIntegrator->IntegralLow(b);
   }
   /**
       define operator() for IntegralLow
    */
   double operator() (double x) {
      return IntegralLow(x);
   }


   /**
      evaluate the Integral over the defined interval (a,b) using the function previously set with Integrator::SetFunction method. The function has known singular points.
      @param pts vector containing both the function singular points and the lower/upper edges of the interval. The vector must have as first element the lower edge of the integration Integral ( \a a) and last element the upper value.

   */
   double Integral( const std::vector<double> & pts) {
      return fIntegrator == 0 ? 0 : fIntegrator->Integral(pts);
   }

   /**
      evaluate the Cauchy principal value of the integral of  a function f over the defined interval (a,b) with a singularity at c

   */
   double IntegralCauchy(double a, double b, double c) {
      return fIntegrator == 0 ? 0 : fIntegrator->IntegralCauchy(a,b,c);
   }

   /**
      return  the Result of the last Integral calculation
   */
   double Result() const { return fIntegrator == 0 ? 0 : fIntegrator->Result(); }

   /**
      return the estimate of the absolute Error of the last Integral calculation
   */
   double Error() const { return fIntegrator == 0 ? 0 : fIntegrator->Error(); }

   /**
      return the Error Status of the last Integral calculation
   */
   int Status() const { return fIntegrator == 0 ? -1 : fIntegrator->Status(); }

   /**
      return number of function evaluations in calculating the integral
      (if integrator do not implement this function returns -1)
   */
   int NEval() const { return fIntegrator == 0 ? -1 : fIntegrator->NEval(); }


   // setter for control Parameters  (getters are not needed so far )

   /**
      set the desired relative Error
   */
   void SetRelTolerance(double relTolerance) { if (fIntegrator) fIntegrator->SetRelTolerance(relTolerance); }


   /**
      set the desired absolute Error
   */
   void SetAbsTolerance(double absTolerance) { if (fIntegrator) fIntegrator->SetAbsTolerance(absTolerance); }

   /**
      return a pointer to integrator object
   */
   VirtualIntegratorOneDim * GetIntegrator() { return fIntegrator; }

   /**
       set the options
   */
   void SetOptions(const ROOT::Math::IntegratorOneDimOptions & opt) { if (fIntegrator) fIntegrator->SetOptions(opt); }

   /**
       retrieve the options
   */
   ROOT::Math::IntegratorOneDimOptions Options() const { return (fIntegrator) ? fIntegrator->Options() : IntegratorOneDimOptions(); }

   /// return name of integrator
   std::string Name() const { return (fIntegrator) ? Options().Integrator() : std::string(""); }

   /// static function to get the enumeration from a string
   static IntegrationOneDim::Type GetType(const char * name);

   /// static function to get a string from the enumeration
   static std::string GetName(IntegrationOneDim::Type);


protected:

   VirtualIntegratorOneDim * CreateIntegrator(IntegrationOneDim::Type type , double absTol, double relTol, unsigned int size, int rule);

private:

   VirtualIntegratorOneDim * fIntegrator;   // pointer to integrator interface class
   IGenFunction            * fFunc;         // pointer to owned function

};


   typedef IntegratorOneDim Integrator;


} // namespace Math
} // namespace ROOT


#ifndef __CINT__


#include "Math/WrappedFunction.h"

template<class Function>
void ROOT::Math::IntegratorOneDim::SetFunction( Function & f) {
  ::ROOT::Math::WrappedFunction<Function &> wf(f);
  // need to copy the wrapper function, the instance created here will be deleted after SetFunction()
  SetFunction(wf, true);
}

template<class Function>
double ROOT::Math::IntegratorOneDim::Integral(Function & f, double a, double b) {
   ::ROOT::Math::WrappedFunction< Function &> wf(f);
   SetFunction(wf,false); // no copy is needed in this case
   return Integral(a,b);
}

// remove because can create ambiguities

// template<class Function>
// double ROOT::Math::IntegratorOneDim::Integral(const  Function & f) {
//    ROOT::Math::WrappedFunction<const Function &> wf(f);
//    SetFunction(wf,false); // no copy is needed in this case
//    return Integral();
// }

// template<class Function>
// double ROOT::Math::IntegratorOneDim::IntegralLow(Function  & f, double x) {
//    ROOT::Math::WrappedFunction< Function &> wf(f);
//    SetFunction(wf,false); // no copy is needed in this case
//    return IntegralLow(x);
// }

// template<class Function>
// double ROOT::Math::IntegratorOneDim::IntegralUp(Function & f, double x) {
//    ROOT::Math::WrappedFunction<Function &> wf(f);
//    SetFunction(wf,false); // no copy is needed in this case
//    return IntegralUp(x);
// }

template<class Function>
double ROOT::Math::IntegratorOneDim::Integral(Function & f, const std::vector<double> & pts) {
   ::ROOT::Math::WrappedFunction<Function &> wf(f);
   SetFunction(wf,false); // no copy is needed in this case
   return Integral(pts);
}

template<class Function>
double ROOT::Math::IntegratorOneDim::IntegralCauchy(Function & f, double a, double b, double c) {
   ::ROOT::Math::WrappedFunction<Function & > wf(f);
   SetFunction(wf,false); // no copy is needed in this case
   return IntegralCauchy(a,b,c);
}


#endif



#endif /* ROOT_Math_Integrator */
