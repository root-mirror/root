// @(#)root/mathcore:$Id$
// Authors: M. Slawinska   08/2007

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2007 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header source file for class IntegratorMultiDim


#ifndef ROOT_Math_IntegratorMultiDim
#define ROOT_Math_IntegratorMultiDim


#include "Math/IFunctionfwd.h"

#include "Math/AllIntegrationTypes.h"

#include "Math/IntegratorOptions.h"

#include "Math/VirtualIntegrator.h"

#ifndef __CINT__

#include "Math/WrappedFunction.h"

#endif

#include <memory>
#include <string>

namespace ROOT {
namespace Math {

//___________________________________________________________________________________________
/**
   User class for performing multidimensional integration

   By default uses adaptive multi-dimensional integration using the algorithm from Genz Mallik
   implemented in the class ROOT::Math::AdaptiveIntegratorMultiDim otherwise it can uses via the
   plug-in manager the MC integration class (ROOT::Math::GSLMCIntegration) from MathMore.

   @ingroup Integration


 */

class IntegratorMultiDim {

public:

   typedef IntegrationMultiDim::Type Type;   // for the enumerations defining the types


    /** Generic constructor of multi dimensional Integrator. By default uses the Adaptive integration method

       @param type   integration type (adaptive, MC methods, etc..)
       @param absTol desired absolute Error
       @param relTol desired relative Error
       @param ncall  number of function calls (apply only to MC integratioon methods)

       In case no parameter  values are passed the default ones used in IntegratorMultiDimOptions are used
    */
   explicit
   IntegratorMultiDim(IntegrationMultiDim::Type type = IntegrationMultiDim::kDEFAULT, double absTol = -1, double relTol = -1, unsigned int ncall = 0) :
      fIntegrator(0)
   {
       fIntegrator = CreateIntegrator(type, absTol, relTol, ncall);
   }

    /** Generic Constructor of multi dimensional Integrator passing a function. By default uses the adaptive integration method

       @param f      integration function (multi-dim interface)
       @param type   integration type (adaptive, MC methods, etc..)
       @param absTol desired absolute Error
       @param relTol desired relative Error
       @param ncall  number of function calls (apply only to MC integratioon methods)
    */
   explicit
   IntegratorMultiDim(const IMultiGenFunction &f, IntegrationMultiDim::Type type = IntegrationMultiDim::kDEFAULT, double absTol = -1, double relTol = -1, unsigned int ncall = 0) :
      fIntegrator(0)
   {
      fIntegrator = CreateIntegrator(type, absTol, relTol, ncall);
      SetFunction(f);
   }

   // remove template constructor since is ambigous

    /** Template Constructor of multi dimensional Integrator passing a generic function. By default uses the adaptive integration method

       @param f      integration function (generic function implementin operator()(const double *)
       @param dim    function dimension
       @param type   integration type (adaptive, MC methods, etc..)
       @param absTol desired absolute Error
       @param relTol desired relative Error
       @param ncall  number of function calls (apply only to MC integratioon methods)
    */
// this is ambigous
//    template<class Function>
//    IntegratorMultiDim(Function & f, unsigned int dim, IntegrationMultiDim::Type type = IntegrationMultiDim::kADAPTIVE, double absTol = 1.E-9, double relTol = 1E-6, unsigned int ncall = 100000) {
//       fIntegrator = CreateIntegrator(type, absTol, relTol, ncall);
//       SetFunction(f, dim);
//    }

   /// destructor
   virtual ~IntegratorMultiDim() {
      if (fIntegrator) delete fIntegrator;
   }


   // disable copy constructur and assignment operator

private:
   IntegratorMultiDim(const IntegratorMultiDim &) : fIntegrator(0), fFunc(nullptr) {}
   IntegratorMultiDim & operator=(const IntegratorMultiDim &) { return *this; }

public:


   /**
      evaluate the integral with the previously given function between xmin[] and xmax[]
   */
   double Integral(const double* xmin, const double * xmax) {
      return fIntegrator == 0 ? 0 : fIntegrator->Integral(xmin,xmax);
   }

   /// evaluate the integral passing a new function
   double Integral(const IMultiGenFunction &f, const double* xmin, const double * xmax) {
      SetFunction(f);
      return Integral(xmin,xmax);
   }

   /// evaluate the integral passing a new generic function
   template<class Function>
   double Integral(Function & f , unsigned int dim, const double* xmin, const double * xmax) {
      SetFunction<Function>(f,dim);
      return Integral(xmin, xmax);
   }


   /**
       set integration function using a generic function implementing the operator()(double *x)
       The dimension of the function is in this case required
   */
   template <class Function>
   void SetFunction(Function & f, unsigned int dim) {
      fFunc.reset(new  WrappedMultiFunction<Function &> (f, dim) );
      fIntegrator->SetFunction(*fFunc);
   }

   // set the function without cloning it
   void SetFunction(const IMultiGenFunction &f) {
      if (fIntegrator)  fIntegrator->SetFunction(f);
   }

   /// return result of last integration
   double Result() const { return fIntegrator == 0 ? 0 : fIntegrator->Result(); }

   /// return integration error
   double Error() const { return fIntegrator == 0 ? 0 : fIntegrator->Error(); }

   ///  return the Error Status of the last Integral calculation
   int Status() const { return fIntegrator == 0 ? -1 : fIntegrator->Status(); }

   // return number of function evaluations in calculating the integral
   //unsigned int NEval() const { return fNEval; }

   /// set the relative tolerance
   void SetRelTolerance(double relTol) { if (fIntegrator) fIntegrator->SetRelTolerance(relTol); }

   /// set absolute tolerance
   void SetAbsTolerance(double absTol)  { if (fIntegrator) fIntegrator->SetAbsTolerance(absTol); }

   /// set the options
   void SetOptions(const ROOT::Math::IntegratorMultiDimOptions & opt) { if (fIntegrator) fIntegrator->SetOptions(opt); }

   /// retrieve the options
   ROOT::Math::IntegratorMultiDimOptions Options() const { return (fIntegrator) ? fIntegrator->Options() : IntegratorMultiDimOptions(); }

   /// return a pointer to integrator object
   VirtualIntegratorMultiDim * GetIntegrator() { return fIntegrator; }

   /// return name of integrator
   std::string Name() const { return (fIntegrator) ? Options().Integrator() : std::string(""); }

   /// static function to get the enumeration from a string
   static IntegrationMultiDim::Type GetType(const char * name);

   /// static function to get a string from the enumeration
   static std::string GetName(IntegrationMultiDim::Type);

protected:

   VirtualIntegratorMultiDim * CreateIntegrator(IntegrationMultiDim::Type type , double absTol, double relTol, unsigned int ncall);

 private:

   VirtualIntegratorMultiDim * fIntegrator;     // pointer to multi-dimensional integrator base class
   std::unique_ptr<IMultiGenFunction> fFunc;    // pointer to owned function


};

}//namespace Math
}//namespace ROOT

#endif /* ROOT_Math_IntegratorMultiDim */
