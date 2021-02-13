// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Nov 14 14:20:07 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class IParamFunction

#ifndef ROOT_Math_IParamFunction
#define ROOT_Math_IParamFunction

#include "Math/IFunction.h"

#include "Math/IParamFunctionfwd.h"

#include "Math/Util.h"

#include <cassert>
#include <string>

/**
   @defgroup ParamFunc Parameteric Function Evaluation Interfaces.
   Interfaces classes for evaluation of parametric functions
   @ingroup CppFunctions
*/


namespace ROOT {

   namespace Math {


//___________________________________________________________________
      /**
          Documentation for the abstract class IBaseParam.
          It defines the interface for dealing with the function parameters
          This is used only for internal convinience, to avoid redefining the Parameter API
          for the one and the multi-dim functions.
          Concrete class should derive from ROOT::Math::IParamFunction and not from this class.

          @ingroup  ParamFunc
      */

      class IBaseParam  {

      public:


         /**
            Virtual Destructor (no operations)
         */
         virtual ~IBaseParam()  {}


         /**
            Access the parameter values
         */
         virtual const double *Parameters() const = 0;

         /**
            Set the parameter values
            @param p vector of doubles containing the parameter values.

            to be defined:  can user change number of params ? At the moment no.

         */
         virtual void SetParameters(const double *p) = 0;


         /**
            Return the number of Parameters
         */
         virtual unsigned int NPar() const = 0;

         /**
            Return the name of the i-th parameter (starting from zero)
            Overwrite if want to avoid the default name ("Par_0, Par_1, ...")
          */
         virtual std::string ParameterName(unsigned int i) const
         {
            assert(i < NPar());
            return "Par_" + Util::ToString(i);
         }


      };

//___________________________________________________________________
      /**
         IParamFunction interface (abstract class) describing multi-dimensional parameteric functions
         It is a derived class from ROOT::Math::IBaseFunctionMultiDim and
         ROOT::Math::IBaseParam

         Provides the interface for evaluating a function passing a coordinate vector and a parameter vector.

         @ingroup  ParamFunc
      */

      template<class T>
      class IParametricFunctionMultiDimTempl: virtual public IBaseFunctionMultiDimTempl<T>,
         virtual public IBaseParam {
      public:

         typedef IBaseFunctionMultiDimTempl<T>  BaseFunc;

         /**
         Evaluate function at a point x and for given parameters p.
         This method does not change the internal status of the function (internal parameter values).
         If for some reason one prefers caching the parameter values, SetParameters(p) and then operator()(x) should be
         called.
         Use the pure virtual function DoEvalPar to implement it
         */

         /* Reimplementation instead of using BaseParamFunc::operator();
         until the bug in VS is fixed */
         T operator()(const T *x, const double   *p) const
         {
            return DoEvalPar(x, p);
         }

         T operator()(const T *x) const
         {
            return DoEval(x);
         }

      private:
         /**
            Implementation of the evaluation function using the x values and the parameters.
            Must be implemented by derived classes
         */
         virtual T DoEvalPar(const T *x, const double *p) const = 0;

         /**
            Implement the ROOT::Math::IBaseFunctionMultiDim interface DoEval(x) using the cached parameter values
         */
         virtual T DoEval(const T *x) const
         {
            return DoEvalPar(x, Parameters());
         }
      };


//___________________________________________________________________
      /**
         Specialized IParamFunction interface (abstract class) for one-dimensional parametric functions
         It is a derived class from ROOT::Math::IBaseFunctionOneDim and
         ROOT::Math::IBaseParam

         @ingroup  ParamFunc
      */

      class IParametricFunctionOneDim :
         virtual public IBaseFunctionOneDim,
         public IBaseParam {


      public:

         typedef IBaseFunctionOneDim   BaseFunc;


         using BaseFunc::operator();

         /**
            Evaluate function at a point x and for given parameters p.
            This method does not change the internal status of the function (internal parameter values).
            If for some reason one prefers caching the parameter values, SetParameters(p) and then operator()(x) should be
            called.
            Use the pure virtual function DoEvalPar to implement it
         */
         double operator()(double x, const double   *p) const
         {
            return DoEvalPar(x, p);
         }


         /**
            multidim-like interface
         */
         double operator()(const double *x, const double   *p) const
         {
            return DoEvalPar(*x, p);
         }

      private:

         /**
            Implementation of the evaluation function using the x value and the parameters.
            Must be implemented by derived classes
         */
         virtual double DoEvalPar(double x, const double *p) const = 0;

         /**
            Implement the ROOT::Math::IBaseFunctionOneDim interface DoEval(x) using the cached parameter values
         */
         virtual double DoEval(double x) const
         {
            return DoEvalPar(x, Parameters());
         }

      };



//_______________________________________________________________________________
      /**
         Interface (abstract class) for parametric gradient multi-dimensional functions providing
         in addition to function evaluation with respect to the coordinates
         also the gradient with respect to the parameters, via the method ParameterGradient.

         It is a derived class from ROOT::Math::IParametricFunctionMultiDim.

         The pure private virtual method DoParameterGradient must be implemented by the derived classes
         in addition to those inherited by the base abstract classes.

         @ingroup  ParamFunc
      */

      template<class T>
      class IParametricGradFunctionMultiDimTempl: virtual public IParametricFunctionMultiDimTempl<T>, virtual public IBaseParam {
      public:

         using  BaseParamFunc = IParametricFunctionMultiDimTempl<T>;
         using BaseGradFunc = IGradientFunctionMultiDimTempl<T>;
         using BaseFunc = typename IParametricFunctionMultiDimTempl<T>::BaseFunc;


         /**
            Virtual Destructor (no operations)
         */
         virtual ~IParametricGradFunctionMultiDimTempl()  {}


         /* Reimplementation instead of using BaseParamFunc::operator();
         until the bug in VS is fixed */
         T operator()(const T *x, const double   *p) const
         {
            return DoEvalPar(x, p);
         }

         T operator()(const T *x) const
         {
            return DoEval(x);
         }

         /**
            Evaluate the all the derivatives (gradient vector) of the function with respect to the parameters at a point x.
            It is optional to be implemented by the derived classes for better efficiency
         */
         virtual void ParameterGradient(const T *x, const double *p, T *grad) const
         {
            unsigned int npar = NPar();
            for (unsigned int ipar  = 0; ipar < npar; ++ipar)
               grad[ipar] = DoParameterDerivative(x, p, ipar);
         }

         /**
            Evaluate the partial derivative w.r.t a parameter ipar from values and parameters
          */
         T ParameterDerivative(const T *x, const double *p, unsigned int ipar = 0) const
         {
            return DoParameterDerivative(x, p, ipar);
         }

         /**
            Evaluate all derivatives using cached parameter values
         */
         void ParameterGradient(const T *x, T *grad) const { return ParameterGradient(x, Parameters(), grad); }
         /**
            Evaluate partial derivative using cached parameter values
         */
         T ParameterDerivative(const T *x, unsigned int ipar = 0) const
         {
            return DoParameterDerivative(x, Parameters() , ipar);
         }

      private:

         /**
            Evaluate the partial derivative w.r.t a parameter ipar , to be implemented by the derived classes
          */
         virtual T DoParameterDerivative(const T *x, const double *p, unsigned int ipar) const = 0;
         virtual T DoEvalPar(const T *x, const double *p) const = 0;
         virtual T DoEval(const T *x) const
         {
            return DoEvalPar(x, Parameters());
         }
      };

//_______________________________________________________________________________
      /**
         Interface (abstract class) for parametric one-dimensional gradient functions providing
         in addition to function evaluation with respect the coordinates
         also the gradient with respect to the parameters, via the method ParameterGradient.

         It is a derived class from ROOT::Math::IParametricFunctionOneDim.

         The pure private virtual method DoParameterGradient must be implemented by the derived classes
         in addition to those inherited by the base abstract classes.

         @ingroup  ParamFunc
      */

      class IParametricGradFunctionOneDim :
         public IParametricFunctionOneDim
//         ,public IGradientFunctionOneDim
      {

      public:

         typedef IParametricFunctionOneDim            BaseParamFunc;
         typedef IGradientFunctionOneDim              BaseGradFunc;
         typedef IParametricFunctionOneDim::BaseFunc  BaseFunc;


         /**
            Virtual Destructor (no operations)
         */
         virtual ~IParametricGradFunctionOneDim()  {}


         using BaseParamFunc::operator();

         /**
            Evaluate the derivatives of the function with respect to the parameters at a point x.
            It is optional to be implemented by the derived classes for better efficiency if needed
         */
         virtual void ParameterGradient(double x , const double *p, double *grad) const
         {
            unsigned int npar = NPar();
            for (unsigned int ipar  = 0; ipar < npar; ++ipar)
               grad[ipar] = DoParameterDerivative(x, p, ipar);
         }

         /**
            Evaluate all derivatives using cached parameter values
         */
         void ParameterGradient(double  x , double *grad) const
         {
            return ParameterGradient(x, Parameters(), grad);
         }

         /**
            Compatibility interface with multi-dimensional functions
         */
         void ParameterGradient(const double *x , const double *p, double *grad) const
         {
            ParameterGradient(*x, p, grad);
         }

         /**
            Evaluate all derivatives using cached parameter values (multi-dim like interface)
         */
         void ParameterGradient(const double *x , double *grad) const
         {
            return ParameterGradient(*x, Parameters(), grad);
         }


         /**
            Partial derivative with respect a parameter
          */
         double ParameterDerivative(double x, const double *p, unsigned int ipar = 0) const
         {
            return DoParameterDerivative(x, p, ipar);
         }

         /**
            Evaluate partial derivative using cached parameter values
         */
         double ParameterDerivative(double x, unsigned int ipar = 0) const
         {
            return DoParameterDerivative(x, Parameters() , ipar);
         }

         /**
            Partial derivative with respect a parameter
            Compatibility interface with multi-dimensional functions
         */
         double ParameterDerivative(const double *x, const double *p, unsigned int ipar = 0) const
         {
            return DoParameterDerivative(*x, p, ipar);
         }


         /**
            Evaluate partial derivative using cached parameter values (multi-dim like interface)
         */
         double ParameterDerivative(const double *x, unsigned int ipar = 0) const
         {
            return DoParameterDerivative(*x, Parameters() , ipar);
         }



      private:


         /**
            Evaluate the gradient, to be implemented by the derived classes
          */
         virtual double DoParameterDerivative(double x, const double *p, unsigned int ipar) const = 0;


      };




   } // end namespace Math

} // end namespace ROOT



#endif /* ROOT_Math_IParamFunction */
