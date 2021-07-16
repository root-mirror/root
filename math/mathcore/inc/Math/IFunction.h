// @(#)root/mathcore:$Id$
// Authors: L. Moneta    11/2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for function interfaces
//
// Generic Interfaces for one or  multi-dimensional functions
//
// Created by: Lorenzo Moneta  : Wed Nov 13 2006
//
//
#ifndef ROOT_Math_IFunction
#define ROOT_Math_IFunction

/**
@defgroup CppFunctions Function Classes and Interfaces

 Interfaces (abstract classes) and Base classes used in MathCore and MathMore numerical methods
 for describing function classes. They define function and gradient evaluation and as well the
 functionality for dealing with parameters in the case of parametric functions which are used for
 fitting and data modeling.
 Included are also adapter classes, such as  functors, to wrap generic callable C++ objects
 in the desired interface.

@ingroup MathCore
*/

//typedefs and tags definitions
#include "Math/IFunctionfwd.h"


namespace ROOT {
   namespace Math {

      /**
         @defgroup GenFunc Generic Function Evaluation Interfaces
         Interface classes for evaluation of function object classes in one or multi-dimensions.
         @ingroup CppFunctions
      */

//___________________________________________________________________________________
      /**
          Documentation for the abstract class IBaseFunctionMultiDim.
          Interface (abstract class) for generic functions objects of multi-dimension
          Provides a method to evaluate the function given a vector of coordinate values,
          by implementing operator() (const double *).
          In addition it defines the interface for copying functions via the pure virtual method Clone()
          and the interface for getting the function dimension via the NDim() method.
          Derived classes must implement the pure private virtual method DoEval(const double *) for the
          function evaluation in addition to NDim() and Clone().

          @ingroup  GenFunc
      */

      template<class T>
      class IBaseFunctionMultiDimTempl {

      public:

         typedef T BackendType;
         typedef  IBaseFunctionMultiDimTempl<T> BaseFunc;


         IBaseFunctionMultiDimTempl() {}

         /**
            virtual destructor
          */
         virtual ~IBaseFunctionMultiDimTempl() {}

         /**
             Clone a function.
             Each derived class must implement their version of the Clone method
         */
         virtual IBaseFunctionMultiDimTempl<T> *Clone() const = 0;

         /**
            Retrieve the dimension of the function
          */
         virtual unsigned int NDim() const = 0;

         /**
             Evaluate the function at a point x[].
             Use the pure virtual private method DoEval which must be implemented by the sub-classes
         */
         T operator()(const T *x) const
         {
            return DoEval(x);
         }

#ifdef LATER
         /**
            Template method to eveluate the function using the begin of an iterator
            User is responsible to provide correct size for the iterator
         */
         template <class Iterator>
         T operator()(const Iterator it) const
         {
            return DoEval(&(*it));
         }
#endif


      private:


         /**
            Implementation of the evaluation function. Must be implemented by derived classes
         */
         virtual T DoEval(const T *x) const = 0;


      };


//___________________________________________________________________________________
      /**
          Interface (abstract class) for generic functions objects of one-dimension
          Provides a method to evaluate the function given a value (simple double)
          by implementing operator() (const double ).
          In addition it defines the interface for copying functions via the pure virtual method Clone().
          Derived classes must implement the pure virtual private method DoEval(double ) for the
          function evaluation in addition to  Clone().
          An interface for evaluating the function passing a vector (like for multidim functions) is also
          provided

          @ingroup  GenFunc
      */
      class IBaseFunctionOneDim {

      public:

         typedef  IBaseFunctionOneDim BaseFunc;

         IBaseFunctionOneDim() {}

         /**
            virtual destructor
          */
         virtual ~IBaseFunctionOneDim() {}

         /**
             Clone a function.
             Each derived class will implement their version of the provate DoClone method
         */
         virtual IBaseFunctionOneDim *Clone() const = 0;

         /**
             Evaluate the function at a point x
             Use the  a pure virtual private method DoEval which must be implemented by sub-classes
         */
         double operator()(double x) const
         {
            return DoEval(x);
         }

         /**
             Evaluate the function at a point x[].
             Compatible method with multi-dimensional functions
         */
         double operator()(const double *x) const
         {
            return DoEval(*x);
         }



      private:

         // use private virtual inheritance

         /// implementation of the evaluation function. Must be implemented by derived classes
         virtual double DoEval(double x) const = 0;

      };


//-------- GRAD  functions---------------------------

//___________________________________________________________________________________
      /**
         Gradient interface (abstract class) defining the signature for calculating the gradient of a
         multi-dimensional function.
         Three methods are provided:
         - Gradient(const double *x, double * grad) evaluate the full gradient vector at the vector value x
         - Derivative(const double * x, int icoord) evaluate the partial derivative for the icoord coordinate
         - FdF(const double *x, double &f, double * g) evaluate at the same time gradient and function/

         Concrete classes should derive from ROOT::Math::IGradientFunctionMultiDim and not from this class.

         @ingroup  GenFunc
       */

      template <class T>
      class IGradientMultiDimTempl {

      public:

         /// virual destructor
         virtual ~IGradientMultiDimTempl() {}

         /**
             Evaluate all the vector of function derivatives (gradient)  at a point x.
             Derived classes must re-implement if it is more efficient than evaluting one at a time
         */
         virtual void Gradient(const T *x, T *grad) const = 0;

         /**
            Return the partial derivative with respect to the passed coordinate
         */
         T Derivative(const T *x, unsigned int icoord = 0) const { return DoDerivative(x, icoord); }
         /// In some cases, the derivative algorithm will use information from the previous step, these can be passed
         /// in with this overload. The `previous_*` arrays can also be used to return second derivative and step size
         /// so that these can be passed forward again as well at the call site, if necessary.
         T Derivative(const T *x, unsigned int icoord, T *previous_grad, T *previous_g2,
                      T *previous_gstep) const
         {
            return DoDerivative(x, icoord, previous_grad, previous_g2, previous_gstep);
         }

         /**
             Optimized method to evaluate at the same time the function value and derivative at a point x.
             Often both value and derivatives are needed and it is often more efficient to compute them at the same time.
             Derived class should implement this method if performances play an important role and if it is faster to
             evaluate value and derivative at the same time

         */
         virtual void FdF(const T *x, T &f, T *df) const = 0;

      private:


         /**
            function to evaluate the derivative with respect each coordinate. To be implemented by the derived class
         */
         virtual T DoDerivative(const T *x, unsigned int icoord) const = 0;
         /// In some cases, the derivative algorithm will use information from the previous step, these can be passed
         /// in with this overload. The `previous_*` arrays can also be used to return second derivative and step size
         /// so that these can be passed forward again as well at the call site, if necessary.
         virtual T DoDerivative(const T *x, unsigned int icoord, T * /*previous_grad*/, T * /*previous_g2*/,
                                T * /*previous_gstep*/) const
         {
            return DoDerivative(x, icoord);
         };
      };

//___________________________________________________________________________________
      /**
         Specialized Gradient interface(abstract class)  for one dimensional functions
         It provides a method to evaluate the derivative of the function, Derivative and a
         method to evaluate at the same time the function and the derivative FdF

         Concrete classes should derive from ROOT::Math::IGradientFunctionOneDim and not from this class.

         @ingroup  GenFunc
       */
      class IGradientOneDim {

      public:

         /// virtual destructor
         virtual ~IGradientOneDim() {}

         /**
            Return the derivative of the function at a point x
            Use the private method DoDerivative
         */
         double Derivative(double x) const
         {
            return DoDerivative(x);
         }


         /**
             Optimized method to evaluate at the same time the function value and derivative at a point x.
             Often both value and derivatives are needed and it is often more efficient to compute them at the same time.
             Derived class should implement this method if performances play an important role and if it is faster to
             evaluate value and derivative at the same time

         */
         virtual void FdF(double x, double &f, double &df) const = 0;


         /**
            Compatibility method with multi-dimensional interface for partial derivative
          */
         double Derivative(const double *x) const
         {
            return DoDerivative(*x);
         }

         /**
            Compatibility method with multi-dimensional interface for Gradient
          */
         void Gradient(const double *x, double *g) const
         {
            g[0] = DoDerivative(*x);
         }

         /**
            Compatibility method with multi-dimensional interface for Gradient and function evaluation
          */
         void FdF(const double *x, double &f, double *df) const
         {
            FdF(*x, f, *df);
         }



      private:


         /**
            function to evaluate the derivative with respect each coordinate. To be implemented by the derived class
         */
         virtual  double  DoDerivative(double x) const = 0;

      };

//___________________________________________________________________________________
      /**
         Interface (abstract class) for multi-dimensional functions providing a gradient calculation.
         It implements both the ROOT::Math::IBaseFunctionMultiDimTempl and
         ROOT::Math::IGradientMultiDimTempl interfaces.
         The method ROOT::Math::IFunction::Gradient calculates the full gradient vector,
         ROOT::Math::IFunction::Derivative calculates the partial derivative for each coordinate and
         ROOT::Math::Fdf calculates the gradient and the function value at the same time.
         The pure private virtual method DoDerivative() must be implemented by the derived classes, while
         Gradient and FdF are by default implemented using DoDerivative, butthey  can be overloaded by the
         derived classes to improve the efficiency in the derivative calculation.

         @ingroup  GenFunc
      */

      template <class T>
      class IGradientFunctionMultiDimTempl : virtual public IBaseFunctionMultiDimTempl<T>,
                                             public IGradientMultiDimTempl<T> {

      public:
         typedef IBaseFunctionMultiDimTempl<T> BaseFunc;
         typedef IGradientMultiDimTempl<T> BaseGrad;

         /**
            Virtual Destructor (no operations)
         */
         virtual ~IGradientFunctionMultiDimTempl() {}

         /**
            Evaluate all the vector of function derivatives (gradient)  at a point x.
            Derived classes must re-implement it if more efficient than evaluting one at a time
         */
         virtual void Gradient(const T *x, T *grad) const
         {
            unsigned int ndim = NDim();
            for (unsigned int icoord  = 0; icoord < ndim; ++icoord)
               grad[icoord] = BaseGrad::Derivative(x, icoord);
         }

         /// In some cases, the gradient algorithm will use information from the previous step, these can be passed
         /// in with this overload. The `previous_*` arrays can also be used to return second derivative and step size
         /// so that these can be passed forward again as well at the call site, if necessary.
         virtual void Gradient(const T *x, T *grad, T *previous_grad, T *previous_g2, T *previous_gstep) const
         {
            unsigned int ndim = NDim();
            for (unsigned int icoord  = 0; icoord < ndim; ++icoord)
               grad[icoord] = BaseGrad::Derivative(x, icoord, previous_grad, previous_g2, previous_gstep);
         }

         using  BaseFunc::NDim;

         /**
            Optimized method to evaluate at the same time the function value and derivative at a point x.
            Often both value and derivatives are needed and it is often more efficient to compute them at the same time.
            Derived class should implement this method if performances play an important role and if it is faster to
            evaluate value and derivative at the same time
         */
         virtual void FdF(const T *x, T &f, T *df) const
         {
            f = BaseFunc::operator()(x);
            Gradient(x, df);
         }

         virtual bool returnsInMinuit2ParameterSpace() const { return false; }
      };

//___________________________________________________________________________________
      /**
         Interface (abstract class) for one-dimensional functions providing a gradient calculation.
         It implements both the ROOT::Math::IBaseFunctionOneDim and
         ROOT::Math::IGradientOneDim interfaces.
         The method  ROOT::Math::IFunction::Derivative calculates the derivative  and
         ROOT::Math::Fdf calculates the derivative and the function values at the same time.
         The pure private virtual method DoDerivative() must be implemented by the derived classes, while
         FdF is by default implemented using DoDerivative, but it can be overloaded by the
         derived classes to improve the efficiency in the derivative calculation.


         @ingroup  GenFunc
      */
      //template <>
      class IGradientFunctionOneDim :
         virtual public IBaseFunctionOneDim ,
         public IGradientOneDim {


      public:

         typedef IBaseFunctionOneDim BaseFunc;
         typedef IGradientOneDim BaseGrad;


         /**
             Virtual Destructor (no operations)
         */
         virtual ~IGradientFunctionOneDim() {}


         /**
             Optimized method to evaluate at the same time the function value and derivative at a point x.
              Often both value and derivatives are needed and it is often more efficient to compute them at the same time.
             Derived class should implement this method if performances play an important role and if it is faster to
             evaluate value and derivative at the same time

         */
         virtual void FdF(double x, double &f, double &df) const
         {
            f = operator()(x);
            df = Derivative(x);
         }



      };



   } // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_IFunction */
