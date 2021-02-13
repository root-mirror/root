// @(#)root/mathmore:$Id$
// Author: L. Moneta Wed Dec  6 11:45:55 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class MultiDimParamFunctionAdapter

#ifndef ROOT_Math_MultiDimParamFunctionAdapter
#define ROOT_Math_MultiDimParamFunctionAdapter

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"

#include "Math/WrappedFunction.h"

namespace ROOT {

   namespace Math {


      /**
         MultiDimParamFunctionAdapter class to wrap a one-dimensional parametric function in
         a multi dimensional parameteric function interface
         This is used typically in fitting where internally the function is stored as multidimension

         To wrap a non-parametric one-dim function in a multi-dim interface one can use simply a
         ROOT::Math::WrappedFunction<ROOT::Math::IGenFunction> or ROOT::Math::Functor
         and ROOT::Math::GradFunctor for gradient functions

         This class differs from WrappedParamFunction in the fact that the parameters are not stored in
         the adapter class and optionally it keeps a cloned and managed copy of the adapter class.

         @ingroup  ParamFunc

      */
      class MultiDimParamFunctionAdapter : public IParametricFunctionMultiDimTempl<double>  {

      public:

         typedef IParamMultiFunction::BaseFunc BaseFunc;


         /**
            Constructor from a parametric one dim function interface from a const reference
            Own the function in this case
         */
         MultiDimParamFunctionAdapter(const IParamFunction &f) :
            fOwn(true)
         {
            fFunc = dynamic_cast<IParamFunction *>(f.Clone());
         }

         /**
            Constructor from a parametric one dim function interface from a non-const reference
            Do not own the function in this case
         */
         MultiDimParamFunctionAdapter(IParamFunction &f) :
            fOwn(false),
            fFunc(&f)
         { }


         /**
            Copy constructor. Different behaviour according if function is owned or not
          */
         MultiDimParamFunctionAdapter(const MultiDimParamFunctionAdapter &rhs) :
            BaseFunc(),
            IParamMultiFunction(),
            fOwn(rhs.fOwn),
            fFunc(0)
         {
            if (fOwn)
               fFunc = dynamic_cast<IParamFunction *>((rhs.fFunc)->Clone());
         }

         /**
            Destructor (no operations)
         */
         virtual ~MultiDimParamFunctionAdapter()
         {
            if (fOwn && fFunc != 0) delete fFunc;
         }


         /**
            Assignment operator
          */
         MultiDimParamFunctionAdapter &operator=(const MultiDimParamFunctionAdapter &rhs)
         {
            fOwn = rhs.fOwn;
            if (fOwn) {
               if (fFunc) delete fFunc; // delete previously existing copy
               fFunc = dynamic_cast<IParamFunction *>((rhs.fFunc)->Clone());
            } else
               fFunc = rhs.fFunc;

            return *this;
         }

         /**
            clone
         */
         virtual BaseFunc *Clone() const
         {
            return new MultiDimParamFunctionAdapter(*this);
         }

      public:

         // methods required by interface
         const double *Parameters() const
         {
            return  fFunc->Parameters();
         }

         void SetParameters(const double *p)
         {
            fFunc->SetParameters(p);
         }

         unsigned int NPar() const
         {
            return fFunc->NPar();
         }

         unsigned int NDim() const
         {
            return 1;
         }


      private:

         /// needed by the interface
         double DoEvalPar(const double *x, const double *p) const
         {
            return (*fFunc)(*x, p);
         }


      private:

         bool fOwn;
         IParamFunction *fFunc;

      };



      /**
         MultiDimParamGradFunctionAdapter class to wrap a one-dimensional parametric gradient function in
         a multi dimensional parameteric gradient function interface
         This is used typically in fitting where internally the function is stored as multidimension

         To wrap a non-parametric one-dim gradient function in a multi-dim interface one can use simply a
           a ROOT::Math::GradFunctor

         The parameters are not stored in the adapter class and by default the pointer to the 1D function is owned.
         This means that deleteing the class deletes also the 1D function and copying the class copies also the
         1D function
         This class differs from WrappedParamFunction in the fact that the parameters are not stored in
         the adapter class and optionally it keeps a cloned and managed copy of the adapter class.

         @ingroup  ParamFunc

      */
      class MultiDimParamGradFunctionAdapter : public IParamMultiGradFunction  {

      public:

         typedef IParamMultiGradFunction::BaseFunc BaseFunc;


         /**
            Constructor from a param one dim function interface from a const reference
            Copy and manage the own function pointer
         */
         MultiDimParamGradFunctionAdapter(const IParamGradFunction &f) :
            fOwn(true)
         {
            fFunc = dynamic_cast<IParamGradFunction *>(f.Clone());
         }

         /**
            Constructor from a param one dim function interface from a non const reference
            Do not  own the function pointer in this case
         */
         MultiDimParamGradFunctionAdapter(IParamGradFunction &f) :
            fOwn(false),
            fFunc(&f)
         { }


         /**
            Copy constructor. Different behaviour according if function is owned or not
          */
         MultiDimParamGradFunctionAdapter(const MultiDimParamGradFunctionAdapter &rhs) :
            BaseFunc(),
            IParamMultiGradFunction(),
            fOwn(rhs.fOwn),
            fFunc(rhs.fFunc)
         {
            if (fOwn)
               fFunc = dynamic_cast<IParamGradFunction *>((rhs.fFunc)->Clone());
         }

         /**
            Destructor (no operations)
         */
         virtual ~MultiDimParamGradFunctionAdapter()
         {
            if (fOwn && fFunc != 0) delete fFunc;
         }


         /**
            Assignment operator
          */
         MultiDimParamGradFunctionAdapter &operator=(const MultiDimParamGradFunctionAdapter &rhs)
         {
            fOwn = rhs.fOwn;
            if (fOwn) {
               if (fFunc) delete fFunc; // delete previously existing copy
               fFunc = dynamic_cast<IParamGradFunction *>((rhs.fFunc)->Clone());
            } else
               fFunc = rhs.fFunc;

            return *this;
         }

         /**
            clone
         */
         virtual BaseFunc *Clone() const
         {
            return new MultiDimParamGradFunctionAdapter(*this);
         }

      public:

         // methods required by interface
         const double *Parameters() const
         {
            return  fFunc->Parameters();
         }

         void SetParameters(const double *p)
         {
            fFunc->SetParameters(p);
         }

         unsigned int NPar() const
         {
            return fFunc->NPar();
         }

         unsigned int NDim() const
         {
            return 1;
         }

//    void Gradient(const double *x, double * grad) const {
//       grad[0] = fFunc->Derivative( *x);
//    }

         void ParameterGradient(const double *x, const double *p, double *grad) const
         {
            fFunc->ParameterGradient(*x, p, grad);
         }

         //  using IParamMultiGradFunction::BaseFunc::operator();

      private:

         /// functions needed by interface
         double DoEvalPar(const double *x, const double *p) const
         {
            return (*fFunc)(*x, p);
         }

//    double DoDerivative(const double * x, unsigned int ) const {
//       return fFunc->Derivative(*x);
//    }

         double DoParameterDerivative(const double *x, const double *p, unsigned int ipar) const
         {
            return fFunc->ParameterDerivative(*x, p, ipar);
         }

      private:

         bool fOwn;
         IParamGradFunction *fFunc;

      };




   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_MultiDimParamFunctionAdapter */
