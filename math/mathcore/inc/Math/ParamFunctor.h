// @(#)root/mathcore:$Id$
// Author: L. Moneta Mon Nov 13 15:58:13 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for Functor classes.
// design is inspired by the Loki Functor

#ifndef ROOT_Math_ParamFunctor
#define ROOT_Math_ParamFunctor

// #ifndef ROOT_Math_IFunction
// #include "Math/IFunction.h"
// #endif

// #ifndef Root_Math_StaticCheck
// #include "Math/StaticCheck.h"
// #endif

//#ifndef __CINT__
//#include <memory>

#include "RtypesCore.h"
#include <functional>
#include <iostream>

namespace ROOT {

namespace Math {


/** class defining the signature for multi-dim parametric functions

   @ingroup  ParamFunctor_int
 */
template<class T>
class ParamFunctionBase {
  public:
   virtual ~ParamFunctionBase() {}
   virtual T operator() (const T * x, const double *p) = 0;
   virtual T operator() (T * x, double *p) = 0;
   virtual ParamFunctionBase * Clone() const = 0;
};



/**
   ParamFunctor Handler class is responsible for wrapping any other functor and pointer to
   free C functions.
   It can be created from any function implementing the correct signature
   corresponding to the requested type

   @ingroup  ParamFunctor_int

*/
#ifndef __CINT__

template<class ParentFunctor, class Func >
class ParamFunctorHandler : public ParentFunctor::Impl {

   typedef typename ParentFunctor::EvalType EvalType;
   typedef typename ParentFunctor::Impl     Base;

public:

   // constructor
   ParamFunctorHandler(const Func & fun) : fFunc(fun) {}


   virtual ~ParamFunctorHandler() {}


   // for 1D functions
   inline EvalType operator() (EvalType x, double *p)  {
      return fFunc(x,p);
   }
//    inline double operator() (double x, const double *p) const {
//       return fFunc(x,p);
//    }
   // for multi-dimensional functions
//    inline double operator() (const double * x, const double *p) const {
//       return fFunc(x,p);
//    }
   inline EvalType operator() (EvalType * x, double *p)  {
      return FuncEvaluator<Func, EvalType>::Eval(fFunc,x,p);
   }

   inline EvalType operator() (const EvalType * x, const double *p)  {
      return FuncEvaluator<Func, EvalType>::EvalConst(fFunc,x,p);
   }

   // clone (use same pointer)
   ParamFunctorHandler  * Clone() const {
      return new ParamFunctorHandler(fFunc);
   }


private :

   Func fFunc;

   // structure to distinguish pointer types
   template <typename F,typename  T> struct FuncEvaluator {
      inline static T Eval( F & f, T *x, double * p) {
         return f(x, p);
      }

      inline static T EvalConst( F & f, const T *x, const double * p) {
         return f((T*)x, (double*)p);
      }
   };

   template <typename F, typename T> struct FuncEvaluator<F*, T> {
      inline static T Eval( F * f, T *x, double * p) {
         return (*f)(x, p);
      }

      inline static T EvalConst( F * f, const T *x, const double * p) {
         return (*f)((T*)x, (double*)p);

      }
   };

   template <typename F,typename  T> struct FuncEvaluator<F* const, T> {
      inline static T Eval( const F * f, T *x, double * p) {
         return (*f)(x, p);
      }

      inline static T EvalConst( const F * f, const T *x, const double * p) {
         return (*f)((T*)x, (double*)p);
      }
   };

   // need maybe also volatile ?
};


#if defined(__MAKECINT__) || defined(G__DICTIONARY)
// needed since CINT initialize it with TRootIOCtor
//class TRootIOCtor;
template<class ParentFunctor>
class ParamFunctorHandler<ParentFunctor,TRootIOCtor *> : public ParentFunctor::Impl
{
public:

   ParamFunctorHandler(TRootIOCtor  *) {}

   double operator() (double *, double * )  { return 0; }

   double operator() (const double *, const double * )  { return 0; }
   // clone (use same pointer)
   ParamFunctorHandler  * Clone() const {
      return 0;
   }

};
#endif


/**
   ParamFunctor Handler to Wrap pointers to member functions

   @ingroup  ParamFunctor_int
*/
template <class ParentFunctor, typename PointerToObj,
          typename PointerToMemFn>
class ParamMemFunHandler : public ParentFunctor::Impl
{
   typedef typename ParentFunctor::Impl Base;


public:

   /// constructor from a pointer to the class and a pointer to the function
   ParamMemFunHandler(const PointerToObj& pObj, PointerToMemFn pMemFn)
      : fObj(pObj), fMemFn(pMemFn)
   {}

   virtual ~ParamMemFunHandler() {}

//    inline double operator() (double x, const double * p) const {
//       return ((*fObj).*fMemFn)(x,p);
//    }

   inline double operator() (double x, double * p)  {
      return ((*fObj).*fMemFn)(x,p);
   }

//    inline double operator() (const double * x, const double * p) const {
//       return ((*fObj).*fMemFn)(x,p);
//    }

   inline double operator() (double * x, double * p)  {
      return MemFuncEvaluator<PointerToObj,PointerToMemFn, double>::Eval(fObj,fMemFn,x,p);
   }

   inline double operator() (const double * x, const double * p)  {
      return MemFuncEvaluator<PointerToObj,PointerToMemFn, double>::EvalConst(fObj,fMemFn,x,p);
   }

   // clone (use same pointer)
   ParamMemFunHandler  * Clone() const {
      return new ParamMemFunHandler(fObj, fMemFn);
   }

private:

   // structure to distinguish pointer types
   template <typename PObj, typename F,typename  T> struct MemFuncEvaluator {
      inline static T Eval(PObj & pobj, F &  f, T *x, double * p) {
         return ((*pobj).*f)(x, p);
      }

      inline static T EvalConst(PObj & pobj, F & f, const T *x, const double * p) {
         return ((*pobj).*f)((T*)x, (double*)p);
      }
   };


   // // these are needed ??
   // template <typename PObj, typename F, typename T> struct MemFuncEvaluator<PObj,F*, T> {
   //    inline static T Eval(PObj & pobj,  F * f, T *x, double * p) {
   //       return ((*pobj).*f)f(x, p);
   //    }

   //    inline static T EvalConst(PObj & pobj,  F * f, const T *x, const double * p) {
   //       return ((*pobj).*f)((T*)x, (double*)p);

   //    }
   // };

   // template <typename PObj, typename F,typename  T> struct FuncEvaluator<PObj,F* const, T> {
   //    inline static T Eval(PObj &, const F * f, T *x, double * p) {
   //       return ((*pobj).*f)f(x, p);
   //    }

   //    inline static T EvalConst(PObj & pobj, const F * f, const T *x, const double * p) {
   //       return ((*pobj).*f)((T*)x, (double*)p);
   //    }
   // };

private :
   ParamMemFunHandler(const ParamMemFunHandler&); // Not implemented
   ParamMemFunHandler& operator=(const ParamMemFunHandler&); // Not implemented

   PointerToObj fObj;
   PointerToMemFn fMemFn;

};

#endif



/**
   Param Functor class for Multidimensional functions.
   It is used to wrap in a very simple and convenient way
   any other C++ callable object (implemention double operator( const double *, const double * ) )
   or a member function with the correct signature,
   like Foo::EvalPar(const double *, const double *)

   @ingroup  ParamFunc

 */


template<class T>
class ParamFunctorTempl   {


public:

   typedef  T                    EvalType;
   typedef  ParamFunctionBase<T> Impl;


   /**
      Default constructor
   */
   ParamFunctorTempl ()  : fImpl(0) {}


   /**
       construct from a pointer to member function (multi-dim type)
    */
   template <class PtrObj, typename MemFn>
   ParamFunctorTempl(const PtrObj& p, MemFn memFn)
      : fImpl(new ParamMemFunHandler<ParamFunctorTempl<T>, PtrObj, MemFn>(p, memFn))
   {}



   /**
      construct from another generic Functor of multi-dimension
    */
   template <typename Func>
   explicit ParamFunctorTempl( const Func & f) :
      fImpl(new ParamFunctorHandler<ParamFunctorTempl<T>,Func>(f) )
   {}



   // specialization used in TF1
   typedef T (* FreeFunc ) (T * , double *);
   ParamFunctorTempl(FreeFunc f) :
      fImpl(new ParamFunctorHandler<ParamFunctorTempl<T>,FreeFunc>(f) )
   {
   }

   // specialization used in TF1
   ParamFunctorTempl(const std::function<T(const T *f, const Double_t *param)> &func) :
      fImpl(new ParamFunctorHandler<ParamFunctorTempl<T>, const std::function<T(const T *f, const Double_t *param)>>(func))
   {
   }

   /**
      Destructor (no operations)
   */
   virtual ~ParamFunctorTempl ()  {
      if (fImpl) delete fImpl;
   }

   /**
      Copy constructor
   */
   ParamFunctorTempl(const ParamFunctorTempl & rhs) :
      fImpl(0)
   {
//       if (rhs.fImpl.get() != 0)
//          fImpl = std::unique_ptr<Impl>( (rhs.fImpl)->Clone() );
      if (rhs.fImpl != 0)  fImpl = rhs.fImpl->Clone();
   }

   /**
      Assignment operator
   */
   ParamFunctorTempl & operator = (const ParamFunctorTempl & rhs)  {
//      ParamFunctor copy(rhs);
      // swap unique_ptr by hand
//       Impl * p = fImpl.release();
//       fImpl.reset(copy.fImpl.release());
//       copy.fImpl.reset(p);

      if(this != &rhs) {
         if (fImpl) delete fImpl;
         fImpl = 0;
         if (rhs.fImpl != 0)
            fImpl = rhs.fImpl->Clone();
      }
      return *this;
   }

   void * GetImpl() { return (void *) fImpl; }


   T operator() ( T * x, double * p)  {
      return (*fImpl)(x,p);
   }

   T operator() (const T * x, const double * p)  {
      return (*fImpl)(x,p);
   }


   bool Empty() const { return fImpl == 0; }


   void SetFunction(Impl * f) {
      fImpl = f;
   }

private :


   //std::unique_ptr<Impl> fImpl;
   Impl * fImpl;


};


using ParamFunctor = ParamFunctorTempl<double>;

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_ParamFunctor */
