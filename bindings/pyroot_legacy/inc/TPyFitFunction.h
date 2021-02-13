// Author: Wim Lavrijsen   November 2010

#ifndef ROOT_TPyFitFunction
#define ROOT_TPyFitFunction

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyFitFunction                                                           //
//                                                                          //
// Python base class to work with Math::IMultiGenFunction                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


//- ROOT
#include "Math/IFunction.h"
#include "Rtypes.h"

// Python
struct _object;
typedef _object PyObject;


class TPyMultiGenFunction : public ROOT::Math::IMultiGenFunction {
public:
// ctor/dtor, and assignment
   TPyMultiGenFunction( PyObject* self = 0 );
   virtual ~TPyMultiGenFunction();

// Math::IMultiGenFunction implementation
   virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const
      { return new TPyMultiGenFunction( fPySelf ); }
   virtual unsigned int NDim() const;
   virtual double DoEval( const double* x ) const;

   ClassDef( TPyMultiGenFunction, 1 );   //Python for Math::IMultiGenFunction equivalent

private:
// to prevent confusion when handing 'self' from python
   TPyMultiGenFunction( const TPyMultiGenFunction& src ) : ROOT::Math::IMultiGenFunction( src ) {}
   TPyMultiGenFunction& operator=( const TPyMultiGenFunction& ) { return *this; }

private:
   PyObject* fPySelf;              //! actual python object
};


class TPyMultiGradFunction : public ROOT::Math::IMultiGradFunction {
public:
// ctor/dtor, and assignment
   TPyMultiGradFunction( PyObject* self = 0 );
   virtual ~TPyMultiGradFunction();

// Math::IMultiGenFunction implementation
   virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const
      { return new TPyMultiGradFunction( fPySelf ); }
   virtual unsigned int NDim() const;
   virtual double DoEval( const double* x ) const;

   virtual void Gradient( const double* x, double* grad ) const;
   virtual void FdF( const double* x, double& f, double* df ) const;
   virtual double DoDerivative( const double * x, unsigned int icoord ) const;

   ClassDef( TPyMultiGradFunction, 1 );   //Python for Math::IMultiGradFunction equivalent

private:
// to prevent confusion when handing 'self' from python
   TPyMultiGradFunction( const TPyMultiGradFunction& src ) :
       ROOT::Math::IMultiGenFunction( src ), ROOT::Math::IMultiGradFunction( src ) {}
   TPyMultiGradFunction& operator=( const TPyMultiGradFunction& ) { return *this; }

private:
   PyObject* fPySelf;              //! actual python object
};

#endif
