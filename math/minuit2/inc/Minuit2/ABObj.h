// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ABObj
#define ROOT_Minuit2_ABObj

#include "Minuit2/ABTypes.h"

namespace ROOT {

namespace Minuit2 {

template <class mtype, class M, class T>
class ABObj {

public:
   typedef mtype Type;

private:
   ABObj() : fObject(M()), fFactor(T(0.)) {}

   ABObj &operator=(const ABObj &) { return *this; }

   template <class a, class b, class c>
   ABObj(const ABObj<a, b, c> &) : fObject(M()), fFactor(T(0.))
   {
   }

   template <class a, class b, class c>
   ABObj &operator=(const ABObj<a, b, c> &)
   {
      return *this;
   }

public:
   ABObj(const M &obj) : fObject(obj), fFactor(T(1.)) {}

   ABObj(const M &obj, T factor) : fObject(obj), fFactor(factor) {}

   ~ABObj() {}

   ABObj(const ABObj &obj) : fObject(obj.fObject), fFactor(obj.fFactor) {}

   template <class b, class c>
   ABObj(const ABObj<mtype, b, c> &obj) : fObject(M(obj.Obj())), fFactor(T(obj.f()))
   {
   }

   const M &Obj() const { return fObject; }

   T f() const { return fFactor; }

private:
   M fObject;
   T fFactor;
};

class LAVector;
template <>
class ABObj<vec, LAVector, double> {

public:
   typedef vec Type;

private:
   ABObj &operator=(const ABObj &) = delete;

public:
   ABObj(const LAVector &obj) : fObject(obj), fFactor(double(1.)) {}

   ABObj(const LAVector &obj, double factor) : fObject(obj), fFactor(factor) {}

   ~ABObj() {}

   // remove copy constructure to Fix a problem in AIX
   // should be able to use the compiler generated one
   //   ABObj(const ABObj& obj) :
   //     fObject(obj.fObject), fFactor(obj.fFactor) {}

   template <class c>
   ABObj(const ABObj<vec, LAVector, c> &obj) : fObject(obj.fObject), fFactor(double(obj.fFactor))
   {
   }

   const LAVector &Obj() const { return fObject; }

   double f() const { return fFactor; }

private:
   const LAVector &fObject;
   double fFactor;
};

class LASymMatrix;
template <>
class ABObj<sym, LASymMatrix, double> {

public:
   typedef sym Type;

private:
   ABObj &operator=(const ABObj &) { return *this; }

public:
   ABObj(const LASymMatrix &obj) : fObject(obj), fFactor(double(1.)) {}

   ABObj(const LASymMatrix &obj, double factor) : fObject(obj), fFactor(factor) {}

   ~ABObj() {}

   ABObj(const ABObj &obj) : fObject(obj.fObject), fFactor(obj.fFactor) {}

   template <class c>
   ABObj(const ABObj<vec, LASymMatrix, c> &obj) : fObject(obj.fObject), fFactor(double(obj.fFactor))
   {
   }

   const LASymMatrix &Obj() const { return fObject; }

   double f() const { return fFactor; }

private:
   const LASymMatrix &fObject;
   double fFactor;
};

// templated scaling operator *
template <class mt, class M, class T>
inline ABObj<mt, M, T> operator*(T f, const M &obj)
{
   return ABObj<mt, M, T>(obj, f);
}

// templated operator /
template <class mt, class M, class T>
inline ABObj<mt, M, T> operator/(const M &obj, T f)
{
   return ABObj<mt, M, T>(obj, T(1.) / f);
}

// templated unary operator -
template <class mt, class M, class T>
inline ABObj<mt, M, T> operator-(const M &obj)
{
   return ABObj<mt, M, T>(obj, T(-1.));
}

/*
// specialization for LAVector

inline ABObj<vec, LAVector, double> operator*(double f, const LAVector& obj) {
  return ABObj<vec, LAVector, double>(obj, f);
}

inline ABObj<vec, LAVector, double> operator/(const LAVector& obj, double f) {
  return ABObj<vec, LAVector, double>(obj, double(1.)/f);
}

inline ABObj<vec,LAVector,double> operator-(const LAVector& obj) {
  return ABObj<vec,LAVector,double>(obj, double(-1.));
}
*/

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ABObj
