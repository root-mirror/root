// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_LAVector
#define ROOT_Minuit2_LAVector

#include "Minuit2/ABSum.h"
#include "Minuit2/ABProd.h"
#include "Minuit2/LASymMatrix.h"

#include <cassert>
#include <memory>

#include "Minuit2/StackAllocator.h"

namespace ROOT {

namespace Minuit2 {

// extern StackAllocator StackAllocatorHolder::Get();

int Mndaxpy(unsigned int, double, const double *, int, double *, int);
int Mndscal(unsigned int, double, double *, int);
int Mndspmv(const char *, unsigned int, double, const double *, const double *, int, double, double *, int);

class LAVector {

private:
   LAVector() : fSize(0), fData(0) {}

public:
   typedef vec Type;

   //   LAVector() : fSize(0), fData(0) {}

   LAVector(unsigned int n) : fSize(n), fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * n))
   {
      //     assert(fSize>0);
      std::memset(fData, 0, size() * sizeof(double));
      //     std::cout<<"LAVector(unsigned int n), n= "<<n<<std::endl;
   }

   ~LAVector()
   {
      //     std::cout<<"~LAVector()"<<std::endl;
      //    if(fData) std::cout<<"deleting "<<fSize<<std::endl;
      //     else std::cout<<"no delete"<<std::endl;
      //     if(fData) delete [] fData;
      if (fData)
         StackAllocatorHolder::Get().Deallocate(fData);
   }

   LAVector(const LAVector &v)
      : fSize(v.size()), fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * v.size()))
   {
      //     std::cout<<"LAVector(const LAVector& v)"<<std::endl;
      std::memcpy(fData, v.Data(), fSize * sizeof(double));
   }

   LAVector &operator=(const LAVector &v)
   {
      //     std::cout<<"LAVector& operator=(const LAVector& v)"<<std::endl;
      //     std::cout<<"fSize= "<<fSize<<std::endl;
      //     std::cout<<"v.size()= "<<v.size()<<std::endl;
      assert(fSize == v.size());
      std::memcpy(fData, v.Data(), fSize * sizeof(double));
      return *this;
   }

   template <class T>
   LAVector(const ABObj<vec, LAVector, T> &v)
      : fSize(v.Obj().size()), fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * v.Obj().size()))
   {
      //     std::cout<<"LAVector(const ABObj<LAVector, T>& v)"<<std::endl;
      //     std::cout<<"allocate "<<fSize<<std::endl;
      std::memcpy(fData, v.Obj().Data(), fSize * sizeof(T));
      (*this) *= T(v.f());
      //     std::cout<<"fData= "<<fData[0]<<" "<<fData[1]<<std::endl;
   }

   template <class A, class B, class T>
   LAVector(const ABObj<vec, ABSum<ABObj<vec, A, T>, ABObj<vec, B, T>>, T> &sum) : fSize(0), fData(0)
   {
      //     std::cout<<"template<class A, class B, class T> LAVector(const ABObj<ABSum<ABObj<A, T>, ABObj<B, T> > >&
      //     sum)"<<std::endl;
      (*this) = sum.Obj().A();
      (*this) += sum.Obj().B();
      (*this) *= double(sum.f());
   }

   template <class A, class T>
   LAVector(const ABObj<vec, ABSum<ABObj<vec, LAVector, T>, ABObj<vec, A, T>>, T> &sum) : fSize(0), fData(0)
   {
      //     std::cout<<"template<class A, class T> LAVector(const ABObj<ABSum<ABObj<LAVector, T>, ABObj<A, T> >,T>&
      //     sum)"<<std::endl;

      // recursive construction
      //     std::cout<<"(*this)=sum.Obj().B();"<<std::endl;
      (*this) = sum.Obj().B();
      //     std::cout<<"(*this)+=sum.Obj().A();"<<std::endl;
      (*this) += sum.Obj().A();
      (*this) *= double(sum.f());
      //     std::cout<<"leaving template<class A, class T> LAVector(const ABObj<ABSum<ABObj<LAVector,.."<<std::endl;
   }

   template <class A, class T>
   LAVector(const ABObj<vec, ABObj<vec, A, T>, T> &something) : fSize(0), fData(0)
   {
      //     std::cout<<"template<class A, class T> LAVector(const ABObj<ABObj<A, T>, T>& something)"<<std::endl;
      (*this) = something.Obj();
      (*this) *= something.f();
   }

   //
   template <class T>
   LAVector(const ABObj<vec, ABProd<ABObj<sym, LASymMatrix, T>, ABObj<vec, LAVector, T>>, T> &prod)
      : fSize(prod.Obj().B().Obj().size()),
        fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * prod.Obj().B().Obj().size()))
   {
      //     std::cout<<"template<class T> LAVector(const ABObj<vec, ABProd<ABObj<sym, LASymMatrix, T>, ABObj<vec,
      //     LAVector, T> >, T>& prod)"<<std::endl;

      Mndspmv("U", fSize, prod.f() * prod.Obj().A().f() * prod.Obj().B().f(), prod.Obj().A().Obj().Data(),
              prod.Obj().B().Obj().Data(), 1, 0., fData, 1);
   }

   //
   template <class T>
   LAVector(const ABObj<
            vec,
            ABSum<ABObj<vec, ABProd<ABObj<sym, LASymMatrix, T>, ABObj<vec, LAVector, T>>, T>, ABObj<vec, LAVector, T>>,
            T> &prod)
      : fSize(0), fData(0)
   {
      (*this) = prod.Obj().B();
      (*this) += prod.Obj().A();
      (*this) *= double(prod.f());
   }

   //
   LAVector &operator+=(const LAVector &m)
   {
      //     std::cout<<"LAVector& operator+=(const LAVector& m)"<<std::endl;
      assert(fSize == m.size());
      Mndaxpy(fSize, 1., m.Data(), 1, fData, 1);
      return *this;
   }

   LAVector &operator-=(const LAVector &m)
   {
      //     std::cout<<"LAVector& operator-=(const LAVector& m)"<<std::endl;
      assert(fSize == m.size());
      Mndaxpy(fSize, -1., m.Data(), 1, fData, 1);
      return *this;
   }

   template <class T>
   LAVector &operator+=(const ABObj<vec, LAVector, T> &m)
   {
      //     std::cout<<"template<class T> LAVector& operator+=(const ABObj<LAVector, T>& m)"<<std::endl;
      assert(fSize == m.Obj().size());
      if (m.Obj().Data() == fData) {
         Mndscal(fSize, 1. + double(m.f()), fData, 1);
      } else {
         Mndaxpy(fSize, double(m.f()), m.Obj().Data(), 1, fData, 1);
      }
      //     std::cout<<"fData= "<<fData[0]<<" "<<fData[1]<<std::endl;
      return *this;
   }

   template <class A, class T>
   LAVector &operator+=(const ABObj<vec, A, T> &m)
   {
      //     std::cout<<"template<class A, class T> LAVector& operator+=(const ABObj<A,T>& m)"<<std::endl;
      (*this) += LAVector(m);
      return *this;
   }

   template <class T>
   LAVector &operator+=(const ABObj<vec, ABProd<ABObj<sym, LASymMatrix, T>, ABObj<vec, LAVector, T>>, T> &prod)
   {
      Mndspmv("U", fSize, prod.f() * prod.Obj().A().f() * prod.Obj().B().f(), prod.Obj().A().Obj().Data(),
              prod.Obj().B().Data(), 1, 1., fData, 1);
      return *this;
   }

   LAVector &operator*=(double scal)
   {
      Mndscal(fSize, scal, fData, 1);
      return *this;
   }

   double operator()(unsigned int i) const
   {
      assert(i < fSize);
      return fData[i];
   }

   double &operator()(unsigned int i)
   {
      assert(i < fSize);
      return fData[i];
   }

   double operator[](unsigned int i) const
   {
      assert(i < fSize);
      return fData[i];
   }

   double &operator[](unsigned int i)
   {
      assert(i < fSize);
      return fData[i];
   }

   const double *Data() const { return fData; }

   double *Data() { return fData; }

   unsigned int size() const { return fSize; }

private:
   unsigned int fSize;
   double *fData;

public:
   template <class T>
   LAVector &operator=(const ABObj<vec, LAVector, T> &v)
   {
      //     std::cout<<"template<class T> LAVector& operator=(ABObj<LAVector, T>& v)"<<std::endl;
      if (fSize == 0 && fData == 0) {
         fSize = v.Obj().size();
         fData = (double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * fSize);
      } else {
         assert(fSize == v.Obj().size());
      }
      std::memcpy(fData, v.Obj().Data(), fSize * sizeof(double));
      (*this) *= T(v.f());
      return *this;
   }

   template <class A, class T>
   LAVector &operator=(const ABObj<vec, ABObj<vec, A, T>, T> &something)
   {
      //     std::cout<<"template<class A, class T> LAVector& operator=(const ABObj<ABObj<A, T>, T>&
      //     something)"<<std::endl;
      if (fSize == 0 && fData == 0) {
         (*this) = something.Obj();
      } else {
         LAVector tmp(something.Obj());
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      (*this) *= something.f();
      return *this;
   }

   template <class A, class B, class T>
   LAVector &operator=(const ABObj<vec, ABSum<ABObj<vec, A, T>, ABObj<vec, B, T>>, T> &sum)
   {
      if (fSize == 0 && fData == 0) {
         (*this) = sum.Obj().A();
         (*this) += sum.Obj().B();
      } else {
         LAVector tmp(sum.Obj().A());
         tmp += sum.Obj().B();
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      (*this) *= sum.f();
      return *this;
   }

   template <class A, class T>
   LAVector &operator=(const ABObj<vec, ABSum<ABObj<vec, LAVector, T>, ABObj<vec, A, T>>, T> &sum)
   {
      if (fSize == 0 && fData == 0) {
         (*this) = sum.Obj().B();
         (*this) += sum.Obj().A();
      } else {
         LAVector tmp(sum.Obj().A());
         tmp += sum.Obj().B();
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      (*this) *= sum.f();
      return *this;
   }

   //
   template <class T>
   LAVector &operator=(const ABObj<vec, ABProd<ABObj<sym, LASymMatrix, T>, ABObj<vec, LAVector, T>>, T> &prod)
   {
      if (fSize == 0 && fData == 0) {
         fSize = prod.Obj().B().Obj().size();
         fData = (double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * fSize);
         Mndspmv("U", fSize, double(prod.f() * prod.Obj().A().f() * prod.Obj().B().f()), prod.Obj().A().Obj().Data(),
                 prod.Obj().B().Obj().Data(), 1, 0., fData, 1);
      } else {
         LAVector tmp(prod.Obj().B());
         assert(fSize == tmp.size());
         Mndspmv("U", fSize, double(prod.f() * prod.Obj().A().f()), prod.Obj().A().Obj().Data(), tmp.Data(), 1, 0.,
                 fData, 1);
      }
      return *this;
   }

   //
   template <class T>
   LAVector &
   operator=(const ABObj<
             vec,
             ABSum<ABObj<vec, ABProd<ABObj<sym, LASymMatrix, T>, ABObj<vec, LAVector, T>>, T>, ABObj<vec, LAVector, T>>,
             T> &prod)
   {
      if (fSize == 0 && fData == 0) {
         (*this) = prod.Obj().B();
         (*this) += prod.Obj().A();
      } else {
         //       std::cout<<"creating tmp variable"<<std::endl;
         LAVector tmp(prod.Obj().B());
         tmp += prod.Obj().A();
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      (*this) *= prod.f();
      return *this;
   }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_LAVector
