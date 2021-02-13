// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_LASymMatrix
#define ROOT_Minuit2_LASymMatrix

#include "Minuit2/MnConfig.h"
#include "Minuit2/ABSum.h"
#include "Minuit2/VectorOuterProduct.h"
#include "Minuit2/MatrixInverse.h"
#include "Minuit2/StackAllocator.h"

#include <cassert>
#include <memory>
#include <cstring> // for memcopy

// extern StackAllocator StackAllocatorHolder::Get();

namespace ROOT {

namespace Minuit2 {

int Mndaxpy(unsigned int, double, const double *, int, double *, int);
int Mndscal(unsigned int, double, double *, int);

class LAVector;

int Invert(LASymMatrix &);

/**
   Class describing a symmetric matrix of size n.
   The size is specified as a run-time argument passed in the
   constructor.
   The class uses expression templates for the operations and functions.
   Only the independent data are kept in the fdata array of size n*(n+1)/2
   containing the lower triangular data
 */

class LASymMatrix {

private:
   LASymMatrix() : fSize(0), fNRow(0), fData(0) {}

public:
   typedef sym Type;

   LASymMatrix(unsigned int n)
      : fSize(n * (n + 1) / 2), fNRow(n),
        fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * n * (n + 1) / 2))
   {
      //     assert(fSize>0);
      std::memset(fData, 0, fSize * sizeof(double));
      //     std::cout<<"LASymMatrix(unsigned int n), n= "<<n<<std::endl;
   }

   ~LASymMatrix()
   {
      //     std::cout<<"~LASymMatrix()"<<std::endl;
      //     if(fData) std::cout<<"deleting "<<fSize<<std::endl;
      //     else std::cout<<"no delete"<<std::endl;
      //     if(fData) delete [] fData;
      if (fData)
         StackAllocatorHolder::Get().Deallocate(fData);
   }

   LASymMatrix(const LASymMatrix &v)
      : fSize(v.size()), fNRow(v.Nrow()),
        fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * v.size()))
   {
      //     std::cout<<"LASymMatrix(const LASymMatrix& v)"<<std::endl;
      std::memcpy(fData, v.Data(), fSize * sizeof(double));
   }

   LASymMatrix &operator=(const LASymMatrix &v)
   {
      //     std::cout<<"LASymMatrix& operator=(const LASymMatrix& v)"<<std::endl;
      //     std::cout<<"fSize= "<<fSize<<std::endl;
      //     std::cout<<"v.size()= "<<v.size()<<std::endl;
      assert(fSize == v.size());
      std::memcpy(fData, v.Data(), fSize * sizeof(double));
      return *this;
   }

   template <class T>
   LASymMatrix(const ABObj<sym, LASymMatrix, T> &v)
      : fSize(v.Obj().size()), fNRow(v.Obj().Nrow()),
        fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * v.Obj().size()))
   {
      //     std::cout<<"LASymMatrix(const ABObj<sym, LASymMatrix, T>& v)"<<std::endl;
      // std::cout<<"allocate "<<fSize<<std::endl;
      std::memcpy(fData, v.Obj().Data(), fSize * sizeof(double));
      Mndscal(fSize, double(v.f()), fData, 1);
      // std::cout<<"fData= "<<fData[0]<<" "<<fData[1]<<std::endl;
   }

   template <class A, class B, class T>
   LASymMatrix(const ABObj<sym, ABSum<ABObj<sym, A, T>, ABObj<sym, B, T>>, T> &sum) : fSize(0), fNRow(0), fData(0)
   {
      //     std::cout<<"template<class A, class B, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym, A, T>,
      //     ABObj<sym, B, T> > >& sum)"<<std::endl; recursive construction
      (*this) = sum.Obj().A();
      (*this) += sum.Obj().B();
      // std::cout<<"leaving template<class A, class B, class T> LASymMatrix(const ABObj..."<<std::endl;
   }

   template <class A, class T>
   LASymMatrix(const ABObj<sym, ABSum<ABObj<sym, LASymMatrix, T>, ABObj<sym, A, T>>, T> &sum)
      : fSize(0), fNRow(0), fData(0)
   {
      //     std::cout<<"template<class A, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym, LASymMatrix, T>,
      //     ABObj<sym, A, T> >,T>& sum)"<<std::endl;

      // recursive construction
      // std::cout<<"(*this)=sum.Obj().B();"<<std::endl;
      (*this) = sum.Obj().B();
      // std::cout<<"(*this)+=sum.Obj().A();"<<std::endl;
      (*this) += sum.Obj().A();
      // std::cout<<"leaving template<class A, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym,
      // LASymMatrix,.."<<std::endl;
   }

   template <class A, class T>
   LASymMatrix(const ABObj<sym, ABObj<sym, A, T>, T> &something) : fSize(0), fNRow(0), fData(0)
   {
      //     std::cout<<"template<class A, class T> LASymMatrix(const ABObj<sym, ABObj<sym, A, T>, T>&
      //     something)"<<std::endl;
      (*this) = something.Obj();
      (*this) *= something.f();
      // std::cout<<"leaving template<class A, class T> LASymMatrix(const ABObj<sym, ABObj<sym, A, T>, T>&
      // something)"<<std::endl;
   }

   template <class T>
   LASymMatrix(const ABObj<sym, MatrixInverse<sym, ABObj<sym, LASymMatrix, T>, T>, T> &inv)
      : fSize(inv.Obj().Obj().Obj().size()), fNRow(inv.Obj().Obj().Obj().Nrow()),
        fData((double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * inv.Obj().Obj().Obj().size()))
   {
      std::memcpy(fData, inv.Obj().Obj().Obj().Data(), fSize * sizeof(double));
      Mndscal(fSize, double(inv.Obj().Obj().f()), fData, 1);
      Invert(*this);
      Mndscal(fSize, double(inv.f()), fData, 1);
   }

   template <class A, class T>
   LASymMatrix(
      const ABObj<sym, ABSum<ABObj<sym, MatrixInverse<sym, ABObj<sym, LASymMatrix, T>, T>, T>, ABObj<sym, A, T>>, T>
         &sum)
      : fSize(0), fNRow(0), fData(0)
   {
      //     std::cout<<"template<class A, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym, MatrixInverse<sym,
      //     ABObj<sym, LASymMatrix, T>, T>, T>, ABObj<sym, A, T> >, T>& sum)"<<std::endl;

      // recursive construction
      (*this) = sum.Obj().B();
      (*this) += sum.Obj().A();
      // std::cout<<"leaving template<class A, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym,
      // LASymMatrix,.."<<std::endl;
   }

   LASymMatrix(const ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, double>, double>, double> &);

   template <class A, class T>
   LASymMatrix(
      const ABObj<sym, ABSum<ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>, ABObj<sym, A, T>>, T> &sum)
      : fSize(0), fNRow(0), fData(0)
   {
      //     std::cout<<"template<class A, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym,
      //     VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> ABObj<sym, A, T> >,T>& sum)"<<std::endl;

      // recursive construction
      (*this) = sum.Obj().B();
      (*this) += sum.Obj().A();
      // std::cout<<"leaving template<class A, class T> LASymMatrix(const ABObj<sym, ABSum<ABObj<sym,
      // LASymMatrix,.."<<std::endl;
   }

   LASymMatrix &operator+=(const LASymMatrix &m)
   {
      //     std::cout<<"LASymMatrix& operator+=(const LASymMatrix& m)"<<std::endl;
      assert(fSize == m.size());
      Mndaxpy(fSize, 1., m.Data(), 1, fData, 1);
      return *this;
   }

   LASymMatrix &operator-=(const LASymMatrix &m)
   {
      //     std::cout<<"LASymMatrix& operator-=(const LASymMatrix& m)"<<std::endl;
      assert(fSize == m.size());
      Mndaxpy(fSize, -1., m.Data(), 1, fData, 1);
      return *this;
   }

   template <class T>
   LASymMatrix &operator+=(const ABObj<sym, LASymMatrix, T> &m)
   {
      //     std::cout<<"template<class T> LASymMatrix& operator+=(const ABObj<sym, LASymMatrix, T>& m)"<<std::endl;
      assert(fSize == m.Obj().size());
      if (m.Obj().Data() == fData) {
         Mndscal(fSize, 1. + double(m.f()), fData, 1);
      } else {
         Mndaxpy(fSize, double(m.f()), m.Obj().Data(), 1, fData, 1);
      }
      // std::cout<<"fData= "<<fData[0]<<" "<<fData[1]<<std::endl;
      return *this;
   }

   template <class A, class T>
   LASymMatrix &operator+=(const ABObj<sym, A, T> &m)
   {
      //     std::cout<<"template<class A, class T> LASymMatrix& operator+=(const ABObj<sym, A,T>& m)"<<std::endl;
      (*this) += LASymMatrix(m);
      return *this;
   }

   template <class T>
   LASymMatrix &operator+=(const ABObj<sym, MatrixInverse<sym, ABObj<sym, LASymMatrix, T>, T>, T> &m)
   {
      //     std::cout<<"template<class T> LASymMatrix& operator+=(const ABObj<sym, MatrixInverse<sym, ABObj<sym,
      //     LASymMatrix, T>, T>, T>& m)"<<std::endl;
      assert(fNRow > 0);
      LASymMatrix tmp(m.Obj().Obj());
      Invert(tmp);
      tmp *= double(m.f());
      (*this) += tmp;
      return *this;
   }

   template <class T>
   LASymMatrix &operator+=(const ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> &m)
   {
      //     std::cout<<"template<class T> LASymMatrix& operator+=(const ABObj<sym, VectorOuterProduct<ABObj<vec,
      //     LAVector, T>, T>, T>&"<<std::endl;
      assert(fNRow > 0);
      Outer_prod(*this, m.Obj().Obj().Obj(), m.f() * m.Obj().Obj().f() * m.Obj().Obj().f());
      return *this;
   }

   LASymMatrix &operator*=(double scal)
   {
      Mndscal(fSize, scal, fData, 1);
      return *this;
   }

   double operator()(unsigned int row, unsigned int col) const
   {
      assert(row < fNRow && col < fNRow);
      if (row > col)
         return fData[col + row * (row + 1) / 2];
      else
         return fData[row + col * (col + 1) / 2];
   }

   double &operator()(unsigned int row, unsigned int col)
   {
      assert(row < fNRow && col < fNRow);
      if (row > col)
         return fData[col + row * (row + 1) / 2];
      else
         return fData[row + col * (col + 1) / 2];
   }

   const double *Data() const { return fData; }

   double *Data() { return fData; }

   unsigned int size() const { return fSize; }

   unsigned int Nrow() const { return fNRow; }

   unsigned int Ncol() const { return Nrow(); }

private:
   unsigned int fSize;
   unsigned int fNRow;
   double *fData;

public:
   template <class T>
   LASymMatrix &operator=(const ABObj<sym, LASymMatrix, T> &v)
   {
      // std::cout<<"template<class T> LASymMatrix& operator=(ABObj<sym, LASymMatrix, T>& v)"<<std::endl;
      if (fSize == 0 && fData == 0) {
         fSize = v.Obj().size();
         fNRow = v.Obj().Nrow();
         fData = (double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * fSize);
      } else {
         assert(fSize == v.Obj().size());
      }
      // std::cout<<"fData= "<<fData[0]<<" "<<fData[1]<<std::endl;
      std::memcpy(fData, v.Obj().Data(), fSize * sizeof(double));
      (*this) *= v.f();
      return *this;
   }

   template <class A, class T>
   LASymMatrix &operator=(const ABObj<sym, ABObj<sym, A, T>, T> &something)
   {
      // std::cout<<"template<class A, class T> LASymMatrix& operator=(const ABObj<sym, ABObj<sym, A, T>, T>&
      // something)"<<std::endl;
      if (fSize == 0 && fData == 0) {
         (*this) = something.Obj();
         (*this) *= something.f();
      } else {
         LASymMatrix tmp(something.Obj());
         tmp *= something.f();
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      // std::cout<<"template<class A, class T> LASymMatrix& operator=(const ABObj<sym, ABObj<sym, A, T>, T>&
      // something)"<<std::endl;
      return *this;
   }

   template <class A, class B, class T>
   LASymMatrix &operator=(const ABObj<sym, ABSum<ABObj<sym, A, T>, ABObj<sym, B, T>>, T> &sum)
   {
      // std::cout<<"template<class A, class B, class T> LASymMatrix& operator=(const ABObj<sym, ABSum<ABObj<sym, A, T>,
      // ABObj<sym, B, T> >,T>& sum)"<<std::endl;
      // recursive construction
      if (fSize == 0 && fData == 0) {
         (*this) = sum.Obj().A();
         (*this) += sum.Obj().B();
         (*this) *= sum.f();
      } else {
         LASymMatrix tmp(sum.Obj().A());
         tmp += sum.Obj().B();
         tmp *= sum.f();
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      return *this;
   }

   template <class A, class T>
   LASymMatrix &operator=(const ABObj<sym, ABSum<ABObj<sym, LASymMatrix, T>, ABObj<sym, A, T>>, T> &sum)
   {
      // std::cout<<"template<class A, class T> LASymMatrix& operator=(const ABObj<sym, ABSum<ABObj<sym, LASymMatrix,
      // T>, ABObj<sym, A, T> >,T>& sum)"<<std::endl;

      if (fSize == 0 && fData == 0) {
         // std::cout<<"fSize == 0 && fData == 0"<<std::endl;
         (*this) = sum.Obj().B();
         (*this) += sum.Obj().A();
         (*this) *= sum.f();
      } else {
         // std::cout<<"creating tmp variable"<<std::endl;
         LASymMatrix tmp(sum.Obj().B());
         tmp += sum.Obj().A();
         tmp *= sum.f();
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      // std::cout<<"leaving LASymMatrix& operator=(const ABObj<sym, ABSum<ABObj<sym, LASymMatrix..."<<std::endl;
      return *this;
   }

   template <class T>
   LASymMatrix &operator=(const ABObj<sym, MatrixInverse<sym, ABObj<sym, LASymMatrix, T>, T>, T> &inv)
   {
      if (fSize == 0 && fData == 0) {
         fSize = inv.Obj().Obj().Obj().size();
         fNRow = inv.Obj().Obj().Obj().Nrow();
         fData = (double *)StackAllocatorHolder::Get().Allocate(sizeof(double) * fSize);
         std::memcpy(fData, inv.Obj().Obj().Obj().Data(), fSize * sizeof(double));
         (*this) *= inv.Obj().Obj().f();
         Invert(*this);
         (*this) *= inv.f();
      } else {
         LASymMatrix tmp(inv.Obj().Obj());
         Invert(tmp);
         tmp *= double(inv.f());
         assert(fSize == tmp.size());
         std::memcpy(fData, tmp.Data(), fSize * sizeof(double));
      }
      return *this;
   }

   LASymMatrix &operator=(const ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, double>, double>, double> &);
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_LASymMatrix
