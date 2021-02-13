// @(#)root/matrix:$Id$
// Authors: Fons Rademakers, Eddy Offermann  Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVectorT
    \ingroup Matrix

 TVectorT

 Template class of Vectors in the linear algebra package.

 See \ref MatrixPage for the documentation of the linear algebra package

 Unless otherwise specified, vector indices always start with 0,
 spanning up to the specified limit-1.

 For (n) vectors where n <= kSizeMax (5 currently) storage space is
 available on the stack, thus avoiding expensive allocation/
 deallocation of heap space . However, this introduces of course
 kSizeMax overhead for each vector object . If this is an issue
 recompile with a new appropriate value (>=0) for kSizeMax

 Another way to assign and store vector data is through Use
 see for instance stressLinear.cxx file .

 Note that Constructors/assignments exists for all different matrix
 views

 For usage examples see $ROOTSYS/test/stressLinear.cxx

*/

#include "TVectorT.h"
#include "TBuffer.h"
#include "TMath.h"
#include "TROOT.h"
#include "Varargs.h"

templateClassImp(TVectorT);


////////////////////////////////////////////////////////////////////////////////
/// Delete data pointer m, if it was assigned on the heap

template<class Element>
void TVectorT<Element>::Delete_m(Int_t size,Element *&m)
{
   if (m) {
      if (size > kSizeMax)
         delete [] m;
      m = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return data pointer . if requested size <= kSizeMax, assign pointer
/// to the stack space

template<class Element>
Element* TVectorT<Element>::New_m(Int_t size)
{
   if (size == 0) return 0;
   else {
      if ( size <= kSizeMax )
         return fDataStack;
      else {
         Element *heap = new Element[size];
         return heap;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add vector v to this vector

template<class Element>
void TVectorT<Element>::Add(const TVectorT<Element> &v)
{
   if (gMatrixCheck && !AreCompatible(*this,v)) {
      Error("Add(TVectorT<Element> &)","vector's not compatible");
      return;
   }

   const Element *sp = v.GetMatrixArray();
         Element *tp = this->GetMatrixArray();
   const Element * const tp_last = tp+fNrows;
   while (tp < tp_last)
      *tp++ += *sp++;
}

////////////////////////////////////////////////////////////////////////////////
/// Set this vector to v1+v2

template<class Element>
void TVectorT<Element>::Add(const TVectorT<Element> &v1,const TVectorT<Element> &v2)
{
   if (gMatrixCheck) {
      if ( !AreCompatible(*this,v1) && !AreCompatible(*this,v2)) {
         Error("Add(TVectorT<Element> &)","vectors not compatible");
         return;
      }
   }

   const Element *sv1 = v1.GetMatrixArray();
   const Element *sv2 = v2.GetMatrixArray();
         Element *tp = this->GetMatrixArray();
   const Element * const tp_last = tp+fNrows;
   while (tp < tp_last)
      *tp++ = *sv1++ + *sv2++;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy copySize doubles from *oldp to *newp . However take care of the
/// situation where both pointers are assigned to the same stack space

template<class Element>
Int_t TVectorT<Element>::Memcpy_m(Element *newp,const Element *oldp,Int_t copySize,
                                  Int_t newSize,Int_t oldSize)
{
   if (copySize == 0 || oldp == newp) return 0;
   else {
      if ( newSize <= kSizeMax && oldSize <= kSizeMax ) {
         // both pointers are inside fDataStack, be careful with copy direction !
         if (newp > oldp) {
            for (Int_t i = copySize-1; i >= 0; i--)
               newp[i] = oldp[i];
         } else {
            for (Int_t i = 0; i < copySize; i++)
               newp[i] = oldp[i];
         }
      }
      else {
         memcpy(newp,oldp,copySize*sizeof(Element));
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate new vector. Arguments are number of rows and row
/// lowerbound (0 default).

template<class Element>
void TVectorT<Element>::Allocate(Int_t nrows,Int_t row_lwb,Int_t init)
{
   fIsOwner  = kTRUE;
   fNrows    = 0;
   fRowLwb   = 0;
   fElements = 0;

   if (nrows < 0)
   {
      Error("Allocate","nrows=%d",nrows);
      return;
   }

   MakeValid();
   fNrows   = nrows;
   fRowLwb  = row_lwb;

   fElements = New_m(fNrows);
   if (init)
      memset(fElements,0,fNrows*sizeof(Element));
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor n-vector

template<class Element>
TVectorT<Element>::TVectorT(Int_t n)
{
   Allocate(n,0,1);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor [lwb..upb]-vector

template<class Element>
TVectorT<Element>::TVectorT(Int_t lwb,Int_t upb)
{
   Allocate(upb-lwb+1,lwb,1);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor n-vector with data copied from array elements

template<class Element>
TVectorT<Element>::TVectorT(Int_t n,const Element *elements)
{
   Allocate(n,0);
   SetElements(elements);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor [lwb..upb]-vector with data copied from array elements

template<class Element>
TVectorT<Element>::TVectorT(Int_t lwb,Int_t upb,const Element *elements)
{
   Allocate(upb-lwb+1,lwb);
   SetElements(elements);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

template<class Element>
TVectorT<Element>::TVectorT(const TVectorT &another) : TObject(another)
{
   R__ASSERT(another.IsValid());
   Allocate(another.GetUpb()-another.GetLwb()+1,another.GetLwb());
   *this = another;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor : create vector from matrix row

template<class Element>
TVectorT<Element>::TVectorT(const TMatrixTRow_const<Element> &mr) : TObject()
{
   const TMatrixTBase<Element> *mt = mr.GetMatrix();
   R__ASSERT(mt->IsValid());
   Allocate(mt->GetColUpb()-mt->GetColLwb()+1,mt->GetColLwb());
   *this = mr;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor : create vector from matrix column

template<class Element>
TVectorT<Element>::TVectorT(const TMatrixTColumn_const<Element> &mc) : TObject()
{
   const TMatrixTBase<Element> *mt = mc.GetMatrix();
   R__ASSERT(mt->IsValid());
   Allocate(mt->GetRowUpb()-mt->GetRowLwb()+1,mt->GetRowLwb());
   *this = mc;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor : create vector from matrix diagonal

template<class Element>
TVectorT<Element>::TVectorT(const TMatrixTDiag_const<Element> &md) : TObject()
{
   const TMatrixTBase<Element> *mt = md.GetMatrix();
   R__ASSERT(mt->IsValid());
   Allocate(TMath::Min(mt->GetNrows(),mt->GetNcols()));
   *this = md;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a vector and assign initial values. Argument list should contain
/// Element values to assign to vector elements. The list must be
/// terminated by the string "END". Example:
/// TVectorT foo(1,3,0.0,1.0,1.5,"END");

template<class Element>
TVectorT<Element>::TVectorT(Int_t lwb,Int_t upb,Double_t iv1, ...)
{
   if (upb < lwb) {
      Error("TVectorT(Int_t, Int_t, ...)","upb(%d) < lwb(%d)",upb,lwb);
      return;
   }

   Allocate(upb-lwb+1,lwb);

   va_list args;
   va_start(args,iv1);             // Init 'args' to the beginning of
                                        // the variable length list of args

   (*this)(lwb) = iv1;
   for (Int_t i = lwb+1; i <= upb; i++)
      (*this)(i) = (Element)va_arg(args,Double_t);

   if (strcmp((char *)va_arg(args,char *),"END"))
      Error("TVectorT(Int_t, Int_t, ...)","argument list must be terminated by \"END\"");

   va_end(args);
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the vector to [lwb:upb] .
/// New dynamic elemenst are created, the overlapping part of the old ones are
/// copied to the new structures, then the old elements are deleted.

template<class Element>
TVectorT<Element> &TVectorT<Element>::ResizeTo(Int_t lwb,Int_t upb)
{
   R__ASSERT(IsValid());
   if (!fIsOwner) {
      Error("ResizeTo(lwb,upb)","Not owner of data array,cannot resize");
      return *this;
   }

   const Int_t new_nrows = upb-lwb+1;

   if (fNrows > 0) {

      if (fNrows == new_nrows && fRowLwb == lwb)
         return *this;
      else if (new_nrows == 0) {
         Clear();
         return *this;
      }

      Element    *elements_old = GetMatrixArray();
      const Int_t  nrows_old    = fNrows;
      const Int_t  rowLwb_old   = fRowLwb;

      Allocate(new_nrows,lwb);
      R__ASSERT(IsValid());
      if (fNrows > kSizeMax || nrows_old > kSizeMax)
         memset(GetMatrixArray(),0,fNrows*sizeof(Element));
      else if (fNrows > nrows_old)
         memset(GetMatrixArray()+nrows_old,0,(fNrows-nrows_old)*sizeof(Element));

    // Copy overlap
      const Int_t rowLwb_copy = TMath::Max(fRowLwb,rowLwb_old);
      const Int_t rowUpb_copy = TMath::Min(fRowLwb+fNrows-1,rowLwb_old+nrows_old-1);
      const Int_t nrows_copy  = rowUpb_copy-rowLwb_copy+1;

      const Int_t nelems_new = fNrows;
      Element *elements_new = GetMatrixArray();
      if (nrows_copy > 0) {
         const Int_t rowOldOff = rowLwb_copy-rowLwb_old;
         const Int_t rowNewOff = rowLwb_copy-fRowLwb;
         Memcpy_m(elements_new+rowNewOff,elements_old+rowOldOff,nrows_copy,nelems_new,nrows_old);
      }

      Delete_m(nrows_old,elements_old);
   } else {
      Allocate(upb-lwb+1,lwb,1);
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Use the array data to fill the vector lwb..upb]

template<class Element>
TVectorT<Element> &TVectorT<Element>::Use(Int_t lwb,Int_t upb,Element *data)
{
   if (upb < lwb) {
      Error("Use","upb(%d) < lwb(%d)",upb,lwb);
      return *this;
   }

   Clear();
   fNrows    = upb-lwb+1;
   fRowLwb   = lwb;
   fElements = data;
   fIsOwner  = kFALSE;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Get subvector [row_lwb..row_upb]; The indexing range of the
/// returned vector depends on the argument option:
///
/// option == "S" : return [0..row_upb-row_lwb+1] (default)
/// else          : return [row_lwb..row_upb]

template<class Element>
TVectorT<Element> &TVectorT<Element>::GetSub(Int_t row_lwb,Int_t row_upb,TVectorT<Element> &target,Option_t *option) const
{
   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      if (row_lwb < fRowLwb || row_lwb > fRowLwb+fNrows-1) {
         Error("GetSub","row_lwb out of bounds");
         return target;
      }
      if (row_upb < fRowLwb || row_upb > fRowLwb+fNrows-1) {
         Error("GetSub","row_upb out of bounds");
         return target;
      }
      if (row_upb < row_lwb) {
         Error("GetSub","row_upb < row_lwb");
         return target;
      }
   }

   TString opt(option);
   opt.ToUpper();
   const Int_t shift = (opt.Contains("S")) ? 1 : 0;

   Int_t row_lwb_sub;
   Int_t row_upb_sub;
   if (shift) {
      row_lwb_sub = 0;
      row_upb_sub = row_upb-row_lwb;
   } else {
      row_lwb_sub = row_lwb;
      row_upb_sub = row_upb;
   }

   target.ResizeTo(row_lwb_sub,row_upb_sub);
   const Int_t nrows_sub = row_upb_sub-row_lwb_sub+1;

   const Element *ap = this->GetMatrixArray()+(row_lwb-fRowLwb);
         Element *bp = target.GetMatrixArray();

   for (Int_t irow = 0; irow < nrows_sub; irow++)
      *bp++ = *ap++;

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert vector source starting at [row_lwb], thereby overwriting the part
/// [row_lwb..row_lwb+nrows_source];

template<class Element>
TVectorT<Element> &TVectorT<Element>::SetSub(Int_t row_lwb,const TVectorT<Element> &source)
{
   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(source.IsValid());

      if (row_lwb < fRowLwb && row_lwb > fRowLwb+fNrows-1) {
         Error("SetSub","row_lwb outof bounds");
         return *this;
      }
      if (row_lwb+source.GetNrows() > fRowLwb+fNrows) {
         Error("SetSub","source vector too large");
         return *this;
      }
   }

   const Int_t nRows_source = source.GetNrows();

   const Element *bp = source.GetMatrixArray();
         Element *ap = this->GetMatrixArray()+(row_lwb-fRowLwb);

   for (Int_t irow = 0; irow < nRows_source; irow++)
      *ap++ = *bp++;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector elements to zero

template<class Element>
TVectorT<Element> &TVectorT<Element>::Zero()
{
   R__ASSERT(IsValid());
   memset(this->GetMatrixArray(),0,fNrows*sizeof(Element));
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Take an absolute value of a vector, i.e. apply Abs() to each element.

template<class Element>
TVectorT<Element> &TVectorT<Element>::Abs()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      *ep = TMath::Abs(*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Square each element of the vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::Sqr()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      *ep = (*ep) * (*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Take square root of all elements.

template<class Element>
TVectorT<Element> &TVectorT<Element>::Sqrt()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      R__ASSERT(*ep >= 0);
      if (*ep >= 0)
         *ep = TMath::Sqrt(*ep);
      else
         Error("Sqrt()","v(%ld) = %g < 0",Long_t(ep-this->GetMatrixArray()),(float)*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// v[i] = 1/v[i]

template<class Element>
TVectorT<Element> &TVectorT<Element>::Invert()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      R__ASSERT(*ep != 0.0);
      if (*ep != 0.0)
         *ep = 1./ *ep;
      else
         Error("Invert()","v(%ld) = %g",Long_t(ep-this->GetMatrixArray()),(float)*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Keep only element as selected through array select non-zero

template<class Element>
TVectorT<Element> &TVectorT<Element>::SelectNonZeros(const TVectorT<Element> &select)
{
   if (gMatrixCheck && !AreCompatible(*this,select)) {
      Error("SelectNonZeros(const TVectorT<Element> &","vector's not compatible");
      return *this;
   }

   const Element *sp = select.GetMatrixArray();
         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      if (*sp == 0.0)
         *ep = 0.0;
      sp++; ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the 1-norm of the vector SUM{ |v[i]| }.

template<class Element>
Element TVectorT<Element>::Norm1() const
{
   R__ASSERT(IsValid());

   Element norm = 0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      norm += TMath::Abs(*ep++);

   return norm;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the square of the 2-norm SUM{ v[i]^2 }.

template<class Element>
Element TVectorT<Element>::Norm2Sqr() const
{
   R__ASSERT(IsValid());

   Element norm = 0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      norm += (*ep) * (*ep);
      ep++;
   }

   return norm;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the infinity-norm of the vector MAX{ |v[i]| }.

template<class Element>
Element TVectorT<Element>::NormInf() const
{
   R__ASSERT(IsValid());

   Element norm = 0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      norm = TMath::Max(norm,TMath::Abs(*ep++));

   return norm;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the number of elements != 0.0

template<class Element>
Int_t TVectorT<Element>::NonZeros() const
{
   R__ASSERT(IsValid());

   Int_t nr_nonzeros = 0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (*ep++) nr_nonzeros++;

   return nr_nonzeros;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute sum of elements

template<class Element>
Element TVectorT<Element>::Sum() const
{
   R__ASSERT(IsValid());

   Element sum = 0.0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      sum += *ep++;

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// return minimum vector element value

template<class Element>
Element TVectorT<Element>::Min() const
{
   R__ASSERT(IsValid());

   const Int_t index = TMath::LocMin(fNrows,fElements);
   return fElements[index];
}

////////////////////////////////////////////////////////////////////////////////
/// return maximum vector element value

template<class Element>
Element TVectorT<Element>::Max() const
{
   R__ASSERT(IsValid());

   const Int_t index = TMath::LocMax(fNrows,fElements);
   return fElements[index];
}

////////////////////////////////////////////////////////////////////////////////
/// Notice that this assignment does NOT change the ownership :
/// if the storage space was adopted, source is copied to
/// this space .

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(const TVectorT<Element> &source)
{
   if (gMatrixCheck && !AreCompatible(*this,source)) {
      Error("operator=(const TVectorT<Element> &)","vectors not compatible");
      return *this;
   }

   if (this->GetMatrixArray() != source.GetMatrixArray()) {
      TObject::operator=(source);
      memcpy(fElements,source.GetMatrixArray(),fNrows*sizeof(Element));
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a matrix row to a vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(const TMatrixTRow_const<Element> &mr)
{
   const TMatrixTBase<Element> *mt = mr.GetMatrix();

   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(mt->IsValid());
      if (mt->GetColLwb() != fRowLwb || mt->GetNcols() != fNrows) {
         Error("operator=(const TMatrixTRow_const &)","vector and row not compatible");
         return *this;
      }
   }

   const Int_t inc   = mr.GetInc();
   const Element *rp = mr.GetPtr();              // Row ptr
         Element *ep = this->GetMatrixArray();   // Vector ptr
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      *ep++ = *rp;
       rp += inc;
   }

   R__ASSERT(rp == mr.GetPtr()+mt->GetNcols());

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a matrix column to a vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(const TMatrixTColumn_const<Element> &mc)
{
   const TMatrixTBase<Element> *mt = mc.GetMatrix();

   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(mt->IsValid());
      if (mt->GetRowLwb() != fRowLwb || mt->GetNrows() != fNrows) {
         Error("operator=(const TMatrixTColumn_const &)","vector and column not compatible");
         return *this;
      }
   }

   const Int_t inc    = mc.GetInc();
   const Element *cp = mc.GetPtr();              // Column ptr
         Element *ep = this->GetMatrixArray();   // Vector ptr
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      *ep++ = *cp;
       cp += inc;
   }

   R__ASSERT(cp == mc.GetPtr()+mt->GetNoElements());

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign the matrix diagonal to a vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(const TMatrixTDiag_const<Element> &md)
{
   const TMatrixTBase<Element> *mt = md.GetMatrix();

   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(mt->IsValid());
      if (md.GetNdiags() != fNrows) {
         Error("operator=(const TMatrixTDiag_const &)","vector and matrix-diagonal not compatible");
        return *this;
      }
   }

   const Int_t    inc = md.GetInc();
   const Element *dp = md.GetPtr();              // Diag ptr
         Element *ep = this->GetMatrixArray();   // Vector ptr
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      *ep++ = *dp;
       dp += inc;
   }

   R__ASSERT(dp < md.GetPtr()+mt->GetNoElements()+inc);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a sparse matrix row to a vector. The matrix row is implicitly transposed
/// to allow the assignment in the strict sense.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(const TMatrixTSparseRow_const<Element> &mr)
{
   const TMatrixTBase<Element> *mt = mr.GetMatrix();

   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(mt->IsValid());
      if (mt->GetColLwb() != fRowLwb || mt->GetNcols() != fNrows) {
         Error("operator=(const TMatrixTSparseRow_const &)","vector and row not compatible");
         return *this;
      }
   }

   const Int_t nIndex = mr.GetNindex();
   const Element * const prData = mr.GetDataPtr();          // Row Data ptr
   const Int_t    * const prCol  = mr.GetColPtr();           // Col ptr
         Element * const pvData = this->GetMatrixArray();   // Vector ptr

   memset(pvData,0,fNrows*sizeof(Element));
   for (Int_t index = 0; index < nIndex; index++) {
      const Int_t icol = prCol[index];
      pvData[icol] = prData[index];
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a sparse matrix diagonal to a vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(const TMatrixTSparseDiag_const<Element> &md)
{
  const TMatrixTBase<Element> *mt = md.GetMatrix();

   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(mt->IsValid());
      if (md.GetNdiags() != fNrows) {
         Error("operator=(const TMatrixTSparseDiag_const &)","vector and matrix-diagonal not compatible");
         return *this;
      }
   }

   Element * const pvData = this->GetMatrixArray();
   for (Int_t idiag = 0; idiag < fNrows; idiag++)
      pvData[idiag] = md(idiag);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign val to every element of the vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator=(Element val)
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      *ep++ = val;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Add val to every element of the vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator+=(Element val)
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      *ep++ += val;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Subtract val from every element of the vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator-=(Element val)
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      *ep++ -= val;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply every element of the vector with val.

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator*=(Element val)
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      *ep++ *= val;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Add vector source

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator+=(const TVectorT<Element> &source)
{
   if (gMatrixCheck && !AreCompatible(*this,source)) {
      Error("operator+=(const TVectorT<Element> &)","vector's not compatible");
      return *this;
   }

   const Element *sp = source.GetMatrixArray();
         Element *tp = this->GetMatrixArray();
   const Element * const tp_last = tp+fNrows;
   while (tp < tp_last)
      *tp++ += *sp++;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Subtract vector source

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator-=(const TVectorT<Element> &source)
{
   if (gMatrixCheck && !AreCompatible(*this,source)) {
      Error("operator-=(const TVectorT<Element> &)","vector's not compatible");
      return *this;
   }

   const Element *sp = source.GetMatrixArray();
         Element *tp = this->GetMatrixArray();
   const Element * const tp_last = tp+fNrows;
   while (tp < tp_last)
      *tp++ -= *sp++;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// "Inplace" multiplication target = A*target. A needn't be a square one
/// If target has to be resized, it should own the storage: fIsOwner = kTRUE

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator*=(const TMatrixT<Element> &a)
{
   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(a.IsValid());
      if (a.GetNcols() != fNrows || a.GetColLwb() != fRowLwb) {
         Error("operator*=(const TMatrixT &)","vector and matrix incompatible");
         return *this;
      }
   }

   const Bool_t doResize = (fNrows != a.GetNrows() || fRowLwb != a.GetRowLwb());
   if (doResize && !fIsOwner) {
      Error("operator*=(const TMatrixT &)","vector has to be resized but not owner");
      return *this;
   }

   Element work[kWorkMax];
   Bool_t isAllocated = kFALSE;
   Element *elements_old = work;
   const Int_t nrows_old = fNrows;
   if (nrows_old > kWorkMax) {
      isAllocated = kTRUE;
      elements_old = new Element[nrows_old];
   }
   memcpy(elements_old,fElements,nrows_old*sizeof(Element));

   if (doResize) {
      const Int_t rowlwb_new = a.GetRowLwb();
      const Int_t nrows_new  = a.GetNrows();
      ResizeTo(rowlwb_new,rowlwb_new+nrows_new-1);
   }
   memset(fElements,0,fNrows*sizeof(Element));

   const Element *mp = a.GetMatrixArray();     // Matrix row ptr
         Element *tp = this->GetMatrixArray(); // Target vector ptr
#ifdef CBLAS
   if (typeid(Element) == typeid(Double_t))
      cblas_dgemv(CblasRowMajor,CblasNoTrans,a.GetNrows(),a.GetNcols(),1.0,mp,
                  a.GetNcols(),elements_old,1,0.0,tp,1);
   else if (typeid(Element) != typeid(Float_t))
      cblas_sgemv(CblasRowMajor,CblasNoTrans,a.GetNrows(),a.GetNcols(),1.0,mp,
                  a.GetNcols(),elements_old,1,0.0,tp,1);
   else
      Error("operator*=","type %s not implemented in BLAS library",typeid(Element));
#else
   const Element * const tp_last = tp+fNrows;
   while (tp < tp_last) {
      Element sum = 0;
      for (const Element *sp = elements_old; sp < elements_old+nrows_old; )
         sum += *sp++ * *mp++;
      *tp++ = sum;
   }
   R__ASSERT(mp == a.GetMatrixArray()+a.GetNoElements());
#endif

   if (isAllocated)
      delete [] elements_old;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// "Inplace" multiplication target = A*target. A needn't be a square one
/// If target has to be resized, it should own the storage: fIsOwner = kTRUE

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator*=(const TMatrixTSparse<Element> &a)
{
   if (gMatrixCheck) {
      R__ASSERT(IsValid());
      R__ASSERT(a.IsValid());
      if (a.GetNcols() != fNrows || a.GetColLwb() != fRowLwb) {
         Error("operator*=(const TMatrixTSparse &)","vector and matrix incompatible");
         return *this;
      }
   }

   const Bool_t doResize = (fNrows != a.GetNrows() || fRowLwb != a.GetRowLwb());
   if (doResize && !fIsOwner) {
      Error("operator*=(const TMatrixTSparse &)","vector has to be resized but not owner");
      return *this;
   }

   Element work[kWorkMax];
   Bool_t isAllocated = kFALSE;
   Element *elements_old = work;
   const Int_t nrows_old = fNrows;
   if (nrows_old > kWorkMax) {
      isAllocated = kTRUE;
      elements_old = new Element[nrows_old];
   }
   memcpy(elements_old,fElements,nrows_old*sizeof(Element));

   if (doResize) {
      const Int_t rowlwb_new = a.GetRowLwb();
      const Int_t nrows_new  = a.GetNrows();
      ResizeTo(rowlwb_new,rowlwb_new+nrows_new-1);
   }
   memset(fElements,0,fNrows*sizeof(Element));

   const Int_t   * const pRowIndex = a.GetRowIndexArray();
   const Int_t   * const pColIndex = a.GetColIndexArray();
   const Element * const mp        = a.GetMatrixArray();     // Matrix row ptr

   const Element * const sp = elements_old;
         Element *       tp = this->GetMatrixArray(); // Target vector ptr

   for (Int_t irow = 0; irow < fNrows; irow++) {
      const Int_t sIndex = pRowIndex[irow];
      const Int_t eIndex = pRowIndex[irow+1];
      Element sum = 0.0;
      for (Int_t index = sIndex; index < eIndex; index++) {
         const Int_t icol = pColIndex[index];
         sum += mp[index]*sp[icol];
      }
      tp[irow] = sum;
   }

   if (isAllocated)
      delete [] elements_old;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// "Inplace" multiplication target = A*target. A is symmetric .
/// vector size will not change

template<class Element>
TVectorT<Element> &TVectorT<Element>::operator*=(const TMatrixTSym<Element> &a)
{
  if (gMatrixCheck) {
     R__ASSERT(IsValid());
     R__ASSERT(a.IsValid());
     if (a.GetNcols() != fNrows || a.GetColLwb() != fRowLwb) {
        Error("operator*=(const TMatrixTSym &)","vector and matrix incompatible");
        return *this;
     }
  }

   Element work[kWorkMax];
   Bool_t isAllocated = kFALSE;
   Element *elements_old = work;
   const Int_t nrows_old = fNrows;
   if (nrows_old > kWorkMax) {
      isAllocated = kTRUE;
      elements_old = new Element[nrows_old];
   }
   memcpy(elements_old,fElements,nrows_old*sizeof(Element));
   memset(fElements,0,fNrows*sizeof(Element));

   const Element *mp = a.GetMatrixArray();     // Matrix row ptr
         Element *tp = this->GetMatrixArray(); // Target vector ptr
#ifdef CBLAS
   if (typeid(Element) == typeid(Double_t))
      cblas_dsymv(CblasRowMajor,CblasUpper,fNrows,1.0,mp,
                  fNrows,elements_old,1,0.0,tp,1);
   else if (typeid(Element) != typeid(Float_t))
      cblas_ssymv(CblasRowMajor,CblasUpper,fNrows,1.0,mp,
                  fNrows,elements_old,1,0.0,tp,1);
   else
      Error("operator*=","type %s not implemented in BLAS library",typeid(Element));
#else
   const Element * const tp_last = tp+fNrows;
   while (tp < tp_last) {
      Element sum = 0;
      for (const Element *sp = elements_old; sp < elements_old+nrows_old; )
         sum += *sp++ * *mp++;
      *tp++ = sum;
   }
   R__ASSERT(mp == a.GetMatrixArray()+a.GetNoElements());
#endif

   if (isAllocated)
      delete [] elements_old;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all vector elements equal to val?

template<class Element>
Bool_t TVectorT<Element>::operator==(Element val) const
{
   R__ASSERT(IsValid());

   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (!(*ep++ == val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all vector elements not equal to val?

template<class Element>
Bool_t TVectorT<Element>::operator!=(Element val) const
{
   R__ASSERT(IsValid());

   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (!(*ep++ != val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all vector elements < val?

template<class Element>
Bool_t TVectorT<Element>::operator<(Element val) const
{
   R__ASSERT(IsValid());

   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (!(*ep++ < val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all vector elements <= val?

template<class Element>
Bool_t TVectorT<Element>::operator<=(Element val) const
{
   R__ASSERT(IsValid());

   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (!(*ep++ <= val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all vector elements > val?

template<class Element>
Bool_t TVectorT<Element>::operator>(Element val) const
{
   R__ASSERT(IsValid());

   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (!(*ep++ > val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all vector elements >= val?

template<class Element>
Bool_t TVectorT<Element>::operator>=(Element val) const
{
   R__ASSERT(IsValid());

   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      if (!(*ep++ >= val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if vector elements as selected through array select are non-zero

template<class Element>
Bool_t TVectorT<Element>::MatchesNonZeroPattern(const TVectorT<Element> &select)
{
   if (gMatrixCheck && !AreCompatible(*this,select)) {
      Error("MatchesNonZeroPattern(const TVectorT&)","vector's not compatible");
      return kFALSE;
   }

   const Element *sp = select.GetMatrixArray();
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      if (*sp == 0.0 && *ep != 0.0)
         return kFALSE;
      sp++; ep++;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if vector elements as selected through array select are all positive

template<class Element>
Bool_t TVectorT<Element>::SomePositive(const TVectorT<Element> &select)
{
   if (gMatrixCheck && !AreCompatible(*this,select)) {
      Error("SomePositive(const TVectorT&)","vector's not compatible");
      return kFALSE;
   }

   const Element *sp = select.GetMatrixArray();
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      if (*sp != 0.0 && *ep <= 0.0)
         return kFALSE;
      sp++; ep++;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add to vector elements as selected through array select the value val

template<class Element>
void TVectorT<Element>::AddSomeConstant(Element val,const TVectorT<Element> &select)
{
   if (gMatrixCheck && !AreCompatible(*this,select))
      Error("AddSomeConstant(Element,const TVectorT&)(const TVectorT&)","vector's not compatible");

   const Element *sp = select.GetMatrixArray();
         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp) {
      if (*sp)
         *ep += val;
      sp++; ep++;
   }
}

extern Double_t Drand(Double_t &ix);

////////////////////////////////////////////////////////////////////////////////
/// randomize vector elements value

template<class Element>
void TVectorT<Element>::Randomize(Element alpha,Element beta,Double_t &seed)
{
   R__ASSERT(IsValid());

   const Element scale = beta-alpha;
   const Element shift = alpha/scale;

         Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNrows;
   while (ep < fp)
      *ep++ = scale*(Drand(seed)+shift);
}

////////////////////////////////////////////////////////////////////////////////
/// Apply action to each element of the vector.

template<class Element>
TVectorT<Element> &TVectorT<Element>::Apply(const TElementActionT<Element> &action)
{
   R__ASSERT(IsValid());
   for (Element *ep = fElements; ep < fElements+fNrows; ep++)
      action.Operation(*ep);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply action to each element of the vector. In action the location
/// of the current element is known.

template<class Element>
TVectorT<Element> &TVectorT<Element>::Apply(const TElementPosActionT<Element> &action)
{
   R__ASSERT(IsValid());

   Element *ep = fElements;
   for (action.fI = fRowLwb; action.fI < fRowLwb+fNrows; action.fI++)
      action.Operation(*ep++);

   R__ASSERT(ep == fElements+fNrows);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this vector
/// The histogram is named "TVectorT" by default and no title

template<class Element>
void TVectorT<Element>::Draw(Option_t *option)
{
   gROOT->ProcessLine(Form("THistPainter::PaintSpecialObjects((TObject*)0x%lx,\"%s\");",
                           (ULong_t)this, option));
}

////////////////////////////////////////////////////////////////////////////////
/// Print the vector as a list of elements.

template<class Element>
void TVectorT<Element>::Print(Option_t *flag) const
{
  if (!IsValid()) {
      Error("Print","Vector is invalid");
      return;
   }

   printf("\nVector (%d) %s is as follows",fNrows,flag);

   printf("\n\n     |   %6d  |", 1);
   printf("\n%s\n", "------------------");
   for (Int_t i = 0; i < fNrows; i++) {
      printf("%4d |",i+fRowLwb);
      //printf("%11.4g \n",(*this)(i+fRowLwb));
      printf("%g \n",(*this)(i+fRowLwb));
   }
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Check to see if two vectors are identical.

template<class Element>
Bool_t operator==(const TVectorT<Element> &v1,const TVectorT<Element> &v2)
{
   if (!AreCompatible(v1,v2)) return kFALSE;
   return (memcmp(v1.GetMatrixArray(),v2.GetMatrixArray(),v1.GetNrows()*sizeof(Element)) == 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the scalar product.

template<class Element>
Element operator*(const TVectorT<Element> &v1,const TVectorT<Element> &v2)
{
   if (gMatrixCheck) {
      if (!AreCompatible(v1,v2)) {
         Error("operator*(const TVectorT<Element> &,const TVectorT<Element> &)","vector's are incompatible");
         return 0.0;
      }
   }

   return Dot(v1,v2);
}

////////////////////////////////////////////////////////////////////////////////
/// Return source1+source2

template<class Element>
TVectorT<Element> operator+(const TVectorT<Element> &source1,const TVectorT<Element> &source2)
{
   TVectorT<Element> target = source1;
   target += source2;
   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Return source1-source2

template<class Element>
TVectorT<Element> operator-(const TVectorT<Element> &source1,const TVectorT<Element> &source2)
{
   TVectorT<Element> target = source1;
   target -= source2;
   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// return A * source

template<class Element>
TVectorT<Element> operator*(const TMatrixT<Element> &a,const TVectorT<Element> &source)
{
   R__ASSERT(a.IsValid());
   TVectorT<Element> target(a.GetRowLwb(),a.GetRowUpb());
   return Add(target,Element(1.0),a,source);
}

////////////////////////////////////////////////////////////////////////////////
/// return A * source

template<class Element>
TVectorT<Element> operator*(const TMatrixTSym<Element> &a,const TVectorT<Element> &source)
{
   R__ASSERT(a.IsValid());
   TVectorT<Element> target(a.GetRowLwb(),a.GetRowUpb());
   return Add(target,Element(1.0),a,source);
}

////////////////////////////////////////////////////////////////////////////////
/// return A * source

template<class Element>
TVectorT<Element> operator*(const TMatrixTSparse<Element> &a,const TVectorT<Element> &source)
{
   R__ASSERT(a.IsValid());
   TVectorT<Element> target(a.GetRowLwb(),a.GetRowUpb());
   return Add(target,Element(1.0),a,source);
}

////////////////////////////////////////////////////////////////////////////////
/// return val * source

template<class Element>
TVectorT<Element> operator*(Element val,const TVectorT<Element> &source)
{
   TVectorT<Element> target = source;
   target *= val;
   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// return inner-produvt v1 . v2

template<class Element>
Element Dot(const TVectorT<Element> &v1,const TVectorT<Element> &v2)
{
   const Element *v1p = v1.GetMatrixArray();
   const Element *v2p = v2.GetMatrixArray();
   Element sum = 0.0;
   const Element * const fv1p = v1p+v1.GetNrows();
   while (v1p < fv1p)
      sum += *v1p++ * *v2p++;

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the matrix M = v1 * v2'

template <class Element1, class Element2>
TMatrixT<Element1>
OuterProduct(const TVectorT<Element1> &v1,const TVectorT<Element2> &v2)
{
   // TMatrixD::GetSub does:
   //   TMatrixT tmp;
   // Doesn't compile here, because we are outside the class?
   // So we'll be explicit:
   TMatrixT<Element1> target;

   return OuterProduct(target,v1,v2);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the matrix M = v1 * v2'

template <class Element1,class Element2,class Element3>
TMatrixT<Element1>
&OuterProduct(TMatrixT<Element1> &target,const TVectorT<Element2> &v1,const TVectorT<Element3> &v2)
{
   target.ResizeTo(v1.GetLwb(), v1.GetUpb(), v2.GetLwb(), v2.GetUpb());

         Element1 *       mp      = target.GetMatrixArray();
   const Element1 * const m_last  = mp + target.GetNoElements();

   const Element2 *       v1p     = v1.GetMatrixArray();
   const Element2 * const v1_last = v1p + v1.GetNrows();

   const Element3 * const v20     = v2.GetMatrixArray();
   const Element3 *       v2p     = v20;
   const Element3 * const v2_last = v2p + v2.GetNrows();

   while (v1p < v1_last) {
      v2p = v20;
      while (v2p < v2_last) {
         *mp++ = *v1p * *v2p++ ;
      }
      v1p++;
  }

  R__ASSERT(v1p == v1_last && mp == m_last && v2p == v2_last);

  return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform v1 * M * v2, a scalar result

template <class Element1, class Element2, class Element3>
Element1 Mult(const TVectorT<Element1> &v1,const TMatrixT<Element2> &m,
              const TVectorT<Element3> &v2)
{
   if (gMatrixCheck) {
      if (!AreCompatible(v1, m)) {
         ::Error("Mult", "Vector v1 and matrix m incompatible");
         return 0;
      }
      if (!AreCompatible(m, v2)) {
         ::Error("Mult", "Matrix m and vector v2 incompatible");
         return 0;
      }
   }

   const Element1 *       v1p     = v1.GetMatrixArray();    // first of v1
   const Element1 * const v1_last = v1p + v1.GetNrows();    // last of  v1

   const Element2 *       mp      = m.GetMatrixArray();     // first of m
   const Element2 * const m_last  = mp + m.GetNoElements(); // last of  m

   const Element3 * const v20     = v2.GetMatrixArray();    // first of v2
   const Element3 *       v2p     = v20;                    // running  v2
   const Element3 * const v2_last = v2p + v2.GetNrows();    // last of  v2

   Element1 sum     = 0;      // scalar result accumulator
   Element3 dot     = 0;      // M_row * v2 dot product accumulator

   while (v1p < v1_last) {
      v2p  = v20;               // at beginning of v2
      while (v2p < v2_last) {   // compute (M[i] * v2) dot product
         dot += *mp++ * *v2p++;
      }
      sum += *v1p++ * dot;      // v1[i] * (M[i] * v2)
      dot  = 0;                 // start next dot product
   }

   R__ASSERT(v1p == v1_last && mp == m_last && v2p == v2_last);

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += scalar * source.

template<class Element>
TVectorT<Element> &Add(TVectorT<Element> &target,Element scalar,const TVectorT<Element> &source)
{
   if (gMatrixCheck && !AreCompatible(target,source)) {
      Error("Add(TVectorT<Element> &,Element,const TVectorT<Element> &)","vector's are incompatible");
      return target;
   }

   const Element *       sp  = source.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();
   if (scalar == 1.0 ) {
      while ( tp < ftp )
         *tp++ += *sp++;
   } else if (scalar == -1.0) {
      while ( tp < ftp )
         *tp++ -= *sp++;
   } else {
      while ( tp < ftp )
         *tp++ += scalar * *sp++;
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += scalar * A * source.
/// NOTE: in case scalar=0, do  target = A * source.

template<class Element>
TVectorT<Element> &Add(TVectorT<Element> &target,Element scalar,
                       const TMatrixT<Element> &a,const TVectorT<Element> &source)
{
   if (gMatrixCheck) {
      R__ASSERT(target.IsValid());
      R__ASSERT(a.IsValid());
      if (a.GetNrows() != target.GetNrows() || a.GetRowLwb() != target.GetLwb()) {
         Error("Add(TVectorT &,const TMatrixT &,const TVectorT &)","target vector and matrix are incompatible");
         return target;
      }

      R__ASSERT(source.IsValid());
      if (a.GetNcols() != source.GetNrows() || a.GetColLwb() != source.GetLwb()) {
         Error("Add(TVectorT &,const TMatrixT &,const TVectorT &)","source vector and matrix are incompatible");
         return target;
      }
   }

   const Element * const sp = source.GetMatrixArray();  // sources vector ptr
   const Element *       mp = a.GetMatrixArray();       // Matrix row ptr
         Element *       tp = target.GetMatrixArray();  // Target vector ptr
#ifdef CBLAS
   if (typeid(Element) == typeid(Double_t))
      cblas_dsymv(CblasRowMajor,CblasUpper,fNrows,scalar,mp,
                  fNrows,sp,1,0.0,tp,1);
   else if (typeid(Element) != typeid(Float_t))
      cblas_ssymv(CblasRowMajor,CblasUpper,fNrows,scalar,mp,
                  fNrows,sp,1,0.0,tp,1);
   else
      Error("operator*=","type %s not implemented in BLAS library",typeid(Element));
#else
   const Element * const sp_last = sp+source.GetNrows();
   const Element * const tp_last = tp+target.GetNrows();
   if (scalar == 1.0) {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++ += sum;
      }
   } else if (scalar == 0.0) {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++  = sum;
      }
   } else if (scalar == -1.0) {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++ -= sum;
      }
   } else {
      while (tp < tp_last) {
        const Element *sp1 = sp;
        Element sum = 0;
        while (sp1 < sp_last)
           sum += *sp1++ * *mp++;
        *tp++ += scalar * sum;
      }
   }

   if (gMatrixCheck) R__ASSERT(mp == a.GetMatrixArray()+a.GetNoElements());
#endif

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += A * source.
/// NOTE: in case scalar=0, do  target = A * source.

template<class Element>
TVectorT<Element> &Add(TVectorT<Element> &target,Element scalar,
                       const TMatrixTSym<Element> &a,const TVectorT<Element> &source)
{
   if (gMatrixCheck) {
      R__ASSERT(target.IsValid());
      R__ASSERT(source.IsValid());
      R__ASSERT(a.IsValid());
      if (a.GetNrows() != target.GetNrows() || a.GetRowLwb() != target.GetLwb()) {
         Error("Add(TVectorT &,const TMatrixT &,const TVectorT &)","target vector and matrix are incompatible");
         return target;
      }
   }

   const Element * const sp = source.GetMatrixArray();  // sources vector ptr
   const Element *       mp = a.GetMatrixArray();       // Matrix row ptr
         Element *       tp = target.GetMatrixArray();  // Target vector ptr
#ifdef CBLAS
   if (typeid(Element) == typeid(Double_t))
      cblas_dsymv(CblasRowMajor,CblasUpper,fNrows,1.0,mp,
                  fNrows,sp,1,0.0,tp,1);
   else if (typeid(Element) != typeid(Float_t))
      cblas_ssymv(CblasRowMajor,CblasUpper,fNrows,1.0,mp,
                  fNrows,sp,1,0.0,tp,1);
   else
      Error("operator*=","type %s not implemented in BLAS library",typeid(Element));
#else
   const Element * const sp_last = sp+source.GetNrows();
   const Element * const tp_last = tp+target.GetNrows();
   if (scalar == 1.0) {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++ += sum;
      }
   } else if (scalar == 0.0) {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++  = sum;
      }
   } else if (scalar == -1.0) {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++ -= sum;
      }
   } else {
      while (tp < tp_last) {
         const Element *sp1 = sp;
         Element sum = 0;
         while (sp1 < sp_last)
            sum += *sp1++ * *mp++;
         *tp++ += scalar * sum;
      }
   }
   R__ASSERT(mp == a.GetMatrixArray()+a.GetNoElements());
#endif

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += A * source.
/// NOTE: in case scalar=0, do  target = A * source.

template<class Element>
TVectorT<Element> &Add(TVectorT<Element> &target,Element scalar,
                       const TMatrixTSparse<Element> &a,const TVectorT<Element> &source)
{
   if (gMatrixCheck) {
      R__ASSERT(target.IsValid());
      R__ASSERT(a.IsValid());
      if (a.GetNrows() != target.GetNrows() || a.GetRowLwb() != target.GetLwb()) {
         Error("Add(TVectorT &,const TMatrixT &,const TVectorT &)","target vector and matrix are incompatible");
         return target;
      }

      R__ASSERT(source.IsValid());
      if (a.GetNcols() != source.GetNrows() || a.GetColLwb() != source.GetLwb()) {
         Error("Add(TVectorT &,const TMatrixT &,const TVectorT &)","source vector and matrix are incompatible");
         return target;
      }
   }

   const Int_t   * const pRowIndex = a.GetRowIndexArray();
   const Int_t   * const pColIndex = a.GetColIndexArray();
   const Element * const mp        = a.GetMatrixArray();     // Matrix row ptr

   const Element * const sp = source.GetMatrixArray(); // Source vector ptr
         Element *       tp = target.GetMatrixArray(); // Target vector ptr

   if (scalar == 1.0) {
      for (Int_t irow = 0; irow < a.GetNrows(); irow++) {
         const Int_t sIndex = pRowIndex[irow];
         const Int_t eIndex = pRowIndex[irow+1];
         Element sum = 0.0;
         for (Int_t index = sIndex; index < eIndex; index++) {
            const Int_t icol = pColIndex[index];
            sum += mp[index]*sp[icol];
         }
         tp[irow] += sum;
      }
   } else if (scalar == 0.0) {
      for (Int_t irow = 0; irow < a.GetNrows(); irow++) {
         const Int_t sIndex = pRowIndex[irow];
         const Int_t eIndex = pRowIndex[irow+1];
         Element sum = 0.0;
         for (Int_t index = sIndex; index < eIndex; index++) {
            const Int_t icol = pColIndex[index];
            sum += mp[index]*sp[icol];
         }
         tp[irow]  = sum;
      }
   } else if (scalar == -1.0) {
     for (Int_t irow = 0; irow < a.GetNrows(); irow++) {
        const Int_t sIndex = pRowIndex[irow];
        const Int_t eIndex = pRowIndex[irow+1];
        Element sum = 0.0;
        for (Int_t index = sIndex; index < eIndex; index++) {
           const Int_t icol = pColIndex[index];
           sum += mp[index]*sp[icol];
        }
        tp[irow] -= sum;
      }
   } else {
      for (Int_t irow = 0; irow < a.GetNrows(); irow++) {
        const Int_t sIndex = pRowIndex[irow];
        const Int_t eIndex = pRowIndex[irow+1];
        Element sum = 0.0;
        for (Int_t index = sIndex; index < eIndex; index++) {
           const Int_t icol = pColIndex[index];
           sum += mp[index]*sp[icol];
        }
        tp[irow] += scalar * sum;
      }
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += scalar * ElementMult(source1,source2) .

template<class Element>
TVectorT<Element> &AddElemMult(TVectorT<Element> &target,Element scalar,
                      const TVectorT<Element> &source1,const TVectorT<Element> &source2)
{
   if (gMatrixCheck && !(AreCompatible(target,source1) && AreCompatible(target,source2))) {
      Error("AddElemMult(TVectorT<Element> &,Element,const TVectorT<Element> &,const TVectorT<Element> &)",
             "vector's are incompatible");
      return target;
   }

   const Element *       sp1 = source1.GetMatrixArray();
   const Element *       sp2 = source2.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();

   if (scalar == 1.0 ) {
      while ( tp < ftp )
         *tp++ += *sp1++ * *sp2++;
   } else if (scalar == -1.0) {
      while ( tp < ftp )
         *tp++ -= *sp1++ * *sp2++;
   } else {
      while ( tp < ftp )
         *tp++ += scalar * *sp1++ * *sp2++;
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += scalar * ElementMult(source1,source2) only for those elements
/// where select[i] != 0.0

template<class Element>
TVectorT<Element> &AddElemMult(TVectorT<Element> &target,Element scalar,
                      const TVectorT<Element> &source1,const TVectorT<Element> &source2,const TVectorT<Element> &select)
{
   if (gMatrixCheck && !( AreCompatible(target,source1) && AreCompatible(target,source2) &&
          AreCompatible(target,select) )) {
      Error("AddElemMult(TVectorT<Element> &,Element,const TVectorT<Element> &,const TVectorT<Element> &,onst TVectorT<Element> &)",
             "vector's are incompatible");
      return target;
   }

   const Element *       sp1 = source1.GetMatrixArray();
   const Element *       sp2 = source2.GetMatrixArray();
   const Element *       mp  = select.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();

   if (scalar == 1.0 ) {
      while ( tp < ftp ) {
         if (*mp) *tp += *sp1 * *sp2;
         mp++; tp++; sp1++; sp2++;
      }
   } else if (scalar == -1.0) {
      while ( tp < ftp ) {
         if (*mp) *tp -= *sp1 * *sp2;
         mp++; tp++; sp1++; sp2++;
      }
   } else {
      while ( tp < ftp ) {
         if (*mp) *tp += scalar * *sp1 * *sp2;
         mp++; tp++; sp1++; sp2++;
      }
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += scalar * ElementDiv(source1,source2) .

template<class Element>
TVectorT<Element> &AddElemDiv(TVectorT<Element> &target,Element scalar,
                     const TVectorT<Element> &source1,const TVectorT<Element> &source2)
{
   if (gMatrixCheck && !(AreCompatible(target,source1) && AreCompatible(target,source2))) {
      Error("AddElemDiv(TVectorT<Element> &,Element,const TVectorT<Element> &,const TVectorT<Element> &)",
             "vector's are incompatible");
      return target;
   }

   const Element *       sp1 = source1.GetMatrixArray();
   const Element *       sp2 = source2.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();

   if (scalar == 1.0 ) {
      while ( tp < ftp ) {
         if (*sp2 != 0.0)
            *tp += *sp1 / *sp2;
         else {
            const Int_t irow = (sp2-source2.GetMatrixArray())/source2.GetNrows();
            Error("AddElemDiv","source2 (%d) is zero",irow);
         }
         tp++; sp1++; sp2++;
      }
   } else if (scalar == -1.0) {
      while ( tp < ftp ) {
         if (*sp2 != 0.0)
            *tp -= *sp1 / *sp2;
         else {
            const Int_t irow = (sp2-source2.GetMatrixArray())/source2.GetNrows();
            Error("AddElemDiv","source2 (%d) is zero",irow);
         }
         tp++; sp1++; sp2++;
      }
   } else {
      while ( tp < ftp ) {
         if (*sp2 != 0.0)
            *tp += scalar * *sp1 / *sp2;
         else {
            const Int_t irow = (sp2-source2.GetMatrixArray())/source2.GetNrows();
            Error("AddElemDiv","source2 (%d) is zero",irow);
         }
         tp++; sp1++; sp2++;
      }
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify addition: target += scalar * ElementDiv(source1,source2) only for those elements
/// where select[i] != 0.0

template<class Element>
TVectorT<Element> &AddElemDiv(TVectorT<Element> &target,Element scalar,
                     const TVectorT<Element> &source1,const TVectorT<Element> &source2,const TVectorT<Element> &select)
{
   if (gMatrixCheck && !( AreCompatible(target,source1) && AreCompatible(target,source2) &&
          AreCompatible(target,select) )) {
      Error("AddElemDiv(TVectorT<Element> &,Element,const TVectorT<Element> &,const TVectorT<Element> &,onst TVectorT<Element> &)",
             "vector's are incompatible");
      return target;
   }

   const Element *       sp1 = source1.GetMatrixArray();
   const Element *       sp2 = source2.GetMatrixArray();
   const Element *       mp  = select.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();

   if (scalar == 1.0 ) {
      while ( tp < ftp ) {
         if (*mp) {
            if (*sp2 != 0.0)
               *tp += *sp1 / *sp2;
            else {
               const Int_t irow = (sp2-source2.GetMatrixArray())/source2.GetNrows();
               Error("AddElemDiv","source2 (%d) is zero",irow);
            }
         }
         mp++; tp++; sp1++; sp2++;
      }
   } else if (scalar == -1.0) {
      while ( tp < ftp ) {
         if (*mp) {
            if (*sp2 != 0.0)
               *tp -= *sp1 / *sp2;
            else {
               const Int_t irow = (sp2-source2.GetMatrixArray())/source2.GetNrows();
               Error("AddElemDiv","source2 (%d) is zero",irow);
            }
         }
         mp++; tp++; sp1++; sp2++;
      }
   } else {
      while ( tp < ftp ) {
         if (*mp) {
            if (*sp2 != 0.0)
               *tp += scalar * *sp1 / *sp2;
            else {
               const Int_t irow = (sp2-source2.GetMatrixArray())/source2.GetNrows();
               Error("AddElemDiv","source2 (%d) is zero",irow);
            }
         }
         mp++; tp++; sp1++; sp2++;
      }
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply target by the source, element-by-element.

template<class Element>
TVectorT<Element> &ElementMult(TVectorT<Element> &target,const TVectorT<Element> &source)
{
   if (gMatrixCheck && !AreCompatible(target,source)) {
      Error("ElementMult(TVectorT<Element> &,const TVectorT<Element> &)","vector's are incompatible");
      return target;
   }

   const Element *       sp  = source.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();
   while ( tp < ftp )
      *tp++ *= *sp++;

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply target by the source, element-by-element only where select[i] != 0.0

template<class Element>
TVectorT<Element> &ElementMult(TVectorT<Element> &target,const TVectorT<Element> &source,const TVectorT<Element> &select)
{
   if (gMatrixCheck && !(AreCompatible(target,source) && AreCompatible(target,select))) {
      Error("ElementMult(TVectorT<Element> &,const TVectorT<Element> &,const TVectorT<Element> &)","vector's are incompatible");
      return target;
   }

   const Element *       sp  = source.GetMatrixArray();
   const Element *       mp  = select.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();
   while ( tp < ftp ) {
      if (*mp) *tp *= *sp;
      mp++; tp++; sp++;
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide target by the source, element-by-element.

template<class Element>
TVectorT<Element> &ElementDiv(TVectorT<Element> &target,const TVectorT<Element> &source)
{
   if (gMatrixCheck && !AreCompatible(target,source)) {
      Error("ElementDiv(TVectorT<Element> &,const TVectorT<Element> &)","vector's are incompatible");
      return target;
   }

   const Element *       sp  = source.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();
   while ( tp < ftp ) {
      if (*sp  != 0.0)
         *tp++ /= *sp++;
      else {
         const Int_t irow = (sp-source.GetMatrixArray())/source.GetNrows();
         Error("ElementDiv","source (%d) is zero",irow);
      }
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide target by the source, element-by-element only where select[i] != 0.0

template<class Element>
TVectorT<Element> &ElementDiv(TVectorT<Element> &target,const TVectorT<Element> &source,const TVectorT<Element> &select)
{
   if (gMatrixCheck && !AreCompatible(target,source)) {
      Error("ElementDiv(TVectorT<Element> &,const TVectorT<Element> &,const TVectorT<Element> &)","vector's are incompatible");
      return target;
   }

   const Element *       sp  = source.GetMatrixArray();
   const Element *       mp  = select.GetMatrixArray();
         Element *       tp  = target.GetMatrixArray();
   const Element * const ftp = tp+target.GetNrows();
   while ( tp < ftp ) {
      if (*mp) {
         if (*sp != 0.0)
            *tp /= *sp;
         else {
            const Int_t irow = (sp-source.GetMatrixArray())/source.GetNrows();
            Error("ElementDiv","source (%d) is zero",irow);
         }
      }
      mp++; tp++; sp++;
   }

   return target;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if v1 and v2 are both valid and have the same shape

template<class Element1,class Element2>
Bool_t AreCompatible(const TVectorT<Element1> &v1,const TVectorT<Element2> &v2,Int_t verbose)
{
   if (!v1.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "vector 1 not valid");
      return kFALSE;
   }
   if (!v2.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "vector 2 not valid");
      return kFALSE;
   }

   if (v1.GetNrows() != v2.GetNrows() || v1.GetLwb() != v2.GetLwb()) {
      if (verbose)
         ::Error("AreCompatible", "matrices 1 and 2 not compatible");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if m and v are both valid and have compatible shapes for M * v

template<class Element1, class Element2>
Bool_t AreCompatible(const TMatrixT<Element1> &m,const TVectorT<Element2> &v,Int_t verbose)
{
   if (!m.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "Matrix not valid");
      return kFALSE;
   }
   if (!v.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "vector not valid");
      return kFALSE;
   }

   if (m.GetNcols() != v.GetNrows() ) {
      if (verbose)
         ::Error("AreCompatible", "matrix and vector not compatible");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if m and v are both valid and have compatible shapes for v * M

template<class Element1, class Element2>
Bool_t AreCompatible(const TVectorT<Element1> &v,const TMatrixT<Element2> &m,Int_t verbose)
{
   if (!m.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "Matrix not valid");
      return kFALSE;
   }
   if (!v.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "vector not valid");
      return kFALSE;
   }

   if (v.GetNrows() != m.GetNrows() ) {
      if (verbose)
         ::Error("AreCompatible", "vector and matrix not compatible");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compare two vectors and print out the result of the comparison.

template<class Element>
void Compare(const TVectorT<Element> &v1,const TVectorT<Element> &v2)
{
   if (!AreCompatible(v1,v2)) {
      Error("Compare(const TVectorT<Element> &,const TVectorT<Element> &)","vectors are incompatible");
      return;
   }

   printf("\n\nComparison of two TVectorTs:\n");

   Element norm1  = 0;       // Norm of the Matrices
   Element norm2  = 0;       // Norm of the Matrices
   Element ndiff  = 0;       // Norm of the difference
   Int_t    imax   = 0;       // For the elements that differ most
   Element difmax = -1;
   const Element *mp1 = v1.GetMatrixArray();    // Vector element pointers
   const Element *mp2 = v2.GetMatrixArray();

   for (Int_t i = 0; i < v1.GetNrows(); i++) {
      const Element mv1  = *mp1++;
      const Element mv2  = *mp2++;
      const Element diff = TMath::Abs(mv1-mv2);

      if (diff > difmax) {
         difmax = diff;
         imax = i;
      }
      norm1 += TMath::Abs(mv1);
      norm2 += TMath::Abs(mv2);
      ndiff += TMath::Abs(diff);
   }

   imax += v1.GetLwb();
   printf("\nMaximal discrepancy    \t\t%g",difmax);
   printf("\n   occured at the point\t\t(%d)",imax);
   const Element mv1 = v1(imax);
   const Element mv2 = v2(imax);
   printf("\n Vector 1 element is    \t\t%g",mv1);
   printf("\n Vector 2 element is    \t\t%g",mv2);
   printf("\n Absolute error v2[i]-v1[i]\t\t%g",mv2-mv1);
   printf("\n Relative error\t\t\t\t%g\n",
          (mv2-mv1)/TMath::Max(TMath::Abs(mv2+mv1)/2,(Element)1e-7));

   printf("\n||Vector 1||   \t\t\t%g",norm1);
   printf("\n||Vector 2||   \t\t\t%g",norm2);
   printf("\n||Vector1-Vector2||\t\t\t\t%g",ndiff);
   printf("\n||Vector1-Vector2||/sqrt(||Vector1|| ||Vector2||)\t%g\n\n",
          ndiff/TMath::Max(TMath::Sqrt(norm1*norm2),1e-7));
}

////////////////////////////////////////////////////////////////////////////////
/// Validate that all elements of vector have value val within maxDevAllow .

template<class Element>
Bool_t VerifyVectorValue(const TVectorT<Element> &v,Element val,
                         Int_t verbose,Element maxDevAllow)
{
   Int_t   imax      = 0;
   Element maxDevObs = 0;

   if (TMath::Abs(maxDevAllow) <= 0.0)
      maxDevAllow = std::numeric_limits<Element>::epsilon();

   for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++) {
      const Element dev = TMath::Abs(v(i)-val);
      if (dev > maxDevObs) {
         imax      = i;
         maxDevObs = dev;
      }
   }

   if (maxDevObs == 0)
      return kTRUE;

   if (verbose) {
      printf("Largest dev for (%d); dev = |%g - %g| = %g\n",imax,v(imax),val,maxDevObs);
      if (maxDevObs > maxDevAllow)
         Error("VerifyVectorValue","Deviation > %g\n",maxDevAllow);
   }

   if (maxDevObs > maxDevAllow)
      return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Verify that elements of the two vectors are equal within maxDevAllow .

template<class Element>
Bool_t VerifyVectorIdentity(const TVectorT<Element> &v1,const TVectorT<Element> &v2,
                            Int_t verbose, Element maxDevAllow)
{
   Int_t   imax      = 0;
   Element maxDevObs = 0;

   if (!AreCompatible(v1,v2))
      return kFALSE;

   if (TMath::Abs(maxDevAllow) <= 0.0)
      maxDevAllow = std::numeric_limits<Element>::epsilon();

   for (Int_t i = v1.GetLwb(); i <= v1.GetUpb(); i++) {
      const Element dev = TMath::Abs(v1(i)-v2(i));
      if (dev > maxDevObs) {
         imax      = i;
         maxDevObs = dev;
      }
   }

   if (maxDevObs == 0)
      return kTRUE;

   if (verbose) {
      printf("Largest dev for (%d); dev = |%g - %g| = %g\n",imax,v1(imax),v2(imax),maxDevObs);
      if (maxDevObs > maxDevAllow)
         Error("VerifyVectorIdentity","Deviation > %g\n",maxDevAllow);
   }

   if (maxDevObs > maxDevAllow) {
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TVectorT.

template<class Element>
void TVectorT<Element>::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s,&R__c);
      if (R__v > 1) {
         Clear();
         R__b.ReadClassBuffer(TVectorT<Element>::Class(),this,R__v,R__s,R__c);
      } else { //====process old versions before automatic schema evolution
         TObject::Streamer(R__b);
         R__b >> fRowLwb;
         fNrows = R__b.ReadArray(fElements);
         R__b.CheckByteCount(R__s, R__c, TVectorT<Element>::IsA());
      }
      if (fNrows > 0 && fNrows <= kSizeMax) {
         memcpy(fDataStack,fElements,fNrows*sizeof(Element));
         delete [] fElements;
         fElements = fDataStack;
      } else if (fNrows < 0)
         Invalidate();

      if (R__v < 3)
         MakeValid();
    } else {
       R__b.WriteClassBuffer(TVectorT<Element>::Class(),this);
   }
}

#include "TMatrixFfwd.h"
#include "TMatrixFSymfwd.h"
#include "TMatrixFSparsefwd.h"

template class TVectorT<Float_t>;

template Bool_t    operator==          <Float_t>(const TVectorF       &source1,const TVectorF &source2);
template TVectorF  operator+           <Float_t>(const TVectorF       &source1,const TVectorF &source2);
template TVectorF  operator-           <Float_t>(const TVectorF       &source1,const TVectorF &source2);
template Float_t   operator*           <Float_t>(const TVectorF       &source1,const TVectorF &source2);
template TVectorF  operator*           <Float_t>(const TMatrixF       &a,      const TVectorF &source);
template TVectorF  operator*           <Float_t>(const TMatrixFSym    &a,      const TVectorF &source);
template TVectorF  operator*           <Float_t>(const TMatrixFSparse &a,      const TVectorF &source);
template TVectorF  operator*           <Float_t>(      Float_t         val,    const TVectorF &source);

template Float_t   Dot                 <Float_t>(const TVectorF       &v1,     const TVectorF &v2);
template TMatrixF  OuterProduct        <Float_t,Float_t>
                                                (const TVectorF       &v1,     const TVectorF &v2);
template TMatrixF &OuterProduct        <Float_t,Float_t,Float_t>
                                                (      TMatrixF       &target, const TVectorF &v1,     const TVectorF       &v2);
template Float_t   Mult                <Float_t,Float_t,Float_t>
                                                (const TVectorF       &v1,     const TMatrixF &m,      const TVectorF       &v2);

template TVectorF &Add                 <Float_t>(      TVectorF       &target,       Float_t   scalar, const TVectorF       &source);
template TVectorF &Add                 <Float_t>(      TVectorF       &target,       Float_t   scalar, const TMatrixF       &a,
                                                                               const TVectorF &source);
template TVectorF &Add                 <Float_t>(      TVectorF       &target,       Float_t   scalar, const TMatrixFSym    &a,
                                                                               const TVectorF &source);
template TVectorF &Add                 <Float_t>(      TVectorF       &target,       Float_t   scalar, const TMatrixFSparse &a,
                                                                               const TVectorF &source);
template TVectorF &AddElemMult         <Float_t>(      TVectorF       &target,       Float_t   scalar, const TVectorF       &source1,
                                                                               const TVectorF &source2);
template TVectorF &AddElemMult         <Float_t>(      TVectorF       &target,       Float_t   scalar, const TVectorF       &source1,
                                                                               const TVectorF &source2,const TVectorF       &select);
template TVectorF &AddElemDiv          <Float_t>(      TVectorF       &target,       Float_t   scalar, const TVectorF       &source1,
                                                                               const TVectorF &source2);
template TVectorF &AddElemDiv          <Float_t>(      TVectorF       &target,       Float_t   scalar, const TVectorF       &source1,
                                                                               const TVectorF &source2,const TVectorF       &select);
template TVectorF &ElementMult         <Float_t>(      TVectorF       &target, const TVectorF &source);
template TVectorF &ElementMult         <Float_t>(      TVectorF       &target, const TVectorF &source, const TVectorF       &select);
template TVectorF &ElementDiv          <Float_t>(      TVectorF       &target, const TVectorF &source);
template TVectorF &ElementDiv          <Float_t>(      TVectorF       &target, const TVectorF &source, const TVectorF       &select);

template Bool_t    AreCompatible       <Float_t,Float_t> (const TVectorF &v1,const TVectorF &v2,Int_t verbose);
template Bool_t    AreCompatible       <Float_t,Double_t>(const TVectorF &v1,const TVectorD &v2,Int_t verbose);
template Bool_t    AreCompatible       <Float_t,Float_t> (const TMatrixF &m, const TVectorF &v, Int_t verbose);
template Bool_t    AreCompatible       <Float_t,Float_t> (const TVectorF &v, const TMatrixF &m, Int_t verbose);

template void      Compare             <Float_t>         (const TVectorF &v1,const TVectorF &v2);
template Bool_t    VerifyVectorValue   <Float_t>         (const TVectorF &m,       Float_t   val,Int_t verbose,Float_t maxDevAllow);
template Bool_t    VerifyVectorIdentity<Float_t>         (const TVectorF &m1,const TVectorF &m2, Int_t verbose,Float_t maxDevAllow);

#include "TMatrixDfwd.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixDSparsefwd.h"

template class TVectorT<Double_t>;

template Bool_t    operator==          <Double_t>(const TVectorD       &source1,const TVectorD &source2);
template TVectorD  operator+           <Double_t>(const TVectorD       &source1,const TVectorD &source2);
template TVectorD  operator-           <Double_t>(const TVectorD       &source1,const TVectorD &source2);
template Double_t  operator*           <Double_t>(const TVectorD       &source1,const TVectorD &source2);
template TVectorD  operator*           <Double_t>(const TMatrixD       &a,      const TVectorD &source);
template TVectorD  operator*           <Double_t>(const TMatrixDSym    &a,      const TVectorD &source);
template TVectorD  operator*           <Double_t>(const TMatrixDSparse &a,      const TVectorD &source);
template TVectorD  operator*           <Double_t>(      Double_t        val,    const TVectorD &source);

template Double_t  Dot                 <Double_t>(const TVectorD       &v1,     const TVectorD &v2);
template TMatrixD  OuterProduct        <Double_t,Double_t>
                                                 (const TVectorD       &v1,     const TVectorD &v2);
template TMatrixD &OuterProduct        <Double_t,Double_t,Double_t>
                                                 (      TMatrixD       &target, const TVectorD &v1,     const TVectorD       &v2);
template Double_t  Mult                <Double_t,Double_t,Double_t>
                                                 (const TVectorD       &v1,     const TMatrixD &m,      const TVectorD       &v2);

template TVectorD &Add                 <Double_t>(      TVectorD       &target,       Double_t  scalar, const TVectorD       &source);
template TVectorD &Add                 <Double_t>(      TVectorD       &target,       Double_t  scalar, const TMatrixD       &a,
                                                                                const TVectorD &source);
template TVectorD &Add                 <Double_t>(      TVectorD       &target,       Double_t  scalar, const TMatrixDSym    &a
                                                                                ,      const TVectorD &source);
template TVectorD &Add                 <Double_t>(      TVectorD       &target,       Double_t  scalar, const TMatrixDSparse &a
                                                                                ,      const TVectorD &source);
template TVectorD &AddElemMult         <Double_t>(      TVectorD       &target,       Double_t  scalar, const TVectorD       &source1,
                                                                                const TVectorD &source2);
template TVectorD &AddElemMult         <Double_t>(      TVectorD       &target,       Double_t  scalar, const TVectorD       &source1,
                                                                                const TVectorD &source2,const TVectorD       &select);
template TVectorD &AddElemDiv          <Double_t>(      TVectorD       &target,       Double_t  scalar, const TVectorD       &source1,
                                                                                const TVectorD &source2);
template TVectorD &AddElemDiv          <Double_t>(      TVectorD       &target,       Double_t  scalar, const TVectorD       &source1,
                                                                                const TVectorD &source2,const TVectorD       &select);
template TVectorD &ElementMult         <Double_t>(      TVectorD       &target, const TVectorD &source);
template TVectorD &ElementMult         <Double_t>(      TVectorD       &target, const TVectorD &source, const TVectorD       &select);
template TVectorD &ElementDiv          <Double_t>(      TVectorD       &target, const TVectorD &source);
template TVectorD &ElementDiv          <Double_t>(      TVectorD       &target, const TVectorD &source, const TVectorD       &select);

template Bool_t    AreCompatible       <Double_t,Double_t>(const TVectorD &v1,const TVectorD &v2,Int_t verbose);
template Bool_t    AreCompatible       <Double_t,Float_t> (const TVectorD &v1,const TVectorF &v2,Int_t verbose);
template Bool_t    AreCompatible       <Double_t,Double_t>(const TMatrixD &m, const TVectorD &v, Int_t verbose);
template Bool_t    AreCompatible       <Double_t,Double_t>(const TVectorD &v, const TMatrixD &m, Int_t verbose);

template void      Compare             <Double_t>         (const TVectorD &v1,const TVectorD &v2);
template Bool_t    VerifyVectorValue   <Double_t>         (const TVectorD &m,       Double_t  val,Int_t verbose,Double_t maxDevAllow);
template Bool_t    VerifyVectorIdentity<Double_t>         (const TVectorD &m1,const TVectorD &m2, Int_t verbose,Double_t maxDevAllow);
