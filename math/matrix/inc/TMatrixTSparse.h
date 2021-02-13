// @(#)root/matrix:$Id$
// Authors: Fons Rademakers, Eddy Offermann   Feb 2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixTSparse
#define ROOT_TMatrixTSparse

#include "TMatrixTBase.h"
#include "TMatrixTUtils.h"

#include <cstring>

#ifdef CBLAS
#include <vecLib/vBLAS.h>
//#include <cblas.h>
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrixTSparse                                                       //
//                                                                      //
// Template class of a general sparse matrix in the Harwell-Boeing      //
// format                                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

template<class Element> class TMatrixT;

template<class Element> class TMatrixTSparse : public TMatrixTBase<Element> {

protected:

   Int_t   *fRowIndex;  //[fNrowIndex] row index
   Int_t   *fColIndex;  //[fNelems]    column index
   Element *fElements;  //[fNelems]

   void Allocate(Int_t nrows,Int_t ncols,Int_t row_lwb = 0,Int_t col_lwb = 0,
                 Int_t init = 0,Int_t nr_nonzeros = 0);

  // Elementary constructors
   void AMultB (const TMatrixTSparse<Element> &a,const TMatrixTSparse<Element> &b,Int_t constr=0) {
                const TMatrixTSparse<Element> bt(TMatrixTSparse::kTransposed,b); AMultBt(a,bt,constr); }
   void AMultB (const TMatrixTSparse<Element> &a,const TMatrixT<Element>       &b,Int_t constr=0) {
                const TMatrixTSparse<Element> bsp = b;
                const TMatrixTSparse<Element> bt(TMatrixTSparse::kTransposed,bsp); AMultBt(a,bt,constr); }
   void AMultB (const TMatrixT<Element>       &a,const TMatrixTSparse<Element> &b,Int_t constr=0) {
                const TMatrixTSparse<Element> bt(TMatrixTSparse::kTransposed,b); AMultBt(a,bt,constr); }

   void AMultBt(const TMatrixTSparse<Element> &a,const TMatrixTSparse<Element> &b,Int_t constr=0);
   void AMultBt(const TMatrixTSparse<Element> &a,const TMatrixT<Element>       &b,Int_t constr=0);
   void AMultBt(const TMatrixT<Element>       &a,const TMatrixTSparse<Element> &b,Int_t constr=0);

   void APlusB (const TMatrixTSparse<Element> &a,const TMatrixTSparse<Element> &b,Int_t constr=0);
   void APlusB (const TMatrixTSparse<Element> &a,const TMatrixT<Element>       &b,Int_t constr=0);
   void APlusB (const TMatrixT<Element>       &a,const TMatrixTSparse<Element> &b,Int_t constr=0) { APlusB(b,a,constr); }

   void AMinusB(const TMatrixTSparse<Element> &a,const TMatrixTSparse<Element> &b,Int_t constr=0);
   void AMinusB(const TMatrixTSparse<Element> &a,const TMatrixT<Element>       &b,Int_t constr=0);
   void AMinusB(const TMatrixT<Element>       &a,const TMatrixTSparse<Element> &b,Int_t constr=0);

public:

   enum EMatrixCreatorsOp1 { kZero,kUnit,kTransposed,kAtA };
   enum EMatrixCreatorsOp2 { kMult,kMultTranspose,kPlus,kMinus };

   TMatrixTSparse() { fElements = 0; fRowIndex = 0; fColIndex = 0; }
   TMatrixTSparse(Int_t nrows,Int_t ncols);
   TMatrixTSparse(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb);
   TMatrixTSparse(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,Int_t nr_nonzeros,
                  Int_t *row, Int_t *col,Element *data);
   TMatrixTSparse(const TMatrixTSparse<Element> &another);
   TMatrixTSparse(const TMatrixT<Element>       &another);

   TMatrixTSparse(EMatrixCreatorsOp1 op,const TMatrixTSparse<Element> &prototype);
   TMatrixTSparse(const TMatrixTSparse<Element> &a,EMatrixCreatorsOp2 op,const TMatrixTSparse<Element> &b);
   TMatrixTSparse(const TMatrixTSparse<Element> &a,EMatrixCreatorsOp2 op,const TMatrixT      <Element> &b);
   TMatrixTSparse(const TMatrixT      <Element> &a,EMatrixCreatorsOp2 op,const TMatrixTSparse<Element> &b);

   virtual ~TMatrixTSparse() { TMatrixTSparse::Clear(); }

   virtual const Element *GetMatrixArray  () const;
   virtual       Element *GetMatrixArray  ();
   virtual const Int_t    *GetRowIndexArray() const;
   virtual       Int_t    *GetRowIndexArray();
   virtual const Int_t    *GetColIndexArray() const;
   virtual       Int_t    *GetColIndexArray();

   virtual TMatrixTBase<Element>   &SetRowIndexArray(Int_t *data) { memmove(fRowIndex,data,(this->fNrows+1)*sizeof(Int_t)); return *this; }
   virtual TMatrixTBase<Element>   &SetColIndexArray(Int_t *data) { memmove(fColIndex,data,this->fNelems*sizeof(Int_t)); return *this; }

           TMatrixTSparse<Element> &SetSparseIndex  (Int_t nelem_new);
           TMatrixTSparse<Element> &SetSparseIndex  (const TMatrixTBase<Element> &another);
           TMatrixTSparse<Element> &SetSparseIndexAB(const TMatrixTSparse<Element> &a,const TMatrixTSparse<Element> &b);
           TMatrixTSparse<Element> &SetSparseIndexAB(const TMatrixT      <Element> &a,const TMatrixTSparse<Element> &b);
           TMatrixTSparse<Element> &SetSparseIndexAB(const TMatrixTSparse<Element> &a,const TMatrixT      <Element> &b)
                                              { return SetSparseIndexAB(b,a); }

   virtual void                     GetMatrix2Array (Element *data,Option_t * /*option*/ ="") const;
   virtual TMatrixTBase<Element>   &SetMatrixArray  (const Element *data,Option_t * /*option*/="")
                                                    { memcpy(fElements,data,this->fNelems*sizeof(Element)); return *this; }
   virtual TMatrixTBase<Element>   &SetMatrixArray  (Int_t nr_nonzeros,Int_t *irow,Int_t *icol,Element *data);
   virtual TMatrixTBase<Element>   &InsertRow       (Int_t row,Int_t col,const Element *v,Int_t n=-1);
   virtual void                     ExtractRow      (Int_t row,Int_t col,      Element *v,Int_t n=-1) const;

   virtual TMatrixTBase<Element>   &ResizeTo(Int_t nrows,Int_t ncols,Int_t nr_nonzeros=-1);
   virtual TMatrixTBase<Element>   &ResizeTo(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,Int_t nr_nonzeros=-1);
   inline  TMatrixTBase<Element>   &ResizeTo(const TMatrixTSparse<Element> &m) {return ResizeTo(m.GetRowLwb(),m.GetRowUpb(),m.GetColLwb(),
                                                                                                m.GetColUpb(),m.GetNoElements()); }

   virtual void Clear(Option_t * /*option*/ ="") { if (this->fIsOwner) {
                                                      if (fElements) { delete [] fElements; fElements = 0; }
                                                      if (fRowIndex) { delete [] fRowIndex; fRowIndex = 0; }
                                                      if (fColIndex) { delete [] fColIndex; fColIndex = 0; }
                                                   }
                                                   this->fNelems    = 0;
                                                   this->fNrowIndex = 0;
                                                 }

           TMatrixTSparse<Element> &Use   (Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,Int_t nr_nonzeros,
                                           Int_t *pRowIndex,Int_t *pColIndex,Element *pData);
   const   TMatrixTSparse<Element> &Use   (Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,Int_t nr_nonzeros,
                                           const Int_t *pRowIndex,const Int_t *pColIndex,const Element *pData) const
                                            { return (const TMatrixTSparse<Element>&)
                                                     ((const_cast<TMatrixTSparse<Element> *>(this))->Use(row_lwb,row_upb,col_lwb,col_upb,nr_nonzeros,
                                                                                             const_cast<Int_t *>(pRowIndex),
                                                                                             const_cast<Int_t *>(pColIndex),
                                                                                             const_cast<Element *>(pData))); }
           TMatrixTSparse<Element> &Use   (Int_t nrows,Int_t ncols,Int_t nr_nonzeros,
                                           Int_t *pRowIndex,Int_t *pColIndex,Element *pData);
   const   TMatrixTSparse<Element> &Use   (Int_t nrows,Int_t ncols,Int_t nr_nonzeros,
                                           const Int_t *pRowIndex,const Int_t *pColIndex,const Element *pData) const;
           TMatrixTSparse<Element> &Use   (TMatrixTSparse<Element> &a);
   const   TMatrixTSparse<Element> &Use   (const TMatrixTSparse<Element> &a) const;

   virtual TMatrixTBase<Element>   &GetSub(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,
                                            TMatrixTBase<Element> &target,Option_t *option="S") const;
           TMatrixTSparse<Element>  GetSub(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,Option_t *option="S") const;
   virtual TMatrixTBase<Element>   &SetSub(Int_t row_lwb,Int_t col_lwb,const TMatrixTBase<Element> &source);

   virtual Bool_t IsSymmetric() const { return (*this == TMatrixTSparse<Element>(kTransposed,*this)); }
   TMatrixTSparse<Element> &Transpose (const TMatrixTSparse<Element> &source);
   inline TMatrixTSparse<Element> &T () { return this->Transpose(*this); }

   inline void Mult(const TMatrixTSparse<Element> &a,const TMatrixTSparse<Element> &b) { AMultB(a,b,0); }

   virtual TMatrixTBase<Element> &Zero       ();
   virtual TMatrixTBase<Element> &UnitMatrix ();

   virtual Element RowNorm () const;
   virtual Element ColNorm () const;
   virtual Int_t   NonZeros() const { return this->fNelems; }

   virtual TMatrixTBase<Element> &NormByDiag(const TVectorT<Element> &/*v*/,Option_t * /*option*/)
                                              { MayNotUse("NormByDiag"); return *this; }

   // Either access a_ij as a(i,j)
   Element  operator()(Int_t rown,Int_t coln) const;
   Element &operator()(Int_t rown,Int_t coln);

   // or as a[i][j]
   inline const TMatrixTSparseRow_const<Element> operator[](Int_t rown) const { return TMatrixTSparseRow_const<Element>(*this,rown); }
   inline       TMatrixTSparseRow      <Element> operator[](Int_t rown)       { return TMatrixTSparseRow      <Element>(*this,rown); }

   TMatrixTSparse<Element> &operator=(const TMatrixT<Element>       &source);
   TMatrixTSparse<Element> &operator=(const TMatrixTSparse<Element> &source);

   TMatrixTSparse<Element> &operator= (Element val);
   TMatrixTSparse<Element> &operator-=(Element val);
   TMatrixTSparse<Element> &operator+=(Element val);
   TMatrixTSparse<Element> &operator*=(Element val);

   TMatrixTSparse<Element> &operator+=(const TMatrixTSparse<Element> &source) { TMatrixTSparse<Element> tmp(*this); Clear();
                                                                                if (this == &source) APlusB (tmp,tmp,1);
                                                                                else                 APlusB (tmp,source,1);
                                                                                return *this; }
   TMatrixTSparse<Element> &operator+=(const TMatrixT<Element>       &source) { TMatrixTSparse<Element> tmp(*this); Clear();
                                                                                APlusB(tmp,source,1); return *this; }
   TMatrixTSparse<Element> &operator-=(const TMatrixTSparse<Element> &source) { TMatrixTSparse<Element> tmp(*this); Clear();
                                                                                if (this == &source) AMinusB (tmp,tmp,1);
                                                                                else                 AMinusB(tmp,source,1);
                                                                                return *this; }
   TMatrixTSparse<Element> &operator-=(const TMatrixT<Element>       &source) { TMatrixTSparse<Element> tmp(*this); Clear();
                                                                                AMinusB(tmp,source,1); return *this; }
   TMatrixTSparse<Element> &operator*=(const TMatrixTSparse<Element> &source) { TMatrixTSparse<Element> tmp(*this); Clear();
                                                                                if (this == &source) AMultB (tmp,tmp,1);
                                                                                else                 AMultB (tmp,source,1);
                                                                                return *this; }
   TMatrixTSparse<Element> &operator*=(const TMatrixT<Element>       &source) { TMatrixTSparse<Element> tmp(*this); Clear();
                                                                                AMultB(tmp,source,1);
                                                                                return *this; }

   virtual TMatrixTBase  <Element> &Randomize  (Element alpha,Element beta,Double_t &seed);
   virtual TMatrixTSparse<Element> &RandomizePD(Element alpha,Element beta,Double_t &seed);

   ClassDef(TMatrixTSparse,3) // Template of Sparse Matrix class
};

#ifndef __CINT__
// When building with -fmodules, it instantiates all pending instantiations,
// instead of delaying them until the end of the translation unit.
// We 'got away with' probably because the use and the definition of the
// explicit specialization do not occur in the same TU.
//
// In case we are building with -fmodules, we need to forward declare the
// specialization in order to compile the dictionary G__Matrix.cxx.
template <> TClass *TMatrixTSparse<double>::Class();
#endif // __CINT__

template <class Element> inline const Element *TMatrixTSparse<Element>::GetMatrixArray  () const { return fElements; }
template <class Element> inline       Element *TMatrixTSparse<Element>::GetMatrixArray  ()       { return fElements; }
template <class Element> inline const Int_t   *TMatrixTSparse<Element>::GetRowIndexArray() const { return fRowIndex; }
template <class Element> inline       Int_t   *TMatrixTSparse<Element>::GetRowIndexArray()       { return fRowIndex; }
template <class Element> inline const Int_t   *TMatrixTSparse<Element>::GetColIndexArray() const { return fColIndex; }
template <class Element> inline       Int_t   *TMatrixTSparse<Element>::GetColIndexArray()       { return fColIndex; }

template <class Element>
inline       TMatrixTSparse<Element> &TMatrixTSparse<Element>::Use   (Int_t nrows,Int_t ncols,Int_t nr_nonzeros,
                                                                      Int_t *pRowIndex,Int_t *pColIndex,Element *pData)
                                                                        { return Use(0,nrows-1,0,ncols-1,nr_nonzeros,pRowIndex,pColIndex,pData); }
template <class Element>
inline const TMatrixTSparse<Element> &TMatrixTSparse<Element>::Use   (Int_t nrows,Int_t ncols,Int_t nr_nonzeros,
                                                                      const Int_t *pRowIndex,const Int_t *pColIndex,const Element *pData) const
                                                                        { return Use(0,nrows-1,0,ncols-1,nr_nonzeros,pRowIndex,pColIndex,pData); }
template <class Element>
inline       TMatrixTSparse<Element> &TMatrixTSparse<Element>::Use   (TMatrixTSparse<Element> &a)
                                                                        { R__ASSERT(a.IsValid());
                                                                           return Use(a.GetRowLwb(),a.GetRowUpb(),a.GetColLwb(),a.GetColUpb(),
                                                                                      a.GetNoElements(),a.GetRowIndexArray(),
                                                                                      a.GetColIndexArray(),a.GetMatrixArray()); }
template <class Element>
inline const TMatrixTSparse<Element> &TMatrixTSparse<Element>::Use   (const TMatrixTSparse<Element> &a) const
                                                                        { R__ASSERT(a.IsValid());
                                                                           return Use(a.GetRowLwb(),a.GetRowUpb(),a.GetColLwb(),a.GetColUpb(),
                                                                                      a.GetNoElements(),a.GetRowIndexArray(),
                                                                                      a.GetColIndexArray(),a.GetMatrixArray()); }

template <class Element>
inline       TMatrixTSparse<Element>  TMatrixTSparse<Element>::GetSub(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,
                                                                      Option_t *option) const
                                                                        {
                                                                          TMatrixTSparse<Element> tmp;
                                                                          this->GetSub(row_lwb,row_upb,col_lwb,col_upb,tmp,option);
                                                                          return tmp;
                                                                        }

template <class Element> TMatrixTSparse<Element>  operator+ (const TMatrixTSparse<Element> &source1,const TMatrixTSparse<Element> &source2);
template <class Element> TMatrixTSparse<Element>  operator+ (const TMatrixTSparse<Element> &source1,const TMatrixT<Element>       &source2);
template <class Element> TMatrixTSparse<Element>  operator+ (const TMatrixT<Element>       &source1,const TMatrixTSparse<Element> &source2);
template <class Element> TMatrixTSparse<Element>  operator+ (const TMatrixTSparse<Element> &source ,      Element                  val    );
template <class Element> TMatrixTSparse<Element>  operator+ (      Element                  val    ,const TMatrixTSparse<Element> &source );
template <class Element> TMatrixTSparse<Element>  operator- (const TMatrixTSparse<Element> &source1,const TMatrixTSparse<Element> &source2);
template <class Element> TMatrixTSparse<Element>  operator- (const TMatrixTSparse<Element> &source1,const TMatrixT<Element>       &source2);
template <class Element> TMatrixTSparse<Element>  operator- (const TMatrixT<Element>       &source1,const TMatrixTSparse<Element> &source2);
template <class Element> TMatrixTSparse<Element>  operator- (const TMatrixTSparse<Element> &source ,      Element                  val    );
template <class Element> TMatrixTSparse<Element>  operator- (      Element                  val    ,const TMatrixTSparse<Element> &source );
template <class Element> TMatrixTSparse<Element>  operator* (const TMatrixTSparse<Element> &source1,const TMatrixTSparse<Element> &source2);
template <class Element> TMatrixTSparse<Element>  operator* (const TMatrixTSparse<Element> &source1,const TMatrixT<Element>       &source2);
template <class Element> TMatrixTSparse<Element>  operator* (const TMatrixT<Element>       &source1,const TMatrixTSparse<Element> &source2);
template <class Element> TMatrixTSparse<Element>  operator* (      Element                  val    ,const TMatrixTSparse<Element> &source );
template <class Element> TMatrixTSparse<Element>  operator* (const TMatrixTSparse<Element> &source,       Element                  val    );

template <class Element> TMatrixTSparse<Element> &Add        (TMatrixTSparse<Element> &target,      Element                   scalar,
                                                              const TMatrixTSparse<Element> &source);
template <class Element> TMatrixTSparse<Element> &ElementMult(TMatrixTSparse<Element> &target,const TMatrixTSparse<Element>  &source);
template <class Element> TMatrixTSparse<Element> &ElementDiv (TMatrixTSparse<Element> &target,const TMatrixTSparse<Element>  &source);

template <class Element> Bool_t AreCompatible(const TMatrixTSparse<Element> &m1,const TMatrixTSparse<Element> &m2,Int_t verbose=0);

#endif
