// @(#)root/matrix:$Id$
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixTBase
#define ROOT_TMatrixTBase

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrixTBase                                                         //
//                                                                      //
// Template of base class in the linear algebra package                 //
//                                                                      //
//  matrix properties are stored here, however the data storage is part //
//  of the derived classes                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//======================================================================//
// Summary of the streamer version history                              //
//======================================================================//
//              3.10/02      4.00/a   4.00/b   4.00/c 4.00-08 5.05-1    //
// TMatrixFBase   -          2        2        2       4      5         //
// TMatrix        2          3        3        3       3      4         //
// TMatrixF       -          3        3        3       3      4         //
// TMatrixFSym    -          1        1        1       1      2         //
// TMatrixDSparse -          -        -        -       -      2         //
//                                                                      //
// TMatrixDBase   -          2        3        3       4      5         //
// TMatrixD       2          3        3        3       3      4         //
// TMatrixDSym    -          1        1        1       1      2         //
// TMatrixDSparse -          -        1        1       1      2         //
//                                                                      //
// TVector        2          3        3        3       3      4         //
// TVectorF       -          2        2        2       3      4         //
//                                                                      //
// TVectorD       2          2        2        2       3      4         //
//======================================================================//
//                                                                      //
// 4.00/a : (Jan 25 2004) introduced new classes/inheritance scheme,    //
//          TMatrix now inherits from TMatrixF                          //
//                                                                      //
//          TMatrixF::TMatrixFBase                                      //
//          TMatrixFSym::TMatrixFBase                                   //
//          TMatrixD::TMatrixDBase                                      //
//          TMatrixDSym::TMatrixDBase                                   //
//                                                                      //
// 4.00/b : (May 12 2004) introduced TMatrixDSparse and added new       //
//          element fNRowIndex to TMatrixFBase and TMatrixDBase         //
//          TMatrixDSparse::TMatrixDBase                                //
//                                                                      //
// 4.00/c : (May 27 2004) Used the TObject::fBits to store validity     //
//           state for vectors and matrices                             //
//                                                                      //
// 5.05-1 :  templates TMatrixTBase,TMatrixT,TMatrixTSym and            //
//           TMatrixTSparse were introduced, all versions were          //
//           increased by 1 .                                           //
//                                                                      //
//======================================================================//

#include "TError.h"
#include "TObject.h"
#include "TMathBase.h"
#include "TMatrixFBasefwd.h"
#include "TMatrixDBasefwd.h"
#include "TVectorFfwd.h"
#include "TVectorDfwd.h"

#include <limits>

template<class Element> class TVectorT;
template<class Element> class TElementActionT;
template<class Element> class TElementPosActionT;

R__EXTERN Int_t gMatrixCheck;

template<class Element> class TMatrixTBase : public TObject {

private:
   Element *GetElements();  // This function is now obsolete (and is not implemented) you should use TMatrix::GetMatrixArray().

protected:
   Int_t    fNrows;               // number of rows
   Int_t    fNcols;               // number of columns
   Int_t    fRowLwb;              // lower bound of the row index
   Int_t    fColLwb;              // lower bound of the col index
   Int_t    fNelems;              // number of elements in matrix
   Int_t    fNrowIndex;           // length of row index array (= fNrows+1) wich is only used for sparse matrices

   Element  fTol;                 // sqrt(epsilon); epsilon is smallest number number so that  1+epsilon > 1
                                  //  fTol is used in matrix decomposition (like in inversion)

   Bool_t   fIsOwner;             //!default kTRUE, when Use array kFALSE

   static  void DoubleLexSort (Int_t n,Int_t *first,Int_t *second,Element *data);
   static  void IndexedLexSort(Int_t n,Int_t *first,Int_t swapFirst,
                               Int_t *second,Int_t swapSecond,Int_t *index);

   enum {kSizeMax = 25};          // size data container on stack, see New_m(),Delete_m()
   enum {kWorkMax = 100};         // size of work array's in several routines

   enum EMatrixStatusBits {
     kStatus = BIT(14) // set if matrix object is valid
   };

public:

   TMatrixTBase():
     fNrows(0), fNcols(0), fRowLwb(0), fColLwb(0), fNelems(0), fNrowIndex(0),
     fTol(0), fIsOwner(kTRUE) { }

   virtual ~TMatrixTBase() {}

           inline       Int_t     GetRowLwb     () const { return fRowLwb; }
           inline       Int_t     GetRowUpb     () const { return fNrows+fRowLwb-1; }
           inline       Int_t     GetNrows      () const { return fNrows; }
           inline       Int_t     GetColLwb     () const { return fColLwb; }
           inline       Int_t     GetColUpb     () const { return fNcols+fColLwb-1; }
           inline       Int_t     GetNcols      () const { return fNcols; }
           inline       Int_t     GetNoElements () const { return fNelems; }
           inline       Element   GetTol        () const { return fTol; }

   virtual        const Element  *GetMatrixArray  () const = 0;
   virtual              Element  *GetMatrixArray  ()       = 0;
   virtual        const Int_t    *GetRowIndexArray() const = 0;
   virtual              Int_t    *GetRowIndexArray()       = 0;
   virtual        const Int_t    *GetColIndexArray() const = 0;
   virtual              Int_t    *GetColIndexArray()       = 0;

   virtual              TMatrixTBase<Element> &SetRowIndexArray(Int_t *data) = 0;
   virtual              TMatrixTBase<Element> &SetColIndexArray(Int_t *data) = 0;
   virtual              TMatrixTBase<Element> &SetMatrixArray  (const Element *data,Option_t *option="");
           inline       Element                SetTol          (Element tol);

   virtual void   Clear      (Option_t *option="") = 0;

   inline  void   Invalidate ()       { SetBit(kStatus); }
   inline  void   MakeValid  ()       { ResetBit(kStatus); }
   inline  Bool_t IsValid    () const { return !TestBit(kStatus); }
   inline  Bool_t IsOwner    () const { return fIsOwner; }
   virtual Bool_t IsSymmetric() const;

   virtual TMatrixTBase<Element> &GetSub(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,
                                         TMatrixTBase<Element> &target,Option_t *option="S") const = 0;
   virtual TMatrixTBase<Element> &SetSub(Int_t row_lwb,Int_t col_lwb,const TMatrixTBase<Element> &source) = 0;

   virtual void                   GetMatrix2Array(Element *data,Option_t *option="") const;
   virtual TMatrixTBase<Element> &InsertRow      (Int_t row,Int_t col,const Element *v,Int_t n = -1);
   virtual void                   ExtractRow     (Int_t row,Int_t col,      Element *v,Int_t n = -1) const;

   virtual TMatrixTBase<Element> &Shift          (Int_t row_shift,Int_t col_shift);
   virtual TMatrixTBase<Element> &ResizeTo       (Int_t nrows,Int_t ncols,Int_t nr_nonzeros=-1) = 0;
   virtual TMatrixTBase<Element> &ResizeTo       (Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,Int_t nr_nonzeros=-1) = 0;

   virtual Double_t Determinant() const                          { AbstractMethod("Determinant()"); return 0.; }
   virtual void     Determinant(Double_t &d1,Double_t &d2) const { AbstractMethod("Determinant()"); d1 = 0.; d2 = 0.; }

   virtual TMatrixTBase<Element> &Zero       ();
   virtual TMatrixTBase<Element> &Abs        ();
   virtual TMatrixTBase<Element> &Sqr        ();
   virtual TMatrixTBase<Element> &Sqrt       ();
   virtual TMatrixTBase<Element> &UnitMatrix ();

   virtual TMatrixTBase<Element> &NormByDiag (const TVectorT<Element> &v,Option_t *option="D");

   virtual Element RowNorm    () const;
   virtual Element ColNorm    () const;
   virtual Element E2Norm     () const;
   inline  Element NormInf    () const { return RowNorm(); }
   inline  Element Norm1      () const { return ColNorm(); }
   virtual Int_t   NonZeros   () const;
   virtual Element Sum        () const;
   virtual Element Min        () const;
   virtual Element Max        () const;

   void Draw (Option_t *option="");       // *MENU*
   void Print(Option_t *name  ="") const; // *MENU*

   virtual Element   operator()(Int_t rown,Int_t coln) const = 0;
   virtual Element  &operator()(Int_t rown,Int_t coln)       = 0;

   Bool_t operator==(Element val) const;
   Bool_t operator!=(Element val) const;
   Bool_t operator< (Element val) const;
   Bool_t operator<=(Element val) const;
   Bool_t operator> (Element val) const;
   Bool_t operator>=(Element val) const;

   virtual TMatrixTBase<Element> &Apply(const TElementActionT<Element>    &action);
   virtual TMatrixTBase<Element> &Apply(const TElementPosActionT<Element> &action);

   virtual TMatrixTBase<Element> &Randomize(Element alpha,Element beta,Double_t &seed);

   // make it public since it can be called by TMatrixTRow
   static Element & NaNValue();

   ClassDef(TMatrixTBase,5) // Matrix base class (template)
};

#ifndef __CLING__
// When building with -fmodules, it instantiates all pending instantiations,
// instead of delaying them until the end of the translation unit.
// We 'got away with' probably because the use and the definition of the
// explicit specialization do not occur in the same TU.
//
// In case we are building with -fmodules, we need to forward declare the
// specialization in order to compile the dictionary G__Matrix.cxx.
template <> TClass *TMatrixTBase<double>::Class();
#endif // __CLING__


template<class Element> Element TMatrixTBase<Element>::SetTol(Element newTol)
{
   const Element  oldTol = fTol;
   if (newTol >= 0.0)
      fTol = newTol;
   return oldTol;
}

template<class Element> Bool_t  operator==   (const TMatrixTBase<Element>  &m1,const TMatrixTBase<Element>  &m2);
template<class Element> Element E2Norm       (const TMatrixTBase<Element>  &m1,const TMatrixTBase<Element>  &m2);
template<class Element1,class Element2>
                        Bool_t  AreCompatible(const TMatrixTBase<Element1> &m1,const TMatrixTBase<Element2> &m2,Int_t verbose=0);
template<class Element> void    Compare      (const TMatrixTBase<Element>  &m1,const TMatrixTBase<Element>  &m2);

// Service functions (useful in the verification code).
// They print some detail info if the validation condition fails

template<class Element> Bool_t VerifyMatrixValue   (const TMatrixTBase<Element> &m,Element val,
                                                    Int_t verbose,Element maxDevAllow);
template<class Element> Bool_t VerifyMatrixValue   (const TMatrixTBase<Element> &m,Element val,Int_t verbose)
                                                                           { return VerifyMatrixValue(m,val,verbose,Element(0.)); }
template<class Element> Bool_t VerifyMatrixValue   (const TMatrixTBase<Element> &m,Element val)
                                                                           { return VerifyMatrixValue(m,val,1,Element(0.)); }
template<class Element> Bool_t VerifyMatrixIdentity(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2,
                                                    Int_t verbose,Element maxDevAllow);
template<class Element> Bool_t VerifyMatrixIdentity(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2,Int_t verbose)
                                                                           { return VerifyMatrixIdentity(m1,m2,verbose,Element(0.)); }
template<class Element> Bool_t VerifyMatrixIdentity(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2)
                                                                           { return VerifyMatrixIdentity(m1,m2,1,Element(0.)); }

#endif
