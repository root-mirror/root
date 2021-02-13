// @(#)root/matrix:$Id: 2d00df45ce4c38c7ea0930d6b520cbf4cfb9152e $
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TMatrixTBase
    \ingroup Matrix

    TMatrixTBase

    Template of base class in the linear algebra package.

    See \ref MatrixPage for the documentation of the linear algebra package.

    Matrix properties are stored here, however the data storage is part
    of the derived classes
*/

#include "TMatrixTBase.h"
#include "TVectorT.h"
#include "TBuffer.h"
#include "TROOT.h"
#include "TMath.h"
#include "snprintf.h"
#include <climits>

Int_t gMatrixCheck = 1;

templateClassImp(TMatrixTBase);


////////////////////////////////////////////////////////////////////////////////
/// Lexical sort on array data using indices first and second

template<class Element>
void TMatrixTBase<Element>::DoubleLexSort(Int_t n,Int_t *first,Int_t *second,Element *data)
{
   const int incs[] = {1,5,19,41,109,209,505,929,2161,3905,8929,16001,INT_MAX};

   Int_t kinc = 0;
   while (incs[kinc] <= n/2)
      kinc++;
   kinc -= 1;

   // incs[kinc] is the greatest value in the sequence that is also <= n/2.
   // If n == {0,1}, kinc == -1 and so no sort will take place.

   for( ; kinc >= 0; kinc--) {
      const Int_t inc = incs[kinc];

      for (Int_t k = inc; k < n; k++) {
         const Element tmp = data[k];
         const Int_t fi = first [k];
         const Int_t se = second[k];
         Int_t j;
         for (j = k; j >= inc; j -= inc) {
            if ( fi < first[j-inc] || (fi == first[j-inc] && se < second[j-inc]) ) {
               data  [j] = data  [j-inc];
               first [j] = first [j-inc];
               second[j] = second[j-inc];
            } else
               break;
         }
         data  [j] = tmp;
         first [j] = fi;
         second[j] = se;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Lexical sort on array data using indices first and second

template<class Element>
void TMatrixTBase<Element>::IndexedLexSort(Int_t n,Int_t *first,Int_t swapFirst,
                                           Int_t *second,Int_t swapSecond,Int_t *index)
{
   const int incs[] = {1,5,19,41,109,209,505,929,2161,3905,8929,16001,INT_MAX};

   Int_t kinc = 0;
   while (incs[kinc] <= n/2)
      kinc++;
   kinc -= 1;

   // incs[kinc] is the greatest value in the sequence that is also less
   // than n/2.

   for( ; kinc >= 0; kinc--) {
      const Int_t inc = incs[kinc];

      if ( !swapFirst && !swapSecond ) {
         for (Int_t k = inc; k < n; k++) {
            // loop over all subarrays defined by the current increment
            const Int_t ktemp = index[k];
            const Int_t fi = first [ktemp];
            const Int_t se = second[ktemp];
            // Insert element k into the sorted subarray
            Int_t j;
            for (j = k; j >= inc; j -= inc) {
               // Loop over the elements in the current subarray
               if (fi < first[index[j-inc]] || (fi == first[index[j-inc]] && se < second[index[j-inc]])) {
                  // Swap elements j and j - inc, implicitly use the fact
                  // that ktemp hold element j to avoid having to assign to
                  // element j-inc
                    index[j] = index[j-inc];
               } else {
                  // There are no more elements in this sorted subarray which
                  // are less than element j
                  break;
               }
            } // End loop over the elements in the current subarray
            // Move index[j] out of temporary storage
            index[j] = ktemp;
            // The element has been inserted into the subarray.
         } // End loop over all subarrays defined by the current increment
      } else if ( swapSecond && !swapFirst ) {
         for (Int_t k = inc; k < n; k++) {
            const Int_t ktemp = index[k];
            const Int_t fi = first [ktemp];
            const Int_t se = second[k];
            Int_t j;
            for (j = k; j >= inc; j -= inc) {
               if (fi < first[index[j-inc]] || (fi == first[index[j-inc]] && se < second[j-inc])) {
                  index [j] = index[j-inc];
                  second[j] = second[j-inc];
               } else {
                  break;
               }
            }
            index[j]  = ktemp;
            second[j] = se;
         }
      } else if (swapFirst  && !swapSecond) {
         for (Int_t k = inc; k < n; k++ ) {
            const Int_t ktemp = index[k];
            const Int_t fi = first[k];
            const Int_t se = second[ktemp];
            Int_t j;
            for (j = k; j >= inc; j -= inc) {
               if ( fi < first[j-inc] || (fi == first[j-inc] && se < second[ index[j-inc]])) {
                  index[j] = index[j-inc];
                  first[j] = first[j-inc];
               } else {
                  break;
               }
            }
            index[j] = ktemp;
            first[j] = fi;
         }
      } else { // Swap both
         for (Int_t k = inc; k < n; k++ ) {
            const Int_t ktemp = index[k];
            const Int_t fi = first [k];
            const Int_t se = second[k];
            Int_t j;
            for (j = k; j >= inc; j -= inc) {
               if ( fi < first[j-inc] || (fi == first[j-inc] && se < second[j-inc])) {
                  index [j] = index [j-inc];
                  first [j] = first [j-inc];
                  second[j] = second[j-inc];
               } else {
                  break;
               }
            }
            index[j]  = ktemp;
            first[j]  = fi;
            second[j] = se;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy array data to matrix . It is assumed that array is of size >= fNelems
/// (=)))) fNrows*fNcols
/// option indicates how the data is stored in the array:
/// option =
///         - 'F'   : column major (Fortran) m[i][j] = array[i+j*fNrows]
///         - else  : row major    (C)       m[i][j] = array[i*fNcols+j] (default)

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::SetMatrixArray(const Element *data,Option_t *option)
{
   R__ASSERT(IsValid());

   TString opt = option;
   opt.ToUpper();

   Element *elem = GetMatrixArray();
   if (opt.Contains("F")) {
      for (Int_t irow = 0; irow < fNrows; irow++) {
         const Int_t off1 = irow*fNcols;
         Int_t off2 = 0;
         for (Int_t icol = 0; icol < fNcols; icol++) {
            elem[off1+icol] = data[off2+irow];
            off2 += fNrows;
         }
      }
   }
   else
      memcpy(elem,data,fNelems*sizeof(Element));

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether matrix is symmetric

template<class Element>
Bool_t TMatrixTBase<Element>::IsSymmetric() const
{
  R__ASSERT(IsValid());

   if ((fNrows != fNcols) || (fRowLwb != fColLwb))
      return kFALSE;

   const Element * const elem = GetMatrixArray();
   for (Int_t irow = 0; irow < fNrows; irow++) {
      const Int_t rowOff = irow*fNcols;
      Int_t colOff = 0;
      for (Int_t icol = 0; icol < irow; icol++) {
         if (elem[rowOff+icol] != elem[colOff+irow])
           return kFALSE;
         colOff += fNrows;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy matrix data to array . It is assumed that array is of size >= fNelems
/// (=)))) fNrows*fNcols
/// option indicates how the data is stored in the array:
/// option =
///         - 'F'   : column major (Fortran) array[i+j*fNrows] = m[i][j]
///         - else  : row major    (C)       array[i*fNcols+j] = m[i][j] (default)

template<class Element>
void TMatrixTBase<Element>::GetMatrix2Array(Element *data,Option_t *option) const
{
   R__ASSERT(IsValid());

   TString opt = option;
   opt.ToUpper();

   const Element * const elem = GetMatrixArray();
   if (opt.Contains("F")) {
      for (Int_t irow = 0; irow < fNrows; irow++) {
         const Int_t off1 = irow*fNcols;
         Int_t off2 = 0;
         for (Int_t icol = 0; icol < fNcols; icol++) {
            data[off2+irow] = elem[off1+icol];
            off2 += fNrows;
         }
      }
   }
   else
      memcpy(data,elem,fNelems*sizeof(Element));
}

////////////////////////////////////////////////////////////////////////////////
/// Copy n elements from array v to row rown starting at column coln

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::InsertRow(Int_t rown,Int_t coln,const Element *v,Int_t n)
{
   const Int_t arown = rown-fRowLwb;
   const Int_t acoln = coln-fColLwb;
   const Int_t nr = (n > 0) ? n : fNcols;

   if (gMatrixCheck) {
      if (arown >= fNrows || arown < 0) {
         Error("InsertRow","row %d out of matrix range",rown);
         return *this;
      }

      if (acoln >= fNcols || acoln < 0) {
         Error("InsertRow","column %d out of matrix range",coln);
         return *this;
      }

      if (acoln+nr > fNcols || nr < 0) {
         Error("InsertRow","row length %d out of range",nr);
         return *this;
      }
   }

   const Int_t off = arown*fNcols+acoln;
   Element * const elem = GetMatrixArray()+off;
   memcpy(elem,v,nr*sizeof(Element));

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Store in array v, n matrix elements of row rown starting at column coln

template<class Element>
void TMatrixTBase<Element>::ExtractRow(Int_t rown,Int_t coln,Element *v,Int_t n) const
{
   const Int_t arown = rown-fRowLwb;
   const Int_t acoln = coln-fColLwb;
   const Int_t nr = (n > 0) ? n : fNcols;

   if (gMatrixCheck) {
      if (arown >= fNrows || arown < 0) {
         Error("ExtractRow","row %d out of matrix range",rown);
         return;
      }

      if (acoln >= fNcols || acoln < 0) {
         Error("ExtractRow","column %d out of matrix range",coln);
         return;
      }

      if (acoln+n >= fNcols || nr < 0) {
         Error("ExtractRow","row length %d out of range",nr);
         return;
      }
   }

   const Int_t off = arown*fNcols+acoln;
   const Element * const elem = GetMatrixArray()+off;
   memcpy(v,elem,nr*sizeof(Element));
}

////////////////////////////////////////////////////////////////////////////////
/// Shift the row index by adding row_shift and the column index by adding
/// col_shift, respectively. So [rowLwb..rowUpb][colLwb..colUpb] becomes
/// [rowLwb+row_shift..rowUpb+row_shift][colLwb+col_shift..colUpb+col_shift]

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Shift(Int_t row_shift,Int_t col_shift)
{
   fRowLwb += row_shift;
   fColLwb += col_shift;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Set matrix elements to zero

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Zero()
{
   R__ASSERT(IsValid());
   memset(this->GetMatrixArray(),0,fNelems*sizeof(Element));

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Take an absolute value of a matrix, i.e. apply Abs() to each element.

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Abs()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNelems;
   while (ep < fp) {
      *ep = TMath::Abs(*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Square each element of the matrix.

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Sqr()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNelems;
   while (ep < fp) {
      *ep = (*ep) * (*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Take square root of all elements.

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Sqrt()
{
   R__ASSERT(IsValid());

         Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNelems;
   while (ep < fp) {
      *ep = TMath::Sqrt(*ep);
      ep++;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a unit matrix (matrix need not be a square one).

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::UnitMatrix()
{
   R__ASSERT(IsValid());

   Element *ep = this->GetMatrixArray();
   memset(ep,0,fNelems*sizeof(Element));
   for (Int_t i = fRowLwb; i <= fRowLwb+fNrows-1; i++)
      for (Int_t j = fColLwb; j <= fColLwb+fNcols-1; j++)
         *ep++ = (i==j ? 1.0 : 0.0);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// option:
///  - "D"   :  b(i,j) = a(i,j)/sqrt(abs*(v(i)*v(j)))  (default)
///  - else  :  b(i,j) = a(i,j)*sqrt(abs*(v(i)*v(j)))  (default)

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::NormByDiag(const TVectorT<Element> &v,Option_t *option)
{
   R__ASSERT(IsValid());
   R__ASSERT(v.IsValid());

   if (gMatrixCheck) {
      const Int_t nMax = TMath::Max(fNrows,fNcols);
      if (v.GetNoElements() < nMax) {
         Error("NormByDiag","vector shorter than matrix diagonal");
         return *this;
      }
   }

   TString opt(option);
   opt.ToUpper();
   const Int_t divide = (opt.Contains("D")) ? 1 : 0;

   const Element *pV = v.GetMatrixArray();
         Element *mp = this->GetMatrixArray();

   if (divide) {
      for (Int_t irow = 0; irow < fNrows; irow++) {
         if (pV[irow] != 0.0) {
            for (Int_t icol = 0; icol < fNcols; icol++) {
               if (pV[icol] != 0.0) {
                  const Element val = TMath::Sqrt(TMath::Abs(pV[irow]*pV[icol]));
                  *mp++ /= val;
               } else {
                  Error("NormbyDiag","vector element %d is zero",icol);
                  mp++;
               }
            }
         } else {
            Error("NormbyDiag","vector element %d is zero",irow);
            mp += fNcols;
         }
      }
   } else {
      for (Int_t irow = 0; irow < fNrows; irow++) {
         for (Int_t icol = 0; icol < fNcols; icol++) {
            const Element val = TMath::Sqrt(TMath::Abs(pV[irow]*pV[icol]));
            *mp++ *= val;
         }
      }
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Row matrix norm, MAX{ SUM{ |M(i,j)|, over j}, over i}.
/// The norm is induced by the infinity vector norm.

template<class Element>
Element TMatrixTBase<Element>::RowNorm() const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
         Element norm = 0;

   // Scan the matrix row-after-row
   while (ep < fp) {
      Element sum = 0;
      // Scan a row to compute the sum
      for (Int_t j = 0; j < fNcols; j++)
         sum += TMath::Abs(*ep++);
      norm = TMath::Max(norm,sum);
   }

   R__ASSERT(ep == fp);

   return norm;
}

////////////////////////////////////////////////////////////////////////////////
/// Column matrix norm, MAX{ SUM{ |M(i,j)|, over i}, over j}.
/// The norm is induced by the 1 vector norm.

template<class Element>
Element TMatrixTBase<Element>::ColNorm() const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNcols;
         Element norm = 0;

   // Scan the matrix col-after-col
   while (ep < fp) {
      Element sum = 0;
      // Scan a col to compute the sum
      for (Int_t i = 0; i < fNrows; i++,ep += fNcols)
         sum += TMath::Abs(*ep);
      ep -= fNelems-1;         // Point ep to the beginning of the next col
      norm = TMath::Max(norm,sum);
   }

   R__ASSERT(ep == fp);

   return norm;
}

////////////////////////////////////////////////////////////////////////////////
/// Square of the Euclidian norm, SUM{ m(i,j)^2 }.

template<class Element>
Element TMatrixTBase<Element>::E2Norm() const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
         Element sum = 0;

   for ( ; ep < fp; ep++)
      sum += (*ep) * (*ep);

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the number of elements != 0.0

template<class Element>
Int_t TMatrixTBase<Element>::NonZeros() const
{
   R__ASSERT(IsValid());

   Int_t nr_nonzeros = 0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNelems;
   while (ep < fp)
      if (*ep++ != 0.0) nr_nonzeros++;

   return nr_nonzeros;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute sum of elements

template<class Element>
Element TMatrixTBase<Element>::Sum() const
{
   R__ASSERT(IsValid());

   Element sum = 0.0;
   const Element *ep = this->GetMatrixArray();
   const Element * const fp = ep+fNelems;
   while (ep < fp)
      sum += *ep++;

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// return minimum matrix element value

template<class Element>
Element TMatrixTBase<Element>::Min() const
{
   R__ASSERT(IsValid());

   const Element * const ep = this->GetMatrixArray();
   const Int_t index = TMath::LocMin(fNelems,ep);
   return ep[index];
}

////////////////////////////////////////////////////////////////////////////////
/// return maximum vector element value

template<class Element>
Element TMatrixTBase<Element>::Max() const
{
   R__ASSERT(IsValid());

   const Element * const ep = this->GetMatrixArray();
   const Int_t index = TMath::LocMax(fNelems,ep);
   return ep[index];
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this matrix
/// The histogram is named "TMatrixT" by default and no title

template<class Element>
void TMatrixTBase<Element>::Draw(Option_t *option)
{
   gROOT->ProcessLine(Form("THistPainter::PaintSpecialObjects((TObject*)0x%lx,\"%s\");",
                           (ULong_t)this, option));
}

////////////////////////////////////////////////////////////////////////////////
/// Print the matrix as a table of elements.
/// By default the format "%11.4g" is used to print one element.
/// One can specify an alternative format with eg
///  option ="f=  %6.2f  "

template<class Element>
void TMatrixTBase<Element>::Print(Option_t *option) const
{
   if (!IsValid()) {
      Error("Print","Matrix is invalid");
      return;
   }

   //build format
   const char *format = "%11.4g ";
   if (option) {
      const char *f = strstr(option,"f=");
      if (f) format = f+2;
   }
   char topbar[100];
   snprintf(topbar,100,format,123.456789);
   Int_t nch = strlen(topbar)+1;
   if (nch > 18) nch = 18;
   char ftopbar[20];
   for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
   Int_t nk = 1 + Int_t(TMath::Log10(fNcols));
   snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
   Int_t nch2 = strlen(ftopbar);
   for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
   ftopbar[nch] = '|';
   ftopbar[nch+1] = 0;

   printf("\n%dx%d matrix is as follows",fNrows,fNcols);

   Int_t cols_per_sheet = 5;
   if (nch <= 8) cols_per_sheet =10;
   const Int_t ncols  = fNcols;
   const Int_t nrows  = fNrows;
   const Int_t collwb = fColLwb;
   const Int_t rowlwb = fRowLwb;
   nk = 5+nch*TMath::Min(cols_per_sheet,fNcols);
   for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
   topbar[nk] = 0;
   for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
      printf("\n\n     |");
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
         printf(ftopbar,j+collwb-1);
      printf("\n%s\n",topbar);
      if (fNelems <= 0) continue;
      for (Int_t i = 1; i <= nrows; i++) {
         printf("%4d |",i+rowlwb-1);
         for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
            printf(format,(*this)(i+rowlwb-1,j+collwb-1));
         printf("\n");
      }
   }
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Are all matrix elements equal to val?

template<class Element>
Bool_t TMatrixTBase<Element>::operator==(Element val) const
{
   R__ASSERT(IsValid());

   if (val == 0. && fNelems == 0)
      return kTRUE;

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   for (; ep < fp; ep++)
      if (!(*ep == val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all matrix elements not equal to val?

template<class Element>
Bool_t TMatrixTBase<Element>::operator!=(Element val) const
{
   R__ASSERT(IsValid());

   if (val == 0. && fNelems == 0)
      return kFALSE;

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   for (; ep < fp; ep++)
      if (!(*ep != val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all matrix elements < val?

template<class Element>
Bool_t TMatrixTBase<Element>::operator<(Element val) const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   for (; ep < fp; ep++)
      if (!(*ep < val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all matrix elements <= val?

template<class Element>
Bool_t TMatrixTBase<Element>::operator<=(Element val) const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   for (; ep < fp; ep++)
      if (!(*ep <= val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all matrix elements > val?

template<class Element>
Bool_t TMatrixTBase<Element>::operator>(Element val) const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   for (; ep < fp; ep++)
      if (!(*ep > val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Are all matrix elements >= val?

template<class Element>
Bool_t TMatrixTBase<Element>::operator>=(Element val) const
{
   R__ASSERT(IsValid());

   const Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   for (; ep < fp; ep++)
      if (!(*ep >= val))
         return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply action to each matrix element

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Apply(const TElementActionT<Element> &action)
{
   R__ASSERT(IsValid());

   Element *ep = this->GetMatrixArray();
   const Element * const ep_last = ep+fNelems;
   while (ep < ep_last)
      action.Operation(*ep++);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply action to each element of the matrix. To action the location
/// of the current element is passed.

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Apply(const TElementPosActionT<Element> &action)
{
   R__ASSERT(IsValid());

   Element *ep = this->GetMatrixArray();
   for (action.fI = fRowLwb; action.fI < fRowLwb+fNrows; action.fI++)
      for (action.fJ = fColLwb; action.fJ < fColLwb+fNcols; action.fJ++)
         action.Operation(*ep++);

   R__ASSERT(ep == this->GetMatrixArray()+fNelems);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Randomize matrix element values

template<class Element>
TMatrixTBase<Element> &TMatrixTBase<Element>::Randomize(Element alpha,Element beta,Double_t &seed)
{
   R__ASSERT(IsValid());

   const Element scale = beta-alpha;
   const Element shift = alpha/scale;

         Element *       ep = GetMatrixArray();
   const Element * const fp = ep+fNelems;
   while (ep < fp)
      *ep++ = scale*(Drand(seed)+shift);

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Check to see if two matrices are identical.

template<class Element>
Bool_t operator==(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2)
{
   if (!AreCompatible(m1,m2)) return kFALSE;
   return (memcmp(m1.GetMatrixArray(),m2.GetMatrixArray(),
                   m1.GetNoElements()*sizeof(Element)) == 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Square of the Euclidian norm of the difference between two matrices.

template<class Element>
Element E2Norm(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2)
{
   if (gMatrixCheck && !AreCompatible(m1,m2)) {
      ::Error("E2Norm","matrices not compatible");
      return -1.0;
   }

   const Element *        mp1 = m1.GetMatrixArray();
   const Element *        mp2 = m2.GetMatrixArray();
   const Element * const fmp1 = mp1+m1.GetNoElements();

   Element sum = 0.0;
   for (; mp1 < fmp1; mp1++, mp2++)
      sum += (*mp1 - *mp2)*(*mp1 - *mp2);

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Check that matrice sm1 and m2 areboth valid and have identical shapes .

template<class Element1,class Element2>
Bool_t AreCompatible(const TMatrixTBase<Element1> &m1,const TMatrixTBase<Element2> &m2,Int_t verbose)
{
   if (!m1.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "matrix 1 not valid");
      return kFALSE;
   }
   if (!m2.IsValid()) {
      if (verbose)
         ::Error("AreCompatible", "matrix 2 not valid");
      return kFALSE;
   }

   if (m1.GetNrows()  != m2.GetNrows()  || m1.GetNcols()  != m2.GetNcols() ||
       m1.GetRowLwb() != m2.GetRowLwb() || m1.GetColLwb() != m2.GetColLwb()) {
      if (verbose)
         ::Error("AreCompatible", "matrices 1 and 2 not compatible");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compare two matrices and print out the result of the comparison.

template<class Element>
void Compare(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2)
{
   if (!AreCompatible(m1,m2)) {
      Error("Compare(const TMatrixTBase<Element> &,const TMatrixTBase<Element> &)","matrices are incompatible");
      return;
   }

   printf("\n\nComparison of two TMatrices:\n");

   Element norm1  = 0;      // Norm of the Matrices
   Element norm2  = 0;
   Element ndiff  = 0;      // Norm of the difference
   Int_t   imax   = 0;      // For the elements that differ most
   Int_t   jmax   = 0;
   Element difmax = -1;

   for (Int_t i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++) {
      for (Int_t j = m1.GetColLwb(); j < m1.GetColUpb(); j++) {
         const Element mv1 = m1(i,j);
         const Element mv2 = m2(i,j);
         const Element diff = TMath::Abs(mv1-mv2);

         if (diff > difmax) {
            difmax = diff;
            imax = i;
            jmax = j;
         }
         norm1 += TMath::Abs(mv1);
         norm2 += TMath::Abs(mv2);
         ndiff += TMath::Abs(diff);
      }
   }

   printf("\nMaximal discrepancy    \t\t%g", difmax);
   printf("\n   occured at the point\t\t(%d,%d)",imax,jmax);
   const Element mv1 = m1(imax,jmax);
   const Element mv2 = m2(imax,jmax);
   printf("\n Matrix 1 element is    \t\t%g", mv1);
   printf("\n Matrix 2 element is    \t\t%g", mv2);
   printf("\n Absolute error v2[i]-v1[i]\t\t%g", mv2-mv1);
   printf("\n Relative error\t\t\t\t%g\n",
          (mv2-mv1)/TMath::Max(TMath::Abs(mv2+mv1)/2,(Element)1e-7));

   printf("\n||Matrix 1||   \t\t\t%g", norm1);
   printf("\n||Matrix 2||   \t\t\t%g", norm2);
   printf("\n||Matrix1-Matrix2||\t\t\t\t%g", ndiff);
   printf("\n||Matrix1-Matrix2||/sqrt(||Matrix1|| ||Matrix2||)\t%g\n\n",
          ndiff/TMath::Max(TMath::Sqrt(norm1*norm2),1e-7));
}

////////////////////////////////////////////////////////////////////////////////
/// Validate that all elements of matrix have value val within maxDevAllow.

template<class Element>
Bool_t VerifyMatrixValue(const TMatrixTBase<Element> &m,Element val,Int_t verbose,Element maxDevAllow)
{
   R__ASSERT(m.IsValid());

   if (m == 0)
      return kTRUE;

   Int_t   imax      = 0;
   Int_t   jmax      = 0;
   Element maxDevObs = 0;

   if (TMath::Abs(maxDevAllow) <= 0.0)
      maxDevAllow = std::numeric_limits<Element>::epsilon();

   for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
         const Element dev = TMath::Abs(m(i,j)-val);
         if (dev > maxDevObs) {
            imax    = i;
            jmax    = j;
            maxDevObs = dev;
         }
      }
   }

   if (maxDevObs == 0)
      return kTRUE;

   if (verbose) {
      printf("Largest dev for (%d,%d); dev = |%g - %g| = %g\n",imax,jmax,m(imax,jmax),val,maxDevObs);
      if(maxDevObs > maxDevAllow)
         Error("VerifyElementValue","Deviation > %g\n",maxDevAllow);
   }

   if(maxDevObs > maxDevAllow)
      return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Verify that elements of the two matrices are equal within MaxDevAllow .

template<class Element>
Bool_t VerifyMatrixIdentity(const TMatrixTBase<Element> &m1,const TMatrixTBase<Element> &m2,Int_t verbose,
                            Element maxDevAllow)
{
   if (!AreCompatible(m1,m2,verbose))
      return kFALSE;

   if (m1 == 0 && m2 == 0)
      return kTRUE;

   Int_t   imax      = 0;
   Int_t   jmax      = 0;
   Element maxDevObs = 0;

   if (TMath::Abs(maxDevAllow) <= 0.0)
      maxDevAllow = std::numeric_limits<Element>::epsilon();

   for (Int_t i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++) {
      for (Int_t j = m1.GetColLwb(); j <= m1.GetColUpb(); j++) {
         const Element dev = TMath::Abs(m1(i,j)-m2(i,j));
         if (dev > maxDevObs) {
            imax = i;
            jmax = j;
            maxDevObs = dev;
         }
      }
   }

   if (maxDevObs == 0)
      return kTRUE;

   if (verbose) {
      printf("Largest dev for (%d,%d); dev = |%g - %g| = %g\n",
              imax,jmax,m1(imax,jmax),m2(imax,jmax),maxDevObs);
      if (maxDevObs > maxDevAllow)
         Error("VerifyMatrixValue","Deviation > %g\n",maxDevAllow);
   }

   if (maxDevObs > maxDevAllow)
      return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TMatrixTBase<Element>.

template<class Element>
void TMatrixTBase<Element>::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TMatrixTBase<Element>::Class(),this,R__v,R__s,R__c);
      } else {
         Error("TMatrixTBase<Element>::Streamer","Unknown version number: %d",R__v);
         R__ASSERT(0);
      }
      if (R__v < 4) MakeValid();
   } else {
      R__b.WriteClassBuffer(TMatrixTBase<Element>::Class(),this);
   }
}

// trick to return a reference to nan in operator(i,j_ when i,j are outside of range
template<class Element>
struct nan_value_t {
   static Element gNanValue;
};
template<>
Double_t nan_value_t<Double_t>::gNanValue = std::numeric_limits<Double_t>::quiet_NaN();
template<>
Float_t nan_value_t<Float_t>::gNanValue = std::numeric_limits<Float_t>::quiet_NaN();

template<class Element>
Element & TMatrixTBase<Element>::NaNValue()
{
   return nan_value_t<Element>::gNanValue;
}


template class TMatrixTBase<Float_t>;

template Bool_t   operator==          <Float_t>(const TMatrixFBase &m1,const TMatrixFBase &m2);
template Float_t  E2Norm              <Float_t>(const TMatrixFBase &m1,const TMatrixFBase &m2);
template Bool_t   AreCompatible<Float_t,Float_t>
                                               (const TMatrixFBase &m1,const TMatrixFBase &m2,Int_t verbose);
template Bool_t   AreCompatible<Float_t,Double_t>
                                               (const TMatrixFBase &m1,const TMatrixDBase &m2,Int_t verbose);
template void     Compare             <Float_t>(const TMatrixFBase &m1,const TMatrixFBase &m2);
template Bool_t   VerifyMatrixValue   <Float_t>(const TMatrixFBase &m,Float_t val,Int_t verbose,Float_t maxDevAllow);
template Bool_t   VerifyMatrixValue   <Float_t>(const TMatrixFBase &m,Float_t val);
template Bool_t   VerifyMatrixIdentity<Float_t>(const TMatrixFBase &m1,const TMatrixFBase &m2,
                                                Int_t verbose,Float_t maxDevAllowN);
template Bool_t   VerifyMatrixIdentity<Float_t>(const TMatrixFBase &m1,const TMatrixFBase &m2);

template class TMatrixTBase<Double_t>;

template Bool_t   operator==          <Double_t>(const TMatrixDBase &m1,const TMatrixDBase &m2);
template Double_t E2Norm              <Double_t>(const TMatrixDBase &m1,const TMatrixDBase &m2);
template Bool_t   AreCompatible<Double_t,Double_t>
                                               (const TMatrixDBase &m1,const TMatrixDBase &m2,Int_t verbose);
template Bool_t   AreCompatible<Double_t,Float_t>
                                               (const TMatrixDBase &m1,const TMatrixFBase &m2,Int_t verbose);
template void     Compare             <Double_t>(const TMatrixDBase &m1,const TMatrixDBase &m2);
template Bool_t   VerifyMatrixValue   <Double_t>(const TMatrixDBase &m,Double_t val,Int_t verbose,Double_t maxDevAllow);
template Bool_t   VerifyMatrixValue   <Double_t>(const TMatrixDBase &m,Double_t val);
template Bool_t   VerifyMatrixIdentity<Double_t>(const TMatrixDBase &m1,const TMatrixDBase &m2,
                                                 Int_t verbose,Double_t maxDevAllow);
template Bool_t   VerifyMatrixIdentity<Double_t>(const TMatrixDBase &m1,const TMatrixDBase &m2);
