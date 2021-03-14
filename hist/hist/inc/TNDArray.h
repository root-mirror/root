// @(#)root/hist:$Id$
// Author: Axel Naumann, Nov 2011

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TNDArray
#define ROOT_TNDArray

#include "TObject.h"
#include "TError.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TNDArray                                                             //
//                                                                      //
// N-Dim array class.                                                   //
//                                                                      //
// Storage layout:                                                      //
// Assume 3 dimensions, array sizes 2, 4 and 3 i.e. 24 bins:            //
// Data is stored as [0,0,0], [0,0,1], [0,0,2], [0,1,0],...             //
//                                                                      //
// fSizes stores the combined size of each bin in a dimension, i.e. in  //
// above example it would contain 24, 12, 3, 1.                         //
//                                                                      //
// Storage is allocated lazily, only when data is written to the array. //
//                                                                      //
// TNDArrayRef gives access to a sub-dimension, e.g. arr[0][1] in above //
// three-dimensional example, up to an element with conversion operator //
// to double: double value = arr[0][1][2];                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

// Array layout:
// nbins[0] = 2, nbins[1] = 4, nbins[2] = 3 => 24 bins
//
// fSizes: 24, 12, 3 [, 1

class TNDArray: public TObject {
public:
   TNDArray(): fNdimPlusOne(), fSizes() {}

   TNDArray(Int_t ndim, const Int_t* nbins, bool addOverflow = false):
   fNdimPlusOne(), fSizes() {
      TNDArray::Init(ndim, nbins, addOverflow);
   }
   ~TNDArray() {}

   virtual void Init(Int_t ndim, const Int_t* nbins, bool addOverflow = false) {
      // Calculate fSize based on ndim dimensions, nbins for each dimension,
      // possibly adding over- and underflow bin to each dimensions' nbins.
      fNdimPlusOne = ndim + 1;
      fSizes.resize(ndim + 1);
      Int_t overBins = addOverflow ? 2 : 0;
      fSizes[ndim] = 1;
      for (Int_t i = 0; i < ndim; ++i) {
         fSizes[ndim - i - 1] = fSizes[ndim - i] * (nbins[ndim - i - 1] + overBins);
      }
   }

   virtual void Reset(Option_t* option = "") = 0;

   Int_t GetNdimensions() const { return fNdimPlusOne - 1; }
   Long64_t GetNbins() const { return fSizes[0]; }
   Long64_t GetCellSize(Int_t dim) const { return fSizes[dim + 1]; }

   Long64_t GetBin(const Int_t* idx) const {
      // Get the linear bin number for each dimension's bin index
      Long64_t bin = idx[fNdimPlusOne - 2];
      for (Int_t d = 0; d < fNdimPlusOne - 2; ++d) {
         bin += fSizes[d + 1] * idx[d];
      }
      return bin;
   }

   virtual Double_t AtAsDouble(ULong64_t linidx) const = 0;
   virtual void SetAsDouble(ULong64_t linidx, Double_t value) = 0;
   virtual void AddAt(ULong64_t linidx, Double_t value) = 0;

protected:
   Int_t  fNdimPlusOne; // Number of dimensions plus one
   std::vector<Long64_t> fSizes; // bin count
   ClassDef(TNDArray, 1); //Base for n-dimensional array
};

template <typename T>
class TNDArrayRef {
public:
   TNDArrayRef(const T* data, const Long64_t* sizes):
   fData(data), fSizes(sizes) {}

   TNDArrayRef<T> operator[] (Int_t idx) const {
      if (!fData) return TNDArrayRef<T>(0, 0);
      R__ASSERT(idx < fSizes[-1] / fSizes[0] && "index out of range!");
      return TNDArrayRef<T>(fData + idx * fSizes[0], (fSizes[0] == 1) ? 0 : (fSizes + 1));
   }
   operator T() const {
      if (!fData) return T();
      R__ASSERT(fSizes == 0 && "Element operator can only be used on non-array element. Missing an operator[] level?");
      return *fData;
   }

private:
   const T* fData; // pointer into TNDArray's fData
   const Long64_t* fSizes; // pointer into TNDArray's fSizes
   ClassDefNV(TNDArrayRef, 0); // subdimension of a TNDArray
};

template <typename T>
class TNDArrayT: public TNDArray {
public:
   TNDArrayT(): fNumData(), fData() {}

   TNDArrayT(Int_t ndim, const Int_t* nbins, bool addOverflow = false):
   TNDArray(ndim, nbins, addOverflow),
   fNumData(), fData() {
      fNumData = fSizes[0];
   }
   ~TNDArrayT() {}

   void Init(Int_t ndim, const Int_t* nbins, bool addOverflow = false) {
      fData.clear();
      TNDArray::Init(ndim, nbins, addOverflow);
      fNumData = fSizes[0];
   }

   void Reset(Option_t* /*option*/ = "") {
      // Reset the content
      fData.clear();
      fData.resize(fNumData, T());
   }

#ifndef __CINT__
   TNDArrayRef<T> operator[](Int_t idx) const {
      if (!fData) return TNDArrayRef<T>(0, 0);
      R__ASSERT(idx < fSizes[0] / fSizes[1] && "index out of range!");
      return TNDArrayRef<T>(fData.data() + idx * fSizes[1], fSizes.data() + 2);
   }
#endif // __CINT__

   T At(const Int_t* idx) const {
      return At(GetBin(idx));
   }
   T& At(const Int_t* idx) {
      return At(GetBin(idx));
   }
   T At(ULong64_t linidx) const {
      if (fData.empty()) return T();
      return fData[linidx];
   }
   T& At(ULong64_t linidx) {
      if (fData.empty()) fData.resize(fNumData, T());
      return fData[linidx];
   }

   Double_t AtAsDouble(ULong64_t linidx) const {
      if (fData.empty()) return 0.;
      return fData[linidx];
   }
   void SetAsDouble(ULong64_t linidx, Double_t value) {
      if (fData.empty()) fData.resize(fNumData, T());
      fData[linidx] = (T) value;
   }
   void AddAt(ULong64_t linidx, Double_t value) {
      if (fData.empty()) fData.resize(fNumData, T());
      fData[linidx] += (T) value;
   }

protected:
   int fNumData; // number of bins, product of fSizes
   std::vector<T> fData;// data
   ClassDef(TNDArrayT, 1); // N-dimensional array
};

// FIXME: Remove once we implement https://sft.its.cern.ch/jira/browse/ROOT-6284
// When building with -fmodules, it instantiates all pending instantiations,
// instead of delaying them until the end of the translation unit.
// We 'got away with' probably because the use and the definition of the
// explicit specialization do not occur in the same TU.
//
// In case we are building with -fmodules, we need to forward declare the
// specialization in order to compile the dictionary G__Hist.cxx.
template<> void TNDArrayT<double>::Streamer(TBuffer &R__b);
template<> TClass *TNDArrayT<double>::Class();


#endif // ROOT_TNDArray
