// @(#)root/mathcore:$Id$
// Authors: B. Rabacal   11/2010

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2010 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// implementation file for class TKDTreeBinning
//

#include <algorithm>
#include <limits>
#include <cmath>

#include "TKDTreeBinning.h"

#include "Fit/BinData.h"
#include "TBuffer.h"

ClassImp(TKDTreeBinning);

/**
\class TKDTreeBinning
<- TKDTreeBinning - A class providing multidimensional binning ->

The class implements multidimensional binning by constructing a TKDTree inner structure from the
data which is used as the bins.
The bins are retrieved as two double*, one for the minimum bin edges,
the other as the maximum bin edges. For one dimension one of these is enough to correctly define the bins.
For the multidimensional case both minimum and maximum ones are necessary for the bins to be well defined.
The bin edges of d-dimensional data is a d-tet of the bin's thresholds. For example if d=3 the minimum bin
edges of bin b is of the form of the following array: {xbmin, ybmin, zbmin}.
You also have the possibility to sort the bins by their density.

Details of usage can be found in `$ROOTSYS/tutorials/math/kdTreeBinning.C` and more information on
the embedded TKDTree documentation.

@ingroup MathCore

*/

struct TKDTreeBinning::CompareAsc {
   // Boolean functor whose predicate depends on the bin's density. Used for ascending sort.
   CompareAsc(const TKDTreeBinning* treebins) : bins(treebins) {}
   Bool_t operator()(UInt_t bin1, UInt_t bin2) {
      return bins->GetBinDensity(bin1) < bins->GetBinDensity(bin2);
   }
   const TKDTreeBinning* bins;
};

struct TKDTreeBinning::CompareDesc {
   // Boolean functor whose predicate depends on the bin's density. Used for descending sort.
   CompareDesc(const TKDTreeBinning* treebins) : bins(treebins) {}
   Bool_t operator()(UInt_t bin1, UInt_t bin2) {
      return bins->GetBinDensity(bin1) > bins->GetBinDensity(bin2);
   }
   const TKDTreeBinning* bins;
};

/// Class's constructor taking the size of the data points, dimension, a data array and the number
/// of bins (default = 100). It is reccomended to have the number of bins as an exact divider of
/// the data size.
/// The data array must be organized with a stride=1 for the points and = N (the dataSize) for the dimension.
///
/// Thus data[] = x1,x2,x3,......xN, y1,y2,y3......yN, z1,z2,...........zN,....
///
/// Note that the passed dataSize is not the size of the array but is the number of points (N)
/// The size of the array must be at least  dataDim*dataSize
///
TKDTreeBinning::TKDTreeBinning(UInt_t dataSize, UInt_t dataDim, Double_t* data, UInt_t nBins, bool adjustBinEdges)
: fData(0), fBinMinEdges(std::vector<Double_t>()), fBinMaxEdges(std::vector<Double_t>()), fDataBins((TKDTreeID*)0), fDim(dataDim),
fDataSize(dataSize), fDataThresholds(std::vector<std::pair<Double_t, Double_t> >(fDim, std::make_pair(0., 0.))),
fIsSorted(kFALSE), fIsSortedAsc(kFALSE), fBinsContent(std::vector<UInt_t>()) {
   if (adjustBinEdges) SetBit(kAdjustBinEdges);
   if (data) {
      SetData(data);
      SetNBins(nBins);
   } else {
      if (fData.empty())
         this->Warning("TKDTreeBinning", "Data is nil. Nothing is built.");
   }
}
/// Class's constructor taking the size of the data points, dimension, a data vector and the number
/// of bins (default = 100). It is reccomended to have the number of bins as an exact divider of
/// the data size.
/// The data array must be organized with a stride=1 for the points and = N (the dataSize) for the dimension.
///
/// Thus data[] = x1,x2,x3,......xN, y1,y2,y3......yN, z1,z2,...........zN,....
///
/// Note that the passed data vector may contains a larger size, in case extra coordinates are associated but not used
/// in building the kdtree
/// The size of thedata vector must be at least  dataDim*dataSize
///
TKDTreeBinning::TKDTreeBinning(UInt_t dataSize, UInt_t dataDim, const std::vector<double> &data, UInt_t nBins, bool adjustBinEdges)
: fData(0), fBinMinEdges(std::vector<Double_t>()), fBinMaxEdges(std::vector<Double_t>()), fDataBins((TKDTreeID*)0), fNBins (nBins), fDim(dataDim),
fDataSize(dataSize), fDataThresholds(std::vector<std::pair<Double_t, Double_t> >(fDim, std::make_pair(0., 0.))),
fIsSorted(kFALSE), fIsSortedAsc(kFALSE), fBinsContent(std::vector<UInt_t>()) {
   if (adjustBinEdges) SetBit(kAdjustBinEdges);
   if (!data.empty()) {
      SetData(data);
      SetNBins(nBins);
   } else {
      if (fData.empty())
         this->Warning("TKDTreeBinning", "Data is nil. Nothing is built.");
   }
}

/// Default constructor (for I/O)
TKDTreeBinning::TKDTreeBinning() :
//   fData(nullptr),
   fDataBins(nullptr),
   fNBins (0),
   fDim(0),
   fDataSize(0),
   fIsSorted(kFALSE), fIsSortedAsc(kFALSE)
{}

/// Class's destructor
TKDTreeBinning::~TKDTreeBinning() {
   // if (fData)     delete[] fData;
   if (fDataBins) delete   fDataBins;
}

/// Sets binning inner structure
void TKDTreeBinning::SetNBins(UInt_t bins) {
   fNBins = bins;
   if (fDim && fNBins && fDataSize) {
      if (fDataSize / fNBins) {
         Bool_t remainingData = fDataSize % fNBins;
         if (remainingData) {
            fNBins += 1;
            this->Info("SetNBins", "Number of bins is not enough to hold the data. Extra bin added.");
         }
         fDataBins = new TKDTreeID(fDataSize, fDim, fDataSize / (fNBins - remainingData)); // TKDTree input is data size, data dimension and the content size of bins ("bucket" size)
         SetTreeData();
         fDataBins->Build();
         SetBinsEdges();
         SetBinsContent();
      } else {
         fDataBins = (TKDTreeID*)0;
         this->Warning("SetNBins", "Number of bins is bigger than data size. Nothing is built.");
      }
   } else {
      fDataBins = (TKDTreeID*)0;
      if (!fDim)
         this->Warning("SetNBins", "Data dimension is nil. Nothing is built.");
      if (!fNBins)
         this->Warning("SetNBins", "Number of bins is nil. Nothing is built.");
      if (!fDataSize)
         this->Warning("SetNBins", "Data size is nil. Nothing is built.");
   }
}

/// Sorts bins by their density
void TKDTreeBinning::SortBinsByDensity(Bool_t sortAsc) {
   if (fDim == 1) {
      // in one dim they are already sorted (no need to do anything)
      return;
   } else {
      std::vector<UInt_t> indices(fNBins);        // vector for indices (for the inverse transformation)
      for (UInt_t i = 0; i < fNBins; ++i)
         indices[i] = i;
      if (sortAsc) {
         std::sort(indices.begin(), indices.end(), CompareAsc(this));
         fIsSortedAsc = kTRUE;
      } else {
         std::sort(indices.begin(), indices.end(), CompareDesc(this));
         fIsSortedAsc = kFALSE;
      }
      std::vector<Double_t> binMinEdges(fNBins * fDim);
      std::vector<Double_t> binMaxEdges(fNBins * fDim);
      std::vector<UInt_t> binContent(fNBins );    // reajust also content (not needed bbut better in general!)
      fIndices.resize(fNBins);
      for (UInt_t i = 0; i < fNBins; ++i) {
         for (UInt_t j = 0; j < fDim; ++j) {
            binMinEdges[i * fDim + j] = fBinMinEdges[indices[i] * fDim + j];
            binMaxEdges[i * fDim + j] = fBinMaxEdges[indices[i] * fDim + j];
         }
         binContent[i] = fBinsContent[indices[i]];
         fIndices[indices[i]] = i;
      }
      fBinMinEdges.swap(binMinEdges);
      fBinMaxEdges.swap(binMaxEdges);
      fBinsContent.swap(binContent);

      // not needed anymore if readjusting bin content all the time
      // re-adjust content of extra bins if exists
      // since it is different than the others
      // if ( fDataSize % fNBins != 0) {
      //    UInt_t k = 0;
      //    Bool_t found = kFALSE;
      //    while (!found) {
      //       if (indices[k] == fNBins - 1) {
      //          found = kTRUE;
      //          break;
      //       }
      //       ++k;
      //    }
      //    fBinsContent[fNBins - 1] = fDataBins->GetBucketSize();
      //    fBinsContent[k] = fDataSize % fNBins-1;
      // }

      fIsSorted = kTRUE;
    }
}

void TKDTreeBinning::SetData(Double_t* data) {
   // Sets the data and finds minimum and maximum by dimensional coordinate
   fData.resize(fDim*fDataSize);
   auto first = fData.begin();
   for (UInt_t i = 0; i < fDim; ++i) {
      for (UInt_t j = 0; j < fDataSize; ++j) {
         fData[i*fDataSize+j] = data[i * fDataSize + j];
      }
      auto end = first+fDataSize;
      fDataThresholds[i] = std::make_pair(*std::min_element(first, end), *std::max_element(first,end));
      first = end;
   }
}
void TKDTreeBinning::SetData(const std::vector<double>& data) {
   // Sets the data and finds minimum and maximum by dimensional coordinate
   fData = data;
   auto first = fData.begin();
   // find min/max
   for (UInt_t i = 0; i < fDim; ++i) {
      auto end = first+fDataSize;
      fDataThresholds[i] = std::make_pair(*std::min_element(first, end), *std::max_element(first,end));
      first = end;
   }
}

void TKDTreeBinning::SetTreeData() {
   // Sets the data for constructing the kD-tree
   for (UInt_t i = 0; i < fDim; ++i)
      fDataBins->SetData(i, &fData[i*fDataSize]);
}

void TKDTreeBinning::SetBinsContent() {
   // Sets the bins' content
   fBinsContent.resize(fNBins);
   for (UInt_t i = 0; i < fNBins; ++i)
      fBinsContent[i] = fDataBins->GetBucketSize();
   if ( fDataSize % fNBins != 0 )
      fBinsContent[fNBins - 1] = fDataSize % (fNBins-1);
}

void TKDTreeBinning::SetBinsEdges() {
   // Sets the bins' edges
   //Double_t* rawBinEdges = fDataBins->GetBoundaryExact(fDataBins->GetNNodes());
   Double_t* rawBinEdges = fDataBins->GetBoundary(fDataBins->GetNNodes());
   fCheckedBinEdges = std::vector<std::vector<std::pair<Bool_t, Bool_t> > >(fDim, std::vector<std::pair<Bool_t, Bool_t> >(fNBins, std::make_pair(kFALSE, kFALSE)));
   fCommonBinEdges = std::vector<std::map<Double_t, std::vector<UInt_t> > >(fDim, std::map<Double_t, std::vector<UInt_t> >());
   SetCommonBinEdges(rawBinEdges);
   if (TestBit(kAdjustBinEdges) ) {
      ReadjustMinBinEdges(rawBinEdges);
      ReadjustMaxBinEdges(rawBinEdges);
   }
   SetBinMinMaxEdges(rawBinEdges);
   fCommonBinEdges.clear();
   fCheckedBinEdges.clear();
}

void TKDTreeBinning::SetBinMinMaxEdges(Double_t* binEdges) {
   // Sets the bins' minimum and maximum edges
   fBinMinEdges.reserve(fNBins * fDim);
   fBinMaxEdges.reserve(fNBins * fDim);
   for (UInt_t i = 0; i < fNBins; ++i) {
      for (UInt_t j = 0; j < fDim; ++j) {
         fBinMinEdges.push_back(binEdges[(i * fDim + j) * 2]);
         fBinMaxEdges.push_back(binEdges[(i * fDim + j) * 2 + 1]);
      }
   }
}

void TKDTreeBinning::SetCommonBinEdges(Double_t* binEdges) {
   // Sets indexing on the bin edges which have common boundaries
   for (UInt_t i = 0; i < fDim; ++i) {
      for (UInt_t j = 0; j < fNBins; ++j) {
         Double_t binEdge = binEdges[(j * fDim + i) * 2];
         if(fCommonBinEdges[i].find(binEdge) == fCommonBinEdges[i].end()) {
            std::vector<UInt_t> commonBinEdges;
            for (UInt_t k = 0; k < fNBins; ++k) {
               UInt_t minBinEdgePos = (k * fDim + i) * 2;
               if (std::fabs(binEdge - binEdges[minBinEdgePos]) < std::numeric_limits<Double_t>::epsilon())
                  commonBinEdges.push_back(minBinEdgePos);
               UInt_t maxBinEdgePos = ++minBinEdgePos;
               if (std::fabs(binEdge - binEdges[maxBinEdgePos]) < std::numeric_limits<Double_t>::epsilon())
                  commonBinEdges.push_back(maxBinEdgePos);
            }
            fCommonBinEdges[i][binEdge] = commonBinEdges;
         }
      }
   }
}

void TKDTreeBinning::ReadjustMinBinEdges(Double_t* binEdges) {
   // Readjusts the bins' minimum edge by shifting it slightly lower
   // to avoid overlapping with the data
   for (UInt_t i = 0; i < fDim; ++i) {
      for (UInt_t j = 0; j < fNBins; ++j) {
         if (!fCheckedBinEdges[i][j].first) {
            Double_t binEdge = binEdges[(j * fDim + i) * 2];
            Double_t adjustedBinEdge = binEdge;
            double eps = -10*std::numeric_limits<Double_t>::epsilon();
            if (adjustedBinEdge != 0)
               adjustedBinEdge *= (1. + eps);
            else
               adjustedBinEdge += eps;

            for (UInt_t k = 0; k < fCommonBinEdges[i][binEdge].size(); ++k) {
               UInt_t binEdgePos = fCommonBinEdges[i][binEdge][k];
               Bool_t isMinBinEdge = binEdgePos % 2 == 0;
               UInt_t bin = isMinBinEdge ? (binEdgePos / 2 - i) / fDim : ((binEdgePos - 1) / 2 - i) / fDim;
               binEdges[binEdgePos] = adjustedBinEdge;
               if (isMinBinEdge)
                  fCheckedBinEdges[i][bin].first = kTRUE;
               else
                  fCheckedBinEdges[i][bin].second = kTRUE;
            }
         }
      }
   }
}

void TKDTreeBinning::ReadjustMaxBinEdges(Double_t* binEdges) {
   // Readjusts the bins' maximum edge
   // and shift it sligtly higher
   for (UInt_t i = 0; i < fDim; ++i) {
      for (UInt_t j = 0; j < fNBins; ++j) {
         if (!fCheckedBinEdges[i][j].second) {
            Double_t& binEdge = binEdges[(j * fDim + i) * 2 + 1];
            double eps = 10*std::numeric_limits<Double_t>::epsilon();
            if (binEdge != 0)
               binEdge *= (1. + eps);
            else
               binEdge += eps;


         }
      }
   }
}

/// Returns an array with all bins' minimum edges
/// The edges are arranges as xmin_1,ymin_1, xmin_2,ymin_2,....xmin_{nbin},ymin_{nbin}
const Double_t* TKDTreeBinning::GetBinsMinEdges() const {
   if (fDataBins)
      return &fBinMinEdges[0];
   this->Warning("GetBinsMinEdges", "Binning kd-tree is nil. No bin edges retrieved.");
   this->Info("GetBinsMinEdges", "Returning null pointer.");
   return (Double_t*)0;
}

/// Returns an array with all bins' maximum edges
/// The edges are arranges as xmax_1,ymax_1, xmax_2,ymax_2,....xmax_{nbin},ymax_{nbin}
const Double_t* TKDTreeBinning::GetBinsMaxEdges() const {
   // Returns the bins' maximum edges
   if (fDataBins)
      return &fBinMaxEdges[0];
   this->Warning("GetBinsMaxEdges", "Binning kd-tree is nil. No bin edges retrieved.");
   this->Info("GetBinsMaxEdges", "Returning null pointer.");
   return (Double_t*)0;
}

/// Returns a pair of an array with all bins minimum and maximum edges
std::pair<const Double_t*, const Double_t*> TKDTreeBinning::GetBinsEdges() const {
   // Returns the bins' edges
   if (fDataBins)
      return std::make_pair(GetBinsMinEdges(), GetBinsMaxEdges());
   this->Warning("GetBinsEdges", "Binning kd-tree is nil. No bin edges retrieved.");
   this->Info("GetBinsEdges", "Returning null pointer pair.");
   return std::make_pair((Double_t*)0, (Double_t*)0);
}

/// Returns the bin's minimum edges. 'bin' is between 0 and fNBins - 1
const Double_t* TKDTreeBinning::GetBinMinEdges(UInt_t bin) const {
   if (fDataBins)
      if (bin < fNBins)
         return &fBinMinEdges[bin * fDim];
      else
         this->Warning("GetBinMinEdges", "No such bin. 'bin' is between 0 and %d", fNBins - 1);
   else
      this->Warning("GetBinMinEdges", "Binning kd-tree is nil. No bin edges retrieved.");
   this->Info("GetBinMinEdges", "Returning null pointer.");
   return (Double_t*)0;
}

/// Returns the bin's maximum edges. 'bin' is between 0 and fNBins - 1
const Double_t* TKDTreeBinning::GetBinMaxEdges(UInt_t bin) const {
   if (fDataBins)
      if (bin < fNBins)
         return &fBinMaxEdges[bin * fDim];
      else
         this->Warning("GetBinMaxEdges", "No such bin. 'bin' is between 0 and %d", fNBins - 1);
   else
      this->Warning("GetBinMaxEdges", "Binning kd-tree is nil. No bin edges retrieved.");
   this->Info("GetBinMaxEdges", "Returning null pointer.");
   return (Double_t*)0;
}

/// Returns a pir with the bin's edges. 'bin' is between 0 and fNBins - 1
std::pair<const Double_t*, const Double_t*> TKDTreeBinning::GetBinEdges(UInt_t bin) const {
   if (fDataBins)
      if (bin < fNBins)
         return std::make_pair(GetBinMinEdges(bin), GetBinMaxEdges(bin));
      else
         this->Warning("GetBinEdges", "No such bin. 'bin' is between 0 and %d", fNBins - 1);
   else
      this->Warning("GetBinEdges", "Binning kd-tree is nil. No bin edges retrieved.");
   this->Info("GetBinEdges", "Returning null pointer pair.");
   return std::make_pair((Double_t*)0, (Double_t*)0);
}

/// Returns the number of bins
UInt_t TKDTreeBinning::GetNBins() const {
   return fNBins;
}

/// Returns the number of dimensions
UInt_t TKDTreeBinning::GetDim() const {
   return fDim;
}

/// Returns the number of points in bin. 'bin' is between 0 and fNBins - 1
UInt_t TKDTreeBinning::GetBinContent(UInt_t bin) const {
   if(bin <= fNBins - 1)
         return fBinsContent[bin];
   this->Warning("GetBinContent", "No such bin. Returning 0.");
   this->Info("GetBinContent", "'bin' is between 0 and %d.", fNBins - 1);
   return 0;
}


/// Returns the kD-Tree structure of the binning
TKDTreeID* TKDTreeBinning::GetTree() const {
   if (fDataBins)
      return fDataBins;
   this->Warning("GetTree", "Binning kd-tree is nil. No embedded kd-tree retrieved. Returning null pointer.");
   return (TKDTreeID*)0;
}

// Returns the data array in the dim coordinate. 'dim' is between 0 and fDim - 1
const Double_t* TKDTreeBinning::GetDimData(UInt_t dim) const {
   if(dim < fDim)
      return &fData[dim*fDataSize];
   this->Warning("GetDimData", "No such dimensional coordinate. No coordinate data retrieved. Returning null pointer.");
   this->Info("GetDimData", "'dim' is between 0 and %d.", fDim - 1);
   return 0;
}

/// Returns the minimum of the data in the dim coordinate. 'dim' is between 0 and fDim - 1
Double_t TKDTreeBinning::GetDataMin(UInt_t dim) const {
   if(dim < fDim)
      return fDataThresholds[dim].first;
   this->Warning("GetDataMin", "No such dimensional coordinate. No coordinate data minimum retrieved. Returning +inf.");
   this->Info("GetDataMin", "'dim' is between 0 and %d.", fDim - 1);
   return std::numeric_limits<Double_t>::infinity();
}

/// Returns the maximum of the data in the dim coordinate. 'dim' is between 0 and fDim - 1
Double_t TKDTreeBinning::GetDataMax(UInt_t dim) const {
   if(dim < fDim)
      return fDataThresholds[dim].second;
   this->Warning("GetDataMax", "No such dimensional coordinate. No coordinate data maximum retrieved. Returning -inf.");
   this->Info("GetDataMax", "'dim' is between 0 and %d.", fDim - 1);
   return -1 * std::numeric_limits<Double_t>::infinity();
}

/// Returns the density in bin. 'bin' is between 0 and fNBins - 1
/// The density is the bin content/ bin volume
Double_t TKDTreeBinning::GetBinDensity(UInt_t bin) const {
   if(bin < fNBins) {
      Double_t volume = GetBinVolume(bin);
      if (!volume)
         this->Warning("GetBinDensity", "Volume is null. Returning -1.");
      return GetBinContent(bin) / volume;
   }
   this->Warning("GetBinDensity", "No such bin. Returning -1.");
   this->Info("GetBinDensity", "'bin' is between 0 and %d.", fNBins - 1);
   return -1.;
}

/// Returns the (hyper)volume of bin. 'bin' is between 0 and fNBins - 1
Double_t TKDTreeBinning::GetBinVolume(UInt_t bin) const {
   if(bin < fNBins) {
      std::pair<const Double_t*, const Double_t*> binEdges = GetBinEdges(bin);
      Double_t volume = 1.;
      for (UInt_t i = 0; i < fDim; ++i) {
         volume *= (binEdges.second[i] - binEdges.first[i]);
      }
      return volume;
   }
   this->Warning("GetBinVolume", "No such bin. Returning 0.");
   this->Info("GetBinVolume", "'bin' is between 0 and %d.", fNBins - 1);
   return 0.;
}

/// Returns a pointer to the vector of the bin edges for one dimensional binning only.
/// size of the vector is fNBins + 1 is the vector has been sorted in increasing bin edges
/// N.B : if one does not call SortOneDimBinEdges the bins are not ordered
const double * TKDTreeBinning::GetOneDimBinEdges() const  {
   if (fDim == 1) {
      // no need to sort here because vector is already sorted in one dim
      return &fBinMinEdges.front();
   }
   this->Warning("GetOneDimBinEdges", "Data is multidimensional. No sorted bin edges retrieved. Returning null pointer.");
   this->Info("GetOneDimBinEdges", "This method can only be invoked if the data is a one dimensional set");
   return 0;
}

/// Sort the one-dimensional bin edges and retuns a pointer to them
const Double_t * TKDTreeBinning::SortOneDimBinEdges(Bool_t sortAsc) {
   if (fDim != 1) {
      this->Warning("SortOneDimBinEdges", "Data is multidimensional. Cannot sorted bin edges. Returning null pointer.");
      this->Info("SortOneDimBinEdges", "This method can only be invoked if the data is a one dimensional set");
      return 0;
   }
   // order bins by increasing (or decreasing ) x positions
   std::vector<UInt_t> indices(fNBins);
   TMath::Sort( fNBins, &fBinMinEdges[0], &indices[0], !sortAsc );

   std::vector<Double_t> binMinEdges(fNBins );
   std::vector<Double_t> binMaxEdges(fNBins );
   std::vector<UInt_t> binContent(fNBins );    // reajust also content (not needed but better in general!)
   fIndices.resize( fNBins );
   for (UInt_t i = 0; i < fNBins; ++i) {
      binMinEdges[i ] = fBinMinEdges[indices[i] ];
      binMaxEdges[i ] = fBinMaxEdges[indices[i] ];
      binContent[i] = fBinsContent[indices[i] ];
      fIndices[indices[i] ] = i;  // for the inverse transformation
   }
   fBinMinEdges.swap(binMinEdges);
   fBinMaxEdges.swap(binMaxEdges);
   fBinsContent.swap(binContent);

   fIsSorted = kTRUE;

   // Add also the upper(lower) edge to the min (max) list
   if (sortAsc) {
      fBinMinEdges.push_back(fBinMaxEdges.back());
      fIsSortedAsc = kTRUE;
      return &fBinMinEdges[0];
   }
   fBinMaxEdges.push_back(fBinMinEdges.back());
   return &fBinMaxEdges[0];

}

/// Returns the geometric center of of the bin. 'bin' is between 0 and fNBins - 1,
/// Array must be deleted as "delete [] res"
const Double_t* TKDTreeBinning::GetBinCenter(UInt_t bin) const {
   if(bin < fNBins) {
      Double_t* result = new Double_t[fDim];
      std::pair<const Double_t*, const Double_t*> binEdges = GetBinEdges(bin);
      for (UInt_t i = 0; i < fDim; ++i) {
         result[i] = (binEdges.second[i] + binEdges.first[i]) / 2.;
      }
      return result;
   }
   this->Warning("GetBinCenter", "No such bin. Returning null pointer.");
   this->Info("GetBinCenter", "'bin' is between 0 and %d.", fNBins - 1);
   return nullptr;
}

/// Returns a pointer to the vector of the bin widths. 'bin' is between 0 and fNBins - 1
/// Array must be deleted as "delete [] res"
const Double_t* TKDTreeBinning::GetBinWidth(UInt_t bin) const {
   if(bin < fNBins) {
      Double_t* result = new Double_t[fDim];
      std::pair<const Double_t*, const Double_t*> binEdges = GetBinEdges(bin);
      for (UInt_t i = 0; i < fDim; ++i) {
         result[i] = (binEdges.second[i] - binEdges.first[i]);
      }
      return result;
   }
   this->Warning("GetBinWidth", "No such bin. Returning null pointer.");
   this->Info("GetBinWidth", "'bin' is between 0 and %d.", fNBins - 1);
   return nullptr;
}

/// Return the bin with maximum density
UInt_t TKDTreeBinning::GetBinMaxDensity() const {
   if (fIsSorted) {
      if (fIsSortedAsc)
         return fNBins - 1;
      else return 0;
   }
   UInt_t* indices = new UInt_t[fNBins];
   for (UInt_t i = 0; i < fNBins; ++i)
      indices[i] = i;
   UInt_t result = *std::max_element(indices, indices + fNBins, CompareAsc(this));
   delete [] indices;
   return  result;
}

/// Return the bin with minimum density
UInt_t TKDTreeBinning::GetBinMinDensity() const {
   if (fIsSorted) {
      if (!fIsSortedAsc)
         return fNBins - 1;
      else return 0;
   }
   UInt_t* indices = new UInt_t[fNBins];
   for (UInt_t i = 0; i < fNBins; ++i)
      indices[i] = i;
   UInt_t result = *std::min_element(indices, indices + fNBins, CompareAsc(this));
   delete [] indices;
   return result;
}

/// Fill the bin data set (class ROOT::Fit::BinData) with the result of the TKDTree binning
void TKDTreeBinning::FillBinData(ROOT::Fit::BinData & data) const {
   if (!fDataBins) return;
   data.Initialize(fNBins, fDim);
   for (unsigned int i = 0; i < fNBins; ++i) {
      data.Add( GetBinMinEdges(i), GetBinDensity(i), std::sqrt(double(GetBinContent(i) ))/ GetBinVolume(i) );
      data.AddBinUpEdge(GetBinMaxEdges(i) );
   }
}

/// find the corresponding bin index given the coordinate of a point
UInt_t TKDTreeBinning::FindBin(const Double_t * point) const {

   Int_t inode = fDataBins->FindNode(point);
   // find node return the index in the total nodes and the bins are only the terminal ones
   // so we subtract all the non-terminal nodes
   inode -= fDataBins->GetNNodes();
   R__ASSERT( inode >= 0);
   UInt_t bin = inode;

   if (!fIsSorted) return  bin;
   //return std::distance(fIndices.begin(), std::find(fIndices.begin(), fIndices.end(), bin ) );
   return fIndices[bin];
}

/// Return the corresponding point belonging to the bin i
std::vector<std::vector<Double_t> > TKDTreeBinning::GetPointsInBin(UInt_t bin) const {
   std::vector<Double_t> point(fDim);
   std::vector< std::vector<Double_t> > thePoints;
   if (fData.size() == 0) {
      Error("GetPointsInBin","Internal data set is not valid");
      return thePoints;
   }
   if (!fDataBins) {
      Error("GetPointsInBin","Internal TKDTree is not valid");
      return thePoints;
   }
   if (bin >= fNBins) {
      Error("GetPointsInBin","Invalid bin number");
      return thePoints;
   }
   // need to find the bin number including the non-terminal node
   int inode = bin + fDataBins->GetNNodes();
   auto indices = fDataBins->GetPointsIndexes(inode);
   int npoints =  fDataBins->GetNPointsNode(inode);
   thePoints.resize(npoints);
   for (int ipoint = 0; ipoint < npoints; ++ipoint) {
      for (unsigned int idim = 0; idim < fDim;  ++idim) {
         point[idim] = fData[idim*fDataSize+indices[ipoint] ];
      }
      thePoints[ipoint] = point;
   }
   return thePoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.
void TKDTreeBinning::Streamer(TBuffer &b) {
   if (b.IsReading() ) {
      UInt_t R__s, R__c;
      Version_t v = b.ReadVersion(&R__s, &R__c);
      b.ReadClassBuffer(TKDTreeBinning::Class(), this, v, R__s, R__c);
      // we need to rebuild the TKDTree
      if (fDataBins) delete fDataBins;
      SetNBins(fNBins);
   }
   else {
      // case of writing
       b.WriteClassBuffer(TKDTreeBinning::Class(),this);
       //std::cout << "writing npar = " << GetNpar() << std::endl;
   }
}
