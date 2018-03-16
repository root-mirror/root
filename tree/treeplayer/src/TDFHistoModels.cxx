// Author: Enrico Guiraud, Danilo Piparo CERN  09/2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/TDFHistoModels.hxx>
#include <ROOT/TSeq.hxx>
#include <TProfile.h>
#include <TProfile2D.h>
#include <stddef.h>
#include <vector>

#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

/**
* \class ROOT::Experimental::TDF::TH1DModel
* \ingroup dataframe
* \brief A struct which stores the parameters of a TH1D
*
* \class ROOT::Experimental::TDF::TH2DModel
* \ingroup dataframe
* \brief A struct which stores the parameters of a TH2D
*
* \class ROOT::Experimental::TDF::TH3DModel
* \ingroup dataframe
* \brief A struct which stores the parameters of a TH3D
*
* \class ROOT::Experimental::TDF::TProfile1DModel
* \ingroup dataframe
* \brief A struct which stores the parameters of a TProfile
*
* \class ROOT::Experimental::TDF::TProfile2DModel
* \ingroup dataframe
* \brief A struct which stores the parameters of a TProfile2D
*/

template <typename T>
inline void FillVector(std::vector<double> &v, int size, T *a)
{
   v.reserve(size);
   for (auto i : ROOT::TSeq<int>(size + 1))
      v.push_back(a[i]);
}

template <>
inline void FillVector<double>(std::vector<double> &v, int size, double *a)
{
   v.assign(a, a + (size_t)(size + 1));
}

inline void SetAxisProperties(const TAxis *axis, double &low, double &up, std::vector<double> &edges)
{
   // Check if this histo has fixed binning
   // Same technique of "Int_t TAxis::FindBin(Double_t)"
   if (!axis->GetXbins()->fN) {
      low = axis->GetXmin();
      up = axis->GetXmax();
   } else {
      // This histo has variable binning
      const auto size = axis->GetNbins() + 1;
      edges.reserve(size);
      for (auto i : ROOT::TSeq<int>(1, size))
         edges.push_back(axis->GetBinLowEdge(i));
      edges.push_back(axis->GetBinUpEdge(size - 1));
   }
}

namespace ROOT {
namespace Experimental {
namespace TDF {

TH1DModel::TH1DModel(const ::TH1D &h) : fName(h.GetName()), fTitle(h.GetTitle()), fNbinsX(h.GetNbinsX())
{
   SetAxisProperties(h.GetXaxis(), fXLow, fXUp, fBinXEdges);
}
TH1DModel::TH1DModel(const char *name, const char *title, int nbinsx, double xlow, double xup)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup)
{
}
TH1DModel::TH1DModel(const char *name, const char *title, int nbinsx, const float *xbins)
   : fName(name), fTitle(title), fNbinsX(nbinsx)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}
TH1DModel::TH1DModel(const char *name, const char *title, int nbinsx, const double *xbins)
   : fName(name), fTitle(title), fNbinsX(nbinsx)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}
std::shared_ptr<::TH1D> TH1DModel::GetHistogram() const
{
   if (fBinXEdges.empty()) {
      return std::make_shared<::TH1D>(fName, fTitle, fNbinsX, fXLow, fXUp);
   } else {
      return std::make_shared<::TH1D>(fName, fTitle, fNbinsX, fBinXEdges.data());
   }
}
TH1DModel::~TH1DModel()
{
}

TH2DModel::TH2DModel(const ::TH2D &h)
   : fName(h.GetName()), fTitle(h.GetTitle()), fNbinsX(h.GetNbinsX()), fNbinsY(h.GetNbinsY())
{
   SetAxisProperties(h.GetXaxis(), fXLow, fXUp, fBinXEdges);
   SetAxisProperties(h.GetYaxis(), fYLow, fYUp, fBinYEdges);
}
TH2DModel::TH2DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, int nbinsy, double ylow,
                     double yup)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fNbinsY(nbinsy), fYLow(ylow), fYUp(yup)
{
}
TH2DModel::TH2DModel(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, double ylow,
                     double yup)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy), fYLow(ylow), fYUp(yup)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}
TH2DModel::TH2DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, int nbinsy,
                     const double *ybins)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fNbinsY(nbinsy)
{
   FillVector(fBinYEdges, nbinsy, ybins);
}
TH2DModel::TH2DModel(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy,
                     const double *ybins)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy)
{
   FillVector(fBinXEdges, nbinsx, xbins);
   FillVector(fBinYEdges, nbinsy, ybins);
}
TH2DModel::TH2DModel(const char *name, const char *title, int nbinsx, const float *xbins, int nbinsy,
                     const float *ybins)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy)
{
   FillVector(fBinXEdges, nbinsx, xbins);
   FillVector(fBinYEdges, nbinsy, ybins);
}
std::shared_ptr<::TH2D> TH2DModel::GetHistogram() const
{
   auto xEdgesEmpty = fBinXEdges.empty();
   auto yEdgesEmpty = fBinYEdges.empty();
   if (xEdgesEmpty && yEdgesEmpty) {
      return std::make_shared<::TH2D>(fName, fTitle, fNbinsX, fXLow, fXUp, fNbinsY, fYLow, fYUp);
   } else if (!xEdgesEmpty && yEdgesEmpty) {
      return std::make_shared<::TH2D>(fName, fTitle, fNbinsX, fBinXEdges.data(), fNbinsY, fYLow, fYUp);
   } else if (xEdgesEmpty && !yEdgesEmpty) {
      return std::make_shared<::TH2D>(fName, fTitle, fNbinsX, fXLow, fXUp, fNbinsY, fBinYEdges.data());
   } else {
      return std::make_shared<::TH2D>(fName, fTitle, fNbinsX, fBinXEdges.data(), fNbinsY, fBinYEdges.data());
   }
}
TH2DModel::~TH2DModel()
{
}

TH3DModel::TH3DModel(const ::TH3D &h)
   : fName(h.GetName()), fTitle(h.GetTitle()), fNbinsX(h.GetNbinsX()), fNbinsY(h.GetNbinsY()), fNbinsZ(h.GetNbinsZ())
{
   SetAxisProperties(h.GetXaxis(), fXLow, fXUp, fBinXEdges);
   SetAxisProperties(h.GetYaxis(), fYLow, fYUp, fBinYEdges);
   SetAxisProperties(h.GetZaxis(), fZLow, fZUp, fBinZEdges);
}
TH3DModel::TH3DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, int nbinsy, double ylow,
                     double yup, int nbinsz, double zlow, double zup)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fNbinsY(nbinsy), fYLow(ylow), fYUp(yup),
     fNbinsZ(nbinsz), fZLow(zlow), fZUp(zup)
{
}
TH3DModel::TH3DModel(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy,
                     const double *ybins, int nbinsz, const double *zbins)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy), fNbinsZ(nbinsz)
{
   FillVector(fBinXEdges, nbinsx, xbins);
   FillVector(fBinYEdges, nbinsy, ybins);
   FillVector(fBinZEdges, nbinsz, zbins);
}
TH3DModel::TH3DModel(const char *name, const char *title, int nbinsx, const float *xbins, int nbinsy,
                     const float *ybins, int nbinsz, const float *zbins)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy), fNbinsZ(nbinsz)
{
   FillVector(fBinXEdges, nbinsx, xbins);
   FillVector(fBinYEdges, nbinsy, ybins);
   FillVector(fBinZEdges, nbinsz, zbins);
}
std::shared_ptr<::TH3D> TH3DModel::GetHistogram() const
{
   if (fBinXEdges.empty() && fBinYEdges.empty() && fBinZEdges.empty()) {
      return std::make_shared<::TH3D>(fName, fTitle, fNbinsX, fXLow, fXUp, fNbinsY, fYLow, fYUp, fNbinsZ, fZLow, fZUp);
   } else {
      return std::make_shared<::TH3D>(fName, fTitle, fNbinsX, fBinXEdges.data(), fNbinsY, fBinYEdges.data(), fNbinsZ,
                                      fBinZEdges.data());
   }
}
TH3DModel::~TH3DModel()
{
}

// Profiles

TProfile1DModel::TProfile1DModel(const ::TProfile &h)
   : fName(h.GetName()), fTitle(h.GetTitle()), fNbinsX(h.GetNbinsX()), fXLow(h.GetXaxis()->GetXmin()),
     fXUp(h.GetXaxis()->GetXmax()), fYLow(h.GetYmin()), fYUp(h.GetYmax()), fOption(h.GetErrorOption())
{
}
TProfile1DModel::TProfile1DModel(const char *name, const char *title, int nbinsx, double xlow, double xup,
                                 const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fOption(option)
{
}

TProfile1DModel::TProfile1DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, double ylow,
                                 double yup, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fYLow(ylow), fYUp(yup), fOption(option)
{
}

TProfile1DModel::TProfile1DModel(const char *name, const char *title, int nbinsx, const float *xbins,
                                 const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fOption(option)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}
TProfile1DModel::TProfile1DModel(const char *name, const char *title, int nbinsx, const double *xbins,
                                 const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fOption(option)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}
TProfile1DModel::TProfile1DModel(const char *name, const char *title, int nbinsx, const double *xbins, double ylow,
                                 double yup, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fYLow(ylow), fYUp(yup), fOption(option)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}
std::shared_ptr<::TProfile> TProfile1DModel::GetProfile() const
{
   return std::make_shared<::TProfile>(fName, fTitle, fNbinsX, fXLow, fXUp, fYLow, fYUp, fOption);
}
TProfile1DModel::~TProfile1DModel()
{
}

TProfile2DModel::TProfile2DModel(const ::TProfile2D &h)
   : fName(h.GetName()), fTitle(h.GetTitle()), fNbinsX(h.GetNbinsX()), fXLow(h.GetXaxis()->GetXmin()),
     fXUp(h.GetXaxis()->GetXmax()), fNbinsY(h.GetNbinsY()), fYLow(h.GetYaxis()->GetXmin()),
     fYUp(h.GetYaxis()->GetXmax()), fZLow(h.GetZmin()), fZUp(h.GetZmax()), fOption(h.GetErrorOption())
{
}
TProfile2DModel::TProfile2DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, int nbinsy,
                                 double ylow, double yup, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fNbinsY(nbinsy), fYLow(ylow), fYUp(yup),
     fOption(option)
{
}

TProfile2DModel::TProfile2DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, int nbinsy,
                                 double ylow, double yup, double zlow, double zup, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fNbinsY(nbinsy), fYLow(ylow), fYUp(yup),
     fZLow(zlow), fZUp(zup), fOption(option)
{
}

TProfile2DModel::TProfile2DModel(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy,
                                 double ylow, double yup, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy), fYLow(ylow), fYUp(yup), fOption(option)
{
   FillVector(fBinXEdges, nbinsx, xbins);
}

TProfile2DModel::TProfile2DModel(const char *name, const char *title, int nbinsx, double xlow, double xup, int nbinsy,
                                 const double *ybins, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fXLow(xlow), fXUp(xup), fNbinsY(nbinsy), fOption(option)
{
   FillVector(fBinYEdges, nbinsy, ybins);
}

TProfile2DModel::TProfile2DModel(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy,
                                 const double *ybins, const char *option)
   : fName(name), fTitle(title), fNbinsX(nbinsx), fNbinsY(nbinsy), fOption(option)
{
   FillVector(fBinXEdges, nbinsx, xbins);
   FillVector(fBinYEdges, nbinsy, ybins);
}

std::shared_ptr<::TProfile2D> TProfile2DModel::GetProfile() const
{
   return std::make_shared<::TProfile2D>(fName, fTitle, fNbinsX, fXLow, fXUp, fNbinsY, fYLow, fYUp, fZLow, fZUp,
                                         fOption);
}

TProfile2DModel::~TProfile2DModel()
{
}

} // ns TDF
} // ns Experimental
} // ns ROOT
