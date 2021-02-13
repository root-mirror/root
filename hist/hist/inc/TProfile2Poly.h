// @(#)root/hist:$Id$
// Author: Filip Ilic

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProfile2Poly
#define ROOT_TProfile2Poly

#include "TH2Poly.h"
#include "TProfile.h"
#include <vector>

class TProfile2PolyBin : public TH2PolyBin {
public:
   friend class TProfile2Poly;

   TProfile2PolyBin();
   TProfile2PolyBin(TObject *poly, Int_t bin_number);
   virtual ~TProfile2PolyBin() {}

   void Merge(const TProfile2PolyBin *toMerge);

   void Update();
   void ClearStats();

   Double_t GetEffectiveEntries() const { return (fSumw * fSumw) / fSumw2; }
   Double_t GetEntries() const { return fSumw; }
   Double_t GetEntriesW2() const { return fSumw2; }
   Double_t GetEntriesVW() const { return fSumvw; }
   Double_t GetEntriesWV2() const { return fSumwv2; }
   Double_t GetError() const { return fError; }


private:
   Double_t fSumw;
   Double_t fSumvw;
   Double_t fSumw2;
   Double_t fSumwv2;
   Double_t fAverage;
   Double_t fError;
   EErrorType fErrorMode = kERRORMEAN;

protected:
   void Fill(Double_t value, Double_t weight);
   void UpdateAverage();
   void UpdateError();
   void SetErrorOption(EErrorType type) { fErrorMode = type; }

   ClassDef(TProfile2PolyBin, 1)
};

class TProfile2Poly : public TH2Poly {
    friend class TProfile2PolyBin;

public:
   friend class TProfileHelper;

   TProfile2Poly() {}
   TProfile2Poly(const char *name, const char *title, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup);
   TProfile2Poly(const char *name, const char *title, Int_t nX, Double_t xlow, Double_t xup, Int_t nY, Double_t ylow,
                 Double_t yup);
   virtual ~TProfile2Poly() {}

   using TH2Poly::Fill;
   virtual Int_t Fill(Double_t xcoord, Double_t ycoord, Double_t value) override;
   virtual Int_t Fill(Double_t xcoord, Double_t ycoord, Double_t value, Double_t weight);

   Long64_t Merge(const std::vector<TProfile2Poly *> &list);
   Long64_t Merge(TCollection *in) override;
   virtual void Reset(Option_t *option = "") override;

   // option to dispay different measures on bins
   void SetContentToAverage(); // this one is used by default
   void SetContentToError();

   void SetErrorOption(EErrorType type);

   Double_t GetBinEffectiveEntries(Int_t bin) const;
   Double_t GetBinEntries(Int_t bin) const;
   Double_t GetBinEntriesW2(Int_t bin) const;
   Double_t GetBinEntriesVW(Int_t bin) const;
   Double_t GetBinEntriesWV2(Int_t bin) const;

   using TH2Poly::GetBinContent;
   virtual Double_t GetBinContent(Int_t bin) const override;

   using TH2Poly::GetBinError;
   virtual Double_t GetBinError(Int_t bin) const override;

   virtual void GetStats(Double_t *stats) const override;


   Double_t GetOverflowContent(Int_t idx) { return fOverflowBins[idx].fSumw; }
   void PrintOverflowRegions();

private:
   TProfile2PolyBin fOverflowBins[kNOverflow];
   EErrorType fErrorMode = kERRORMEAN;
   Double_t fTsumwz;
   Double_t fTsumwz2;

protected:
   virtual TProfile2PolyBin *CreateBin(TObject *poly) override;

   Int_t GetOverflowRegionFromCoordinates(Double_t x, Double_t y);
   Int_t OverflowIdxToArrayIdx(Int_t val) { return -val - 1; }


   ClassDefOverride(TProfile2Poly, 2)
};
#endif
