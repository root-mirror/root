// @(#)root/hist:$Id$
// Author: Filip Ilic

#ifndef ROOT_TProfile2Poly
#define ROOT_TProfile2Poly

#include "TH2Poly.h"
#include "TProfile.h"

class TProfile2PolyBin : public TH2PolyBin {
public:
   friend class TProfile2Poly;

   TProfile2PolyBin();
   TProfile2PolyBin(TObject *poly, Int_t bin_number);
   virtual ~TProfile2PolyBin() {}

   void Merge(const TProfile2PolyBin *toMerge);

   void Update();
   void ClearStats();

   Double_t GetEffectiveEntries() { return (fSumw * fSumw) / fSumw2; }
   Double_t GetEntries() { return fSumw; }
   Double_t GetBinError() { return fError; }


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

   Long64_t Merge(std::vector<TProfile2Poly *> list);
   Long64_t Merge(TCollection *in) override;
   virtual void Reset(Option_t *option = "") override;

   Int_t GetOverflowRegionFromCoordinates(Double_t x, Double_t y);
   Int_t OverflowIdxToArrayIdx(Int_t val) { return -val - 1; }

   // option to dispay different measures on bins
   void SetContentToAverage(); // this one is used by default
   void SetContentToError();

   void SetErrorOption(EErrorType type);

   Double_t GetBinEffectiveEntries(Int_t bin) const;
   Double_t GetBinEntries(Int_t bin) const;

   using TH2Poly::GetBinError;
   virtual Double_t GetBinError(Int_t bin) const override;

   Double_t GetOverflowContent(Int_t idx) { return regions[idx].fSumw; }
   void printOverflowRegions();

private:
   TProfile2PolyBin regions[kNOverflow];
   EErrorType fErrorMode = kERRORMEAN;
   Double_t fTsumwz;
   Double_t fTsumwz2;

protected:
   virtual TProfile2PolyBin *CreateBin(TObject *poly) override;

   ClassDefOverride(TProfile2Poly, 1)
};
#endif
