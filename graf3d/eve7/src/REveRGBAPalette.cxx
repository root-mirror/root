// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/REveRGBAPalette.hxx"

#include "TColor.h"
#include "TStyle.h"
#include "TMath.h"


using namespace ROOT::Experimental;

/** \class REveRGBAPalette
\ingroup REve
A generic, speed-optimised mapping from value to RGBA color
supporting different wrapping and range truncation modes.

Flag fFixColorRange: specifies how the palette is mapped to signal values:
  - true  - LowLimit -> HighLimit
  - false - MinValue -> MaxValue
*/

ClassImp(REveRGBAPalette);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveRGBAPalette::REveRGBAPalette() :
   TObject(),
   REveRefCnt(),

   fUIf(1), fUIc(0),

   fLowLimit(0), fHighLimit(0), fMinVal(0), fMaxVal(0),

   fUIDoubleRep     (kFALSE),
   fInterpolate     (kTRUE),
   fShowDefValue    (kTRUE),
   fFixColorRange   (kFALSE),
   fUnderflowAction (kLA_Cut),
   fOverflowAction  (kLA_Clip),

   fDefaultColor(-1),
   fUnderColor  (-1),
   fOverColor   (-1),

   fNBins(0), fCAMin(0), fCAMax(0), fColorArray(0)
{
   SetLimits(0, 1024);
   SetMinMax(0,  512);

   SetDefaultColor(0);
   SetUnderColor(1);
   SetOverColor(2);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveRGBAPalette::REveRGBAPalette(Int_t min, Int_t max, Bool_t interp,
                                 Bool_t showdef, Bool_t fixcolrng) :
   TObject(),
   REveRefCnt(),

   fUIf(1), fUIc(0),

   fLowLimit(0), fHighLimit(0), fMinVal(0), fMaxVal(0),

   fUIDoubleRep     (kFALSE),
   fInterpolate     (interp),
   fShowDefValue    (showdef),
   fFixColorRange   (fixcolrng),
   fUnderflowAction (kLA_Cut),
   fOverflowAction  (kLA_Clip),

   fDefaultColor(-1),
   fUnderColor  (-1),
   fOverColor   (-1),

   fNBins(0), fCAMin(0), fCAMax(0), fColorArray(0)
{
   SetLimits(min, max);
   SetMinMax(min, max);

   SetDefaultColor(0);
   SetUnderColor(1);
   SetOverColor(2);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveRGBAPalette::~REveRGBAPalette()
{
   delete [] fColorArray;
}

////////////////////////////////////////////////////////////////////////////////
/// Set RGBA color 'pixel' for signal-value 'val'.

void REveRGBAPalette::SetupColor(Int_t val, UChar_t* pixel) const
{
   using namespace TMath;
   Float_t div  = Max(1, fCAMax - fCAMin);
   Int_t   nCol = gStyle->GetNumberOfColors();

   Float_t f;
   if      (val >= fCAMax) f = nCol - 1;
   else if (val <= fCAMin) f = 0;
   else                    f = (val - fCAMin)/div*(nCol - 1);

   if (fInterpolate) {
      Int_t  bin = (Int_t) f;
      Float_t f2 = f - bin, f1 = 1.0f - f2;
      REveUtil::ColorFromIdx(f1, gStyle->GetColorPalette(bin),
                             f2, gStyle->GetColorPalette(Min(bin + 1, nCol - 1)),
                             pixel);
   } else {
      REveUtil::ColorFromIdx(gStyle->GetColorPalette((Int_t) Nint(f)), pixel);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Construct internal color array that maps signal value to RGBA color.

void REveRGBAPalette::SetupColorArray() const
{
   if (fColorArray)
      delete [] fColorArray;

   if (fFixColorRange) {
      fCAMin = fLowLimit; fCAMax = fHighLimit;
   } else {
      fCAMin = fMinVal;   fCAMax = fMaxVal;
   }
   fNBins = fCAMax - fCAMin + 1;

   fColorArray = new UChar_t [4 * fNBins];
   UChar_t* p = fColorArray;
   for(Int_t v = fCAMin; v <= fCAMax; ++v, p+=4)
      SetupColor(v, p);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear internal color array.

void REveRGBAPalette::ClearColorArray()
{
   if (fColorArray) {
      delete [] fColorArray;
      fColorArray = 0;
      fNBins = fCAMin = fCAMax = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set low/high limits on signal value. Current min/max values are
/// clamped into the new limits.

void REveRGBAPalette::SetLimits(Int_t low, Int_t high)
{
   fLowLimit  = low;
   fHighLimit = high;

   if (fMaxVal < fLowLimit)  SetMax(fLowLimit);
   if (fMinVal < fLowLimit)  SetMin(fLowLimit);
   if (fMinVal > fHighLimit) SetMin(fHighLimit);
   if (fMaxVal > fHighLimit) SetMax(fHighLimit);

   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set low/high limits and rescale current min/max values.

void REveRGBAPalette::SetLimitsScaleMinMax(Int_t low, Int_t high)
{
   Float_t rng_old = fHighLimit - fLowLimit;
   Float_t rng_new = high - low;

   fMinVal = TMath::Nint(low + (fMinVal - fLowLimit)*rng_new/rng_old);
   fMaxVal = TMath::Nint(low + (fMaxVal - fLowLimit)*rng_new/rng_old);
   fLowLimit  = low;
   fHighLimit = high;

   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set current min value.

void REveRGBAPalette::SetMin(Int_t min)
{
   fMinVal = TMath::Min(min, fMaxVal);
   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set current max value.

void REveRGBAPalette::SetMax(Int_t max)
{
   fMaxVal = TMath::Max(max, fMinVal);
   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set current min/max values.

void REveRGBAPalette::SetMinMax(Int_t min, Int_t max)
{
   fMinVal = min;
   fMaxVal = max;
   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set flag determining whether GUI editor and overlays should show limits
/// and axis values as real values with mapping from integer value i to real
/// value d as: d = f*i + fc

void REveRGBAPalette::SetUIDoubleRep(Bool_t b, Double_t f, Double_t c)
{
   fUIDoubleRep = b;
   if (fUIDoubleRep) {
      fUIf = f;  fUIc = c;
   } else {
      fUIf = 1;  fUIc = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set interpolation flag. This determines how colors from ROOT's
/// palette are mapped into RGBA values for given signal.

void REveRGBAPalette::SetInterpolate(Bool_t b)
{
   fInterpolate = b;
   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set flag specifying how the palette is mapped to signal values:
///  true  - LowLimit -> HighLimit
///  false - MinValue -> MaxValue

void REveRGBAPalette::SetFixColorRange(Bool_t v)
{
   fFixColorRange = v;
   ClearColorArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Set default color.

void REveRGBAPalette::SetDefaultColor(Color_t ci)
{
   fDefaultColor = ci;
   REveUtil::ColorFromIdx(ci, fDefaultRGBA, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set default color.

void REveRGBAPalette::SetDefaultColorPixel(Pixel_t pix)
{
   SetDefaultColor(Color_t(TColor::GetColor(pix)));
}

////////////////////////////////////////////////////////////////////////////////
/// Set default color.

void REveRGBAPalette::SetDefaultColorRGBA(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
   fDefaultColor = Color_t(TColor::GetColor(r, g, b));
   fDefaultRGBA[0] = r;
   fDefaultRGBA[1] = g;
   fDefaultRGBA[2] = b;
   fDefaultRGBA[3] = a;
}

////////////////////////////////////////////////////////////////////////////////
/// Set underflow color.

void REveRGBAPalette::SetUnderColor(Color_t ci)
{
   fUnderColor = ci;
   REveUtil::ColorFromIdx(ci, fUnderRGBA, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set underflow color.

void REveRGBAPalette::SetUnderColorPixel(Pixel_t pix)
{
   SetUnderColor(Color_t(TColor::GetColor(pix)));
}

////////////////////////////////////////////////////////////////////////////////
/// Set underflow color.

void REveRGBAPalette::SetUnderColorRGBA(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
   fUnderColor = Color_t(TColor::GetColor(r, g, b));
   fUnderRGBA[0] = r;
   fUnderRGBA[1] = g;
   fUnderRGBA[2] = b;
   fUnderRGBA[3] = a;
}

////////////////////////////////////////////////////////////////////////////////
/// Set overflow color.

void REveRGBAPalette::SetOverColor(Color_t ci)
{
   fOverColor = ci;
   REveUtil::ColorFromIdx(ci, fOverRGBA, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set overflow color.

void REveRGBAPalette::SetOverColorPixel(Pixel_t pix)
{
   SetOverColor(Color_t(TColor::GetColor(pix)));
}

////////////////////////////////////////////////////////////////////////////////
/// Set overflow color.

void REveRGBAPalette::SetOverColorRGBA(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
   fOverColor = Color_t(TColor::GetColor(r, g, b));
   fOverRGBA[0] = r;
   fOverRGBA[1] = g;
   fOverRGBA[2] = b;
   fOverRGBA[3] = a;
}
