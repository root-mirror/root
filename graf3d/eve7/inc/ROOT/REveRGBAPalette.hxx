// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_REveRGBAPalette
#define ROOT_REveRGBAPalette

#include "ROOT/REveUtil.hxx"

#include "TObject.h"

#include "TMath.h"

namespace ROOT {
namespace Experimental {

class REveRGBAPalette : public TObject,
                        public REveRefCnt
{
   friend class REveRGBAPaletteEditor;
   friend class REveRGBAPaletteSubEditor;

   friend class REveRGBAPaletteOverlay;

public:
   enum ELimitAction_e { kLA_Cut, kLA_Mark, kLA_Clip, kLA_Wrap };

private:
   REveRGBAPalette(const REveRGBAPalette&);            // Not implemented
   REveRGBAPalette& operator=(const REveRGBAPalette&); // Not implemented

protected:
   Double_t  fUIf;       // UI representation calculated as: d = fUIf*i + fUIc
   Double_t  fUIc;       // UI representation calculated as: d = fUIf*i + fUIc

   Int_t     fLowLimit;  // Low  limit for Min/Max values (used by editor)
   Int_t     fHighLimit; // High limit for Min/Max values (used by editor)
   Int_t     fMinVal;
   Int_t     fMaxVal;

   Bool_t    fUIDoubleRep;    // Represent UI parts with real values.
   Bool_t    fInterpolate;    // Interpolate colors for signal values.
   Bool_t    fShowDefValue;   // Flags whether signals with default value should be shown.
   Bool_t    fFixColorRange;  // If true, map palette to low/high limit otherwise to min/max value.
   Int_t     fUnderflowAction;
   Int_t     fOverflowAction;

   Color_t   fDefaultColor;   // Color for when value is not specified
   UChar_t   fDefaultRGBA[4];
   Color_t   fUnderColor;     // Underflow color
   UChar_t   fUnderRGBA[4];
   Color_t   fOverColor;      // Overflow color
   UChar_t   fOverRGBA[4];

   mutable Int_t    fNBins;      // Number of signal-color entries.
   mutable Int_t    fCAMin;      // Minimal signal in color-array.
   mutable Int_t    fCAMax;      // Maximal signal in color-array.
   mutable UChar_t* fColorArray; //[4*fNBins]

   void SetupColor(Int_t val, UChar_t* pix) const;

   Double_t IntToDouble(Int_t i)    const { return fUIf*i + fUIc; }
   Int_t    DoubleToInt(Double_t d) const { return TMath::Nint((d - fUIc) / fUIf); }

   Double_t GetCAMinAsDouble() const { return IntToDouble(fCAMin); }
   Double_t GetCAMaxAsDouble() const { return IntToDouble(fCAMax); }

   static REveRGBAPalette* fgDefaultPalette;

public:
   REveRGBAPalette();
   REveRGBAPalette(Int_t min, Int_t max, Bool_t interp=kTRUE,
                   Bool_t showdef=kTRUE, Bool_t fixcolrng=kFALSE);
   virtual ~REveRGBAPalette();

   void SetupColorArray() const;
   void ClearColorArray();

   Bool_t   WithinVisibleRange(Int_t val) const;
   const UChar_t* ColorFromValue(Int_t val) const;
   void     ColorFromValue(Int_t val, UChar_t* pix, Bool_t alpha=kTRUE) const;
   Bool_t   ColorFromValue(Int_t val, Int_t defVal, UChar_t* pix, Bool_t alpha=kTRUE) const;

   Int_t  GetMinVal() const { return fMinVal; }
   Int_t  GetMaxVal() const { return fMaxVal; }

   void   SetLimits(Int_t low, Int_t high);
   void   SetLimitsScaleMinMax(Int_t low, Int_t high);
   void   SetMinMax(Int_t min, Int_t max);
   void   SetMin(Int_t min);
   void   SetMax(Int_t max);

   Int_t  GetLowLimit()  const { return fLowLimit;  }
   Int_t  GetHighLimit() const { return fHighLimit; }

   // ================================================================

   Bool_t GetUIDoubleRep() const { return fUIDoubleRep; }
   void   SetUIDoubleRep(Bool_t b, Double_t f=1, Double_t c=0);

   Bool_t GetInterpolate() const { return fInterpolate; }
   void   SetInterpolate(Bool_t b);

   Bool_t GetShowDefValue() const { return fShowDefValue; }
   void   SetShowDefValue(Bool_t v) { fShowDefValue = v; }

   Bool_t GetFixColorRange() const { return fFixColorRange; }
   void   SetFixColorRange(Bool_t v);

   Int_t GetUnderflowAction() const  { return fUnderflowAction; }
   Int_t GetOverflowAction()  const  { return fOverflowAction;  }
   void  SetUnderflowAction(Int_t a) { fUnderflowAction = a;    }
   void  SetOverflowAction(Int_t a)  { fOverflowAction  = a;    }

   // ================================================================

   Color_t  GetDefaultColor() const { return fDefaultColor; }
   Color_t* PtrDefaultColor() { return &fDefaultColor; }
   UChar_t* GetDefaultRGBA()  { return fDefaultRGBA;  }
   const UChar_t* GetDefaultRGBA() const { return fDefaultRGBA;  }

   void   SetDefaultColor(Color_t ci);
   void   SetDefaultColorPixel(Pixel_t pix);
   void   SetDefaultColorRGBA(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

   // ----------------------------------------------------------------

   Color_t  GetUnderColor() const { return fUnderColor; }
   Color_t* PtrUnderColor() { return &fUnderColor; }
   UChar_t* GetUnderRGBA()  { return fUnderRGBA;  }
   const UChar_t* GetUnderRGBA() const { return fUnderRGBA;  }

   void   SetUnderColor(Color_t ci);
   void   SetUnderColorPixel(Pixel_t pix);
   void   SetUnderColorRGBA(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

   // ----------------------------------------------------------------

   Color_t  GetOverColor() const { return fOverColor; }
   Color_t* PtrOverColor() { return &fOverColor; }
   UChar_t* GetOverRGBA()  { return fOverRGBA;  }
   const UChar_t* GetOverRGBA() const { return fOverRGBA;  }

   void   SetOverColor(Color_t ci);
   void   SetOverColorPixel(Pixel_t pix);
   void   SetOverColorRGBA(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

   // ================================================================

   ClassDef(REveRGBAPalette, 0); // A generic, speed-optimised mapping from value to RGBA color supporting different wrapping and range truncation modes.
};


/******************************************************************************/
// Inlines for REveRGBAPalette
/******************************************************************************/

//______________________________________________________________________________
inline Bool_t REveRGBAPalette::WithinVisibleRange(Int_t val) const
{
   if ((val < fMinVal && fUnderflowAction == kLA_Cut) ||
       (val > fMaxVal && fOverflowAction  == kLA_Cut))
      return kFALSE;
   else
      return kTRUE;
}

//______________________________________________________________________________
inline const UChar_t* REveRGBAPalette::ColorFromValue(Int_t val) const
{
   // Here we expect that kLA_Cut has been checked; we further check
   // for kLA_Wrap and kLA_Clip otherwise we proceed as for kLA_Mark.

   if (!fColorArray)  SetupColorArray();

   if (val < fMinVal)
   {
      if (fUnderflowAction == kLA_Wrap)
         val = (val+1-fCAMin)%fNBins + fCAMax;
      else if (fUnderflowAction == kLA_Clip)
         val = fMinVal;
      else
         return fUnderRGBA;
   }
   else if(val > fMaxVal)
   {
      if (fOverflowAction == kLA_Wrap)
         val = (val-1-fCAMax)%fNBins + fCAMin;
      else if (fOverflowAction == kLA_Clip)
         val = fMaxVal;
      else
         return fOverRGBA;
   }

   return fColorArray + 4 * (val - fCAMin);
}

//______________________________________________________________________________
inline void REveRGBAPalette::ColorFromValue(Int_t val, UChar_t* pix, Bool_t alpha) const
{
   const UChar_t* c = ColorFromValue(val);
   pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2];
   if (alpha) pix[3] = c[3];
}

//______________________________________________________________________________
inline Bool_t REveRGBAPalette::ColorFromValue(Int_t val, Int_t defVal, UChar_t* pix, Bool_t alpha) const
{
   if (val == defVal) {
      if (fShowDefValue) {
         pix[0] = fDefaultRGBA[0];
         pix[1] = fDefaultRGBA[1];
         pix[2] = fDefaultRGBA[2];
         if (alpha) pix[3] = fDefaultRGBA[3];
         return kTRUE;
      } else {
         return kFALSE;
      }
   }

   if (WithinVisibleRange(val)) {
      ColorFromValue(val, pix, alpha);
      return kTRUE;
   } else {
      return kFALSE;
   }
}

} // namespace Experimental
} // namespace ROOT

#endif
