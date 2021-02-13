// @(#)root/base:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TColor
#define ROOT_TColor

#include "TNamed.h"

class TArrayI;

class TColor : public TNamed {
protected:
   Int_t          fNumber;        ///< Color number identifier
private:
   Float_t        fRed;           ///< Fraction of Red
   Float_t        fGreen;         ///< Fraction of Green
   Float_t        fBlue;          ///< Fraction of Blue
   Float_t        fHue;           ///< Hue
   Float_t        fLight;         ///< Light
   Float_t        fSaturation;    ///< Saturation
   Float_t        fAlpha;         ///< Alpha (transparency)

   void           Allocate();
   static Float_t HLStoRGB1(Float_t rn1, Float_t rn2, Float_t huei);

public:
   TColor();
   TColor(Int_t color, Float_t r, Float_t g, Float_t b, const char *name="", Float_t a = 1);
   TColor(Float_t r, Float_t g, Float_t b, Float_t a = 1);
   TColor(const TColor &color);
   TColor &operator=(const TColor &color);
   virtual ~TColor();
   const char   *AsHexString() const;
   void          Copy(TObject &color) const;
   static void   CreateColorWheel();
   static void   CreateColorsGray();
   static void   CreateColorsCircle(Int_t offset, const char *name, UChar_t *rgb);
   static void   CreateColorsRectangle(Int_t offset, const char *name, UChar_t *rgb);
   static Int_t  CreateGradientColorTable(UInt_t Number, Double_t* Stops,
                    Double_t* Red, Double_t* Green, Double_t* Blue, UInt_t NColors, Float_t alpha=1.);
   static Int_t  GetColorPalette(Int_t i);
   static Int_t  GetNumberOfColors();
   virtual void  GetRGB(Float_t &r, Float_t &g, Float_t &b) const
                    { r=GetRed(); g=GetGreen(); b=GetBlue(); }
   virtual void  GetHLS(Float_t &h, Float_t &l, Float_t &s) const
                    { h=GetHue(); l=GetLight(); s=GetSaturation(); }
   Int_t         GetNumber() const { return fNumber; }
   ULong_t       GetPixel() const;
   Float_t       GetRed() const { return IsGrayscale() ? GetGrayscale() : fRed; }
   Float_t       GetGreen() const { return IsGrayscale() ? GetGrayscale() : fGreen; }
   Float_t       GetBlue() const { return IsGrayscale() ? GetGrayscale() : fBlue; }
   Float_t       GetHue() const { return fHue; }
   Float_t       GetLight() const { return fLight; }
   Float_t       GetSaturation() const { return IsGrayscale() ? 0 : fSaturation; }
   Float_t       GetAlpha() const { return fAlpha; }
   virtual Float_t GetGrayscale() const { /*ITU*/ return 0.299f*fRed + 0.587f*fGreen + 0.114f*fBlue; }
   virtual void  ls(Option_t *option="") const;
   virtual void  Print(Option_t *option="") const;
   virtual void  SetAlpha(Float_t a) { fAlpha = a; }
   virtual void  SetRGB(Float_t r, Float_t g, Float_t b);

   static void    InitializeColors();
   static void    HLS2RGB(Float_t h, Float_t l, Float_t s, Float_t &r, Float_t &g, Float_t &b);
   static void    HLS2RGB(Int_t h, Int_t l, Int_t s, Int_t &r, Int_t &g, Int_t &b);
   static void    HLStoRGB(Float_t h, Float_t l, Float_t s, Float_t &r, Float_t &g, Float_t &b)
                     { TColor::HLS2RGB(h, l, s, r, g, b); } // backward compatible
   static void    HSV2RGB(Float_t h, Float_t s, Float_t v, Float_t &r, Float_t &g, Float_t &b);
   static void    RGB2HLS(Float_t r, Float_t g, Float_t b, Float_t &h, Float_t &l, Float_t &s);
   static void    RGB2HLS(Int_t r, Int_t g, Int_t b, Int_t &h, Int_t &l, Int_t &s);
   static void    RGBtoHLS(Float_t r, Float_t g, Float_t b, Float_t &h, Float_t &l, Float_t &s)
                     { TColor::RGB2HLS(r, g, b, h, l, s); } // backward compatible
   static void    RGB2HSV(Float_t r, Float_t g, Float_t b, Float_t &h, Float_t &s, Float_t &v);
   static Int_t   GetColor(const char *hexcolor);
   static Int_t   GetColor(Float_t r, Float_t g, Float_t b);
   static Int_t   GetColor(Int_t r, Int_t g, Int_t b);
   static Int_t   GetColor(ULong_t pixel);
   static Int_t   GetColorBright(Int_t color);
   static Int_t   GetColorDark(Int_t color);
   static Int_t   GetColorTransparent(Int_t color, Float_t a);
   static Int_t   GetFreeColorIndex();
   static const TArrayI& GetPalette();
   static ULong_t Number2Pixel(Int_t ci);
   static ULong_t RGB2Pixel(Int_t r, Int_t g, Int_t b);
   static ULong_t RGB2Pixel(Float_t r, Float_t g, Float_t b);
   static void    Pixel2RGB(ULong_t pixel, Int_t &r, Int_t &g, Int_t &b);
   static void    Pixel2RGB(ULong_t pixel, Float_t &r, Float_t &g, Float_t &b);
   static const char *PixelAsHexString(ULong_t pixel);
   static void    SaveColor(std::ostream &out, Int_t ci);
   static void    SetColorThreshold(Float_t t);
   static Bool_t  DefinedColors();
   static void    InvertPalette();
   static Bool_t  IsGrayscale();
   static void    SetGrayscale(Bool_t set = kTRUE);
   static void    SetPalette(Int_t ncolors, Int_t *colors,Float_t alpha=1.);

   ClassDef(TColor,2)  //Color defined by RGB or HLS
};

   enum EColorPalette {kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
                       kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
                       kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
                       kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
                       kAlpine=63,           kAquamarine=64,   kArmy=65,
                       kAtlantic=66,         kAurora=67,       kAvocado=68,
                       kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
                       kBrownCyan=72,        kCMYK=73,         kCandy=74,
                       kCherry=75,           kCoffee=76,       kDarkRainBow=77,
                       kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
                       kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
                       kGreenPink=84,        kIsland=85,       kLake=86,
                       kLightTemperature=87, kLightTerrain=88, kMint=89,
                       kNeon=90,             kPastel=91,       kPearl=92,
                       kPigeon=93,           kPlum=94,         kRedBlue=95,
                       kRose=96,             kRust=97,         kSandyTerrain=98,
                       kSienna=99,           kSolar=100,       kSouthWest=101,
                       kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
                       kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
                       kWaterMelon=108,      kCool=109,        kCopper=110,
                       kGistEarth=111,       kViridis=112,     kCividis=113};
#endif

