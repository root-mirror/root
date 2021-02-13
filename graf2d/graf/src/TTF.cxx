// @(#)root/graf:$Id$
// Author: Olivier Couet     01/10/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTF
\ingroup BasicGraphics

Interface to the freetype 2 library.
*/

#  include <ft2build.h>
#  include FT_FREETYPE_H
#  include FT_GLYPH_H
#include "TROOT.h"
#include "TTF.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TMath.h"
#include "TError.h"

// to scale fonts to the same size as the old TT version
const Float_t kScale = 0.93376068;

TTF gCleanupTTF; // Allows to call "Cleanup" at the end of the session

Bool_t         TTF::fgInit           = kFALSE;
Bool_t         TTF::fgSmoothing      = kTRUE;
Bool_t         TTF::fgKerning        = kTRUE;
Bool_t         TTF::fgHinting        = kFALSE;
Int_t          TTF::fgTBlankW        = 0;
Int_t          TTF::fgWidth          = 0;
Int_t          TTF::fgAscent         = 0;
Int_t          TTF::fgCurFontIdx     = -1;
Int_t          TTF::fgSymbItaFontIdx = -1;
Int_t          TTF::fgFontCount      = 0;
Int_t          TTF::fgNumGlyphs      = 0;
char          *TTF::fgFontName[kTTMaxFonts];
FT_Matrix     *TTF::fgRotMatrix      = nullptr;
FT_Library     TTF::fgLibrary;
FT_BBox        TTF::fgCBox;
FT_Face        TTF::fgFace[kTTMaxFonts];
FT_CharMap     TTF::fgCharMap[kTTMaxFonts];
TTF::TTGlyph   TTF::fgGlyphs[kMaxGlyphs];

ClassImp(TTF);

////////////////////////////////////////////////////////////////////////////////
/// Cleanup TTF environment.

TTF::~TTF()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialise the TrueType fonts interface.

void TTF::Init()
{
   fgInit = kTRUE;

   // initialize FTF library
   if (FT_Init_FreeType(&fgLibrary)) {
      Error("TTF::Init", "error initializing FreeType");
      return;
   }

   // load default font (arialbd)
   SetTextFont(62);
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup. Is called by the gCleanupTTF destructor.

void TTF::Cleanup()
{
   if (!fgInit) return;

   for (int i = 0; i < fgFontCount; i++) {
      delete [] fgFontName[i];
      FT_Done_Face(fgFace[i]);
   }
   if (fgRotMatrix) delete fgRotMatrix;
   FT_Done_FreeType(fgLibrary);

   fgInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Map char to unicode. Returns 0 in case no mapping exists.

Short_t TTF::CharToUnicode(UInt_t code)
{
   if (!fgCharMap[fgCurFontIdx]) {
      UShort_t i, platform, encoding;
      FT_CharMap  charmap;

      if (!fgFace[fgCurFontIdx]) return 0;
      Int_t n = fgFace[fgCurFontIdx]->num_charmaps;
      for (i = 0; i < n; i++) {
         if (!fgFace[fgCurFontIdx]) continue;
         charmap  = fgFace[fgCurFontIdx]->charmaps[i];
         platform = charmap->platform_id;
         encoding = charmap->encoding_id;
         if ((platform == 3 && encoding == 1) ||
             (platform == 0 && encoding == 0) ||
             (platform == 1 && encoding == 0 &&
              !strcmp(fgFontName[fgCurFontIdx], "wingding.ttf")) ||
             (platform == 1 && encoding == 0 &&
              !strcmp(fgFontName[fgCurFontIdx], "symbol.ttf")))
         {
            fgCharMap[fgCurFontIdx] = charmap;
            if (FT_Set_Charmap(fgFace[fgCurFontIdx], fgCharMap[fgCurFontIdx]))
                Error("TTF::CharToUnicode", "error in FT_Set_CharMap");
            return FT_Get_Char_Index(fgFace[fgCurFontIdx], (FT_ULong)code);
         }
      }
   }
   return FT_Get_Char_Index(fgFace[fgCurFontIdx], (FT_ULong)code);
}

////////////////////////////////////////////////////////////////////////////////
/// Get width (w) and height (h) when text is horizontal.

void TTF::GetTextExtent(UInt_t &w, UInt_t &h, char *text)
{
   if (!fgInit) Init();

   SetRotationMatrix(0);
   PrepareString(text);
   LayoutGlyphs();
   Int_t Xoff = 0; if (fgCBox.xMin < 0) Xoff = -fgCBox.xMin;
   Int_t Yoff = 0; if (fgCBox.yMin < 0) Yoff = -fgCBox.yMin;
   w = fgCBox.xMax + Xoff + fgTBlankW;
   h = fgCBox.yMax + Yoff;
}

////////////////////////////////////////////////////////////////////////////////
/// Get advance (a) when text is horizontal.

void TTF::GetTextAdvance(UInt_t &a, char *text)
{
   if (!fgInit) Init();

   SetRotationMatrix(0);
   PrepareString(text);
   LayoutGlyphs();
   a = GetWidth()>>6;
}

////////////////////////////////////////////////////////////////////////////////
/// Get width (w) and height (h) when text is horizontal.

void TTF::GetTextExtent(UInt_t &w, UInt_t &h, wchar_t *text)
{
   if (!fgInit) Init();

   SetRotationMatrix(0);
   PrepareString(text);
   LayoutGlyphs();
   Int_t Xoff = 0; if (fgCBox.xMin < 0) Xoff = -fgCBox.xMin;
   Int_t Yoff = 0; if (fgCBox.yMin < 0) Yoff = -fgCBox.yMin;
   w = fgCBox.xMax + Xoff + fgTBlankW;
   h = fgCBox.yMax + Yoff;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the glyphs positions, fgAscent and fgWidth (needed for alignment).
/// Perform the Glyphs transformation.
/// Compute the string control box.
/// If required take the "kerning" into account.
/// SetRotation and PrepareString should have been called before.

void TTF::LayoutGlyphs()
{
   TTGlyph*  glyph = fgGlyphs;
   FT_Vector origin;
   FT_UInt   load_flags;
   FT_UInt   prev_index = 0;

   fgAscent = 0;
   fgWidth  = 0;

   load_flags = FT_LOAD_DEFAULT;
   if (!fgHinting) load_flags |= FT_LOAD_NO_HINTING;

   fgCBox.xMin = fgCBox.yMin =  32000;
   fgCBox.xMax = fgCBox.yMax = -32000;

   for (int n = 0; n < fgNumGlyphs; n++, glyph++) {

      // compute glyph origin
      if (fgKerning) {
         if (prev_index) {
            FT_Vector  kern;
            FT_Get_Kerning(fgFace[fgCurFontIdx], prev_index, glyph->fIndex,
                           fgHinting ? ft_kerning_default : ft_kerning_unfitted,
                           &kern);
            fgWidth += kern.x;
         }
         prev_index = glyph->fIndex;
      }

      origin.x = fgWidth;
      origin.y = 0;

      // clear existing image if there is one
      if (glyph->fImage) {
         FT_Done_Glyph(glyph->fImage);
         glyph->fImage = nullptr;
      }

      // load the glyph image (in its native format)
      if (FT_Load_Glyph(fgFace[fgCurFontIdx], glyph->fIndex, load_flags))
         continue;

      // extract the glyph image
      if (FT_Get_Glyph (fgFace[fgCurFontIdx]->glyph, &glyph->fImage))
         continue;

      glyph->fPos = origin;
      fgWidth    += fgFace[fgCurFontIdx]->glyph->advance.x;
      fgAscent    = TMath::Max((Int_t)(fgFace[fgCurFontIdx]->glyph->metrics.horiBearingY), fgAscent);

      // transform the glyphs
      FT_Vector_Transform(&glyph->fPos, fgRotMatrix);
      if (FT_Glyph_Transform(glyph->fImage, fgRotMatrix, &glyph->fPos))
         continue;

      // compute the string control box
      FT_BBox  bbox;
      FT_Glyph_Get_CBox(glyph->fImage, ft_glyph_bbox_pixels, &bbox);
      if (bbox.xMin < fgCBox.xMin) fgCBox.xMin = bbox.xMin;
      if (bbox.yMin < fgCBox.yMin) fgCBox.yMin = bbox.yMin;
      if (bbox.xMax > fgCBox.xMax) fgCBox.xMax = bbox.xMax;
      if (bbox.yMax > fgCBox.yMax) fgCBox.yMax = bbox.yMax;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Put the characters in "string" in the "glyphs" array.

void TTF::PrepareString(const char *string)
{
   const unsigned char *p = (const unsigned char*) string;
   TTGlyph *glyph = fgGlyphs;
   UInt_t index;       // Unicode value
   Int_t NbTBlank = 0; // number of trailing blanks

   fgTBlankW   = 0;
   fgNumGlyphs = 0;
   while (*p) {
      index = CharToUnicode((FT_ULong)*p);
      if (index != 0) {
         glyph->fIndex = index;
         glyph++;
         fgNumGlyphs++;
      }
      if (*p == ' ') {
         NbTBlank++;
      } else {
         NbTBlank = 0;
      }
      if (fgNumGlyphs >= kMaxGlyphs) break;
      p++;
   }

   // compute the trailing blanks width. It is use to compute the text
   // width in GetTextExtent
   if (NbTBlank) {
      FT_UInt load_flags = FT_LOAD_DEFAULT;
      if (!fgHinting) load_flags |= FT_LOAD_NO_HINTING;
      if (FT_Load_Glyph(fgFace[fgCurFontIdx], 3, load_flags)) return;
      fgTBlankW = (Int_t)((fgFace[fgCurFontIdx]->glyph->advance.x)>>6)*NbTBlank;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Put the characters in "string" in the "glyphs" array.

void TTF::PrepareString(const wchar_t *string)
{
   const wchar_t *p = string;
   TTGlyph *glyph = fgGlyphs;
   UInt_t index;       // Unicode value
   Int_t NbTBlank = 0; // number of trailing blanks

   fgTBlankW   = 0;
   fgNumGlyphs = 0;
   while (*p) {
      index = FT_Get_Char_Index(fgFace[fgCurFontIdx], (FT_ULong)*p);
      if (index != 0) {
         glyph->fIndex = index;
         glyph++;
         fgNumGlyphs++;
      }
      if (*p == ' ') {
         NbTBlank++;
      } else {
         NbTBlank = 0;
      }
      if (fgNumGlyphs >= kMaxGlyphs) break;
      p++;
   }

   // compute the trailing blanks width. It is use to compute the text
   // width in GetTextExtent
   if (NbTBlank) {
      FT_UInt load_flags = FT_LOAD_DEFAULT;
      if (!fgHinting) load_flags |= FT_LOAD_NO_HINTING;
      if (FT_Load_Glyph(fgFace[fgCurFontIdx], 3, load_flags)) return;
      fgTBlankW = (Int_t)((fgFace[fgCurFontIdx]->glyph->advance.x)>>6)*NbTBlank;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set hinting flag.

void TTF::SetHinting(Bool_t state)
{
   fgHinting = state;
}

////////////////////////////////////////////////////////////////////////////////
/// Set kerning flag.

void TTF::SetKerning(Bool_t state)
{
   fgKerning = state;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the rotation matrix used to rotate the font outlines.

void TTF::SetRotationMatrix(Float_t angle)
{
   Float_t rangle = Float_t(angle * TMath::Pi() / 180.); // Angle in radian
#if defined(FREETYPE_PATCH) && \
    (FREETYPE_MAJOR == 2) && (FREETYPE_MINOR == 1) && (FREETYPE_PATCH == 2)
   Float_t sin    = TMath::Sin(rangle);
   Float_t cos    = TMath::Cos(rangle);
#else
   Float_t sin    = TMath::Sin(-rangle);
   Float_t cos    = TMath::Cos(-rangle);
#endif

   if (!fgRotMatrix) fgRotMatrix = new FT_Matrix;

   fgRotMatrix->xx = (FT_Fixed) (cos * (1<<16));
   fgRotMatrix->xy = (FT_Fixed) (sin * (1<<16));
   fgRotMatrix->yx = -fgRotMatrix->xy;
   fgRotMatrix->yy =  fgRotMatrix->xx;
}

////////////////////////////////////////////////////////////////////////////////
/// Set smoothing (anti-aliasing) flag.

void TTF::SetSmoothing(Bool_t state)
{
   fgSmoothing = state;
}

////////////////////////////////////////////////////////////////////////////////
/// Set text font to specified name.
///  - font       : font name
///  - italic     : the fonts should be slanted. Used for symbol font.
///
/// Set text font to specified name. This function returns 0 if
/// the specified font is found, 1 if not.

Int_t TTF::SetTextFont(const char *fontname, Int_t italic)
{
   if (!fgInit) Init();

   if (!fontname || !fontname[0]) {
      Warning("TTF::SetTextFont",
              "no font name specified, using default font %s", fgFontName[0]);
      fgCurFontIdx = 0;
      return 0;
   }
   const char *basename = gSystem->BaseName(fontname);

   // check if font is in cache
   int i;
   for (i = 0; i < fgFontCount; i++) {
      if (!strcmp(fgFontName[i], basename)) {
         if (italic) {
            if (i==fgSymbItaFontIdx) {
               fgCurFontIdx = i;
               return 0;
            }
         } else {
            if (i!=fgSymbItaFontIdx) {
               fgCurFontIdx = i;
               return 0;
            }
         }
      }
   }

   // enough space in cache to load font?
   if (fgFontCount >= kTTMaxFonts) {
      Error("TTF::SetTextFont", "too many fonts opened (increase kTTMaxFont = %d)",
            kTTMaxFonts);
      Warning("TTF::SetTextFont", "using default font %s", fgFontName[0]);
      fgCurFontIdx = 0;    // use font 0 (default font, set in ctor)
      return 0;
   }

   // try to load font (font must be in Root.TTFontPath resource)
   const char *ttpath = gEnv->GetValue("Root.TTFontPath",
                                       TROOT::GetTTFFontDir());
   char *ttfont = gSystem->Which(ttpath, fontname, kReadPermission);

   if (!ttfont) {
      Error("TTF::SetTextFont", "font file %s not found in path", fontname);
      if (fgFontCount) {
         Warning("TTF::SetTextFont", "using default font %s", fgFontName[0]);
         fgCurFontIdx = 0;    // use font 0 (default font, set in ctor)
         return 0;
      } else {
         return 1;
      }
   }

   FT_Face  tface = 0;

   if (FT_New_Face(fgLibrary, ttfont, 0, &tface)) {
      Error("TTF::SetTextFont", "error loading font %s", ttfont);
      delete [] ttfont;
      if (tface) FT_Done_Face(tface);
      if (fgFontCount) {
         Warning("TTF::SetTextFont", "using default font %s", fgFontName[0]);
         fgCurFontIdx = 0;
         return 0;
      } else {
         return 1;
      }
   }

   delete [] ttfont;

   fgFontName[fgFontCount] = StrDup(basename);
   fgCurFontIdx            = fgFontCount;
   fgFace[fgCurFontIdx]    = tface;
   fgCharMap[fgCurFontIdx] = 0;
   fgFontCount++;

   if (italic) {
      fgSymbItaFontIdx = fgCurFontIdx;
      FT_Matrix slantMat;
      slantMat.xx = (1 << 16);
      slantMat.xy = ((1 << 16) >> 2);
      slantMat.yx = 0;
      slantMat.yy = (1 << 16);
      FT_Set_Transform( fgFace[fgSymbItaFontIdx], &slantMat, NULL );
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set specified font.
/// List of the currently supported fonts (screen and PostScript)
///
/// | Font ID |   X11                     |     TTF          |
/// |---------|---------------------------|------------------|
/// |      1  | times-medium-i-normal     | timesi.ttf       |
/// |      2  | times-bold-r-normal       | timesbd.ttf      |
/// |      3  | times-bold-i-normal       | timesbi.ttf      |
/// |      4  | helvetica-medium-r-normal | arial.ttf        |
/// |      5  | helvetica-medium-o-normal | ariali.ttf       |
/// |      6  | helvetica-bold-r-normal   | arialbd.ttf      |
/// |      7  | helvetica-bold-o-normal   | arialbi.ttf      |
/// |      8  | courier-medium-r-normal   | cour.ttf         |
/// |      9  | courier-medium-o-normal   | couri.ttf        |
/// |     10  | courier-bold-r-normal     | courbd.ttf       |
/// |     11  | courier-bold-o-normal     | courbi.ttf       |
/// |     12  | symbol-medium-r-normal    | symbol.ttf       |
/// |     13  | times-medium-r-normal     | times.ttf        |
/// |     14  |                           | wingding.ttf     |
/// |     15  | symbol oblique is emulated from symbol.ttf | |

void TTF::SetTextFont(Font_t fontnumber)
{
   // Added by cholm for use of DFSG - fonts - based on Kevins fix.
   // Table of Microsoft and (for non-MSFT operating systems) backup
   // FreeFont TTF fonts.
   static const char *fonttable[][2] = {
     { "Root.TTFont.0", "FreeSansBold.otf" },
     { "Root.TTFont.1", "FreeSerifItalic.otf" },
     { "Root.TTFont.2", "FreeSerifBold.otf" },
     { "Root.TTFont.3", "FreeSerifBoldItalic.otf" },
     { "Root.TTFont.4", "FreeSans.otf" },
     { "Root.TTFont.5", "FreeSansOblique.otf" },
     { "Root.TTFont.6", "FreeSansBold.otf" },
     { "Root.TTFont.7", "FreeSansBoldOblique.otf" },
     { "Root.TTFont.8", "FreeMono.otf" },
     { "Root.TTFont.9", "FreeMonoOblique.otf" },
     { "Root.TTFont.10", "FreeMonoBold.otf" },
     { "Root.TTFont.11", "FreeMonoBoldOblique.otf" },
     { "Root.TTFont.12", "symbol.ttf" },
     { "Root.TTFont.13", "FreeSerif.otf" },
     { "Root.TTFont.14", "wingding.ttf" },
     { "Root.TTFont.15", "symbol.ttf" },
     { "Root.TTFont.STIXGen", "STIXGeneral.otf" },
     { "Root.TTFont.STIXGenIt", "STIXGeneralItalic.otf" },
     { "Root.TTFont.STIXGenBd", "STIXGeneralBol.otf" },
     { "Root.TTFont.STIXGenBdIt", "STIXGeneralBolIta.otf" },
     { "Root.TTFont.STIXSiz1Sym", "STIXSiz1Sym.otf" },
     { "Root.TTFont.STIXSiz1SymBd", "STIXSiz1SymBol.otf" },
     { "Root.TTFont.STIXSiz2Sym", "STIXSiz2Sym.otf" },
     { "Root.TTFont.STIXSiz2SymBd", "STIXSiz2SymBol.otf" },
     { "Root.TTFont.STIXSiz3Sym", "STIXSiz3Sym.otf" },
     { "Root.TTFont.STIXSiz3SymBd", "STIXSiz3SymBol.otf" },
     { "Root.TTFont.STIXSiz4Sym", "STIXSiz4Sym.otf" },
     { "Root.TTFont.STIXSiz4SymBd", "STIXSiz4SymBol.otf" },
     { "Root.TTFont.STIXSiz5Sym", "STIXSiz5Sym.otf" },
     { "Root.TTFont.ME", "DroidSansFallback.ttf" },
     { "Root.TTFont.CJKMing", "DroidSansFallback.ttf" },
     { "Root.TTFont.CJKGothic", "DroidSansFallback.ttf" }
   };

   static int fontset = -1;
   int        thisset = fontset;

   int fontid = fontnumber / 10;
   if (fontid < 0 || fontid > 31) fontid = 0;

   if (thisset == -1) {
      // try to load font (font must be in Root.TTFontPath resource)
      // to see which fontset we have available
      const char *ttpath = gEnv->GetValue("Root.TTFontPath",
                                          TROOT::GetTTFFontDir());
      char *ttfont = gSystem->Which(ttpath, gEnv->GetValue(fonttable[fontid][0], fonttable[fontid][1]), kReadPermission);
      if (ttfont) {
         delete [] ttfont;
         thisset = 0;
      } else {
         // try backup free font
         thisset = 1;
      }
   }
   Int_t italic = 0;
   if (fontid==15) italic = 1;
   int ret = SetTextFont(gEnv->GetValue(fonttable[fontid][thisset], fonttable[fontid][1]), italic);
   // Do not define font set is we're loading the symbol.ttf - it's
   // the same in both cases.
   if (ret == 0 && fontid != 12) fontset = thisset;
}

////////////////////////////////////////////////////////////////////////////////
/// Set current text size.

void TTF::SetTextSize(Float_t textsize)
{
   if (!fgInit) Init();
   if (textsize < 0) return;

   if (fgCurFontIdx < 0 || fgFontCount <= fgCurFontIdx) {
      Error("TTF::SetTextSize", "current font index out of bounds");
      fgCurFontIdx = 0;
      return;
   }

   Int_t tsize = (Int_t)(textsize*kScale+0.5) << 6;
   if (FT_Set_Char_Size(fgFace[fgCurFontIdx], tsize, tsize, 72, 72))
      Error("TTF::SetTextSize", "error in FT_Set_Char_Size");
}

////////////////////////////////////////////////////////////////////////////////

void TTF::Version(Int_t &major, Int_t &minor, Int_t &patch)
{
   FT_Library_Version(fgLibrary, &major, &minor, &patch);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTF::GetHinting()
{
    return fgHinting;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTF::GetKerning()
{
    return fgKerning;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTF::GetSmoothing()
{
    return fgSmoothing;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTF::IsInitialized()
{
    return fgInit;
}

////////////////////////////////////////////////////////////////////////////////

Int_t  TTF::GetWidth()
{
    return fgWidth;
}

////////////////////////////////////////////////////////////////////////////////

Int_t  TTF::GetAscent()
{
    return fgAscent;
}

////////////////////////////////////////////////////////////////////////////////

Int_t  TTF::GetNumGlyphs()
{
    return fgNumGlyphs;
}

////////////////////////////////////////////////////////////////////////////////

FT_Matrix *TTF::GetRotMatrix()
{
    return fgRotMatrix;
}

////////////////////////////////////////////////////////////////////////////////

const FT_BBox &TTF::GetBox()
{
    return fgCBox;
}

////////////////////////////////////////////////////////////////////////////////

TTF::TTGlyph *TTF::GetGlyphs()
{
    return fgGlyphs;
}
