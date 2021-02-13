// @(#)root/gl:$Id$
// Author: Alja Mrak-Tadel 2008

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RConfigure.h"
#include "TGLFontManager.h"


#include "TROOT.h"
#include "TVirtualX.h"
#include "TMath.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TObjString.h"
#include "TGLUtil.h"
#include "TGLIncludes.h"

// Direct inclusion of FTGL headers is deprecated in ftgl-2.1.3 while
// ftgl-2.1.2 shipped with ROOT requires manual inclusion.
#ifndef BUILTIN_FTGL
# include <FTGL/ftgl.h>
#else
# include "FTFont.h"
# include "FTGLExtrdFont.h"
# include "FTGLOutlineFont.h"
# include "FTGLPolygonFont.h"
# include "FTGLTextureFont.h"
# include "FTGLPixmapFont.h"
# include "FTGLBitmapFont.h"
#endif


/** \class TGLFont
\ingroup opengl
A wrapper class for FTFont.
Holds pointer to FTFont object and its description: face size, font file
and class ID. It  wraps Render and BBox functions.
*/

ClassImp(TGLFont);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLFont::TGLFont():
   fFont(0), fManager(0), fDepth(0),
   fSize(0), fFile(0), fMode(kUndef),
   fTrashCount(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLFont::TGLFont(Int_t size, Int_t font, EMode mode, FTFont* f, TGLFontManager* mng):
   fFont(f), fManager(mng), fDepth(0),
   fSize(size), fFile(font), fMode(mode),
   fTrashCount(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TGLFont::TGLFont(const TGLFont &o):
   fFont(0), fManager(0), fDepth(0), fTrashCount(0)
{
   fFont = (FTFont*)o.GetFont();

   fSize  = o.fSize;
   fFile  = o.fFile;
   fMode  = o.fMode;

   fTrashCount = o.fTrashCount;
}

////////////////////////////////////////////////////////////////////////////////
///Destructor

TGLFont::~TGLFont()
{
   if (fManager) fManager->ReleaseFont(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

void TGLFont::CopyAttributes(const TGLFont &o)
{
   SetFont(o.fFont);
   SetManager(o.fManager);

   SetDepth(o.fDepth);

   fSize  = o.fSize;
   fFile  = o.fFile;
   fMode  = o.fMode;

   fTrashCount = o.fTrashCount;
}


/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Get font's ascent.

Float_t TGLFont::GetAscent() const
{
   return fFont->Ascender();
}

////////////////////////////////////////////////////////////////////////////////
/// Get font's descent. The returned value is positive.

Float_t TGLFont::GetDescent() const
{
   return -fFont->Descender();
}

////////////////////////////////////////////////////////////////////////////////
/// Get font's line-height.

Float_t TGLFont::GetLineHeight() const
{
   return fFont->LineHeight();
}

////////////////////////////////////////////////////////////////////////////////
/// Measure font's base-line parameters from the passed text.
/// Note that the measured parameters are not the same as the ones
/// returned by get-functions - those were set by the font designer.

void TGLFont::MeasureBaseLineParams(Float_t& ascent, Float_t& descent, Float_t& line_height,
                                    const char* txt) const
{
   Float_t dum, lly, ury;
   const_cast<FTFont*>(fFont)->BBox(txt, dum, lly, dum, dum, ury, dum);
   ascent      =  ury;
   descent     = -lly;
   line_height =  ury - lly;
}

////////////////////////////////////////////////////////////////////////////////
/// Get bounding box.

void TGLFont::BBox(const char* txt,
                   Float_t& llx, Float_t& lly, Float_t& llz,
                   Float_t& urx, Float_t& ury, Float_t& urz) const
{
   // FTGL is not const correct.
   const_cast<FTFont*>(fFont)->BBox(txt, llx, lly, llz, urx, ury, urz);
}

////////////////////////////////////////////////////////////////////////////////
/// Get bounding box.

void TGLFont::BBox(const wchar_t* txt,
                   Float_t& llx, Float_t& lly, Float_t& llz,
                   Float_t& urx, Float_t& ury, Float_t& urz) const
{
   // FTGL is not const correct.
   const_cast<FTFont*>(fFont)->BBox(txt, llx, lly, llz, urx, ury, urz);
}

////////////////////////////////////////////////////////////////////////////////
///mgn is simply ignored, because ROOT's TVirtualX TGX11 are complete mess with
///painting attributes.

template<class Char>
void TGLFont::RenderHelper(const Char *txt, Double_t x, Double_t y, Double_t angle, Double_t /*mgn*/) const
{
   glPushMatrix();
   //glLoadIdentity();

   // FTGL is not const correct.
   Float_t llx = 0.f, lly = 0.f, llz = 0.f, urx = 0.f, ury = 0.f, urz = 0.f;
   BBox(txt, llx, lly, llz, urx, ury, urz);

   /*
    V\H   | left | center | right
   _______________________________
   bottom |  7   |   8    |   9
   _______________________________
   center |  4   |   5    |   6
   _______________________________
    top   |  1   |   2    |   3
   */
   const Double_t dx = urx - llx, dy = ury - lly;
   Double_t xc = 0., yc = 0.;
   const UInt_t align = gVirtualX->GetTextAlign();

   //Here's the nice X11 bullshido: you call gVirtualX->SetTextAlign(11),
   //later gVirtualX->GetTextAling() will give you 7. Brilliant!
   //But with Cocoa you'll have 11. As it should be, of course.

   if (gVirtualX->InheritsFrom("TGCocoa")) {
      const UInt_t hAlign = UInt_t(align / 10);
      switch (hAlign) {
      case 1:
         xc = 0.5 * dx;
         break;
      case 2:
         break;
      case 3:
         xc = -0.5 * dy;
         break;
      }

      const UInt_t vAlign = UInt_t(align % 10);
      switch (vAlign) {
      case 1:
         yc = 0.5 * dy;
         break;
      case 2:
         break;
      case 3:
         yc = -0.5 * dy;
         break;
      }
   } else {
      switch (align) {
      case 7:
         xc += 0.5 * dx;
         yc += 0.5 * dy;
         break;
      case 8:
         yc += 0.5 * dy;
         break;
      case 9:
         xc -= 0.5 * dx;
         yc += 0.5 * dy;
         break;
      case 4:
         xc += 0.5 * dx;
         break;
      case 5:
         break;
      case 6:
         xc = -0.5 * dx;
         break;
      case 1:
         xc += 0.5 * dx;
         yc -= 0.5 * dy;
         break;
      case 2:
         yc -= 0.5 * dy;
         break;
      case 3:
         xc -= 0.5 * dx;
         yc -= 0.5 * dy;
         break;
      }
   }

   glTranslated(x, y, 0.);
   glRotated(angle, 0., 0., 1.);
   glTranslated(xc, yc, 0.);
   glTranslated(-0.5 * dx, -0.5 * dy, 0.);
   //glScaled(mgn, mgn, 1.);

   const_cast<FTFont*>(fFont)->Render(txt);

   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLFont::Render(const wchar_t* txt, Double_t x, Double_t y, Double_t angle, Double_t mgn) const
{
   RenderHelper(txt, x, y, angle, mgn);
}

////////////////////////////////////////////////////////////////////////////////

void TGLFont::Render(const char* txt, Double_t x, Double_t y, Double_t angle, Double_t mgn) const
{
   RenderHelper(txt, x, y, angle, mgn);
}

////////////////////////////////////////////////////////////////////////////////
/// Render text.

void TGLFont::Render(const TString &txt) const
{
   Bool_t scaleDepth = (fMode == kExtrude && fDepth != 1.0f);

   if (scaleDepth) {
      glPushMatrix();
      // !!! 0.2*fSize is hard-coded in TGLFontManager::GetFont(), too.
      glTranslatef(0.0f, 0.0f, 0.5f*fDepth * 0.2f*fSize);
      glScalef(1.0f, 1.0f, fDepth);
   }

   // FTGL is not const correct.
   const_cast<FTFont*>(fFont)->Render(txt);

   if (scaleDepth) {
      glPopMatrix();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Render text with given alignmentrepl and at given position.

void  TGLFont:: Render(const TString &txt, Float_t x, Float_t y, Float_t z,
             ETextAlignH_e alignH, ETextAlignV_e alignV) const
{
   glPushMatrix();

   glTranslatef(x, y, z);

   x=0, y=0;
   Float_t llx, lly, llz, urx, ury, urz;
   BBox(txt, llx, lly, llz, urx, ury, urz);

   switch (alignH)
   {
      case TGLFont::kRight:
         x = -urx;
         break;

      case  TGLFont::kCenterH:
         x = -urx*0.5;
         break;
      default:
         break;
   };

   switch (alignV)
   {
      case TGLFont::kBottom:
         y = -ury;
         break;
      case  TGLFont::kCenterV:
         y = -ury*0.5;
         break;
      default:
         break;
   };

   if (fMode == TGLFont::kPixmap || fMode ==  TGLFont::kBitmap)
   {
      glRasterPos2i(0, 0);
      glBitmap(0, 0, 0, 0, x, y, 0);
   }
   else
   {
      glTranslatef(x, y, 0);
   }
   Render(txt);
   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Set-up GL state before FTFont rendering.

void TGLFont::PreRender(Bool_t autoLight, Bool_t lightOn) const
{
   switch (fMode)
   {
      case kBitmap:
      case kPixmap:
         glPushAttrib(GL_CURRENT_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
         glEnable(GL_ALPHA_TEST);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
         glAlphaFunc(GL_GEQUAL, 0.0625);
         break;
      case kTexture:
         glPushAttrib(GL_POLYGON_BIT | GL_ENABLE_BIT);
         glEnable(GL_TEXTURE_2D);
         glDisable(GL_CULL_FACE);
         glEnable(GL_ALPHA_TEST);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
         glAlphaFunc(GL_GEQUAL, 0.0625);
         break;
      case kOutline:
      case kPolygon:
      case kExtrude:
         glPushAttrib(GL_POLYGON_BIT | GL_ENABLE_BIT);
         glEnable(GL_NORMALIZE);
         glDisable(GL_CULL_FACE);
         break;
      default:
         Warning("TGLFont::PreRender", "Font mode undefined.");
         glPushAttrib(GL_LIGHTING_BIT);
         break;
   }

   if ((autoLight && fMode > TGLFont::kOutline) || (!autoLight && lightOn))
      glEnable(GL_LIGHTING);
   else
      glDisable(GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset GL state after FTFont rendering.

void TGLFont::PostRender() const
{
   glPopAttrib();
}

/** \class TGLFontManager
\ingroup opengl
A FreeType GL font manager.

Each GL rendering context has an instance of FTGLManager.
This enables FTGL fonts to be shared same way as textures and display lists.
*/

ClassImp(TGLFontManager);

TObjArray   TGLFontManager::fgFontFileArray;
Int_t TGLFontManager::fgExtendedFontStart;
TGLFontManager::FontSizeVec_t TGLFontManager::fgFontSizeArray;
Bool_t  TGLFontManager::fgStaticInitDone = kFALSE;

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLFontManager::~TGLFontManager()
{
   FontMap_i it = fFontMap.begin();
   while (it != fFontMap.end()) {
      delete it->first.GetFont();
      it++;
   }
   fFontMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Provide font with given size, file and FTGL class.

void TGLFontManager::RegisterFont(Int_t sizeIn, Int_t fileID, TGLFont::EMode mode, TGLFont &out)
{
   if (fgStaticInitDone == kFALSE) InitStatics();

   Int_t  size = GetFontSize(sizeIn);
   if (mode == out.GetMode() && fileID == out.GetFile() && size == out.GetSize())
      return;

   FontMap_i it = fFontMap.find(TGLFont(size, fileID, mode));
   if (it == fFontMap.end())
   {
      TString ttpath, file;
      ttpath = gEnv->GetValue("Root.TTGLFontPath", TROOT::GetTTFFontDir());
      {
         //For extenede we have both ttf and otf.
         char *fp = gSystem->Which(ttpath, fileID < fgExtendedFontStart ?
                                   ((TObjString*)fgFontFileArray[fileID])->String() + ".ttf" :
                                   ((TObjString*)fgFontFileArray[fileID])->String());
         file = fp;
         delete [] fp;
      }

      FTFont* ftfont = 0;
      switch (mode)
      {
         case TGLFont::kBitmap:
            ftfont = new FTGLBitmapFont(file);
            break;
         case TGLFont::kPixmap:
            ftfont = new FTGLPixmapFont(file);
            break;
         case TGLFont::kOutline:
            ftfont = new FTGLOutlineFont(file);
            break;
         case TGLFont::kPolygon:
            ftfont = new FTGLPolygonFont(file);
            break;
         case TGLFont::kExtrude:
            ftfont = new FTGLExtrdFont(file);
            ftfont->Depth(0.2*size);
            break;
         case TGLFont::kTexture:
            ftfont = new FTGLTextureFont(file);
            break;
         default:
            Error("TGLFontManager::RegisterFont", "invalid FTGL type");
            return;
            break;
      }
      ftfont->FaceSize(size);
      const TGLFont &mf = fFontMap.insert(std::make_pair(TGLFont(size, fileID, mode, ftfont, 0), 1)).first->first;
      out.CopyAttributes(mf);
   }
   else
   {
      if (it->first.GetTrashCount() > 0) {
         fFontTrash.remove(&(it->first));
         it->first.SetTrashCount(0);
      }
      ++(it->second);
      out.CopyAttributes(it->first);
   }
   out.SetManager(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Get mapping from ttf id to font names. Table taken from TTF.cxx.

void TGLFontManager::RegisterFont(Int_t size, const char* name, TGLFont::EMode mode, TGLFont &out)
{
   TObjArray* farr = GetFontFileArray();
   TIter next(farr);
   TObjString* os;
   Int_t cnt = 0;
   while ((os = (TObjString*) next()) != 0)
   {
      if (os->String() == name)
         break;
      cnt++;
   }

   if (cnt < farr->GetEntries())
      RegisterFont(size, cnt, mode, out);
   else
      Error("TGLFontManager::RegisterFont", "unknown font name %s", name);
}

////////////////////////////////////////////////////////////////////////////////
/// Release font with given attributes. Returns false if font has
/// not been found in the managers font set.

void TGLFontManager::ReleaseFont(TGLFont& font)
{
   FontMap_i it = fFontMap.find(font);

   if (it != fFontMap.end())
   {
      --(it->second);
      if (it->second == 0)
      {
         assert(it->first.GetTrashCount() == 0);
         it->first.IncTrashCount();
         fFontTrash.push_back(&it->first);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get id to file name map.

TObjArray* TGLFontManager::GetFontFileArray()
{
   if (fgStaticInitDone == kFALSE) InitStatics();
   return &fgFontFileArray;
}

////////////////////////////////////////////////////////////////////////////////
/// Get valid font size vector.

TGLFontManager::FontSizeVec_t* TGLFontManager::GetFontSizeArray()
{
   if (fgStaticInitDone == kFALSE) InitStatics();
   return &fgFontSizeArray;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TGLFontManager::GetExtendedFontStartIndex()
{
   if (fgStaticInitDone == kFALSE) InitStatics();

   assert(fgExtendedFontStart > 0 && "GetExtendedFontStartIndex, invalid index");

   return fgExtendedFontStart;
}

////////////////////////////////////////////////////////////////////////////////
/// Get availabe font size.

Int_t TGLFontManager::GetFontSize(Int_t ds)
{
   if (fgStaticInitDone == kFALSE) InitStatics();

   Int_t idx = TMath::BinarySearch(fgFontSizeArray.size(), &fgFontSizeArray[0],
                                   TMath::CeilNint(ds));

   if (idx < 0) idx = 0;
   return fgFontSizeArray[idx];
}

////////////////////////////////////////////////////////////////////////////////
/// Get availabe font size.

Int_t TGLFontManager::GetFontSize(Int_t ds, Int_t min, Int_t max)
{
   if (ds < min) ds = min;
   if (ds > max) ds = max;
   return GetFontSize(ds);
}

////////////////////////////////////////////////////////////////////////////////
/// Get font name from TAttAxis font id.

const char* TGLFontManager::GetFontNameFromId(Int_t id)
{
   if (fgStaticInitDone == kFALSE) InitStatics();

   Int_t fontIndex = id / 10;

   if (fontIndex > fgFontFileArray.GetEntries() || !fontIndex)
      fontIndex = 5;//arialbd
   else
      fontIndex -= 1;

   TObjString* os = (TObjString*)fgFontFileArray[fontIndex];
   return os->String().Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a list of available font files and allowed font sizes.

void TGLFontManager::InitStatics()
{
   fgFontFileArray.Add(new TObjString("timesi"));   //  10
   fgFontFileArray.Add(new TObjString("timesbd"));  //  20
   fgFontFileArray.Add(new TObjString("timesbi"));  //  30

   fgFontFileArray.Add(new TObjString("arial"));    //  40
   fgFontFileArray.Add(new TObjString("ariali"));   //  50
   fgFontFileArray.Add(new TObjString("arialbd"));  //  60
   fgFontFileArray.Add(new TObjString("arialbi"));  //  70

   fgFontFileArray.Add(new TObjString("cour"));     //  80
   fgFontFileArray.Add(new TObjString("couri"));    //  90
   fgFontFileArray.Add(new TObjString("courbd"));   // 100
   fgFontFileArray.Add(new TObjString("courbi"));   // 110

   fgFontFileArray.Add(new TObjString("symbol"));   // 120
   fgFontFileArray.Add(new TObjString("times"));    // 130
   fgFontFileArray.Add(new TObjString("wingding")); // 140
   fgFontFileArray.Add(new TObjString("symbol"));   // 150

   fgExtendedFontStart = fgFontFileArray.GetEntries();
   //"Extended" fonts for gl-pad.
   //fgPadFontStart + ...
   fgFontFileArray.Add(new TObjString("FreeSerifItalic.otf"));      //  10 (160)
   fgFontFileArray.Add(new TObjString("FreeSerifBold.otf"));        //  20 (170)
   fgFontFileArray.Add(new TObjString("FreeSerifBoldItalic.otf"));  //  30

   fgFontFileArray.Add(new TObjString("FreeSans.otf"));             //  40
   fgFontFileArray.Add(new TObjString("FreeSansOblique.otf"));      //  50
   fgFontFileArray.Add(new TObjString("FreeSansBold.otf"));         //  60
   fgFontFileArray.Add(new TObjString("FreeSansBoldOblique.otf"));  //  70

   fgFontFileArray.Add(new TObjString("FreeMono.otf"));             //  80
   fgFontFileArray.Add(new TObjString("FreeMonoOblique.otf"));      //  90
   fgFontFileArray.Add(new TObjString("FreeMonoBold.otf"));         // 100
   fgFontFileArray.Add(new TObjString("FreeMonoBoldOblique.otf"));  // 110

   fgFontFileArray.Add(new TObjString("symbol.ttf"));               // 120
   fgFontFileArray.Add(new TObjString("FreeSerif.otf"));            // 130
   fgFontFileArray.Add(new TObjString("wingding.ttf"));             // 140
   fgFontFileArray.Add(new TObjString("symbol.ttf"));               // 150

   fgFontFileArray.Add(new TObjString("STIXGeneral.otf"));          // 200
   fgFontFileArray.Add(new TObjString("STIXGeneralItalic.otf"));    // 210
   fgFontFileArray.Add(new TObjString("STIXGeneralBol.otf"));       // 220
   fgFontFileArray.Add(new TObjString("STIXGeneralBolIta.otf"));    // 230

   fgFontFileArray.Add(new TObjString("STIXSiz1Sym.otf"));          // 240
   fgFontFileArray.Add(new TObjString("STIXSiz1SymBol.otf"));       // 250
   fgFontFileArray.Add(new TObjString("STIXSiz2Sym.otf"));          // 260
   fgFontFileArray.Add(new TObjString("STIXSiz2SymBol.otf"));       // 270

   fgFontFileArray.Add(new TObjString("STIXSiz3Sym.otf"));          // 280
   fgFontFileArray.Add(new TObjString("STIXSiz3SymBol.otf"));       // 290
   fgFontFileArray.Add(new TObjString("STIXSiz4Sym.otf"));          // 300
   fgFontFileArray.Add(new TObjString("STIXSiz4SymBol.otf"));       // 310

   fgFontFileArray.Add(new TObjString("STIXSiz5Sym.otf"));          // 320
   fgFontFileArray.Add(new TObjString("DroidSansFallback.ttf"));    // 330
   fgFontFileArray.Add(new TObjString("DroidSansFallback.ttf"));    // 340
   fgFontFileArray.Add(new TObjString("DroidSansFallback.ttf"));    // 350

   for (Int_t i = 10; i <= 20; i+=2)
      fgFontSizeArray.push_back(i);
   for (Int_t i = 24; i <= 64; i+=4)
      fgFontSizeArray.push_back(i);
   for (Int_t i = 72; i <= 128; i+=8)
      fgFontSizeArray.push_back(i);

   fgStaticInitDone = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete FTFFont objects registered for destruction.

void TGLFontManager::ClearFontTrash()
{
   FontList_i it = fFontTrash.begin();
   while (it != fFontTrash.end())
   {
      if ((*it)->IncTrashCount() > 10000)
      {
         FontMap_i mi = fFontMap.find(**it);
         assert(mi != fFontMap.end());
         fFontMap.erase(mi);
         delete (*it)->GetFont();

         FontList_i li = it++;
         fFontTrash.erase(li);
      }
      else
      {
         ++it;
      }
   }
}
