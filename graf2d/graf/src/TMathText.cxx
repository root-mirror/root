// @(#)root/graf:$Id: TMathText.cxx  $
// Author: Yue Shi Lai  16/10/12

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TROOT.h"
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H
#include "TTF.h"
#include "TMathText.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TVirtualPS.h"
#include "TVirtualX.h"
#include "TText.h"

#include "../../../graf2d/mathtext/inc/mathtext.h"
#include "../../../graf2d/mathtext/inc/mathrender.h"

/** \class TMathText
\ingroup BasicGraphics

To draw TeX Mathematical Formula

TMathText's purpose is to write mathematical equations, exactly as TeX would
do it. The syntax is the same as the TeX's one.

The following example demonstrate how to use TMathText:

Begin_Macro(source)
../../../tutorials/graphics/tmathtext.C
End_Macro

The list of all available symbols is given in the following example:

Begin_Macro
../../../tutorials/graphics/tmathtext2.C
End_Macro

#### Limitation:
TMathText rendering is not implemented for the PDF output.
PostScript output should be used instead.
*/

const Double_t kPI      = TMath::Pi();

class TMathTextRenderer : public TText, public TAttFill,
                          public mathtext::math_text_renderer_t {
private:
   TMathText *_parent;
   float _font_size;
   float _x0;
   float _y0;
   float _angle_degree;
   float _pad_pixel_transform[6];
   float _pad_scale;
   float _pad_scale_x;
   float _pad_scale_y;
   float _pad_scale_x_relative;
   float _pad_scale_y_relative;
   float _current_font_size[mathtext::math_text_renderer_t::NFAMILY];
   inline size_t root_face_number(
      const unsigned int family, const bool serif = false) const
   {
      static const int precision = 2;

      if (family >= mathtext::math_text_renderer_t::
         FAMILY_REGULAR &&
         family <= mathtext::math_text_renderer_t::
         FAMILY_BOLD_ITALIC) {
         const unsigned int offset = family -
            mathtext::math_text_renderer_t::FAMILY_REGULAR;
         return serif ?
            ((offset == 0 ? 13 : offset) * 10 + precision) :
            ((offset + 4) * 10 + precision);
      } else if (family >= mathtext::math_text_renderer_t::
            FAMILY_STIX_REGULAR) {
         const unsigned int offset = family -
            mathtext::math_text_renderer_t::FAMILY_STIX_REGULAR;
         return (offset + 16) * 10 + precision;
      }

      return precision;
   }
   inline bool is_cyrillic_or_cjk(const wchar_t c) const
   {
      return mathtext::math_text_renderer_t::is_cyrillic(c) ||
         mathtext::math_text_renderer_t::is_cjk(c);
   }
   inline size_t root_cjk_face_number(
      const bool serif = false) const
   {
      return (serif ? 28 : 29) * 10 + 2;
   }
protected:
   inline mathtext::affine_transform_t
   transform_logical_to_pixel(void) const
   {
      return mathtext::affine_transform_t::identity;
   }
   inline mathtext::affine_transform_t
   transform_pixel_to_logical(void) const
   {
      return mathtext::affine_transform_t::identity;
   }
public:
   inline TMathTextRenderer(TMathText *parent)
      : TText(), TAttFill(0, 1001),
        _parent(parent), _font_size(0), _angle_degree(0)
   {
      int i;
      _font_size = 0;
      _x0 = 0;
      _y0 = 0;
      _angle_degree = 0;
      for (i = 0; i<6; i++) _pad_pixel_transform[i] = 0;
      _pad_scale = 0;
      _pad_scale_x = 0;
      _pad_scale_y = 0;
      _pad_scale_x_relative = 0;
      _pad_scale_y_relative = 0;
      for (i = 0; i < mathtext::math_text_renderer_t::NFAMILY; i++) _current_font_size[i] = 0;
   }
   inline float
   font_size(const unsigned int family = FAMILY_PLAIN) const
   {
      return _current_font_size[family];
   }
   inline void
   point(const float /*x*/, const float /*y*/)
   {
   }
   inline void
   set_font_size(const float size, const unsigned int family)
   {
      _current_font_size[family] = size;
   }
   inline void
   set_font_size(const float size)
   {
      _font_size = size;
      std::fill(_current_font_size,
              _current_font_size + NFAMILY, size);
   }
   inline void
   reset_font_size(const unsigned int /*family*/)
   {
   }
   inline void
   set_parameter(const float x, const float y, const float size,
              const float angle_degree)
   {
      _x0 = gPad->XtoAbsPixel(x);
      _y0 = gPad->YtoAbsPixel(y);
      _pad_scale_x =
         gPad->XtoPixel(gPad->GetX2()) -
         gPad->XtoPixel(gPad->GetX1());
      _pad_scale_y =
         gPad->YtoPixel(gPad->GetY1()) -
         gPad->YtoPixel(gPad->GetY2());
      _pad_scale = std::min(_pad_scale_x, _pad_scale_y);

      _angle_degree = angle_degree;

      const float angle_radiant = _angle_degree * (kPI / 180.0);

      // Initialize the affine transform
      _pad_pixel_transform[0] = _pad_scale * cosf(angle_radiant);
      _pad_pixel_transform[1] = -_pad_scale * sinf(angle_radiant);
      _pad_pixel_transform[2] = _x0;
      _pad_pixel_transform[3] = _pad_pixel_transform[1];
      _pad_pixel_transform[4] = -_pad_pixel_transform[0];
      _pad_pixel_transform[5] = _y0;

      set_font_size(size);
      SetTextAngle(_angle_degree);
      SetTextColor(_parent->fTextColor);
   }
   inline void
   transform_pad(double &xt, double &yt,
              const float x, const float y) const
   {
      xt = gPad->AbsPixeltoX(Int_t(
         x * _pad_pixel_transform[0] +
         y * _pad_pixel_transform[1] + _pad_pixel_transform[2]));
      yt = gPad->AbsPixeltoY(Int_t(
         x * _pad_pixel_transform[3] +
         y * _pad_pixel_transform[4] + _pad_pixel_transform[5]));
   }
   inline void
   filled_rectangle(const mathtext::bounding_box_t &bounding_box_0)
   {
      SetFillColor(_parent->fTextColor);
      SetFillStyle(1001);
      TAttFill::Modify();

      double xt[4];
      double yt[4];

      transform_pad(xt[0], yt[0],
                 bounding_box_0.left(),
                 bounding_box_0.bottom());
      transform_pad(xt[1], yt[1],
                 bounding_box_0.right(),
                 bounding_box_0.bottom());
      transform_pad(xt[2], yt[2],
                 bounding_box_0.right(),
                 bounding_box_0.top());
      transform_pad(xt[3], yt[3],
                 bounding_box_0.left(),
                 bounding_box_0.top());
      gPad->PaintFillArea(4, xt, yt);
   }
   inline void
   rectangle(const mathtext::bounding_box_t &/*bounding_box*/)
   {
   }
   inline mathtext::bounding_box_t
   bounding_box(const wchar_t character, float &current_x,
             const unsigned int family)
   {
      const size_t old_font_index = TTF::fgCurFontIdx;
      const bool cyrillic_or_cjk = is_cyrillic_or_cjk(character);

      if (cyrillic_or_cjk) {
         TTF::SetTextFont((Font_t) root_cjk_face_number());
      } else {
         TTF::SetTextFont((Font_t) root_face_number(family));
      }
      FT_Load_Glyph(
         TTF::fgFace[TTF::fgCurFontIdx],
         FT_Get_Char_Index(
            TTF::fgFace[TTF::fgCurFontIdx], character),
         FT_LOAD_NO_SCALE);

      const float scale = _current_font_size[family] /
         TTF::fgFace[TTF::fgCurFontIdx]->units_per_EM;
      const FT_Glyph_Metrics metrics =
         TTF::fgFace[TTF::fgCurFontIdx]->glyph->metrics;
      const float lower_left_x = metrics.horiBearingX;
      const float lower_left_y =
         metrics.horiBearingY - metrics.height;
      const float upper_right_x =
         metrics.horiBearingX + metrics.width;
      const float upper_right_y = metrics.horiBearingY;
      const float advance = metrics.horiAdvance;
      const float margin = std::max(0.0F, lower_left_x);
      const float italic_correction =
         upper_right_x <= advance ? 0.0F :
         std::max(0.0F, upper_right_x + margin - advance);
      const mathtext::bounding_box_t ret =
         mathtext::bounding_box_t(
            lower_left_x, lower_left_y,
            upper_right_x, upper_right_y,
            advance, italic_correction) * scale;

      current_x += ret.advance();
      TTF::fgCurFontIdx = old_font_index;

      return ret;
   }
   inline mathtext::bounding_box_t
   bounding_box(const std::wstring string,
             const unsigned int family = FAMILY_PLAIN)
   {
      if (TTF::fgCurFontIdx<0) return mathtext::bounding_box_t(0, 0, 0, 0, 0, 0);
      if (string.empty() || TTF::fgFace[TTF::fgCurFontIdx] == NULL ||
         TTF::fgFace[TTF::fgCurFontIdx]->units_per_EM == 0) {
         return mathtext::bounding_box_t(0, 0, 0, 0, 0, 0);
      }

      std::wstring::const_iterator iterator = string.begin();
      float current_x = 0;
      mathtext::bounding_box_t ret =
         bounding_box(*iterator, current_x, family);

      ++iterator;
      for (; iterator != string.end(); ++iterator) {
         const mathtext::point_t position =
            mathtext::point_t(current_x, 0);
         const mathtext::bounding_box_t glyph_bounding_box =
            bounding_box(*iterator, current_x, family);
         ret = ret.merge(position + glyph_bounding_box);
      }

      return ret;
   }
   inline void
   text_raw(const float x, const float y,
          const std::wstring string,
          const unsigned int family = FAMILY_PLAIN)
   {
      SetTextFont((Font_t) root_face_number(family));
      SetTextSize(_current_font_size[family]);
      TAttText::Modify();

      wchar_t buf[2];
      float advance = 0;

      buf[1] = L'\0';
      for (std::wstring::const_iterator iterator = string.begin(); iterator != string.end(); ++iterator) {
         buf[0] = *iterator;
         const bool cyrillic_or_cjk = is_cyrillic_or_cjk(buf[0]);

         if (cyrillic_or_cjk) {
            SetTextFont((Font_t) root_cjk_face_number());
            TAttText::Modify();
         }

         const mathtext::bounding_box_t b =
            bounding_box(buf, family);
         double xt;
         double yt;

         transform_pad(xt, yt, x + advance, y);
         gPad->PaintText(xt, yt, buf);
         advance += b.advance();
         if (cyrillic_or_cjk) {
            SetTextFont((Font_t) root_face_number(family));
            TAttText::Modify();
         }
      }
   }
   inline void
   text_with_bounding_box(const float /*x*/, const float /*y*/,
                     const std::wstring /*string*/,
                     const unsigned int /*family = FAMILY_PLAIN*/)
   {
   }
   using mathtext::math_text_renderer_t::bounding_box;
};

ClassImp(TMathText);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TMathText::TMathText(void)
   : TAttFill(0, 1001)
{
   fRenderer = new TMathTextRenderer(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TMathText::TMathText(Double_t x, Double_t y, const char *text)
   : TText(x, y, text), TAttFill(0, 1001)
{
   fRenderer = new TMathTextRenderer(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TMathText::~TMathText(void)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TMathText::TMathText(const TMathText &text)
   : TText(text), TAttFill(text)
{
   ((TMathText &)text).Copy(*this);
   fRenderer = new TMathTextRenderer(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TMathText &TMathText::operator=(const TMathText &rhs)
{
   if (this != &rhs) {
      TText::operator    = (rhs);
      TAttFill::operator = (rhs);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TMathText::Copy(TObject &obj) const
{
   ((TMathText &)obj).fRenderer = fRenderer;
   TText::Copy(obj);
   TAttFill::Copy((TAttFill &)obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Render the text.

void TMathText::
Render(const Double_t x, const Double_t y, const Double_t size,
      const Double_t angle, const Char_t *t, const Int_t /*length*/)
{
   const mathtext::math_text_t math_text(t);
   TMathTextRenderer *renderer = (TMathTextRenderer *)fRenderer;

   renderer->set_parameter(x, y, size, angle);
   renderer->text(0, 0, math_text);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the text bounding box.

void TMathText::
GetSize(Double_t &x0, Double_t &y0, Double_t &x1, Double_t &y1,
      const Double_t size, const Double_t angle, const Char_t *t,
      const Int_t /*length*/)
{
   const mathtext::math_text_t math_text(t);
   TMathTextRenderer *renderer = (TMathTextRenderer *)fRenderer;

   renderer->set_parameter(0, 0, size, angle);

   const mathtext::bounding_box_t bounding_box =
      renderer->bounding_box(math_text);
   double x[4];
   double y[4];

   renderer->transform_pad(
      x[0], y[0], bounding_box.left(), bounding_box.bottom());
   renderer->transform_pad(
      x[1], y[1], bounding_box.right(), bounding_box.bottom());
   renderer->transform_pad(
      x[2], y[2], bounding_box.right(), bounding_box.top());
   renderer->transform_pad(
      x[3], y[3], bounding_box.left(), bounding_box.top());

   x0 = std::min(std::min(x[0], x[1]), std::min(x[2], x[3]));
   y0 = std::min(std::min(y[0], y[1]), std::min(y[2], y[3]));
   x1 = std::max(std::max(x[0], x[1]), std::max(x[2], x[3]));
   y1 = std::max(std::max(y[0], y[1]), std::max(y[2], y[3]));
}

////////////////////////////////////////////////////////////////////////////////
/// Alignment.

void TMathText::
GetAlignPoint(Double_t &x0, Double_t &y0,
           const Double_t size, const Double_t angle,
           const Char_t *t, const Int_t /*length*/,
           const Short_t align)
{
   const mathtext::math_text_t math_text(t);
   TMathTextRenderer *renderer = (TMathTextRenderer *)fRenderer;

   renderer->set_parameter(0, 0, size, angle);

   const mathtext::bounding_box_t bounding_box =
      renderer->bounding_box(math_text);
   float x = 0;
   float y = 0;

   Short_t halign = align / 10;
   Short_t valign = align - 10 * halign;

   switch(halign) {
      case 0:   x = bounding_box.left();              break;
      case 1:   x = 0;                                break;
      case 2:   x = bounding_box.horizontal_center(); break;
      case 3:   x = bounding_box.right();             break;
   }
   switch(valign) {
      case 0:   y = bounding_box.bottom();            break;
      case 1:   y = 0;                                break;
      case 2:   y = bounding_box.vertical_center();   break;
      case 3:   y = bounding_box.top();               break;
   }
   renderer->transform_pad(x0, y0, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the text width and height.

void TMathText::GetBoundingBox(UInt_t &w, UInt_t &h, Bool_t /*angle*/)
{
   const TString newText = GetTitle();
   const Int_t length = newText.Length();
   const Char_t *text = newText.Data();
   const Double_t size = GetTextSize();

   Double_t x0;
   Double_t y0;
   Double_t x1;
   Double_t y1;

   GetSize(x0, y0, x1, y1, size, 0, text, length);
   w = (UInt_t)(TMath::Abs(gPad->XtoAbsPixel(x1) - gPad->XtoAbsPixel(x0)));
   h = (UInt_t)(TMath::Abs(gPad->YtoAbsPixel(y0) - gPad->YtoAbsPixel(y1)));
}

////////////////////////////////////////////////////////////////////////////////
/// Get X size.

Double_t TMathText::GetXsize(void)
{
   const TString newText = GetTitle();
   const Int_t length    = newText.Length();
   const Char_t *text    = newText.Data();
   const Double_t size   = GetTextSize();
   const Double_t angle  = GetTextAngle();

   Double_t x0;
   Double_t y0;
   Double_t x1;
   Double_t y1;

   GetSize(x0, y0, x1, y1, size, angle, text, length);

   return TMath::Abs(x1 - x0);
}

////////////////////////////////////////////////////////////////////////////////
/// Get Y size.

Double_t TMathText::GetYsize(void)
{
   const TString newText = GetTitle();
   const Int_t length    = newText.Length();
   const Char_t *text    = newText.Data();
   const Double_t size   = GetTextSize();
   const Double_t angle  = GetTextAngle();

   Double_t x0;
   Double_t y0;
   Double_t x1;
   Double_t y1;

   GetSize(x0, y0, x1, y1, size, angle, text, length);

   return TMath::Abs(y0 - y1);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this object with the new parameters
/// and copy object attributes.

TMathText *TMathText::DrawMathText(Double_t x, Double_t y, const char *text)
{
   TMathText *newtext = new TMathText(x, y, text);
   TAttText::Copy(*newtext);

   newtext->SetBit(kCanDelete);
   if (TestBit(kTextNDC)) newtext->SetNDC();
   newtext->AppendPad();

   return newtext;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint text.

void TMathText::Paint(Option_t *)
{
   Double_t xsave = fX;
   Double_t ysave = fY;

   if (TestBit(kTextNDC)) {
      fX = gPad->GetX1() + xsave * (gPad->GetX2() - gPad->GetX1());
      fY = gPad->GetY1() + ysave * (gPad->GetY2() - gPad->GetY1());
      PaintMathText(fX, fY, GetTextAngle(), GetTextSize(), GetTitle());
   } else {
      PaintMathText(gPad->XtoPad(fX), gPad->YtoPad(fY),
                 GetTextAngle(), GetTextSize(), GetTitle());
   }
   fX = xsave;
   fY = ysave;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint text (used by Paint()).

void TMathText::PaintMathText(Double_t x, Double_t y, Double_t angle,
                              Double_t size, const Char_t *text1)
{
   Double_t saveSize = size;
   Int_t saveFont    = fTextFont;
   Short_t saveAlign = fTextAlign;

   TAttText::Modify();
   if (gVirtualPS) { // Initialise TMathTextRenderer
      if (gPad->IsBatch()) {
         if (gVirtualPS->InheritsFrom("TImageDump")) gPad->PaintText(0, 0, "");
      }
   }

   // Do not use Latex if font is low precision.
   if (fTextFont % 10 < 2) {
      if (gVirtualX) {
         gVirtualX->SetTextAngle(angle);
      }
      if (gVirtualPS) {
         gVirtualPS->SetTextAngle(angle);
      }
      gPad->PaintText(x, y, text1);
      return;
   }

   if (fTextFont % 10 > 2) {
      UInt_t w = TMath::Abs(gPad->XtoAbsPixel(gPad->GetX2()) -
                       gPad->XtoAbsPixel(gPad->GetX1()));
      UInt_t h = TMath::Abs(gPad->YtoAbsPixel(gPad->GetY2()) -
                       gPad->YtoAbsPixel(gPad->GetY1()));
      size = size / std::min(w, h);
      SetTextFont(10 * (saveFont / 10) + 2);
   }

   TString newText = text1;

   if (newText.Length() == 0) return;

   // Compatibility with TLatex and Latex
   newText.ReplaceAll("\\omicron","o");
   newText.ReplaceAll("\\Alpha","A");
   newText.ReplaceAll("\\Beta","B");
   newText.ReplaceAll("\\Epsilon","E");
   newText.ReplaceAll("\\Zeta","Z");
   newText.ReplaceAll("\\Eta","H");
   newText.ReplaceAll("\\Iota","I");
   newText.ReplaceAll("\\Kappa","K");
   newText.ReplaceAll("\\Mu","M");
   newText.ReplaceAll("\\Nu","N");
   newText.ReplaceAll("\\Omicron","O");
   newText.ReplaceAll("\\Rho","P");
   newText.ReplaceAll("\\Tau","T");
   newText.ReplaceAll("\\Chi","X");
   newText.ReplaceAll("\\varomega","\\varpi");
   newText.ReplaceAll("\\mbox","\\hbox");
   newText.ReplaceAll("\\bar","\\wwbar");
   if (newText.Contains("\\frac")) {
      Int_t len,i1,i2;
      TString str;
      while (newText.Contains("\\frac")) {
         len = newText.Length();
         i1  = newText.Index("\\frac");
         str = newText(i1,len).Data();
         i2  = str.Index("}{");
         newText.Replace(i1+i2,2," \\over ");
         newText.Remove(i1,5);
      }
   }
   if (newText.Contains("\\splitline")) {
      Int_t len,i1,i2;
      TString str;
      while (newText.Contains("\\splitline")) {
         len = newText.Length();
         i1  = newText.Index("\\splitline");
         str = newText(i1,len).Data();
         i2  = str.Index("}{");
         newText.Replace(i1+i2,2," \\atop ");
         newText.Remove(i1,10);
      }
   }

   const Int_t length = newText.Length();
   const Char_t *text = newText.Data();
   Double_t x0;
   Double_t y0;
   GetAlignPoint(x0, y0, size, angle, text, length, fTextAlign);

   Render(x - x0, y - y0, size, angle, text, length);

   SetTextSize(saveSize);
   SetTextFont(saveFont);
   SetTextAlign(saveAlign);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TMathText::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   const char quote = '"';

   if (gROOT->ClassSaved(TMathText::Class())) {
      out << "   ";
   } else {
      out << "   TMathText *";
   }

   TString s = GetTitle();

   s.ReplaceAll("\\","\\\\");
   s.ReplaceAll("\"","\\\"");
   out << "mathtex = new TMathText("<< fX << "," << fY << ","
      << quote << s.Data() << quote << ");" << std::endl;
   if (TestBit(kTextNDC)) {
      out << "mathtex->SetNDC();" << std::endl;
   }

   SaveTextAttributes(out, "mathtex", 11, 0, 1, 42, 0.05);
   SaveFillAttributes(out, "mathtex", 0, 1001);

   out<<"   mathtex->Draw();" << std::endl;
}
