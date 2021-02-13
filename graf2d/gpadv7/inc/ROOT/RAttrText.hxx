/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrText
#define ROOT7_RAttrText

#include <ROOT/RAttrBase.hxx>
#include <ROOT/RAttrValue.hxx>
#include <ROOT/RAttrColor.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrText
\ingroup GpadROOT7
\brief A text attributes.
\authors Axel Naumann <axel@cern.ch> Sergey Linev <s.linev@gsi.de>
\date 2018-10-12
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrText : public RAttrBase {

   RAttrColor fColor{this, "color"};                          ///<! text color
   RAttrValue<double> fSize{this, "size", 12.};               ///<! text size
   RAttrValue<double> fAngle{this, "angle", 0.};              ///<! text angle
   RAttrValue<int> fAlign{this, "align", 22};                 ///<! text align
   RAttrValue<std::string> fFontFamily{this, "font_familty"}; ///<! font family, corresponds to css font-familty attribute
   RAttrValue<std::string> fFontStyle{this, "font_style"};    ///<! font style, corresponds to css font-style attribute
   RAttrValue<std::string> fFontWeight{this, "font_weight"};  ///<! font weight, corresponds to css font-weight attribute

   R__ATTR_CLASS(RAttrText, "text");

   ///The text size
   RAttrText &SetSize(double sz) { fSize = sz; return *this; }
   double GetSize() const { return fSize; }

   ///The text angle
   RAttrText &SetAngle(double angle) { fAngle = angle; return *this; }
   double GetAngle() const { return fAngle; }

   ///The text alignment
   RAttrText &SetAlign(int align) { fAlign = align; return *this; }
   int GetAlign() const { return fAlign; }

   ///Set text font by id as usually handled in the ROOT, set number between 1 and 15
   RAttrText &SetFont(int font)
   {
      std::string family, style, weight;
      switch(font) {
      case 1: family = "Times New Roman"; style = "italic"; break;
      case 2: family = "Times New Roman"; weight = "bold"; break;
      case 3: family = "Times New Roman"; style = "italic"; weight = "bold"; break;
      case 4: family = "Arial"; break;
      case 5: family = "Arial"; style = "oblique"; break;
      case 6: family = "Arial"; weight = "bold"; break;
      case 7: family = "Arial"; style = "oblique"; weight = "bold"; break;
      case 8: family = "Courier New"; break;
      case 9: family = "Courier New"; style = "oblique"; break;
      case 10: family = "Courier New"; weight = "bold"; break;
      case 11: family = "Courier New"; style = "oblique"; weight = "bold"; break;
      case 12: family = "Symbol"; break;
      case 13: family = "Times New Roman"; break;
      case 14: family = "Wingdings"; break;
      case 15: family = "Symbol"; style = "italic"; break;
      }

      SetFontFamily(family);
      SetFontWeight(style);
      SetFontStyle(weight);

      return *this;
   }

   ///The text font-family attribute
   RAttrText &SetFontFamily(const std::string &family)
   {
      if (family.empty())
         fFontFamily.Clear();
      else
         fFontFamily = family;
      return *this;
   }
   std::string GetFontFamily() const { return fFontFamily; }

   ///The text font-style attribute
   RAttrText &SetFontStyle(const std::string &style)
   {
      if (style.empty())
         fFontStyle.Clear();
      else
         fFontStyle = style;
      return *this;
   }
   std::string GetFontStyle() const { return fFontStyle; }

   ///The text font-weight attribute
   RAttrText &SetFontWeight(const std::string &weight)
   {
      if (weight.empty())
         fFontWeight.Clear();
      else
         fFontWeight = weight;
      return *this;
   }
   std::string GetFontWeight() const { return fFontWeight; }


   ///The color of the text.
   RAttrText &SetColor(const RColor &color) { fColor = color; return *this; }
   RColor GetColor() const { return fColor.GetColor(); }
   RAttrColor &AttrColor() { return fColor; }

};



} // namespace Experimental
} // namespace ROOT

#endif
