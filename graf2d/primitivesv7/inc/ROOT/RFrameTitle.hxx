/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RFrameTitle
#define ROOT7_RFrameTitle

#include <ROOT/RDrawable.hxx>
#include <ROOT/RAttrText.hxx>
#include <ROOT/RAttrValue.hxx>
#include <ROOT/RPadPos.hxx>

#include <string>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::RFrameTitle
\ingroup GrafROOT7
\brief A title for the RFrame.
\author Sergey Linev <s.linev@gsi.de>
\date 2020-02-26
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is
welcome!
*/

class RFrameTitle final : public RDrawable {

   std::string fText;                                           ///< title to display
   RAttrText fAttrText{this, "text"};                           ///<! title text attributes
   RAttrValue<RPadLength> fMargin{this, "margin", 0.02_normal}; ///<! title margin
   RAttrValue<RPadLength> fHeight{this, "height", 0.05_normal}; ///<! title height

protected:
   bool IsFrameRequired() const final { return true; }

public:
   RFrameTitle() : RDrawable("title") {}

   RFrameTitle(const std::string &txt) : RFrameTitle() { fText = txt; }

   RFrameTitle &SetText(const std::string &t)
   {
      fText = t;
      return *this;
   }
   const std::string &GetText() const { return fText; }

   RFrameTitle &SetMargin(const RPadLength &pos)
   {
      fMargin = pos;
      return *this;
   }
   RPadLength GetMargin() const { return fMargin; }

   RFrameTitle &SetHeight(const RPadLength &pos)
   {
      fHeight = pos;
      return *this;
   }
   RPadLength GetHeight() const { return fHeight; }

   const RAttrText &GetAttrText() const { return fAttrText; }
   RFrameTitle &SetAttrText(const RAttrText &attr)
   {
      fAttrText = attr;
      return *this;
   }
   RAttrText &AttrText() { return fAttrText; }
};

} // namespace Experimental
} // namespace ROOT

#endif
