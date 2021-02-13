/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RText
#define ROOT7_RText

#include <ROOT/RDrawable.hxx>
#include <ROOT/RAttrText.hxx>
#include <ROOT/RPadPos.hxx>

#include <string>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::RText
\ingroup GrafROOT7
\brief A text.
\author Olivier Couet <Olivier.Couet@cern.ch>
\date 2017-10-16
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is
welcome!
*/

class RText : public RDrawable {

   std::string fText;                 ///< text to display
   RPadPos fPos;                      ///< position
   RAttrText fAttrText{this, "text"}; ///<! text attributes

public:
   RText() : RDrawable("text") {}

   RText(const std::string &txt) : RText() { fText = txt; }

   RText(const RPadPos &p, const std::string &txt) : RText()
   {
      fText = txt;
      fPos = p;
   }

   RText &SetText(const std::string &t)
   {
      fText = t;
      return *this;
   }
   const std::string &GetText() const { return fText; }

   RText &SetPos(const RPadPos &p)
   {
      fPos = p;
      return *this;
   }
   const RPadPos &GetPos() const { return fPos; }

   const RAttrText &GetAttrText() const { return fAttrText; }
   RText &SetAttrText(const RAttrText &attr)
   {
      fAttrText = attr;
      return *this;
   }
   RAttrText &AttrText() { return fAttrText; }
};

} // namespace Experimental
} // namespace ROOT

#endif
