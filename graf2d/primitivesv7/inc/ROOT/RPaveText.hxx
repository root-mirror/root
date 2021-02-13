/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPaveText
#define ROOT7_RPaveText

#include <ROOT/RPave.hxx>

namespace ROOT {
namespace Experimental {

/** \class RPaveText
\ingroup GrafROOT7
\brief A RPave with text content
\author Sergey Linev <S.Linev@gsi.de>
\date 2020-06-19
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is
welcome!
*/

class RPaveText : public RPave {

   std::vector<std::string> fText; ///< list of text entries

public:
   RPaveText() : RPave("pavetext") {}

   void AddLine(const std::string &txt) { fText.emplace_back(txt); }

   auto NumLines() const { return fText.size(); }

   const std::string &GetLine(int n) const { return fText[n]; }

   void ClearLines() { fText.clear(); }
};

} // namespace Experimental
} // namespace ROOT

#endif
