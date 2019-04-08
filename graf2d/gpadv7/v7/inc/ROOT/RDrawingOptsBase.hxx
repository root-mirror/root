/// \file ROOT/RDrawingOptionsBase.hxx
/// \ingroup Gpad ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2018-02-12
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RDrawingOptsBase
#define ROOT7_RDrawingOptsBase

#include <ROOT/RDrawingAttr.hxx>

#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {

class RDrawingOptsBase: public RDrawingAttrHolderBase {
   /// Attribute style classes of these options.
   std::vector<std::string> fStyleClasses;

public:
   RDrawingOptsBase() = default;

   /// Initialize the options with a (possibly empty) set of style classes.
   RDrawingOptsBase(std::vector<std::string> &styleClasses): fStyleClasses(styleClasses) {}

   /// Get the attribute style classes of these options.
   const std::vector<std::string> &GetStyleClasses() const { return fStyleClasses; }

   std::string GetAttrFromStyle(const Name_t &attrName) override;
};

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RDrawingOptsBase
