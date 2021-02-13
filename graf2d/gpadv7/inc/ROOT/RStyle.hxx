/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RStyle
#define ROOT7_RStyle

#include <ROOT/RAttrMap.hxx>

#include <string>
#include <list>
#include <memory>

namespace ROOT {
namespace Experimental {

class RDrawable;

/** \class RStyle
\ingroup GpadROOT7
\brief A set of defaults for graphics attributes, e.g. for histogram fill color, line width, frame offsets etc. Uses CSS syntax.
\author Axel Naumann <axel@cern.ch>
\author Sergey Linev <s.linev@gsi.de>
\date 2017-10-10
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RStyle {

public:

   struct Block_t {
      std::string selector;       ///<  css selector
      RAttrMap map;               ///<  container with attributes
      Block_t() = default;
      Block_t(const std::string &_selector) : selector(_selector) {}

      Block_t(const Block_t &src) : selector(src.selector), map(src.map) {} // not required
      Block_t& operator=(const Block_t &) = delete;
   };

   const RAttrMap::Value_t *Eval(const std::string &field, const RDrawable &drawable) const;

   const RAttrMap::Value_t *Eval(const std::string &field, const std::string &selector) const;

   RAttrMap &AddBlock(const std::string &selector)
   {
      fBlocks.emplace_back(selector);
      return fBlocks.back().map;
   }

   void Clear();

   static std::shared_ptr<RStyle> Parse(const std::string &css_code);

   bool ParseString(const std::string &css_code);

private:

   std::list<Block_t> fBlocks;  // use std::list to avoid calling of Block_t copy constructor

};

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RStyle
