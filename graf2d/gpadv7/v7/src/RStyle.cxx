/// \file RStyle.cxx
/// \ingroup Gpad ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \author Sergey Linev <s.linev@gsi.de>
/// \date 2017-10-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RStyle.hxx>

#include <ROOT/RDrawable.hxx>

// #include "RStyleReader.hxx" // in src/


using namespace std::string_literals;


///////////////////////////////////////////////////////////////////////////////
/// Evaluate style

const ROOT::Experimental::RAttrMap::Value_t *ROOT::Experimental::RStyle::Eval(const std::string &field, const RDrawable *drawable) const
{
   for (const auto &block : fBlocks) {

      bool match = !drawable || drawable->MatchSelector(block.selector);

      if (match) {
         auto res = block.map.Find(field);
         if (res) return res;
      }
   }

   return nullptr;
}
