// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveSecondarySelectable
#define ROOT7_REveSecondarySelectable

#include "Rtypes.h"

#include <set>

// XXXX class TGLSelectRecord;

namespace ROOT {
namespace Experimental {

class REveSecondarySelectable
{
private:
   REveSecondarySelectable(const REveSecondarySelectable &) = delete;
   REveSecondarySelectable &operator=(const REveSecondarySelectable &) = delete;

public:
   typedef std::set<Int_t>                SelectionSet_t;

protected:
   Bool_t fAlwaysSecSelect{kFALSE}; // Always do secondary-selection in GL.

   SelectionSet_t fSelectedSet;    // Selected indices.
   SelectionSet_t fHighlightedSet; // Highlighted indices.

   // XXXX
   // void ProcessGLSelectionInternal(TGLSelectRecord& rec, SelectionSet_t& sset);

public:
   REveSecondarySelectable() = default;
   virtual ~REveSecondarySelectable() {}

   Bool_t GetAlwaysSecSelect()   const { return fAlwaysSecSelect; }
   void   SetAlwaysSecSelect(Bool_t f) { fAlwaysSecSelect = f; }

   SelectionSet_t &RefSelectedSet()    { return fSelectedSet; }
   SelectionSet_t &RefHighlightedSet() { return fHighlightedSet; }

   // XXXX
   // void   ProcessGLSelection(TGLSelectRecord& rec);
};

} // namespace Experimental
} // namespace ROOT

#endif
