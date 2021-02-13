/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RLevelIter.hxx>

#include <ROOT/Browsable/RElement.hxx>
#include <ROOT/Browsable/RItem.hxx>

using namespace ROOT::Experimental::Browsable;

/////////////////////////////////////////////////////////////////////
/// Find item with specified name
/// Default implementation, should work for all
/// If index specified, not only name but also index should match

bool RLevelIter::Find(const std::string &name, int indx)
{
   int i = -1;

   while (Next()) {
      if (indx >= 0) {
         i++;
         if (i > indx) return false;
         if (i < indx) continue;
      }

      if (GetItemName() == name)
         return true;
   }

   return false;
}

/////////////////////////////////////////////////////////////////////
/// Create generic description item for RBrowser

std::unique_ptr<RItem> RLevelIter::CreateItem()
{
   auto have_childs = CanItemHaveChilds();

   return std::make_unique<RItem>(GetItemName(), have_childs ? -1 : 0, have_childs ? "sap-icon://folder-blank" : "sap-icon://document");
}

