/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RGroup.hxx>

#include <ROOT/Browsable/RLevelIter.hxx>

#include <ROOT/RLogger.hxx>

using namespace ROOT::Experimental::Browsable;

class RGroupIter : public RLevelIter {
   int fIndx{-1};
   RGroup &fComp;

public:

   explicit RGroupIter(RGroup &comp) : fComp(comp) {}
   virtual ~RGroupIter() = default;

   /** Shift to next element */
   bool Next() override { return ++fIndx < (int) fComp.GetChilds().size(); }

   /** Returns current element name  */
   std::string GetItemName() const override { return fComp.GetChilds()[fIndx]->GetName(); }

   /** Returns true if item can have childs */
   bool CanItemHaveChilds() const override { return true; }

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override { return fComp.GetChilds()[fIndx]; }

   /** Find item with specified name, use item MatchName() functionality */
   bool Find(const std::string &name, int indx = -1) override
   {
      if ((indx >= 0) && (indx <= (int) fComp.GetChilds().size()))
         if (fComp.GetChilds()[indx]->MatchName(name)) {
            fIndx = indx;
            return true;
         }

      while (Next()) {
         if (fComp.GetChilds()[fIndx]->MatchName(name))
            return true;
      }

      return false;
   }

};

/////////////////////////////////////////////////////////////////////
/// Create iterator for childs of composite

std::unique_ptr<RLevelIter> RGroup::GetChildsIter()
{
   return std::make_unique<RGroupIter>(*this);
}

