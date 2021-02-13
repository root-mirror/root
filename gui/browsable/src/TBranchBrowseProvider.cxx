/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/TObjectElement.hxx>
#include <ROOT/Browsable/RProvider.hxx>
#include <ROOT/Browsable/RLevelIter.hxx>

#include "TBranchElement.h"
#include "TBranchBrowsable.h"

using namespace ROOT::Experimental::Browsable;


////////////////////////////////////////////////////////////
/// Representing TBranchElement in browsables
/// Kept here only for a demo - default TObject-based API is enough for handling TBranchElement

class TBrElement : public TObjectElement {

public:
   TBrElement(std::unique_ptr<RHolder> &br) : TObjectElement(br) {}

   virtual ~TBrElement() = default;

   int GetNumChilds() override
   {
      auto br = fObject->Get<TBranchElement>();
      return br && br->IsFolder() ? TObjectElement::GetNumChilds() : 0;
   }

   /** Create iterator for childs elements if any */
   std::unique_ptr<RLevelIter> GetChildsIter() override
   {
      auto br = fObject->Get<TBranchElement>();
      if (br && br->IsFolder())
         return TObjectElement::GetChildsIter();
      return nullptr;
   }
};


////////////////////////////////////////////////////////////
/// Representing TVirtualBranchBrowsable in browsables

class TBrBrowsableElement : public TObjectElement {

public:
   TBrBrowsableElement(std::unique_ptr<RHolder> &br) : TObjectElement(br) {}

   virtual ~TBrBrowsableElement() = default;

   int GetNumChilds() override
   {
      auto br = fObject->Get<TVirtualBranchBrowsable>();
      return br && br->GetLeaves() ? br->GetLeaves()->GetSize() : 0;
   }

   /** Create iterator for childs elements if any */
   std::unique_ptr<RLevelIter> GetChildsIter() override
   {
      auto br = fObject->Get<TVirtualBranchBrowsable>();
      if (br && br->GetLeaves())
         return GetCollectionIter(br->GetLeaves());
      return nullptr;
   }
};


// ==============================================================================================

class TBranchBrowseProvider : public RProvider {

public:
   TBranchBrowseProvider()
   {
      RegisterBrowse(TBranchElement::Class(), [](std::unique_ptr<RHolder> &object) -> std::shared_ptr<RElement> {
         return std::make_shared<TBrElement>(object);
      });
   }

} newTBranchBrowseProvider;
