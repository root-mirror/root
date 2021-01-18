/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/TObjectElement.hxx>
#include <ROOT/Browsable/RProvider.hxx>
#include <ROOT/Browsable/TObjectHolder.hxx>
#include <ROOT/Browsable/TObjectItem.hxx>
#include <ROOT/Browsable/RLevelIter.hxx>

#include <ROOT/RLogger.hxx>

#include "TBrowser.h"
#include "TBrowserImp.h"
#include "TFolder.h"
#include "TList.h"
#include "TDirectory.h"
#include "TBufferJSON.h"

#include <sstream>

using namespace std::string_literals;

using namespace ROOT::Experimental::Browsable;

/** \class TObjectLevelIter
\ingroup rbrowser

Iterator over list of elements, designed for support TBrowser usage
*/

class TObjectLevelIter : public RLevelIter {

   std::vector<std::shared_ptr<RElement>> fElements;

   int fCounter{-1};

public:
   explicit TObjectLevelIter() {}

   virtual ~TObjectLevelIter() = default;

   void AddElement(std::shared_ptr<RElement> &&elem)
   {
      fElements.emplace_back(std::move(elem));
   }

   auto NumElements() const { return fElements.size(); }

   bool Next() override { return ++fCounter < (int) fElements.size(); }

   // use default implementation for now
   // bool Find(const std::string &name) override { return FindDirEntry(name); }

   std::string GetItemName() const override { return fElements[fCounter]->GetName(); }

   int GetNumItemChilds() const override
   {
      std::shared_ptr<TObjectElement> telem = std::dynamic_pointer_cast<TObjectElement>(fElements[fCounter]);
      if (!telem) return 0;
      return telem->IsFolder() ? -1 : 0;
   }

   /** Create element for the browser */
   std::unique_ptr<RItem> CreateItem() override;

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override
   {
      return fElements[fCounter];
   }

};

// ===============================================================================================================

class TMyBrowserImp : public TBrowserImp {
   TObjectLevelIter *fIter{nullptr};   ///<!  back-reference on iterator
   const TObject *fBrowseObj{nullptr}; ///<!  object which wil be browsed
   bool fDuplicated{false};            ///<! is object was duplicated?

public:

   TMyBrowserImp(TObjectLevelIter *iter, TObject *obj) : TBrowserImp(nullptr), fIter(iter), fBrowseObj(obj) {}
   virtual ~TMyBrowserImp() = default;

   void Add(TObject* obj, const char* name, Int_t) override;

   bool IsDuplicated() const { return fDuplicated; }
};

// ===============================================================================================================


TObjectElement::TObjectElement(TObject *obj, const std::string &name) : fObj(obj), fName(name)
{
   fObject = std::make_unique<TObjectHolder>(fObj);
   if (fName.empty())
      fName = fObj->GetName();
}

TObjectElement::TObjectElement(std::unique_ptr<RHolder> &obj, const std::string &name)
{
   fObject = std::move(obj); // take responsibility
   fObj = const_cast<TObject *>(fObject->Get<TObject>()); // try to cast into TObject

   fName = name;
   if (!fObj)
      fObject.reset();
   else if (fName.empty())
      fName = fObj->GetName();
}

TObjectElement::~TObjectElement()
{
}

std::string TObjectElement::GetName() const
{
   if (!fName.empty()) return fName;
   return fObj ? fObj->GetName() : "";
}

/** Title of TObject */
std::string TObjectElement::GetTitle() const
{
   return fObj ? fObj->GetTitle() : "";
}

/** Returns IsFolder of contained object */
bool TObjectElement::IsFolder()
{
   return fObj ? fObj->IsFolder() : false;
}

/** Create iterator for childs elements if any */
std::unique_ptr<RLevelIter> TObjectElement::GetChildsIter()
{
   if (!IsFolder()) return nullptr;

   auto iter = std::make_unique<TObjectLevelIter>();

   TMyBrowserImp *imp = new TMyBrowserImp(iter.get(), fObj);

   // must be new, otherwise TBrowser constructor ignores imp
   TBrowser *br = new TBrowser("name", "title", imp);

   fObj->Browse(br);

   auto dupl = imp->IsDuplicated();

   delete br; // also will destroy implementaion

   if (dupl || (iter->NumElements() == 0)) return nullptr;

   return iter;
}

/** Return copy of TObject holder - if possible */
std::unique_ptr<RHolder> TObjectElement::GetObject()
{
   if (!fObject)
      return nullptr;

   return fObject->Copy();
}

std::string TObjectElement::ClassName() const
{
   return fObj ? fObj->ClassName() : "";
}

// ==============================================================================================


void TMyBrowserImp::Add(TObject *obj, const char *name, Int_t)
{
   // prevent duplication of object itself - ignore such browsing
   if (fBrowseObj == obj) fDuplicated = true;
   if (fDuplicated) return;

   std::unique_ptr<RHolder> holder = std::make_unique<TObjectHolder>(obj);

   std::shared_ptr<RElement> elem = RProvider::Browse(holder);

   if (name && *name) {
      std::shared_ptr<TObjectElement> telem = std::dynamic_pointer_cast<TObjectElement>(elem);
      if (telem) telem->SetName(name);
   }

   fIter->AddElement(std::move(elem));
}


// ==============================================================================================


///////////////////////////////////////////////////////////////
/// Create element for the browser

std::unique_ptr<RItem> TObjectLevelIter::CreateItem()
{
   std::shared_ptr<TObjectElement> elem = std::dynamic_pointer_cast<TObjectElement>(fElements[fCounter]);
   // should never happen
   if (!elem) return nullptr;

   bool can_have_childs = false;
   if (GetNumItemChilds() != 0) {
      // TODO: make via RProvider methods
      std::string clname = elem->ClassName();
      can_have_childs = (clname.find("TDirectory") == 0) || (clname.find("TTree") == 0) ||
                        (clname.find("TNtuple") == 0) || (clname.find("TBranchElement") == 0) ||
                        (clname.find("TGeoManager") == 0) || (clname.find("TGeoVolume") == 0) || (clname.find("TGeoNode") == 0);
   }

   auto item = std::make_unique<TObjectItem>(elem->GetName(), can_have_childs ? -1 : 0);

   item->SetClassName(elem->ClassName());

   item->SetIcon(RProvider::GetClassIcon(elem->ClassName()));

   return item;
}


// ==============================================================================================

class TFolderElement : public TObjectElement {

public:

   TFolderElement(std::unique_ptr<RHolder> &obj) : TObjectElement(obj) {}

   std::unique_ptr<RLevelIter> GetChildsIter() override;
};

class TCollectionElement : public TObjectElement {
public:

   TCollectionElement(std::unique_ptr<RHolder> &obj) : TObjectElement(obj) {}

   std::unique_ptr<RLevelIter> GetChildsIter() override;
};



class TCollectionIter : public RLevelIter {

   TIter  fIter;

public:
   explicit TCollectionIter(const TFolder *f) : RLevelIter(), fIter(f->GetListOfFolders()) {};

   explicit TCollectionIter(const TCollection *coll) : RLevelIter(), fIter(coll) {};

   virtual ~TCollectionIter() = default;

   bool Next() override { return fIter.Next() != nullptr; }

   // use default implementation for now
   // bool Find(const std::string &name) override { return FindDirEntry(name); }

   std::string GetItemName() const override { return (*fIter)->GetName(); }

   int GetNumItemChilds() const override
   {
      // TODO: add test based on class name
      return -1;
   }

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override
   {
      std::unique_ptr<RHolder> holder = std::make_unique<TObjectHolder>(*fIter, kFALSE);

      return RProvider::Browse(holder);
   }
};

//////////////////////////////////////////////////////////////////////////////////////
/// Provides iterator for TFolder

std::unique_ptr<RLevelIter> TFolderElement::GetChildsIter()
{
  auto folder = fObject->Get<TFolder>();
  if (folder)
     return std::make_unique<TCollectionIter>(folder->GetListOfFolders());

  return TObjectElement::GetChildsIter();
}

//////////////////////////////////////////////////////////////////////////////////////
/// Provides iterator for generic TCollecion

std::unique_ptr<RLevelIter> TCollectionElement::GetChildsIter()
{
   auto coll = fObject->Get<TCollection>();
   if (coll)
      return std::make_unique<TCollectionIter>(coll);

   return TObjectElement::GetChildsIter();
}

// ==============================================================================================

class RTObjectProvider : public RProvider {

public:
   RTObjectProvider()
   {
      RegisterClass("TTree", "sap-icon://tree");
      RegisterClass("TNtuple", "sap-icon://tree");

      RegisterBrowse(TFolder::Class(), [](std::unique_ptr<RHolder> &object) -> std::shared_ptr<RElement> {
         return std::make_shared<TFolderElement>(object);
      });

      RegisterBrowse(TCollection::Class(), [](std::unique_ptr<RHolder> &object) -> std::shared_ptr<RElement> {
         return std::make_shared<TCollectionElement>(object);
      });

      RegisterBrowse(nullptr, [](std::unique_ptr<RHolder> &object) -> std::shared_ptr<RElement> {
         if (object->CanCastTo<TObject>())
            return std::make_shared<TObjectElement>(object);
         return nullptr;
      });
   }

} newRTObjectProvider;
