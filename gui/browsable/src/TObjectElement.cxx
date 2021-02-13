/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
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

   bool CanItemHaveChilds() const override
   {
      std::shared_ptr<TObjectElement> telem = std::dynamic_pointer_cast<TObjectElement>(fElements[fCounter]);
      return telem ? telem->IsFolder() : false;
   }

   /** Create element for the browser */
   std::unique_ptr<RItem> CreateItem() override
   {
      std::shared_ptr<TObjectElement> elem = std::dynamic_pointer_cast<TObjectElement>(fElements[fCounter]);
      // should never happen
      if (!elem) return nullptr;

      auto cl = elem->GetClass();

      auto nchilds = elem->GetNumChilds();

      auto item = std::make_unique<TObjectItem>(elem->GetName(), nchilds);

      item->SetClassName(cl ? cl->GetName() : "");

      item->SetIcon(RProvider::GetClassIcon(cl, nchilds > 0));

      item->SetTitle(elem->GetTitle());

      return item;
   }

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override
   {
      return fElements[fCounter];
   }

   bool Find(const std::string &name, int indx = -1) override
   {
      if ((indx >= 0) && (indx < (int) fElements.size()) && (name == fElements[indx]->GetName())) {
         fCounter = indx;
         return true;
      }

      return RLevelIter::Find(name, -1);
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

   bool IsDuplicated() const { return fDuplicated; }

   void Add(TObject* obj, const char* name, Int_t) override
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

};

// ===============================================================================================================

/** \class TCollectionIter
\ingroup rbrowser

Iterator over elements in TCollection
*/

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

   bool CanItemHaveChilds() const override
   {
      TObject *obj = *fIter;
      return obj ? obj->IsFolder() : false;
   }

   /** Create item for current TObject */
   std::unique_ptr<RItem> CreateItem() override
   {
      TObject *obj = *fIter;
      if (!obj) return nullptr;

      auto item = std::make_unique<TObjectItem>(obj->GetName(), obj->IsFolder() ? -1 : 0);

      item->SetClassName(obj->ClassName());

      item->SetIcon(RProvider::GetClassIcon(obj->IsA(), obj->IsFolder()));

      item->SetTitle(obj->GetTitle());

      return item;
   }

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override
   {
      std::unique_ptr<RHolder> holder = std::make_unique<TObjectHolder>(*fIter, kFALSE);

      return RProvider::Browse(holder);
   }

};


// ==============================================================================================

/** \class TFolderElement
\ingroup rbrowser

Browsable element for TFolder
*/


class TFolderElement : public TObjectElement {

public:

   TFolderElement(std::unique_ptr<RHolder> &obj) : TObjectElement(obj) {}

   std::unique_ptr<RLevelIter> GetChildsIter() override
   {
      auto folder = fObject->Get<TFolder>();
      if (folder)
         return std::make_unique<TCollectionIter>(folder->GetListOfFolders());

      return TObjectElement::GetChildsIter();
   }

   int GetNumChilds() override
   {
      auto folder = fObject->Get<TFolder>();
      return folder && folder->GetListOfFolders() ? folder->GetListOfFolders()->GetEntries() : 0;
   }
};

// ==============================================================================================

/** \class TCollectionElement
\ingroup rbrowser

Browsable element for TCollection
*/


class TCollectionElement : public TObjectElement {
public:

   TCollectionElement(std::unique_ptr<RHolder> &obj) : TObjectElement(obj) {}

   std::unique_ptr<RLevelIter> GetChildsIter() override
   {
      auto coll = fObject->Get<TCollection>();
      if (coll && (coll->GetSize() > 0))
         return std::make_unique<TCollectionIter>(coll);

      return TObjectElement::GetChildsIter();
   }

   int GetNumChilds() override
   {
      auto coll = fObject->Get<TCollection>();
      return coll ? coll->GetSize() : 0;
   }
};

// =================================================================================

////////////////////////////////////////////////////////////////////////////////
/// Constructor with plain TObject* as argument - ownership is not defined

TObjectElement::TObjectElement(TObject *obj, const std::string &name) : fObj(obj), fName(name)
{
   fObject = std::make_unique<TObjectHolder>(fObj);
   if (fName.empty())
      fName = fObj->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with std::unique_ptr<RHolder> as argument

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

////////////////////////////////////////////////////////////////////////////////
/// Returns name of the TObject

std::string TObjectElement::GetName() const
{
   if (!fName.empty()) return fName;
   return fObj ? fObj->GetName() : "";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns title of the TObject

std::string TObjectElement::GetTitle() const
{
   return fObj ? fObj->GetTitle() : "";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns IsFolder of contained TObject

bool TObjectElement::IsFolder()
{
   return fObj ? fObj->IsFolder() : false;
}

////////////////////////////////////////////////////////////////////////////////
/// Create iterator for childs elements if any

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

////////////////////////////////////////////////////////////////////////////////
/// Returns copy of TObject holder - if possible

std::unique_ptr<RHolder> TObjectElement::GetObject()
{
   if (!fObject)
      return nullptr;

   return fObject->Copy();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns class for contained object

const TClass *TObjectElement::GetClass() const
{
   return fObj ? fObj->IsA() : nullptr;
}


////////////////////////////////////////////////////////////////////////////////
/// Provides default action which can be performed with the object

RElement::EActionKind TObjectElement::GetDefaultAction() const
{
   auto cl = GetClass();
   if (!cl) return kActNone;
   if (("TGeoManager"s == cl->GetName()) || ("TGeoVolume"s == cl->GetName())) return kActGeom;
   if (RProvider::CanDraw6(cl)) return kActDraw6;
   if (RProvider::CanDraw7(cl)) return kActDraw7;
   if (RProvider::CanHaveChilds(cl)) return kActBrowse;
   return kActNone;
}

////////////////////////////////////////////////////////////////////////////////
/// Check object capability

bool TObjectElement::IsCapable(RElement::EActionKind action) const
{
   auto cl = GetClass();
   if (!cl) return false;

   switch(action) {
      case kActBrowse: return RProvider::CanHaveChilds(cl);
      case kActEdit: return true;
      case kActImage:
      case kActDraw6: return RProvider::CanDraw6(cl); // if can draw in TCanvas, can produce image
      case kActDraw7: return RProvider::CanDraw7(cl);
      case kActGeom: return ("TGeoManager"s == cl->GetName()) || ("TGeoVolume"s == cl->GetName());
      default: return false;
   }

   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates iterator for TCollection object

std::unique_ptr<RLevelIter> TObjectElement::GetCollectionIter(const TCollection *coll)
{
   return std::make_unique<TCollectionIter>(coll);
}

// ==============================================================================================

/** \class RTObjectProvider
\ingroup rbrowser

Provider for all known TObject-based classes
*/

class RTObjectProvider : public RProvider {

public:
   //////////////////////////////////////////////////////////////////////////////////
   // Register TObject-based class with standard browsing/drawing libs

   void RegisterTObject(const std::string &clname, const std::string &iconname, bool can_browse = false, int can_draw = 3)
   {
      RegisterClass(clname, iconname, can_browse ? "dflt"s : ""s,
                                      can_draw & 1 ? "libROOTObjectDraw6Provider"s : ""s,
                                      can_draw & 2 ? "libROOTObjectDraw7Provider"s : ""s);
   }

   RTObjectProvider()
   {
      RegisterTObject("TTree", "sap-icon://tree", true, 0);
      RegisterTObject("TNtuple", "sap-icon://tree", true, 0);
      RegisterClass("TBranchElement", "sap-icon://e-care", "libROOTBranchBrowseProvider", "libROOTLeafDraw6Provider", "libROOTLeafDraw7Provider");
      RegisterClass("TLeaf", "sap-icon://e-care", ""s, "libROOTLeafDraw6Provider", "libROOTLeafDraw7Provider");
      RegisterClass("TVirtualBranchBrowsable", "sap-icon://e-care", ""s, "libROOTLeafDraw6Provider", "libROOTLeafDraw7Provider");

      RegisterTObject("TDirectory", "sap-icon://folder-blank", true, 0);
      RegisterTObject("TH1", "sap-icon://bar-chart");
      RegisterTObject("TH2", "sap-icon://pixelate");
      RegisterTObject("TH3", "sap-icon://product");
      RegisterTObject("TProfile", "sap-icon://vertical-bar-chart");
      RegisterTObject("TGraph", "sap-icon://line-chart");
      RegisterTObject("TCanvas", "sap-icon://business-objects-experience", false, 1); // only can use TWebCanvas

      RegisterTObject("THStack", "sap-icon://multiple-bar-chart");
      RegisterTObject("TMultiGraph", "sap-icon://multiple-line-chart");

      RegisterTObject("TCollection", "sap-icon://list", true, 0);
      RegisterTObject("TGeoManager", "sap-icon://overview-chart", true, 0);
      RegisterTObject("TGeoVolume", "sap-icon://product", true, 0);
      RegisterTObject("TGeoNode", "sap-icon://product", true, 0);

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
