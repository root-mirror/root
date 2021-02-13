/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include <ROOT/Browsable/RAnyObjectHolder.hxx>
#include <ROOT/Browsable/RLevelIter.hxx>
#include <ROOT/Browsable/RProvider.hxx>
#include <ROOT/Browsable/TObjectHolder.hxx>
#include <ROOT/Browsable/TKeyItem.hxx>

#include <ROOT/RLogger.hxx>

#include "TKey.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TFile.h"

using namespace std::string_literals;

using namespace ROOT::Experimental::Browsable;


/** \class TDirectoryLevelIter
\ingroup rbrowser

Iterator over keys in TDirectory
*/


class TDirectoryLevelIter : public RLevelIter {
   TDirectory *fDir{nullptr};         ///<! current directory handle
   std::unique_ptr<TIterator>  fIter; ///<! created iterator
   TKey *fKey{nullptr};               ///<! currently selected key
   std::string fCurrentName;          ///<! current key name

   bool CreateIter()
   {
      if (!fDir) return false;
      fIter.reset(fDir->GetListOfKeys()->MakeIterator());
      fKey = nullptr;
      return true;
   }

   bool NextDirEntry()
   {
      fCurrentName.clear();
      if (!fIter) return false;

      fKey = dynamic_cast<TKey *>(fIter->Next());

      if (!fKey) {
         fIter.reset();
         return false;
      }

      fCurrentName = fKey->GetName();
      fCurrentName.append(";");
      fCurrentName.append(std::to_string(fKey->GetCycle()));

      return true;
   }

public:
   explicit TDirectoryLevelIter(TDirectory *dir) : fDir(dir) { CreateIter(); }

   virtual ~TDirectoryLevelIter() = default;

   bool Next() override { return NextDirEntry(); }

   // use default implementation for now
   // bool Find(const std::string &name) override { return FindDirEntry(name); }

   std::string GetItemName() const override { return fCurrentName; }

   bool CanItemHaveChilds() const override
   {
      return RProvider::CanHaveChilds(fKey->GetClassName());
   }

   /** Create element for the browser */
   std::unique_ptr<RItem> CreateItem() override
   {
      auto item = std::make_unique<TKeyItem>(GetItemName(), CanItemHaveChilds() ? -1 : 0);
      item->SetClassName(fKey->GetClassName());
      item->SetIcon(RProvider::GetClassIcon(fKey->GetClassName()));
      item->SetTitle(fKey->GetTitle());
      return item;
   }

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override;

};

// ===============================================================================================================


/** \class TKeyElement
\ingroup rbrowser

Element representing TKey from TDirectory
*/


class TKeyElement : public RElement {
   TDirectory *fDir{nullptr};
   TKey *fKey{nullptr};
   std::shared_ptr<RElement> fElement; ///<! holder of read object

public:
   TKeyElement(TDirectory *dir, TKey *key) : fDir(dir), fKey(key) {}

   virtual ~TKeyElement() = default;

   /** Name of TKeyElement, includes key cycle */
   std::string GetName() const override
   {
      if (fElement)
         return fElement->GetName();
      std::string name = fKey->GetName();
      name.append(";");
      name.append(std::to_string(fKey->GetCycle()));
      return name;
   }

   /** Title of TKeyElement (optional) */
   std::string GetTitle() const override
   {
      if (fElement)
         return fElement->GetTitle();

      return fKey->GetTitle();
   }

   /** Create iterator for childs elements if any
     * Means we should try to browse inside.
     * Either it is directory or some complex object */
   std::unique_ptr<RLevelIter> GetChildsIter() override
   {
      if (fElement)
         return fElement->GetChildsIter();

      std::string clname = fKey->GetClassName();

      if (clname.find("TDirectory") == 0) {
          auto subdir = fDir->GetDirectory(GetName().c_str());
          if (!subdir) return nullptr;
          return std::make_unique<TDirectoryLevelIter>(subdir);
      }

      auto obj = GetObject();

      if (obj)
         fElement = RProvider::Browse(obj);

      if (fElement)
         return fElement->GetChildsIter();

      return nullptr;
   }

   /** Return object associated with TKey, if TDirectory has object of that name it will be returned */
   std::unique_ptr<RHolder> GetObject() override
   {
      if (fElement)
         return fElement->GetObject();

      std::string clname = fKey->GetClassName();

      auto obj_class = TClass::GetClass(clname.c_str());
      if (!obj_class)
         return nullptr;

      void *obj = fDir->GetObjectChecked(fKey->GetName(), obj_class);
      if (!obj)
         return nullptr;

      TObject *tobj = (TObject *) obj_class->DynamicCast(TObject::Class(), obj);

      if (tobj) {
         bool owned_by_dir = (fDir->FindObject(tobj) == tobj) || (clname == "TGeoManager");

         return std::make_unique<TObjectHolder>(tobj, !owned_by_dir);
      }

      return std::make_unique<RAnyObjectHolder>(obj_class, obj, true);
   }


   EActionKind GetDefaultAction() const override
   {
      if (fElement)
         return fElement->GetDefaultAction();

      std::string clname = fKey->GetClassName();
      if (clname.empty()) return kActNone;
      if (clname == "TGeoManager") return kActGeom;
      if (RProvider::CanDraw6(clname)) return kActDraw6;
      if (RProvider::CanDraw7(clname)) return kActDraw7;
      if (RProvider::CanHaveChilds(clname)) return kActBrowse;
      return kActNone;

   }

   bool IsCapable(EActionKind action) const override
   {
      if (fElement)
         return fElement->IsCapable(action);

      std::string clname = fKey->GetClassName();
      if (clname.empty()) return false;

      switch(action) {
         case kActBrowse: return RProvider::CanHaveChilds(clname);
         case kActEdit: return true;
         case kActImage:
         case kActDraw6: return RProvider::CanDraw6(clname); // if can draw in TCanvas, can produce image
         case kActDraw7: return RProvider::CanDraw7(clname);
         case kActGeom: return (clname == "TGeoManager");
         default: return false;
      }

      return false;
   }

};

// ==============================================================================================


/////////////////////////////////////////////////////////////////////////////////
/// Return element for current TKey object in TDirectory

std::shared_ptr<RElement> TDirectoryLevelIter::GetElement()
{
   return std::make_shared<TKeyElement>(fDir, fKey);
}

// ==============================================================================================

/** \class TDirectoryElement
\ingroup rbrowser

Element representing TDirectory
*/


class TDirectoryElement : public RElement {
   std::string fFileName;       ///<!   file name
   TDirectory *fDir{nullptr};   ///<!   subdirectory (ifany)

   ///////////////////////////////////////////////////////////////////
   /// Get TDirectory. Checks if parent file is still there. If not, means it was closed outside ROOT

   TDirectory *GetDir()
   {
      if (fDir) {
         if (!gROOT->GetListOfFiles()->FindObject(fDir->GetFile()))
            fDir = nullptr;
      } else if (!fFileName.empty()) {
         fDir = TFile::Open(fFileName.c_str());
      }

      return fDir;
   }

public:

   TDirectoryElement(const std::string &fname, TDirectory *dir = nullptr)
   {
      fFileName = fname;
      fDir = dir;
   }

   virtual ~TDirectoryElement() = default;

   /** Name of TDirectoryElement */
   std::string GetName() const override
   {
      if (fDir)
         return fDir->GetName();

      if (!fFileName.empty()) {
         auto pos = fFileName.rfind("/");
         return ((pos == std::string::npos) || (pos > fFileName.length() - 2)) ? fFileName : fFileName.substr(pos + 1);
      }

      return ""s;
   }

   /** Title of TDirectoryElement */
   std::string GetTitle() const override
   {
      if (fDir)
         return fDir->GetTitle();

      return "ROOT file "s + fFileName;
   }

   std::unique_ptr<RLevelIter> GetChildsIter() override
   {
      auto dir = GetDir();

      return dir ? std::make_unique<TDirectoryLevelIter>(dir) : nullptr;
   }
};


// ==============================================================================================


/** \class RTFileProvider
\ingroup rbrowser

Provides access to ROOT files with extension "root"
Other extensions can be registered
*/


class RTFileProvider : public RProvider {

public:
   RTFileProvider()
   {
      RegisterFile("root", [] (const std::string &fullname) -> std::shared_ptr<RElement> {
         auto f = dynamic_cast<TFile *> (gROOT->GetListOfFiles()->FindObject(fullname.c_str()));
         if (!f) f = TFile::Open(fullname.c_str());
         if (!f) return nullptr;
         return std::make_shared<TDirectoryElement>(fullname, f);
      });

      RegisterBrowse(TFile::Class(), [](std::unique_ptr<RHolder> &object) -> std::shared_ptr<RElement> {
         return std::make_shared<TDirectoryElement>("", const_cast<TFile*>(object->Get<TFile>()));
      });

      RegisterBrowse(TDirectory::Class(), [](std::unique_ptr<RHolder> &object) -> std::shared_ptr<RElement> {
         return std::make_shared<TDirectoryElement>("", const_cast<TDirectory*>(object->Get<TDirectory>()));
      });
   }

} newRTFileProvider ;


