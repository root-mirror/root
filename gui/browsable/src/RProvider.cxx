/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RProvider.hxx>

#include <ROOT/RLogger.hxx>

#include "TBaseClass.h"
#include "TList.h"
#include "TSystem.h"

using namespace ROOT::Experimental::Browsable;
using namespace std::string_literals;

//////////////////////////////////////////////////////////////////////////////////
// Provide map of browsing for different classes

RProvider::BrowseMap_t &RProvider::GetBrowseMap()
{
   static RProvider::BrowseMap_t sMap;
   return sMap;
}

//////////////////////////////////////////////////////////////////////////////////
// Provide map of files opening

RProvider::FileMap_t &RProvider::GetFileMap()
{
   static RProvider::FileMap_t sMap;
   return sMap;
}

//////////////////////////////////////////////////////////////////////////////////
// Returns map of registered drawing functions for v6 canvas

RProvider::Draw6Map_t &RProvider::GetDraw6Map()
{
   static RProvider::Draw6Map_t sMap;
   return sMap;
}

//////////////////////////////////////////////////////////////////////////////////
// Returns map of registered drawing functions for v7 canvas

RProvider::Draw7Map_t &RProvider::GetDraw7Map()
{
   static RProvider::Draw7Map_t sMap;
   return sMap;
}

//////////////////////////////////////////////////////////////////////////////////
// Destructor
/// Automatically unregister provider from global lists

RProvider::~RProvider()
{
   // here to remove all correspondent entries
   auto &fmap = GetFileMap();

   for (auto fiter = fmap.begin();fiter != fmap.end();) {
      if (fiter->second.provider == this)
         fiter = fmap.erase(fiter);
      else
         fiter++;
   }

   auto &bmap = GetBrowseMap();
   for (auto biter = bmap.begin(); biter != bmap.end();) {
      if (biter->second.provider == this)
         biter = bmap.erase(biter);
      else
         biter++;
   }

   auto &map6 = GetDraw6Map();
   for (auto iter6 = map6.begin(); iter6 != map6.end();) {
      if (iter6->second.provider == this)
         iter6 = map6.erase(iter6);
      else
         iter6++;
   }

   auto &map7 = GetDraw7Map();
   for (auto iter7 = map7.begin(); iter7 != map7.end();) {
      if (iter7->second.provider == this)
         iter7 = map7.erase(iter7);
      else
         iter7++;
   }

}

//////////////////////////////////////////////////////////////////////////////////
// Register file open function for specified extension

void RProvider::RegisterFile(const std::string &extension, FileFunc_t func)
{
    auto &fmap = GetFileMap();

    if ((extension != "*") && (fmap.find(extension) != fmap.end()))
       R__ERROR_HERE("Browserv7") << "Provider for file extension  " << extension << " already exists";

    fmap.emplace(extension, StructFile{this,func});
}

//////////////////////////////////////////////////////////////////////////////////
// Register browse function for specified class

void RProvider::RegisterBrowse(const TClass *cl, BrowseFunc_t func)
{
    auto &bmap = GetBrowseMap();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__ERROR_HERE("Browserv7") << "Browse provider for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructBrowse{this,func});
}


//////////////////////////////////////////////////////////////////////////////////
// Register drawing function for v6 canvas

void RProvider::RegisterDraw6(const TClass *cl, Draw6Func_t func)
{
    auto &bmap = GetDraw6Map();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__ERROR_HERE("Browserv7") << "Draw v6 handler for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructDraw6{this, func});
}

//////////////////////////////////////////////////////////////////////////////////
// Register drawing function for v6 canvas

void RProvider::RegisterDraw7(const TClass *cl, Draw7Func_t func)
{
    auto &bmap = GetDraw7Map();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__ERROR_HERE("Browserv7") << "Draw v7 handler for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructDraw7{this, func});
}


//////////////////////////////////////////////////////////////////////////////////
// remove provider from all registered lists

std::shared_ptr<RElement> RProvider::OpenFile(const std::string &extension, const std::string &fullname)
{
   auto &fmap = GetFileMap();

   auto iter = fmap.find(extension);

   if (iter != fmap.end()) {
      auto res = iter->second.func(fullname);
      if (res) return res;
   }

   for (auto &pair : fmap)
      if ((pair.first == "*") || (pair.first == extension)) {
         auto res = pair.second.func(fullname);
         if (res) return res;
      }

   return nullptr;
}

/////////////////////////////////////////////////////////////////////////
/// Create browsable element for the object
/// Created element may take ownership over the object

std::shared_ptr<RElement> RProvider::Browse(std::unique_ptr<RHolder> &object)
{
   if (!object)
      return nullptr;

   auto &bmap = GetBrowseMap();

   auto cl = object->GetClass();
   if (!cl)
      return nullptr;

   auto iter = bmap.find(cl);

   if (iter != bmap.end()) {
      auto res = iter->second.func(object);
      if (res || !object) return res;
   }

   for (auto &pair : bmap)
      if ((pair.first == nullptr) || (cl == pair.first)) {
         auto res = pair.second.func(object);
         if (res || !object) return res;
      }

   return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////
/// Invoke drawing of object on TCanvas sub-pad
/// All existing providers are checked, first checked are class matches (including direct parents)

bool RProvider::Draw6(TVirtualPad *subpad, std::unique_ptr<Browsable::RHolder> &obj, const std::string &opt)
{
   if (!obj || !obj->GetClass())
      return false;

   auto &map6 = GetDraw6Map();

   TClass *cl = const_cast<TClass *>(obj->GetClass());
   while (cl) {
      auto iter6 = map6.find(cl);

      if (iter6 != map6.end()) {
         if (iter6->second.func(subpad, obj, opt))
            return true;
      }

      auto bases = cl->GetListOfBases();

      cl = bases && (bases->GetSize() > 0) ? dynamic_cast<TBaseClass *>(bases->First())->GetClassPointer() : nullptr;
   }

   for (auto &pair : map6)
      if ((pair.first == obj->GetClass()) || !pair.first)
         if (pair.second.func(subpad, obj, opt))
            return true;

   // try to load necessary library and repeat action again
   // TODO: need factory methods for that

   if (obj->GetClass()->InheritsFrom("TLeaf"))
      gSystem->Load("libROOTTreeDrawProvider");
   else if (obj->GetClass()->InheritsFrom(TObject::Class()))
      gSystem->Load("libROOTObjectDrawProvider");
   else
      return false;

   cl = const_cast<TClass *>(obj->GetClass());
   while (cl) {
      auto iter6 = map6.find(cl);

      if (iter6 != map6.end()) {
         if (iter6->second.func(subpad, obj, opt))
            return true;
      }

      auto bases = cl->GetListOfBases();

      cl = bases && (bases->GetSize() > 0) ? dynamic_cast<TBaseClass *>(bases->First())->GetClassPointer() : nullptr;
   }

   for (auto &pair : map6)
      if ((pair.first == obj->GetClass()) || !pair.first)
         if (pair.second.func(subpad, obj, opt))
            return true;

   return false;
}

/////////////////////////////////////////////////////////////////////////////////
/// Invoke drawing of object on RCanvas sub-pad
/// All existing providers are checked, first checked are class matches (including direct parents)

bool RProvider::Draw7(std::shared_ptr<RPadBase> &subpad, std::unique_ptr<Browsable::RHolder> &obj, const std::string &opt)
{
   if (!obj || !obj->GetClass())
      return false;

   auto &map7 = GetDraw7Map();

   TClass *cl = const_cast<TClass *>(obj->GetClass());
   while (cl) {
      auto iter7 = map7.find(cl);

      if (iter7 != map7.end()) {
         if (iter7->second.func(subpad, obj, opt))
            return true;
      }

      auto bases = cl->GetListOfBases();

      cl = bases && (bases->GetSize() > 0) ? dynamic_cast<TBaseClass *>(bases->First())->GetClassPointer() : nullptr;
   }

   for (auto &pair : map7)
      if ((pair.first == obj->GetClass()) || !pair.first)
         if (pair.second.func(subpad, obj, opt))
            return true;

   // try to load necessary library and repeat action again
   // TODO: need factory methods for that

   if (obj->GetClass()->InheritsFrom("TLeaf"))
      gSystem->Load("libROOTTreeDrawProvider");
   else if (obj->GetClass()->InheritsFrom(TObject::Class()))
      gSystem->Load("libROOTObjectDrawProvider");
   else if (obj->GetClass()->InheritsFrom("ROOT::Experimental::RH1D") || obj->GetClass()->InheritsFrom("ROOT::Experimental::RH2D") || obj->GetClass()->InheritsFrom("ROOT::Experimental::RH2D"))
      gSystem->Load("libROOTHistDrawProvider");
   else
      return false;

   cl = const_cast<TClass *>(obj->GetClass());
   while (cl) {
      auto iter7 = map7.find(cl);

      if (iter7 != map7.end()) {
         if (iter7->second.func(subpad, obj, opt))
            return true;
      }

      auto bases = cl->GetListOfBases();

      cl = bases && (bases->GetSize() > 0) ? dynamic_cast<TBaseClass *>(bases->First())->GetClassPointer() : nullptr;
   }

   for (auto &pair : map7)
      if ((pair.first == obj->GetClass()) || !pair.first)
         if (pair.second.func(subpad, obj, opt))
            return true;

   return false;
}


/////////////////////////////////////////////////////////////////////
/// Return icon name for the given class
/// TODO: should be factorized out from here

std::string RProvider::GetClassIcon(const std::string &classname)
{
   if (classname == "TTree" || classname == "TNtuple")
      return "sap-icon://tree"s;
   else if (classname == "TDirectory" || classname == "TDirectoryFile")
      return "sap-icon://folder-blank"s;
   else if (classname.find("TLeaf") == 0)
      return "sap-icon://e-care"s;

   return "sap-icon://electronic-medical-record"s;
}
