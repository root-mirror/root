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
#include "TClass.h"
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
// Returns map of registered icons base on class pointer

RProvider::ClassMap_t &RProvider::GetClassMap()
{
   static RProvider::ClassMap_t sMap;
   return sMap;
}

//////////////////////////////////////////////////////////////////////////////////
// Returns map of registered icons base on class name

//////////////////////////////////////////////////////////////////////////////////
// Destructor
/// Automatically unregister provider from all maps

RProvider::~RProvider()
{
   // here to remove all correspondent entries
   CleanThis(GetClassMap());
   CleanThis(GetFileMap());
   CleanThis(GetBrowseMap());
   CleanThis(GetDraw6Map());
   CleanThis(GetDraw7Map());
}

//////////////////////////////////////////////////////////////////////////////////
// Register file open function for specified extension

void RProvider::RegisterFile(const std::string &extension, FileFunc_t func)
{
    auto &fmap = GetFileMap();

    if ((extension != "*") && (fmap.find(extension) != fmap.end()))
       R__LOG_ERROR(BrowsableLog()) << "Provider for file extension  " << extension << " already exists";

    fmap.emplace(extension, StructFile{this,func});
}

//////////////////////////////////////////////////////////////////////////////////
// Register browse function for specified class

void RProvider::RegisterBrowse(const TClass *cl, BrowseFunc_t func)
{
    auto &bmap = GetBrowseMap();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__LOG_ERROR(BrowsableLog()) << "Browse provider for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructBrowse{this,func});
}


//////////////////////////////////////////////////////////////////////////////////
// Register drawing function for v6 canvas

void RProvider::RegisterDraw6(const TClass *cl, Draw6Func_t func)
{
    auto &bmap = GetDraw6Map();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__LOG_ERROR(BrowsableLog()) << "Draw v6 handler for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructDraw6{this, func});
}

//////////////////////////////////////////////////////////////////////////////////
// Register drawing function for v7 canvas

void RProvider::RegisterDraw7(const TClass *cl, Draw7Func_t func)
{
    auto &bmap = GetDraw7Map();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__LOG_ERROR(BrowsableLog()) << "Draw v7 handler for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructDraw7{this, func});
}

//////////////////////////////////////////////////////////////////////////////////
// Register class with supported libs (if any)

void RProvider::RegisterClass(const std::string &clname, const std::string &iconname,
                              const std::string &browselib, const std::string &draw6lib, const std::string &draw7lib)
{
   auto &bmap = GetClassMap();

   if (!clname.empty() && (bmap.find(clname) != bmap.end()))
      R__LOG_ERROR(BrowsableLog()) << "Entry for class " << clname << " already exists";

   std::string blib = browselib;
   bool can_have_childs = !browselib.empty();
   if ((blib == "dflt") || (blib == "TObject")) blib = ""; // just use as indicator that browsing is possible

   bmap.emplace(clname, StructClass{this, can_have_childs, iconname, blib, draw6lib, draw7lib});
}

//////////////////////////////////////////////////////////////////////////////////
// Returns entry for the requested class
const RProvider::StructClass &RProvider::GetClassEntry(const ClassArg &cl)
{
   if (!cl.empty()) {
      auto &bmap = GetClassMap();
      auto iter = bmap.find(cl.cl ? cl.cl->GetName() : cl.name.c_str());
      if (iter != bmap.end())
         return iter->second;

      if (!cl.name.empty()) {
         for (auto &elem : bmap)
            if (cl.name.compare(0, elem.first.length(), elem.first) == 0)
               return elem.second;
      } else {
         auto bases = const_cast<TClass *>(cl.cl)->GetListOfBases();
         const TClass *basecl = bases && (bases->GetSize() > 0) ? dynamic_cast<TBaseClass *>(bases->First())->GetClassPointer() : nullptr;
         if (basecl) return RProvider::GetClassEntry(basecl);
      }
   }

   static StructClass dummy;
   return dummy;
}

//////////////////////////////////////////////////////////////////////////////////
// Returns true if file extension is supported

bool RProvider::IsFileFormatSupported(const std::string &extension)
{
   if (extension.empty())
      return false;

   auto &fmap = GetFileMap();

   return fmap.find(extension) != fmap.end();
}

//////////////////////////////////////////////////////////////////////////////////
// Try to open file using provided extension.

std::shared_ptr<RElement> RProvider::OpenFile(const std::string &extension, const std::string &fullname)
{
   auto &fmap = GetFileMap();

   auto iter = fmap.find(extension);

   if (iter != fmap.end())
      return iter->second.func(fullname);

   return nullptr;
}

//////////////////////////////////////////////////////////////////////////////////
// Template function to scan class entries, including parent object classes

template<class Map_t, class Func_t>
bool ScanProviderMap(Map_t &fmap, const RProvider::ClassArg &cl, bool test_all = false, std::function<bool(Func_t &)> check_func = nullptr)
{
   if (cl.empty())
      return false;

   if (cl.GetClass()) {
      TClass *testcl = const_cast<TClass *>(cl.GetClass());
      while (testcl) {
         auto iter = fmap.find(testcl);
         if (iter != fmap.end())
            if (!check_func || check_func(iter->second.func))
               return true;

         auto bases = testcl->GetListOfBases();

         testcl = bases && (bases->GetSize() > 0) ? dynamic_cast<TBaseClass *>(bases->First())->GetClassPointer() : nullptr;
      }
   } else {
      for (auto &entry : fmap) {
         if (!entry.first) continue;
         std::string name = entry.first->GetName();
         if (!check_func) {
            // when check_func not specified, just try to guess if class can match
            if ((cl.GetName() == name) || (cl.GetName().compare(0, name.length(), name) == 0))
               return true;
         } else if (cl.GetName() == name) {
            if (check_func(entry.second.func))
               return true;
         }
      }
   }

   if (test_all && check_func) {
      for (auto &entry : fmap)
         if (!entry.first && check_func(entry.second.func))
            return true;
   }

   return false;
}


/////////////////////////////////////////////////////////////////////////
/// Create browsable element for the object
/// Created element may take ownership over the object

std::shared_ptr<RElement> RProvider::Browse(std::unique_ptr<RHolder> &object)
{
   std::shared_ptr<RElement> res;

   if (!object) return res;

   auto browse_func = [&res, &object] (BrowseFunc_t &func) -> bool {
      res = func(object);
      return (res || !object) ? true : false;
   };

   // check only class entries
   if (ScanProviderMap<BrowseMap_t,BrowseFunc_t>(GetBrowseMap(), object->GetClass(), false, browse_func))
      return res;

   auto &entry = GetClassEntry(object->GetClass());
   if (!entry.dummy() && !entry.browselib.empty())
      gSystem->Load(entry.browselib.c_str());

   // let call also generic browse functions (multicast)
   ScanProviderMap<BrowseMap_t,BrowseFunc_t>(GetBrowseMap(), object->GetClass(), true, browse_func);

   return res;
}

/////////////////////////////////////////////////////////////////////////////////
/// Invoke drawing of object on TCanvas sub-pad
/// All existing providers are checked, first checked are class matches (including direct parents)

bool RProvider::Draw6(TVirtualPad *subpad, std::unique_ptr<RHolder> &object, const std::string &opt)
{
   if (!object || !object->GetClass())
      return false;

   auto draw_func = [subpad, &object, &opt](Draw6Func_t &func) -> bool {
      return func(subpad, object, opt);
   };

   if (ScanProviderMap<Draw6Map_t,Draw6Func_t>(GetDraw6Map(), object->GetClass(), false, draw_func))
      return true;

   auto &entry = GetClassEntry(object->GetClass());
   if (!entry.dummy() && !entry.draw6lib.empty())
      gSystem->Load(entry.draw6lib.c_str());

   return ScanProviderMap<Draw6Map_t,Draw6Func_t>(GetDraw6Map(), object->GetClass(), true, draw_func);
}

/////////////////////////////////////////////////////////////////////////////////
/// Invoke drawing of object on RCanvas sub-pad
/// All existing providers are checked, first checked are class matches (including direct parents)

bool RProvider::Draw7(std::shared_ptr<ROOT::Experimental::RPadBase> &subpad, std::unique_ptr<RHolder> &object, const std::string &opt)
{
   if (!object || !object->GetClass())
      return false;

   auto draw_func = [&subpad, &object, &opt](Draw7Func_t &func) -> bool {
      return func(subpad, object, opt);
   };

   if (ScanProviderMap<Draw7Map_t,Draw7Func_t>(GetDraw7Map(), object->GetClass(), false, draw_func))
      return true;

   auto &entry = GetClassEntry(object->GetClass());
   if (!entry.dummy() && !entry.draw7lib.empty())
      gSystem->Load(entry.draw7lib.c_str());

   return ScanProviderMap<Draw7Map_t,Draw7Func_t>(GetDraw7Map(), object->GetClass(), true, draw_func);
}

/////////////////////////////////////////////////////////////////////
/// Return icon name for the given class - either class name or TClass *

std::string RProvider::GetClassIcon(const ClassArg &arg, bool is_folder)
{
   auto &entry = GetClassEntry(arg);
   if (!entry.iconname.empty())
      return entry.iconname;

   return is_folder ? "sap-icon://folder-blank"s : "sap-icon://electronic-medical-record"s;
}


/////////////////////////////////////////////////////////////////////
/// Return true if provided class can have childs

bool RProvider::CanHaveChilds(const ClassArg &arg)
{
   return GetClassEntry(arg).can_have_childs;
}

/////////////////////////////////////////////////////////////////////
/// Return true if provided class can be drawn on the TCanvas

bool RProvider::CanDraw6(const ClassArg &arg)
{
   if (ScanProviderMap<Draw6Map_t,Draw6Func_t>(GetDraw6Map(), arg))
      return true;

   if (!GetClassEntry(arg).draw6lib.empty())
      return true;

   return false;
}

/////////////////////////////////////////////////////////////////////
/// Return true if provided class can be drawn on the RCanvas

bool RProvider::CanDraw7(const ClassArg &arg)
{
   if (ScanProviderMap<Draw7Map_t,Draw7Func_t>(GetDraw7Map(), arg))
      return true;

   if (!GetClassEntry(arg).draw7lib.empty())
      return true;

   return false;
}


// ==============================================================================================

class RDefaultProvider : public RProvider {

public:
   RDefaultProvider()
   {
      // TODO: let read from rootrc or any other files
      RegisterClass("ROOT::Experimental::RH1D", "sap-icon://bar-chart", "", "", "libROOTHistDrawProvider");
      RegisterClass("ROOT::Experimental::RH2D", "sap-icon://pixelate", "", "", "libROOTHistDrawProvider");
      RegisterClass("ROOT::Experimental::RH3D", "sap-icon://product", "", "", "libROOTHistDrawProvider");
   }

} newRDefaultProvider;

