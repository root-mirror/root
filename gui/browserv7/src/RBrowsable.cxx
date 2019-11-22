/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RBrowsable.hxx"

#include "ROOT/RLogger.hxx"

#include "TClass.h"
#include "TBufferJSON.h"

#include <algorithm>

using namespace ROOT::Experimental;
using namespace ROOT::Experimental::Browsable;
using namespace std::string_literals;


/////////////////////////////////////////////////////////////////////
/// Find item with specified name
/// Default implementation, should work for all

RElement::EContentKind RElement::GetContentKind(const std::string &kind)
{
   if (kind == "text") return kText;
   if ((kind == "image") || (kind == "image64")) return kImage;
   if (kind == "png") return kPng;
   if ((kind == "jpg") || (kind == "jpeg")) return kJpeg;
   return kNone;
}


std::shared_ptr<RElement> RElement::GetSubElement(std::shared_ptr<RElement> &elem, const RElementPath_t &path)
{
   auto curr = elem;

   for (auto &itemname : path) {
      if (!curr)
         return nullptr;

      auto iter = curr->GetChildsIter();
      if (!iter || !iter->Find(itemname))
         return nullptr;

      curr = iter->GetElement();
   }

   return curr;
}


/////////////////////////////////////////////////////////////////////
/// Find item with specified name
/// Default implementation, should work for all

bool RLevelIter::Find(const std::string &name)
{
   if (!Reset()) return false;

   while (Next()) {
      if (GetName() == name)
         return true;
   }

   return false;
}



std::string RProvider::GetClassIcon(const std::string &classname)
{
   if (classname == "TTree" || classname == "TNtuple")
      return "sap-icon://tree"s;
   if (classname == "TDirectory" || classname == "TDirectoryFile")
      return "sap-icon://folder-blank"s;
   if (classname.find("TLeaf") == 0)
      return "sap-icon://e-care"s;

   return "sap-icon://electronic-medical-record"s;
}


RProvider::BrowseMap_t &RProvider::GetBrowseMap()
{
   static RProvider::BrowseMap_t sMap;
   return sMap;
}

RProvider::FileMap_t &RProvider::GetFileMap()
{
   static RProvider::FileMap_t sMap;
   return sMap;
}


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
}


void RProvider::RegisterFile(const std::string &extension, FileFunc_t func)
{
    auto &fmap = GetFileMap();

    if ((extension != "*") && (fmap.find(extension) != fmap.end()))
       R__ERROR_HERE("Browserv7") << "Provider for file extension  " << extension << " already exists";

    fmap.emplace(extension, StructFile{this,func});
}

void RProvider::RegisterBrowse(const TClass *cl, BrowseFunc_t func)
{
    auto &bmap = GetBrowseMap();

    if (cl && (bmap.find(cl) != bmap.end()))
       R__ERROR_HERE("Browserv7") << "Browse provider for class " << cl->GetName() << " already exists";

    bmap.emplace(cl, StructBrowse{this,func});
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

std::shared_ptr<RElement> RProvider::Browse(std::unique_ptr<Browsable::RHolder> &object)
{
   auto &bmap = GetBrowseMap();

   auto cl = object->GetClass();

   auto iter = bmap.find(cl);

   if (iter != bmap.end()) {
      auto res = iter->second.func(object);
      if (res || !object) return res;
   }

   for (auto &pair : bmap)
      if ((pair.first == nullptr) || (pair.first == cl)) {
         auto res = pair.second.func(object);
         if (res || !object) return res;
      }

   return nullptr;
}

/////////////////////////////////////////////////////////////////////
/// set top element for browsing

void RBrowsable::SetTopElement(std::shared_ptr<Browsable::RElement> elem)
{
   fTopElement = elem;

   SetWorkingDirectory("");
}

/////////////////////////////////////////////////////////////////////
/// set working directory relative to top element

void RBrowsable::SetWorkingDirectory(const std::string &strpath)
{
   auto path = DecomposePath(strpath, false);

   SetWorkingPath(path);
}

/////////////////////////////////////////////////////////////////////
/// set working directory relative to top element

void RBrowsable::SetWorkingPath(const RElementPath_t &path)
{
   fWorkingPath = path;
   fWorkElement = RElement::GetSubElement(fTopElement, path);

   ResetLastRequest();

   fLevels.clear();
   fLevels.emplace_back("");
   fLevels.front().fElement = elem;
}

/////////////////////////////////////////////////////////////////////
/// Navigate to the top level

bool RBrowsable::ResetLevels()
{
   while (fLevels.size() > 1)
      fLevels.pop_back();

   if (fLevels.size() != 1)
      return false;

   return true;
}


/////////////////////////////////////////////////////////////////////
/// Reset all data correspondent to last request

void RBrowsable::ResetLastRequest()
{
   fLastAllChilds = false;
   fLastSortedItems.clear();
   fLastSortMethod.clear();
   fLastItems.clear();
   fLastPath.clear();
   fLastElement.reset();
}

/////////////////////////////////////////////////////////////////////
/// Direct navigate to specified path without building levels
/// Do not support ".." - navigation to level up

std::shared_ptr<Browsable::RElement> RBrowsable::DirectNavigate(std::shared_ptr<Browsable::RElement> item, const RElementPath_t &paths, int indx)
{
   while (item && (indx < (int) paths.size())) {
      auto subdir = paths[indx++];

      if (subdir == ".") continue;

      // not supported with direct navigate
      if (subdir == "..") return nullptr;

      auto iter = item->GetChildsIter();

      if (!iter || !iter->Find(subdir))
         return nullptr;

      item = iter->GetElement();
   }
   return item;
}


/////////////////////////////////////////////////////////////////////
/// Navigate to specified path
/// If specified paths array empty - just do nothing
/// If first element is "/", start from the most-top element
/// If first element ".", start from current position
/// Any other will be ignored
/// level_indx returns found index in the levels. If pointer not specified, levels will not be touched

std::shared_ptr<Browsable::RElement> RBrowsable::Navigate(const RElementPath_t &paths, int *level_indx)
{
   if (fLevels.empty())
      return nullptr;

   if (paths.empty())
      return fLevels.back().fElement;

   int lindx = 0;
   if (paths[0] == "/")
      lindx = 0;
   else if (paths[0] == ".")
      lindx = (int)fLevels.size() - 1;
   else
      return nullptr;

   for (int pindx = 1; pindx < (int) paths.size(); pindx ++) {

      auto subdir = paths[pindx];
      if (subdir == ".") continue; // do nothing

      if (subdir == "..") { // one level up
         if (--lindx < 0)
            return nullptr;
         continue;
      }

      // path is matching
      if ((lindx + 1 < (int) fLevels.size()) && (fLevels[lindx + 1].fName == subdir)) {
         ++lindx;
         continue;
      }

      auto &level = fLevels[lindx];

      // up to here we could follow existing paths
      // but here , any next elements should be created directly
      if (!level_indx)
         return DirectNavigate(level.fElement, paths, pindx);

      auto iter = level.fElement->GetChildsIter();

      if (!iter || !iter->Find(subdir))
         return nullptr;

      auto subitem = iter->GetElement();
      if (!subitem)
         return nullptr;

      fLevels.resize(++lindx); // remove all other items which may be there before

      fLevels.emplace_back(subdir, subitem);
   }

   if (level_indx) *level_indx = lindx;

   return fLevels[lindx].fElement;
}

/////////////////////////////////////////////////////////////////////////
/// Decompose path to elements
/// Returns array of names for each element in the path, first element either "/" or "."
/// If returned array empty - it is error

RElementPath_t RBrowsable::DecomposePath(const std::string &strpath, bool relative_path)
{
   RElementPath_t arr;

   if (strpath.empty()) {
      if (relative_path)
         arr.emplace_back("."s);
      return arr;
   }

   std::string slash = "/";

   if (relative_path) arr.emplace_back(slash);

   std::string::size_type previous = 0;
   if (strpath[0] == slash[0]) previous++;

   auto current = strpath.find(slash, previous);
   while (current != std::string::npos) {
      if (current > previous)
         arr.emplace_back(strpath.substr(previous, current - previous));
      previous = current + 1;
      current = strpath.find(slash, previous);
   }

   if (previous < strpath.length())
      arr.emplace_back(strpath.substr(previous));

   return arr;
}


/////////////////////////////////////////////////////////////////////////
/// Process browser request

bool RBrowsable::SamePath(const RElementPath_t &p1, const RElementPath_t &p2) const
{
   if (p1.size() != p2.size())
      return false;

   for (unsigned n = 0; n < p1.size(); ++n)
      if (p1[n] != p2[n])
         return false;

   return true;
}


/////////////////////////////////////////////////////////////////////////
/// Process browser request

bool RBrowsable::ProcessRequest(const RBrowserRequest &request, RBrowserReply &reply)
{
   if (gDebug > 0)
      printf("REQ: Do decompose path '%s'\n",request.path.c_str());

   auto arr = DecomposePath(request.path, false);

   if (!SamePath(arr, fLastPath) || !fLastElement) {

      auto elem = RElement::GetSubElement(fWorkElement, arr);
      if (!elem) return false;

      ResetLastRequest();

      fLastPath = arr;
      fLastElement = elem;
   }

   if (gDebug > 0) {
      printf("REQ:Try to navigate %d\n", (int) arr.size());
      for (auto & subdir : arr) printf("   %s\n", subdir.c_str());
   }

   // when request childs, always try to make elements
   if (fLastItems.empty()) {
      auto iter = fLastElement->GetChildsIter();
      if (!iter) return false;
      int id = 0;
      fLastAllChilds = true;

      while (iter->Next() && fLastAllChilds) {
         fLastItems.emplace_back(iter->CreateBrowserItem());
         if (id++ > 10000)
            fLastAllChilds = false;
      }

      fLastSortedItems.clear();
      fLastSortMethod.clear();
   }

   // create sorted array
   if ((fLastSortedItems.size() != fLastItems.size()) || (fLastSortMethod != request.sort)) {
      fLastSortedItems.resize(fLastItems.size(), nullptr);
      int id = 0;
      if (request.sort.empty()) {
         // no sorting, just move all folders up
         for (auto &item : fLastItems)
            if (item->IsFolder())
               fLastSortedItems[id++] = item.get();
         for (auto &item : fLastItems)
            if (!item->IsFolder())
               fLastSortedItems[id++] = item.get();
      } else {
         // copy items
         for (auto &item : fLastItems)
            fLastSortedItems[id++] = item.get();

         if (request.sort != "unsorted")
            std::sort(fLastSortedItems.begin(), fLastSortedItems.end(),
                      [request](const RBrowserItem *a, const RBrowserItem *b) { return a->Compare(b, request.sort); });
      }
      fLastSortMethod = request.sort;
   }

   int id = 0;
   for (auto &item : fLastSortedItems) {
      if (!request.filter.empty() && (item->GetName().compare(0, request.filter.length(), request.filter) != 0))
         continue;

      if ((id >= request.first) && ((request.number == 0) || (id < request.first + request.number)))
         reply.nodes.emplace_back(item);
      id++;
   }

   reply.first = request.first;
   reply.nchilds = id; // total number of childs

   return true;
}


std::string RBrowsable::ProcessRequest(const RBrowserRequest &request)
{
   RBrowserReply reply;

   reply.path = request.path;
   reply.first = 0;
   reply.nchilds = 0;

   ProcessRequest(request, reply);

   return TBufferJSON::ToJSON(&reply, TBufferJSON::kSkipTypeInfo + TBufferJSON::kNoSpaces).Data();
}


std::shared_ptr<RElement> RBrowsable::GetElement(const std::string &path)
{
   auto arr = DecomposePath(path);
   if (arr.empty()) return nullptr;

   // find element using current levels when possible
   return Navigate(arr);
}


