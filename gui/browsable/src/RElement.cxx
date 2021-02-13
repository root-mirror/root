/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RElement.hxx>

#include <ROOT/Browsable/RLevelIter.hxx>
#include <ROOT/RLogger.hxx>
#include "TBufferJSON.h"

using namespace ROOT::Experimental::Browsable;
using namespace std::string_literals;

ROOT::Experimental::RLogChannel &ROOT::Experimental::BrowsableLog() {
   static RLogChannel sLog("ROOT.Browsable");
   return sLog;
}


/////////////////////////////////////////////////////////////////////
/// Returns child iterator (if any)

std::unique_ptr<RLevelIter> RElement::GetChildsIter()
{
   return nullptr;
}

/////////////////////////////////////////////////////////////////////
/// Returns number of childs
/// By default creates iterator and iterates over all items

int RElement::GetNumChilds()
{
   auto iter = GetChildsIter();
   if (!iter) return 0;
   int cnt = 0;
   while (iter->Next()) cnt++;
   return cnt;
}

/////////////////////////////////////////////////////////////////////
/// Find item with specified name
/// Default implementation, should work for all

RElement::EContentKind RElement::GetContentKind(const std::string &kind)
{
    std::string lkind = kind;
    std::transform(lkind.begin(), lkind.end(), lkind.begin(), ::tolower);

   if (lkind == "text") return kText;
   if ((lkind == "image") || (lkind == "image64")) return kImage;
   if (lkind == "png") return kPng;
   if ((lkind == "jpg") || (lkind == "jpeg")) return kJpeg;
   if (lkind == "json") return kJson;
   if (lkind == "filename") return kFileName;
   return kNone;
}

/////////////////////////////////////////////////////////////////////
/// Returns sub element

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
/// Returns string content like text file content or json representation

std::string RElement::GetContent(const std::string &kind)
{
   if (GetContentKind(kind) == kJson) {
      auto obj = GetObject();
      if (obj)
         return TBufferJSON::ConvertToJSON(obj->GetObject(), obj->GetClass()).Data();
   }

   return ""s;
}

/////////////////////////////////////////////////////////////////////
/// Parse string path to produce RElementPath_t
/// One should avoid to use string pathes as much as possible

RElementPath_t RElement::ParsePath(const std::string &strpath)
{
   RElementPath_t arr;
   if (strpath.empty())
      return arr;

   std::string slash = "/";

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

/////////////////////////////////////////////////////////////////////
/// Compare two paths,
/// Returns number of elements matches in both paths

int RElement::ComparePaths(const RElementPath_t &path1, const RElementPath_t &path2)
{
   int sz = path1.size();
   if (sz > (int) path2.size()) sz = path2.size();

   for (int n = 0; n < sz; ++n)
      if (path1[n] != path2[n])
         return n;

   return sz;
}

/////////////////////////////////////////////////////////////////////
/// Converts element path back to string

std::string RElement::GetPathAsString(const RElementPath_t &path)
{
   std::string res;
   for (auto &elem : path) {
      res.append("/");
      std::string subname = elem;
      ExtractItemIndex(subname);
      res.append(subname);
   }

   return res;
}

/////////////////////////////////////////////////////////////////////
/// Extract index from name
/// Index coded by client with ###<indx>$$$ suffix
/// Such coding used by browser to identify element by index

int RElement::ExtractItemIndex(std::string &name)
{
   auto p1 = name.rfind("###"), p2 = name.rfind("$$$");
   if ((p1 == std::string::npos) || (p2 == std::string::npos) || (p1 >= p2) || (p2 != name.length()-3)) return -1;

   int indx = std::stoi(name.substr(p1+3,p2-p1-3));
   name.resize(p1);
   return indx;
}
