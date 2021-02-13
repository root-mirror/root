// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2019-10-14
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RBrowserData.hxx>

#include <ROOT/Browsable/RProvider.hxx>
#include <ROOT/Browsable/RLevelIter.hxx>
#include <ROOT/RLogger.hxx>

#include "TBufferJSON.h"

#include <algorithm>
#include <regex>

using namespace ROOT::Experimental;
using namespace std::string_literals;

ROOT::Experimental::RLogChannel &ROOT::Experimental::BrowserLog() {
   static RLogChannel sLog("ROOT.Browser");
   return sLog;
}


/** \class ROOT::Experimental::RBrowserData
\ingroup rbrowser
\brief Way to browse (hopefully) everything in %ROOT
*/


/////////////////////////////////////////////////////////////////////
/// set top element for browsing

void RBrowserData::SetTopElement(std::shared_ptr<Browsable::RElement> elem)
{
   fTopElement = elem;

   SetWorkingPath({});
}

/////////////////////////////////////////////////////////////////////
/// set working directory relative to top element

void RBrowserData::SetWorkingPath(const Browsable::RElementPath_t &path)
{
   fWorkingPath = path;

   ResetLastRequestData(true);
}

/////////////////////////////////////////////////////////////////////
/// Reset all data correspondent to last request

void RBrowserData::ResetLastRequestData(bool with_element)
{
   fLastAllChilds = false;
   fLastSortedItems.clear();
   fLastSortMethod.clear();
   fLastItems.clear();
   if (with_element) {
      fLastPath.clear();
      fLastElement.reset();
   }
}

/////////////////////////////////////////////////////////////////////////
/// Decompose path to elements
/// Returns array of names for each element in the path, first element either "/" or "."
/// If returned array empty - it is error

Browsable::RElementPath_t RBrowserData::DecomposePath(const std::string &strpath, bool relative_to_work_element)
{
   Browsable::RElementPath_t arr;
   if (relative_to_work_element) arr = fWorkingPath;

   if (strpath.empty())
      return arr;

   auto arr2 = Browsable::RElement::ParsePath(strpath);
   arr.insert(arr.end(), arr2.begin(), arr2.end());
   return arr;
}

/////////////////////////////////////////////////////////////////////////
/// Process browser request

bool RBrowserData::ProcessBrowserRequest(const RBrowserRequest &request, RBrowserReply &reply)
{
   auto path = fWorkingPath;
   path.insert(path.end(), request.path.begin(), request.path.end());

   if ((path != fLastPath) || !fLastElement) {

      auto elem = GetSubElement(path);
      if (!elem) return false;

      ResetLastRequestData(true);

      fLastPath = path;
      fLastElement = std::move(elem);
   } else if (request.reload) {
      // only reload items from element, not need to reset element itself
      ResetLastRequestData(false);
   }

   // when request childs, always try to make elements
   if (fLastItems.empty()) {

      auto iter = fLastElement->GetChildsIter();

      if (!iter) return false;
      int id = 0;
      fLastAllChilds = true;

      while (iter->Next() && fLastAllChilds) {
         fLastItems.emplace_back(iter->CreateItem());
         if (id++ > 10000)
            fLastAllChilds = false;
      }

      fLastSortedItems.clear();
      fLastSortMethod.clear();
   }

   // create sorted array
   if ((fLastSortedItems.size() != fLastItems.size()) ||
       (fLastSortMethod != request.sort) ||
       (fLastSortReverse != request.reverse)) {
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
                      [request](const Browsable::RItem *a, const Browsable::RItem *b) { return a->Compare(b, request.sort); });
      }

      if (request.reverse)
         std::reverse(fLastSortedItems.begin(), fLastSortedItems.end());

      fLastSortMethod = request.sort;
      fLastSortReverse = request.reverse;
   }

   const std::regex expr(request.regex);

   int id = 0;
   for (auto &item : fLastSortedItems) {

      // check if element is hidden
      if (!request.hidden && item->IsHidden())
         continue;

      if (!request.regex.empty() && !item->IsFolder() && !std::regex_match(item->GetName(), expr))
         continue;

      if ((id >= request.first) && ((request.number == 0) || (id < request.first + request.number)))
         reply.nodes.emplace_back(item);

      id++;
   }

   reply.first = request.first;
   reply.nchilds = id; // total number of childs

   return true;
}

/////////////////////////////////////////////////////////////////////////
/// Process browser request, returns string with JSON of RBrowserReply data

std::string RBrowserData::ProcessRequest(const RBrowserRequest &request)
{
   RBrowserReply reply;

   reply.path = request.path;
   reply.first = 0;
   reply.nchilds = 0;

   ProcessBrowserRequest(request, reply);

   return TBufferJSON::ToJSON(&reply, TBufferJSON::kSkipTypeInfo + TBufferJSON::kNoSpaces).Data();
}

/////////////////////////////////////////////////////////////////////////
/// Returns element with path, specified as string

std::shared_ptr<Browsable::RElement> RBrowserData::GetElement(const std::string &str)
{
   auto path = DecomposePath(str, true);

   return GetSubElement(path);
}

/////////////////////////////////////////////////////////////////////////
/// Returns element with path, specified as Browsable::RElementPath_t

std::shared_ptr<Browsable::RElement> RBrowserData::GetElementFromTop(const Browsable::RElementPath_t &path)
{
   return GetSubElement(path);
}

/////////////////////////////////////////////////////////////////////////
/// Returns sub-element starting from top, using cached data

std::shared_ptr<Browsable::RElement> RBrowserData::GetSubElement(const Browsable::RElementPath_t &path)
{
   if (path.empty())
      return fTopElement;

   // find best possible entry in cache
   int pos = 0;
   auto elem = fTopElement;

   for (auto &entry : fCache) {
      auto comp = Browsable::RElement::ComparePaths(path, entry.first);
      if (comp > pos) { pos = comp; elem = entry.second; }
   }

   while (pos < (int) path.size()) {
      std::string subname = path[pos];
      int indx = Browsable::RElement::ExtractItemIndex(subname);

      auto iter = elem->GetChildsIter();
      if (!iter)
         return nullptr;

      if (!iter->Find(subname, indx)) {
         if (indx < 0)
            return nullptr;
         iter = elem->GetChildsIter();
         if (!iter || !iter->Find(subname))
            return nullptr;
      }

      elem = iter->GetElement();

      if (!elem)
         return nullptr;

      auto subpath = path;
      subpath.resize(pos+1);
      fCache[subpath] = elem;

      pos++; // switch to next element
   }

   return elem;
}

