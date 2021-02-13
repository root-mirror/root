// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2021-01-22
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RBrowserWidget.hxx"

#include <algorithm>

#include "TSystem.h"


using namespace std::string_literals;
using namespace ROOT::Experimental;

RBrowserWidgetProvider::RBrowserWidgetProvider(const std::string &kind)
{
   GetMap().emplace(kind, this);
}

RBrowserWidgetProvider::~RBrowserWidgetProvider()
{
   auto &map = GetMap();
   auto iter = std::find_if(map.begin(), map.end(),
         [this](const ProvidersMap_t::value_type &pair) { return this == pair.second; });

   if (iter != map.end())
      map.erase(iter);
}

RBrowserWidgetProvider::ProvidersMap_t& RBrowserWidgetProvider::GetMap()
{
   static RBrowserWidgetProvider::ProvidersMap_t mMap;
   return mMap;
}

///////////////////////////////////////////////////////////////
/// Create specified widget

std::shared_ptr<RBrowserWidget> RBrowserWidgetProvider::CreateWidget(const std::string &kind, const std::string &name)
{
   auto &map = GetMap();
   auto iter = map.find(kind);
   if (iter == map.end()) {
      // try to load necessary libraries
      if (kind == "geom")
         gSystem->Load("libROOTBrowserGeomWidget");
      else if (kind == "tcanvas")
         gSystem->Load("libROOTBrowserTCanvasWidget");
      else if (kind == "rcanvas")
         gSystem->Load("libROOTBrowserRCanvasWidget");
      iter = map.find(kind);
      if (iter == map.end())
         return nullptr;
   }
   return iter->second->Create(name);
}
