// @(#)root/io:$Id$
// Author: Jakob Blomer

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RConfig.h>
#include <ROOT/RRawFile.hxx>

#include "TPluginManager.h"
#include "TROOT.h"

#include <memory>

std::unique_ptr<ROOT::Internal::RRawFile> ROOT::Internal::RRawFile::CreateHttpRRawFile(std::string_view url, ROptions options)
{  
   std::string transport = GetTransport(url);
   std::unique_ptr<RRawFile>  helper;
   if (transport == "http" || transport == "https") {
      if (TPluginHandler *h = gROOT->GetPluginManager()->FindHandler("ROOT::Internal::RRawFile")) {
         if (h->LoadPlugin() == 0) {
            //return std::unique_ptr<RRawFile>(reinterpret_cast<RRawFile *>(h->ExecPlugin(2, &url, &options)));
            helper = std::unique_ptr<RRawFile>(reinterpret_cast<RRawFile *>(h->ExecPlugin(2, &url, &options)));
            ROOT::Internal::RRawFile::InitHelper(std::move(helper));
            return helper;
         }
         throw std::runtime_error("Cannot load plugin handler for RRawFileDavix");
      }
      throw std::runtime_error("Cannot find plugin handler for RRawFileDavix");
   } else {
         return nullptr;
   }
}
