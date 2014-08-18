// @(#)root/utils:$Id$
// Author: Axel Naumann, 2014-04-07

/*************************************************************************
 * Copyright (C) 1995-2014, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Provides bindings to TCling (compiled with rtti) from rootcling (compiled
// without rtti).

#include "rootclingTCling.h"

#include "TClass.h"
#include "TCling.h"
#include "TEnum.h"
#include "TFile.h"
#include "TProtoClass.h"
#include "TROOT.h"
#include "TStreamerInfo.h"
#include <iostream>

std::string gPCMFilename;
std::vector<std::string> gClassesToStore;
std::vector<std::string> gTypedefsToStore;
std::vector<std::string> gEnumsToStore;
std::vector<std::string> gAncestorPCMsNames;

extern "C"
const char ** *TROOT__GetExtraInterpreterArgs()
{
   return &TROOT::GetExtraInterpreterArgs();
}

extern "C"
cling::Interpreter *TCling__GetInterpreter()
{
   static bool sInitialized = false;
   gROOT; // trigger initialization
   if (!sInitialized) {
      gCling->SetClassAutoloading(false);
      sInitialized = true;
   }
   return ((TCling *)gCling)->GetInterpreter();
}

extern "C"
void InitializeStreamerInfoROOTFile(const char *filename)
{
   gPCMFilename = filename;
}

extern "C"
void AddStreamerInfoToROOTFile(const char *normName)
{
   gClassesToStore.emplace_back(normName);
}

extern "C"
void AddTypedefToROOTFile(const char *tdname)
{
   gTypedefsToStore.push_back(tdname);
}

extern "C"
void AddEnumToROOTFile(const char *enumname)
{
   gEnumsToStore.push_back(enumname);
}

extern "C"
void AddAncestorPCMROOTFile(const char *pcmName)
{
   gAncestorPCMsNames.emplace_back(pcmName);
}

extern "C"
bool CloseStreamerInfoROOTFile()
{
   // Write all persistent TClasses.

   // Avoid plugins.
   TVirtualStreamerInfo::SetFactory(new TStreamerInfo());

   TObjArray protoClasses;
   for (const auto & normName : gClassesToStore) {
      TClass *cl = TClass::GetClass(normName.c_str(), kTRUE /*load*/);
      if (!cl) {
         std::cerr << "ERROR in CloseStreamerInfoROOTFile(): cannot find class "
                   << normName << '\n';
         return false;
      }
      // We include transient classes as they could be used by a derived
      // class which may have rules setting the member of the transient class.
      // (And the derived class RealData *do* contain member from the transient
      // base classes.
//      if (cl->GetClassVersion() == 0)
//         continue;

      // Let's include also proxied collections in order to delay parsing as long as possible.
      // In the first implementations, proxied collections did not result in a protoclass.
      // If this is a proxied collection then offsets are not needed.
//       if (cl->GetCollectionProxy())
//          continue;
      cl->Property(); // Force initialization of the bits and property fields.

      protoClasses.AddLast(new TProtoClass(cl));
   }

   TObjArray typedefs;
   for (const auto & dtname : gTypedefsToStore) {
      TDataType *dt = (TDataType *)gROOT->GetListOfTypes()->FindObject(dtname.c_str());
      if (!dt) {
         std::cerr << "ERROR in CloseStreamerInfoROOTFile(): cannot find typedef "
                   << dtname << '\n';
         return false;
      }
      if (dt->GetType() == -1) {
         dt->Property(); // Force initialization of the bits and property fields.
         dt->GetTypeName(); // Force caching of type name.
         typedefs.AddLast(dt);
      }
   }


   TObjArray enums;
   for (const auto & enumname : gEnumsToStore) {
      TEnum *en = nullptr;
      const size_t lastSepPos = enumname.find_last_of("::");
      if (lastSepPos != std::string::npos) {
         const std::string nsName = enumname.substr(0, lastSepPos - 1);
         TClass *tclassInstance = TClass::GetClass(nsName.c_str());
         if (!tclassInstance) {
            std::cerr << "ERROR in CloseStreamerInfoROOTFile(): cannot find TClass instance for namespace "
                      << nsName << '\n';
            return false;
         }
         auto enumListPtr = tclassInstance->GetListOfEnums();
         if (!enumListPtr) {
            std::cerr << "ERROR in CloseStreamerInfoROOTFile(): TClass instance for namespace "
                      << nsName << " does not have any enum associated. This is an inconsistency." << '\n';
            return false;
         }
         const std::string unqualifiedEnumName = enumname.substr(lastSepPos + 1);
         en = (TEnum *)enumListPtr->FindObject(unqualifiedEnumName.c_str());
         en->SetTitle(nsName.c_str());
      } else {
         en = (TEnum *)gROOT->GetListOfEnums()->FindObject(enumname.c_str());
         en->SetTitle("");
      }
      if (!en) {
         std::cerr << "ERROR in CloseStreamerInfoROOTFile(): cannot find enum "
                   << enumname << '\n';
         return false;
      }
      en->Property(); // Force initialization of the bits and property fields.
      enums.AddLast(en);
   }

   // Don't use TFile::Open(); we don't need plugins.
   TFile dictFile(gPCMFilename.c_str(), "RECREATE");
   if (dictFile.IsZombie())
      return false;
   // Instead of plugins:
   protoClasses.Write("__ProtoClasses", TObject::kSingleKey);
   protoClasses.Delete();
   typedefs.Write("__Typedefs", TObject::kSingleKey);
   enums.Write("__Enums", TObject::kSingleKey);

   dictFile.WriteObjectAny(&gAncestorPCMsNames, "std::vector<std::string>", "__AncestorPCMsNames");


   return true;
}
