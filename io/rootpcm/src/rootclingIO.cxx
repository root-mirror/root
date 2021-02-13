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

#include "TClass.h"
#include "TEnum.h"
#include "TError.h"
#include "TFile.h"
#include "TProtoClass.h"
#include "TDataMember.h"
#include "TROOT.h"
#include "TStreamerInfo.h"
#include "TClassEdit.h"
#include <memory>

std::string gPCMFilename;
std::vector<std::string> gClassesToStore;
std::vector<std::string> gTypedefsToStore;
std::vector<std::string> gEnumsToStore;

extern "C"
void InitializeStreamerInfoROOTFile(const char *filename)
{
   gPCMFilename = filename;
}

extern "C"
void AddStreamerInfoToROOTFile(const char *normName)
{
   // Filter unnamed and (anonymous) classes.
   if (normName && normName[0] && normName[0] != '(')
      gClassesToStore.emplace_back(normName);
}

extern "C"
void AddTypedefToROOTFile(const char *tdname)
{
   gTypedefsToStore.emplace_back(tdname);
}

extern "C"
void AddEnumToROOTFile(const char *enumname)
{
   gEnumsToStore.emplace_back(enumname);
}

static bool IsUniquePtrOffsetZero()
{
   auto regularPtr = (long *)0x42;
   std::unique_ptr<long> uniquePtr(regularPtr);
   auto regularPtr_2 = reinterpret_cast<long **>(&uniquePtr);
   bool isZero = uniquePtr.get() == *regularPtr_2;
   uniquePtr.release();
   if (!isZero) {
      Error("CloseStreamerInfoROOTFile",
            "UniquePtr points to %p, reinterpreting it gives %p and should have %p", uniquePtr.get(), *(regularPtr_2), regularPtr);
   }
   return isZero;
}

static bool IsUnsupportedUniquePointer(const char *normName, TDataMember *dm)
{
   auto dmTypeName = dm->GetTypeName();
   static bool isUniquePtrOffsetZero = IsUniquePtrOffsetZero(); // call this only once
   if (TClassEdit::IsUniquePtr(dmTypeName)) {

      if (!isUniquePtrOffsetZero) return true;

      auto clm = TClass::GetClass(dmTypeName);
      if (!clm) {
         Error("CloseStreamerInfoROOTFile", "Class %s is not available.", dmTypeName);
         return true;
      }

      clm->BuildRealData();
      auto upDms = clm->GetListOfRealData();
      if (!upDms) {
         Error("CloseStreamerInfoROOTFile", "Cannot determine unique pointer %s data members.", dmTypeName);
         return true;
      }

      if (0 == upDms->GetSize()) {
         Error("CloseStreamerInfoROOTFile", "Unique pointer %s has zero data members.", dmTypeName);
         return true;
      }

      // We check if the unique_ptr has a default deleter
      std::vector<std::string> out;
      int i;
      TClassEdit::GetSplit(dmTypeName, out, i);
      // out[2] contains the deleter type.
      if (0 != out[2].find("default_delete<")) {
         Error("CloseStreamerInfoROOTFile", "I/O is supported only for unique_ptrs with a default deleter. %s::%s  appears to have a custom one, %s.", normName, dm->GetName(), out[2].c_str());
         return true;
      }
   }
   return false;
}

// static bool IsSupportedClass(TClass *cl)
// {
//    // Check if the Class is of an unsupported type
//    using namespace ROOT::TMetaUtils;
//
//    // Check if this is a collection of unique_ptrs
//    if (ROOT::ESTLType::kNotSTL != cl->GetCollectionType()) {
//       std::vector<std::string> out;
//       int i;
//       TClassEdit::GetSplit(cl->GetName(), out, i);
//       std::string_view containedObjectTypeName(out[1].c_str());
//       if (TClassEdit::IsUniquePtr(containedObjectTypeName)) {
//          auto clName = cl->GetName();
//          // Here we can use the new name for the error message
//          Error("CloseStreamerInfoROOTFile", "A collection of unique pointers was selected: %s. These are not supported. If you wish to perform I/O operations with %s, just select the same collection of raw C pointers.\n", clName, clName);
//          return false;
//       }
//    }
//    return true;
//
// }

extern "C"
bool CloseStreamerInfoROOTFile(bool writeEmptyRootPCM)
{
   // Write all persistent TClasses.

   // Avoid plugins.
   TVirtualStreamerInfo::SetFactory(new TStreamerInfo());

   // Don't use TFile::Open(); we don't need plugins.
   TFile dictFile((gPCMFilename + "?filetype=pcm").c_str(), "RECREATE");

   // Reset the content of the pcm
   if (writeEmptyRootPCM) {
      TObject obj;
      obj.Write("EMPTY");
      return true;
   };

   TObjArray protoClasses(gClassesToStore.size());
   for (const auto & normName : gClassesToStore) {
      TClass *cl = TClass::GetClass(normName.c_str(), kTRUE /*load*/);
      if (!cl) {
         Error("CloseStreamerInfoROOTFile", "Cannot find class %s.", normName.c_str());
         return false;
      }

//       if (!IsSupportedClass(cl)) return false;

      // Check if a datamember is a unique_ptr and if yes that it has a default
      // deleter.
      auto dms = cl->GetListOfDataMembers();
      if (!dms) {
         Error("CloseStreamerInfoROOTFile", "Cannot find data members for %s.", normName.c_str());
         return false;
      }

      for (auto dmObj : *dms) {
         auto dm = (TDataMember *) dmObj;
         if (!dm->IsPersistent() || cl->GetClassVersion()==0) continue;
         if (IsUnsupportedUniquePointer(normName.c_str(), dm)) return false;
      }

      // Never store a proto class for a class which rootcling already has
      // an 'official' TClass (i.e. the dictionary is in libCore or libRIO).
      if (cl->IsLoaded()) continue;

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

   TObjArray typedefs(gTypedefsToStore.size());

   for (const auto & dtname : gTypedefsToStore) {
      TDataType *dt = (TDataType *)gROOT->GetListOfTypes()->FindObject(dtname.c_str());
      if (!dt) {
         Error("CloseStreamerInfoROOTFile", "Cannot find typedef %s.", dtname.c_str());
         return false;
      }
      if (dt->GetType() == -1) {
         dt->Property(); // Force initialization of the bits and property fields.
         dt->GetTypeName(); // Force caching of type name.
         typedefs.AddLast(dt);
      }
   }


   TObjArray enums(gEnumsToStore.size());
   for (const auto & enumname : gEnumsToStore) {
      TEnum *en = nullptr;
      const size_t lastSepPos = enumname.find_last_of("::");
      if (lastSepPos != std::string::npos) {
         const std::string nsName = enumname.substr(0, lastSepPos - 1);
         TClass *tclassInstance = TClass::GetClass(nsName.c_str());
         if (!tclassInstance) {
            Error("CloseStreamerInfoROOTFile", "Cannot find TClass instance for namespace %s.", nsName.c_str());
            return false;
         }
         auto enumListPtr = tclassInstance->GetListOfEnums();
         if (!enumListPtr) {
            Error("CloseStreamerInfoROOTFile", "TClass instance for namespace %s does not have any enum associated. This is an inconsistency.", nsName.c_str());
            return false;
         }
         const std::string unqualifiedEnumName = enumname.substr(lastSepPos + 1);
         en = (TEnum *)enumListPtr->FindObject(unqualifiedEnumName.c_str());
         if (en) en->SetTitle(nsName.c_str());
      } else {
         en = (TEnum *)gROOT->GetListOfEnums()->FindObject(enumname.c_str());
         if (en) en->SetTitle("");
      }
      if (!en) {
         Error("CloseStreamerInfoROOTFile", "Cannot find enum %s.", enumname.c_str());
         return false;
      }
      en->Property(); // Force initialization of the bits and property fields.
      enums.AddLast(en);
   }

   if (dictFile.IsZombie())
      return false;
// Instead of plugins:
   protoClasses.Write("__ProtoClasses", TObject::kSingleKey);
   protoClasses.Delete();
   typedefs.Write("__Typedefs", TObject::kSingleKey);
   enums.Write("__Enums", TObject::kSingleKey);

   return true;
}
