// @(#)root/thread:$Id$
// // Author: Enric Tejedor Saavedra   03/12/15
//
/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TImplicitMT                                                          //
//                                                                      //
// This file implements the methods to enable, disable and check the    //
// status of the global implicit multi-threading in ROOT.               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include "TThread.h"
#include "ROOT/TPoolManager.hxx"
#include <atomic>
#include <string>

static std::shared_ptr<ROOT::Internal::TPoolManager> &R__GetPoolManagerMT()
{
   static std::shared_ptr<ROOT::Internal::TPoolManager> schedMT;
   return schedMT;
}

extern "C" UInt_t ROOT_MT_GetThreadPoolSize()
{
   return ROOT::Internal::TPoolManager::GetPoolSize();
};

static bool &GetImplicitMTFlag()
{
   static bool enabled = false;
   return enabled;
}

static std::string &GetImplicitMTBackend()
{
   static std::string backend = "TBB";
   return backend;
}

static std::atomic_int &GetParBranchProcessingCount()
{
   static std::atomic_int count(0);
   return count;
}

extern "C" void ROOT_TImplicitMT_EnableImplicitMT(UInt_t numthreads, std::string backend)
{
   if (backend == "TBB") {
      GetImplicitMTBackend() = backend;
#ifdef R__HAS_HPX
   } else if (backend == "HPX") {
      GetImplicitMTBackend() = backend;
#endif
   } else {
      ::Warning("ROOT_TImplicitMT_EnableImplicitMT", "Backend %s not supported, proceeding with TBB as backend.", backend.c_str());
      GetImplicitMTBackend() = "TBB";
   }

   if (!GetImplicitMTFlag()) {
      if (ROOT::Internal::TPoolManager::GetPoolSize() == 0) {
         TThread::Initialize();
      }
      R__GetPoolManagerMT() = ROOT::Internal::GetPoolManager(numthreads);
      GetImplicitMTFlag() = true;
   } else {
      ::Warning("ROOT_TImplicitMT_EnableImplicitMT", "Implicit multi-threading is already enabled");
   }
};

extern "C" void ROOT_TImplicitMT_DisableImplicitMT()
{
   if (GetImplicitMTFlag()) {
      GetImplicitMTFlag() = false;
      R__GetPoolManagerMT().reset();
   } else {
      ::Warning("ROOT_TImplicitMT_DisableImplicitMT", "Implicit multi-threading is already disabled");
   }
};

extern "C" std::string ROOT_TImplicitMT_ImplicitMTBackend()
{
   return GetImplicitMTBackend();
};

extern "C" void ROOT_TImplicitMT_EnableParBranchProcessing()
{
   ++GetParBranchProcessingCount();
};

extern "C" void ROOT_TImplicitMT_DisableParBranchProcessing()
{
   --GetParBranchProcessingCount();
};

extern "C" bool ROOT_TImplicitMT_IsParBranchProcessingEnabled()
{
   return GetParBranchProcessingCount() > 0;
};
