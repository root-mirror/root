// @(#)root/memstat:$Id$
// Author: Anar Manafov (A.Manafov@gsi.de) 2008-03-02

/*************************************************************************
* Copyright (C) 1995-2010, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/
#ifndef ROOT_TMemStatMng
#define ROOT_TMemStatMng

// STD
#include <map>
#include <cstring>
// ROOT
#include "TObject.h"
#include "TTimeStamp.h"
// Memstat
#include "TMemStatHook.h"
#include "TMemStatDef.h"

class TTree;
class TFile;
class TH1I;
class TObjArray;

namespace Memstat {

   class TMemStatFAddrContainer {
      typedef std::map<ULong_t, Int_t> Container_t;
      typedef Container_t::iterator pos_type;
      typedef Container_t::value_type value_type;

   public:
      bool add(ULong_t addr, Int_t idx) {
         std::pair<pos_type, bool> ret = fContainer.insert(value_type(addr, idx));
         return (ret.second);
      }

      Int_t find(ULong_t addr) {
         pos_type iter = fContainer.find(addr);
         if(fContainer.end() == iter)
            return -1;

         return iter->second;
      }

   private:
      Container_t fContainer;
   };

   const UShort_t g_digestSize = 16;
   struct SCustomDigest {
      SCustomDigest() {
         memset(fValue, 0, g_digestSize);
      }
      SCustomDigest(UChar_t _val[g_digestSize]) {
         memcpy(fValue, _val, g_digestSize);
      }

      UChar_t fValue[g_digestSize];
   };
   inline bool operator< (const SCustomDigest &a, const SCustomDigest &b)
   {
      for(int i = 0; i < g_digestSize; ++i) {
         if(a.fValue[i] != b.fValue[i])
            return (a.fValue[i] < b.fValue[i]);
      }
      return false;
   }


   class TMemStatMng: public TObject {
      typedef std::map<SCustomDigest, Int_t> CRCSet_t;

   private:
      TMemStatMng();
      virtual ~TMemStatMng();

   public:
      void Enable();                       //enable memory statistic
      void Disable();                      //Disable memory statistic
      static TMemStatMng* GetInstance();   //get instance of class - ONLY ONE INSTANCE
      static void Close();                 //close MemStatManager
      void SetBufferSize(Int_t buffersize);
      void SetMaxCalls(Int_t maxcalls);

   public:
      //stack data members
      void SetUseGNUBuiltinBacktrace(Bool_t newVal) {
         fUseGNUBuiltinBacktrace = newVal;
      }

   protected:
#if !defined(__APPLE__)
      TMemStatHook::MallocHookFunc_t fPreviousMallocHook;    //!old malloc function
      TMemStatHook::FreeHookFunc_t fPreviousFreeHook;        //!old free function
#endif
      void Init();
      void AddPointer(void *ptr, Int_t size);    //add pointer to the table
      void FillTree();
      static void *AllocHook(size_t size, const void* /*caller*/);
      static void FreeHook(void* ptr, const void* /*caller*/);
      static void MacAllocHook(void *ptr, size_t size);
      static void MacFreeHook(void *ptr);
      Int_t generateBTID(UChar_t *CRCdigest, Int_t stackEntries,
                         void **stackPointers);


      //  memory information
      TFile* fDumpFile;               //!file to dump current information
      TTree *fDumpTree;               //!tree to dump information
      static TMemStatMng *fgInstance; // pointer to instance
      static void *fgStackTop;        // stack top pointer

      Bool_t fUseGNUBuiltinBacktrace;
      TTimeStamp fTimeStamp;
      Double_t  fBeginTime;   //time when monitoring starts
      ULong64_t fPos;         //position in memory where alloc/free happens
      Int_t     fTimems;      //10000*(current time - begin time)
      Int_t     fNBytes;      //number of bytes allocated/freed
      Int_t     fBtID;        //back trace identifier
      Int_t     fMaxCalls;    //max number of malloc/frees to register in the output Tree
      Int_t     fBufferSize;  //max number of malloc/free to keep in the buffer
      Int_t     fBufN;        //current number of alloc or free in the buffer
      ULong64_t *fBufPos;     //position in memory where alloc/free happens
      Int_t     *fBufTimems;  //10000*(current time - begin time)
      Int_t     *fBufNBytes;  //number of bytes allocated/freed
      Int_t     *fBufBtID;    //back trace identifier
      Int_t     *fIndex;      //array to sort fBufPos
      Bool_t    *fMustWrite;  //flag to write or not the entry

   private:
      TMemStatFAddrContainer fFAddrs;
      TObjArray *fFAddrsList;
      TH1I *fHbtids;
      CRCSet_t fBTChecksums;
      Int_t fBTCount;
      // for Debug. A counter of all (de)allacations.
      UInt_t  fBTIDCount;
      TNamed *fSysInfo;

      ClassDef(TMemStatMng, 0)   // a manager of memstat sessions.
   };

}

#endif
