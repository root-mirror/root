// @(#)root/cont:$Id$
// Author: Fons Rademakers   11/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClassTable
#define ROOT_TClassTable


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClassTable                                                          //
//                                                                      //
// This class registers for all classes their name, id and dictionary   //
// function in a hash table.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include <string>

class TProtoClass;

namespace ROOT {
   class TClassAlt;
   class TClassRec;
   class TMapTypeToClassRec;
}

class TClassTable : public TObject {

friend  void ROOT::ResetClassVersion(TClass*, const char*, Short_t);
friend  class TROOT;

private:
   typedef ROOT::TMapTypeToClassRec IdMap_t;

   static ROOT::TClassAlt **fgAlternate;
   static ROOT::TClassRec **fgTable;
   static ROOT::TClassRec **fgSortedTable;
   static IdMap_t     *fgIdMap;
   static UInt_t       fgSize;
   static UInt_t       fgTally;
   static Bool_t       fgSorted;
   static UInt_t       fgCursor;

   TClassTable();

   static ROOT::TClassRec   *FindElementImpl(const char *cname, Bool_t insert);
   static ROOT::TClassRec   *FindElement(const char *cname, Bool_t insert=kFALSE);
   static void         SortTable();

   static Bool_t CheckClassTableInit();

public:
   // bits that can be set in pragmabits
   enum {
      kNoStreamer = 0x01, kNoInputOperator = 0x02, kAutoStreamer = 0x04,
      kHasVersion = 0x08, kHasCustomStreamerMember = 0x10
   };

   ~TClassTable();

   static void          Add(const char *cname, Version_t id,
                            const std::type_info &info, DictFuncPtr_t dict,
                            Int_t pragmabits);
   static void          Add(TProtoClass *protoClass);
   static void          AddAlternate(const char *normname, const char *alternate);
   static char         *At(UInt_t index);
   int                  Classes();
   static Bool_t        Check(const char *cname, std::string &normname);
   static Version_t     GetID(const char *cname);
   static Int_t         GetPragmaBits(const char *name);
   static DictFuncPtr_t GetDict(const char *cname);
   static DictFuncPtr_t GetDict(const std::type_info& info);
   static DictFuncPtr_t GetDictNorm(const char *cname);
   static TProtoClass  *GetProto(const char *cname);
   static TProtoClass  *GetProtoNorm(const char *cname);
   static void          Init();
   static char         *Next();
   void                 Print(Option_t *option="") const;
   static void          PrintTable();
   static void          Remove(const char *cname);
   static void          Terminate();

   ClassDef(TClassTable,0)  //Table of known classes
};

R__EXTERN TClassTable *gClassTable;

namespace ROOT {
   extern void AddClass(const char *cname, Version_t id, DictFuncPtr_t dict,
                        Int_t pragmabits);
   extern void RemoveClass(const char *cname);
}

#endif
