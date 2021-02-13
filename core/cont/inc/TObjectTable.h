// @(#)root/cont:$Id$
// Author: Fons Rademakers   11/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TObjectTable
#define ROOT_TObjectTable


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TObjectTable                                                         //
//                                                                      //
// This class registers all instances of TObject and its derived        //
// classes in a hash table. The Add() and Remove() members are called   //
// from the TObject ctor and dtor, repectively. Using the Print()       //
// member one can see all currently active objects in the system.       //
// Using the runtime flag: Root.ObjectStat one can toggle this feature  //
// on or off.                                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TObject.h"

class TClass;


class TObjectTable : public TObject {

private:
   TObject  **fTable;
   Int_t      fSize;
   Int_t      fTally;

   Bool_t     HighWaterMark();
   void       Expand(Int_t newsize);
   Int_t      FindElement(TObject *obj);
   void       FixCollisions(Int_t index);

private:
   TObjectTable(const TObjectTable&) = delete;
   TObjectTable& operator=(const TObjectTable&) = delete;

public:
   TObjectTable(Int_t tableSize = 100);
   ~TObjectTable();

   void      Add(TObject *obj);
   void     *CheckPtrAndWarn(const char *msg, void *vp);
   void      Delete(Option_t *opt = "");
   Int_t     GetSize() const { return fSize; }
   Int_t     Instances() const { return fTally; }
   void      InstanceStatistics() const;
   void      Print(Option_t *option="") const;
   Bool_t    PtrIsValid(TObject *obj);
   void      Remove(TObject *obj);
   void      RemoveQuietly(TObject *obj);
   void      Statistics() { Print(); }
   void      Terminate();
   void      UpdateInstCount() const;

   static void AddObj(TObject *obj);

   ClassDef(TObjectTable,0)  //Table of active objects
};


inline Bool_t TObjectTable::HighWaterMark()
   { return (Bool_t) (fTally >= ((3*fSize)/4)); }

inline Bool_t TObjectTable::PtrIsValid(TObject *op)
   { return fTable[FindElement(op)] != 0; }


R__EXTERN TObjectTable *gObjectTable;

#endif
