// @(#)root/cont:$Id$
// Author: Fons Rademakers   11/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TObjectTable
\ingroup Containers
This class registers all instances of TObject and its derived
classes in a hash table. The Add() and Remove() members are called
from the TObject ctor and dtor, respectively. Using the Print()
member one can see all currently active objects in the system.
Using the resource (in .rootrc): Root.ObjectStat one can toggle this
feature on or off.

Using the compile option R__NOSTATS one can de-active this feature
for the entire system (for maximum performance in highly time
critical applications).

The following output has been produced in a ROOT interactive session
via the command gObjectTable->Print()
~~~ {.cpp}
  class                     cnt    on heap     size    total size    heap size
  ============================================================================
  TKey                        4          4       72           288          288
  TClass                     84         84       80          6720         6720
  TDataMember               276        276       24          6624         6624
  TObject                    11         11       12           132          132
  TMethod                  1974       1974       64        126336       126336
  TDataType                  34         34       56          1904         1904
  TList                    2328       2328       36         83808        83808
  TH1F                        1          1      448           448          448
  TText                    2688       2688       56        150528       150528
  TGaxis                      1          0      120           120            0
  TAxis                       6          3       88           528          264
  TBox                       57         57       52          2964         2964
  TLine                     118        118       40          4720         4720
  TWbox                       1          1       56            56           56
  TArrow                      1          1       64            64           64
  TPaveText                  59         59      124          7316         7316
  TPave                       1          1       92            92           92
  TFile                       1          1      136           136          136
  TCanvas                     3          3      444          1332         1332
  TPad                        1          1      312           312          312
  TContextMenu                3          3       48           144          144
  TMethodArg               2166       2166       44         95304        95304
  TPaveLabel                  1          1      120           120          120
  THtml                       1          1       32            32           32
  TROOT                       1          0      208           208            0
  TApplication                1          1       28            28           28
  TFileHandler                1          1       20            20           20
  TColor                    163        163       40          6520         6520
  TStyle                      1          1      364           364          364
  TRealData                 117        117       28          3276         3276
  TBaseClass                 88         88       36          3168         3168
  THashList                   5          5       40           200          200
  THashTable                  5          5       36           180          180
  TGeometry                   1          1       64            64           64
  TLink                       7          7       60           420          420
  TPostScript                 1          1      764           764          764
  TMinuit                     1          1      792           792          792
  TStopwatch                  1          0       56            56            0
  TRootGuiFactory             1          1       28            28           28
  TGX11                       1          1      172           172          172
  TUnixSystem                 1          1      252           252          252
  TSignalHandler              1          1       20            20           20
  TOrdCollection              3          3       40           120          120
  TEnv                        1          1       24            24           24
  TCling                      1          1      208           208          208
  TBenchmark                  1          1       52            52           52
  TClassTable                 1          1       12            12           12
  TObjectTable                1          1       12            12           12
  ----------------------------------------------------------------------------
  Total:                  10225      10219     5976        506988       506340
  ============================================================================
~~~
*/

#include "TObjectTable.h"
#include "TList.h"
#include "TROOT.h"
#include "TClass.h"
#include "TError.h"


TObjectTable *gObjectTable;


ClassImp(TObjectTable);

////////////////////////////////////////////////////////////////////////////////
/// Create an object table.

TObjectTable::TObjectTable(Int_t tableSize)
{
   fSize  = (Int_t)TMath::NextPrime(tableSize);
   fTable = new TObject* [fSize];
   memset(fTable, 0, fSize*sizeof(TObject*));
   fTally = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete TObjectTable.

TObjectTable::~TObjectTable()
{
   delete [] fTable; fTable = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the object table.
/// If option ="all" prints the list of all objects with the format
/// object number, pointer, class name, object name

void TObjectTable::Print(Option_t *option) const
{
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("all")) {
      TObject *obj;
      int i, num = 0;
      Printf("\nList of all objects");
      Printf("object   address            class                    name");
      Printf("================================================================================");
      for (i = 0; i < fSize; i++) {
         if (!fTable[i]) continue;
         num++;
         obj = fTable[i];
         printf("%-8d 0x%-16lx %-24s %s\n", num, (Long_t)obj, obj->ClassName(),
                obj->GetName());
      }
      Printf("================================================================================\n");
   }

   //print the number of instances per class
   InstanceStatistics();
}

////////////////////////////////////////////////////////////////////////////////
/// Add an object to the object table.

void TObjectTable::Add(TObject *op)
{
   if (!op) {
      Error("Add", "op is 0");
      return;
   }
   if (!fTable)
      return;

   Int_t slot = FindElement(op);
   if (fTable[slot] == 0) {
      fTable[slot] = op;
      fTally++;
      if (HighWaterMark())
         Expand(2 * fSize);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add an object to the global object table gObjectTable. If the global
/// table does not exist create it first. This member function may only
/// be used by TObject::TObject. Use Add() to add objects to any other
/// TObjectTable object. This is a static function.

void TObjectTable::AddObj(TObject *op)
{
   static Bool_t olock = kFALSE;

   if (!op) {
      ::Error("TObjectTable::AddObj", "op is 0");
      return;
   }
   if (olock)
      return;

   if (!gObjectTable) {
      olock = kTRUE;
      gObjectTable = new TObjectTable(10000);
      olock = kFALSE;
      gObjectTable->Add(gObjectTable);
   }

   gObjectTable->Add(op);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all objects stored in the TObjectTable.

void TObjectTable::Delete(Option_t *)
{
   for (int i = 0; i < fSize; i++) {
      if (fTable[i]) {
         delete fTable[i];
         fTable[i] = 0;
      }
   }
   fTally = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an object from the object table.

void TObjectTable::Remove(TObject *op)
{
   if (op == 0) {
      Error("Remove", "remove 0 from TObjectTable");
      return;
   }

   if (!fTable)
      return;

   Int_t i = FindElement(op);
   if (fTable[i] == 0) {
      Warning("Remove", "0x%lx not found at %d", (Long_t)op, i);
      for (int j = 0; j < fSize; j++) {
         if (fTable[j] == op) {
            Error("Remove", "0x%lx found at %d !!!", (Long_t)op, j);
            i = j;
         }
      }
   }

   if (fTable[i]) {
      fTable[i] = 0;
      FixCollisions(i);
      fTally--;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an object from the object table. If op is 0 or not in the table
/// don't complain. Currently only used by the TClonesArray dtor. Should not
/// be used anywhere else, except in places where "special" allocation and
/// de-allocation tricks are performed.

void TObjectTable::RemoveQuietly(TObject *op)
{
   if (op == 0) return;

   if (!fTable)
      return;

   Int_t i = FindElement(op);
   if (fTable[i] == 0)
      for (int j = 0; j < fSize; j++)
         if (fTable[j] == op)
            i = j;

   fTable[i] = 0;
   FixCollisions(i);
   fTally--;
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the object table (this static class function calls the dtor).

void TObjectTable::Terminate()
{
   InstanceStatistics();
   delete [] fTable; fTable = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in the object table. Returns the slot where to put
/// the object. To test if the object is actually already in the table
/// use PtrIsValid().

Int_t TObjectTable::FindElement(TObject *op)
{
   Int_t    slot, n;
   TObject *slotOp;

   if (!fTable)
      return 0;

   //slot = Int_t(((ULong_t) op >> 2) % fSize);
   slot = Int_t(TString::Hash(&op, sizeof(TObject*)) % fSize);
   for (n = 0; n < fSize; n++) {
      if ((slotOp = fTable[slot]) == 0)
         break;
      if (op == slotOp)
         break;
      if (++slot == fSize)
         slot = 0;
   }
   return slot;
}

////////////////////////////////////////////////////////////////////////////////
/// Rehash the object table in case an object has been removed.

void TObjectTable::FixCollisions(Int_t index)
{
   Int_t oldIndex, nextIndex;
   TObject *nextObject;

   for (oldIndex = index+1; ;oldIndex++) {
      if (oldIndex >= fSize)
         oldIndex = 0;
      nextObject = fTable[oldIndex];
      if (nextObject == 0)
         break;
      nextIndex = FindElement(nextObject);
      if (nextIndex != oldIndex) {
         fTable[nextIndex] = nextObject;
         fTable[oldIndex] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand the object table.

void TObjectTable::Expand(Int_t newSize)
{
   TObject **oldTable = fTable, *op;
   int oldsize = fSize;
   newSize = (Int_t)TMath::NextPrime(newSize);
   fTable  = new TObject* [newSize];
   memset(fTable, 0, newSize*sizeof(TObject*));
   fSize   = newSize;
   fTally  = 0;
   for (int i = 0; i < oldsize; i++)
      if ((op = oldTable[i]))
         Add(op);
   delete [] oldTable;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the object table.

void TObjectTable::InstanceStatistics() const
{
   int n, h, s, ncum = 0, hcum = 0, scum = 0, tcum = 0, thcum = 0;

   if (fTally == 0 || !fTable)
      return;

   UpdateInstCount();

   Printf("\nObject statistics");
   Printf("class                         cnt    on heap     size    total size    heap size");
   Printf("================================================================================");
   TIter next(gROOT->GetListOfClasses());
   TClass *cl;
   while ((cl = (TClass*) next())) {
      n = cl->GetInstanceCount();
      h = cl->GetHeapInstanceCount();
      s = cl->Size();
      if (n > 0) {
         Printf("%-24s %8d%11d%9d%14d%13d", cl->GetName(), n, h, s, n*s, h*s);
         ncum  += n;
         hcum  += h;
         scum  += s;
         tcum  += n*s;
         thcum += h*s;
      }
   }
   Printf("--------------------------------------------------------------------------------");
   Printf("Total:                   %8d%11d%9d%14d%13d", ncum, hcum, scum, tcum, thcum);
   Printf("================================================================================\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Histogram all objects according to their classes.

void TObjectTable::UpdateInstCount() const
{
   TObject *op;

   if (!fTable || !TROOT::Initialized())
      return;

   gROOT->GetListOfClasses()->R__FOR_EACH(TClass,ResetInstanceCount)();

   for (int i = 0; i < fSize; i++)
      if ((op = fTable[i])) {                // attention: no ==
         if (op->TestBit(TObject::kNotDeleted))
            op->IsA()->AddInstance(op->IsOnHeap());
         else
            Error("UpdateInstCount", "oops 0x%lx\n", (Long_t)op);
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Issue a warning in case an object still appears in the table
/// while it should not.

void *TObjectTable::CheckPtrAndWarn(const char *msg, void *vp)
{
   if (fTable && vp && fTable[FindElement((TObject*)vp)]) {
      Remove((TObject*)vp);
      Warning("CheckPtrAndWarn", "%s (0x%lx)\n", msg, (Long_t)vp);
   }
   return vp;
}
