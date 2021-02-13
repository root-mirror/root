// @(#)root/test:$Id$
// Author: Nikolay Root   05/07/98

#include <cstdlib>
#include <iostream>
#include "snprintf.h"
#include "TCollection.h"
#include "TSortedList.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TOrdCollection.h"
#include "THashList.h"
#include "THashTable.h"
#include "TBtree.h"

#include "TStopwatch.h"
#include "TRandom.h"
//
// This program benchmarks access time to objects by name or index
// for TObjArray,TOrdCollection,TList,TSortedList,THashList,TBtree,
// TClonesArray and THashTable collections.
//
// Usage: tcollbm -h                               - to print a usage info
//        tcollbm [-n|-i|-m] [nobjects] [ntimes]   - to run tests
//
// switches:
//       -n            - benchmark access by name  (default)
//       -i            - benchmark access by index
//       -m            - benchmark of objects allocation
//
// parameters:
//       nobjects      - number of objects to be inserted into collections
//       ntimes        - number of random lookups in the collection
//
//  default values  ( adjusted in 'main' according 'moda' value )

int nobjects = 100;       // Number of objects to be inserted.
int ntimes   = 100000;    // Number of random lookups.
int moda     = 0;         // Default - access by name

//_____________________________________________________________

class HitNames {          // Utility class to handle the table of names
  Int_t  fSize;           // Array size
  char** fNames;          // Array for objects names
public:
  void InitNames();
  void DeleteNames();
  void Resize(Int_t s) {
    if(s <= fSize) return;
    DeleteNames(); fSize=s; InitNames();
  }
  HitNames() : fSize(0), fNames(0) {};
  HitNames(Int_t s) : fSize(s), fNames(0) { InitNames(); };
  virtual ~HitNames()                     { DeleteNames(); };
  const char *operator[](Int_t i)       { return fNames[i]; };
};

void HitNames::InitNames() {        // Create and fill names
  //  std::cout << "HitNames::InitNames() are called" << std::endl;
  fNames = new char* [fSize];
  for (int i=0;i<fSize;i++) {
    fNames[i] = new char [20];
    snprintf(fNames[i],20,"W%d",i+1);
  }
}
void HitNames::DeleteNames() {      // Clean up
  if(!fNames) return;
  //  std::cout << "HitNames::DeleteNames() are called" << std::endl;
  for (int i=0;i<fSize;i++) {
    delete[] fNames[i];
  }
  delete[] fNames;
  fNames = 0;
}

HitNames names;   // We needs only one static object of this class

//_______________________________________________________________
//
// Just for fun, I introduce the class 'Tester' to collect
// a 'TNamed' objects. I use a objects of this class to perform a real
// tests.

class Tester : public TObject {
  enum EClass { Clones, Array, BTree, Other };
  EClass        fWhat;
  Int_t         fNobj;       // Number of objects to be stored in collection
  Int_t         fNtimes;     // How many times to perform a test
  Int_t         fModa;       // What we need to test
  TCollection  *fColl;       // Collection under test
public:

  void Fill();               // Fill the collection
  Double_t TestAllocation(); // Memory allocation test
  Double_t TestByName();     // benchmark by name
  Double_t TestByIndex();    // benchmark by index
  Double_t DoTest();         // Tests multiplexsor

  void        CleanUp()    { fColl->Delete(); }
  void        Dump() const { fColl->Dump(); }

  virtual const char* GetName() const
  { return fColl->ClassName(); }

  Tester() :
    fWhat(Other),fNobj(10),fNtimes(100000),fModa(0),fColl(0) {}
  Tester(Int_t no,Int_t nt,Int_t m,TCollection *p) :
    fWhat(Other),fNobj(no),fNtimes(nt),fModa(m),fColl(p) {
      if(!strcmp(GetName(),"TClonesArray")) fWhat = Clones;
      if(!strcmp(GetName(),"TObjArray"))    fWhat = Array;
      if(!strcmp(GetName(),"TBtree"))       fWhat = BTree;
  }
  virtual ~Tester() { if(fColl) delete fColl; }
};

void Tester::Fill() {
  if (fWhat == Clones) {
    TClonesArray &base = *((TClonesArray*)fColl);
    for (Int_t i=0;i<fNobj;i++) new (base[i]) TNamed(names[i],GetName());
  } else {
    for (Int_t i=0;i<fNobj;i++) fColl->Add(new TNamed(names[i],GetName()));
  }
}

Double_t Tester::TestAllocation() {   // benchmark memory allocation
  TStopwatch timer;
  timer.Start();
  for (Int_t i=0;i<fNtimes;i++) {
    Fill();
    CleanUp();
  }
  timer.Stop();
  return timer.CpuTime();
}

Double_t Tester::TestByIndex() {      // benchmark access by index
  Fill();
  TStopwatch timer;
  Int_t i;
  if (fWhat == Clones) {
    TClonesArray *o = (TClonesArray*)fColl;
    timer.Start();
    for (Int_t j=0;j<fNtimes;j++) {
      i=Int_t(fNobj*gRandom->Rndm());
      if(!((*o)[i])) Printf("Object %d not found !!!",i);
    }
    timer.Stop();
  } else if(fWhat == Array) {
    TObjArray    *o = (TObjArray*)fColl;
    timer.Start();
    for (Int_t j=0;j<fNtimes;j++) {
      i=Int_t(fNobj*gRandom->Rndm());
      if(!((*o)[i])) Printf("Object %d not found !!!",i);
    }
    timer.Stop();
  } else if(fWhat == BTree) {
    TBtree       *o = (TBtree*)fColl;
    timer.Start();
    for (Int_t j=0;j<fNtimes;j++) {
      i=Int_t(fNobj*gRandom->Rndm());
      if(!((*o)[i])) Printf("Object %d not found !!!",i);
    }
    timer.Stop();
  } else {
     std::cout << "Class " << GetName()
     << " doesn't support access by index" << std::endl;
  }
  CleanUp();
  return timer.CpuTime();
};

Double_t Tester::TestByName() {           // benchmark access by name
  Fill();
  TStopwatch timer;
  Int_t i;
  timer.Start();
  for (Int_t j=0;j<fNtimes;j++) {
    i=Int_t(fNobj*gRandom->Rndm());
    if(!((*fColl)(names[i]))) Printf("Object %5s not found !!!",names[i]);
  }
  timer.Stop();
  CleanUp();
  return timer.CpuTime();
};

Double_t Tester::DoTest() {
  // Return the average time in msec.
  Double_t v;
  if(fModa==2) {
    printf("Memory allocation test for %-20s", GetName());
    v=TestAllocation();
  } else if(fModa==1) {
    printf("Random lookups by Index for %-20s", GetName());
    v=TestByIndex();
  } else {
    printf("Random lookups by Name for %-20s", GetName());
    v=TestByName();
  }
  Printf(" CpuTime=%7.2f seconds", v);
  return 1000*v/fNtimes;
};
////////////////////////////////////////////////////////////////////////////////

const int ntests = 8;     // Number of classes to be tested
Double_t deltas[ntests];  // benchmark results

int main(int argc,char **argv)
{
  // Initialize the ROOT framework
  if(argc == 2 && !strcmp(argv[1],"-h")) {
    Printf("Usage: tcollbm [-n|-i] [nobjects] [ntimes]");
    Printf("  -n        - benchmark access by name");
    Printf("  -i        - benchmark access by index");
    Printf("  -m        - benchmark memory allocation");
    Printf("  nobjects  - number of objects to be inserted into collections");
    Printf("  ntimes    - number of random lookups in the collection");
    return 1;
  }
  if(argc > 1 && !strcmp(argv[1],"-n")) { moda = 0; argc--; argv++; };
  if(argc > 1 && !strcmp(argv[1],"-i")) { moda = 1; argc--; argv++; };
  if(argc > 1 && !strcmp(argv[1],"-m")) { moda = 2; argc--; argv++; };
  //
  // Set defaults values for selected test
  //
  if(moda == 0) { nobjects = 100;  ntimes = 10000; }
  if(moda == 1) { nobjects = 100;  ntimes = 1000000; }
  if(moda == 2) { nobjects = 1000; ntimes = 100; }

  Int_t no = nobjects;
  Int_t nt = ntimes;
  if(argc > 1) { no = atoi(argv[1]); }
  if(argc > 2) { nt = atoi(argv[2]); }
  if(no > 1) {
    nobjects = no;
  } else {
    nobjects = 2;
    Printf("Reset nobjects to %d",nobjects);
  }
  if(nt > 99) {
    ntimes   = nt;
  } else {
    ntimes = 100;
    Printf("Reset ntimes   to %d",ntimes);
  }

  Printf("Nobjects = %d , Ntimes = %d",nobjects,ntimes);

  if(ntimes > 5000000 || (moda==2 && nobjects*ntimes > 500000)) {
    Printf("This test takes some time. Be partient ...");
  }

  names.Resize(nobjects);
  Double_t smin = 10000000;
  int j=-1,idx=-1;
  TObjArray array(ntests);                         // Array of collections

  array.Add(new Tester(nobjects,ntimes,moda,
                       new TObjArray()));         // Add TObjArray
  j++;
  deltas[j] = ((Tester*)(array[j]))->DoTest();
  if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

  array.Add(new Tester(nobjects,ntimes,moda,
                       new TClonesArray("TNamed",nobjects)));  // Add TClonesArray
  j++;
  deltas[j] = ((Tester*)(array[j]))->DoTest();
  if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

  array.Add(new Tester(nobjects,ntimes,moda,
                       new TBtree()));                       // Add TBTree
  j++;
  deltas[j] = ((Tester*)(array[j]))->DoTest();
  if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

  if(moda != 1) {   // Skip classes without operator[] from index test

    array.Add(new Tester(nobjects,ntimes,moda,
                         new TOrdCollection()));            // Add TOrdCollection
    j++;
    deltas[j] = ((Tester*)(array[j]))->DoTest();
    if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

    array.Add(new Tester(nobjects,ntimes,moda,
                         new TList()));                     // Add TList
    j++;
    deltas[j] = ((Tester*)(array[j]))->DoTest();
    if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

    array.Add(new Tester(nobjects,ntimes,moda,
                         new TSortedList()));               // Add TSortedList
    j++;
    deltas[j] = ((Tester*)(array[j]))->DoTest();
    if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

    array.Add(new Tester(nobjects,ntimes,moda,
                         new THashList()));                 // Add THashList
    j++;
    deltas[j] = ((Tester*)(array[j]))->DoTest();
    if(deltas[j] < smin) { idx=j; smin=deltas[j]; }

    array.Add(new Tester(nobjects,ntimes,moda,
                         new THashTable()));                // Add THashTable
    j++;
    deltas[j] = ((Tester*)(array[j]))->DoTest();
    if(deltas[j] < smin) { idx=j; smin=deltas[j]; }
  }

  //          Results info ...
  j++;
  if(idx<0) {
    Printf("Can not find the winner. Sorry...\n");
    array.Delete();
    return 1;
  }
  if(smin==0) {
    Printf(
    "Input parameters nobjects=%d and ntimes=%d too small for Your computer.",
    nobjects,ntimes);
    Printf("Please, increase its and try again...");
    array.Delete();
    return 1;
  }

  Printf("\n\t\tBenchmark results\n");
  Printf("\tAbsolute winner - %s",array[idx]->GetName());

  Printf("\tAverage (msec)\tRatio\t\tClassName");

  for (int i=0;i<j;i++) {
     Printf("\t%f\t%f\t%s",
            deltas[i],deltas[i]/smin,array[i]->GetName());
  }
  array.Delete();
  return 0;
}
