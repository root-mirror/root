// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**
  \defgroup tree Tree Library

  In order to store columnar datasets, ROOT provides the TTree, TChain,
  TNtuple and TNtupleD classes.
  The TTree class represents a columnar dataset. Any C++ type can be stored in the
  columns. The TTree has allowed to store about **1 EB** of data coming from the LHC alone:
  it is demonstrated to scale and it's battle tested. It has been optimized during the years
  to reduce dataset sizes on disk and to deliver excellent runtime performance.
  It allows to access only part of the columns of the datasets, too.
  The TNtuple and TNtupleD classes are specialisations of the TTree class which can
  only hold single precision and double precision floating-point numbers respectively;
  The TChain is a collection of TTrees, which can be located also in different files.

*/

/** \class TTree
\ingroup tree

A TTree represents a columnar dataset. Any C++ type can be stored in its columns.

A TTree, often called in jargon *tree*, consists of a list of independent columns or *branches*,
represented by the TBranch class.
Behind each branch, buffers are allocated automatically by ROOT.
Such buffers are automatically written to disk or kept in memory until the size stored in the
attribute fMaxVirtualSize is reached.
Variables of one branch are written to the same buffer. A branch buffer is
automatically compressed if the file compression attribute is set (default).
Branches may be written to different files (see TBranch::SetFile).

The ROOT user can decide to make one single branch and serialize one object into
one single I/O buffer or to make several branches.
Making several branches is particularly interesting in the data analysis phase,
when it is desirable to have a high reading rate and not all columns are equally interesting

## Table of contents:
- [Creating a TTree](#creatingattree)
- [Add a Column of Fundamental Types and Arrays thereof](#addcolumnoffundamentaltypes)
- [Add a Column of a STL Collection instances](#addingacolumnofstl)
- [Add a column holding an object](#addingacolumnofobjs)
- [Add a column holding a TObjectArray](#addingacolumnofobjs)
- [Fill the tree](#fillthetree)
- [Add a column to an already existing Tree](#addcoltoexistingtree)
- [An Example](#fullexample)

## <a name="creatingattree"></a>Creating a TTree

~~~ {.cpp}
    TTree tree(name, title)
~~~
Creates a Tree with name and title.

Various kinds of branches can be added to a tree:
- Variables representing fundamental types, simple classes/structures or list of variables: for example for C or Fortran
structures.
- Any C++ object or collection, provided by the STL or ROOT.

In the following, the details about the creation of different types of branches are given.

## <a name="addcolumnoffundamentaltypes"></a>Add a column (`branch`) of fundamental types and arrays thereof
This strategy works also for lists of variables, e.g. to describe simple structures.
It is strongly recommended to persistify those as objects rather than lists of leaves.

~~~ {.cpp}
    auto branch = tree.Branch(branchname, address, leaflist, bufsize)
~~~
- address is the address of the first item of a structure
- leaflist is the concatenation of all the variable names and types
  separated by a colon character :
  The variable name and the variable type are separated by a
  slash (/). The variable type must be 1 character. (Characters
  after the first are legal and will be appended to the visible
  name of the leaf, but have no effect.) If no type is given, the
  type of the variable is assumed to be the same as the previous
  variable. If the first variable does not have a type, it is
  assumed of type F by default. The list of currently supported
  types is given below:
   - `C` : a character string terminated by the 0 character
   - `B` : an 8 bit signed integer (`Char_t`)
   - `b` : an 8 bit unsigned integer (`UChar_t`)
   - `S` : a 16 bit signed integer (`Short_t`)
   - `s` : a 16 bit unsigned integer (`UShort_t`)
   - `I` : a 32 bit signed integer (`Int_t`)
   - `i` : a 32 bit unsigned integer (`UInt_t`)
   - `F` : a 32 bit floating point (`Float_t`)
   - `f` : a 24 bit floating point with truncated mantissa (`Float16_t`)
   - `D` : a 64 bit floating point (`Double_t`)
   - `d` : a 24 bit truncated floating point (`Double32_t`)
   - `L` : a 64 bit signed integer (`Long64_t`)
   - `l` : a 64 bit unsigned integer (`ULong64_t`)
   - `G` : a long signed integer, stored as 64 bit (`Long_t`)
   - `g` : a long unsigned integer, stored as 64 bit (`ULong_t`)
   - `O` : [the letter `o`, not a zero] a boolean (`Bool_t`)

  Examples:
   - A int: "myVar/I"
   - A float array with fixed size: "myArrfloat[42]/F"
   - An double array with variable size, held by the `myvar` column: "myArrdouble[myvar]/D"
   - An Double32_t array with variable size, held by the `myvar` column , with values between 0 and 16: "myArr[myvar]/d[0,10]"

- If the address points to a single numerical variable, the leaflist is optional:
~~~ {.cpp}
  int value;
  `tree->Branch(branchname, &value);`
~~~
- If the address points to more than one numerical variable, we strongly recommend
  that the variable be sorted in decreasing order of size.  Any other order will
  result in a non-portable TTree (i.e. you will not be able to read it back on a
  platform with a different padding strategy).
  We recommend to persistify objects rather than composite leaflists.
- In case of the truncated floating point types (Float16_t and Double32_t) you can
  furthermore specify the range in the style [xmin,xmax] or [xmin,xmax,nbits] after
  the type character. For example, for storing a variable size array `myArr` of
  `Double32_t` with values within a range of `[0, 2*pi]` and the size of which is
  stored in a branch called `myArrSize`, the syntax for the `leaflist` string would
  be: `myArr[myArrSize]/d[0,twopi]`. Of course the number of bits could be specified,
  the standard rules of opaque typedefs annotation are valid. For example, if only
  18 bits were sufficient, the syntax would become: `myArr[myArrSize]/d[0,twopi,18]`

## <a name="addingacolumnofstl"></a>Adding a column of STL collection instances (e.g. std::vector, std::list, std::unordered_map)

~~~ {.cpp}
    auto branch = tree.Branch( branchname, STLcollection, buffsize, splitlevel);
~~~
STLcollection is the address of a pointer to std::vector, std::list,
std::deque, std::set or std::multiset containing pointers to objects.
If the splitlevel is a value bigger than 100 (TTree::kSplitCollectionOfPointers)
then the collection will be written in split mode, e.g. if it contains objects of
any types deriving from TTrack this function will sort the objects
based on their type and store them in separate branches in split
mode.

~~~ {.cpp}
    branch->SetAddress(void *address)
~~~
In case of dynamic structures changing with each entry for example, one must
redefine the branch address before filling the branch again.
This is done via the TBranch::SetAddress member function.

## <a name="addingacolumnofobjs">Add a column of objects

~~~ {.cpp}
    MyClass object;
    auto branch = tree.Branch(branchname, &object, bufsize, splitlevel)
~~~
Note: The 2nd parameter must be the address of a valid object.
      The object must not be destroyed (i.e. be deleted) until the TTree
      is deleted or TTree::ResetBranchAddress is called.

- if splitlevel=0, the object is serialized in the branch buffer.
- if splitlevel=1 (default), this branch will automatically be split
  into subbranches, with one subbranch for each data member or object
  of the object itself. In case the object member is a TClonesArray,
  the mechanism described in case C is applied to this array.
- if splitlevel=2 ,this branch will automatically be split
  into subbranches, with one subbranch for each data member or object
  of the object itself. In case the object member is a TClonesArray,
  it is processed as a TObject*, only one branch.

Another available syntax is the following:

~~~ {.cpp}
    auto branch = tree.Branch(branchname, &p_object, bufsize, splitlevel)
    auto branch = tree.Branch(branchname, className, &p_object, bufsize, splitlevel)
~~~
- p_object is a pointer to an object.
- If className is not specified, Branch uses the type of p_object to determine the
  type of the object.
- If className is used to specify explicitly the object type, the className must
  be of a type related to the one pointed to by the pointer.  It should be either
  a parent or derived class.

Note: The pointer whose address is passed to TTree::Branch must not
      be destroyed (i.e. go out of scope) until the TTree is deleted or
      TTree::ResetBranchAddress is called.

Note: The pointer p_object must be initialized before calling TTree::Branch
- Do either:
~~~ {.cpp}
    MyDataClass* p_object = nullptr;
    tree.Branch(branchname, &p_object);
~~~
- Or:
~~~ {.cpp}
    auto p_object = new MyDataClass;
    tree.Branch(branchname, &p_object);
~~~
Whether the pointer is set to zero or not, the ownership of the object
is not taken over by the TTree.  I.e. even though an object will be allocated
by TTree::Branch if the pointer p_object is zero, the object will <b>not</b>
be deleted when the TTree is deleted.

## <a name="addingacolumnoftclonesarray">Add a column of TClonesArray instances

*It is recommended to use STL containers instead of TClonesArrays*.

~~~ {.cpp}
    // clonesarray is the address of a pointer to a TClonesArray.
    auto branch = tree.Branch(branchname,clonesarray, bufsize, splitlevel)
~~~
The TClonesArray is a direct access list of objects of the same class.
For example, if the TClonesArray is an array of TTrack objects,
this function will create one subbranch for each data member of
the object TTrack.

## <a name="fillthetree">Fill the Tree:

A TTree instance is filled with the invocation of the TTree::Fill method:
~~~ {.cpp}
    tree.Fill()
~~~
Upon its invocation, a loop on all defined branches takes place that for each branch invokes
the TBranch::Fill method.

## <a name="addcoltoexistingtree">Add a column to an already existing Tree

You may want to add a branch to an existing tree. For example,
if one variable in the tree was computed with a certain algorithm,
you may want to try another algorithm and compare the results.
One solution is to add a new branch, fill it, and save the tree.
The code below adds a simple branch to an existing tree.
Note the kOverwrite option in the Write method, it overwrites the
existing tree. If it is not specified, two copies of the tree headers
are saved.
~~~ {.cpp}
    void tree3AddBranch() {
        TFile f("tree3.root", "update");

        Float_t new_v;
        auto t3 = f->Get<TTree>("t3");
        auto newBranch = t3->Branch("new_v", &new_v, "new_v/F");

        Long64_t nentries = t3->GetEntries(); // read the number of entries in the t3

        for (Long64_t i = 0; i < nentries; i++) {
            new_v = gRandom->Gaus(0, 1);
            newBranch->Fill();
        }

        t3->Write("", TObject::kOverwrite); // save only the new version of the tree
    }
~~~
It is not always possible to add branches to existing datasets stored in TFiles: for example,
these files might not be writeable, just readable. In addition, modifying in place a TTree
causes a new TTree instance to be written and the previous one to be deleted.
For this reasons, ROOT offers the concept of friends for TTree and TChain:
if is good practice to rely on friend trees rather than adding a branch manually.

## <a name="fullexample">An Example

Begin_Macro
../../../tutorials/tree/tree.C
End_Macro

~~~ {.cpp}
    // A simple example with histograms and a tree
    //
    // This program creates :
    //    - a one dimensional histogram
    //    - a two dimensional histogram
    //    - a profile histogram
    //    - a tree
    //
    // These objects are filled with some random numbers and saved on a file.

    #include "TFile.h"
    #include "TH1.h"
    #include "TH2.h"
    #include "TProfile.h"
    #include "TRandom.h"
    #include "TTree.h"

    //__________________________________________________________________________
    main(int argc, char **argv)
    {
    // Create a new ROOT binary machine independent file.
    // Note that this file may contain any kind of ROOT objects, histograms,trees
    // pictures, graphics objects, detector geometries, tracks, events, etc..
    // This file is now becoming the current directory.
    TFile hfile("htree.root","RECREATE","Demo ROOT file with histograms & trees");

    // Create some histograms and a profile histogram
    TH1F hpx("hpx","This is the px distribution",100,-4,4);
    TH2F hpxpy("hpxpy","py ps px",40,-4,4,40,-4,4);
    TProfile hprof("hprof","Profile of pz versus px",100,-4,4,0,20);

    // Define some simple structures
    typedef struct {Float_t x,y,z;} POINT;
    typedef struct {
       Int_t ntrack,nseg,nvertex;
       UInt_t flag;
       Float_t temperature;
    } EVENTN;
    POINT point;
    EVENTN eventn;

    // Create a ROOT Tree
    TTree tree("T","An example of ROOT tree with a few branches");
    tree.Branch("point",&point,"x:y:z");
    tree.Branch("eventn",&eventn,"ntrack/I:nseg:nvertex:flag/i:temperature/F");
    tree.Branch("hpx","TH1F",&hpx,128000,0);

    Float_t px,py,pz;

    // Here we start a loop on 1000 events
    for ( Int_t i=0; i<1000; i++) {
       gRandom->Rannor(px,py);
       pz = px*px + py*py;
       const auto random = gRandom->::Rndm(1);

       // Fill histograms
       hpx.Fill(px);
       hpxpy.Fill(px,py,1);
       hprof.Fill(px,pz,1);

       // Fill structures
       point.x = 10*(random-1);
       point.y = 5*random;
       point.z = 20*random;
       eventn.ntrack  = Int_t(100*random);
       eventn.nseg    = Int_t(2*eventn.ntrack);
       eventn.nvertex = 1;
       eventn.flag    = Int_t(random+0.5);
       eventn.temperature = 20+random;

       // Fill the tree. For each event, save the 2 structures and 3 objects
       // In this simple example, the objects hpx, hprof and hpxpy are slightly
       // different from event to event. We expect a big compression factor!
       tree->Fill();
    }
    // End of the loop

    tree.Print();

    // Save all objects in this file
    hfile.Write();

    // Close the file. Note that this is automatically done when you leave
    // the application upon file destruction.
    hfile.Close();

    return 0;
}
~~~
*/

#include <ROOT/RConfig.hxx>
#include "TTree.h"

#include "ROOT/TIOFeatures.hxx"
#include "TArrayC.h"
#include "TBufferFile.h"
#include "TBaseClass.h"
#include "TBasket.h"
#include "TBranchClones.h"
#include "TBranchElement.h"
#include "TBranchObject.h"
#include "TBranchRef.h"
#include "TBrowser.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClonesArray.h"
#include "TCut.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TDirectory.h"
#include "TError.h"
#include "TEntryList.h"
#include "TEnv.h"
#include "TEventList.h"
#include "TFile.h"
#include "TFolder.h"
#include "TFriendElement.h"
#include "TInterpreter.h"
#include "TLeaf.h"
#include "TLeafB.h"
#include "TLeafC.h"
#include "TLeafD.h"
#include "TLeafElement.h"
#include "TLeafF.h"
#include "TLeafI.h"
#include "TLeafL.h"
#include "TLeafObject.h"
#include "TLeafS.h"
#include "TList.h"
#include "TMath.h"
#include "TMemFile.h"
#include "TROOT.h"
#include "TRealData.h"
#include "TRegexp.h"
#include "TRefTable.h"
#include "TStreamerElement.h"
#include "TStreamerInfo.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTreeCloner.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TVirtualCollectionProxy.h"
#include "TEmulatedCollectionProxy.h"
#include "TVirtualIndex.h"
#include "TVirtualPerfStats.h"
#include "TVirtualPad.h"
#include "TBranchSTL.h"
#include "TSchemaRuleSet.h"
#include "TFileMergeInfo.h"
#include "ROOT/StringConv.hxx"
#include "TVirtualMutex.h"
#include "strlcpy.h"
#include "snprintf.h"

#include "TBranchIMTHelper.h"
#include "TNotifyLink.h"

#include <chrono>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <climits>
#include <algorithm>
#include <set>

#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#include <thread>
#endif

constexpr Int_t   kNEntriesResort    = 100;
constexpr Float_t kNEntriesResortInv = 1.f/kNEntriesResort;

Int_t    TTree::fgBranchStyle = 1;  // Use new TBranch style with TBranchElement.
Long64_t TTree::fgMaxTreeSize = 100000000000LL;

ClassImp(TTree);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static char DataTypeToChar(EDataType datatype)
{
   // Return the leaflist 'char' for a given datatype.

   switch(datatype) {
   case kChar_t:     return 'B';
   case kUChar_t:    return 'b';
   case kBool_t:     return 'O';
   case kShort_t:    return 'S';
   case kUShort_t:   return 's';
   case kCounter:
   case kInt_t:      return 'I';
   case kUInt_t:     return 'i';
   case kDouble_t:   return 'D';
   case kDouble32_t: return 'd';
   case kFloat_t:    return 'F';
   case kFloat16_t:  return 'f';
   case kLong_t:     return 'G';
   case kULong_t:    return 'g';
   case kchar:       return 0; // unsupported
   case kLong64_t:   return 'L';
   case kULong64_t:  return 'l';

   case kCharStar:   return 'C';
   case kBits:       return 0; //unsupported

   case kOther_t:
   case kNoType_t:
   default:
      return 0;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// \class TTree::TFriendLock
/// Helper class to prevent infinite recursion in the usage of TTree Friends.

////////////////////////////////////////////////////////////////////////////////
/// Record in tree that it has been used while recursively looks through the friends.

TTree::TFriendLock::TFriendLock(TTree* tree, UInt_t methodbit)
: fTree(tree)
{
   // We could also add some code to acquire an actual
   // lock to prevent multi-thread issues
   fMethodBit = methodbit;
   if (fTree) {
      fPrevious = fTree->fFriendLockStatus & fMethodBit;
      fTree->fFriendLockStatus |= fMethodBit;
   } else {
      fPrevious = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TTree::TFriendLock::TFriendLock(const TFriendLock& tfl) :
  fTree(tfl.fTree),
  fMethodBit(tfl.fMethodBit),
  fPrevious(tfl.fPrevious)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TTree::TFriendLock& TTree::TFriendLock::operator=(const TTree::TFriendLock& tfl)
{
   if(this!=&tfl) {
      fTree=tfl.fTree;
      fMethodBit=tfl.fMethodBit;
      fPrevious=tfl.fPrevious;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Restore the state of tree the same as before we set the lock.

TTree::TFriendLock::~TFriendLock()
{
   if (fTree) {
      if (!fPrevious) {
         fTree->fFriendLockStatus &= ~(fMethodBit & kBitMask);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// \class TTree::TClusterIterator
/// Helper class to iterate over cluster of baskets.

////////////////////////////////////////////////////////////////////////////////
/// Regular constructor.
/// TTree is not set as const, since we might modify if it is a TChain.

TTree::TClusterIterator::TClusterIterator(TTree *tree, Long64_t firstEntry) : fTree(tree), fClusterRange(0), fStartEntry(0), fNextEntry(0), fEstimatedSize(-1)
{
   if (fTree->fNClusterRange) {
      // Find the correct cluster range.
      //
      // Since fClusterRangeEnd contains the inclusive upper end of the range, we need to search for the
      // range that was containing the previous entry and add 1 (because BinarySearch consider the values
      // to be the inclusive start of the bucket).
      fClusterRange = TMath::BinarySearch(fTree->fNClusterRange, fTree->fClusterRangeEnd, firstEntry - 1) + 1;

      Long64_t entryInRange;
      Long64_t pedestal;
      if (fClusterRange == 0) {
         pedestal = 0;
         entryInRange = firstEntry;
      } else {
         pedestal = fTree->fClusterRangeEnd[fClusterRange-1] + 1;
         entryInRange = firstEntry - pedestal;
      }
      Long64_t autoflush;
      if (fClusterRange == fTree->fNClusterRange) {
         autoflush = fTree->fAutoFlush;
      } else {
         autoflush = fTree->fClusterSize[fClusterRange];
      }
      if (autoflush <= 0) {
         autoflush = GetEstimatedClusterSize();
      }
      fStartEntry = pedestal + entryInRange - entryInRange%autoflush;
   } else if ( fTree->GetAutoFlush() <= 0 ) {
      // Case of old files before November 9 2009 *or* small tree where AutoFlush was never set.
      fStartEntry = firstEntry;
   } else {
      fStartEntry = firstEntry - firstEntry%fTree->GetAutoFlush();
   }
   fNextEntry = fStartEntry; // Position correctly for the first call to Next()
}

////////////////////////////////////////////////////////////////////////////////
/// Estimate the cluster size.
///
/// In almost all cases, this quickly returns the size of the auto-flush
/// in the TTree.
///
/// However, in the case where the cluster size was not fixed (old files and
/// case where autoflush was explicitly set to zero), we need estimate
/// a cluster size in relation to the size of the cache.
///
/// After this value is calculated once for the TClusterIterator, it is
/// cached and reused in future calls.

Long64_t TTree::TClusterIterator::GetEstimatedClusterSize()
{
   auto autoFlush = fTree->GetAutoFlush();
   if (autoFlush > 0) return autoFlush;
   if (fEstimatedSize > 0) return fEstimatedSize;

   Long64_t zipBytes = fTree->GetZipBytes();
   if (zipBytes == 0) {
      fEstimatedSize = fTree->GetEntries() - 1;
      if (fEstimatedSize <= 0)
         fEstimatedSize = 1;
   } else {
      Long64_t clusterEstimate = 1;
      Long64_t cacheSize = fTree->GetCacheSize();
      if (cacheSize == 0) {
         // Humm ... let's double check on the file.
         TFile *file = fTree->GetCurrentFile();
         if (file) {
            TFileCacheRead *cache = fTree->GetReadCache(file);
            if (cache) {
               cacheSize = cache->GetBufferSize();
            }
         }
      }
      // If neither file nor tree has a cache, use the current default.
      if (cacheSize <= 0) {
         cacheSize = 30000000;
      }
      clusterEstimate = fTree->GetEntries() * cacheSize / zipBytes;
      // If there are no entries, then just default to 1.
      fEstimatedSize = clusterEstimate ? clusterEstimate : 1;
   }
   return fEstimatedSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Move on to the next cluster and return the starting entry
/// of this next cluster

Long64_t TTree::TClusterIterator::Next()
{
   fStartEntry = fNextEntry;
   if (fTree->fNClusterRange || fTree->GetAutoFlush() > 0) {
      if (fClusterRange == fTree->fNClusterRange) {
         // We are looking at a range which size
         // is defined by AutoFlush itself and goes to the GetEntries.
         fNextEntry += GetEstimatedClusterSize();
      } else {
         if (fStartEntry > fTree->fClusterRangeEnd[fClusterRange]) {
            ++fClusterRange;
         }
         if (fClusterRange == fTree->fNClusterRange) {
            // We are looking at the last range which size
            // is defined by AutoFlush itself and goes to the GetEntries.
            fNextEntry += GetEstimatedClusterSize();
         } else {
            Long64_t clusterSize = fTree->fClusterSize[fClusterRange];
            if (clusterSize == 0) {
               clusterSize = GetEstimatedClusterSize();
            }
            fNextEntry += clusterSize;
            if (fNextEntry > fTree->fClusterRangeEnd[fClusterRange]) {
               // The last cluster of the range was a partial cluster,
               // so the next cluster starts at the beginning of the
               // next range.
               fNextEntry = fTree->fClusterRangeEnd[fClusterRange] + 1;
            }
         }
      }
   } else {
      // Case of old files before November 9 2009
      fNextEntry = fStartEntry + GetEstimatedClusterSize();
   }
   if (fNextEntry > fTree->GetEntries()) {
      fNextEntry = fTree->GetEntries();
   }
   return fStartEntry;
}

////////////////////////////////////////////////////////////////////////////////
/// Move on to the previous cluster and return the starting entry
/// of this previous cluster

Long64_t TTree::TClusterIterator::Previous()
{
   fNextEntry = fStartEntry;
   if (fTree->fNClusterRange || fTree->GetAutoFlush() > 0) {
      if (fClusterRange == 0 || fTree->fNClusterRange == 0) {
         // We are looking at a range which size
         // is defined by AutoFlush itself.
         fStartEntry -= GetEstimatedClusterSize();
      } else {
         if (fNextEntry <= fTree->fClusterRangeEnd[fClusterRange]) {
            --fClusterRange;
         }
         if (fClusterRange == 0) {
            // We are looking at the first range.
            fStartEntry = 0;
         } else {
            Long64_t clusterSize = fTree->fClusterSize[fClusterRange];
            if (clusterSize == 0) {
               clusterSize = GetEstimatedClusterSize();
            }
            fStartEntry -= clusterSize;
         }
      }
   } else {
      // Case of old files before November 9 2009 or trees that never auto-flushed.
      fStartEntry = fNextEntry - GetEstimatedClusterSize();
   }
   if (fStartEntry < 0) {
      fStartEntry = 0;
   }
   return fStartEntry;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Default constructor and I/O constructor.
///
/// Note: We do *not* insert ourself into the current directory.
///

TTree::TTree()
: TNamed()
, TAttLine()
, TAttFill()
, TAttMarker()
, fEntries(0)
, fTotBytes(0)
, fZipBytes(0)
, fSavedBytes(0)
, fFlushedBytes(0)
, fWeight(1)
, fTimerInterval(0)
, fScanField(25)
, fUpdate(0)
, fDefaultEntryOffsetLen(1000)
, fNClusterRange(0)
, fMaxClusterRange(0)
, fMaxEntries(0)
, fMaxEntryLoop(0)
, fMaxVirtualSize(0)
, fAutoSave( -300000000)
, fAutoFlush(-30000000)
, fEstimate(1000000)
, fClusterRangeEnd(0)
, fClusterSize(0)
, fCacheSize(0)
, fChainOffset(0)
, fReadEntry(-1)
, fTotalBuffers(0)
, fPacketSize(100)
, fNfill(0)
, fDebug(0)
, fDebugMin(0)
, fDebugMax(9999999)
, fMakeClass(0)
, fFileNumber(0)
, fNotify(0)
, fDirectory(0)
, fBranches()
, fLeaves()
, fAliases(0)
, fEventList(0)
, fEntryList(0)
, fIndexValues()
, fIndex()
, fTreeIndex(0)
, fFriends(0)
, fExternalFriends(0)
, fPerfStats(0)
, fUserInfo(0)
, fPlayer(0)
, fClones(0)
, fBranchRef(0)
, fFriendLockStatus(0)
, fTransientBuffer(0)
, fCacheDoAutoInit(kTRUE)
, fCacheDoClusterPrefetch(kFALSE)
, fCacheUserSet(kFALSE)
, fIMTEnabled(ROOT::IsImplicitMTEnabled())
, fNEntriesSinceSorting(0)
{
   fMaxEntries = 1000000000;
   fMaxEntries *= 1000;

   fMaxEntryLoop = 1000000000;
   fMaxEntryLoop *= 1000;

   fBranches.SetOwner(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal tree constructor.
///
/// The tree is created in the current directory.
/// Use the various functions Branch below to add branches to this tree.
///
/// If the first character of title is a "/", the function assumes a folder name.
/// In this case, it creates automatically branches following the folder hierarchy.
/// splitlevel may be used in this case to control the split level.

TTree::TTree(const char* name, const char* title, Int_t splitlevel /* = 99 */,
             TDirectory* dir /* = gDirectory*/)
: TNamed(name, title)
, TAttLine()
, TAttFill()
, TAttMarker()
, fEntries(0)
, fTotBytes(0)
, fZipBytes(0)
, fSavedBytes(0)
, fFlushedBytes(0)
, fWeight(1)
, fTimerInterval(0)
, fScanField(25)
, fUpdate(0)
, fDefaultEntryOffsetLen(1000)
, fNClusterRange(0)
, fMaxClusterRange(0)
, fMaxEntries(0)
, fMaxEntryLoop(0)
, fMaxVirtualSize(0)
, fAutoSave( -300000000)
, fAutoFlush(-30000000)
, fEstimate(1000000)
, fClusterRangeEnd(0)
, fClusterSize(0)
, fCacheSize(0)
, fChainOffset(0)
, fReadEntry(-1)
, fTotalBuffers(0)
, fPacketSize(100)
, fNfill(0)
, fDebug(0)
, fDebugMin(0)
, fDebugMax(9999999)
, fMakeClass(0)
, fFileNumber(0)
, fNotify(0)
, fDirectory(dir)
, fBranches()
, fLeaves()
, fAliases(0)
, fEventList(0)
, fEntryList(0)
, fIndexValues()
, fIndex()
, fTreeIndex(0)
, fFriends(0)
, fExternalFriends(0)
, fPerfStats(0)
, fUserInfo(0)
, fPlayer(0)
, fClones(0)
, fBranchRef(0)
, fFriendLockStatus(0)
, fTransientBuffer(0)
, fCacheDoAutoInit(kTRUE)
, fCacheDoClusterPrefetch(kFALSE)
, fCacheUserSet(kFALSE)
, fIMTEnabled(ROOT::IsImplicitMTEnabled())
, fNEntriesSinceSorting(0)
{
   // TAttLine state.
   SetLineColor(gStyle->GetHistLineColor());
   SetLineStyle(gStyle->GetHistLineStyle());
   SetLineWidth(gStyle->GetHistLineWidth());

   // TAttFill state.
   SetFillColor(gStyle->GetHistFillColor());
   SetFillStyle(gStyle->GetHistFillStyle());

   // TAttMarkerState.
   SetMarkerColor(gStyle->GetMarkerColor());
   SetMarkerStyle(gStyle->GetMarkerStyle());
   SetMarkerSize(gStyle->GetMarkerSize());

   fMaxEntries = 1000000000;
   fMaxEntries *= 1000;

   fMaxEntryLoop = 1000000000;
   fMaxEntryLoop *= 1000;

   // Insert ourself into the current directory.
   // FIXME: This is very annoying behaviour, we should
   //        be able to choose to not do this like we
   //        can with a histogram.
   if (fDirectory) fDirectory->Append(this);

   fBranches.SetOwner(kTRUE);

   // If title starts with "/" and is a valid folder name, a superbranch
   // is created.
   // FIXME: Why?
   if (strlen(title) > 2) {
      if (title[0] == '/') {
         Branch(title+1,32000,splitlevel);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TTree::~TTree()
{
   if (auto link = dynamic_cast<TNotifyLinkBase*>(fNotify)) {
      link->Clear();
   }
   if (fAllocationCount && (gDebug > 0)) {
      Info("TTree::~TTree", "For tree %s, allocation count is %u.", GetName(), fAllocationCount.load());
#ifdef R__TRACK_BASKET_ALLOC_TIME
      Info("TTree::~TTree", "For tree %s, allocation time is %lluus.", GetName(), fAllocationTime.load());
#endif
   }

   if (fDirectory) {
      // We are in a directory, which may possibly be a file.
      if (fDirectory->GetList()) {
         // Remove us from the directory listing.
         fDirectory->Remove(this);
      }
      //delete the file cache if it points to this Tree
      TFile *file = fDirectory->GetFile();
      MoveReadCache(file,0);
   }
   // We don't own the leaves in fLeaves, the branches do.
   fLeaves.Clear();
   // I'm ready to destroy any objects allocated by
   // SetAddress() by my branches.  If I have clones,
   // tell them to zero their pointers to this shared
   // memory.
   if (fClones && fClones->GetEntries()) {
      // I have clones.
      // I am about to delete the objects created by
      // SetAddress() which we are sharing, so tell
      // the clones to release their pointers to them.
      for (TObjLink* lnk = fClones->FirstLink(); lnk; lnk = lnk->Next()) {
         TTree* clone = (TTree*) lnk->GetObject();
         // clone->ResetBranchAddresses();

         // Reset only the branch we have set the address of.
         CopyAddresses(clone,kTRUE);
      }
   }
   // Get rid of our branches, note that this will also release
   // any memory allocated by TBranchElement::SetAddress().
   fBranches.Delete();
   // FIXME: We must consider what to do with the reset of these if we are a clone.
   delete fPlayer;
   fPlayer = 0;
   if (fExternalFriends) {
      using namespace ROOT::Detail;
      for(auto fetree : TRangeStaticCast<TFriendElement>(*fExternalFriends))
         fetree->Reset();
      fExternalFriends->Clear("nodelete");
      SafeDelete(fExternalFriends);
   }
   if (fFriends) {
      fFriends->Delete();
      delete fFriends;
      fFriends = 0;
   }
   if (fAliases) {
      fAliases->Delete();
      delete fAliases;
      fAliases = 0;
   }
   if (fUserInfo) {
      fUserInfo->Delete();
      delete fUserInfo;
      fUserInfo = 0;
   }
   if (fClones) {
      // Clone trees should no longer be removed from fClones when they are deleted.
     {
        R__LOCKGUARD(gROOTMutex);
        gROOT->GetListOfCleanups()->Remove(fClones);
     }
      // Note: fClones does not own its content.
      delete fClones;
      fClones = 0;
   }
   if (fEntryList) {
      if (fEntryList->TestBit(kCanDelete) && fEntryList->GetDirectory()==0) {
         // Delete the entry list if it is marked to be deleted and it is not also
         // owned by a directory.  (Otherwise we would need to make sure that a
         // TDirectoryFile that has a TTree in it does a 'slow' TList::Delete.
         delete fEntryList;
         fEntryList=0;
      }
   }
   delete fTreeIndex;
   fTreeIndex = 0;
   delete fBranchRef;
   fBranchRef = 0;
   delete [] fClusterRangeEnd;
   fClusterRangeEnd = 0;
   delete [] fClusterSize;
   fClusterSize = 0;
   // Must be done after the destruction of friends.
   // Note: We do *not* own our directory.
   fDirectory = 0;

   if (fTransientBuffer) {
      delete fTransientBuffer;
      fTransientBuffer = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the transient buffer currently used by this TTree for reading/writing baskets.

TBuffer* TTree::GetTransientBuffer(Int_t size)
{
   if (fTransientBuffer) {
      if (fTransientBuffer->BufferSize() < size) {
         fTransientBuffer->Expand(size);
      }
      return fTransientBuffer;
   }
   fTransientBuffer = new TBufferFile(TBuffer::kRead, size);
   return fTransientBuffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Add branch with name bname to the Tree cache.
/// If bname="*" all branches are added to the cache.
/// if subbranches is true all the branches of the subbranches are
/// also put to the cache.
///
/// Returns:
/// - 0 branch added or already included
/// - -1 on error

Int_t TTree::AddBranchToCache(const char*bname, Bool_t subbranches)
{
   if (!GetTree()) {
      if (LoadTree(0)<0) {
         Error("AddBranchToCache","Could not load a tree");
         return -1;
      }
   }
   if (GetTree()) {
      if (GetTree() != this) {
         return GetTree()->AddBranchToCache(bname, subbranches);
      }
   } else {
      Error("AddBranchToCache", "No tree is available. Branch was not added to the cache");
      return -1;
   }

   TFile *f = GetCurrentFile();
   if (!f) {
      Error("AddBranchToCache", "No file is available. Branch was not added to the cache");
      return -1;
   }
   TTreeCache *tc = GetReadCache(f,kTRUE);
   if (!tc) {
      Error("AddBranchToCache", "No cache is available, branch not added");
      return -1;
   }
   return tc->AddBranch(bname,subbranches);
}

////////////////////////////////////////////////////////////////////////////////
/// Add branch b to the Tree cache.
/// if subbranches is true all the branches of the subbranches are
/// also put to the cache.
///
/// Returns:
/// - 0 branch added or already included
/// - -1 on error

Int_t TTree::AddBranchToCache(TBranch *b, Bool_t subbranches)
{
   if (!GetTree()) {
      if (LoadTree(0)<0) {
         Error("AddBranchToCache","Could not load a tree");
         return -1;
      }
   }
   if (GetTree()) {
      if (GetTree() != this) {
         Int_t res = GetTree()->AddBranchToCache(b, subbranches);
         if (res<0) {
             Error("AddBranchToCache", "Error adding branch");
         }
         return res;
      }
   } else {
      Error("AddBranchToCache", "No tree is available. Branch was not added to the cache");
      return -1;
   }

   TFile *f = GetCurrentFile();
   if (!f) {
      Error("AddBranchToCache", "No file is available. Branch was not added to the cache");
      return -1;
   }
   TTreeCache *tc = GetReadCache(f,kTRUE);
   if (!tc) {
      Error("AddBranchToCache", "No cache is available, branch not added");
      return -1;
   }
   return tc->AddBranch(b,subbranches);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the branch with name 'bname' from the Tree cache.
/// If bname="*" all branches are removed from the cache.
/// if subbranches is true all the branches of the subbranches are
/// also removed from the cache.
///
/// Returns:
/// - 0 branch dropped or not in cache
/// - -1 on error

Int_t TTree::DropBranchFromCache(const char*bname, Bool_t subbranches)
{
   if (!GetTree()) {
      if (LoadTree(0)<0) {
         Error("DropBranchFromCache","Could not load a tree");
         return -1;
      }
   }
   if (GetTree()) {
      if (GetTree() != this) {
         return GetTree()->DropBranchFromCache(bname, subbranches);
      }
   } else {
      Error("DropBranchFromCache", "No tree is available. Branch was not dropped from the cache");
      return -1;
   }

   TFile *f = GetCurrentFile();
   if (!f) {
      Error("DropBranchFromCache", "No file is available. Branch was not dropped from the cache");
      return -1;
   }
   TTreeCache *tc = GetReadCache(f,kTRUE);
   if (!tc) {
      Error("DropBranchFromCache", "No cache is available, branch not dropped");
      return -1;
   }
   return tc->DropBranch(bname,subbranches);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the branch b from the Tree cache.
/// if subbranches is true all the branches of the subbranches are
/// also removed from the cache.
///
/// Returns:
/// - 0 branch dropped or not in cache
/// - -1 on error

Int_t TTree::DropBranchFromCache(TBranch *b, Bool_t subbranches)
{
   if (!GetTree()) {
      if (LoadTree(0)<0) {
         Error("DropBranchFromCache","Could not load a tree");
         return -1;
      }
   }
   if (GetTree()) {
      if (GetTree() != this) {
         Int_t res = GetTree()->DropBranchFromCache(b, subbranches);
         if (res<0) {
             Error("DropBranchFromCache", "Error dropping branch");
         }
         return res;
      }
   } else {
      Error("DropBranchFromCache", "No tree is available. Branch was not dropped from the cache");
      return -1;
   }

   TFile *f = GetCurrentFile();
   if (!f) {
      Error("DropBranchFromCache", "No file is available. Branch was not dropped from the cache");
      return -1;
   }
   TTreeCache *tc = GetReadCache(f,kTRUE);
   if (!tc) {
      Error("DropBranchFromCache", "No cache is available, branch not dropped");
      return -1;
   }
   return tc->DropBranch(b,subbranches);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a cloned tree to our list of trees to be notified whenever we change
/// our branch addresses or when we are deleted.

void TTree::AddClone(TTree* clone)
{
   if (!fClones) {
      fClones = new TList();
      fClones->SetOwner(false);
      // So that the clones are automatically removed from the list when
      // they are deleted.
      {
         R__LOCKGUARD(gROOTMutex);
         gROOT->GetListOfCleanups()->Add(fClones);
      }
   }
   if (!fClones->FindObject(clone)) {
      fClones->Add(clone);
   }
}

// Check whether mainTree and friendTree can be friends w.r.t. the kEntriesReshuffled bit.
// In particular, if any has the bit set, then friendTree must have a TTreeIndex and the
// branches used for indexing must be present in mainTree.
// Return true if the trees can be friends, false otherwise.
bool CheckReshuffling(TTree &mainTree, TTree &friendTree)
{
   const auto isMainReshuffled = mainTree.TestBit(TTree::kEntriesReshuffled);
   const auto isFriendReshuffled = friendTree.TestBit(TTree::kEntriesReshuffled);
   const auto friendHasValidIndex = [&] {
      auto idx = friendTree.GetTreeIndex();
      return idx ? idx->IsValidFor(&mainTree) : kFALSE;
   }();

   if ((isMainReshuffled || isFriendReshuffled) && !friendHasValidIndex) {
      const auto reshuffledTreeName = isMainReshuffled ? mainTree.GetName() : friendTree.GetName();
      const auto msg = "Tree '%s' has the kEntriesReshuffled bit set, and cannot be used as friend nor can be added as "
                       "a friend unless the main tree has a TTreeIndex on the friend tree '%s'. You can also unset the "
                       "bit manually if you know what you are doing.";
      Error("AddFriend", msg, reshuffledTreeName, friendTree.GetName());
      return false;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TFriendElement to the list of friends.
///
/// This function:
/// - opens a file if filename is specified
/// - reads a Tree with name treename from the file (current directory)
/// - adds the Tree to the list of friends
/// see other AddFriend functions
///
/// A TFriendElement TF describes a TTree object TF in a file.
/// When a TFriendElement TF is added to the the list of friends of an
/// existing TTree T, any variable from TF can be referenced in a query
/// to T.
///
///   A tree keeps a list of friends. In the context of a tree (or a chain),
/// friendship means unrestricted access to the friends data. In this way
/// it is much like adding another branch to the tree without taking the risk
/// of damaging it. To add a friend to the list, you can use the TTree::AddFriend
/// method.  The tree in the diagram below has two friends (friend_tree1 and
/// friend_tree2) and now has access to the variables a,b,c,i,j,k,l and m.
///
/// \image html ttree_friend1.png
///
/// The AddFriend method has two parameters, the first is the tree name and the
/// second is the name of the ROOT file where the friend tree is saved.
/// AddFriend automatically opens the friend file. If no file name is given,
/// the tree called ft1 is assumed to be in the same file as the original tree.
///
/// tree.AddFriend("ft1","friendfile1.root");
/// If the friend tree has the same name as the original tree, you can give it
/// an alias in the context of the friendship:
///
/// tree.AddFriend("tree1 = tree","friendfile1.root");
/// Once the tree has friends, we can use TTree::Draw as if the friend's
/// variables were in the original tree. To specify which tree to use in
/// the Draw method, use the syntax:
/// ~~~ {.cpp}
///     <treeName>.<branchname>.<varname>
/// ~~~
/// If the variablename is enough to uniquely identify the variable, you can
/// leave out the tree and/or branch name.
/// For example, these commands generate a 3-d scatter plot of variable "var"
/// in the TTree tree versus variable v1 in TTree ft1 versus variable v2 in
/// TTree ft2.
/// ~~~ {.cpp}
///     tree.AddFriend("ft1","friendfile1.root");
///     tree.AddFriend("ft2","friendfile2.root");
///     tree.Draw("var:ft1.v1:ft2.v2");
/// ~~~
/// \image html ttree_friend2.png
///
/// The picture illustrates the access of the tree and its friends with a
/// Draw command.
/// When AddFriend is called, the ROOT file is automatically opened and the
/// friend tree (ft1) is read into memory. The new friend (ft1) is added to
/// the list of friends of tree.
/// The number of entries in the friend must be equal or greater to the number
/// of entries of the original tree. If the friend tree has fewer entries a
/// warning is given and the missing entries are not included in the histogram.
/// To retrieve the list of friends from a tree use TTree::GetListOfFriends.
/// When the tree is written to file (TTree::Write), the friends list is saved
/// with it. And when the tree is retrieved, the trees on the friends list are
/// also retrieved and the friendship restored.
/// When a tree is deleted, the elements of the friend list are also deleted.
/// It is possible to declare a friend tree that has the same internal
/// structure (same branches and leaves) as the original tree, and compare the
/// same values by specifying the tree.
/// ~~~ {.cpp}
///     tree.Draw("var:ft1.var:ft2.var")
/// ~~~

TFriendElement *TTree::AddFriend(const char *treename, const char *filename)
{
   if (!fFriends) {
      fFriends = new TList();
   }
   TFriendElement *fe = new TFriendElement(this, treename, filename);

   TTree *t = fe->GetTree();
   bool canAddFriend = true;
   if (t) {
      canAddFriend = CheckReshuffling(*this, *t);
      if (!t->GetTreeIndex() && (t->GetEntries() < fEntries)) {
         Warning("AddFriend", "FriendElement %s in file %s has less entries %lld than its parent Tree: %lld", treename,
                 filename, t->GetEntries(), fEntries);
      }
   } else {
      Error("AddFriend", "Cannot find tree '%s' in file '%s', friend not added", treename, filename);
      canAddFriend = false;
   }

   if (canAddFriend)
      fFriends->Add(fe);
   return fe;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TFriendElement to the list of friends.
///
/// The TFile is managed by the user (e.g. the user must delete the file).
/// For complete description see AddFriend(const char *, const char *).
/// This function:
/// - reads a Tree with name treename from the file
/// - adds the Tree to the list of friends

TFriendElement *TTree::AddFriend(const char *treename, TFile *file)
{
   if (!fFriends) {
      fFriends = new TList();
   }
   TFriendElement *fe = new TFriendElement(this, treename, file);
   R__ASSERT(fe);
   TTree *t = fe->GetTree();
   bool canAddFriend = true;
   if (t) {
      canAddFriend = CheckReshuffling(*this, *t);
      if (!t->GetTreeIndex() && (t->GetEntries() < fEntries)) {
         Warning("AddFriend", "FriendElement %s in file %s has less entries %lld than its parent tree: %lld", treename,
                 file->GetName(), t->GetEntries(), fEntries);
      }
   } else {
      Error("AddFriend", "Cannot find tree '%s' in file '%s', friend not added", treename, file->GetName());
      canAddFriend = false;
   }

   if (canAddFriend)
      fFriends->Add(fe);
   return fe;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TFriendElement to the list of friends.
///
/// The TTree is managed by the user (e.g., the user must delete the file).
/// For a complete description see AddFriend(const char *, const char *).

TFriendElement *TTree::AddFriend(TTree *tree, const char *alias, Bool_t warn)
{
   if (!tree) {
      return 0;
   }
   if (!fFriends) {
      fFriends = new TList();
   }
   TFriendElement *fe = new TFriendElement(this, tree, alias);
   R__ASSERT(fe); // this assert is for historical reasons. Don't remove it unless you understand all the consequences.
   TTree *t = fe->GetTree();
   if (warn && (t->GetEntries() < fEntries)) {
      Warning("AddFriend", "FriendElement '%s' in file '%s' has less entries %lld than its parent tree: %lld",
              tree->GetName(), fe->GetFile() ? fe->GetFile()->GetName() : "(memory resident)", t->GetEntries(),
              fEntries);
   }
   if (CheckReshuffling(*this, *t)) {
      fFriends->Add(fe);
      tree->RegisterExternalFriend(fe);
   }
   return fe;
}

////////////////////////////////////////////////////////////////////////////////
/// AutoSave tree header every fAutoSave bytes.
///
/// When large Trees are produced, it is safe to activate the AutoSave
/// procedure. Some branches may have buffers holding many entries.
/// If fAutoSave is negative, AutoSave is automatically called by
/// TTree::Fill when the number of bytes generated since the previous
/// AutoSave is greater than -fAutoSave bytes.
/// If fAutoSave is positive, AutoSave is automatically called by
/// TTree::Fill every N entries.
/// This function may also be invoked by the user.
/// Each AutoSave generates a new key on the file.
/// Once the key with the tree header has been written, the previous cycle
/// (if any) is deleted.
///
/// Note that calling TTree::AutoSave too frequently (or similarly calling
/// TTree::SetAutoSave with a small value) is an expensive operation.
/// You should make tests for your own application to find a compromise
/// between speed and the quantity of information you may loose in case of
/// a job crash.
///
/// In case your program crashes before closing the file holding this tree,
/// the file will be automatically recovered when you will connect the file
/// in UPDATE mode.
/// The Tree will be recovered at the status corresponding to the last AutoSave.
///
/// if option contains "SaveSelf", gDirectory->SaveSelf() is called.
/// This allows another process to analyze the Tree while the Tree is being filled.
///
/// if option contains "FlushBaskets", TTree::FlushBaskets is called and all
/// the current basket are closed-out and written to disk individually.
///
/// By default the previous header is deleted after having written the new header.
/// if option contains "Overwrite", the previous Tree header is deleted
/// before written the new header. This option is slightly faster, but
/// the default option is safer in case of a problem (disk quota exceeded)
/// when writing the new header.
///
/// The function returns the number of bytes written to the file.
/// if the number of bytes is null, an error has occurred while writing
/// the header to the file.
///
/// ## How to write a Tree in one process and view it from another process
///
/// The following two scripts illustrate how to do this.
/// The script treew.C is executed by process1, treer.C by process2
///
/// script treew.C:
/// ~~~ {.cpp}
///     void treew() {
///        TFile f("test.root","recreate");
///        TNtuple *ntuple = new TNtuple("ntuple","Demo","px:py:pz:random:i");
///        Float_t px, py, pz;
///        for ( Int_t i=0; i<10000000; i++) {
///           gRandom->Rannor(px,py);
///           pz = px*px + py*py;
///           Float_t random = gRandom->Rndm(1);
///           ntuple->Fill(px,py,pz,random,i);
///           if (i%1000 == 1) ntuple->AutoSave("SaveSelf");
///        }
///     }
/// ~~~
/// script treer.C:
/// ~~~ {.cpp}
///     void treer() {
///        TFile f("test.root");
///        TTree *ntuple = (TTree*)f.Get("ntuple");
///        TCanvas c1;
///        Int_t first = 0;
///        while(1) {
///           if (first == 0) ntuple->Draw("px>>hpx", "","",10000000,first);
///           else            ntuple->Draw("px>>+hpx","","",10000000,first);
///           first = (Int_t)ntuple->GetEntries();
///           c1.Update();
///           gSystem->Sleep(1000); //sleep 1 second
///           ntuple->Refresh();
///        }
///     }
/// ~~~

Long64_t TTree::AutoSave(Option_t* option)
{
   if (!fDirectory || fDirectory == gROOT || !fDirectory->IsWritable()) return 0;
   if (gDebug > 0) {
      Info("AutoSave", "Tree:%s after %lld bytes written\n",GetName(),GetTotBytes());
   }
   TString opt = option;
   opt.ToLower();

   if (opt.Contains("flushbaskets")) {
      if (gDebug > 0) Info("AutoSave", "calling FlushBaskets \n");
      FlushBasketsImpl();
   }

   fSavedBytes = GetZipBytes();

   TKey *key = (TKey*)fDirectory->GetListOfKeys()->FindObject(GetName());
   Long64_t nbytes;
   if (opt.Contains("overwrite")) {
      nbytes = fDirectory->WriteTObject(this,"","overwrite");
   } else {
      nbytes = fDirectory->WriteTObject(this); //nbytes will be 0 if Write failed (disk space exceeded)
      if (nbytes && key && strcmp(ClassName(), key->GetClassName()) == 0) {
         key->Delete();
         delete key;
      }
   }
   // save StreamerInfo
   TFile *file = fDirectory->GetFile();
   if (file) file->WriteStreamerInfo();

   if (opt.Contains("saveself")) {
      fDirectory->SaveSelf();
      //the following line is required in case GetUserInfo contains a user class
      //for which the StreamerInfo must be written. One could probably be a bit faster (Rene)
      if (file) file->WriteHeader();
   }

   return nbytes;
}

namespace {
   // This error message is repeated several times in the code. We write it once.
   const char* writeStlWithoutProxyMsg = "The class requested (%s) for the branch \"%s\""
                                      " is an instance of an stl collection and does not have a compiled CollectionProxy."
                                      " Please generate the dictionary for this collection (%s) to avoid to write corrupted data.";
}

////////////////////////////////////////////////////////////////////////////////
/// Same as TTree::Branch() with added check that addobj matches className.
///
/// See TTree::Branch() for other details.
///

TBranch* TTree::BranchImp(const char* branchname, const char* classname, TClass* ptrClass, void* addobj, Int_t bufsize, Int_t splitlevel)
{
   TClass* claim = TClass::GetClass(classname);
   if (!ptrClass) {
      if (claim && claim->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(claim->GetCollectionProxy())) {
         Error("Branch", writeStlWithoutProxyMsg,
               claim->GetName(), branchname, claim->GetName());
         return 0;
      }
      return Branch(branchname, classname, (void*) addobj, bufsize, splitlevel);
   }
   TClass* actualClass = 0;
   void** addr = (void**) addobj;
   if (addr) {
      actualClass = ptrClass->GetActualClass(*addr);
   }
   if (ptrClass && claim) {
      if (!(claim->InheritsFrom(ptrClass) || ptrClass->InheritsFrom(claim))) {
         // Note we currently do not warn in case of splicing or over-expectation).
         if (claim->IsLoaded() && ptrClass->IsLoaded() && strcmp( claim->GetTypeInfo()->name(), ptrClass->GetTypeInfo()->name() ) == 0) {
            // The type is the same according to the C++ type_info, we must be in the case of
            // a template of Double32_t.  This is actually a correct case.
         } else {
            Error("Branch", "The class requested (%s) for \"%s\" is different from the type of the pointer passed (%s)",
                  claim->GetName(), branchname, ptrClass->GetName());
         }
      } else if (actualClass && (claim != actualClass) && !actualClass->InheritsFrom(claim)) {
         if (claim->IsLoaded() && actualClass->IsLoaded() && strcmp( claim->GetTypeInfo()->name(), actualClass->GetTypeInfo()->name() ) == 0) {
            // The type is the same according to the C++ type_info, we must be in the case of
            // a template of Double32_t.  This is actually a correct case.
         } else {
            Error("Branch", "The actual class (%s) of the object provided for the definition of the branch \"%s\" does not inherit from %s",
                  actualClass->GetName(), branchname, claim->GetName());
         }
      }
   }
   if (claim && claim->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(claim->GetCollectionProxy())) {
      Error("Branch", writeStlWithoutProxyMsg,
            claim->GetName(), branchname, claim->GetName());
      return 0;
   }
   return Branch(branchname, classname, (void*) addobj, bufsize, splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Same as TTree::Branch but automatic detection of the class name.
/// See TTree::Branch for other details.

TBranch* TTree::BranchImp(const char* branchname, TClass* ptrClass, void* addobj, Int_t bufsize, Int_t splitlevel)
{
   if (!ptrClass) {
      Error("Branch", "The pointer specified for %s is not of a class known to ROOT", branchname);
      return 0;
   }
   TClass* actualClass = 0;
   void** addr = (void**) addobj;
   if (addr && *addr) {
      actualClass = ptrClass->GetActualClass(*addr);
      if (!actualClass) {
         Warning("Branch", "The actual TClass corresponding to the object provided for the definition of the branch \"%s\" is missing.\n\tThe object will be truncated down to its %s part",
                 branchname, ptrClass->GetName());
         actualClass = ptrClass;
      } else if ((ptrClass != actualClass) && !actualClass->InheritsFrom(ptrClass)) {
         Error("Branch", "The actual class (%s) of the object provided for the definition of the branch \"%s\" does not inherit from %s", actualClass->GetName(), branchname, ptrClass->GetName());
         return 0;
      }
   } else {
      actualClass = ptrClass;
   }
   if (actualClass && actualClass->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(actualClass->GetCollectionProxy())) {
      Error("Branch", writeStlWithoutProxyMsg,
            actualClass->GetName(), branchname, actualClass->GetName());
      return 0;
   }
   return Branch(branchname, actualClass->GetName(), (void*) addobj, bufsize, splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Same as TTree::Branch but automatic detection of the class name.
/// See TTree::Branch for other details.

TBranch* TTree::BranchImpRef(const char* branchname, const char *classname, TClass* ptrClass, void *addobj, Int_t bufsize, Int_t splitlevel)
{
   TClass* claim = TClass::GetClass(classname);
   if (!ptrClass) {
      if (claim && claim->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(claim->GetCollectionProxy())) {
         Error("Branch", writeStlWithoutProxyMsg,
               claim->GetName(), branchname, claim->GetName());
         return 0;
      } else if (claim == 0) {
         Error("Branch", "The pointer specified for %s is not of a class known to ROOT and %s is not a known class", branchname, classname);
         return 0;
      }
      ptrClass = claim;
   }
   TClass* actualClass = 0;
   if (!addobj) {
      Error("Branch", "Reference interface requires a valid object (for branch: %s)!", branchname);
      return 0;
   }
   actualClass = ptrClass->GetActualClass(addobj);
   if (ptrClass && claim) {
      if (!(claim->InheritsFrom(ptrClass) || ptrClass->InheritsFrom(claim))) {
         // Note we currently do not warn in case of splicing or over-expectation).
         if (claim->IsLoaded() && ptrClass->IsLoaded() && strcmp( claim->GetTypeInfo()->name(), ptrClass->GetTypeInfo()->name() ) == 0) {
            // The type is the same according to the C++ type_info, we must be in the case of
            // a template of Double32_t.  This is actually a correct case.
         } else {
            Error("Branch", "The class requested (%s) for \"%s\" is different from the type of the object passed (%s)",
                  claim->GetName(), branchname, ptrClass->GetName());
         }
      } else if (actualClass && (claim != actualClass) && !actualClass->InheritsFrom(claim)) {
         if (claim->IsLoaded() && actualClass->IsLoaded() && strcmp( claim->GetTypeInfo()->name(), actualClass->GetTypeInfo()->name() ) == 0) {
            // The type is the same according to the C++ type_info, we must be in the case of
            // a template of Double32_t.  This is actually a correct case.
         } else {
            Error("Branch", "The actual class (%s) of the object provided for the definition of the branch \"%s\" does not inherit from %s",
                  actualClass->GetName(), branchname, claim->GetName());
         }
      }
   }
   if (!actualClass) {
      Warning("Branch", "The actual TClass corresponding to the object provided for the definition of the branch \"%s\" is missing.\n\tThe object will be truncated down to its %s part",
              branchname, ptrClass->GetName());
      actualClass = ptrClass;
   } else if ((ptrClass != actualClass) && !actualClass->InheritsFrom(ptrClass)) {
      Error("Branch", "The actual class (%s) of the object provided for the definition of the branch \"%s\" does not inherit from %s", actualClass->GetName(), branchname, ptrClass->GetName());
      return 0;
   }
   if (actualClass && actualClass->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(actualClass->GetCollectionProxy())) {
      Error("Branch", writeStlWithoutProxyMsg,
            actualClass->GetName(), branchname, actualClass->GetName());
      return 0;
   }
   return BronchExec(branchname, actualClass->GetName(), (void*) addobj, kFALSE, bufsize, splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Same as TTree::Branch but automatic detection of the class name.
/// See TTree::Branch for other details.

TBranch* TTree::BranchImpRef(const char* branchname, TClass* ptrClass, EDataType datatype, void* addobj, Int_t bufsize, Int_t splitlevel)
{
   if (!ptrClass) {
      if (datatype == kOther_t || datatype == kNoType_t) {
         Error("Branch", "The pointer specified for %s is not of a class or type known to ROOT", branchname);
      } else {
         TString varname; varname.Form("%s/%c",branchname,DataTypeToChar(datatype));
         return Branch(branchname,addobj,varname.Data(),bufsize);
      }
      return 0;
   }
   TClass* actualClass = 0;
   if (!addobj) {
      Error("Branch", "Reference interface requires a valid object (for branch: %s)!", branchname);
      return 0;
   }
   actualClass = ptrClass->GetActualClass(addobj);
   if (!actualClass) {
      Warning("Branch", "The actual TClass corresponding to the object provided for the definition of the branch \"%s\" is missing.\n\tThe object will be truncated down to its %s part",
              branchname, ptrClass->GetName());
      actualClass = ptrClass;
   } else if ((ptrClass != actualClass) && !actualClass->InheritsFrom(ptrClass)) {
      Error("Branch", "The actual class (%s) of the object provided for the definition of the branch \"%s\" does not inherit from %s", actualClass->GetName(), branchname, ptrClass->GetName());
      return 0;
   }
   if (actualClass && actualClass->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(actualClass->GetCollectionProxy())) {
      Error("Branch", writeStlWithoutProxyMsg,
            actualClass->GetName(), branchname, actualClass->GetName());
      return 0;
   }
   return BronchExec(branchname, actualClass->GetName(), (void*) addobj, kFALSE, bufsize, splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper to turn Branch call with an std::array into the relevant leaf list
// call
TBranch *TTree::BranchImpArr(const char *branchname, EDataType datatype, std::size_t N, void *addobj, Int_t bufsize,
                             Int_t /* splitlevel */)
{
   if (datatype == kOther_t || datatype == kNoType_t) {
      Error("Branch",
            "The inner type of the std::array passed specified for %s is not of a class or type known to ROOT",
            branchname);
   } else {
      TString varname;
      varname.Form("%s[%d]/%c", branchname, (int)N, DataTypeToChar(datatype));
      return Branch(branchname, addobj, varname.Data(), bufsize);
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Deprecated function. Use next function instead.

Int_t TTree::Branch(TList* li, Int_t bufsize /* = 32000 */ , Int_t splitlevel /* = 99 */)
{
   return Branch((TCollection*) li, bufsize, splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Create one branch for each element in the collection.
///
/// Each entry in the collection becomes a top level branch if the
/// corresponding class is not a collection. If it is a collection, the entry
/// in the collection becomes in turn top level branches, etc.
/// The splitlevel is decreased by 1 every time a new collection is found.
/// For example if list is a TObjArray*
///   - if splitlevel = 1, one top level branch is created for each element
///      of the TObjArray.
///   - if splitlevel = 2, one top level branch is created for each array element.
///     if, in turn, one of the array elements is a TCollection, one top level
///     branch will be created for each element of this collection.
///
/// In case a collection element is a TClonesArray, the special Tree constructor
/// for TClonesArray is called.
/// The collection itself cannot be a TClonesArray.
///
/// The function returns the total number of branches created.
///
/// If name is given, all branch names will be prefixed with name_.
///
/// IMPORTANT NOTE1: This function should not be called with splitlevel < 1.
///
/// IMPORTANT NOTE2: The branches created by this function will have names
/// corresponding to the collection or object names. It is important
/// to give names to collections to avoid misleading branch names or
/// identical branch names. By default collections have a name equal to
/// the corresponding class name, e.g. the default name for a TList is "TList".
///
/// And in general, in case two or more master branches contain subbranches
/// with identical names, one must add a "." (dot) character at the end
/// of the master branch name. This will force the name of the subbranches
/// to be of the form `master.subbranch` instead of simply `subbranch`.
/// This situation happens when the top level object
/// has two or more members referencing the same class.
/// For example, if a Tree has two branches B1 and B2 corresponding
/// to objects of the same class MyClass, one can do:
/// ~~~ {.cpp}
///     tree.Branch("B1.","MyClass",&b1,8000,1);
///     tree.Branch("B2.","MyClass",&b2,8000,1);
/// ~~~
/// if MyClass has 3 members a,b,c, the two instructions above will generate
/// subbranches called B1.a, B1.b ,B1.c, B2.a, B2.b, B2.c
///
/// Example:
/// ~~~ {.cpp}
///     {
///           TTree T("T","test list");
///           TList *list = new TList();
///
///           TObjArray *a1 = new TObjArray();
///           a1->SetName("a1");
///           list->Add(a1);
///           TH1F *ha1a = new TH1F("ha1a","ha1",100,0,1);
///           TH1F *ha1b = new TH1F("ha1b","ha1",100,0,1);
///           a1->Add(ha1a);
///           a1->Add(ha1b);
///           TObjArray *b1 = new TObjArray();
///           b1->SetName("b1");
///           list->Add(b1);
///           TH1F *hb1a = new TH1F("hb1a","hb1",100,0,1);
///           TH1F *hb1b = new TH1F("hb1b","hb1",100,0,1);
///           b1->Add(hb1a);
///           b1->Add(hb1b);
///
///           TObjArray *a2 = new TObjArray();
///           a2->SetName("a2");
///           list->Add(a2);
///           TH1S *ha2a = new TH1S("ha2a","ha2",100,0,1);
///           TH1S *ha2b = new TH1S("ha2b","ha2",100,0,1);
///           a2->Add(ha2a);
///           a2->Add(ha2b);
///
///           T.Branch(list,16000,2);
///           T.Print();
///     }
/// ~~~

Int_t TTree::Branch(TCollection* li, Int_t bufsize /* = 32000 */, Int_t splitlevel /* = 99 */, const char* name /* = "" */)
{

   if (!li) {
      return 0;
   }
   TObject* obj = 0;
   Int_t nbranches = GetListOfBranches()->GetEntries();
   if (li->InheritsFrom(TClonesArray::Class())) {
      Error("Branch", "Cannot call this constructor for a TClonesArray");
      return 0;
   }
   Int_t nch = strlen(name);
   TString branchname;
   TIter next(li);
   while ((obj = next())) {
      if ((splitlevel > 1) &&  obj->InheritsFrom(TCollection::Class()) && !obj->InheritsFrom(TClonesArray::Class())) {
         TCollection* col = (TCollection*) obj;
         if (nch) {
            branchname.Form("%s_%s_", name, col->GetName());
         } else {
            branchname.Form("%s_", col->GetName());
         }
         Branch(col, bufsize, splitlevel - 1, branchname);
      } else {
         if (nch && (name[nch-1] == '_')) {
            branchname.Form("%s%s", name, obj->GetName());
         } else {
            if (nch) {
               branchname.Form("%s_%s", name, obj->GetName());
            } else {
               branchname.Form("%s", obj->GetName());
            }
         }
         if (splitlevel > 99) {
            branchname += ".";
         }
         Bronch(branchname, obj->ClassName(), li->GetObjectRef(obj), bufsize, splitlevel - 1);
      }
   }
   return GetListOfBranches()->GetEntries() - nbranches;
}

////////////////////////////////////////////////////////////////////////////////
/// Create one branch for each element in the folder.
/// Returns the total number of branches created.

Int_t TTree::Branch(const char* foldername, Int_t bufsize /* = 32000 */, Int_t splitlevel /* = 99 */)
{
   TObject* ob = gROOT->FindObjectAny(foldername);
   if (!ob) {
      return 0;
   }
   if (ob->IsA() != TFolder::Class()) {
      return 0;
   }
   Int_t nbranches = GetListOfBranches()->GetEntries();
   TFolder* folder = (TFolder*) ob;
   TIter next(folder->GetListOfFolders());
   TObject* obj = 0;
   char* curname = new char[1000];
   char occur[20];
   while ((obj = next())) {
      snprintf(curname,1000, "%s/%s", foldername, obj->GetName());
      if (obj->IsA() == TFolder::Class()) {
         Branch(curname, bufsize, splitlevel - 1);
      } else {
         void* add = (void*) folder->GetListOfFolders()->GetObjectRef(obj);
         for (Int_t i = 0; i < 1000; ++i) {
            if (curname[i] == 0) {
               break;
            }
            if (curname[i] == '/') {
               curname[i] = '.';
            }
         }
         Int_t noccur = folder->Occurence(obj);
         if (noccur > 0) {
            snprintf(occur,20, "_%d", noccur);
            strlcat(curname, occur,1000);
         }
         TBranchElement* br = (TBranchElement*) Bronch(curname, obj->ClassName(), add, bufsize, splitlevel - 1);
         if (br) br->SetBranchFolder();
      }
   }
   delete[] curname;
   return GetListOfBranches()->GetEntries() - nbranches;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new TTree Branch.
///
/// This Branch constructor is provided to support non-objects in
/// a Tree. The variables described in leaflist may be simple
/// variables or structures.  // See the two following
/// constructors for writing objects in a Tree.
///
/// By default the branch buffers are stored in the same file as the Tree.
/// use TBranch::SetFile to specify a different file
///
///    * address is the address of the first item of a structure.
///    * leaflist is the concatenation of all the variable names and types
///      separated by a colon character :
///      The variable name and the variable type are separated by a slash (/).
///      The variable type may be 0,1 or 2 characters. If no type is given,
///      the type of the variable is assumed to be the same as the previous
///      variable. If the first variable does not have a type, it is assumed
///      of type F by default. The list of currently supported types is given below:
///         - `C` : a character string terminated by the 0 character
///         - `B` : an 8 bit signed integer (`Char_t`)
///         - `b` : an 8 bit unsigned integer (`UChar_t`)
///         - `S` : a 16 bit signed integer (`Short_t`)
///         - `s` : a 16 bit unsigned integer (`UShort_t`)
///         - `I` : a 32 bit signed integer (`Int_t`)
///         - `i` : a 32 bit unsigned integer (`UInt_t`)
///         - `F` : a 32 bit floating point (`Float_t`)
///         - `f` : a 24 bit floating point with truncated mantissa (`Float16_t`)
///         - `D` : a 64 bit floating point (`Double_t`)
///         - `d` : a 24 bit truncated floating point (`Double32_t`)
///         - `L` : a 64 bit signed integer (`Long64_t`)
///         - `l` : a 64 bit unsigned integer (`ULong64_t`)
///         - `G` : a long signed integer, stored as 64 bit (`Long_t`)
///         - `g` : a long unsigned integer, stored as 64 bit (`ULong_t`)
///         - `O` : [the letter `o`, not a zero] a boolean (`Bool_t`)
///
///      Arrays of values are supported with the following syntax:
///         - If leaf name has the form var[nelem], where nelem is alphanumeric, then
///           if nelem is a leaf name, it is used as the variable size of the array,
///           otherwise return 0.
///         - If leaf name has the form var[nelem], where nelem is a non-negative integer, then
///           it is used as the fixed size of the array.
///         - If leaf name has the form of a multi-dimensional array (e.g. var[nelem][nelem2])
///           where nelem and nelem2 are non-negative integer) then
///           it is used as a 2 dimensional array of fixed size.
///         - In case of the truncated floating point types (Float16_t and Double32_t) you can
///           furthermore specify the range in the style [xmin,xmax] or [xmin,xmax,nbits] after
///           the type character. See `TStreamerElement::GetRange()` for further information.
///
///      Any of other form is not supported.
///
/// Note that the TTree will assume that all the item are contiguous in memory.
/// On some platform, this is not always true of the member of a struct or a class,
/// due to padding and alignment.  Sorting your data member in order of decreasing
/// sizeof usually leads to their being contiguous in memory.
///
///    * bufsize is the buffer size in bytes for this branch
///      The default value is 32000 bytes and should be ok for most cases.
///      You can specify a larger value (e.g. 256000) if your Tree is not split
///      and each entry is large (Megabytes)
///      A small value for bufsize is optimum if you intend to access
///      the entries in the Tree randomly and your Tree is in split mode.

TBranch* TTree::Branch(const char* name, void* address, const char* leaflist, Int_t bufsize /* = 32000 */)
{
   TBranch* branch = new TBranch(this, name, address, leaflist, bufsize);
   if (branch->IsZombie()) {
      delete branch;
      branch = 0;
      return 0;
   }
   fBranches.Add(branch);
   return branch;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new branch with the object of class classname at address addobj.
///
/// WARNING:
///
/// Starting with Root version 3.01, the Branch function uses the new style
/// branches (TBranchElement). To get the old behaviour, you can:
///   - call BranchOld or
///   - call TTree::SetBranchStyle(0)
///
/// Note that with the new style, classname does not need to derive from TObject.
/// It must derived from TObject if the branch style has been set to 0 (old)
///
/// Note: See the comments in TBranchElement::SetAddress() for a more
///       detailed discussion of the meaning of the addobj parameter in
///       the case of new-style branches.
///
/// Use splitlevel < 0 instead of splitlevel=0 when the class
/// has a custom Streamer
///
/// Note: if the split level is set to the default (99),  TTree::Branch will
/// not issue a warning if the class can not be split.

TBranch* TTree::Branch(const char* name, const char* classname, void* addobj, Int_t bufsize /* = 32000 */, Int_t splitlevel /* = 99 */)
{
   if (fgBranchStyle == 1) {
      return Bronch(name, classname, addobj, bufsize, splitlevel);
   } else {
      if (splitlevel < 0) {
         splitlevel = 0;
      }
      return BranchOld(name, classname, addobj, bufsize, splitlevel);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new TTree BranchObject.
///
/// Build a TBranchObject for an object of class classname.
/// addobj is the address of a pointer to an object of class classname.
/// IMPORTANT: classname must derive from TObject.
/// The class dictionary must be available (ClassDef in class header).
///
/// This option requires access to the library where the corresponding class
/// is defined. Accessing one single data member in the object implies
/// reading the full object.
/// See the next Branch constructor for a more efficient storage
/// in case the entry consists of arrays of identical objects.
///
/// By default the branch buffers are stored in the same file as the Tree.
/// use TBranch::SetFile to specify a different file
///
/// IMPORTANT NOTE about branch names:
///
/// And in general, in case two or more master branches contain subbranches
/// with identical names, one must add a "." (dot) character at the end
/// of the master branch name. This will force the name of the subbranches
/// to be of the form `master.subbranch` instead of simply `subbranch`.
/// This situation happens when the top level object
/// has two or more members referencing the same class.
/// For example, if a Tree has two branches B1 and B2 corresponding
/// to objects of the same class MyClass, one can do:
/// ~~~ {.cpp}
///     tree.Branch("B1.","MyClass",&b1,8000,1);
///     tree.Branch("B2.","MyClass",&b2,8000,1);
/// ~~~
/// if MyClass has 3 members a,b,c, the two instructions above will generate
/// subbranches called B1.a, B1.b ,B1.c, B2.a, B2.b, B2.c
///
/// bufsize is the buffer size in bytes for this branch
/// The default value is 32000 bytes and should be ok for most cases.
/// You can specify a larger value (e.g. 256000) if your Tree is not split
/// and each entry is large (Megabytes)
/// A small value for bufsize is optimum if you intend to access
/// the entries in the Tree randomly and your Tree is in split mode.

TBranch* TTree::BranchOld(const char* name, const char* classname, void* addobj, Int_t bufsize /* = 32000 */, Int_t splitlevel /* = 1 */)
{
   TClass* cl = TClass::GetClass(classname);
   if (!cl) {
      Error("BranchOld", "Cannot find class: '%s'", classname);
      return 0;
   }
   if (!cl->IsTObject()) {
      if (fgBranchStyle == 0) {
        Fatal("BranchOld", "The requested class ('%s') does not inherit from TObject.\n"
              "\tfgBranchStyle is set to zero requesting by default to use BranchOld.\n"
              "\tIf this is intentional use Bronch instead of Branch or BranchOld.", classname);
      } else {
        Fatal("BranchOld", "The requested class ('%s') does not inherit from TObject.\n"
              "\tYou can not use BranchOld to store objects of this type.",classname);
      }
      return 0;
   }
   TBranch* branch = new TBranchObject(this, name, classname, addobj, bufsize, splitlevel);
   fBranches.Add(branch);
   if (!splitlevel) {
      return branch;
   }
   // We are going to fully split the class now.
   TObjArray* blist = branch->GetListOfBranches();
   const char* rdname = 0;
   const char* dname = 0;
   TString branchname;
   char** apointer = (char**) addobj;
   TObject* obj = (TObject*) *apointer;
   Bool_t delobj = kFALSE;
   if (!obj) {
      obj = (TObject*) cl->New();
      delobj = kTRUE;
   }
   // Build the StreamerInfo if first time for the class.
   BuildStreamerInfo(cl, obj);
   // Loop on all public data members of the class and its base classes.
   Int_t lenName = strlen(name);
   Int_t isDot = 0;
   if (name[lenName-1] == '.') {
      isDot = 1;
   }
   TBranch* branch1 = 0;
   TRealData* rd = 0;
   TRealData* rdi = 0;
   TIter nexti(cl->GetListOfRealData());
   TIter next(cl->GetListOfRealData());
   // Note: This loop results in a full split because the
   //       real data list includes all data members of
   //       data members.
   while ((rd = (TRealData*) next())) {
      if (rd->TestBit(TRealData::kTransient)) continue;

      // Loop over all data members creating branches for each one.
      TDataMember* dm = rd->GetDataMember();
      if (!dm->IsPersistent()) {
         // Do not process members with an "!" as the first character in the comment field.
         continue;
      }
      if (rd->IsObject()) {
         // We skip data members of class type.
         // But we do build their real data, their
         // streamer info, and write their streamer
         // info to the current directory's file.
         // Oh yes, and we also do this for all of
         // their base classes.
         TClass* clm = TClass::GetClass(dm->GetFullTypeName());
         if (clm) {
            BuildStreamerInfo(clm, (char*) obj + rd->GetThisOffset());
         }
         continue;
      }
      rdname = rd->GetName();
      dname = dm->GetName();
      if (cl->CanIgnoreTObjectStreamer()) {
         // Skip the TObject base class data members.
         // FIXME: This prevents a user from ever
         //        using these names themself!
         if (!strcmp(dname, "fBits")) {
            continue;
         }
         if (!strcmp(dname, "fUniqueID")) {
            continue;
         }
      }
      TDataType* dtype = dm->GetDataType();
      Int_t code = 0;
      if (dtype) {
         code = dm->GetDataType()->GetType();
      }
      // Encode branch name. Use real data member name
      branchname = rdname;
      if (isDot) {
         if (dm->IsaPointer()) {
            // FIXME: This is wrong!  The asterisk is not usually in the front!
            branchname.Form("%s%s", name, &rdname[1]);
         } else {
            branchname.Form("%s%s", name, &rdname[0]);
         }
      }
      // FIXME: Change this to a string stream.
      TString leaflist;
      Int_t offset = rd->GetThisOffset();
      char* pointer = ((char*) obj) + offset;
      if (dm->IsaPointer()) {
         // We have a pointer to an object or a pointer to an array of basic types.
         TClass* clobj = 0;
         if (!dm->IsBasic()) {
            clobj = TClass::GetClass(dm->GetTypeName());
         }
         if (clobj && clobj->InheritsFrom(TClonesArray::Class())) {
            // We have a pointer to a clones array.
            char* cpointer = (char*) pointer;
            char** ppointer = (char**) cpointer;
            TClonesArray* li = (TClonesArray*) *ppointer;
            if (splitlevel != 2) {
               if (isDot) {
                  branch1 = new TBranchClones(branch,branchname, pointer, bufsize);
               } else {
                  // FIXME: This is wrong!  The asterisk is not usually in the front!
                  branch1 = new TBranchClones(branch,&branchname.Data()[1], pointer, bufsize);
               }
               blist->Add(branch1);
            } else {
               if (isDot) {
                  branch1 = new TBranchObject(branch, branchname, li->ClassName(), pointer, bufsize);
               } else {
                  // FIXME: This is wrong!  The asterisk is not usually in the front!
                  branch1 = new TBranchObject(branch, &branchname.Data()[1], li->ClassName(), pointer, bufsize);
               }
               blist->Add(branch1);
            }
         } else if (clobj) {
            // We have a pointer to an object.
            //
            // It must be a TObject object.
            if (!clobj->IsTObject()) {
               continue;
            }
            branch1 = new TBranchObject(branch, dname, clobj->GetName(), pointer, bufsize, 0);
            if (isDot) {
               branch1->SetName(branchname);
            } else {
               // FIXME: This is wrong!  The asterisk is not usually in the front!
               // Do not use the first character (*).
               branch1->SetName(&branchname.Data()[1]);
            }
            blist->Add(branch1);
         } else {
            // We have a pointer to an array of basic types.
            //
            // Check the comments in the text of the code for an index specification.
            const char* index = dm->GetArrayIndex();
            if (index[0]) {
               // We are a pointer to a varying length array of basic types.
               //check that index is a valid data member name
               //if member is part of an object (e.g. fA and index=fN)
               //index must be changed from fN to fA.fN
               TString aindex (rd->GetName());
               Ssiz_t rdot = aindex.Last('.');
               if (rdot>=0) {
                  aindex.Remove(rdot+1);
                  aindex.Append(index);
               }
               nexti.Reset();
               while ((rdi = (TRealData*) nexti())) {
                  if (rdi->TestBit(TRealData::kTransient)) continue;

                  if (!strcmp(rdi->GetName(), index)) {
                     break;
                  }
                  if (!strcmp(rdi->GetName(), aindex)) {
                     index = rdi->GetName();
                     break;
                  }
               }

               char vcode = DataTypeToChar((EDataType)code);
               // Note that we differentiate between strings and
               // char array by the fact that there is NO specified
               // size for a string (see next if (code == 1)

               if (vcode) {
                  leaflist.Form("%s[%s]/%c", &rdname[0], index, vcode);
               } else {
                  Error("BranchOld", "Cannot create branch for rdname: %s code: %d", branchname.Data(), code);
                  leaflist = "";
               }
            } else {
               // We are possibly a character string.
               if (code == 1) {
                  // We are a character string.
                  leaflist.Form("%s/%s", dname, "C");
               } else {
                  // Invalid array specification.
                  // FIXME: We need an error message here.
                  continue;
               }
            }
            // There are '*' in both the branchname and leaflist, remove them.
            TString bname( branchname );
            bname.ReplaceAll("*","");
            leaflist.ReplaceAll("*","");
            // Add the branch to the tree and indicate that the address
            // is that of a pointer to be dereferenced before using.
            branch1 = new TBranch(branch, bname, *((void**) pointer), leaflist, bufsize);
            TLeaf* leaf = (TLeaf*) branch1->GetListOfLeaves()->At(0);
            leaf->SetBit(TLeaf::kIndirectAddress);
            leaf->SetAddress((void**) pointer);
            blist->Add(branch1);
         }
      } else if (dm->IsBasic()) {
         // We have a basic type.

         char vcode = DataTypeToChar((EDataType)code);
         if (vcode) {
            leaflist.Form("%s/%c", rdname, vcode);
         } else {
            Error("BranchOld", "Cannot create branch for rdname: %s code: %d", branchname.Data(), code);
            leaflist = "";
         }
         branch1 = new TBranch(branch, branchname, pointer, leaflist, bufsize);
         branch1->SetTitle(rdname);
         blist->Add(branch1);
      } else {
         // We have a class type.
         // Note: This cannot happen due to the rd->IsObject() test above.
         // FIXME: Put an error message here just in case.
      }
      if (branch1) {
         branch1->SetOffset(offset);
      } else {
         Warning("BranchOld", "Cannot process member: '%s'", rdname);
      }
   }
   if (delobj) {
      delete obj;
      obj = 0;
   }
   return branch;
}

////////////////////////////////////////////////////////////////////////////////
/// Build the optional branch supporting the TRefTable.
/// This branch will keep all the information to find the branches
/// containing referenced objects.
///
/// At each Tree::Fill, the branch numbers containing the
/// referenced objects are saved to the TBranchRef basket.
/// When the Tree header is saved (via TTree::Write), the branch
/// is saved keeping the information with the pointers to the branches
/// having referenced objects.

TBranch* TTree::BranchRef()
{
   if (!fBranchRef) {
      fBranchRef = new TBranchRef(this);
   }
   return fBranchRef;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new TTree BranchElement.
///
/// ## WARNING about this new function
///
/// This function is designed to replace the internal
/// implementation of the old TTree::Branch (whose implementation
/// has been moved to BranchOld).
///
/// NOTE: The 'Bronch' method supports only one possible calls
/// signature (where the object type has to be specified
/// explicitly and the address must be the address of a pointer).
/// For more flexibility use 'Branch'.  Use Bronch only in (rare)
/// cases (likely to be legacy cases) where both the new and old
/// implementation of Branch needs to be used at the same time.
///
/// This function is far more powerful than the old Branch
/// function.  It supports the full C++, including STL and has
/// the same behaviour in split or non-split mode. classname does
/// not have to derive from TObject.  The function is based on
/// the new TStreamerInfo.
///
/// Build a TBranchElement for an object of class classname.
///
/// addr is the address of a pointer to an object of class
/// classname.  The class dictionary must be available (ClassDef
/// in class header).
///
/// Note: See the comments in TBranchElement::SetAddress() for a more
///       detailed discussion of the meaning of the addr parameter.
///
/// This option requires access to the library where the
/// corresponding class is defined. Accessing one single data
/// member in the object implies reading the full object.
///
/// By default the branch buffers are stored in the same file as the Tree.
/// use TBranch::SetFile to specify a different file
///
/// IMPORTANT NOTE about branch names:
///
/// And in general, in case two or more master branches contain subbranches
/// with identical names, one must add a "." (dot) character at the end
/// of the master branch name. This will force the name of the subbranches
/// to be of the form `master.subbranch` instead of simply `subbranch`.
/// This situation happens when the top level object
/// has two or more members referencing the same class.
/// For example, if a Tree has two branches B1 and B2 corresponding
/// to objects of the same class MyClass, one can do:
/// ~~~ {.cpp}
///     tree.Branch("B1.","MyClass",&b1,8000,1);
///     tree.Branch("B2.","MyClass",&b2,8000,1);
/// ~~~
/// if MyClass has 3 members a,b,c, the two instructions above will generate
/// subbranches called B1.a, B1.b ,B1.c, B2.a, B2.b, B2.c
///
/// bufsize is the buffer size in bytes for this branch
/// The default value is 32000 bytes and should be ok for most cases.
/// You can specify a larger value (e.g. 256000) if your Tree is not split
/// and each entry is large (Megabytes)
/// A small value for bufsize is optimum if you intend to access
/// the entries in the Tree randomly and your Tree is in split mode.
///
/// Use splitlevel < 0 instead of splitlevel=0 when the class
/// has a custom Streamer
///
/// Note: if the split level is set to the default (99),  TTree::Branch will
/// not issue a warning if the class can not be split.

TBranch* TTree::Bronch(const char* name, const char* classname, void* addr, Int_t bufsize /* = 32000 */, Int_t splitlevel /* = 99 */)
{
   return BronchExec(name, classname, addr, kTRUE, bufsize, splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function implementing TTree::Bronch and TTree::Branch(const char *name, T &obj);

TBranch* TTree::BronchExec(const char* name, const char* classname, void* addr, Bool_t isptrptr, Int_t bufsize /* = 32000 */, Int_t splitlevel /* = 99 */)
{
   TClass* cl = TClass::GetClass(classname);
   if (!cl) {
      Error("Bronch", "Cannot find class:%s", classname);
      return 0;
   }

   //if splitlevel <= 0 and class has a custom Streamer, we must create
   //a TBranchObject. We cannot assume that TClass::ReadBuffer is consistent
   //with the custom Streamer. The penalty is that one cannot process
   //this Tree without the class library containing the class.

   char* objptr = 0;
   if (!isptrptr) {
      objptr = (char*)addr;
   } else if (addr) {
      objptr = *((char**) addr);
   }

   if (cl == TClonesArray::Class()) {
      TClonesArray* clones = (TClonesArray*) objptr;
      if (!clones) {
         Error("Bronch", "Pointer to TClonesArray is null");
         return 0;
      }
      if (!clones->GetClass()) {
         Error("Bronch", "TClonesArray with no class defined in branch: %s", name);
         return 0;
      }
      if (!clones->GetClass()->HasDataMemberInfo()) {
         Error("Bronch", "TClonesArray with no dictionary defined in branch: %s", name);
         return 0;
      }
      bool hasCustomStreamer = clones->GetClass()->TestBit(TClass::kHasCustomStreamerMember);
      if (splitlevel > 0) {
         if (hasCustomStreamer)
            Warning("Bronch", "Using split mode on a class: %s with a custom Streamer", clones->GetClass()->GetName());
      } else {
         if (hasCustomStreamer) clones->BypassStreamer(kFALSE);
         TBranchObject *branch = new TBranchObject(this,name,classname,addr,bufsize,0,/*compress=*/ -1,isptrptr);
         fBranches.Add(branch);
         return branch;
      }
   }

   if (cl->GetCollectionProxy()) {
      TVirtualCollectionProxy* collProxy = cl->GetCollectionProxy();
      //if (!collProxy) {
      //   Error("Bronch", "%s is missing its CollectionProxy (for branch %s)", classname, name);
      //}
      TClass* inklass = collProxy->GetValueClass();
      if (!inklass && (collProxy->GetType() == 0)) {
         Error("Bronch", "%s with no class defined in branch: %s", classname, name);
         return 0;
      }
      if ((splitlevel > 0) && inklass && (inklass->GetCollectionProxy() == 0)) {
         ROOT::ESTLType stl = cl->GetCollectionType();
         if ((stl != ROOT::kSTLmap) && (stl != ROOT::kSTLmultimap)) {
            if (!inklass->HasDataMemberInfo()) {
               Error("Bronch", "Container with no dictionary defined in branch: %s", name);
               return 0;
            }
            if (inklass->TestBit(TClass::kHasCustomStreamerMember)) {
               Warning("Bronch", "Using split mode on a class: %s with a custom Streamer", inklass->GetName());
            }
         }
      }
      //-------------------------------------------------------------------------
      // If the splitting switch is enabled, the split level is big enough and
      // the collection contains pointers we can split it
      //////////////////////////////////////////////////////////////////////////

      TBranch *branch;
      if( splitlevel > kSplitCollectionOfPointers && collProxy->HasPointers() )
         branch = new TBranchSTL( this, name, collProxy, bufsize, splitlevel );
      else
         branch = new TBranchElement(this, name, collProxy, bufsize, splitlevel);
      fBranches.Add(branch);
      if (isptrptr) {
         branch->SetAddress(addr);
      } else {
         branch->SetObject(addr);
      }
      return branch;
   }

   Bool_t hasCustomStreamer = kFALSE;
   if (!cl->HasDataMemberInfo() && !cl->GetCollectionProxy()) {
      Error("Bronch", "Cannot find dictionary for class: %s", classname);
      return 0;
   }

   if (!cl->GetCollectionProxy() && cl->TestBit(TClass::kHasCustomStreamerMember)) {
      // Not an STL container and the linkdef file had a "-" after the class name.
      hasCustomStreamer = kTRUE;
   }

   if (splitlevel < 0 || ((splitlevel == 0) && hasCustomStreamer && cl->IsTObject())) {
      TBranchObject* branch = new TBranchObject(this, name, classname, addr, bufsize, 0, /*compress=*/ ROOT::RCompressionSetting::EAlgorithm::kInherit, isptrptr);
      fBranches.Add(branch);
      return branch;
   }

   if (cl == TClonesArray::Class()) {
      // Special case of TClonesArray.
      // No dummy object is created.
      // The streamer info is not rebuilt unoptimized.
      // No dummy top-level branch is created.
      // No splitting is attempted.
      TBranchElement* branch = new TBranchElement(this, name, (TClonesArray*) objptr, bufsize, splitlevel%kSplitCollectionOfPointers);
      fBranches.Add(branch);
      if (isptrptr) {
         branch->SetAddress(addr);
      } else {
         branch->SetObject(addr);
      }
      return branch;
   }

   //
   // If we are not given an object to use as an i/o buffer
   // then create a temporary one which we will delete just
   // before returning.
   //

   Bool_t delobj = kFALSE;

   if (!objptr) {
      objptr = (char*) cl->New();
      delobj = kTRUE;
   }

   //
   // Avoid splitting unsplittable classes.
   //

   if ((splitlevel > 0) && !cl->CanSplit()) {
      if (splitlevel != 99) {
         Warning("Bronch", "%s cannot be split, resetting splitlevel to 0", cl->GetName());
      }
      splitlevel = 0;
   }

   //
   // Make sure the streamer info is built and fetch it.
   //
   // If we are splitting, then make sure the streamer info
   // is built unoptimized (data members are not combined).
   //

   TStreamerInfo* sinfo = BuildStreamerInfo(cl, objptr, splitlevel==0);
   if (!sinfo) {
      Error("Bronch", "Cannot build the StreamerInfo for class: %s", cl->GetName());
      return 0;
   }

   //
   // Create a dummy top level branch object.
   //

   Int_t id = -1;
   if (splitlevel > 0) {
      id = -2;
   }
   TBranchElement* branch = new TBranchElement(this, name, sinfo, id, objptr, bufsize, splitlevel);
   fBranches.Add(branch);

   //
   // Do splitting, if requested.
   //

   if (splitlevel%kSplitCollectionOfPointers > 0) {
      branch->Unroll(name, cl, sinfo, objptr, bufsize, splitlevel);
   }

   //
   // Setup our offsets into the user's i/o buffer.
   //

   if (isptrptr) {
      branch->SetAddress(addr);
   } else {
      branch->SetObject(addr);
   }

   if (delobj) {
      cl->Destructor(objptr);
      objptr = 0;
   }

   return branch;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse content of the TTree.

void TTree::Browse(TBrowser* b)
{
   fBranches.Browse(b);
   if (fUserInfo) {
      if (strcmp("TList",fUserInfo->GetName())==0) {
         fUserInfo->SetName("UserInfo");
         b->Add(fUserInfo);
         fUserInfo->SetName("TList");
      } else {
         b->Add(fUserInfo);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Build a Tree Index (default is TTreeIndex).
/// See a description of the parameters and functionality in
/// TTreeIndex::TTreeIndex().
///
/// The return value is the number of entries in the Index (< 0 indicates failure).
///
/// A TTreeIndex object pointed by fTreeIndex is created.
/// This object will be automatically deleted by the TTree destructor.
/// If an index is already existing, this is replaced by the new one without being
/// deleted. This behaviour prevents the deletion of a previously external index
/// assigned to the TTree via the TTree::SetTreeIndex() method.
/// See also comments in TTree::SetTreeIndex().

Int_t TTree::BuildIndex(const char* majorname, const char* minorname /* = "0" */)
{
   fTreeIndex = GetPlayer()->BuildIndex(this, majorname, minorname);
   if (fTreeIndex->IsZombie()) {
      delete fTreeIndex;
      fTreeIndex = 0;
      return 0;
   }
   return fTreeIndex->GetN();
}

////////////////////////////////////////////////////////////////////////////////
/// Build StreamerInfo for class cl.
/// pointer is an optional argument that may contain a pointer to an object of cl.

TStreamerInfo* TTree::BuildStreamerInfo(TClass* cl, void* pointer /* = 0 */, Bool_t canOptimize /* = kTRUE */ )
{
   if (!cl) {
      return 0;
   }
   cl->BuildRealData(pointer);
   TStreamerInfo* sinfo = (TStreamerInfo*)cl->GetStreamerInfo(cl->GetClassVersion());

   // Create StreamerInfo for all base classes.
   TBaseClass* base = 0;
   TIter nextb(cl->GetListOfBases());
   while((base = (TBaseClass*) nextb())) {
      if (base->IsSTLContainer()) {
         continue;
      }
      TClass* clm = TClass::GetClass(base->GetName());
      BuildStreamerInfo(clm, pointer, canOptimize);
   }
   if (sinfo && fDirectory) {
      sinfo->ForceWriteInfo(fDirectory->GetFile());
   }
   return sinfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Called by TTree::Fill() when file has reached its maximum fgMaxTreeSize.
/// Create a new file. If the original file is named "myfile.root",
/// subsequent files are named "myfile_1.root", "myfile_2.root", etc.
///
/// Returns a pointer to the new file.
///
/// Currently, the automatic change of file is restricted
/// to the case where the tree is in the top level directory.
/// The file should not contain sub-directories.
///
/// Before switching to a new file, the tree header is written
/// to the current file, then the current file is closed.
///
/// To process the multiple files created by ChangeFile, one must use
/// a TChain.
///
/// The new file name has a suffix "_N" where N is equal to fFileNumber+1.
/// By default a Root session starts with fFileNumber=0. One can set
/// fFileNumber to a different value via TTree::SetFileNumber.
/// In case a file named "_N" already exists, the function will try
/// a file named "__N", then "___N", etc.
///
/// fgMaxTreeSize can be set via the static function TTree::SetMaxTreeSize.
/// The default value of fgMaxTreeSize is 100 Gigabytes.
///
/// If the current file contains other objects like TH1 and TTree,
/// these objects are automatically moved to the new file.
///
/// \warning Be careful when writing the final Tree header to the file!
///      Don't do:
/// ~~~ {.cpp}
///     TFile *file = new TFile("myfile.root","recreate");
///     TTree *T = new TTree("T","title");
///     T->Fill(); // Loop
///     file->Write();
///     file->Close();
/// ~~~
/// \warning but do the following:
/// ~~~ {.cpp}
///     TFile *file = new TFile("myfile.root","recreate");
///     TTree *T = new TTree("T","title");
///     T->Fill(); // Loop
///     file = T->GetCurrentFile(); // To get the pointer to the current file
///     file->Write();
///     file->Close();
/// ~~~
///
/// \note This method is never called if the input file is a `TMemFile` or derivate.

TFile* TTree::ChangeFile(TFile* file)
{
   file->cd();
   Write();
   Reset();
   constexpr auto kBufSize = 2000;
   char* fname = new char[kBufSize];
   ++fFileNumber;
   char uscore[10];
   for (Int_t i = 0; i < 10; ++i) {
      uscore[i] = 0;
   }
   Int_t nus = 0;
   // Try to find a suitable file name that does not already exist.
   while (nus < 10) {
      uscore[nus] = '_';
      fname[0] = 0;
      strlcpy(fname, file->GetName(), kBufSize);

      if (fFileNumber > 1) {
         char* cunder = strrchr(fname, '_');
         if (cunder) {
            snprintf(cunder, kBufSize - Int_t(cunder - fname), "%s%d", uscore, fFileNumber);
            const char* cdot = strrchr(file->GetName(), '.');
            if (cdot) {
               strlcat(fname, cdot, kBufSize);
            }
         } else {
            char fcount[21];
            snprintf(fcount,21, "%s%d", uscore, fFileNumber);
            strlcat(fname, fcount, kBufSize);
         }
      } else {
         char* cdot = strrchr(fname, '.');
         if (cdot) {
            snprintf(cdot, kBufSize - Int_t(fname-cdot), "%s%d", uscore, fFileNumber);
            strlcat(fname, strrchr(file->GetName(), '.'), kBufSize);
         } else {
            char fcount[21];
            snprintf(fcount,21, "%s%d", uscore, fFileNumber);
            strlcat(fname, fcount, kBufSize);
         }
      }
      if (gSystem->AccessPathName(fname)) {
         break;
      }
      ++nus;
      Warning("ChangeFile", "file %s already exist, trying with %d underscores", fname, nus+1);
   }
   Int_t compress = file->GetCompressionSettings();
   TFile* newfile = TFile::Open(fname, "recreate", "chain files", compress);
   if (newfile == 0) {
      Error("Fill","Failed to open new file %s, continuing as a memory tree.",fname);
   } else {
      Printf("Fill: Switching to new file: %s", fname);
   }
   // The current directory may contain histograms and trees.
   // These objects must be moved to the new file.
   TBranch* branch = 0;
   TObject* obj = 0;
   while ((obj = file->GetList()->First())) {
      file->Remove(obj);
      // Histogram: just change the directory.
      if (obj->InheritsFrom("TH1")) {
         gROOT->ProcessLine(TString::Format("((%s*)0x%lx)->SetDirectory((TDirectory*)0x%lx);", obj->ClassName(), (Long_t) obj, (Long_t) newfile));
         continue;
      }
      // Tree: must save all trees in the old file, reset them.
      if (obj->InheritsFrom(TTree::Class())) {
         TTree* t = (TTree*) obj;
         if (t != this) {
            t->AutoSave();
            t->Reset();
            t->fFileNumber = fFileNumber;
         }
         t->SetDirectory(newfile);
         TIter nextb(t->GetListOfBranches());
         while ((branch = (TBranch*)nextb())) {
            branch->SetFile(newfile);
         }
         if (t->GetBranchRef()) {
            t->GetBranchRef()->SetFile(newfile);
         }
         continue;
      }
      // Not a TH1 or a TTree, move object to new file.
      if (newfile) newfile->Append(obj);
      file->Remove(obj);
   }
   file->TObject::Delete();
   file = 0;
   delete[] fname;
   fname = 0;
   return newfile;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether or not the address described by the last 3 parameters
/// matches the content of the branch. If a Data Model Evolution conversion
/// is involved, reset the fInfo of the branch.
/// The return values are:
//
/// - kMissingBranch (-5) : Missing branch
/// - kInternalError (-4) : Internal error (could not find the type corresponding to a data type number)
/// - kMissingCompiledCollectionProxy (-3) : Missing compiled collection proxy for a compiled collection
/// - kMismatch (-2) : Non-Class Pointer type given does not match the type expected by the branch
/// - kClassMismatch (-1) : Class Pointer type given does not match the type expected by the branch
/// - kMatch (0) : perfect match
/// - kMatchConversion (1) : match with (I/O) conversion
/// - kMatchConversionCollection (2) : match with (I/O) conversion of the content of a collection
/// - kMakeClass (3) : MakeClass mode so we can not check.
/// - kVoidPtr (4) : void* passed so no check was made.
/// - kNoCheck (5) : Underlying TBranch not yet available so no check was made.
/// In addition this can be multiplexed with the two bits:
/// - kNeedEnableDecomposedObj : in order for the address (type) to be 'usable' the branch needs to be in Decomposed Object (aka MakeClass) mode.
/// - kNeedDisableDecomposedObj : in order for the address (type) to be 'usable' the branch needs to not be in Decomposed Object (aka MakeClass) mode.
/// This bits can be masked out by using kDecomposedObjMask

Int_t TTree::CheckBranchAddressType(TBranch* branch, TClass* ptrClass, EDataType datatype, Bool_t isptr)
{
   if (GetMakeClass()) {
      // If we are in MakeClass mode so we do not really use classes.
      return kMakeClass;
   }

   // Let's determine what we need!
   TClass* expectedClass = 0;
   EDataType expectedType = kOther_t;
   if (0 != branch->GetExpectedType(expectedClass,expectedType) ) {
      // Something went wrong, the warning message has already been issued.
      return kInternalError;
   }
   bool isBranchElement = branch->InheritsFrom( TBranchElement::Class() );
   if (expectedClass && datatype == kOther_t && ptrClass == 0) {
      if (isBranchElement) {
         TBranchElement* bEl = (TBranchElement*)branch;
         bEl->SetTargetClass( expectedClass->GetName() );
      }
      if (expectedClass && expectedClass->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(expectedClass->GetCollectionProxy())) {
         Error("SetBranchAddress", "Unable to determine the type given for the address for \"%s\". "
               "The class expected (%s) refers to an stl collection and do not have a compiled CollectionProxy.  "
               "Please generate the dictionary for this class (%s)",
               branch->GetName(), expectedClass->GetName(), expectedClass->GetName());
         return kMissingCompiledCollectionProxy;
      }
      if (!expectedClass->IsLoaded()) {
         // The originally expected class does not have a dictionary, it is then plausible that the pointer being passed is the right type
         // (we really don't know).  So let's express that.
         Error("SetBranchAddress", "Unable to determine the type given for the address for \"%s\". "
               "The class expected (%s) does not have a dictionary and needs to be emulated for I/O purposes but is being passed a compiled object."
               "Please generate the dictionary for this class (%s)",
               branch->GetName(), expectedClass->GetName(), expectedClass->GetName());
      } else {
         Error("SetBranchAddress", "Unable to determine the type given for the address for \"%s\". "
               "This is probably due to a missing dictionary, the original data class for this branch is %s.", branch->GetName(), expectedClass->GetName());
      }
      return kClassMismatch;
   }
   if (expectedClass && ptrClass && (branch->GetMother() == branch)) {
      // Top Level branch
      if (!isptr) {
         Error("SetBranchAddress", "The address for \"%s\" should be the address of a pointer!", branch->GetName());
      }
   }
   if (expectedType == kFloat16_t) {
      expectedType = kFloat_t;
   }
   if (expectedType == kDouble32_t) {
      expectedType = kDouble_t;
   }
   if (datatype == kFloat16_t) {
      datatype = kFloat_t;
   }
   if (datatype == kDouble32_t) {
      datatype = kDouble_t;
   }

   /////////////////////////////////////////////////////////////////////////////
   // Deal with the class renaming
   /////////////////////////////////////////////////////////////////////////////

   if( expectedClass && ptrClass &&
       expectedClass != ptrClass &&
       isBranchElement &&
       ptrClass->GetSchemaRules() &&
       ptrClass->GetSchemaRules()->HasRuleWithSourceClass( expectedClass->GetName() ) ) {
      TBranchElement* bEl = (TBranchElement*)branch;

      if ( ptrClass->GetCollectionProxy() && expectedClass->GetCollectionProxy() ) {
         if (gDebug > 7)
            Info("SetBranchAddress", "Matching STL collection (at least according to the SchemaRuleSet when "
               "reading a %s into a %s",expectedClass->GetName(),ptrClass->GetName());

         bEl->SetTargetClass( ptrClass->GetName() );
         return kMatchConversion;

      } else if ( !ptrClass->GetConversionStreamerInfo( expectedClass, bEl->GetClassVersion() ) &&
          !ptrClass->FindConversionStreamerInfo( expectedClass, bEl->GetCheckSum() ) ) {
         Error("SetBranchAddress", "The pointer type given \"%s\" does not correspond to the type needed \"%s\" by the branch: %s", ptrClass->GetName(), bEl->GetClassName(), branch->GetName());

         bEl->SetTargetClass( expectedClass->GetName() );
         return kClassMismatch;
      }
      else {

         bEl->SetTargetClass( ptrClass->GetName() );
         return kMatchConversion;
      }

   } else if (expectedClass && ptrClass && !expectedClass->InheritsFrom(ptrClass)) {

      if (expectedClass->GetCollectionProxy() && ptrClass->GetCollectionProxy() &&
          isBranchElement &&
          expectedClass->GetCollectionProxy()->GetValueClass() &&
          ptrClass->GetCollectionProxy()->GetValueClass() )
      {
         // In case of collection, we know how to convert them, if we know how to convert their content.
         // NOTE: we need to extend this to std::pair ...

         TClass *onfileValueClass = expectedClass->GetCollectionProxy()->GetValueClass();
         TClass *inmemValueClass = ptrClass->GetCollectionProxy()->GetValueClass();

         if (inmemValueClass->GetSchemaRules() &&
             inmemValueClass->GetSchemaRules()->HasRuleWithSourceClass(onfileValueClass->GetName() ) )
         {
            TBranchElement* bEl = (TBranchElement*)branch;
            bEl->SetTargetClass( ptrClass->GetName() );
            return kMatchConversionCollection;
         }
      }

      Error("SetBranchAddress", "The pointer type given (%s) does not correspond to the class needed (%s) by the branch: %s", ptrClass->GetName(), expectedClass->GetName(), branch->GetName());
      if (isBranchElement) {
         TBranchElement* bEl = (TBranchElement*)branch;
         bEl->SetTargetClass( expectedClass->GetName() );
      }
      return kClassMismatch;

   } else if ((expectedType != kOther_t) && (datatype != kOther_t) && (expectedType != kNoType_t) && (datatype != kNoType_t) && (expectedType != datatype)) {
      if (datatype != kChar_t) {
         // For backward compatibility we assume that (char*) was just a cast and/or a generic address
         Error("SetBranchAddress", "The pointer type given \"%s\" (%d) does not correspond to the type needed \"%s\" (%d) by the branch: %s",
               TDataType::GetTypeName(datatype), datatype, TDataType::GetTypeName(expectedType), expectedType, branch->GetName());
         return kMismatch;
      }
   } else if ((expectedClass && (datatype != kOther_t && datatype != kNoType_t && datatype != kInt_t)) ||
              (ptrClass && (expectedType != kOther_t && expectedType != kNoType_t && datatype != kInt_t)) ) {
      // Sometime a null pointer can look an int, avoid complaining in that case.
      if (expectedClass) {
         Error("SetBranchAddress", "The pointer type given \"%s\" (%d) does not correspond to the type needed \"%s\" by the branch: %s",
               TDataType::GetTypeName(datatype), datatype, expectedClass->GetName(), branch->GetName());
         if (isBranchElement) {
            TBranchElement* bEl = (TBranchElement*)branch;
            bEl->SetTargetClass( expectedClass->GetName() );
         }
      } else {
         // In this case, it is okay if the first data member is of the right type (to support the case where we are being passed
         // a struct).
         bool found = false;
         if (ptrClass->IsLoaded()) {
            TIter next(ptrClass->GetListOfRealData());
            TRealData *rdm;
            while ((rdm = (TRealData*)next())) {
               if (rdm->GetThisOffset() == 0) {
                  TDataType *dmtype = rdm->GetDataMember()->GetDataType();
                  if (dmtype) {
                     EDataType etype = (EDataType)dmtype->GetType();
                     if (etype == expectedType) {
                        found = true;
                     }
                  }
                  break;
               }
            }
         } else {
            TIter next(ptrClass->GetListOfDataMembers());
            TDataMember *dm;
            while ((dm = (TDataMember*)next())) {
               if (dm->GetOffset() == 0) {
                  TDataType *dmtype = dm->GetDataType();
                  if (dmtype) {
                     EDataType etype = (EDataType)dmtype->GetType();
                     if (etype == expectedType) {
                        found = true;
                     }
                  }
                  break;
               }
            }
         }
         if (found) {
            // let's check the size.
            TLeaf *last = (TLeaf*)branch->GetListOfLeaves()->Last();
            long len = last->GetOffset() + last->GetLenType() * last->GetLen();
            if (len <= ptrClass->Size()) {
               return kMatch;
            }
         }
         Error("SetBranchAddress", "The pointer type given \"%s\" does not correspond to the type needed \"%s\" (%d) by the branch: %s",
               ptrClass->GetName(), TDataType::GetTypeName(expectedType), expectedType, branch->GetName());
      }
      return kMismatch;
   }
   if (expectedClass && expectedClass->GetCollectionProxy() && dynamic_cast<TEmulatedCollectionProxy*>(expectedClass->GetCollectionProxy())) {
      Error("SetBranchAddress", writeStlWithoutProxyMsg,
            expectedClass->GetName(), branch->GetName(), expectedClass->GetName());
      if (isBranchElement) {
         TBranchElement* bEl = (TBranchElement*)branch;
         bEl->SetTargetClass( expectedClass->GetName() );
      }
      return kMissingCompiledCollectionProxy;
   }
   if (isBranchElement) {
      if (expectedClass) {
         TBranchElement* bEl = (TBranchElement*)branch;
         bEl->SetTargetClass( expectedClass->GetName() );
      } else if (expectedType != kNoType_t && expectedType != kOther_t) {
         return kMatch | kNeedEnableDecomposedObj;
      }
   }
   return kMatch;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a clone of this tree and copy nentries.
///
/// By default copy all entries.
/// The compression level of the cloned tree is set to the destination
/// file's compression level.
///
/// NOTE: Only active branches are copied.
/// NOTE: If the TTree is a TChain, the structure of the first TTree
///       is used for the copy.
///
/// IMPORTANT: The cloned tree stays connected with this tree until
///            this tree is deleted. In particular, any changes in
///            branch addresses in this tree are forwarded to the
///            clone trees, unless a branch in a clone tree has had
///            its address changed, in which case that change stays in
///            effect. When this tree is deleted, all the addresses of
///            the cloned tree are reset to their default values.
///
/// If 'option' contains the word 'fast' and nentries is -1, the
/// cloning will be done without unzipping or unstreaming the baskets
/// (i.e., a direct copy of the raw bytes on disk).
///
/// When 'fast' is specified, 'option' can also contain a sorting
/// order for the baskets in the output file.
///
/// There are currently 3 supported sorting order:
///
/// - SortBasketsByOffset (the default)
/// - SortBasketsByBranch
/// - SortBasketsByEntry
///
/// When using SortBasketsByOffset the baskets are written in the
/// output file in the same order as in the original file (i.e. the
/// baskets are sorted by their offset in the original file; Usually
/// this also means that the baskets are sorted by the index/number of
/// the _last_ entry they contain)
///
/// When using SortBasketsByBranch all the baskets of each individual
/// branches are stored contiguously. This tends to optimize reading
/// speed when reading a small number (1->5) of branches, since all
/// their baskets will be clustered together instead of being spread
/// across the file. However it might decrease the performance when
/// reading more branches (or the full entry).
///
/// When using SortBasketsByEntry the baskets with the lowest starting
/// entry are written first. (i.e. the baskets are sorted by the
/// index/number of the first entry they contain). This means that on
/// the file the baskets will be in the order in which they will be
/// needed when reading the whole tree sequentially.
///
/// For examples of CloneTree, see tutorials:
///
/// - copytree.C:
///     A macro to copy a subset of a TTree to a new TTree.
///     The input file has been generated by the program in
///     $ROOTSYS/test/Event with: Event 1000 1 1 1
///
/// - copytree2.C:
///     A macro to copy a subset of a TTree to a new TTree.
///     One branch of the new Tree is written to a separate file.
///     The input file has been generated by the program in
///     $ROOTSYS/test/Event with: Event 1000 1 1 1

TTree* TTree::CloneTree(Long64_t nentries /* = -1 */, Option_t* option /* = "" */)
{
   // Options
   Bool_t fastClone = kFALSE;

   TString opt = option;
   opt.ToLower();
   if (opt.Contains("fast")) {
      fastClone = kTRUE;
   }

   // If we are a chain, switch to the first tree.
   if ((fEntries > 0) && (LoadTree(0) < 0)) {
         // FIXME: We need an error message here.
         return 0;
   }

   // Note: For a tree we get the this pointer, for
   //       a chain we get the chain's current tree.
   TTree* thistree = GetTree();

   // We will use this to override the IO features on the cloned branches.
   ROOT::TIOFeatures features = this->GetIOFeatures();
   ;

   // Note: For a chain, the returned clone will be
   //       a clone of the chain's first tree.
   TTree* newtree = (TTree*) thistree->Clone();
   if (!newtree) {
      return 0;
   }

   // The clone should not delete any objects allocated by SetAddress().
   TObjArray* branches = newtree->GetListOfBranches();
   Int_t nb = branches->GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i) {
      TBranch* br = (TBranch*) branches->UncheckedAt(i);
      if (br->InheritsFrom(TBranchElement::Class())) {
         ((TBranchElement*) br)->ResetDeleteObject();
      }
   }

   // Add the new tree to the list of clones so that
   // we can later inform it of changes to branch addresses.
   thistree->AddClone(newtree);
   if (thistree != this) {
      // In case this object is a TChain, add the clone
      // also to the TChain's list of clones.
      AddClone(newtree);
   }

   newtree->Reset();

   TDirectory* ndir = newtree->GetDirectory();
   TFile* nfile = 0;
   if (ndir) {
      nfile = ndir->GetFile();
   }
   Int_t newcomp = -1;
   if (nfile) {
      newcomp = nfile->GetCompressionSettings();
   }

   //
   // Delete non-active branches from the clone.
   //
   // Note: If we are a chain, this does nothing
   //       since chains have no leaves.
   TObjArray* leaves = newtree->GetListOfLeaves();
   Int_t nleaves = leaves->GetEntriesFast();
   for (Int_t lndx = 0; lndx < nleaves; ++lndx) {
      TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(lndx);
      if (!leaf) {
         continue;
      }
      TBranch* branch = leaf->GetBranch();
      if (branch && (newcomp > -1)) {
         branch->SetCompressionSettings(newcomp);
      }
      if (branch) branch->SetIOFeatures(features);
      if (!branch || !branch->TestBit(kDoNotProcess)) {
         continue;
      }
      // size might change at each iteration of the loop over the leaves.
      nb = branches->GetEntriesFast();
      for (Long64_t i = 0; i < nb; ++i) {
         TBranch* br = (TBranch*) branches->UncheckedAt(i);
         if (br == branch) {
            branches->RemoveAt(i);
            delete br;
            br = 0;
            branches->Compress();
            break;
         }
         TObjArray* lb = br->GetListOfBranches();
         Int_t nb1 = lb->GetEntriesFast();
         for (Int_t j = 0; j < nb1; ++j) {
            TBranch* b1 = (TBranch*) lb->UncheckedAt(j);
            if (!b1) {
               continue;
            }
            if (b1 == branch) {
               lb->RemoveAt(j);
               delete b1;
               b1 = 0;
               lb->Compress();
               break;
            }
            TObjArray* lb1 = b1->GetListOfBranches();
            Int_t nb2 = lb1->GetEntriesFast();
            for (Int_t k = 0; k < nb2; ++k) {
               TBranch* b2 = (TBranch*) lb1->UncheckedAt(k);
               if (!b2) {
                  continue;
               }
               if (b2 == branch) {
                  lb1->RemoveAt(k);
                  delete b2;
                  b2 = 0;
                  lb1->Compress();
                  break;
               }
            }
         }
      }
   }
   leaves->Compress();

   // Copy MakeClass status.
   newtree->SetMakeClass(fMakeClass);

   // Copy branch addresses.
   CopyAddresses(newtree);

   //
   // Copy entries if requested.
   //

   if (nentries != 0) {
      if (fastClone && (nentries < 0)) {
         if ( newtree->CopyEntries( this, -1, option ) < 0 ) {
            // There was a problem!
            Error("CloneTTree", "TTree has not been cloned\n");
            delete newtree;
            newtree = 0;
            return 0;
         }
      } else {
         newtree->CopyEntries( this, nentries, option );
      }
   }

   return newtree;
}

////////////////////////////////////////////////////////////////////////////////
/// Set branch addresses of passed tree equal to ours.
/// If undo is true, reset the branch address instead of copying them.
/// This insures 'separation' of a cloned tree from its original

void TTree::CopyAddresses(TTree* tree, Bool_t undo)
{
   // Copy branch addresses starting from branches.
   TObjArray* branches = GetListOfBranches();
   Int_t nbranches = branches->GetEntriesFast();
   for (Int_t i = 0; i < nbranches; ++i) {
      TBranch* branch = (TBranch*) branches->UncheckedAt(i);
      if (branch->TestBit(kDoNotProcess)) {
         continue;
      }
      if (undo) {
         TBranch* br = tree->GetBranch(branch->GetName());
         tree->ResetBranchAddress(br);
      } else {
         char* addr = branch->GetAddress();
         if (!addr) {
            if (branch->IsA() == TBranch::Class()) {
               // If the branch was created using a leaflist, the branch itself may not have
               // an address but the leaf might already.
               TLeaf *firstleaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
               if (!firstleaf || firstleaf->GetValuePointer()) {
                  // Either there is no leaf (and thus no point in copying the address)
                  // or the leaf has an address but we can not copy it via the branche
                  // this will be copied via the next loop (over the leaf).
                  continue;
               }
            }
            // Note: This may cause an object to be allocated.
            branch->SetAddress(0);
            addr = branch->GetAddress();
         }
         TBranch* br = tree->GetBranch(branch->GetFullName());
         if (br) {
            if (br->GetMakeClass() != branch->GetMakeClass())
               br->SetMakeClass(branch->GetMakeClass());
            br->SetAddress(addr);
            // The copy does not own any object allocated by SetAddress().
            if (br->InheritsFrom(TBranchElement::Class())) {
               ((TBranchElement*) br)->ResetDeleteObject();
            }
         } else {
            Warning("CopyAddresses", "Could not find branch named '%s' in tree named '%s'", branch->GetName(), tree->GetName());
         }
      }
   }

   // Copy branch addresses starting from leaves.
   TObjArray* tleaves = tree->GetListOfLeaves();
   Int_t ntleaves = tleaves->GetEntriesFast();
   std::set<TLeaf*> updatedLeafCount;
   for (Int_t i = 0; i < ntleaves; ++i) {
      TLeaf* tleaf = (TLeaf*) tleaves->UncheckedAt(i);
      TBranch* tbranch = tleaf->GetBranch();
      TBranch* branch = GetBranch(tbranch->GetName());
      if (!branch) {
         continue;
      }
      TLeaf* leaf = branch->GetLeaf(tleaf->GetName());
      if (!leaf) {
         continue;
      }
      if (branch->TestBit(kDoNotProcess)) {
         continue;
      }
      if (undo) {
         // Now we know whether the address has been transfered
         tree->ResetBranchAddress(tbranch);
      } else {
         TBranchElement *mother = dynamic_cast<TBranchElement*>(leaf->GetBranch()->GetMother());
         bool needAddressReset = false;
         if (leaf->GetLeafCount() && (leaf->TestBit(TLeaf::kNewValue) || !leaf->GetValuePointer() || (mother && mother->IsObjectOwner())) && tleaf->GetLeafCount())
         {
            // If it is an array and it was allocated by the leaf itself,
            // let's make sure it is large enough for the incoming data.
            if (leaf->GetLeafCount()->GetMaximum() < tleaf->GetLeafCount()->GetMaximum()) {
               leaf->GetLeafCount()->IncludeRange( tleaf->GetLeafCount() );
               updatedLeafCount.insert(leaf->GetLeafCount());
               needAddressReset = true;
             } else {
               needAddressReset = (updatedLeafCount.find(leaf->GetLeafCount()) != updatedLeafCount.end());
             }
         }
         if (needAddressReset && leaf->GetValuePointer()) {
            if (leaf->IsA() == TLeafElement::Class() && mother)
               mother->ResetAddress();
            else
               leaf->SetAddress(nullptr);
         }
         if (!branch->GetAddress() && !leaf->GetValuePointer()) {
            // We should attempts to set the address of the branch.
            // something like:
            //(TBranchElement*)branch->GetMother()->SetAddress(0)
            //plus a few more subtleties (see TBranchElement::GetEntry).
            //but for now we go the simplest route:
            //
            // Note: This may result in the allocation of an object.
            branch->SetupAddresses();
         }
         if (branch->GetAddress()) {
            tree->SetBranchAddress(branch->GetName(), (void*) branch->GetAddress());
            TBranch* br = tree->GetBranch(branch->GetName());
            if (br) {
               // The copy does not own any object allocated by SetAddress().
               // FIXME: We do too much here, br may not be a top-level branch.
               if (br->InheritsFrom(TBranchElement::Class())) {
                  ((TBranchElement*) br)->ResetDeleteObject();
               }
            } else {
               Warning("CopyAddresses", "Could not find branch named '%s' in tree named '%s'", branch->GetName(), tree->GetName());
            }
         } else {
            tleaf->SetAddress(leaf->GetValuePointer());
         }
      }
   }

   if (undo &&
       ( tree->IsA()->InheritsFrom("TNtuple") || tree->IsA()->InheritsFrom("TNtupleD") )
       ) {
      tree->ResetBranchAddresses();
   }
}

namespace {

   enum EOnIndexError { kDrop, kKeep, kBuild };

   static Bool_t R__HandleIndex(EOnIndexError onIndexError, TTree *newtree, TTree *oldtree)
   {
      // Return true if we should continue to handle indices, false otherwise.

      Bool_t withIndex = kTRUE;

      if ( newtree->GetTreeIndex() ) {
         if ( oldtree->GetTree()->GetTreeIndex() == 0 ) {
            switch (onIndexError) {
               case kDrop:
                  delete newtree->GetTreeIndex();
                  newtree->SetTreeIndex(0);
                  withIndex = kFALSE;
                  break;
               case kKeep:
                  // Nothing to do really.
                  break;
               case kBuild:
                  // Build the index then copy it
                  if (oldtree->GetTree()->BuildIndex(newtree->GetTreeIndex()->GetMajorName(), newtree->GetTreeIndex()->GetMinorName())) {
                     newtree->GetTreeIndex()->Append(oldtree->GetTree()->GetTreeIndex(), kTRUE);
                     // Clean up
                     delete oldtree->GetTree()->GetTreeIndex();
                     oldtree->GetTree()->SetTreeIndex(0);
                  }
                  break;
            }
         } else {
            newtree->GetTreeIndex()->Append(oldtree->GetTree()->GetTreeIndex(), kTRUE);
         }
      } else if ( oldtree->GetTree()->GetTreeIndex() != 0 ) {
         // We discover the first index in the middle of the chain.
         switch (onIndexError) {
            case kDrop:
               // Nothing to do really.
               break;
            case kKeep: {
               TVirtualIndex *index = (TVirtualIndex*) oldtree->GetTree()->GetTreeIndex()->Clone();
               index->SetTree(newtree);
               newtree->SetTreeIndex(index);
               break;
            }
            case kBuild:
               if (newtree->GetEntries() == 0) {
                  // Start an index.
                  TVirtualIndex *index = (TVirtualIndex*) oldtree->GetTree()->GetTreeIndex()->Clone();
                  index->SetTree(newtree);
                  newtree->SetTreeIndex(index);
               } else {
                  // Build the index so far.
                  if (newtree->BuildIndex(oldtree->GetTree()->GetTreeIndex()->GetMajorName(), oldtree->GetTree()->GetTreeIndex()->GetMinorName())) {
                     newtree->GetTreeIndex()->Append(oldtree->GetTree()->GetTreeIndex(), kTRUE);
                  }
               }
               break;
         }
      } else if ( onIndexError == kDrop ) {
         // There is no index on this or on tree->GetTree(), we know we have to ignore any further
         // index
         withIndex = kFALSE;
      }
      return withIndex;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy nentries from given tree to this tree.
/// This routines assumes that the branches that intended to be copied are
/// already connected.   The typical case is that this tree was created using
/// tree->CloneTree(0).
///
/// By default copy all entries.
///
/// Returns number of bytes copied to this tree.
///
/// If 'option' contains the word 'fast' and nentries is -1, the cloning will be
/// done without unzipping or unstreaming the baskets (i.e., a direct copy of the
/// raw bytes on disk).
///
/// When 'fast' is specified, 'option' can also contains a sorting order for the
/// baskets in the output file.
///
/// There are currently 3 supported sorting order:
///
/// - SortBasketsByOffset (the default)
/// - SortBasketsByBranch
/// - SortBasketsByEntry
///
/// See TTree::CloneTree for a detailed explanation of the semantics of these 3 options.
///
/// If the tree or any of the underlying tree of the chain has an index, that index and any
/// index in the subsequent underlying TTree objects will be merged.
///
/// There are currently three 'options' to control this merging:
/// - NoIndex             : all the TTreeIndex object are dropped.
/// - DropIndexOnError    : if any of the underlying TTree object do no have a TTreeIndex,
///                          they are all dropped.
/// - AsIsIndexOnError [default]: In case of missing TTreeIndex, the resulting TTree index has gaps.
/// - BuildIndexOnError : If any of the underlying TTree objects do not have a TTreeIndex,
///                          all TTreeIndex are 'ignored' and the missing piece are rebuilt.

Long64_t TTree::CopyEntries(TTree* tree, Long64_t nentries /* = -1 */, Option_t* option /* = "" */)
{
   if (!tree) {
      return 0;
   }
   // Options
   TString opt = option;
   opt.ToLower();
   Bool_t fastClone = opt.Contains("fast");
   Bool_t withIndex = !opt.Contains("noindex");
   EOnIndexError onIndexError;
   if (opt.Contains("asisindex")) {
      onIndexError = kKeep;
   } else if (opt.Contains("buildindex")) {
      onIndexError = kBuild;
   } else if (opt.Contains("dropindex")) {
      onIndexError = kDrop;
   } else {
      onIndexError = kBuild;
   }
   Ssiz_t cacheSizeLoc = opt.Index("cachesize=");
   Int_t cacheSize = -1;
   if (cacheSizeLoc != TString::kNPOS) {
      // If the parse faile, cacheSize stays at -1.
      Ssiz_t cacheSizeEnd = opt.Index(" ",cacheSizeLoc+10) - (cacheSizeLoc+10);
      TSubString cacheSizeStr( opt(cacheSizeLoc+10,cacheSizeEnd) );
      auto parseResult = ROOT::FromHumanReadableSize(cacheSizeStr,cacheSize);
      if (parseResult == ROOT::EFromHumanReadableSize::kParseFail) {
         Warning("CopyEntries","The cachesize option can not be parsed: %s. The default size will be used.",cacheSizeStr.String().Data());
      } else if (parseResult == ROOT::EFromHumanReadableSize::kOverflow) {
         double m;
         const char *munit = nullptr;
         ROOT::ToHumanReadableSize(std::numeric_limits<decltype(cacheSize)>::max(),false,&m,&munit);

         Warning("CopyEntries","The cachesize option is too large: %s (%g%s max). The default size will be used.",cacheSizeStr.String().Data(),m,munit);
      }
   }
   if (gDebug > 0 && cacheSize != -1) Info("CopyEntries","Using Cache size: %d\n",cacheSize);

   Long64_t nbytes = 0;
   Long64_t treeEntries = tree->GetEntriesFast();
   if (nentries < 0) {
      nentries = treeEntries;
   } else if (nentries > treeEntries) {
      nentries = treeEntries;
   }

   if (fastClone && (nentries < 0 || nentries == tree->GetEntriesFast())) {
      // Quickly copy the basket without decompression and streaming.
      Long64_t totbytes = GetTotBytes();
      for (Long64_t i = 0; i < nentries; i += tree->GetTree()->GetEntries()) {
         if (tree->LoadTree(i) < 0) {
            break;
         }
         if ( withIndex ) {
            withIndex = R__HandleIndex( onIndexError, this, tree );
         }
         if (this->GetDirectory()) {
            TFile* file2 = this->GetDirectory()->GetFile();
            if (file2 && (file2->GetEND() > TTree::GetMaxTreeSize())) {
               if (this->GetDirectory() == (TDirectory*) file2) {
                  this->ChangeFile(file2);
               }
            }
         }
         TTreeCloner cloner(tree->GetTree(), this, option, TTreeCloner::kNoWarnings);
         if (cloner.IsValid()) {
            this->SetEntries(this->GetEntries() + tree->GetTree()->GetEntries());
            if (cacheSize != -1) cloner.SetCacheSize(cacheSize);
            cloner.Exec();
         } else {
            if (i == 0) {
               Warning("CopyEntries","%s",cloner.GetWarning());
               // If the first cloning does not work, something is really wrong
               // (since apriori the source and target are exactly the same structure!)
               return -1;
            } else {
               if (cloner.NeedConversion()) {
                  TTree *localtree = tree->GetTree();
                  Long64_t tentries = localtree->GetEntries();
                  for (Long64_t ii = 0; ii < tentries; ii++) {
                     if (localtree->GetEntry(ii) <= 0) {
                        break;
                     }
                     this->Fill();
                  }
                  if (this->GetTreeIndex()) {
                     this->GetTreeIndex()->Append(tree->GetTree()->GetTreeIndex(), kTRUE);
                  }
               } else {
                  Warning("CopyEntries","%s",cloner.GetWarning());
                  if (tree->GetDirectory() && tree->GetDirectory()->GetFile()) {
                     Warning("CopyEntries", "Skipped file %s\n", tree->GetDirectory()->GetFile()->GetName());
                  } else {
                     Warning("CopyEntries", "Skipped file number %d\n", tree->GetTreeNumber());
                  }
               }
            }
         }

      }
      if (this->GetTreeIndex()) {
         this->GetTreeIndex()->Append(0,kFALSE); // Force the sorting
      }
      nbytes = GetTotBytes() - totbytes;
   } else {
      if (nentries < 0) {
         nentries = treeEntries;
      } else if (nentries > treeEntries) {
         nentries = treeEntries;
      }
      Int_t treenumber = -1;
      for (Long64_t i = 0; i < nentries; i++) {
         if (tree->LoadTree(i) < 0) {
            break;
         }
         if (treenumber != tree->GetTreeNumber()) {
            if ( withIndex ) {
               withIndex = R__HandleIndex( onIndexError, this, tree );
            }
            treenumber = tree->GetTreeNumber();
         }
         if (tree->GetEntry(i) <= 0) {
            break;
         }
         nbytes += this->Fill();
      }
      if (this->GetTreeIndex()) {
         this->GetTreeIndex()->Append(0,kFALSE); // Force the sorting
      }
   }
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a tree with selection.
///
/// ### Important:
///
/// The returned copied tree stays connected with the original tree
/// until the original tree is deleted.  In particular, any changes
/// to the branch addresses in the original tree are also made to
/// the copied tree.  Any changes made to the branch addresses of the
/// copied tree are overridden anytime the original tree changes its
/// branch addresses.  When the original tree is deleted, all the
/// branch addresses of the copied tree are set to zero.
///
/// For examples of CopyTree, see the tutorials:
///
/// - copytree.C:
/// Example macro to copy a subset of a tree to a new tree.
/// The input file was generated by running the program in
/// $ROOTSYS/test/Event in this way:
/// ~~~ {.cpp}
///     ./Event 1000 1 1 1
/// ~~~
/// - copytree2.C
/// Example macro to copy a subset of a tree to a new tree.
/// One branch of the new tree is written to a separate file.
/// The input file was generated by running the program in
/// $ROOTSYS/test/Event in this way:
/// ~~~ {.cpp}
///     ./Event 1000 1 1 1
/// ~~~
/// - copytree3.C
/// Example macro to copy a subset of a tree to a new tree.
/// Only selected entries are copied to the new tree.
/// NOTE that only the active branches are copied.

TTree* TTree::CopyTree(const char* selection, Option_t* option /* = 0 */, Long64_t nentries /* = TTree::kMaxEntries */, Long64_t firstentry /* = 0 */)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->CopyTree(selection, option, nentries, firstentry);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a basket for this tree and given branch.

TBasket* TTree::CreateBasket(TBranch* branch)
{
   if (!branch) {
      return 0;
   }
   return new TBasket(branch->GetName(), GetName(), branch);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete this tree from memory or/and disk.
///
/// - if option == "all" delete Tree object from memory AND from disk
///                     all baskets on disk are deleted. All keys with same name
///                     are deleted.
/// - if option =="" only Tree object in memory is deleted.

void TTree::Delete(Option_t* option /* = "" */)
{
   TFile *file = GetCurrentFile();

   // delete all baskets and header from file
   if (file && !strcmp(option,"all")) {
      if (!file->IsWritable()) {
         Error("Delete","File : %s is not writable, cannot delete Tree:%s", file->GetName(),GetName());
         return;
      }

      //find key and import Tree header in memory
      TKey *key = fDirectory->GetKey(GetName());
      if (!key) return;

      TDirectory *dirsav = gDirectory;
      file->cd();

      //get list of leaves and loop on all the branches baskets
      TIter next(GetListOfLeaves());
      TLeaf *leaf;
      char header[16];
      Int_t ntot  = 0;
      Int_t nbask = 0;
      Int_t nbytes,objlen,keylen;
      while ((leaf = (TLeaf*)next())) {
         TBranch *branch = leaf->GetBranch();
         Int_t nbaskets = branch->GetMaxBaskets();
         for (Int_t i=0;i<nbaskets;i++) {
            Long64_t pos = branch->GetBasketSeek(i);
            if (!pos) continue;
            TFile *branchFile = branch->GetFile();
            if (!branchFile) continue;
            branchFile->GetRecordHeader(header,pos,16,nbytes,objlen,keylen);
            if (nbytes <= 0) continue;
            branchFile->MakeFree(pos,pos+nbytes-1);
            ntot += nbytes;
            nbask++;
         }
      }

      // delete Tree header key and all keys with the same name
      // A Tree may have been saved many times. Previous cycles are invalid.
      while (key) {
         ntot += key->GetNbytes();
         key->Delete();
         delete key;
         key = fDirectory->GetKey(GetName());
      }
      if (dirsav) dirsav->cd();
      if (gDebug) Info("TTree::Delete", "Deleting Tree: %s: %d baskets deleted. Total space freed = %d bytes\n",GetName(),nbask,ntot);
   }

   if (fDirectory) {
      fDirectory->Remove(this);
      //delete the file cache if it points to this Tree
      MoveReadCache(file,0);
      fDirectory = 0;
      ResetBit(kMustCleanup);
   }

   // Delete object from CINT symbol table so it can not be used anymore.
   gCling->DeleteGlobal(this);

   // Warning: We have intentional invalidated this object while inside a member function!
   delete this;
}

 ///////////////////////////////////////////////////////////////////////////////
 /// Called by TKey and TObject::Clone to automatically add us to a directory
 /// when we are read from a file.

void TTree::DirectoryAutoAdd(TDirectory* dir)
{
   if (fDirectory == dir) return;
   if (fDirectory) {
      fDirectory->Remove(this);
      // Delete or move the file cache if it points to this Tree
      TFile *file = fDirectory->GetFile();
      MoveReadCache(file,dir);
   }
   fDirectory = dir;
   TBranch* b = 0;
   TIter next(GetListOfBranches());
   while((b = (TBranch*) next())) {
      b->UpdateFile();
   }
   if (fBranchRef) {
      fBranchRef->UpdateFile();
   }
   if (fDirectory) fDirectory->Append(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expression varexp for specified entries.
///
/// \return -1 in case of error or number of selected events in case of success.
///
/// This function accepts TCut objects as arguments.
/// Useful to use the string operator +
///
/// Example:
///
/// ~~~ {.cpp}
///     ntuple.Draw("x",cut1+cut2+cut3);
/// ~~~


Long64_t TTree::Draw(const char* varexp, const TCut& selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   return TTree::Draw(varexp, selection.GetTitle(), option, nentries, firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expression varexp for specified entries.
///
/// \return -1 in case of error or number of selected events in case of success.
///
/// \param [in] varexp is an expression of the general form
///  - "e1"           produces a 1-d histogram (TH1F) of expression "e1"
///  - "e1:e2"        produces an unbinned 2-d scatter-plot (TGraph) of "e1"
///                   on the y-axis versus "e2" on the x-axis
///  - "e1:e2:e3"     produces an unbinned 3-d scatter-plot (TPolyMarker3D) of "e1"
///                   vs "e2" vs "e3" on the x-, y-, z-axis, respectively.
///  - "e1:e2:e3:e4"  produces an unbinned 3-d scatter-plot (TPolyMarker3D) of "e1"
///                   vs "e2" vs "e3" and "e4" mapped on the current color palette.
///                   (to create histograms in the 2, 3, and 4 dimensional case,
///                   see section "Saving the result of Draw to an histogram")
///
///   Example:
///    -  varexp = x     simplest case: draw a 1-Dim distribution of column named x
///    -  varexp = sqrt(x)            : draw distribution of sqrt(x)
///    -  varexp = x*y/z
///    -  varexp = y:sqrt(x) 2-Dim distribution of y versus sqrt(x)
///    -  varexp = px:py:pz:2.5*E  produces a 3-d scatter-plot of px vs py ps pz
///               and the color number of each marker will be 2.5*E.
///               If the color number is negative it is set to 0.
///               If the color number is greater than the current number of colors
///               it is set to the highest color number.The default number of
///               colors is 50. see TStyle::SetPalette for setting a new color palette.
///
///   Note that the variables e1, e2 or e3 may contain a selection.
///   example, if e1= x*(y<0), the value histogrammed will be x if y<0
///   and will be 0 otherwise.
///
///   The expressions can use all the operations and build-in functions
///   supported by TFormula (See TFormula::Analyze), including free
///   standing function taking numerical arguments (TMath::Bessel).
///   In addition, you can call member functions taking numerical
///   arguments. For example:
///   ~~~ {.cpp}
///       TMath::BreitWigner(fPx,3,2)
///       event.GetHistogram()->GetXaxis()->GetXmax()
///   ~~~
///   Note: You can only pass expression that depend on the TTree's data
///   to static functions and you can only call non-static member function
///   with 'fixed' parameters.
///
/// \param [in] selection is an expression with a combination of the columns.
///   In a selection all the C++ operators are authorized.
///   The value corresponding to the selection expression is used as a weight
///   to fill the histogram.
///   If the expression includes only boolean operations, the result
///   is 0 or 1. If the result is 0, the histogram is not filled.
///   In general, the expression may be of the form:
///   ~~~ {.cpp}
///       value*(boolean expression)
///   ~~~
///   if boolean expression is true, the histogram is filled with
///   a `weight = value`.
///   Examples:
///    -  selection1 = "x<y && sqrt(z)>3.2"
///    -  selection2 = "(x+y)*(sqrt(z)>3.2)"
///    -  selection1 returns a weight = 0 or 1
///    -  selection2 returns a weight = x+y if sqrt(z)>3.2
///                  returns a weight = 0 otherwise.
///
/// \param [in] option is the drawing option.
///    - When an histogram is produced it can be any histogram drawing option
///      listed in THistPainter.
///    - when no option is specified:
///        - the default histogram drawing option is used
///          if the expression is of the form "e1".
///        - if the expression is of the form "e1:e2"or "e1:e2:e3" a cloud of
///          unbinned 2D or 3D points is drawn respectively.
///        - if the expression  has four fields "e1:e2:e3:e4" a cloud of unbinned 3D
///          points is produced with e1 vs e2 vs e3, and e4 is mapped on the current color
///          palette.
///    - If option COL is specified when varexp has three fields:
///   ~~~ {.cpp}
///        tree.Draw("e1:e2:e3","","col");
///   ~~~
///      a 2D scatter is produced with e1 vs e2, and e3 is mapped on the current
///      color palette. The colors for e3 are evaluated once in linear scale before
///      painting. Therefore changing the pad to log scale along Z as no effect
///      on the colors.
///    - if expression has more than four fields the option "PARA"or "CANDLE"
///      can be used.
///    - If option contains the string "goff", no graphics is generated.
///
/// \param [in] nentries is the number of entries to process (default is all)
///
/// \param [in] firstentry is the first entry to process (default is 0)
///
/// ### Drawing expressions using arrays and array elements
///
/// Let assumes, a leaf fMatrix, on the branch fEvent, which is a 3 by 3 array,
/// or a TClonesArray.
/// In a TTree::Draw expression you can now access fMatrix using the following
/// syntaxes:
///
/// | String passed   | What is used for each entry of the tree
/// |-----------------|--------------------------------------------------------|
/// | `fMatrix`       | the 9 elements of fMatrix |
/// | `fMatrix[][]`   | the 9 elements of fMatrix |
/// | `fMatrix[2][2]` | only the elements fMatrix[2][2] |
/// | `fMatrix[1]`    | the 3 elements fMatrix[1][0], fMatrix[1][1] and fMatrix[1][2] |
/// | `fMatrix[1][]`  | the 3 elements fMatrix[1][0], fMatrix[1][1] and fMatrix[1][2] |
/// | `fMatrix[][0]`  | the 3 elements fMatrix[0][0], fMatrix[1][0] and fMatrix[2][0] |
///
/// "fEvent.fMatrix...." same as "fMatrix..." (unless there is more than one leaf named fMatrix!).
///
/// In summary, if a specific index is not specified for a dimension, TTree::Draw
/// will loop through all the indices along this dimension.  Leaving off the
/// last (right most) dimension of specifying then with the two characters '[]'
/// is equivalent.  For variable size arrays (and TClonesArray) the range
/// of the first dimension is recalculated for each entry of the tree.
/// You can also specify the index as an expression of any other variables from the
/// tree.
///
/// TTree::Draw also now properly handling operations involving 2 or more arrays.
///
/// Let assume a second matrix fResults[5][2], here are a sample of some
/// of the possible combinations, the number of elements they produce and
/// the loop used:
///
/// | expression                       | element(s) | Loop                     |
/// |----------------------------------|------------|--------------------------|
/// | `fMatrix[2][1] - fResults[5][2]` |  one       | no loop |
/// | `fMatrix[2][]  - fResults[5][2]` |  three     | on 2nd dim fMatrix |
/// | `fMatrix[2][]  - fResults[5][]`  |  two       | on both 2nd dimensions |
/// | `fMatrix[][2]  - fResults[][1]`  |  three     | on both 1st dimensions |
/// | `fMatrix[][2]  - fResults[][]`   |  six       | on both 1st and 2nd dimensions of fResults |
/// | `fMatrix[][2]  - fResults[3][]`  |  two       | on 1st dim of fMatrix and 2nd of fResults (at the same time) |
/// | `fMatrix[][]   - fResults[][]`   |  six       | on 1st dim then on  2nd dim |
/// | `fMatrix[][fResult[][]]`         |  30        | on 1st dim of fMatrix then on both dimensions of fResults.  The value if fResults[j][k] is used as the second index of fMatrix.|
///
///
/// In summary, TTree::Draw loops through all unspecified dimensions.  To
/// figure out the range of each loop, we match each unspecified dimension
/// from left to right (ignoring ALL dimensions for which an index has been
/// specified), in the equivalent loop matched dimensions use the same index
/// and are restricted to the smallest range (of only the matched dimensions).
/// When involving variable arrays, the range can of course be different
/// for each entry of the tree.
///
/// So the loop equivalent to "fMatrix[][2] - fResults[3][]" is:
/// ~~~ {.cpp}
///     for (Int_t i0; i < min(3,2); i++) {
///        use the value of (fMatrix[i0][2] - fMatrix[3][i0])
///     }
/// ~~~
/// So the loop equivalent to "fMatrix[][2] - fResults[][]" is:
/// ~~~ {.cpp}
///     for (Int_t i0; i < min(3,5); i++) {
///        for (Int_t i1; i1 < 2; i1++) {
///           use the value of (fMatrix[i0][2] - fMatrix[i0][i1])
///        }
///     }
/// ~~~
/// So the loop equivalent to "fMatrix[][] - fResults[][]" is:
/// ~~~ {.cpp}
///     for (Int_t i0; i < min(3,5); i++) {
///        for (Int_t i1; i1 < min(3,2); i1++) {
///           use the value of (fMatrix[i0][i1] - fMatrix[i0][i1])
///        }
///     }
/// ~~~
/// So the loop equivalent to "fMatrix[][fResults[][]]" is:
/// ~~~ {.cpp}
///     for (Int_t i0; i0 < 3; i0++) {
///        for (Int_t j2; j2 < 5; j2++) {
///           for (Int_t j3; j3 < 2; j3++) {
///              i1 = fResults[j2][j3];
///              use the value of fMatrix[i0][i1]
///        }
///     }
/// ~~~
/// ### Retrieving the result of Draw
///
/// By default a temporary histogram called `htemp` is created. It will be:
///
///  - A TH1F* in case of a mono-dimensional distribution: `Draw("e1")`,
///  - A TH2F* in case of a bi-dimensional distribution: `Draw("e1:e2")`,
///  - A TH3F* in case of a three-dimensional distribution: `Draw("e1:e2:e3")`.
///
/// In the one dimensional case the `htemp` is filled and drawn whatever the drawing
/// option is.
///
/// In the two and three dimensional cases, with the default drawing option (`""`),
/// a cloud of points is drawn and the histogram `htemp` is not filled. For all the other
/// drawing options `htemp` will be filled.
///
/// In all cases `htemp` can be retrieved by calling:
///
/// ~~~ {.cpp}
///     auto htemp = (TH1F*)gPad->GetPrimitive("htemp"); // 1D
///     auto htemp = (TH2F*)gPad->GetPrimitive("htemp"); // 2D
///     auto htemp = (TH3F*)gPad->GetPrimitive("htemp"); // 3D
/// ~~~
///
/// In the two dimensional case (`Draw("e1;e2")`), with the default drawing option, the
/// data is filled into a TGraph named `Graph`. This TGraph can be retrieved by
/// calling
///
/// ~~~ {.cpp}
///     auto graph = (TGraph*)gPad->GetPrimitive("Graph");
/// ~~~
///
/// For the three and four dimensional cases, with the default drawing option, an unnamed
/// TPolyMarker3D is produced, and therefore cannot be retrieved.
///
/// In all cases `htemp` can be used to access the axes. For instance in the 2D case:
///
/// ~~~ {.cpp}
///     auto htemp = (TH2F*)gPad->GetPrimitive("htemp");
///     auto xaxis = htemp->GetXaxis();
/// ~~~
///
/// When the option `"A"` is used (with TGraph painting option) to draw a 2D
/// distribution:
/// ~~~ {.cpp}
///     tree.Draw("e1:e2","","A*");
/// ~~~
/// a scatter plot is produced (with stars in that case) but the axis creation is
/// delegated to TGraph and `htemp` is not created.
///
/// ### Saving the result of Draw to a histogram
///
/// If `varexp` contains `>>hnew` (following the variable(s) name(s)),
/// the new histogram called `hnew` is created and it is kept in the current
/// directory (and also the current pad). This works for all dimensions.
///
/// Example:
/// ~~~ {.cpp}
///     tree.Draw("sqrt(x)>>hsqrt","y>0")
/// ~~~
/// will draw `sqrt(x)` and save the histogram as "hsqrt" in the current
/// directory. To retrieve it do:
/// ~~~ {.cpp}
///     TH1F *hsqrt = (TH1F*)gDirectory->Get("hsqrt");
/// ~~~
/// The binning information is taken from the environment variables
/// ~~~ {.cpp}
///     Hist.Binning.?D.?
/// ~~~
/// In addition, the name of the histogram can be followed by up to 9
/// numbers between '(' and ')', where the numbers describe the
/// following:
///
/// -  1 - bins in x-direction
/// -  2 - lower limit in x-direction
/// -  3 - upper limit in x-direction
/// -  4-6 same for y-direction
/// -  7-9 same for z-direction
///
/// When a new binning is used the new value will become the default.
/// Values can be skipped.
///
/// Example:
/// ~~~ {.cpp}
///     tree.Draw("sqrt(x)>>hsqrt(500,10,20)")
///          // plot sqrt(x) between 10 and 20 using 500 bins
///     tree.Draw("sqrt(x):sin(y)>>hsqrt(100,10,60,50,.1,.5)")
///          // plot sqrt(x) against sin(y)
///          // 100 bins in x-direction; lower limit on x-axis is 10; upper limit is 60
///          //  50 bins in y-direction; lower limit on y-axis is .1; upper limit is .5
/// ~~~
/// By default, the specified histogram is reset.
/// To continue to append data to an existing histogram, use "+" in front
/// of the histogram name.
///
/// A '+' in front of the histogram name is ignored, when the name is followed by
/// binning information as described in the previous paragraph.
/// ~~~ {.cpp}
///     tree.Draw("sqrt(x)>>+hsqrt","y>0")
/// ~~~
/// will not reset `hsqrt`, but will continue filling. This works for 1-D, 2-D
/// and 3-D histograms.
///
/// ### Accessing collection objects
///
/// TTree::Draw default's handling of collections is to assume that any
/// request on a collection pertain to it content.  For example, if fTracks
/// is a collection of Track objects, the following:
/// ~~~ {.cpp}
///     tree->Draw("event.fTracks.fPx");
/// ~~~
/// will plot the value of fPx for each Track objects inside the collection.
/// Also
/// ~~~ {.cpp}
///     tree->Draw("event.fTracks.size()");
/// ~~~
/// would plot the result of the member function Track::size() for each
/// Track object inside the collection.
/// To access information about the collection itself, TTree::Draw support
/// the '@' notation.  If a variable which points to a collection is prefixed
/// or postfixed with '@', the next part of the expression will pertain to
/// the collection object.  For example:
/// ~~~ {.cpp}
///     tree->Draw("event.@fTracks.size()");
/// ~~~
/// will plot the size of the collection referred to by `fTracks` (i.e the number
/// of Track objects).
///
/// ### Drawing 'objects'
///
/// When a class has a member function named AsDouble or AsString, requesting
/// to directly draw the object will imply a call to one of the 2 functions.
/// If both AsDouble and AsString are present, AsDouble will be used.
/// AsString can return either a char*, a std::string or a TString.s
/// For example, the following
/// ~~~ {.cpp}
///     tree->Draw("event.myTTimeStamp");
/// ~~~
/// will draw the same histogram as
/// ~~~ {.cpp}
///     tree->Draw("event.myTTimeStamp.AsDouble()");
/// ~~~
/// In addition, when the object is a type TString or std::string, TTree::Draw
/// will call respectively `TString::Data` and `std::string::c_str()`
///
/// If the object is a TBits, the histogram will contain the index of the bit
/// that are turned on.
///
/// ### Retrieving  information about the tree itself.
///
/// You can refer to the tree (or chain) containing the data by using the
/// string 'This'.
/// You can then could any TTree methods. For example:
/// ~~~ {.cpp}
///     tree->Draw("This->GetReadEntry()");
/// ~~~
/// will display the local entry numbers be read.
/// ~~~ {.cpp}
///     tree->Draw("This->GetUserInfo()->At(0)->GetName()");
/// ~~~
///  will display the name of the first 'user info' object.
///
/// ### Special functions and variables
///
/// `Entry$`:  A TTree::Draw formula can use the special variable `Entry$`
/// to access the entry number being read. For example to draw every
/// other entry use:
/// ~~~ {.cpp}
///     tree.Draw("myvar","Entry$%2==0");
/// ~~~
/// -  `Entry$`      : return the current entry number (`== TTree::GetReadEntry()`)
/// -  `LocalEntry$` : return the current entry number in the current tree of a
///     chain (`== GetTree()->GetReadEntry()`)
/// -  `Entries$`    : return the total number of entries (== TTree::GetEntries())
/// -  `LocalEntries$` : return the total number of entries in the current tree
///     of a chain (== GetTree()->TTree::GetEntries())
/// -  `Length$`     : return the total number of element of this formula for this
///     entry (`==TTreeFormula::GetNdata()`)
/// -  `Iteration$`  : return the current iteration over this formula for this
///     entry (i.e. varies from 0 to `Length$`).
/// -  `Length$(formula )`  : return the total number of element of the formula
///     given as a parameter.
/// -  `Sum$(formula )`  : return the sum of the value of the elements of the
///     formula given as a parameter.  For example the mean for all the elements in
///     one entry can be calculated with: `Sum$(formula )/Length$(formula )`
/// -  `Min$(formula )` : return the minimun (within one TTree entry) of the value of the
///     elements of the formula given as a parameter.
/// -  `Max$(formula )` : return the maximum (within one TTree entry) of the value of the
///     elements of the formula given as a parameter.
/// -  `MinIf$(formula,condition)`
/// -  `MaxIf$(formula,condition)` : return the minimum (maximum) (within one TTree entry)
///     of the value of the elements of the formula given as a parameter
///     if they match the condition. If no element matches the condition,
///     the result is zero.  To avoid the resulting peak at zero, use the
///     pattern:
/// ~~~ {.cpp}
///        tree->Draw("MinIf$(formula,condition)","condition");
/// ~~~
///     which will avoid calculation `MinIf$` for the entries that have no match
///     for the condition.
/// -  `Alt$(primary,alternate)` : return the value of "primary" if it is available
///     for the current iteration otherwise return the value of "alternate".
///     For example, with arr1[3] and arr2[2]
/// ~~~ {.cpp}
///        tree->Draw("arr1+Alt$(arr2,0)");
/// ~~~
///     will draw arr1[0]+arr2[0] ; arr1[1]+arr2[1] and arr1[2]+0
///     Or with a variable size array arr3
/// ~~~ {.cpp}
///        tree->Draw("Alt$(arr3[0],0)+Alt$(arr3[1],0)+Alt$(arr3[2],0)");
/// ~~~
///     will draw the sum arr3 for the index 0 to min(2,actual_size_of_arr3-1)
///     As a comparison
/// ~~~ {.cpp}
///        tree->Draw("arr3[0]+arr3[1]+arr3[2]");
/// ~~~
///     will draw the sum arr3 for the index 0 to 2 only if the
///     actual_size_of_arr3 is greater or equal to 3.
///     Note that the array in 'primary' is flattened/linearized thus using
///     `Alt$` with multi-dimensional arrays of different dimensions in unlikely
///     to yield the expected results.  To visualize a bit more what elements
///     would be matched by TTree::Draw, TTree::Scan can be used:
/// ~~~ {.cpp}
///        tree->Scan("arr1:Alt$(arr2,0)");
/// ~~~
///     will print on one line the value of arr1 and (arr2,0) that will be
///     matched by
/// ~~~ {.cpp}
///        tree->Draw("arr1-Alt$(arr2,0)");
/// ~~~
/// The ternary operator is not directly supported in TTree::Draw however, to plot the
/// equivalent of `var2<20 ? -99 : var1`, you can use:
/// ~~~ {.cpp}
///     tree->Draw("(var2<20)*99+(var2>=20)*var1","");
/// ~~~
///
/// ### Drawing a user function accessing the TTree data directly
///
/// If the formula contains  a file name, TTree::MakeProxy will be used
/// to load and execute this file.   In particular it will draw the
/// result of a function with the same name as the file.  The function
/// will be executed in a context where the name of the branches can
/// be used as a C++ variable.
///
/// For example draw px using the file hsimple.root (generated by the
/// hsimple.C tutorial), we need a file named hsimple.cxx:
/// ~~~ {.cpp}
///     double hsimple() {
///        return px;
///     }
/// ~~~
/// MakeProxy can then be used indirectly via the TTree::Draw interface
/// as follow:
/// ~~~ {.cpp}
///     new TFile("hsimple.root")
///     ntuple->Draw("hsimple.cxx");
/// ~~~
/// A more complete example is available in the tutorials directory:
/// `h1analysisProxy.cxx`, `h1analysProxy.h` and `h1analysisProxyCut.C`
/// which reimplement the selector found in `h1analysis.C`
///
/// The main features of this facility are:
///
///  * on-demand loading of branches
///  * ability to use the 'branchname' as if it was a data member
///  * protection against array out-of-bound
///  * ability to use the branch data as object (when the user code is available)
///
///  See TTree::MakeProxy for more details.
///
/// ### Making a Profile histogram
///
///  In case of a 2-Dim expression, one can generate a TProfile histogram
///  instead of a TH2F histogram by specifying option=prof or option=profs
///  or option=profi or option=profg ; the trailing letter select the way
///  the bin error are computed, See TProfile2D::SetErrorOption for
///  details on the differences.
///  The option=prof is automatically selected in case of y:x>>pf
///  where pf is an existing TProfile histogram.
///
/// ### Making a 2D Profile histogram
///
/// In case of a 3-Dim expression, one can generate a TProfile2D histogram
/// instead of a TH3F histogram by specifying option=prof or option=profs.
/// or option=profi or option=profg ; the trailing letter select the way
/// the bin error are computed, See TProfile2D::SetErrorOption for
/// details on the differences.
/// The option=prof is automatically selected in case of z:y:x>>pf
/// where pf is an existing TProfile2D histogram.
///
/// ### Making a 5D plot using GL
///
/// If option GL5D is specified together with 5 variables, a 5D plot is drawn
/// using OpenGL. See $ROOTSYS/tutorials/tree/staff.C as example.
///
/// ### Making a parallel coordinates plot
///
/// In case of a 2-Dim or more expression with the option=para, one can generate
/// a parallel coordinates plot. With that option, the number of dimensions is
/// arbitrary. Giving more than 4 variables without the option=para or
/// option=candle or option=goff will produce an error.
///
/// ### Making a candle sticks chart
///
/// In case of a 2-Dim or more expression with the option=candle, one can generate
/// a candle sticks chart. With that option, the number of dimensions is
/// arbitrary. Giving more than 4 variables without the option=para or
/// option=candle or option=goff will produce an error.
///
/// ### Normalizing the output histogram to 1
///
/// When option contains "norm" the output histogram is normalized to 1.
///
/// ### Saving the result of Draw to a TEventList, a TEntryList or a TEntryListArray
///
/// TTree::Draw can be used to fill a TEventList object (list of entry numbers)
/// instead of histogramming one variable.
/// If varexp0 has the form >>elist , a TEventList object named "elist"
/// is created in the current directory. elist will contain the list
/// of entry numbers satisfying the current selection.
/// If option "entrylist" is used, a TEntryList object is created
/// If the selection contains arrays, vectors or any container class and option
/// "entrylistarray" is used, a TEntryListArray object is created
/// containing also the subentries satisfying the selection, i.e. the indices of
/// the branches which hold containers classes.
/// Example:
/// ~~~ {.cpp}
///     tree.Draw(">>yplus","y>0")
/// ~~~
/// will create a TEventList object named "yplus" in the current directory.
/// In an interactive session, one can type (after TTree::Draw)
/// ~~~ {.cpp}
///     yplus.Print("all")
/// ~~~
/// to print the list of entry numbers in the list.
/// ~~~ {.cpp}
///     tree.Draw(">>yplus", "y>0", "entrylist")
/// ~~~
/// will create a TEntryList object names "yplus" in the current directory
/// ~~~ {.cpp}
///     tree.Draw(">>yplus", "y>0", "entrylistarray")
/// ~~~
/// will create a TEntryListArray object names "yplus" in the current directory
///
/// By default, the specified entry list is reset.
/// To continue to append data to an existing list, use "+" in front
/// of the list name;
/// ~~~ {.cpp}
///     tree.Draw(">>+yplus","y>0")
/// ~~~
/// will not reset yplus, but will enter the selected entries at the end
/// of the existing list.
///
/// ### Using a TEventList, TEntryList or TEntryListArray as Input
///
/// Once a TEventList or a TEntryList object has been generated, it can be used as input
/// for TTree::Draw. Use TTree::SetEventList or TTree::SetEntryList to set the
/// current event list
///
/// Example 1:
/// ~~~ {.cpp}
///     TEventList *elist = (TEventList*)gDirectory->Get("yplus");
///     tree->SetEventList(elist);
///     tree->Draw("py");
/// ~~~
/// Example 2:
/// ~~~ {.cpp}
///     TEntryList *elist = (TEntryList*)gDirectory->Get("yplus");
///     tree->SetEntryList(elist);
///     tree->Draw("py");
/// ~~~
/// If a TEventList object is used as input, a new TEntryList object is created
/// inside the SetEventList function. In case of a TChain, all tree headers are loaded
/// for this transformation. This new object is owned by the chain and is deleted
/// with it, unless the user extracts it by calling GetEntryList() function.
/// See also comments to SetEventList() function of TTree and TChain.
///
/// If arrays are used in the selection criteria and TEntryListArray is not used,
/// all the entries that have at least one element of the array that satisfy the selection
/// are entered in the list.
///
/// Example:
/// ~~~ {.cpp}
///     tree.Draw(">>pyplus","fTracks.fPy>0");
///     tree->SetEventList(pyplus);
///     tree->Draw("fTracks.fPy");
/// ~~~
///  will draw the fPy of ALL tracks in event with at least one track with
///  a positive fPy.
///
/// To select only the elements that did match the original selection
/// use TEventList::SetReapplyCut or TEntryList::SetReapplyCut.
///
/// Example:
/// ~~~ {.cpp}
///     tree.Draw(">>pyplus","fTracks.fPy>0");
///     pyplus->SetReapplyCut(kTRUE);
///     tree->SetEventList(pyplus);
///     tree->Draw("fTracks.fPy");
/// ~~~
/// will draw the fPy of only the tracks that have a positive fPy.
///
/// To draw only the elements that match a selection in case of arrays,
/// you can also use TEntryListArray (faster in case of a more general selection).
///
/// Example:
/// ~~~ {.cpp}
///     tree.Draw(">>pyplus","fTracks.fPy>0", "entrylistarray");
///     tree->SetEntryList(pyplus);
///     tree->Draw("fTracks.fPy");
/// ~~~
/// will draw the fPy of only the tracks that have a positive fPy,
/// but without redoing the selection.
///
///  Note: Use tree->SetEventList(0) if you do not want use the list as input.
///
/// ### How to obtain more info from TTree::Draw
///
///  Once TTree::Draw has been called, it is possible to access useful
///  information still stored in the TTree object via the following functions:
///
/// - GetSelectedRows() // return the number of values accepted by the selection expression. In case where no selection was specified, returns the number of values processed.
/// - GetV1()           // returns a pointer to the double array of V1
/// - GetV2()           // returns a pointer to the double array of V2
/// - GetV3()           // returns a pointer to the double array of V3
/// - GetV4()           // returns a pointer to the double array of V4
/// - GetW()            // returns a pointer to the double array of Weights where weight equal the result of the selection expression.
///
/// where V1,V2,V3 correspond to the expressions in
/// ~~~ {.cpp}
///     TTree::Draw("V1:V2:V3:V4",selection);
/// ~~~
/// If the expression has more than 4 component use GetVal(index)
///
/// Example:
/// ~~~ {.cpp}
///     Root > ntuple->Draw("py:px","pz>4");
///     Root > TGraph *gr = new TGraph(ntuple->GetSelectedRows(),
///                                   ntuple->GetV2(), ntuple->GetV1());
///     Root > gr->Draw("ap"); //draw graph in current pad
/// ~~~
///
/// A more complete complete tutorial (treegetval.C) shows how to use the
/// GetVal() method.
///
/// creates a TGraph object with a number of points corresponding to the
/// number of entries selected by the expression "pz>4", the x points of the graph
/// being the px values of the Tree and the y points the py values.
///
/// Important note: By default TTree::Draw creates the arrays obtained
/// with GetW, GetV1, GetV2, GetV3, GetV4, GetVal with a length corresponding
/// to the parameter fEstimate.  The content will be the last `GetSelectedRows() % GetEstimate()`
/// values calculated.
/// By default fEstimate=1000000 and can be modified
/// via TTree::SetEstimate. To keep in memory all the results (in case
/// where there is only one result per entry), use
/// ~~~ {.cpp}
///     tree->SetEstimate(tree->GetEntries()+1); // same as tree->SetEstimate(-1);
/// ~~~
/// You must call SetEstimate if the expected number of selected rows
/// you need to look at is greater than 1000000.
///
/// You can use the option "goff" to turn off the graphics output
/// of TTree::Draw in the above example.
///
/// ### Automatic interface to TTree::Draw via the TTreeViewer
///
/// A complete graphical interface to this function is implemented
/// in the class TTreeViewer.
/// To start the TTreeViewer, three possibilities:
/// - select TTree context menu item "StartViewer"
/// - type the command  "TTreeViewer TV(treeName)"
/// - execute statement "tree->StartViewer();"

Long64_t TTree::Draw(const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer)
      return fPlayer->DrawSelect(varexp,selection,option,nentries,firstentry);
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove some baskets from memory.

void TTree::DropBaskets()
{
   TBranch* branch = 0;
   Int_t nb = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i) {
      branch = (TBranch*) fBranches.UncheckedAt(i);
      branch->DropBaskets("all");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Drop branch buffers to accommodate nbytes below MaxVirtualsize.

void TTree::DropBuffers(Int_t)
{
   // Be careful not to remove current read/write buffers.
   Int_t ndrop = 0;
   Int_t nleaves = fLeaves.GetEntriesFast();
   for (Int_t i = 0; i < nleaves; ++i)  {
      TLeaf* leaf = (TLeaf*) fLeaves.UncheckedAt(i);
      TBranch* branch = (TBranch*) leaf->GetBranch();
      Int_t nbaskets = branch->GetListOfBaskets()->GetEntries();
      for (Int_t j = 0; j < nbaskets - 1; ++j) {
         if ((j == branch->GetReadBasket()) || (j == branch->GetWriteBasket())) {
            continue;
         }
         TBasket* basket = (TBasket*)branch->GetListOfBaskets()->UncheckedAt(j);
         if (basket) {
            ndrop += basket->DropBuffers();
            if (fTotalBuffers < fMaxVirtualSize) {
               return;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill all branches.
///
/// This function loops on all the branches of this tree.  For
/// each branch, it copies to the branch buffer (basket) the current
/// values of the leaves data types. If a leaf is a simple data type,
/// a simple conversion to a machine independent format has to be done.
///
/// This machine independent version of the data is copied into a
/// basket (each branch has its own basket).  When a basket is full
/// (32k worth of data by default), it is then optionally compressed
/// and written to disk (this operation is also called committing or
/// 'flushing' the basket).  The committed baskets are then
/// immediately removed from memory.
///
/// The function returns the number of bytes committed to the
/// individual branches.
///
/// If a write error occurs, the number of bytes returned is -1.
///
/// If no data are written, because, e.g., the branch is disabled,
/// the number of bytes returned is 0.
///
/// __The baskets are flushed and the Tree header saved at regular intervals__
///
/// At regular intervals, when the amount of data written so far is
/// greater than fAutoFlush (see SetAutoFlush) all the baskets are flushed to disk.
/// This makes future reading faster as it guarantees that baskets belonging to nearby
/// entries will be on the same disk region.
/// When the first call to flush the baskets happen, we also take this opportunity
/// to optimize the baskets buffers.
/// We also check if the amount of data written is greater than fAutoSave (see SetAutoSave).
/// In this case we also write the Tree header. This makes the Tree recoverable up to this point
/// in case the program writing the Tree crashes.
/// The decisions to FlushBaskets and Auto Save can be made based either on the number
/// of bytes written (fAutoFlush and fAutoSave negative) or on the number of entries
/// written (fAutoFlush and fAutoSave positive).
/// Note that the user can decide to call FlushBaskets and AutoSave in her event loop
/// base on the number of events written instead of the number of bytes written.
///
/// \note Calling `TTree::FlushBaskets` too often increases the IO time. 
/// 
/// \note Calling `TTree::AutoSave` too often increases the IO time and also the
///       file size.
///
/// \note This method calls `TTree::ChangeFile` when the tree reaches a size
///       greater than `TTree::fgMaxTreeSize`. This doesn't happen if the tree is
///       attached to a `TMemFile` or derivate.

Int_t TTree::Fill()
{
   Int_t nbytes = 0;
   Int_t nwrite = 0;
   Int_t nerror = 0;
   Int_t nbranches = fBranches.GetEntriesFast();

   // Case of one single super branch. Automatically update
   // all the branch addresses if a new object was created.
   if (nbranches == 1)
      ((TBranch *)fBranches.UncheckedAt(0))->UpdateAddress();

   if (fBranchRef)
      fBranchRef->Clear();

#ifdef R__USE_IMT
   const auto useIMT = ROOT::IsImplicitMTEnabled() && fIMTEnabled;
   ROOT::Internal::TBranchIMTHelper imtHelper;
   if (useIMT) {
      fIMTFlush = true;
      fIMTZipBytes.store(0);
      fIMTTotBytes.store(0);
   }
#endif

   for (Int_t i = 0; i < nbranches; ++i) {
      // Loop over all branches, filling and accumulating bytes written and error counts.
      TBranch *branch = (TBranch *)fBranches.UncheckedAt(i);

      if (branch->TestBit(kDoNotProcess))
         continue;

#ifndef R__USE_IMT
      nwrite = branch->FillImpl(nullptr);
#else
      nwrite = branch->FillImpl(useIMT ? &imtHelper : nullptr);
#endif
      if (nwrite < 0) {
         if (nerror < 2) {
            Error("Fill", "Failed filling branch:%s.%s, nbytes=%d, entry=%lld\n"
                          " This error is symptomatic of a Tree created as a memory-resident Tree\n"
                          " Instead of doing:\n"
                          "    TTree *T = new TTree(...)\n"
                          "    TFile *f = new TFile(...)\n"
                          " you should do:\n"
                          "    TFile *f = new TFile(...)\n"
                          "    TTree *T = new TTree(...)\n\n",
                  GetName(), branch->GetName(), nwrite, fEntries + 1);
         } else {
            Error("Fill", "Failed filling branch:%s.%s, nbytes=%d, entry=%lld", GetName(), branch->GetName(), nwrite,
                  fEntries + 1);
         }
         ++nerror;
      } else {
         nbytes += nwrite;
      }
   }

#ifdef R__USE_IMT
   if (fIMTFlush) {
      imtHelper.Wait();
      fIMTFlush = false;
      const_cast<TTree *>(this)->AddTotBytes(fIMTTotBytes);
      const_cast<TTree *>(this)->AddZipBytes(fIMTZipBytes);
      nbytes += imtHelper.GetNbytes();
      nerror += imtHelper.GetNerrors();
   }
#endif

   if (fBranchRef)
      fBranchRef->Fill();

   ++fEntries;

   if (fEntries > fMaxEntries)
      KeepCircular();

   if (gDebug > 0)
      Info("TTree::Fill", " - A: %d %lld %lld %lld %lld %lld %lld \n", nbytes, fEntries, fAutoFlush, fAutoSave,
           GetZipBytes(), fFlushedBytes, fSavedBytes);

   bool autoFlush = false;
   bool autoSave = false;

   if (fAutoFlush != 0 || fAutoSave != 0) {
      // Is it time to flush or autosave baskets?
      if (fFlushedBytes == 0) {
         // If fFlushedBytes == 0, it means we never flushed or saved, so
         // we need to check if it's time to do it and recompute the values
         // of fAutoFlush and fAutoSave in terms of the number of entries.
         // Decision can be based initially either on the number of bytes
         // or the number of entries written.
         Long64_t zipBytes = GetZipBytes();

         if (fAutoFlush)
            autoFlush = fAutoFlush < 0 ? (zipBytes > -fAutoFlush) : fEntries % fAutoFlush == 0;

         if (fAutoSave)
            autoSave = fAutoSave < 0 ? (zipBytes > -fAutoSave) : fEntries % fAutoSave == 0;

         if (autoFlush || autoSave) {
            // First call FlushBasket to make sure that fTotBytes is up to date.
            FlushBasketsImpl();
            autoFlush = false; // avoid auto flushing again later

            // When we are in one-basket-per-cluster mode, there is no need to optimize basket:
            // they will automatically grow to the size needed for an event cluster (with the basket
            // shrinking preventing them from growing too much larger than the actually-used space).
            if (!TestBit(TTree::kOnlyFlushAtCluster)) {
               OptimizeBaskets(GetTotBytes(), 1, "");
               if (gDebug > 0)
                  Info("TTree::Fill", "OptimizeBaskets called at entry %lld, fZipBytes=%lld, fFlushedBytes=%lld\n",
                       fEntries, GetZipBytes(), fFlushedBytes);
            }
            fFlushedBytes = GetZipBytes();
            fAutoFlush = fEntries; // Use test on entries rather than bytes

            // subsequently in run
            if (fAutoSave < 0) {
               // Set fAutoSave to the largest integer multiple of
               // fAutoFlush events such that fAutoSave*fFlushedBytes
               // < (minus the input value of fAutoSave)
               Long64_t totBytes = GetTotBytes();
               if (zipBytes != 0) {
                  fAutoSave = TMath::Max(fAutoFlush, fEntries * ((-fAutoSave / zipBytes) / fEntries));
               } else if (totBytes != 0) {
                  fAutoSave = TMath::Max(fAutoFlush, fEntries * ((-fAutoSave / totBytes) / fEntries));
               } else {
                  TBufferFile b(TBuffer::kWrite, 10000);
                  TTree::Class()->WriteBuffer(b, (TTree *)this);
                  Long64_t total = b.Length();
                  fAutoSave = TMath::Max(fAutoFlush, fEntries * ((-fAutoSave / total) / fEntries));
               }
            } else if (fAutoSave > 0) {
               fAutoSave = fAutoFlush * (fAutoSave / fAutoFlush);
            }

            if (fAutoSave != 0 && fEntries >= fAutoSave)
               autoSave = true;

            if (gDebug > 0)
               Info("TTree::Fill", "First AutoFlush.  fAutoFlush = %lld, fAutoSave = %lld\n", fAutoFlush, fAutoSave);
         }
      } else {
         // Check if we need to auto flush
         if (fAutoFlush) {
            if (fNClusterRange == 0)
               autoFlush = fEntries > 1 && fEntries % fAutoFlush == 0;
            else
               autoFlush = (fEntries - (fClusterRangeEnd[fNClusterRange - 1] + 1)) % fAutoFlush == 0;
         }
         // Check if we need to auto save
         if (fAutoSave)
            autoSave = fEntries % fAutoSave == 0;
      }
   }

   if (autoFlush) {
      FlushBasketsImpl();
      if (gDebug > 0)
         Info("TTree::Fill", "FlushBaskets() called at entry %lld, fZipBytes=%lld, fFlushedBytes=%lld\n", fEntries,
              GetZipBytes(), fFlushedBytes);
      fFlushedBytes = GetZipBytes();
   }

   if (autoSave) {
      AutoSave(); // does not call FlushBasketsImpl() again
      if (gDebug > 0)
         Info("TTree::Fill", "AutoSave called at entry %lld, fZipBytes=%lld, fSavedBytes=%lld\n", fEntries,
              GetZipBytes(), fSavedBytes);
   }

   // Check that output file is still below the maximum size.
   // If above, close the current file and continue on a new file.
   // Currently, the automatic change of file is restricted
   // to the case where the tree is in the top level directory.
   if (fDirectory)
      if (TFile *file = fDirectory->GetFile())
         if (static_cast<TDirectory *>(file) == fDirectory && (file->GetEND() > fgMaxTreeSize))
            // Changing file clashes with the design of TMemFile and derivates, see #6523.
            if (!(dynamic_cast<TMemFile *>(file)))
               ChangeFile(file);

   return nerror == 0 ? nbytes : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Search in the array for a branch matching the branch name,
/// with the branch possibly expressed as a 'full' path name (with dots).

static TBranch *R__FindBranchHelper(TObjArray *list, const char *branchname) {
   if (list==0 || branchname == 0 || branchname[0] == '\0') return 0;

   Int_t nbranches = list->GetEntries();

   UInt_t brlen = strlen(branchname);

   for(Int_t index = 0; index < nbranches; ++index) {
      TBranch *where = (TBranch*)list->UncheckedAt(index);

      const char *name = where->GetName();
      UInt_t len = strlen(name);
      if (len && name[len-1]==']') {
         const  char *dim = strchr(name,'[');
         if (dim) {
            len = dim - name;
         }
      }
      if (brlen == len && strncmp(branchname,name,len)==0) {
         return where;
      }
      TBranch *next = 0;
      if ((brlen >= len) && (branchname[len] == '.')
          && strncmp(name, branchname, len) == 0) {
         // The prefix subbranch name match the branch name.

         next = where->FindBranch(branchname);
         if (!next) {
            next = where->FindBranch(branchname+len+1);
         }
         if (next) return next;
      }
      const char *dot = strchr((char*)branchname,'.');
      if (dot) {
         if (len==(size_t)(dot-branchname) &&
             strncmp(branchname,name,dot-branchname)==0 ) {
            return R__FindBranchHelper(where->GetListOfBranches(),dot+1);
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the branch that correspond to the path 'branchname', which can
/// include the name of the tree or the omitted name of the parent branches.
/// In case of ambiguity, returns the first match.

TBranch* TTree::FindBranch(const char* branchname)
{
   // We already have been visited while recursively looking
   // through the friends tree, let return
   if (kFindBranch & fFriendLockStatus) {
      return 0;
   }

   TBranch* branch = 0;
   // If the first part of the name match the TTree name, look for the right part in the
   // list of branches.
   // This will allow the branchname to be preceded by
   // the name of this tree.
   if (strncmp(fName.Data(),branchname,fName.Length())==0 && branchname[fName.Length()]=='.') {
      branch = R__FindBranchHelper( GetListOfBranches(), branchname + fName.Length() + 1);
      if (branch) return branch;
   }
   // If we did not find it, let's try to find the full name in the list of branches.
   branch = R__FindBranchHelper(GetListOfBranches(), branchname);
   if (branch) return branch;

   // If we still did not find, let's try to find it within each branch assuming it does not the branch name.
   TIter next(GetListOfBranches());
   while ((branch = (TBranch*) next())) {
      TBranch* nestedbranch = branch->FindBranch(branchname);
      if (nestedbranch) {
         return nestedbranch;
      }
   }

   // Search in list of friends.
   if (!fFriends) {
      return 0;
   }
   TFriendLock lock(this, kFindBranch);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree* t = fe->GetTree();
      if (!t) {
         continue;
      }
      // If the alias is present replace it with the real name.
      const char *subbranch = strstr(branchname, fe->GetName());
      if (subbranch != branchname) {
         subbranch = 0;
      }
      if (subbranch) {
         subbranch += strlen(fe->GetName());
         if (*subbranch != '.') {
            subbranch = 0;
         } else {
            ++subbranch;
         }
      }
      std::ostringstream name;
      if (subbranch) {
         name << t->GetName() << "." << subbranch;
      } else {
         name << branchname;
      }
      branch = t->FindBranch(name.str().c_str());
      if (branch) {
         return branch;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find leaf..

TLeaf* TTree::FindLeaf(const char* searchname)
{
   // We already have been visited while recursively looking
   // through the friends tree, let's return.
   if (kFindLeaf & fFriendLockStatus) {
      return 0;
   }

   // This will allow the branchname to be preceded by
   // the name of this tree.
   char* subsearchname = (char*) strstr(searchname, GetName());
   if (subsearchname != searchname) {
      subsearchname = 0;
   }
   if (subsearchname) {
      subsearchname += strlen(GetName());
      if (*subsearchname != '.') {
         subsearchname = 0;
      } else {
         ++subsearchname;
         if (subsearchname[0]==0) {
            subsearchname = 0;
         }
      }
   }

   TString leafname;
   TString leaftitle;
   TString longname;
   TString longtitle;

   const bool searchnameHasDot = strchr(searchname, '.') != nullptr;

   // For leaves we allow for one level up to be prefixed to the name.
   TIter next(GetListOfLeaves());
   TLeaf* leaf = 0;
   while ((leaf = (TLeaf*) next())) {
      leafname = leaf->GetName();
      Ssiz_t dim = leafname.First('[');
      if (dim >= 0) leafname.Remove(dim);

      if (leafname == searchname) {
         return leaf;
      }
      if (subsearchname && leafname == subsearchname) {
         return leaf;
      }
      // The TLeafElement contains the branch name
      // in its name, let's use the title.
      leaftitle = leaf->GetTitle();
      dim = leaftitle.First('[');
      if (dim >= 0) leaftitle.Remove(dim);

      if (leaftitle == searchname) {
         return leaf;
      }
      if (subsearchname && leaftitle == subsearchname) {
         return leaf;
      }
      if (!searchnameHasDot)
         continue;
      TBranch* branch = leaf->GetBranch();
      if (branch) {
         longname.Form("%s.%s",branch->GetName(),leafname.Data());
         dim = longname.First('[');
         if (dim>=0) longname.Remove(dim);
         if (longname == searchname) {
            return leaf;
         }
         if (subsearchname && longname == subsearchname) {
            return leaf;
         }
         longtitle.Form("%s.%s",branch->GetName(),leaftitle.Data());
         dim = longtitle.First('[');
         if (dim>=0) longtitle.Remove(dim);
         if (longtitle == searchname) {
            return leaf;
         }
         if (subsearchname && longtitle == subsearchname) {
            return leaf;
         }
         // The following is for the case where the branch is only
         // a sub-branch.  Since we do not see it through
         // TTree::GetListOfBranches, we need to see it indirectly.
         // This is the less sturdy part of this search ... it may
         // need refining ...
         if (strstr(searchname, ".") && !strcmp(searchname, branch->GetName())) {
            return leaf;
         }
         if (subsearchname && strstr(subsearchname, ".") && !strcmp(subsearchname, branch->GetName())) {
            return leaf;
         }
      }
   }
   // Search in list of friends.
   if (!fFriends) {
      return 0;
   }
   TFriendLock lock(this, kFindLeaf);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree* t = fe->GetTree();
      if (!t) {
         continue;
      }
      // If the alias is present replace it with the real name.
      subsearchname = (char*) strstr(searchname, fe->GetName());
      if (subsearchname != searchname) {
         subsearchname = 0;
      }
      if (subsearchname) {
         subsearchname += strlen(fe->GetName());
         if (*subsearchname != '.') {
            subsearchname = 0;
         } else {
            ++subsearchname;
         }
      }
      if (subsearchname) {
         leafname.Form("%s.%s",t->GetName(),subsearchname);
      } else {
         leafname = searchname;
      }
      leaf = t->FindLeaf(leafname);
      if (leaf) {
         return leaf;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fit  a projected item(s) from a tree.
///
/// funcname is a TF1 function.
///
/// See TTree::Draw() for explanations of the other parameters.
///
/// By default the temporary histogram created is called htemp.
/// If varexp contains >>hnew , the new histogram created is called hnew
/// and it is kept in the current directory.
///
/// The function returns the number of selected entries.
///
/// Example:
/// ~~~ {.cpp}
///     tree.Fit(pol4,"sqrt(x)>>hsqrt","y>0")
/// ~~~
/// will fit sqrt(x) and save the histogram as "hsqrt" in the current
/// directory.
///
/// See also TTree::UnbinnedFit
///
/// ## Return status
///
///  The function returns the status of the histogram fit (see TH1::Fit)
///  If no entries were selected, the function returns -1;
///   (i.e. fitResult is null if the fit is OK)

Int_t TTree::Fit(const char* funcname, const char* varexp, const char* selection, Option_t* option, Option_t* goption, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->Fit(funcname, varexp, selection, option, goption, nentries, firstentry);
   }
   return -1;
}

namespace {
struct BoolRAIIToggle {
   Bool_t &m_val;

   BoolRAIIToggle(Bool_t &val) : m_val(val) { m_val = true; }
   ~BoolRAIIToggle() { m_val = false; }
};
}

////////////////////////////////////////////////////////////////////////////////
/// Write to disk all the basket that have not yet been individually written and
/// create an event cluster boundary (by default).
///
/// If the caller wishes to flush the baskets but not create an event cluster,
/// then set create_cluster to false.
///
/// If ROOT has IMT-mode enabled, this will launch multiple TBB tasks in parallel
/// via TThreadExecutor to do this operation; one per basket compression.  If the
///  caller utilizes TBB also, care must be taken to prevent deadlocks.
///
/// For example, let's say the caller holds mutex A and calls FlushBaskets; while
/// TBB is waiting for the ROOT compression tasks to complete, it may decide to
/// run another one of the user's tasks in this thread.  If the second user task
/// tries to acquire A, then a deadlock will occur.  The example call sequence
/// looks like this:
///
/// - User acquires mutex A
/// - User calls FlushBaskets.
/// - ROOT launches N tasks and calls wait.
/// - TBB schedules another user task, T2.
/// - T2 tries to acquire mutex A.
///
/// At this point, the thread will deadlock: the code may function with IMT-mode
/// disabled if the user assumed the legacy code never would run their own TBB
/// tasks.
///
/// SO: users of TBB who want to enable IMT-mode should carefully review their
/// locking patterns and make sure they hold no coarse-grained application
/// locks when they invoke ROOT.
///
/// Return the number of bytes written or -1 in case of write error.
Int_t TTree::FlushBaskets(Bool_t create_cluster) const
{
    Int_t retval = FlushBasketsImpl();
    if (retval == -1) return retval;

    if (create_cluster) const_cast<TTree *>(this)->MarkEventCluster();
    return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal implementation of the FlushBaskets algorithm.
/// Unlike the public interface, this does NOT create an explicit event cluster
/// boundary; it is up to the (internal) caller to determine whether that should
/// done.
///
/// Otherwise, the comments for FlushBaskets applies.
///
Int_t TTree::FlushBasketsImpl() const
{
   if (!fDirectory) return 0;
   Int_t nbytes = 0;
   Int_t nerror = 0;
   TObjArray *lb = const_cast<TTree*>(this)->GetListOfBranches();
   Int_t nb = lb->GetEntriesFast();

#ifdef R__USE_IMT
   const auto useIMT = ROOT::IsImplicitMTEnabled() && fIMTEnabled;
   if (useIMT) {
      // ROOT-9668: here we need to check if the size of fSortedBranches is different from the
      // size of the list of branches before triggering the initialisation of the fSortedBranches
      // container to cover two cases:
      // 1. This is the first time we flush. fSortedBranches is empty and we need to fill it.
      // 2. We flushed at least once already but a branch has been be added to the tree since then
      if (fSortedBranches.size() != unsigned(nb)) { const_cast<TTree*>(this)->InitializeBranchLists(false); }

      BoolRAIIToggle sentry(fIMTFlush);
      fIMTZipBytes.store(0);
      fIMTTotBytes.store(0);
      std::atomic<Int_t> nerrpar(0);
      std::atomic<Int_t> nbpar(0);
      std::atomic<Int_t> pos(0);

      auto mapFunction  = [&]() {
        // The branch to process is obtained when the task starts to run.
        // This way, since branches are sorted, we make sure that branches
        // leading to big tasks are processed first. If we assigned the
        // branch at task creation time, the scheduler would not necessarily
        // respect our sorting.
        Int_t j = pos.fetch_add(1);

        auto branch = fSortedBranches[j].second;
        if (R__unlikely(!branch)) { return; }

        if (R__unlikely(gDebug > 0)) {
            std::stringstream ss;
            ss << std::this_thread::get_id();
            Info("FlushBaskets", "[IMT] Thread %s", ss.str().c_str());
            Info("FlushBaskets", "[IMT] Running task for branch #%d: %s", j, branch->GetName());
        }

        Int_t nbtask = branch->FlushBaskets();

        if (nbtask < 0) { nerrpar++; }
        else            { nbpar += nbtask; }
      };

      ROOT::TThreadExecutor pool;
      pool.Foreach(mapFunction, nb);

      fIMTFlush = false;
      const_cast<TTree*>(this)->AddTotBytes(fIMTTotBytes);
      const_cast<TTree*>(this)->AddZipBytes(fIMTZipBytes);

      return nerrpar ? -1 : nbpar.load();
   }
#endif
   for (Int_t j = 0; j < nb; j++) {
      TBranch* branch = (TBranch*) lb->UncheckedAt(j);
      if (branch) {
         Int_t nwrite = branch->FlushBaskets();
         if (nwrite<0) {
            ++nerror;
         } else {
            nbytes += nwrite;
         }
      }
   }
   if (nerror) {
      return -1;
   } else {
      return nbytes;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the expanded value of the alias.  Search in the friends if any.

const char* TTree::GetAlias(const char* aliasName) const
{
   // We already have been visited while recursively looking
   // through the friends tree, let's return.
   if (kGetAlias & fFriendLockStatus) {
      return 0;
   }
   if (fAliases) {
      TObject* alias = fAliases->FindObject(aliasName);
      if (alias) {
         return alias->GetTitle();
      }
   }
   if (!fFriends) {
      return 0;
   }
   TFriendLock lock(const_cast<TTree*>(this), kGetAlias);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree* t = fe->GetTree();
      if (t) {
         const char* alias = t->GetAlias(aliasName);
         if (alias) {
            return alias;
         }
         const char* subAliasName = strstr(aliasName, fe->GetName());
         if (subAliasName && (subAliasName[strlen(fe->GetName())] == '.')) {
            alias = t->GetAlias(aliasName + strlen(fe->GetName()) + 1);
            if (alias) {
               return alias;
            }
         }
      }
   }
   return 0;
}

namespace {
/// Do a breadth first search through the implied hierarchy
/// of branches.
/// To avoid scanning through the list multiple time
/// we also remember the 'depth-first' match.
TBranch *R__GetBranch(const TObjArray &branches, const char *name)
{
   TBranch *result = nullptr;
   Int_t nb = branches.GetEntriesFast();
   for (Int_t i = 0; i < nb; i++) {
      TBranch* b = (TBranch*)branches.UncheckedAt(i);
      if (!b)
          continue;
      if (!strcmp(b->GetName(), name)) {
         return b;
      }
      if (!strcmp(b->GetFullName(), name)) {
         return b;
      }
      if (!result)
         result = R__GetBranch(*(b->GetListOfBranches()), name);
   }
   return result;
}
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the branch with the given name in this tree or its friends.
/// The search is done breadth first.

TBranch* TTree::GetBranch(const char* name)
{
   if (name == 0) return 0;

   // We already have been visited while recursively
   // looking through the friends tree, let's return.
   if (kGetBranch & fFriendLockStatus) {
      return 0;
   }

   // Search using branches, breadth first.
   TBranch *result = R__GetBranch(fBranches, name);
   if (result)
     return result;

   // Search using leaves.
   TObjArray* leaves = GetListOfLeaves();
   Int_t nleaves = leaves->GetEntriesFast();
   for (Int_t i = 0; i < nleaves; i++) {
      TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(i);
      TBranch* branch = leaf->GetBranch();
      if (!strcmp(branch->GetName(), name)) {
         return branch;
      }
      if (!strcmp(branch->GetFullName(), name)) {
         return branch;
      }
   }

   if (!fFriends) {
      return 0;
   }

   // Search in list of friends.
   TFriendLock lock(this, kGetBranch);
   TIter next(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) next())) {
      TTree* t = fe->GetTree();
      if (t) {
         TBranch* branch = t->GetBranch(name);
         if (branch) {
            return branch;
         }
      }
   }

   // Second pass in the list of friends when
   // the branch name is prefixed by the tree name.
   next.Reset();
   while ((fe = (TFriendElement*) next())) {
      TTree* t = fe->GetTree();
      if (!t) {
         continue;
      }
      char* subname = (char*) strstr(name, fe->GetName());
      if (subname != name) {
         continue;
      }
      Int_t l = strlen(fe->GetName());
      subname += l;
      if (*subname != '.') {
         continue;
      }
      subname++;
      TBranch* branch = t->GetBranch(subname);
      if (branch) {
         return branch;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return status of branch with name branchname.
///
/// - 0 if branch is not activated
/// - 1 if branch is activated

Bool_t TTree::GetBranchStatus(const char* branchname) const
{
   TBranch* br = const_cast<TTree*>(this)->GetBranch(branchname);
   if (br) {
      return br->TestBit(kDoNotProcess) == 0;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning the current branch style.
///
/// - style = 0 old Branch
/// - style = 1 new Bronch

Int_t TTree::GetBranchStyle()
{
   return fgBranchStyle;
}

////////////////////////////////////////////////////////////////////////////////
/// Used for automatic sizing of the cache.
///
/// Estimates a suitable size for the tree cache based on AutoFlush.
/// A cache sizing factor is taken from the configuration. If this yields zero
/// and withDefault is true the historical algorithm for default size is used.

Long64_t TTree::GetCacheAutoSize(Bool_t withDefault /* = kFALSE */ ) const
{
   const char *stcs;
   Double_t cacheFactor = 0.0;
   if (!(stcs = gSystem->Getenv("ROOT_TTREECACHE_SIZE")) || !*stcs) {
      cacheFactor = gEnv->GetValue("TTreeCache.Size", 1.0);
   } else {
      cacheFactor = TString(stcs).Atof();
   }

   if (cacheFactor < 0.0) {
     // ignore negative factors
     cacheFactor = 0.0;
   }

   Long64_t cacheSize = 0;

   if (fAutoFlush < 0) cacheSize = Long64_t(-cacheFactor*fAutoFlush);
   else if (fAutoFlush == 0) cacheSize = 0;
   else cacheSize = Long64_t(cacheFactor*1.5*fAutoFlush*GetZipBytes()/(fEntries+1));

   if (cacheSize >= (INT_MAX / 4)) {
      cacheSize = INT_MAX / 4;
   }

   if (cacheSize < 0) {
      cacheSize = 0;
   }

   if (cacheSize == 0 && withDefault) {
      if (fAutoFlush < 0) cacheSize = -fAutoFlush;
      else if (fAutoFlush == 0) cacheSize = 0;
      else cacheSize = Long64_t(1.5*fAutoFlush*GetZipBytes()/(fEntries+1));
   }

   return cacheSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return an iterator over the cluster of baskets starting at firstentry.
///
/// This iterator is not yet supported for TChain object.
/// ~~~ {.cpp}
///      TTree::TClusterIterator clusterIter = tree->GetClusterIterator(entry);
///      Long64_t clusterStart;
///      while( (clusterStart = clusterIter()) < tree->GetEntries() ) {
///         printf("The cluster starts at %lld and ends at %lld (inclusive)\n",clusterStart,clusterIter.GetNextEntry()-1);
///      }
/// ~~~

TTree::TClusterIterator TTree::GetClusterIterator(Long64_t firstentry)
{
   // create cache if wanted
   if (fCacheDoAutoInit)
      SetCacheSizeAux();

   return TClusterIterator(this,firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the current file.

TFile* TTree::GetCurrentFile() const
{
   if (!fDirectory || fDirectory==gROOT) {
      return 0;
   }
   return fDirectory->GetFile();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of entries matching the selection.
/// Return -1 in case of errors.
///
/// If the selection uses any arrays or containers, we return the number
/// of entries where at least one element match the selection.
/// GetEntries is implemented using the selector class TSelectorEntries,
/// which can be used directly (see code in TTreePlayer::GetEntries) for
/// additional option.
/// If SetEventList was used on the TTree or TChain, only that subset
/// of entries will be considered.

Long64_t TTree::GetEntries(const char *selection)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->GetEntries(selection);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the 1st Leaf named name in any Branch of this Tree or
/// any branch in the list of friend trees.

Long64_t TTree::GetEntriesFriend() const
{
   if (fEntries) return fEntries;
   if (!fFriends) return 0;
   TFriendElement *fr = (TFriendElement*)fFriends->At(0);
   if (!fr) return 0;
   TTree *t = fr->GetTree();
   if (t==0) return 0;
   return t->GetEntriesFriend();
}

////////////////////////////////////////////////////////////////////////////////
/// Read all branches of entry and return total number of bytes read.
///
/// - `getall = 0` : get only active branches
/// - `getall = 1` : get all branches
///
/// The function returns the number of bytes read from the input buffer.
/// If entry does not exist the function returns 0.
/// If an I/O error occurs, the function returns -1.
///
/// If the Tree has friends, also read the friends entry.
///
/// To activate/deactivate one or more branches, use TBranch::SetBranchStatus
/// For example, if you have a Tree with several hundred branches, and you
/// are interested only by branches named "a" and "b", do
/// ~~~ {.cpp}
///     mytree.SetBranchStatus("*",0); //disable all branches
///     mytree.SetBranchStatus("a",1);
///     mytree.SetBranchStatus("b",1);
/// ~~~
/// when calling mytree.GetEntry(i); only branches "a" and "b" will be read.
///
/// __WARNING!!__
/// If your Tree has been created in split mode with a parent branch "parent.",
/// ~~~ {.cpp}
///     mytree.SetBranchStatus("parent",1);
/// ~~~
/// will not activate the sub-branches of "parent". You should do:
/// ~~~ {.cpp}
///     mytree.SetBranchStatus("parent*",1);
/// ~~~
/// Without the trailing dot in the branch creation you have no choice but to
/// call SetBranchStatus explicitly for each of the sub branches.
///
/// An alternative is to call directly
/// ~~~ {.cpp}
///     brancha.GetEntry(i)
///     branchb.GetEntry(i);
/// ~~~
/// ## IMPORTANT NOTE
///
/// By default, GetEntry reuses the space allocated by the previous object
/// for each branch. You can force the previous object to be automatically
/// deleted if you call mybranch.SetAutoDelete(kTRUE) (default is kFALSE).
///
/// Example:
///
/// Consider the example in $ROOTSYS/test/Event.h
/// The top level branch in the tree T is declared with:
/// ~~~ {.cpp}
///     Event *event = 0;  //event must be null or point to a valid object
///                        //it must be initialized
///     T.SetBranchAddress("event",&event);
/// ~~~
/// When reading the Tree, one can choose one of these 3 options:
///
/// ## OPTION 1
///
/// ~~~ {.cpp}
///     for (Long64_t i=0;i<nentries;i++) {
///        T.GetEntry(i);
///        // the object event has been filled at this point
///     }
/// ~~~
/// The default (recommended). At the first entry an object of the class
/// Event will be created and pointed by event. At the following entries,
/// event will be overwritten by the new data. All internal members that are
/// TObject* are automatically deleted. It is important that these members
/// be in a valid state when GetEntry is called. Pointers must be correctly
/// initialized. However these internal members will not be deleted if the
/// characters "->" are specified as the first characters in the comment
/// field of the data member declaration.
///
/// If "->" is specified, the pointer member is read via pointer->Streamer(buf).
/// In this case, it is assumed that the pointer is never null (case of
/// pointer TClonesArray *fTracks in the Event example). If "->" is not
/// specified, the pointer member is read via buf >> pointer. In this case
/// the pointer may be null. Note that the option with "->" is faster to
/// read or write and it also consumes less space in the file.
///
/// ## OPTION 2
///
/// The option AutoDelete is set
/// ~~~ {.cpp}
///     TBranch *branch = T.GetBranch("event");
///     branch->SetAddress(&event);
///     branch->SetAutoDelete(kTRUE);
///     for (Long64_t i=0;i<nentries;i++) {
///        T.GetEntry(i);
///        // the object event has been filled at this point
///     }
/// ~~~
/// In this case, at each iteration, the object event is deleted by GetEntry
/// and a new instance of Event is created and filled.
///
/// ## OPTION 3
///
/// ~~~ {.cpp}
/// Same as option 1, but you delete yourself the event.
///
///     for (Long64_t i=0;i<nentries;i++) {
///        delete event;
///        event = 0;  // EXTREMELY IMPORTANT
///        T.GetEntry(i);
///        // the object event has been filled at this point
///     }
/// ~~~
/// It is strongly recommended to use the default option 1. It has the
/// additional advantage that functions like TTree::Draw (internally calling
/// TTree::GetEntry) will be functional even when the classes in the file are
/// not available.
///
/// Note: See the comments in TBranchElement::SetAddress() for the
/// object ownership policy of the underlying (user) data.

Int_t TTree::GetEntry(Long64_t entry, Int_t getall)
{

   // We already have been visited while recursively looking
   // through the friends tree, let return
   if (kGetEntry & fFriendLockStatus) return 0;

   if (entry < 0 || entry >= fEntries) return 0;
   Int_t i;
   Int_t nbytes = 0;
   fReadEntry = entry;

   // create cache if wanted
   if (fCacheDoAutoInit)
      SetCacheSizeAux();

   Int_t nbranches = fBranches.GetEntriesUnsafe();
   Int_t nb=0;

   auto seqprocessing = [&]() {
      TBranch *branch;
      for (i=0;i<nbranches;i++)  {
         branch = (TBranch*)fBranches.UncheckedAt(i);
         nb = branch->GetEntry(entry, getall);
         if (nb < 0) break;
         nbytes += nb;
      }
   };

#ifdef R__USE_IMT
   if (nbranches > 1 && ROOT::IsImplicitMTEnabled() && fIMTEnabled && !TTreeCacheUnzip::IsParallelUnzip()) {
      if (fSortedBranches.empty())
         InitializeBranchLists(true);

      // Count branches are processed first and sequentially
      for (auto branch : fSeqBranches) {
         nb = branch->GetEntry(entry, getall);
         if (nb < 0) break;
         nbytes += nb;
      }
      if (nb < 0) return nb;

      // Enable this IMT use case (activate its locks)
      ROOT::Internal::TParBranchProcessingRAII pbpRAII;

      Int_t errnb = 0;
      std::atomic<Int_t> pos(0);
      std::atomic<Int_t> nbpar(0);

      auto mapFunction = [&]() {
            // The branch to process is obtained when the task starts to run.
            // This way, since branches are sorted, we make sure that branches
            // leading to big tasks are processed first. If we assigned the
            // branch at task creation time, the scheduler would not necessarily
            // respect our sorting.
            Int_t j = pos.fetch_add(1);

            Int_t nbtask = 0;
            auto branch = fSortedBranches[j].second;

            if (gDebug > 0) {
               std::stringstream ss;
               ss << std::this_thread::get_id();
               Info("GetEntry", "[IMT] Thread %s", ss.str().c_str());
               Info("GetEntry", "[IMT] Running task for branch #%d: %s", j, branch->GetName());
            }

            std::chrono::time_point<std::chrono::system_clock> start, end;

            start = std::chrono::system_clock::now();
            nbtask = branch->GetEntry(entry, getall);
            end = std::chrono::system_clock::now();

            Long64_t tasktime = (Long64_t)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            fSortedBranches[j].first += tasktime;

            if (nbtask < 0) errnb = nbtask;
            else            nbpar += nbtask;
         };

      ROOT::TThreadExecutor pool;
      pool.Foreach(mapFunction, fSortedBranches.size());

      if (errnb < 0) {
         nb = errnb;
      }
      else {
         // Save the number of bytes read by the tasks
         nbytes += nbpar;

         // Re-sort branches if necessary
         if (++fNEntriesSinceSorting == kNEntriesResort) {
            SortBranchesByTime();
            fNEntriesSinceSorting = 0;
         }
      }
   }
   else {
      seqprocessing();
   }
#else
   seqprocessing();
#endif
   if (nb < 0) return nb;

   // GetEntry in list of friends
   if (!fFriends) return nbytes;
   TFriendLock lock(this,kGetEntry);
   TIter nextf(fFriends);
   TFriendElement *fe;
   while ((fe = (TFriendElement*)nextf())) {
      TTree *t = fe->GetTree();
      if (t) {
         if (fe->TestBit(TFriendElement::kFromChain)) {
            nb = t->GetEntry(t->GetReadEntry(),getall);
         } else {
            if ( t->LoadTreeFriend(entry,this) >= 0 ) {
               nb = t->GetEntry(t->GetReadEntry(),getall);
            } else nb = 0;
         }
         if (nb < 0) return nb;
         nbytes += nb;
      }
   }
   return nbytes;
}


////////////////////////////////////////////////////////////////////////////////
/// Divides the top-level branches into two vectors: (i) branches to be
/// processed sequentially and (ii) branches to be processed in parallel.
/// Even if IMT is on, some branches might need to be processed first and in a
/// sequential fashion: in the parallelization of GetEntry, those are the
/// branches that store the size of another branch for every entry
/// (e.g. the size of an array branch). If such branches were processed
/// in parallel with the rest, there could be two threads invoking
/// TBranch::GetEntry on one of them at the same time, since a branch that
/// depends on a size (or count) branch will also invoke GetEntry on the latter.
/// This method can be invoked several times during the event loop if the TTree
/// is being written, for example when adding new branches. In these cases, the
/// `checkLeafCount` parameter is false.
/// \param[in] checkLeafCount True if we need to check whether some branches are
///                           count leaves.

void TTree::InitializeBranchLists(bool checkLeafCount)
{
   Int_t nbranches = fBranches.GetEntriesFast();

   // The special branch fBranchRef needs to be processed sequentially:
   // we add it once only.
   if (fBranchRef && fBranchRef != fSeqBranches[0]) {
      fSeqBranches.push_back(fBranchRef);
   }

   // The branches to be processed sequentially are those that are the leaf count of another branch
   if (checkLeafCount) {
      for (Int_t i = 0; i < nbranches; i++)  {
         TBranch* branch = (TBranch*)fBranches.UncheckedAt(i);
         auto leafCount = ((TLeaf*)branch->GetListOfLeaves()->At(0))->GetLeafCount();
         if (leafCount) {
            auto countBranch = leafCount->GetBranch();
            if (std::find(fSeqBranches.begin(), fSeqBranches.end(), countBranch) == fSeqBranches.end()) {
               fSeqBranches.push_back(countBranch);
            }
         }
      }
   }

   // Any branch that is not a leaf count can be safely processed in parallel when reading
   // We need to reset the vector to make sure we do not re-add several times the same branch.
   if (!checkLeafCount) {
      fSortedBranches.clear();
   }
   for (Int_t i = 0; i < nbranches; i++)  {
      Long64_t bbytes = 0;
      TBranch* branch = (TBranch*)fBranches.UncheckedAt(i);
      if (std::find(fSeqBranches.begin(), fSeqBranches.end(), branch) == fSeqBranches.end()) {
         bbytes = branch->GetTotBytes("*");
         fSortedBranches.emplace_back(bbytes, branch);
      }
   }

   // Initially sort parallel branches by size
   std::sort(fSortedBranches.begin(),
             fSortedBranches.end(),
             [](std::pair<Long64_t,TBranch*> a, std::pair<Long64_t,TBranch*> b) {
                return a.first > b.first;
             });

   for (size_t i = 0; i < fSortedBranches.size(); i++)  {
      fSortedBranches[i].first = 0LL;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sorts top-level branches by the last average task time recorded per branch.

void TTree::SortBranchesByTime()
{
   for (size_t i = 0; i < fSortedBranches.size(); i++)  {
      fSortedBranches[i].first *= kNEntriesResortInv;
   }

   std::sort(fSortedBranches.begin(),
             fSortedBranches.end(),
             [](std::pair<Long64_t,TBranch*> a, std::pair<Long64_t,TBranch*> b) {
                return a.first > b.first;
             });

   for (size_t i = 0; i < fSortedBranches.size(); i++)  {
      fSortedBranches[i].first = 0LL;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Returns the entry list assigned to this tree

TEntryList* TTree::GetEntryList()
{
   return fEntryList;
}

////////////////////////////////////////////////////////////////////////////////
/// Return entry number corresponding to entry.
///
/// if no TEntryList set returns entry
/// else returns the entry number corresponding to the list index=entry

Long64_t TTree::GetEntryNumber(Long64_t entry) const
{
   if (!fEntryList) {
      return entry;
   }

   return fEntryList->GetEntry(entry);
}

////////////////////////////////////////////////////////////////////////////////
/// Return entry number corresponding to major and minor number.
/// Note that this function returns only the entry number, not the data
/// To read the data corresponding to an entry number, use TTree::GetEntryWithIndex
/// the BuildIndex function has created a table of Long64_t* of sorted values
/// corresponding to val = major<<31 + minor;
/// The function performs binary search in this sorted table.
/// If it finds a pair that matches val, it returns directly the
/// index in the table.
/// If an entry corresponding to major and minor is not found, the function
/// returns the index of the major,minor pair immediately lower than the
/// requested value, ie it will return -1 if the pair is lower than
/// the first entry in the index.
///
/// See also GetEntryNumberWithIndex

Long64_t TTree::GetEntryNumberWithBestIndex(Long64_t major, Long64_t minor) const
{
   if (!fTreeIndex) {
      return -1;
   }
   return fTreeIndex->GetEntryNumberWithBestIndex(major, minor);
}

////////////////////////////////////////////////////////////////////////////////
/// Return entry number corresponding to major and minor number.
/// Note that this function returns only the entry number, not the data
/// To read the data corresponding to an entry number, use TTree::GetEntryWithIndex
/// the BuildIndex function has created a table of Long64_t* of sorted values
/// corresponding to val = major<<31 + minor;
/// The function performs binary search in this sorted table.
/// If it finds a pair that matches val, it returns directly the
/// index in the table, otherwise it returns -1.
///
/// See also GetEntryNumberWithBestIndex

Long64_t TTree::GetEntryNumberWithIndex(Long64_t major, Long64_t minor) const
{
   if (!fTreeIndex) {
      return -1;
   }
   return fTreeIndex->GetEntryNumberWithIndex(major, minor);
}

////////////////////////////////////////////////////////////////////////////////
/// Read entry corresponding to major and minor number.
///
///  The function returns the total number of bytes read.
///  If the Tree has friend trees, the corresponding entry with
///  the index values (major,minor) is read. Note that the master Tree
///  and its friend may have different entry serial numbers corresponding
///  to (major,minor).

Int_t TTree::GetEntryWithIndex(Int_t major, Int_t minor)
{
   // We already have been visited while recursively looking
   // through the friends tree, let's return.
   if (kGetEntryWithIndex & fFriendLockStatus) {
      return 0;
   }
   Long64_t serial = GetEntryNumberWithIndex(major, minor);
   if (serial < 0) {
      return -1;
   }
   // create cache if wanted
   if (fCacheDoAutoInit)
      SetCacheSizeAux();

   Int_t i;
   Int_t nbytes = 0;
   fReadEntry = serial;
   TBranch *branch;
   Int_t nbranches = fBranches.GetEntriesFast();
   Int_t nb;
   for (i = 0; i < nbranches; ++i) {
      branch = (TBranch*)fBranches.UncheckedAt(i);
      nb = branch->GetEntry(serial);
      if (nb < 0) return nb;
      nbytes += nb;
   }
   // GetEntry in list of friends
   if (!fFriends) return nbytes;
   TFriendLock lock(this,kGetEntryWithIndex);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree *t = fe->GetTree();
      if (t) {
         serial = t->GetEntryNumberWithIndex(major,minor);
         if (serial <0) return -nbytes;
         nb = t->GetEntry(serial);
         if (nb < 0) return nb;
         nbytes += nb;
      }
   }
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TTree friend whose name or alias is `friendname`.

TTree* TTree::GetFriend(const char *friendname) const
{

   // We already have been visited while recursively
   // looking through the friends tree, let's return.
   if (kGetFriend & fFriendLockStatus) {
      return 0;
   }
   if (!fFriends) {
      return 0;
   }
   TFriendLock lock(const_cast<TTree*>(this), kGetFriend);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      if (strcmp(friendname,fe->GetName())==0
          || strcmp(friendname,fe->GetTreeName())==0) {
         return fe->GetTree();
      }
   }
   // After looking at the first level,
   // let's see if it is a friend of friends.
   nextf.Reset();
   fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree *res = fe->GetTree()->GetFriend(friendname);
      if (res) {
         return res;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// If the 'tree' is a friend, this method returns its alias name.
///
/// This alias is an alternate name for the tree.
///
/// It can be used in conjunction with a branch or leaf name in a TTreeFormula,
/// to specify in which particular tree the branch or leaf can be found if
/// the friend trees have branches or leaves with the same name as the master
/// tree.
///
/// It can also be used in conjunction with an alias created using
/// TTree::SetAlias in a TTreeFormula, e.g.:
/// ~~~ {.cpp}
///      maintree->Draw("treealias.fPx - treealias.myAlias");
/// ~~~
/// where fPx is a branch of the friend tree aliased as 'treealias' and 'myAlias'
/// was created using TTree::SetAlias on the friend tree.
///
/// However, note that 'treealias.myAlias' will be expanded literally,
/// without remembering that it comes from the aliased friend and thus
/// the branch name might not be disambiguated properly, which means
/// that you may not be able to take advantage of this feature.
///

const char* TTree::GetFriendAlias(TTree* tree) const
{
   if ((tree == this) || (tree == GetTree())) {
      return 0;
   }

   // We already have been visited while recursively
   // looking through the friends tree, let's return.
   if (kGetFriendAlias & fFriendLockStatus) {
      return 0;
   }
   if (!fFriends) {
      return 0;
   }
   TFriendLock lock(const_cast<TTree*>(this), kGetFriendAlias);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree* t = fe->GetTree();
      if (t == tree) {
         return fe->GetName();
      }
      // Case of a chain:
      if (t && t->GetTree() == tree) {
         return fe->GetName();
      }
   }
   // After looking at the first level,
   // let's see if it is a friend of friends.
   nextf.Reset();
   fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      const char* res = fe->GetTree()->GetFriendAlias(tree);
      if (res) {
         return res;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current set of IO settings
ROOT::TIOFeatures TTree::GetIOFeatures() const
{
   return fIOFeatures;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a new iterator that will go through all the leaves on the tree itself and its friend.

TIterator* TTree::GetIteratorOnAllLeaves(Bool_t dir)
{
   return new TTreeFriendLeafIter(this, dir);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the 1st Leaf named name in any Branch of this
/// Tree or any branch in the list of friend trees.
///
/// The leaf name can contain the name of a friend tree with the
/// syntax: friend_dir_and_tree.full_leaf_name
/// the friend_dir_and_tree can be of the form:
/// ~~~ {.cpp}
///     TDirectoryName/TreeName
/// ~~~

TLeaf* TTree::GetLeafImpl(const char* branchname, const char *leafname)
{
   TLeaf *leaf = 0;
   if (branchname) {
      TBranch *branch = FindBranch(branchname);
      if (branch) {
         leaf = branch->GetLeaf(leafname);
         if (leaf) {
            return leaf;
         }
      }
   }
   TIter nextl(GetListOfLeaves());
   while ((leaf = (TLeaf*)nextl())) {
      if (strcmp(leaf->GetFullName(), leafname) != 0 && strcmp(leaf->GetName(), leafname) != 0)
         continue; // leafname does not match GetName() nor GetFullName(), this is not the right leaf
      if (branchname) {
         // check the branchname is also a match
         TBranch *br = leaf->GetBranch();
         // if a quick comparison with the branch full name is a match, we are done
         if (!strcmp(br->GetFullName(), branchname))
            return leaf;
         UInt_t nbch = strlen(branchname);
         const char* brname = br->GetName();
         TBranch *mother = br->GetMother();
         if (strncmp(brname,branchname,nbch)) {
            if (mother != br) {
               const char *mothername = mother->GetName();
               UInt_t motherlen = strlen(mothername);
               if (!strcmp(mothername, branchname)) {
                  return leaf;
               } else if (nbch > motherlen && strncmp(mothername,branchname,motherlen)==0 && (mothername[motherlen-1]=='.' || branchname[motherlen]=='.')) {
                  // The left part of the requested name match the name of the mother, let's see if the right part match the name of the branch.
                  if (strncmp(brname,branchname+motherlen+1,nbch-motherlen-1)) {
                     // No it does not
                     continue;
                  } // else we have match so we can proceed.
               } else {
                  // no match
                  continue;
               }
            } else {
               continue;
            }
         }
         // The start of the branch name is identical to the content
         // of 'aname' before the first '/'.
         // Let's make sure that it is not longer (we are trying
         // to avoid having jet2/value match the branch jet23
         if ((strlen(brname) > nbch) && (brname[nbch] != '.') && (brname[nbch] != '[')) {
            continue;
         }
      }
      return leaf;
   }
   if (!fFriends) return 0;
   TFriendLock lock(this,kGetLeaf);
   TIter next(fFriends);
   TFriendElement *fe;
   while ((fe = (TFriendElement*)next())) {
      TTree *t = fe->GetTree();
      if (t) {
         leaf = t->GetLeaf(branchname, leafname);
         if (leaf) return leaf;
      }
   }

   //second pass in the list of friends when the leaf name
   //is prefixed by the tree name
   TString strippedArg;
   next.Reset();
   while ((fe = (TFriendElement*)next())) {
      TTree *t = fe->GetTree();
      if (t==0) continue;
      char *subname = (char*)strstr(leafname,fe->GetName());
      if (subname != leafname) continue;
      Int_t l = strlen(fe->GetName());
      subname += l;
      if (*subname != '.') continue;
      subname++;
      strippedArg += subname;
      leaf = t->GetLeaf(branchname,subname);
      if (leaf) return leaf;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the 1st Leaf named name in any Branch of this
/// Tree or any branch in the list of friend trees.
///
/// The leaf name can contain the name of a friend tree with the
/// syntax: friend_dir_and_tree.full_leaf_name
/// the friend_dir_and_tree can be of the form:
///
///     TDirectoryName/TreeName

TLeaf* TTree::GetLeaf(const char* branchname, const char *leafname)
{
   if (leafname == 0) return 0;

   // We already have been visited while recursively looking
   // through the friends tree, let return
   if (kGetLeaf & fFriendLockStatus) {
      return 0;
   }

   return GetLeafImpl(branchname,leafname);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to first leaf named \param[name] in any branch of this
/// tree or its friend trees.
///
/// \param[name] may be in the form 'branch/leaf'
///

TLeaf* TTree::GetLeaf(const char *name)
{
   // Return nullptr if name is invalid or if we have
   // already been visited while searching friend trees
   if (!name || (kGetLeaf & fFriendLockStatus))
      return nullptr;

   std::string path(name);
   const auto sep = path.find_last_of("/");
   if (sep != std::string::npos)
      return GetLeafImpl(path.substr(0, sep).c_str(), name+sep+1);

   return GetLeafImpl(nullptr, name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum of column with name columname.
/// if the Tree has an associated TEventList or TEntryList, the maximum
/// is computed for the entries in this list.

Double_t TTree::GetMaximum(const char* columname)
{
   TLeaf* leaf = this->GetLeaf(columname);
   if (!leaf) {
      return 0;
   }

   // create cache if wanted
   if (fCacheDoAutoInit)
      SetCacheSizeAux();

   TBranch* branch = leaf->GetBranch();
   Double_t cmax = -DBL_MAX;
   for (Long64_t i = 0; i < fEntries; ++i) {
      Long64_t entryNumber = this->GetEntryNumber(i);
      if (entryNumber < 0) break;
      branch->GetEntry(entryNumber);
      for (Int_t j = 0; j < leaf->GetLen(); ++j) {
         Double_t val = leaf->GetValue(j);
         if (val > cmax) {
            cmax = val;
         }
      }
   }
   return cmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function which returns the tree file size limit in bytes.

Long64_t TTree::GetMaxTreeSize()
{
   return fgMaxTreeSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return minimum of column with name columname.
/// if the Tree has an associated TEventList or TEntryList, the minimum
/// is computed for the entries in this list.

Double_t TTree::GetMinimum(const char* columname)
{
   TLeaf* leaf = this->GetLeaf(columname);
   if (!leaf) {
      return 0;
   }

   // create cache if wanted
   if (fCacheDoAutoInit)
      SetCacheSizeAux();

   TBranch* branch = leaf->GetBranch();
   Double_t cmin = DBL_MAX;
   for (Long64_t i = 0; i < fEntries; ++i) {
      Long64_t entryNumber = this->GetEntryNumber(i);
      if (entryNumber < 0) break;
      branch->GetEntry(entryNumber);
      for (Int_t j = 0;j < leaf->GetLen(); ++j) {
         Double_t val = leaf->GetValue(j);
         if (val < cmin) {
            cmin = val;
         }
      }
   }
   return cmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the TTreePlayer (if not already done).

TVirtualTreePlayer* TTree::GetPlayer()
{
   if (fPlayer) {
      return fPlayer;
   }
   fPlayer = TVirtualTreePlayer::TreePlayer(this);
   return fPlayer;
}

////////////////////////////////////////////////////////////////////////////////
/// Find and return the TTreeCache registered with the file and which may
/// contain branches for us.

TTreeCache *TTree::GetReadCache(TFile *file) const
{
   TTreeCache *pe = dynamic_cast<TTreeCache*>(file->GetCacheRead(GetTree()));
   if (pe && pe->GetTree() != GetTree())
      pe = nullptr;
   return pe;
}

////////////////////////////////////////////////////////////////////////////////
/// Find and return the TTreeCache registered with the file and which may
/// contain branches for us. If create is true and there is no cache
/// a new cache is created with default size.

TTreeCache *TTree::GetReadCache(TFile *file, Bool_t create)
{
   TTreeCache *pe = GetReadCache(file);
   if (create && !pe) {
      if (fCacheDoAutoInit)
         SetCacheSizeAux(kTRUE, -1);
      pe = dynamic_cast<TTreeCache*>(file->GetCacheRead(GetTree()));
      if (pe && pe->GetTree() != GetTree()) pe = 0;
   }
   return pe;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the list containing user objects associated to this tree.
///
/// The list is automatically created if it does not exist.
///
/// WARNING: By default the TTree destructor will delete all objects added
/// to this list. If you do not want these objects to be deleted,
/// call:
///
///     mytree->GetUserInfo()->Clear();
///
/// before deleting the tree.

TList* TTree::GetUserInfo()
{
   if (!fUserInfo) {
      fUserInfo = new TList();
      fUserInfo->SetName("UserInfo");
   }
   return fUserInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Appends the cluster range information stored in 'fromtree' to this tree,
/// including the value of fAutoFlush.
///
/// This is used when doing a fast cloning (by TTreeCloner).
/// See also fAutoFlush and fAutoSave if needed.

void TTree::ImportClusterRanges(TTree *fromtree)
{
   Long64_t autoflush = fromtree->GetAutoFlush();
   if (fromtree->fNClusterRange == 0 && fromtree->fAutoFlush == fAutoFlush) {
      // nothing to do
   } else if (fNClusterRange || fromtree->fNClusterRange) {
      Int_t newsize = fNClusterRange + 1 + fromtree->fNClusterRange;
      if (newsize > fMaxClusterRange) {
         if (fMaxClusterRange) {
            fClusterRangeEnd = (Long64_t*)TStorage::ReAlloc(fClusterRangeEnd,
                                                            newsize*sizeof(Long64_t),fMaxClusterRange*sizeof(Long64_t));
            fClusterSize = (Long64_t*)TStorage::ReAlloc(fClusterSize,
                                                        newsize*sizeof(Long64_t),fMaxClusterRange*sizeof(Long64_t));
            fMaxClusterRange = newsize;
         } else {
            fMaxClusterRange = newsize;
            fClusterRangeEnd = new Long64_t[fMaxClusterRange];
            fClusterSize = new Long64_t[fMaxClusterRange];
         }
      }
      if (fEntries) {
         fClusterRangeEnd[fNClusterRange] = fEntries - 1;
         fClusterSize[fNClusterRange] = fAutoFlush<0 ? 0 : fAutoFlush;
         ++fNClusterRange;
      }
      for (Int_t i = 0 ; i < fromtree->fNClusterRange; ++i) {
         fClusterRangeEnd[fNClusterRange] = fEntries + fromtree->fClusterRangeEnd[i];
         fClusterSize[fNClusterRange] = fromtree->fClusterSize[i];
         ++fNClusterRange;
      }
      fAutoFlush = autoflush;
   } else {
      SetAutoFlush( autoflush );
   }
   Long64_t autosave = GetAutoSave();
   if (autoflush > 0 && autosave > 0) {
      SetAutoSave( autoflush*(autosave/autoflush) );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Keep a maximum of fMaxEntries in memory.

void TTree::KeepCircular()
{
   Int_t nb = fBranches.GetEntriesFast();
   Long64_t maxEntries = fMaxEntries - (fMaxEntries / 10);
   for (Int_t i = 0; i < nb; ++i)  {
      TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
      branch->KeepCircular(maxEntries);
   }
   if (fNClusterRange) {
      Long64_t entriesOffset = fEntries - maxEntries;
      Int_t oldsize = fNClusterRange;
      for(Int_t i = 0, j = 0; j < oldsize; ++j) {
         if (fClusterRangeEnd[j] > entriesOffset) {
            fClusterRangeEnd[i] =  fClusterRangeEnd[j] - entriesOffset;
            ++i;
         } else {
            --fNClusterRange;
         }
      }
   }
   fEntries = maxEntries;
   fReadEntry = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Read in memory all baskets from all branches up to the limit of maxmemory bytes.
///
/// If maxmemory is non null and positive SetMaxVirtualSize is called
/// with this value. Default for maxmemory is 2000000000 (2 Gigabytes).
/// The function returns the total number of baskets read into memory
/// if negative an error occurred while loading the branches.
/// This method may be called to force branch baskets in memory
/// when random access to branch entries is required.
/// If random access to only a few branches is required, you should
/// call directly TBranch::LoadBaskets.

Int_t TTree::LoadBaskets(Long64_t maxmemory)
{
   if (maxmemory > 0) SetMaxVirtualSize(maxmemory);

   TIter next(GetListOfLeaves());
   TLeaf *leaf;
   Int_t nimported = 0;
   while ((leaf=(TLeaf*)next())) {
      nimported += leaf->GetBranch()->LoadBaskets();//break;
   }
   return nimported;
}

////////////////////////////////////////////////////////////////////////////////
/// Set current entry.
///
/// Returns -2 if entry does not exist (just as TChain::LoadTree()).
/// Returns -6 if an error occurs in the notification callback (just as TChain::LoadTree()).
///
/// Note: This function is overloaded in TChain.
///

Long64_t TTree::LoadTree(Long64_t entry)
{
   // We already have been visited while recursively looking
   // through the friends tree, let return
   if (kLoadTree & fFriendLockStatus) {
      // We need to return a negative value to avoid a circular list of friend
      // to think that there is always an entry somewhere in the list.
      return -1;
   }

   // create cache if wanted
   if (fCacheDoAutoInit && entry >=0)
      SetCacheSizeAux();

   if (fNotify) {
      if (fReadEntry < 0) {
         fNotify->Notify();
      }
   }
   fReadEntry = entry;

   Bool_t friendHasEntry = kFALSE;
   if (fFriends) {
      // Set current entry in friends as well.
      //
      // An alternative would move this code to each of the
      // functions calling LoadTree (and to overload a few more).
      Bool_t needUpdate = kFALSE;
      {
         // This scope is need to insure the lock is released at the right time
         TIter nextf(fFriends);
         TFriendLock lock(this, kLoadTree);
         TFriendElement* fe = 0;
         while ((fe = (TFriendElement*) nextf())) {
            if (fe->TestBit(TFriendElement::kFromChain)) {
               // This friend element was added by the chain that owns this
               // tree, the chain will deal with loading the correct entry.
               continue;
            }
            TTree* friendTree = fe->GetTree();
            if (friendTree) {
               if (friendTree->LoadTreeFriend(entry, this) >= 0) {
                  friendHasEntry = kTRUE;
               }
            }
            if (fe->IsUpdated()) {
               needUpdate = kTRUE;
               fe->ResetUpdated();
            }
         } // for each friend
      }
      if (needUpdate) {
         //update list of leaves in all TTreeFormula of the TTreePlayer (if any)
         if (fPlayer) {
            fPlayer->UpdateFormulaLeaves();
         }
         //Notify user if requested
         if (fNotify) {
            if(!fNotify->Notify()) return -6;
         }
      }
   }

   if ((fReadEntry >= fEntries) && !friendHasEntry) {
      fReadEntry = -1;
      return -2;
   }
   return fReadEntry;
}

////////////////////////////////////////////////////////////////////////////////
/// Load entry on behalf of our master tree, we may use an index.
///
/// Called by LoadTree() when the masterTree looks for the entry
/// number in a friend tree (us) corresponding to the passed entry
/// number in the masterTree.
///
/// If we have no index, our entry number and the masterTree entry
/// number are the same.
///
/// If we *do* have an index, we must find the (major, minor) value pair
/// in masterTree to locate our corresponding entry.
///

Long64_t TTree::LoadTreeFriend(Long64_t entry, TTree* masterTree)
{
   if (!fTreeIndex) {
      return LoadTree(entry);
   }
   return LoadTree(fTreeIndex->GetEntryNumberFriend(masterTree));
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a skeleton analysis class for this tree.
///
/// The following files are produced: classname.h and classname.C.
/// If classname is 0, classname will be called "nameoftree".
///
/// The generated code in classname.h includes the following:
///
/// - Identification of the original tree and the input file name.
/// - Definition of an analysis class (data members and member functions).
/// - The following member functions:
///   - constructor (by default opening the tree file),
///   - GetEntry(Long64_t entry),
///   - Init(TTree* tree) to initialize a new TTree,
///   - Show(Long64_t entry) to read and dump entry.
///
/// The generated code in classname.C includes only the main
/// analysis function Loop.
///
/// To use this function:
///
/// - Open your tree file (eg: TFile f("myfile.root");)
/// - T->MakeClass("MyClass");
///
/// where T is the name of the TTree in file myfile.root,
/// and MyClass.h, MyClass.C the name of the files created by this function.
/// In a ROOT session, you can do:
/// ~~~ {.cpp}
///     root > .L MyClass.C
///     root > MyClass* t = new MyClass;
///     root > t->GetEntry(12); // Fill data members of t with entry number 12.
///     root > t->Show();       // Show values of entry 12.
///     root > t->Show(16);     // Read and show values of entry 16.
///     root > t->Loop();       // Loop on all entries.
/// ~~~
/// NOTE: Do not use the code generated for a single TTree which is part
/// of a TChain to process that entire TChain.  The maximum dimensions
/// calculated for arrays on the basis of a single TTree from the TChain
/// might be (will be!) too small when processing all of the TTrees in
/// the TChain.  You must use myChain.MakeClass() to generate the code,
/// not myTree.MakeClass(...).

Int_t TTree::MakeClass(const char* classname, Option_t* option)
{
   GetPlayer();
   if (!fPlayer) {
      return 0;
   }
   return fPlayer->MakeClass(classname, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a skeleton function for this tree.
///
/// The function code is written on filename.
/// If filename is 0, filename will be called nameoftree.C
///
/// The generated code includes the following:
/// - Identification of the original Tree and Input file name,
/// - Opening the Tree file,
/// - Declaration of Tree variables,
/// - Setting of branches addresses,
/// - A skeleton for the entry loop.
///
/// To use this function:
///
/// - Open your Tree file (eg: TFile f("myfile.root");)
/// - T->MakeCode("MyAnalysis.C");
///
/// where T is the name of the TTree in file myfile.root
/// and MyAnalysis.C the name of the file created by this function.
///
/// NOTE: Since the implementation of this function, a new and better
/// function TTree::MakeClass() has been developed.

Int_t TTree::MakeCode(const char* filename)
{
   Warning("MakeCode", "MakeCode is obsolete. Use MakeClass or MakeSelector instead");

   GetPlayer();
   if (!fPlayer) return 0;
   return fPlayer->MakeCode(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a skeleton analysis class for this Tree using TBranchProxy.
///
/// TBranchProxy is the base of a class hierarchy implementing an
/// indirect access to the content of the branches of a TTree.
///
/// "proxyClassname" is expected to be of the form:
/// ~~~ {.cpp}
///     [path/]fileprefix
/// ~~~
/// The skeleton will then be generated in the file:
/// ~~~ {.cpp}
///     fileprefix.h
/// ~~~
/// located in the current directory or in 'path/' if it is specified.
/// The class generated will be named 'fileprefix'
///
/// "macrofilename" and optionally "cutfilename" are expected to point
/// to source files which will be included by the generated skeleton.
/// Method of the same name as the file(minus the extension and path)
/// will be called by the generated skeleton's Process method as follow:
/// ~~~ {.cpp}
///     [if (cutfilename())] htemp->Fill(macrofilename());
/// ~~~
/// "option" can be used select some of the optional features during
/// the code generation.  The possible options are:
///
/// - nohist : indicates that the generated ProcessFill should not fill the histogram.
///
/// 'maxUnrolling' controls how deep in the class hierarchy does the
/// system 'unroll' classes that are not split.  Unrolling a class
/// allows direct access to its data members (this emulates the behavior
/// of TTreeFormula).
///
/// The main features of this skeleton are:
///
/// * on-demand loading of branches
/// * ability to use the 'branchname' as if it was a data member
/// * protection against array out-of-bounds errors
/// * ability to use the branch data as an object (when the user code is available)
///
/// For example with Event.root, if
/// ~~~ {.cpp}
///     Double_t somePx = fTracks.fPx[2];
/// ~~~
/// is executed by one of the method of the skeleton,
/// somePx will updated with the current value of fPx of the 3rd track.
///
/// Both macrofilename and the optional cutfilename are expected to be
/// the name of source files which contain at least a free standing
/// function with the signature:
/// ~~~ {.cpp}
///     x_t macrofilename(); // i.e function with the same name as the file
/// ~~~
/// and
/// ~~~ {.cpp}
///     y_t cutfilename();   // i.e function with the same name as the file
/// ~~~
/// x_t and y_t needs to be types that can convert respectively to a double
/// and a bool (because the skeleton uses:
///
///     if (cutfilename()) htemp->Fill(macrofilename());
///
/// These two functions are run in a context such that the branch names are
/// available as local variables of the correct (read-only) type.
///
/// Note that if you use the same 'variable' twice, it is more efficient
/// to 'cache' the value. For example:
/// ~~~ {.cpp}
///     Int_t n = fEventNumber; // Read fEventNumber
///     if (n<10 || n>10) { ... }
/// ~~~
/// is more efficient than
/// ~~~ {.cpp}
///     if (fEventNumber<10 || fEventNumber>10)
/// ~~~
/// Also, optionally, the generated selector will also call methods named
/// macrofilename_methodname in each of 6 main selector methods if the method
/// macrofilename_methodname exist (Where macrofilename is stripped of its
/// extension).
///
/// Concretely, with the script named h1analysisProxy.C,
///
/// - The method         calls the method (if it exist)
/// - Begin           -> void h1analysisProxy_Begin(TTree*);
/// - SlaveBegin      -> void h1analysisProxy_SlaveBegin(TTree*);
/// - Notify          -> Bool_t h1analysisProxy_Notify();
/// - Process         -> Bool_t h1analysisProxy_Process(Long64_t);
/// - SlaveTerminate  -> void h1analysisProxy_SlaveTerminate();
/// - Terminate       -> void h1analysisProxy_Terminate();
///
/// If a file name macrofilename.h (or .hh, .hpp, .hxx, .hPP, .hXX) exist
/// it is included before the declaration of the proxy class.  This can
/// be used in particular to insure that the include files needed by
/// the macro file are properly loaded.
///
/// The default histogram is accessible via the variable named 'htemp'.
///
/// If the library of the classes describing the data in the branch is
/// loaded, the skeleton will add the needed `include` statements and
/// give the ability to access the object stored in the branches.
///
/// To draw px using the file hsimple.root (generated by the
/// hsimple.C tutorial), we need a file named hsimple.cxx:
/// ~~~ {.cpp}
///     double hsimple() {
///        return px;
///     }
/// ~~~
/// MakeProxy can then be used indirectly via the TTree::Draw interface
/// as follow:
/// ~~~ {.cpp}
///     new TFile("hsimple.root")
///     ntuple->Draw("hsimple.cxx");
/// ~~~
/// A more complete example is available in the tutorials directory:
/// h1analysisProxy.cxx , h1analysProxy.h and h1analysisProxyCut.C
/// which reimplement the selector found in h1analysis.C

Int_t TTree::MakeProxy(const char* proxyClassname, const char* macrofilename, const char* cutfilename, const char* option, Int_t maxUnrolling)
{
   GetPlayer();
   if (!fPlayer) return 0;
   return fPlayer->MakeProxy(proxyClassname,macrofilename,cutfilename,option,maxUnrolling);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate skeleton selector class for this tree.
///
/// The following files are produced: selector.h and selector.C.
/// If selector is 0, the selector will be called "nameoftree".
/// The option can be used to specify the branches that will have a data member.
///    - If option is "=legacy", a pre-ROOT6 selector will be generated (data
///      members and branch pointers instead of TTreeReaders).
///    - If option is empty, readers will be generated for each leaf.
///    - If option is "@", readers will be generated for the topmost branches.
///    - Individual branches can also be picked by their name:
///       - "X" generates readers for leaves of X.
///       - "@X" generates a reader for X as a whole.
///       - "@X;Y" generates a reader for X as a whole and also readers for the
///         leaves of Y.
///    - For further examples see the figure below.
///
/// \image html ttree_makeselector_option_examples.png
///
/// The generated code in selector.h includes the following:
///    - Identification of the original Tree and Input file name
///    - Definition of selector class (data and functions)
///    - The following class functions:
///       - constructor and destructor
///       - void    Begin(TTree *tree)
///       - void    SlaveBegin(TTree *tree)
///       - void    Init(TTree *tree)
///       - Bool_t  Notify()
///       - Bool_t  Process(Long64_t entry)
///       - void    Terminate()
///       - void    SlaveTerminate()
///
/// The class selector derives from TSelector.
/// The generated code in selector.C includes empty functions defined above.
///
/// To use this function:
///
///    - connect your Tree file (eg: `TFile f("myfile.root");`)
///    - `T->MakeSelector("myselect");`
///
/// where T is the name of the Tree in file myfile.root
/// and myselect.h, myselect.C the name of the files created by this function.
/// In a ROOT session, you can do:
/// ~~~ {.cpp}
///     root > T->Process("myselect.C")
/// ~~~

Int_t TTree::MakeSelector(const char* selector, Option_t* option)
{
   TString opt(option);
   if(opt.EqualTo("=legacy", TString::ECaseCompare::kIgnoreCase)) {
      return MakeClass(selector, "selector");
   } else {
      GetPlayer();
      if (!fPlayer) return 0;
      return fPlayer->MakeReader(selector, option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if adding nbytes to memory we are still below MaxVirtualsize.

Bool_t TTree::MemoryFull(Int_t nbytes)
{
   if ((fTotalBuffers + nbytes) < fMaxVirtualSize) {
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function merging the trees in the TList into a new tree.
///
/// Trees in the list can be memory or disk-resident trees.
/// The new tree is created in the current directory (memory if gROOT).

TTree* TTree::MergeTrees(TList* li, Option_t* options)
{
   if (!li) return 0;
   TIter next(li);
   TTree *newtree = 0;
   TObject *obj;

   while ((obj=next())) {
      if (!obj->InheritsFrom(TTree::Class())) continue;
      TTree *tree = (TTree*)obj;
      Long64_t nentries = tree->GetEntries();
      if (nentries == 0) continue;
      if (!newtree) {
         newtree = (TTree*)tree->CloneTree();
         if (!newtree) continue;

         // Once the cloning is done, separate the trees,
         // to avoid as many side-effects as possible
         // The list of clones is guaranteed to exist since we
         // just cloned the tree.
         tree->GetListOfClones()->Remove(newtree);
         tree->ResetBranchAddresses();
         newtree->ResetBranchAddresses();
         continue;
      }

      newtree->CopyAddresses(tree);

      newtree->CopyEntries(tree,-1,options);

      tree->ResetBranchAddresses(); // Disconnect from new tree.
   }
   if (newtree && newtree->GetTreeIndex()) {
      newtree->GetTreeIndex()->Append(0,kFALSE); // Force the sorting
   }
   return newtree;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge the trees in the TList into this tree.
///
/// Returns the total number of entries in the merged tree.

Long64_t TTree::Merge(TCollection* li, Option_t *options)
{
   if (!li) return 0;
   Long64_t storeAutoSave = fAutoSave;
   // Disable the autosave as the TFileMerge keeps a list of key and deleting the underlying
   // key would invalidate its iteration (or require costly measure to not use the deleted keys).
   // Also since this is part of a merging operation, the output file is not as precious as in
   // the general case since the input file should still be around.
   fAutoSave = 0;
   TIter next(li);
   TTree *tree;
   while ((tree = (TTree*)next())) {
      if (tree==this) continue;
      if (!tree->InheritsFrom(TTree::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s", tree->ClassName(), ClassName());
         fAutoSave = storeAutoSave;
         return -1;
      }

      Long64_t nentries = tree->GetEntries();
      if (nentries == 0) continue;

      CopyAddresses(tree);

      CopyEntries(tree,-1,options);

      tree->ResetBranchAddresses();
   }
   fAutoSave = storeAutoSave;
   return GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Merge the trees in the TList into this tree.
/// If info->fIsFirst is true, first we clone this TTree info the directory
/// info->fOutputDirectory and then overlay the new TTree information onto
/// this TTree object (so that this TTree object is now the appropriate to
/// use for further merging).
///
/// Returns the total number of entries in the merged tree.

Long64_t TTree::Merge(TCollection* li, TFileMergeInfo *info)
{
   const char *options = info ? info->fOptions.Data() : "";
   if (info && info->fIsFirst && info->fOutputDirectory && info->fOutputDirectory->GetFile() != GetCurrentFile()) {
      TDirectory::TContext ctxt(info->fOutputDirectory);
      TIOFeatures saved_features = fIOFeatures;
      if (info->fIOFeatures) {
         fIOFeatures = *(info->fIOFeatures);
      }
      TTree *newtree = CloneTree(-1, options);
      fIOFeatures = saved_features;
      if (newtree) {
         newtree->Write();
         delete newtree;
      }
      // Make sure things are really written out to disk before attempting any reading.
      info->fOutputDirectory->GetFile()->Flush();
      info->fOutputDirectory->ReadTObject(this,this->GetName());
   }
   if (!li) return 0;
   Long64_t storeAutoSave = fAutoSave;
   // Disable the autosave as the TFileMerge keeps a list of key and deleting the underlying
   // key would invalidate its iteration (or require costly measure to not use the deleted keys).
   // Also since this is part of a merging operation, the output file is not as precious as in
   // the general case since the input file should still be around.
   fAutoSave = 0;
   TIter next(li);
   TTree *tree;
   while ((tree = (TTree*)next())) {
      if (tree==this) continue;
      if (!tree->InheritsFrom(TTree::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s", tree->ClassName(), ClassName());
         fAutoSave = storeAutoSave;
         return -1;
      }
      // Copy MakeClass status.
      tree->SetMakeClass(fMakeClass);

      // Copy branch addresses.
      CopyAddresses(tree);

      CopyEntries(tree,-1,options);

      tree->ResetBranchAddresses();
   }
   fAutoSave = storeAutoSave;
   return GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Move a cache from a file to the current file in dir.
/// if src is null no operation is done, if dir is null or there is no
/// current file the cache is deleted.

void TTree::MoveReadCache(TFile *src, TDirectory *dir)
{
   if (!src) return;
   TFile *dst = (dir && dir != gROOT) ? dir->GetFile() : 0;
   if (src == dst) return;

   TTreeCache *pf = GetReadCache(src);
   if (dst) {
      src->SetCacheRead(0,this);
      dst->SetCacheRead(pf, this);
   } else {
      if (pf) {
         pf->WaitFinishPrefetch();
      }
      src->SetCacheRead(0,this);
      delete pf;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Function called when loading a new class library.

Bool_t TTree::Notify()
{
   TIter next(GetListOfLeaves());
   TLeaf* leaf = 0;
   while ((leaf = (TLeaf*) next())) {
      leaf->Notify();
      leaf->GetBranch()->Notify();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This function may be called after having filled some entries in a Tree.
/// Using the information in the existing branch buffers, it will reassign
/// new branch buffer sizes to optimize time and memory.
///
/// The function computes the best values for branch buffer sizes such that
/// the total buffer sizes is less than maxMemory and nearby entries written
/// at the same time.
/// In case the branch compression factor for the data written so far is less
/// than compMin, the compression is disabled.
///
/// if option ="d" an analysis report is printed.

void TTree::OptimizeBaskets(ULong64_t maxMemory, Float_t minComp, Option_t *option)
{
   //Flush existing baskets if the file is writable
   if (this->GetDirectory()->IsWritable()) this->FlushBasketsImpl();

   TString opt( option );
   opt.ToLower();
   Bool_t pDebug = opt.Contains("d");
   TObjArray *leaves = this->GetListOfLeaves();
   Int_t nleaves = leaves->GetEntries();
   Double_t treeSize = (Double_t)this->GetTotBytes();

   if (nleaves == 0 || treeSize == 0) {
      // We're being called too early, we really have nothing to do ...
      return;
   }
   Double_t aveSize = treeSize/nleaves;
   UInt_t bmin = 512;
   UInt_t bmax = 256000;
   Double_t memFactor = 1;
   Int_t i, oldMemsize,newMemsize,oldBaskets,newBaskets;
   i = oldMemsize = newMemsize = oldBaskets = newBaskets = 0;

   //we make two passes
   //one pass to compute the relative branch buffer sizes
   //a second pass to compute the absolute values
   for (Int_t pass =0;pass<2;pass++) {
      oldMemsize = 0;  //to count size of baskets in memory with old buffer size
      newMemsize = 0;  //to count size of baskets in memory with new buffer size
      oldBaskets = 0;  //to count number of baskets with old buffer size
      newBaskets = 0;  //to count number of baskets with new buffer size
      for (i=0;i<nleaves;i++) {
         TLeaf *leaf = (TLeaf*)leaves->At(i);
         TBranch *branch = leaf->GetBranch();
         Double_t totBytes = (Double_t)branch->GetTotBytes();
         Double_t idealFactor = totBytes/aveSize;
         UInt_t sizeOfOneEntry;
         if (branch->GetEntries() == 0) {
            // There is no data, so let's make a guess ...
            sizeOfOneEntry = aveSize;
         } else {
            sizeOfOneEntry = 1+(UInt_t)(totBytes / (Double_t)branch->GetEntries());
         }
         Int_t oldBsize = branch->GetBasketSize();
         oldMemsize += oldBsize;
         oldBaskets += 1+Int_t(totBytes/oldBsize);
         Int_t nb = branch->GetListOfBranches()->GetEntries();
         if (nb > 0) {
            newBaskets += 1+Int_t(totBytes/oldBsize);
            continue;
         }
         Double_t bsize = oldBsize*idealFactor*memFactor; //bsize can be very large !
         if (bsize < 0) bsize = bmax;
         if (bsize > bmax) bsize = bmax;
         UInt_t newBsize = UInt_t(bsize);
         if (pass) { // only on the second pass so that it doesn't interfere with scaling
            // If there is an entry offset, it will be stored in the same buffer as the object data; hence,
            // we must bump up the size of the branch to account for this extra footprint.
            // If fAutoFlush is not set yet, let's assume that it is 'in the process of being set' to
            // the value of GetEntries().
            Long64_t clusterSize = (fAutoFlush > 0) ? fAutoFlush : branch->GetEntries();
            if (branch->GetEntryOffsetLen()) {
               newBsize = newBsize + (clusterSize * sizeof(Int_t) * 2);
            }
            // We used ATLAS fully-split xAOD for testing, which is a rather unbalanced TTree, 10K branches,
            // with 8K having baskets smaller than 512 bytes. To achieve good I/O performance ATLAS uses auto-flush 100,
            // resulting in the smallest baskets being ~300-400 bytes, so this change increases their memory by about 8k*150B =~ 1MB,
            // at the same time it significantly reduces the number of total baskets because it ensures that all 100 entries can be
            // stored in a single basket (the old optimization tended to make baskets too small). In a toy example with fixed sized
            // structures we found a factor of 2 fewer baskets needed in the new scheme.
            // rounds up, increases basket size to ensure all entries fit into single basket as intended
            newBsize = newBsize - newBsize%512 + 512;
         }
         if (newBsize < sizeOfOneEntry) newBsize = sizeOfOneEntry;
         if (newBsize < bmin) newBsize = bmin;
         if (newBsize > 10000000) newBsize = bmax;
         if (pass) {
            if (pDebug) Info("OptimizeBaskets", "Changing buffer size from %6d to %6d bytes for %s\n",oldBsize,newBsize,branch->GetName());
            branch->SetBasketSize(newBsize);
         }
         newMemsize += newBsize;
         // For this number to be somewhat accurate when newBsize is 'low'
         // we do not include any space for meta data in the requested size (newBsize) even-though SetBasketSize will
         // not let it be lower than 100+TBranch::fEntryOffsetLen.
         newBaskets += 1+Int_t(totBytes/newBsize);
         if (pass == 0) continue;
         //Reset the compression level in case the compression factor is small
         Double_t comp = 1;
         if (branch->GetZipBytes() > 0) comp = totBytes/Double_t(branch->GetZipBytes());
         if (comp > 1 && comp < minComp) {
            if (pDebug) Info("OptimizeBaskets", "Disabling compression for branch : %s\n",branch->GetName());
            branch->SetCompressionSettings(ROOT::RCompressionSetting::EAlgorithm::kUseGlobal);
         }
      }
      // coverity[divide_by_zero] newMemsize can not be zero as there is at least one leaf
      memFactor = Double_t(maxMemory)/Double_t(newMemsize);
      if (memFactor > 100) memFactor = 100;
      Double_t bmin_new = bmin*memFactor;
      Double_t bmax_new = bmax*memFactor;
      static const UInt_t hardmax = 1*1024*1024*1024; // Really, really never give more than 1Gb to a single buffer.

      // Really, really never go lower than 8 bytes (we use this number
      // so that the calculation of the number of basket is consistent
      // but in fact SetBasketSize will not let the size go below
      // TBranch::fEntryOffsetLen + (100 + strlen(branch->GetName())
      // (The 2nd part being a slight over estimate of the key length.
      static const UInt_t hardmin = 8;
      bmin = (bmin_new > hardmax) ? hardmax : ( bmin_new < hardmin ? hardmin : (UInt_t)bmin_new );
      bmax = (bmax_new > hardmax) ? bmin : (UInt_t)bmax_new;
   }
   if (pDebug) {
      Info("OptimizeBaskets", "oldMemsize = %d,  newMemsize = %d\n",oldMemsize, newMemsize);
      Info("OptimizeBaskets", "oldBaskets = %d,  newBaskets = %d\n",oldBaskets, newBaskets);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to the Principal Components Analysis class.
///
/// Create an instance of TPrincipal
///
/// Fill it with the selected variables
///
/// - if option "n" is specified, the TPrincipal object is filled with
///                 normalized variables.
/// - If option "p" is specified, compute the principal components
/// - If option "p" and "d" print results of analysis
/// - If option "p" and "h" generate standard histograms
/// - If option "p" and "c" generate code of conversion functions
/// - return a pointer to the TPrincipal object. It is the user responsibility
/// - to delete this object.
/// - The option default value is "np"
///
/// see TTree::Draw for explanation of the other parameters.
///
/// The created object is  named "principal" and a reference to it
/// is added to the list of specials Root objects.
/// you can retrieve a pointer to the created object via:
/// ~~~ {.cpp}
///     TPrincipal *principal =
///     (TPrincipal*)gROOT->GetListOfSpecials()->FindObject("principal");
/// ~~~

TPrincipal* TTree::Principal(const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->Principal(varexp, selection, option, nentries, firstentry);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print a summary of the tree contents.
///
/// -  If option contains "all" friend trees are also printed.
/// -  If option contains "toponly" only the top level branches are printed.
/// -  If option contains "clusters" information about the cluster of baskets is printed.
///
/// Wildcarding can be used to print only a subset of the branches, e.g.,
/// `T.Print("Elec*")` will print all branches with name starting with "Elec".

void TTree::Print(Option_t* option) const
{
   // We already have been visited while recursively looking
   // through the friends tree, let's return.
   if (kPrint & fFriendLockStatus) {
      return;
   }
   Int_t s = 0;
   Int_t skey = 0;
   if (fDirectory) {
      TKey* key = fDirectory->GetKey(GetName());
      if (key) {
         skey = key->GetKeylen();
         s = key->GetNbytes();
      }
   }
   Long64_t total = skey;
   Long64_t zipBytes = GetZipBytes();
   if (zipBytes > 0) {
      total += GetTotBytes();
   }
   TBufferFile b(TBuffer::kWrite, 10000);
   TTree::Class()->WriteBuffer(b, (TTree*) this);
   total += b.Length();
   Long64_t file = zipBytes + s;
   Float_t cx = 1;
   if (zipBytes) {
      cx = (GetTotBytes() + 0.00001) / zipBytes;
   }
   Printf("******************************************************************************");
   Printf("*Tree    :%-10s: %-54s *", GetName(), GetTitle());
   Printf("*Entries : %8lld : Total = %15lld bytes  File  Size = %10lld *", fEntries, total, file);
   Printf("*        :          : Tree compression factor = %6.2f                       *", cx);
   Printf("******************************************************************************");

   // Avoid many check of option validity
   if (option == nullptr)
      option = "";

   if (strncmp(option,"clusters",strlen("clusters"))==0) {
      Printf("%-16s %-16s %-16s %8s %20s",
             "Cluster Range #", "Entry Start", "Last Entry", "Size", "Number of clusters");
      Int_t index= 0;
      Long64_t clusterRangeStart = 0;
      Long64_t totalClusters = 0;
      bool estimated = false;
      bool unknown = false;
      auto printer = [this, &totalClusters, &estimated, &unknown](Int_t ind, Long64_t start, Long64_t end, Long64_t recordedSize) {
            Long64_t nclusters = 0;
            if (recordedSize > 0) {
               nclusters = (1 + end - start) / recordedSize;
               Printf("%-16d %-16lld %-16lld %8lld %10lld",
                      ind, start, end, recordedSize, nclusters);
            } else {
               // NOTE: const_cast ... DO NOT Merge for now
               TClusterIterator iter((TTree*)this, start);
               iter.Next();
               auto estimated_size = iter.GetNextEntry() - start;
               if (estimated_size > 0) {
                  nclusters = (1 + end - start) / estimated_size;
                  Printf("%-16d %-16lld %-16lld %8lld %10lld (estimated)",
                      ind, start, end, recordedSize, nclusters);
                  estimated = true;
               } else {
                  Printf("%-16d %-16lld %-16lld %8lld    (unknown)",
                        ind, start, end, recordedSize);
                  unknown = true;
               }
            }
            start = end + 1;
            totalClusters += nclusters;
      };
      if (fNClusterRange) {
         for( ; index < fNClusterRange; ++index) {
            printer(index, clusterRangeStart, fClusterRangeEnd[index], fClusterSize[index]);
            clusterRangeStart = fClusterRangeEnd[index] + 1;
         }
      }
      printer(index, clusterRangeStart, fEntries - 1, fAutoFlush);
      if (unknown) {
         Printf("Total number of clusters: (unknown)");
      } else  {
         Printf("Total number of clusters: %lld %s", totalClusters, estimated ? "(estimated)" : "");
      }
      return;
   }

   Int_t nl = const_cast<TTree*>(this)->GetListOfLeaves()->GetEntries();
   Int_t l;
   TBranch* br = 0;
   TLeaf* leaf = 0;
   if (strstr(option, "toponly")) {
      Long64_t *count = new Long64_t[nl];
      Int_t keep =0;
      for (l=0;l<nl;l++) {
         leaf = (TLeaf *)const_cast<TTree*>(this)->GetListOfLeaves()->At(l);
         br   = leaf->GetBranch();
         if (strchr(br->GetName(),'.')) {
            count[l] = -1;
            count[keep] += br->GetZipBytes();
         } else {
            keep = l;
            count[keep]  = br->GetZipBytes();
         }
      }
      for (l=0;l<nl;l++) {
         if (count[l] < 0) continue;
         leaf = (TLeaf *)const_cast<TTree*>(this)->GetListOfLeaves()->At(l);
         br   = leaf->GetBranch();
         Printf("branch: %-20s %9lld\n",br->GetName(),count[l]);
      }
      delete [] count;
   } else {
      TString reg = "*";
      if (strlen(option) && strchr(option,'*')) reg = option;
      TRegexp re(reg,kTRUE);
      TIter next(const_cast<TTree*>(this)->GetListOfBranches());
      TBranch::ResetCount();
      while ((br= (TBranch*)next())) {
         TString st = br->GetName();
         st.ReplaceAll("/","_");
         if (st.Index(re) == kNPOS) continue;
         br->Print(option);
      }
   }

   //print TRefTable (if one)
   if (fBranchRef) fBranchRef->Print(option);

   //print friends if option "all"
   if (!fFriends || !strstr(option,"all")) return;
   TIter nextf(fFriends);
   TFriendLock lock(const_cast<TTree*>(this),kPrint);
   TFriendElement *fr;
   while ((fr = (TFriendElement*)nextf())) {
      TTree * t = fr->GetTree();
      if (t) t->Print(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print statistics about the TreeCache for this tree.
/// Like:
/// ~~~ {.cpp}
///     ******TreeCache statistics for file: cms2.root ******
///     Reading 73921562 bytes in 716 transactions
///     Average transaction = 103.242405 Kbytes
///     Number of blocks in current cache: 202, total size : 6001193
/// ~~~
/// if option = "a" the list of blocks in the cache is printed

void TTree::PrintCacheStats(Option_t* option) const
{
   TFile *f = GetCurrentFile();
   if (!f) return;
   TTreeCache *tc = GetReadCache(f);
   if (tc) tc->Print(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Process this tree executing the TSelector code in the specified filename.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.
///
/// The code in filename is loaded (interpreted or compiled, see below),
/// filename must contain a valid class implementation derived from TSelector,
/// where TSelector has the following member functions:
///
/// - `Begin()`:         called every time a loop on the tree starts,
///                      a convenient place to create your histograms.
/// - `SlaveBegin()`:    called after Begin(), when on PROOF called only on the
///                      slave servers.
/// - `Process()`:       called for each event, in this function you decide what
///                      to read and fill your histograms.
/// - `SlaveTerminate`:  called at the end of the loop on the tree, when on PROOF
///                      called only on the slave servers.
/// - `Terminate()`:     called at the end of the loop on the tree,
///                      a convenient place to draw/fit your histograms.
///
/// If filename is of the form file.C, the file will be interpreted.
///
/// If filename is of the form file.C++, the file file.C will be compiled
/// and dynamically loaded.
///
/// If filename is of the form file.C+, the file file.C will be compiled
/// and dynamically loaded. At next call, if file.C is older than file.o
/// and file.so, the file.C is not compiled, only file.so is loaded.
///
/// ## NOTE1
///
/// It may be more interesting to invoke directly the other Process function
/// accepting a TSelector* as argument.eg
/// ~~~ {.cpp}
///     MySelector *selector = (MySelector*)TSelector::GetSelector(filename);
///     selector->CallSomeFunction(..);
///     mytree.Process(selector,..);
/// ~~~
/// ## NOTE2
//
/// One should not call this function twice with the same selector file
/// in the same script. If this is required, proceed as indicated in NOTE1,
/// by getting a pointer to the corresponding TSelector,eg
///
/// ### Workaround 1
///
/// ~~~ {.cpp}
///     void stubs1() {
///        TSelector *selector = TSelector::GetSelector("h1test.C");
///        TFile *f1 = new TFile("stubs_nood_le1.root");
///        TTree *h1 = (TTree*)f1->Get("h1");
///        h1->Process(selector);
///        TFile *f2 = new TFile("stubs_nood_le1_coarse.root");
///        TTree *h2 = (TTree*)f2->Get("h1");
///        h2->Process(selector);
///     }
/// ~~~
/// or use ACLIC to compile the selector
///
/// ### Workaround 2
///
/// ~~~ {.cpp}
///     void stubs2() {
///        TFile *f1 = new TFile("stubs_nood_le1.root");
///        TTree *h1 = (TTree*)f1->Get("h1");
///        h1->Process("h1test.C+");
///        TFile *f2 = new TFile("stubs_nood_le1_coarse.root");
///        TTree *h2 = (TTree*)f2->Get("h1");
///        h2->Process("h1test.C+");
///     }
/// ~~~

Long64_t TTree::Process(const char* filename, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->Process(filename, option, nentries, firstentry);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Process this tree executing the code in the specified selector.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.
///
///   The TSelector class has the following member functions:
///
/// - `Begin()`:        called every time a loop on the tree starts,
///                     a convenient place to create your histograms.
/// - `SlaveBegin()`:   called after Begin(), when on PROOF called only on the
///                     slave servers.
/// - `Process()`:      called for each event, in this function you decide what
///                     to read and fill your histograms.
/// - `SlaveTerminate`: called at the end of the loop on the tree, when on PROOF
///                     called only on the slave servers.
/// - `Terminate()`:    called at the end of the loop on the tree,
///                     a convenient place to draw/fit your histograms.
///
///  If the Tree (Chain) has an associated EventList, the loop is on the nentries
///  of the EventList, starting at firstentry, otherwise the loop is on the
///  specified Tree entries.

Long64_t TTree::Process(TSelector* selector, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->Process(selector, option, nentries, firstentry);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a projection of a tree using selections.
///
/// Depending on the value of varexp (described in Draw) a 1-D, 2-D, etc.,
/// projection of the tree will be filled in histogram hname.
/// Note that the dimension of hname must match with the dimension of varexp.
///

Long64_t TTree::Project(const char* hname, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TString var;
   var.Form("%s>>%s", varexp, hname);
   TString opt("goff");
   if (option) {
      opt.Form("%sgoff", option);
   }
   Long64_t nsel = Draw(var, selection, opt, nentries, firstentry);
   return nsel;
}

////////////////////////////////////////////////////////////////////////////////
/// Loop over entries and return a TSQLResult object containing entries following selection.

TSQLResult* TTree::Query(const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->Query(varexp, selection, option, nentries, firstentry);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create or simply read branches from filename.
///
/// if branchDescriptor = "" (default), it is assumed that the Tree descriptor
/// is given in the first line of the file with a syntax like
/// ~~~ {.cpp}
///     A/D:Table[2]/F:Ntracks/I:astring/C
/// ~~~
/// otherwise branchDescriptor must be specified with the above syntax.
///
/// - If the type of the first variable is not specified, it is assumed to be "/F"
/// - If the type of any other variable is not specified, the type of the previous
///   variable is assumed. eg
///     - `x:y:z`      (all variables are assumed of type "F"
///     - `x/D:y:z`    (all variables are of type "D"
///     - `x:y/D:z`    (x is type "F", y and z of type "D"
///
/// delimiter allows for the use of another delimiter besides whitespace.
/// This provides support for direct import of common data file formats
/// like csv.  If delimiter != ' ' and branchDescriptor == "", then the
/// branch description is taken from the first line in the file, but
/// delimiter is used for the branch names tokenization rather than ':'.
/// Note however that if the values in the first line do not use the
/// /[type] syntax, all variables are assumed to be of type "F".
/// If the filename ends with extensions .csv or .CSV and a delimiter is
/// not specified (besides ' '), the delimiter is automatically set to ','.
///
/// Lines in the input file starting with "#" are ignored. Leading whitespace
/// for each column data is skipped. Empty lines are skipped.
///
/// A TBranch object is created for each variable in the expression.
/// The total number of rows read from the file is returned.
///
/// ## FILLING a TTree WITH MULTIPLE INPUT TEXT FILES
///
/// To fill a TTree with multiple input text files, proceed as indicated above
/// for the first input file and omit the second argument for subsequent calls
/// ~~~ {.cpp}
///      T.ReadFile("file1.dat","branch descriptor");
///      T.ReadFile("file2.dat");
/// ~~~

Long64_t TTree::ReadFile(const char* filename, const char* branchDescriptor, char delimiter)
{
   std::ifstream in;
   in.open(filename);
   if (!in.good()) {
      Error("ReadFile","Cannot open file: %s",filename);
      return 0;
   }
   const char* ext = strrchr(filename, '.');
   if(ext != NULL && ((strcmp(ext, ".csv") == 0) || (strcmp(ext, ".CSV") == 0)) && delimiter == ' ') {
      delimiter = ',';
   }
   return ReadStream(in, branchDescriptor, delimiter);
}

////////////////////////////////////////////////////////////////////////////////
/// Determine which newline this file is using.
/// Return '\\r' for Windows '\\r\\n' as that already terminates.

char TTree::GetNewlineValue(std::istream &inputStream)
{
   Long_t inPos = inputStream.tellg();
   char newline = '\n';
   while(1) {
      char c = 0;
      inputStream.get(c);
      if(!inputStream.good()) {
         Error("ReadStream","Error reading stream: no newline found.");
         return 0;
      }
      if(c == newline) break;
      if(c == '\r') {
         newline = '\r';
         break;
      }
   }
   inputStream.clear();
   inputStream.seekg(inPos);
   return newline;
}

////////////////////////////////////////////////////////////////////////////////
/// Create or simply read branches from an input stream.
///
/// See reference information for TTree::ReadFile

Long64_t TTree::ReadStream(std::istream& inputStream, const char *branchDescriptor, char delimiter)
{
   char newline = 0;
   std::stringstream ss;
   std::istream *inTemp;
   Long_t inPos = inputStream.tellg();
   if (!inputStream.good()) {
      Error("ReadStream","Error reading stream");
      return 0;
   }
   if (inPos == -1) {
      ss << std::cin.rdbuf();
      newline = GetNewlineValue(ss);
      inTemp = &ss;
   } else {
      newline = GetNewlineValue(inputStream);
      inTemp = &inputStream;
   }
   std::istream& in = *inTemp;
   Long64_t nlines = 0;

   TBranch *branch = 0;
   Int_t nbranches = fBranches.GetEntries();
   if (nbranches == 0) {
      char *bdname = new char[4000];
      char *bd = new char[100000];
      Int_t nch = 0;
      if (branchDescriptor) nch = strlen(branchDescriptor);
      // branch Descriptor is null, read its definition from the first line in the file
      if (!nch) {
         do {
            in.getline(bd, 100000, newline);
            if (!in.good()) {
               delete [] bdname;
               delete [] bd;
               Error("ReadStream","Error reading stream");
               return 0;
            }
            char *cursor = bd;
            while( isspace(*cursor) && *cursor != '\n' && *cursor != '\0') {
               ++cursor;
            }
            if (*cursor != '#' && *cursor != '\n' && *cursor != '\0') {
               break;
            }
         } while (true);
         ++nlines;
         nch = strlen(bd);
      } else {
         strlcpy(bd,branchDescriptor,100000);
      }

      //parse the branch descriptor and create a branch for each element
      //separated by ":"
      void *address = &bd[90000];
      char *bdcur = bd;
      TString desc="", olddesc="F";
      char bdelim = ':';
      if(delimiter != ' ') {
         bdelim = delimiter;
         if (strchr(bdcur,bdelim)==0 && strchr(bdcur,':') != 0) {
            // revert to the default
            bdelim = ':';
         }
      }
      while (bdcur) {
         char *colon = strchr(bdcur,bdelim);
         if (colon) *colon = 0;
         strlcpy(bdname,bdcur,4000);
         char *slash = strchr(bdname,'/');
         if (slash) {
            *slash = 0;
            desc = bdcur;
            olddesc = slash+1;
         } else {
            desc.Form("%s/%s",bdname,olddesc.Data());
         }
         char *bracket = strchr(bdname,'[');
         if (bracket) {
            *bracket = 0;
         }
         branch = new TBranch(this,bdname,address,desc.Data(),32000);
         if (branch->IsZombie()) {
            delete branch;
            Warning("ReadStream","Illegal branch definition: %s",bdcur);
         } else {
            fBranches.Add(branch);
            branch->SetAddress(0);
         }
         if (!colon)break;
         bdcur = colon+1;
      }
      delete [] bdname;
      delete [] bd;
   }

   nbranches = fBranches.GetEntries();

   if (gDebug > 1) {
      Info("ReadStream", "Will use branches:");
      for (int i = 0 ; i < nbranches; ++i) {
         TBranch* br = (TBranch*) fBranches.At(i);
         Info("ReadStream", "  %s: %s [%s]", br->GetName(),
              br->GetTitle(), br->GetListOfLeaves()->At(0)->IsA()->GetName());
      }
      if (gDebug > 3) {
         Info("ReadStream", "Dumping read tokens, format:");
         Info("ReadStream", "LLLLL:BBB:gfbe:GFBE:T");
         Info("ReadStream", "   L: line number");
         Info("ReadStream", "   B: branch number");
         Info("ReadStream", "   gfbe: good / fail / bad / eof of token");
         Info("ReadStream", "   GFBE: good / fail / bad / eof of file");
         Info("ReadStream", "   T: Token being read");
      }
   }

   //loop on all lines in the file
   Long64_t nGoodLines = 0;
   std::string line;
   const char sDelimBuf[2] = { delimiter, 0 };
   const char* sDelim = sDelimBuf;
   if (delimiter == ' ') {
      // ' ' really means whitespace
      sDelim = "[ \t]";
   }
   while(in.good()) {
      if (newline == '\r' && in.peek() == '\n') {
         // Windows, skip '\n':
         in.get();
      }
      std::getline(in, line, newline);
      ++nlines;

      TString sLine(line);
      sLine = sLine.Strip(TString::kLeading); // skip leading whitespace
      if (sLine.IsNull()) {
         if (gDebug > 2) {
            Info("ReadStream", "Skipping empty line number %lld", nlines);
         }
         continue; // silently skip empty lines
      }
      if (sLine[0] == '#') {
         if (gDebug > 2) {
            Info("ReadStream", "Skipping comment line number %lld: '%s'",
                 nlines, line.c_str());
         }
         continue;
      }
      if (gDebug > 2) {
         Info("ReadStream", "Parsing line number %lld: '%s'",
              nlines, line.c_str());
      }

      // Loop on branches and read the branch values into their buffer
      branch = 0;
      TString tok; // one column's data
      TString leafData; // leaf data, possibly multiple tokens for e.g. /I[2]
      std::stringstream sToken; // string stream feeding leafData into leaves
      Ssiz_t pos = 0;
      Int_t iBranch = 0;
      Bool_t goodLine = kTRUE; // whether the row can be filled into the tree
      Int_t remainingLeafLen = 0; // remaining columns for the current leaf
      while (goodLine && iBranch < nbranches
             && sLine.Tokenize(tok, pos, sDelim)) {
         tok = tok.Strip(TString::kLeading); // skip leading whitespace
         if (tok.IsNull() && delimiter == ' ') {
            // 1   2 should not be interpreted as 1,,,2 but 1, 2.
            // Thus continue until we have a non-empty token.
            continue;
         }

         if (!remainingLeafLen) {
            // next branch!
            branch = (TBranch*)fBranches.At(iBranch);
         }
         TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
         if (!remainingLeafLen) {
            remainingLeafLen = leaf->GetLen();
            if (leaf->GetMaximum() > 0) {
               // This is a dynamic leaf length, i.e. most likely a TLeafC's
               // string size. This still translates into one token:
               remainingLeafLen = 1;
            }

            leafData = tok;
         } else {
            // append token to laf data:
            leafData += " ";
            leafData += tok;
         }
         --remainingLeafLen;
         if (remainingLeafLen) {
            // need more columns for this branch:
            continue;
         }
         ++iBranch;

         // initialize stringstream with token
         sToken.clear();
         sToken.seekp(0, std::ios_base::beg);
         sToken.str(leafData.Data());
         sToken.seekg(0, std::ios_base::beg);
         leaf->ReadValue(sToken, 0 /* 0 = "all" */);
         if (gDebug > 3) {
            Info("ReadStream", "%5lld:%3d:%d%d%d%d:%d%d%d%d:%s",
                 nlines, iBranch,
                 (int)sToken.good(), (int)sToken.fail(),
                 (int)sToken.bad(), (int)sToken.eof(),
                 (int)in.good(), (int)in.fail(),
                 (int)in.bad(), (int)in.eof(),
                 sToken.str().c_str());
         }

         // Error handling
         if (sToken.bad()) {
            // How could that happen for a stringstream?
            Warning("ReadStream",
                    "Buffer error while reading data for branch %s on line %lld",
                    branch->GetName(), nlines);
         } else if (!sToken.eof()) {
            if (sToken.fail()) {
               Warning("ReadStream",
                       "Couldn't read formatted data in \"%s\" for branch %s on line %lld; ignoring line",
                       tok.Data(), branch->GetName(), nlines);
               goodLine = kFALSE;
            } else {
               std::string remainder;
               std::getline(sToken, remainder, newline);
               if (!remainder.empty()) {
                  Warning("ReadStream",
                          "Ignoring trailing \"%s\" while reading data for branch %s on line %lld",
                          remainder.c_str(), branch->GetName(), nlines);
               }
            }
         }
      } // tokenizer loop

      if (iBranch < nbranches) {
         Warning("ReadStream",
                 "Read too few columns (%d < %d) in line %lld; ignoring line",
                 iBranch, nbranches, nlines);
         goodLine = kFALSE;
      } else if (pos != kNPOS) {
         sLine = sLine.Strip(TString::kTrailing);
         if (pos < sLine.Length()) {
            Warning("ReadStream",
                    "Ignoring trailing \"%s\" while reading line %lld",
                    sLine.Data() + pos - 1 /* also print delimiter */,
                    nlines);
         }
      }

      //we are now ready to fill the tree
      if (goodLine) {
         Fill();
         ++nGoodLines;
      }
   }

   return nGoodLines;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that obj (which is being deleted or will soon be) is no
/// longer referenced by this TTree.

void TTree::RecursiveRemove(TObject *obj)
{
   if (obj == fEventList) {
      fEventList = 0;
   }
   if (obj == fEntryList) {
      fEntryList = 0;
   }
   if (fUserInfo) {
      fUserInfo->RecursiveRemove(obj);
   }
   if (fPlayer == obj) {
      fPlayer = 0;
   }
   if (fTreeIndex == obj) {
      fTreeIndex = 0;
   }
   if (fAliases) {
      fAliases->RecursiveRemove(obj);
   }
   if (fFriends) {
      fFriends->RecursiveRemove(obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Refresh contents of this tree and its branches from the current status on disk.
///
///  One can call this function in case the tree file is being
///  updated by another process.

void TTree::Refresh()
{
   if (!fDirectory->GetFile()) {
      return;
   }
   fDirectory->ReadKeys();
   fDirectory->Remove(this);
   TTree* tree; fDirectory->GetObject(GetName(),tree);
   if (!tree) {
      return;
   }
   //copy info from tree header into this Tree
   fEntries = 0;
   fNClusterRange = 0;
   ImportClusterRanges(tree);

   fAutoSave = tree->fAutoSave;
   fEntries = tree->fEntries;
   fTotBytes = tree->GetTotBytes();
   fZipBytes = tree->GetZipBytes();
   fSavedBytes = tree->fSavedBytes;
   fTotalBuffers = tree->fTotalBuffers.load();

   //loop on all branches and update them
   Int_t nleaves = fLeaves.GetEntriesFast();
   for (Int_t i = 0; i < nleaves; i++)  {
      TLeaf* leaf = (TLeaf*) fLeaves.UncheckedAt(i);
      TBranch* branch = (TBranch*) leaf->GetBranch();
      branch->Refresh(tree->GetBranch(branch->GetName()));
   }
   fDirectory->Remove(tree);
   fDirectory->Append(this);
   delete tree;
   tree = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Record a TFriendElement that we need to warn when the chain switches to
/// a new file (typically this is because this chain is a friend of another
/// TChain)

void TTree::RegisterExternalFriend(TFriendElement *fe)
{
   if (!fExternalFriends)
      fExternalFriends = new TList();
   fExternalFriends->Add(fe);
}


////////////////////////////////////////////////////////////////////////////////
/// Removes external friend

void TTree::RemoveExternalFriend(TFriendElement *fe)
{
   if (fExternalFriends) fExternalFriends->Remove((TObject*)fe);
}


////////////////////////////////////////////////////////////////////////////////
/// Remove a friend from the list of friends.

void TTree::RemoveFriend(TTree* oldFriend)
{
   // We already have been visited while recursively looking
   // through the friends tree, let return
   if (kRemoveFriend & fFriendLockStatus) {
      return;
   }
   if (!fFriends) {
      return;
   }
   TFriendLock lock(this, kRemoveFriend);
   TIter nextf(fFriends);
   TFriendElement* fe = 0;
   while ((fe = (TFriendElement*) nextf())) {
      TTree* friend_t = fe->GetTree();
      if (friend_t == oldFriend) {
         fFriends->Remove(fe);
         delete fe;
         fe = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset baskets, buffers and entries count in all branches and leaves.

void TTree::Reset(Option_t* option)
{
   fNotify        = 0;
   fEntries       = 0;
   fNClusterRange = 0;
   fTotBytes      = 0;
   fZipBytes      = 0;
   fFlushedBytes  = 0;
   fSavedBytes    = 0;
   fTotalBuffers  = 0;
   fChainOffset   = 0;
   fReadEntry     = -1;

   delete fTreeIndex;
   fTreeIndex = 0;

   Int_t nb = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i)  {
      TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
      branch->Reset(option);
   }

   if (fBranchRef) {
      fBranchRef->Reset();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resets the state of this TTree after a merge (keep the customization but
/// forget the data).

void TTree::ResetAfterMerge(TFileMergeInfo *info)
{
   fEntries       = 0;
   fNClusterRange = 0;
   fTotBytes      = 0;
   fZipBytes      = 0;
   fSavedBytes    = 0;
   fFlushedBytes  = 0;
   fTotalBuffers  = 0;
   fChainOffset   = 0;
   fReadEntry     = -1;

   delete fTreeIndex;
   fTreeIndex     = 0;

   Int_t nb = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i)  {
      TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
      branch->ResetAfterMerge(info);
   }

   if (fBranchRef) {
      fBranchRef->ResetAfterMerge(info);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Tell all of our branches to set their addresses to zero.
///
/// Note: If any of our branches own any objects, they are deleted.

void TTree::ResetBranchAddress(TBranch *br)
{
   if (br && br->GetTree()) {
      br->ResetAddress();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Tell all of our branches to drop their current objects and allocate new ones.

void TTree::ResetBranchAddresses()
{
   TObjArray* branches = GetListOfBranches();
   Int_t nbranches = branches->GetEntriesFast();
   for (Int_t i = 0; i < nbranches; ++i) {
      TBranch* branch = (TBranch*) branches->UncheckedAt(i);
      branch->ResetAddress();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Loop over tree entries and print entries passing selection.
///
/// - If varexp is 0 (or "") then print only first 8 columns.
/// - If varexp = "*" print all columns.
///
/// Otherwise a columns selection can be made using "var1:var2:var3".
/// See TTreePlayer::Scan for more information

Long64_t TTree::Scan(const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->Scan(varexp, selection, option, nentries, firstentry);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a tree variable alias.
///
/// Set an alias for an expression/formula based on the tree 'variables'.
///
/// The content of 'aliasName' can be used in TTreeFormula (i.e. TTree::Draw,
/// TTree::Scan, TTreeViewer) and will be evaluated as the content of
/// 'aliasFormula'.
///
/// If the content of 'aliasFormula' only contains symbol names, periods and
/// array index specification (for example event.fTracks[3]), then
/// the content of 'aliasName' can be used as the start of symbol.
///
/// If the alias 'aliasName' already existed, it is replaced by the new
/// value.
///
/// When being used, the alias can be preceded by an eventual 'Friend Alias'
/// (see TTree::GetFriendAlias)
///
/// Return true if it was added properly.
///
/// For example:
/// ~~~ {.cpp}
///     tree->SetAlias("x1","(tdc1[1]-tdc1[0])/49");
///     tree->SetAlias("y1","(tdc1[3]-tdc1[2])/47");
///     tree->SetAlias("x2","(tdc2[1]-tdc2[0])/49");
///     tree->SetAlias("y2","(tdc2[3]-tdc2[2])/47");
///     tree->Draw("y2-y1:x2-x1");
///
///     tree->SetAlias("theGoodTrack","event.fTracks[3]");
///     tree->Draw("theGoodTrack.fPx"); // same as "event.fTracks[3].fPx"
/// ~~~

Bool_t TTree::SetAlias(const char* aliasName, const char* aliasFormula)
{
   if (!aliasName || !aliasFormula) {
      return kFALSE;
   }
   if (!aliasName[0] || !aliasFormula[0]) {
      return kFALSE;
   }
   if (!fAliases) {
      fAliases = new TList;
   } else {
      TNamed* oldHolder = (TNamed*) fAliases->FindObject(aliasName);
      if (oldHolder) {
         oldHolder->SetTitle(aliasFormula);
         return kTRUE;
      }
   }
   TNamed* holder = new TNamed(aliasName, aliasFormula);
   fAliases->Add(holder);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This function may be called at the start of a program to change
/// the default value for fAutoFlush.
///
/// ### CASE 1 : autof > 0
///
/// autof is the number of consecutive entries after which TTree::Fill will
/// flush all branch buffers to disk.
///
/// ### CASE 2 : autof < 0
///
/// When filling the Tree the branch buffers will be flushed to disk when
/// more than autof bytes have been written to the file. At the first FlushBaskets
/// TTree::Fill will replace fAutoFlush by the current value of fEntries.
///
/// Calling this function with autof<0 is interesting when it is hard to estimate
/// the size of one entry. This value is also independent of the Tree.
///
/// The Tree is initialized with fAutoFlush=-30000000, ie that, by default,
/// the first AutoFlush will be done when 30 MBytes of data are written to the file.
///
/// ### CASE 3 : autof = 0
///
/// The AutoFlush mechanism is disabled.
///
/// Flushing the buffers at regular intervals optimize the location of
/// consecutive entries on the disk by creating clusters of baskets.
///
/// A cluster of baskets is a set of baskets that contains all
/// the data for a (consecutive) set of entries and that is stored
/// consecutively on the disk.   When reading all the branches, this
/// is the minimum set of baskets that the TTreeCache will read.

void TTree::SetAutoFlush(Long64_t autof /* = -30000000 */ )
{
   // Implementation note:
   //
   // A positive value of autoflush determines the size (in number of entries) of
   // a cluster of baskets.
   //
   // If the value of autoflush is changed over time (this happens in
   // particular when the TTree results from fast merging many trees),
   // we record the values of fAutoFlush in the data members:
   //     fClusterRangeEnd and fClusterSize.
   // In the code we refer to a range of entries where the size of the
   // cluster of baskets is the same (i.e the value of AutoFlush was
   // constant) is called a ClusterRange.
   //
   // The 2 arrays (fClusterRangeEnd and fClusterSize) have fNClusterRange
   // active (used) values and have fMaxClusterRange allocated entries.
   //
   // fClusterRangeEnd contains the last entries number of a cluster range.
   // In particular this means that the 'next' cluster starts at fClusterRangeEnd[]+1
   // fClusterSize contains the size in number of entries of all the cluster
   // within the given range.
   // The last range (and the only one if fNClusterRange is zero) start at
   // fNClusterRange[fNClusterRange-1]+1 and ends at the end of the TTree.  The
   // size of the cluster in this range is given by the value of fAutoFlush.
   //
   // For example printing the beginning and end of each the ranges can be done by:
   //
   //   Printf("%-16s %-16s %-16s %5s",
   //          "Cluster Range #", "Entry Start", "Last Entry", "Size");
   //   Int_t index= 0;
   //   Long64_t clusterRangeStart = 0;
   //   if (fNClusterRange) {
   //      for( ; index < fNClusterRange; ++index) {
   //         Printf("%-16d %-16lld %-16lld %5lld",
   //                index, clusterRangeStart, fClusterRangeEnd[index], fClusterSize[index]);
   //         clusterRangeStart = fClusterRangeEnd[index] + 1;
   //      }
   //   }
   //   Printf("%-16d %-16lld %-16lld %5lld",
   //          index, prevEntry, fEntries - 1, fAutoFlush);
   //

   // Note:  We store the entry number corresponding to the end of the cluster
   // rather than its start in order to avoid using the array if the cluster
   // size never varies (If there is only one value of AutoFlush for the whole TTree).

   if( fAutoFlush != autof) {
      if ((fAutoFlush > 0 || autof > 0) && fFlushedBytes) {
         // The mechanism was already enabled, let's record the previous
         // cluster if needed.
         MarkEventCluster();
      }
      fAutoFlush = autof;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Mark the previous event as being at the end of the event cluster.
///
/// So, if fEntries is set to 10 (and this is the first cluster) when MarkEventCluster
/// is called, then the first cluster has 9 events.
void TTree::MarkEventCluster()
{
    if (!fEntries) return;

    if ( (fNClusterRange+1) > fMaxClusterRange ) {
        if (fMaxClusterRange) {
            // Resize arrays to hold a larger event cluster.
            Int_t newsize = TMath::Max(10,Int_t(2*fMaxClusterRange));
            fClusterRangeEnd = (Long64_t*)TStorage::ReAlloc(fClusterRangeEnd,
                                                            newsize*sizeof(Long64_t),fMaxClusterRange*sizeof(Long64_t));
            fClusterSize = (Long64_t*)TStorage::ReAlloc(fClusterSize,
                                                        newsize*sizeof(Long64_t),fMaxClusterRange*sizeof(Long64_t));
            fMaxClusterRange = newsize;
        } else {
            // Cluster ranges have never been initialized; create them now.
            fMaxClusterRange = 2;
            fClusterRangeEnd = new Long64_t[fMaxClusterRange];
            fClusterSize = new Long64_t[fMaxClusterRange];
        }
    }
    fClusterRangeEnd[fNClusterRange] = fEntries - 1;
    // If we are auto-flushing, then the cluster size is the same as the current auto-flush setting.
    if (fAutoFlush > 0) {
        // Even if the user triggers MarkEventRange prior to fAutoFlush being present, the TClusterIterator
        // will appropriately go to the next event range.
        fClusterSize[fNClusterRange] = fAutoFlush;
    // Otherwise, assume there is one cluster per event range (e.g., user is manually controlling the flush).
    } else if (fNClusterRange == 0) {
        fClusterSize[fNClusterRange] = fEntries;
    } else {
        fClusterSize[fNClusterRange] = fClusterRangeEnd[fNClusterRange] - fClusterRangeEnd[fNClusterRange-1];
    }
    ++fNClusterRange;
}

////////////////////////////////////////////////////////////////////////////////
/// This function may be called at the start of a program to change
/// the default value for fAutoSave (and for SetAutoSave) is -300000000, ie 300 MBytes.
/// When filling the Tree the branch buffers as well as the Tree header
/// will be flushed to disk when the watermark is reached.
/// If fAutoSave is positive the watermark is reached when a multiple of fAutoSave
/// entries have been written.
/// If fAutoSave is negative the watermark is reached when -fAutoSave bytes
/// have been written to the file.
/// In case of a program crash, it will be possible to recover the data in the Tree
/// up to the last AutoSave point.

void TTree::SetAutoSave(Long64_t autos)
{
   fAutoSave = autos;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a branch's basket size.
///
/// bname is the name of a branch.
///
/// - if bname="*", apply to all branches.
/// - if bname="xxx*", apply to all branches with name starting with xxx
///
/// see TRegexp for wildcarding options
/// buffsize = branc basket size

void TTree::SetBasketSize(const char* bname, Int_t buffsize)
{
   Int_t nleaves = fLeaves.GetEntriesFast();
   TRegexp re(bname, kTRUE);
   Int_t nb = 0;
   for (Int_t i = 0; i < nleaves; i++)  {
      TLeaf* leaf = (TLeaf*) fLeaves.UncheckedAt(i);
      TBranch* branch = (TBranch*) leaf->GetBranch();
      TString s = branch->GetName();
      if (strcmp(bname, branch->GetName()) && (s.Index(re) == kNPOS)) {
         continue;
      }
      nb++;
      branch->SetBasketSize(buffsize);
   }
   if (!nb) {
      Error("SetBasketSize", "unknown branch -> '%s'", bname);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change branch address, dealing with clone trees properly.
/// See TTree::CheckBranchAddressType for the semantic of the return value.
///
/// Note: See the comments in TBranchElement::SetAddress() for the
/// meaning of the addr parameter and the object ownership policy.

Int_t TTree::SetBranchAddress(const char* bname, void* addr, TBranch** ptr)
{
   TBranch* branch = GetBranch(bname);
   if (!branch) {
      if (ptr) *ptr = 0;
      Error("SetBranchAddress", "unknown branch -> %s", bname);
      return kMissingBranch;
   }
   return SetBranchAddressImp(branch,addr,ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Verify the validity of the type of addr before calling SetBranchAddress.
/// See TTree::CheckBranchAddressType for the semantic of the return value.
///
/// Note: See the comments in TBranchElement::SetAddress() for the
/// meaning of the addr parameter and the object ownership policy.

Int_t TTree::SetBranchAddress(const char* bname, void* addr, TClass* ptrClass, EDataType datatype, Bool_t isptr)
{
   return SetBranchAddress(bname, addr, 0, ptrClass, datatype, isptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Verify the validity of the type of addr before calling SetBranchAddress.
/// See TTree::CheckBranchAddressType for the semantic of the return value.
///
/// Note: See the comments in TBranchElement::SetAddress() for the
/// meaning of the addr parameter and the object ownership policy.

Int_t TTree::SetBranchAddress(const char* bname, void* addr, TBranch** ptr, TClass* ptrClass, EDataType datatype, Bool_t isptr)
{
   TBranch* branch = GetBranch(bname);
   if (!branch) {
      if (ptr) *ptr = 0;
      Error("SetBranchAddress", "unknown branch -> %s", bname);
      return kMissingBranch;
   }

   Int_t res = CheckBranchAddressType(branch, ptrClass, datatype, isptr);

   // This will set the value of *ptr to branch.
   if (res >= 0) {
      // The check succeeded.
      if ((res & kNeedEnableDecomposedObj) && !branch->GetMakeClass())
         branch->SetMakeClass(kTRUE);
      SetBranchAddressImp(branch,addr,ptr);
   } else {
      if (ptr) *ptr = 0;
   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Change branch address, dealing with clone trees properly.
/// See TTree::CheckBranchAddressType for the semantic of the return value.
///
/// Note: See the comments in TBranchElement::SetAddress() for the
/// meaning of the addr parameter and the object ownership policy.

Int_t TTree::SetBranchAddressImp(TBranch *branch, void* addr, TBranch** ptr)
{
   if (ptr) {
      *ptr = branch;
   }
   if (fClones) {
      void* oldAddr = branch->GetAddress();
      TIter next(fClones);
      TTree* clone = 0;
      const char *bname = branch->GetName();
      while ((clone = (TTree*) next())) {
         TBranch* cloneBr = clone->GetBranch(bname);
         if (cloneBr && (cloneBr->GetAddress() == oldAddr)) {
            cloneBr->SetAddress(addr);
         }
      }
   }
   branch->SetAddress(addr);
   return kVoidPtr;
}

////////////////////////////////////////////////////////////////////////////////
/// Set branch status to Process or DoNotProcess.
///
/// When reading a Tree, by default, all branches are read.
/// One can speed up considerably the analysis phase by activating
/// only the branches that hold variables involved in a query.
///
/// bname is the name of a branch.
///
/// - if bname="*", apply to all branches.
/// - if bname="xxx*", apply to all branches with name starting with xxx
///
/// see TRegexp for wildcarding options
///
/// - status = 1  branch will be processed
/// - = 0  branch will not be processed
///
/// Example:
///
/// Assume a tree T with sub-branches a,b,c,d,e,f,g,etc..
/// when doing T.GetEntry(i) all branches are read for entry i.
/// to read only the branches c and e, one can do
/// ~~~ {.cpp}
///     T.SetBranchStatus("*",0); //disable all branches
///     T.SetBranchStatus("c",1);
///     T.setBranchStatus("e",1);
///     T.GetEntry(i);
/// ~~~
/// bname is interpreted as a wild-carded TRegexp (see TRegexp::MakeWildcard).
/// Thus, "a*b" or "a.*b" matches branches starting with "a" and ending with
/// "b", but not any other branch with an "a" followed at some point by a
/// "b". For this second behavior, use "*a*b*". Note that TRegExp does not
/// support '|', and so you cannot select, e.g. track and shower branches
/// with "track|shower".
///
/// __WARNING! WARNING! WARNING!__
///
/// SetBranchStatus is matching the branch based on match of the branch
/// 'name' and not on the branch hierarchy! In order to be able to
/// selectively enable a top level object that is 'split' you need to make
/// sure the name of the top level branch is prefixed to the sub-branches'
/// name (by adding a dot ('.') at the end of the Branch creation and use the
/// corresponding bname.
///
/// I.e If your Tree has been created in split mode with a parent branch "parent."
/// (note the trailing dot).
/// ~~~ {.cpp}
///     T.SetBranchStatus("parent",1);
/// ~~~
/// will not activate the sub-branches of "parent". You should do:
/// ~~~ {.cpp}
///     T.SetBranchStatus("parent*",1);
/// ~~~
/// Without the trailing dot in the branch creation you have no choice but to
/// call SetBranchStatus explicitly for each of the sub branches.
///
/// An alternative to this function is to read directly and only
/// the interesting branches. Example:
/// ~~~ {.cpp}
///     TBranch *brc = T.GetBranch("c");
///     TBranch *bre = T.GetBranch("e");
///     brc->GetEntry(i);
///     bre->GetEntry(i);
/// ~~~
/// If found is not 0, the number of branch(es) found matching the regular
/// expression is returned in *found AND the error message 'unknown branch'
/// is suppressed.

void TTree::SetBranchStatus(const char* bname, Bool_t status, UInt_t* found)
{
   // We already have been visited while recursively looking
   // through the friends tree, let return
   if (kSetBranchStatus & fFriendLockStatus) {
      return;
   }

   if (0 == strcmp(bname, "")) {
      Error("SetBranchStatus", "Input regexp is an empty string: no match against branch names will be attempted.");
      return;
   }

   TBranch *branch, *bcount, *bson;
   TLeaf *leaf, *leafcount;

   Int_t i,j;
   Int_t nleaves = fLeaves.GetEntriesFast();
   TRegexp re(bname,kTRUE);
   Int_t nb = 0;

   // first pass, loop on all branches
   // for leafcount branches activate/deactivate in function of status
   for (i=0;i<nleaves;i++)  {
      leaf = (TLeaf*)fLeaves.UncheckedAt(i);
      branch = (TBranch*)leaf->GetBranch();
      TString s = branch->GetName();
      if (strcmp(bname,"*")) { //Regexp gives wrong result for [] in name
         TString longname;
         longname.Form("%s.%s",GetName(),branch->GetName());
         if (strcmp(bname,branch->GetName())
             && longname != bname
             && s.Index(re) == kNPOS) continue;
      }
      nb++;
      if (status) branch->ResetBit(kDoNotProcess);
      else        branch->SetBit(kDoNotProcess);
      leafcount = leaf->GetLeafCount();
      if (leafcount) {
         bcount = leafcount->GetBranch();
         if (status) bcount->ResetBit(kDoNotProcess);
         else        bcount->SetBit(kDoNotProcess);
      }
   }
   if (nb==0 && strchr(bname,'*')==0) {
      branch = GetBranch(bname);
      if (branch) {
         if (status) branch->ResetBit(kDoNotProcess);
         else        branch->SetBit(kDoNotProcess);
         ++nb;
      }
   }

   //search in list of friends
   UInt_t foundInFriend = 0;
   if (fFriends) {
      TFriendLock lock(this,kSetBranchStatus);
      TIter nextf(fFriends);
      TFriendElement *fe;
      TString name;
      while ((fe = (TFriendElement*)nextf())) {
         TTree *t = fe->GetTree();
         if (t==0) continue;

         // If the alias is present replace it with the real name.
         char *subbranch = (char*)strstr(bname,fe->GetName());
         if (subbranch!=bname) subbranch = 0;
         if (subbranch) {
            subbranch += strlen(fe->GetName());
            if ( *subbranch != '.' ) subbranch = 0;
            else subbranch ++;
         }
         if (subbranch) {
            name.Form("%s.%s",t->GetName(),subbranch);
         } else {
            name = bname;
         }
         t->SetBranchStatus(name,status, &foundInFriend);
      }
   }
   if (!nb && !foundInFriend) {
      if (found==0) {
         if (status) {
            if (strchr(bname,'*') != 0)
               Error("SetBranchStatus", "No branch name is matching wildcard -> %s", bname);
            else
               Error("SetBranchStatus", "unknown branch -> %s", bname);
         } else {
            if (strchr(bname,'*') != 0)
               Warning("SetBranchStatus", "No branch name is matching wildcard -> %s", bname);
            else
               Warning("SetBranchStatus", "unknown branch -> %s", bname);
         }
      }
      return;
   }
   if (found) *found = nb + foundInFriend;

   // second pass, loop again on all branches
   // activate leafcount branches for active branches only
   for (i = 0; i < nleaves; i++) {
      leaf = (TLeaf*)fLeaves.UncheckedAt(i);
      branch = (TBranch*)leaf->GetBranch();
      if (!branch->TestBit(kDoNotProcess)) {
         leafcount = leaf->GetLeafCount();
         if (leafcount) {
            bcount = leafcount->GetBranch();
            bcount->ResetBit(kDoNotProcess);
         }
      } else {
         //Int_t nbranches = branch->GetListOfBranches()->GetEntriesFast();
         Int_t nbranches = branch->GetListOfBranches()->GetEntries();
         for (j=0;j<nbranches;j++) {
            bson = (TBranch*)branch->GetListOfBranches()->UncheckedAt(j);
            if (!bson) continue;
            if (!bson->TestBit(kDoNotProcess)) {
               if (bson->GetNleaves() <= 0) continue;
               branch->ResetBit(kDoNotProcess);
               break;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current branch style.  (static function)
///
/// - style = 0 old Branch
/// - style = 1 new Bronch

void TTree::SetBranchStyle(Int_t style)
{
   fgBranchStyle = style;
}

////////////////////////////////////////////////////////////////////////////////
/// Set maximum size of the file cache .
//
/// - if cachesize = 0 the existing cache (if any) is deleted.
/// - if cachesize = -1 (default) it is set to the AutoFlush value when writing
///    the Tree (default is 30 MBytes).
///
/// Returns:
/// - 0 size set, cache was created if possible
/// - -1 on error

Int_t TTree::SetCacheSize(Long64_t cacheSize)
{
   // remember that the user has requested an explicit cache setup
   fCacheUserSet = kTRUE;

   return SetCacheSizeAux(kFALSE, cacheSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the size of the file cache and create it if possible.
///
/// If autocache is true:
/// this may be an autocreated cache, possibly enlarging an existing
/// autocreated cache. The size is calculated. The value passed in cacheSize:
/// - cacheSize =  0  make cache if default cache creation is enabled
/// - cacheSize = -1  make a default sized cache in any case
///
/// If autocache is false:
/// this is a user requested cache. cacheSize is used to size the cache.
/// This cache should never be automatically adjusted.
///
/// Returns:
/// - 0 size set, or existing autosized cache almost large enough.
///   (cache was created if possible)
/// - -1 on error

Int_t TTree::SetCacheSizeAux(Bool_t autocache /* = kTRUE */, Long64_t cacheSize /* = 0 */ )
{
   if (autocache) {
      // used as a once only control for automatic cache setup
      fCacheDoAutoInit = kFALSE;
   }

   if (!autocache) {
      // negative size means the user requests the default
      if (cacheSize < 0) {
         cacheSize = GetCacheAutoSize(kTRUE);
      }
   } else {
      if (cacheSize == 0) {
         cacheSize = GetCacheAutoSize();
      } else if (cacheSize < 0) {
         cacheSize = GetCacheAutoSize(kTRUE);
      }
   }

   TFile* file = GetCurrentFile();
   if (!file || GetTree() != this) {
      // if there's no file or we are not a plain tree (e.g. if we're a TChain)
      // do not create a cache, only record the size if one was given
      if (!autocache) {
         fCacheSize = cacheSize;
      }
      if (GetTree() != this) {
         return 0;
      }
      if (!autocache && cacheSize>0) {
         Warning("SetCacheSizeAux", "A TTreeCache could not be created because the TTree has no file");
      }
      return 0;
   }

   // Check for an existing cache
   TTreeCache* pf = GetReadCache(file);
   if (pf) {
      if (autocache) {
         // reset our cache status tracking in case existing cache was added
         // by the user without using one of the TTree methods
         fCacheSize = pf->GetBufferSize();
         fCacheUserSet = !pf->IsAutoCreated();

         if (fCacheUserSet) {
            // existing cache was created by the user, don't change it
            return 0;
         }
      } else {
         // update the cache to ensure it records the user has explicitly
         // requested it
         pf->SetAutoCreated(kFALSE);
      }

      // if we're using an automatically calculated size and the existing
      // cache is already almost large enough don't resize
      if (autocache && Long64_t(0.80*cacheSize) < fCacheSize) {
         // already large enough
         return 0;
      }

      if (cacheSize == fCacheSize) {
         return 0;
      }

      if (cacheSize == 0) {
         // delete existing cache
         pf->WaitFinishPrefetch();
         file->SetCacheRead(0,this);
         delete pf;
         pf = 0;
      } else {
         // resize
         Int_t res = pf->SetBufferSize(cacheSize);
         if (res < 0) {
            return -1;
         }
      }
   } else {
      // no existing cache
      if (autocache) {
         if (fCacheUserSet) {
            // value was already set manually.
            if (fCacheSize == 0) return 0;
            // Expected a cache should exist; perhaps the user moved it
            // Do nothing more here.
            if (cacheSize) {
               Error("SetCacheSizeAux", "Not setting up an automatically sized TTreeCache because of missing cache previously set");
            }
            return -1;
         }
      }
   }

   fCacheSize = cacheSize;
   if (cacheSize == 0 || pf) {
      return 0;
   }

#ifdef R__USE_IMT
   if(TTreeCacheUnzip::IsParallelUnzip() && file->GetCompressionLevel() > 0)
      pf = new TTreeCacheUnzip(this, cacheSize);
   else
#endif
      pf = new TTreeCache(this, cacheSize);

   pf->SetAutoCreated(autocache);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///interface to TTreeCache to set the cache entry range
///
/// Returns:
/// - 0 entry range set
/// - -1 on error

Int_t TTree::SetCacheEntryRange(Long64_t first, Long64_t last)
{
   if (!GetTree()) {
      if (LoadTree(0)<0) {
         Error("SetCacheEntryRange","Could not load a tree");
         return -1;
      }
   }
   if (GetTree()) {
      if (GetTree() != this) {
         return GetTree()->SetCacheEntryRange(first, last);
      }
   } else {
      Error("SetCacheEntryRange", "No tree is available. Could not set cache entry range");
      return -1;
   }

   TFile *f = GetCurrentFile();
   if (!f) {
      Error("SetCacheEntryRange", "No file is available. Could not set cache entry range");
      return -1;
   }
   TTreeCache *tc = GetReadCache(f,kTRUE);
   if (!tc) {
      Error("SetCacheEntryRange", "No cache is available. Could not set entry range");
      return -1;
   }
   tc->SetEntryRange(first,last);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to TTreeCache to set the number of entries for the learning phase

void TTree::SetCacheLearnEntries(Int_t n)
{
   TTreeCache::SetLearnEntries(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable circularity for this tree.
///
/// if maxEntries > 0 a maximum of maxEntries is kept in one buffer/basket
/// per branch in memory.
///   Note that when this function is called (maxEntries>0) the Tree
///   must be empty or having only one basket per branch.
/// if maxEntries <= 0 the tree circularity is disabled.
///
/// #### NOTE 1:
///  Circular Trees are interesting in online real time environments
///  to store the results of the last maxEntries events.
/// #### NOTE 2:
///  Calling SetCircular with maxEntries <= 0 is necessary before
///  merging circular Trees that have been saved on files.
/// #### NOTE 3:
///  SetCircular with maxEntries <= 0 is automatically called
///  by TChain::Merge
/// #### NOTE 4:
///  A circular Tree can still be saved in a file. When read back,
///  it is still a circular Tree and can be filled again.

void TTree::SetCircular(Long64_t maxEntries)
{
   if (maxEntries <= 0) {
      // Disable circularity.
      fMaxEntries = 1000000000;
      fMaxEntries *= 1000;
      ResetBit(kCircular);
      //in case the Tree was originally created in gROOT, the branch
      //compression level was set to -1. If the Tree is now associated to
      //a file, reset the compression level to the file compression level
      if (fDirectory) {
         TFile* bfile = fDirectory->GetFile();
         Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault;
         if (bfile) {
            compress = bfile->GetCompressionSettings();
         }
         Int_t nb = fBranches.GetEntriesFast();
         for (Int_t i = 0; i < nb; i++) {
            TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
            branch->SetCompressionSettings(compress);
         }
      }
   } else {
      // Enable circularity.
      fMaxEntries = maxEntries;
      SetBit(kCircular);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the debug level and the debug range.
///
/// For entries in the debug range, the functions TBranchElement::Fill
/// and TBranchElement::GetEntry will print the number of bytes filled
/// or read for each branch.

void TTree::SetDebug(Int_t level, Long64_t min, Long64_t max)
{
   fDebug = level;
   fDebugMin = min;
   fDebugMax = max;
}

////////////////////////////////////////////////////////////////////////////////
/// Update the default value for the branch's fEntryOffsetLen.
/// If updateExisting is true, also update all the existing branches.
/// If newdefault is less than 10, the new default value will be 10.

void TTree::SetDefaultEntryOffsetLen(Int_t newdefault, Bool_t updateExisting)
{
   if (newdefault < 10) {
      newdefault = 10;
   }
   fDefaultEntryOffsetLen = newdefault;
   if (updateExisting) {
      TIter next( GetListOfBranches() );
      TBranch *b;
      while ( ( b = (TBranch*)next() ) ) {
         b->SetEntryOffsetLen( newdefault, kTRUE );
      }
      if (fBranchRef) {
         fBranchRef->SetEntryOffsetLen( newdefault, kTRUE );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change the tree's directory.
///
/// Remove reference to this tree from current directory and
/// add reference to new directory dir.  The dir parameter can
/// be 0 in which case the tree does not belong to any directory.
///

void TTree::SetDirectory(TDirectory* dir)
{
   if (fDirectory == dir) {
      return;
   }
   if (fDirectory) {
      fDirectory->Remove(this);

      // Delete or move the file cache if it points to this Tree
      TFile *file = fDirectory->GetFile();
      MoveReadCache(file,dir);
   }
   fDirectory = dir;
   if (fDirectory) {
      fDirectory->Append(this);
   }
   TFile* file = 0;
   if (fDirectory) {
      file = fDirectory->GetFile();
   }
   if (fBranchRef) {
      fBranchRef->SetFile(file);
   }
   TBranch* b = 0;
   TIter next(GetListOfBranches());
   while((b = (TBranch*) next())) {
      b->SetFile(file);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change number of entries in the tree.
///
/// If n >= 0, set number of entries in the tree = n.
///
/// If n < 0, set number of entries in the tree to match the
/// number of entries in each branch. (default for n is -1)
///
/// This function should be called only when one fills each branch
/// independently via TBranch::Fill without calling TTree::Fill.
/// Calling TTree::SetEntries() make sense only if the number of entries
/// in each branch is identical, a warning is issued otherwise.
/// The function returns the number of entries.
///

Long64_t TTree::SetEntries(Long64_t n)
{
   // case 1 : force number of entries to n
   if (n >= 0) {
      fEntries = n;
      return n;
   }

   // case 2; compute the number of entries from the number of entries in the branches
   TBranch* b(nullptr), *bMin(nullptr), *bMax(nullptr);
   Long64_t nMin = kMaxEntries;
   Long64_t nMax = 0;
   TIter next(GetListOfBranches());
   while((b = (TBranch*) next())){
      Long64_t n2 = b->GetEntries();
      if (!bMin || n2 < nMin) {
         nMin = n2;
         bMin = b;
      }
      if (!bMax || n2 > nMax) {
         nMax = n2;
         bMax = b;
      }
   }
   if (bMin && nMin != nMax) {
      Warning("SetEntries", "Tree branches have different numbers of entries, eg %s has %lld entries while %s has %lld entries.",
              bMin->GetName(), nMin, bMax->GetName(), nMax);
   }
   fEntries = nMax;
   return fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Set an EntryList

void TTree::SetEntryList(TEntryList *enlist, Option_t * /*opt*/)
{
   if (fEntryList) {
      //check if the previous entry list is owned by the tree
      if (fEntryList->TestBit(kCanDelete)){
         delete fEntryList;
      }
   }
   fEventList = 0;
   if (!enlist) {
      fEntryList = 0;
      return;
   }
   fEntryList = enlist;
   fEntryList->SetTree(this);

}

////////////////////////////////////////////////////////////////////////////////
/// This function transfroms the given TEventList into a TEntryList
/// The new TEntryList is owned by the TTree and gets deleted when the tree
/// is deleted. This TEntryList can be returned by GetEntryList() function.

void TTree::SetEventList(TEventList *evlist)
{
   fEventList = evlist;
   if (fEntryList){
      if (fEntryList->TestBit(kCanDelete)) {
         TEntryList *tmp = fEntryList;
         fEntryList = 0; // Avoid problem with RecursiveRemove.
         delete tmp;
      } else {
         fEntryList = 0;
      }
   }

   if (!evlist) {
      fEntryList = 0;
      fEventList = 0;
      return;
   }

   fEventList = evlist;
   char enlistname[100];
   snprintf(enlistname,100, "%s_%s", evlist->GetName(), "entrylist");
   fEntryList = new TEntryList(enlistname, evlist->GetTitle());
   fEntryList->SetDirectory(0); // We own this.
   Int_t nsel = evlist->GetN();
   fEntryList->SetTree(this);
   Long64_t entry;
   for (Int_t i=0; i<nsel; i++){
      entry = evlist->GetEntry(i);
      fEntryList->Enter(entry);
   }
   fEntryList->SetReapplyCut(evlist->GetReapplyCut());
   fEntryList->SetBit(kCanDelete, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of entries to estimate variable limits.
/// If n is -1, the estimate is set to be the current maximum
/// for the tree (i.e. GetEntries() + 1)
/// If n is less than -1, the behavior is undefined.

void TTree::SetEstimate(Long64_t n /* = 1000000 */)
{
   if (n == 0) {
      n = 10000;
   } else if (n < 0) {
      n = fEntries - n;
   }
   fEstimate = n;
   GetPlayer();
   if (fPlayer) {
      fPlayer->SetEstimate(n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Provide the end-user with the ability to enable/disable various experimental
/// IO features for this TTree.
///
/// Returns all the newly-set IO settings.

ROOT::TIOFeatures TTree::SetIOFeatures(const ROOT::TIOFeatures &features)
{
   // Purposely ignore all unsupported bits; TIOFeatures implementation already warned the user about the
   // error of their ways; this is just a safety check.
   UChar_t featuresRequested = features.GetFeatures() & static_cast<UChar_t>(TBasket::EIOBits::kSupported);

   UChar_t curFeatures = fIOFeatures.GetFeatures();
   UChar_t newFeatures = ~curFeatures & featuresRequested;
   curFeatures |= newFeatures;
   fIOFeatures.Set(curFeatures);

   ROOT::TIOFeatures newSettings(newFeatures);
   return newSettings;
}

////////////////////////////////////////////////////////////////////////////////
/// Set fFileNumber to number.
/// fFileNumber is used by TTree::Fill to set the file name
/// for a new file to be created when the current file exceeds fgTreeMaxSize.
///    (see TTree::ChangeFile)
/// if fFileNumber=10, the new file name will have a suffix "_11",
/// ie, fFileNumber is incremented before setting the file name

void TTree::SetFileNumber(Int_t number)
{
   if (fFileNumber < 0) {
      Warning("SetFileNumber", "file number must be positive. Set to 0");
      fFileNumber = 0;
      return;
   }
   fFileNumber = number;
}

////////////////////////////////////////////////////////////////////////////////
/// Set all the branches in this TTree to be in decomposed object mode
/// (also known as MakeClass mode).
///
/// For MakeClass mode 0, the TTree expects the address where the data is stored
/// to be set by either the user or the TTree to the address of a full object
/// through the top level branch.
/// For MakeClass mode 1, this address is expected to point to a numerical type
/// or C-style array (variable or not) of numerical type, representing the
/// primitive data members.
/// The function's primary purpose is to allow the user to access the data
/// directly with numerical type variable rather than having to have the original
/// set of classes (or a reproduction thereof).

void TTree::SetMakeClass(Int_t make)
{
   fMakeClass = make;

   Int_t nb = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i)  {
      TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
      branch->SetMakeClass(make);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum size in bytes of a Tree file (static function).
/// The default size is 100000000000LL, ie 100 Gigabytes.
///
/// In TTree::Fill, when the file has a size > fgMaxTreeSize,
/// the function closes the current file and starts writing into
/// a new file with a name of the style "file_1.root" if the original
/// requested file name was "file.root".

void TTree::SetMaxTreeSize(Long64_t maxsize)
{
   fgMaxTreeSize = maxsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Change the name of this tree.

void TTree::SetName(const char* name)
{
   if (gPad) {
      gPad->Modified();
   }
   // Trees are named objects in a THashList.
   // We must update hashlists if we change the name.
   TFile *file = 0;
   TTreeCache *pf = 0;
   if (fDirectory) {
      fDirectory->Remove(this);
      if ((file = GetCurrentFile())) {
         pf = GetReadCache(file);
         file->SetCacheRead(0,this,TFile::kDoNotDisconnect);
      }
   }
   // This changes our hash value.
   fName = name;
   if (fDirectory) {
      fDirectory->Append(this);
      if (pf) {
         file->SetCacheRead(pf,this,TFile::kDoNotDisconnect);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change the name and title of this tree.

void TTree::SetObject(const char* name, const char* title)
{
   if (gPad) {
      gPad->Modified();
   }

   //  Trees are named objects in a THashList.
   //  We must update hashlists if we change the name
   TFile *file = 0;
   TTreeCache *pf = 0;
   if (fDirectory) {
      fDirectory->Remove(this);
      if ((file = GetCurrentFile())) {
         pf = GetReadCache(file);
         file->SetCacheRead(0,this,TFile::kDoNotDisconnect);
      }
   }
   // This changes our hash value.
   fName = name;
   fTitle = title;
   if (fDirectory) {
      fDirectory->Append(this);
      if (pf) {
         file->SetCacheRead(pf,this,TFile::kDoNotDisconnect);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enable or disable parallel unzipping of Tree buffers.

void TTree::SetParallelUnzip(Bool_t opt, Float_t RelSize)
{
#ifdef R__USE_IMT
   if (GetTree() == 0) {
      LoadTree(GetReadEntry());
      if (!GetTree())
         return;
   }
   if (GetTree() != this) {
      GetTree()->SetParallelUnzip(opt, RelSize);
      return;
   }
   TFile* file = GetCurrentFile();
   if (!file)
      return;

   TTreeCache* pf = GetReadCache(file);
   if (pf && !( opt ^ (nullptr != dynamic_cast<TTreeCacheUnzip*>(pf)))) {
      // done with opt and type are in agreement.
      return;
   }
   delete pf;
   auto cacheSize = GetCacheAutoSize(kTRUE);
   if (opt) {
      auto unzip = new TTreeCacheUnzip(this, cacheSize);
      unzip->SetUnzipBufferSize( Long64_t(cacheSize * RelSize) );
   } else {
      pf = new TTreeCache(this, cacheSize);
   }
#else
   (void)opt;
   (void)RelSize;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Set perf stats

void TTree::SetPerfStats(TVirtualPerfStats *perf)
{
   fPerfStats = perf;
}

////////////////////////////////////////////////////////////////////////////////
/// The current TreeIndex is replaced by the new index.
/// Note that this function does not delete the previous index.
/// This gives the possibility to play with more than one index, e.g.,
/// ~~~ {.cpp}
///     TVirtualIndex* oldIndex = tree.GetTreeIndex();
///     tree.SetTreeIndex(newIndex);
///     tree.Draw();
///     tree.SetTreeIndex(oldIndex);
///     tree.Draw(); etc
/// ~~~

void TTree::SetTreeIndex(TVirtualIndex* index)
{
   if (fTreeIndex) {
      fTreeIndex->SetTree(0);
   }
   fTreeIndex = index;
}

////////////////////////////////////////////////////////////////////////////////
/// Set tree weight.
///
/// The weight is used by TTree::Draw to automatically weight each
/// selected entry in the resulting histogram.
///
/// For example the equivalent of:
/// ~~~ {.cpp}
///     T.Draw("x", "w")
/// ~~~
/// is:
/// ~~~ {.cpp}
///     T.SetWeight(w);
///     T.Draw("x");
/// ~~~
/// This function is redefined by TChain::SetWeight. In case of a
/// TChain, an option "global" may be specified to set the same weight
/// for all trees in the TChain instead of the default behaviour
/// using the weights of each tree in the chain (see TChain::SetWeight).

void TTree::SetWeight(Double_t w, Option_t*)
{
   fWeight = w;
}

////////////////////////////////////////////////////////////////////////////////
/// Print values of all active leaves for entry.
///
/// - if entry==-1, print current entry (default)
/// - if a leaf is an array, a maximum of lenmax elements is printed.

void TTree::Show(Long64_t entry, Int_t lenmax)
{
   if (entry != -1) {
      Int_t ret = LoadTree(entry);
      if (ret == -2) {
         Error("Show()", "Cannot read entry %lld (entry does not exist)", entry);
         return;
      } else if (ret == -1) {
         Error("Show()", "Cannot read entry %lld (I/O error)", entry);
         return;
      }
      ret = GetEntry(entry);
      if (ret == -1) {
         Error("Show()", "Cannot read entry %lld (I/O error)", entry);
         return;
      } else if (ret == 0) {
         Error("Show()", "Cannot read entry %lld (no data read)", entry);
         return;
      }
   }
   printf("======> EVENT:%lld\n", fReadEntry);
   TObjArray* leaves  = GetListOfLeaves();
   Int_t nleaves = leaves->GetEntriesFast();
   Int_t ltype;
   for (Int_t i = 0; i < nleaves; i++) {
      TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(i);
      TBranch* branch = leaf->GetBranch();
      if (branch->TestBit(kDoNotProcess)) {
         continue;
      }
      Int_t len = leaf->GetLen();
      if (len <= 0) {
         continue;
      }
      len = TMath::Min(len, lenmax);
      if (leaf->IsA() == TLeafElement::Class()) {
         leaf->PrintValue(lenmax);
         continue;
      }
      if (branch->GetListOfBranches()->GetEntriesFast() > 0) {
         continue;
      }
      ltype = 10;
      if (leaf->IsA() == TLeafF::Class()) {
         ltype = 5;
      }
      if (leaf->IsA() == TLeafD::Class()) {
         ltype = 5;
      }
      if (leaf->IsA() == TLeafC::Class()) {
         len = 1;
         ltype = 5;
      };
      printf(" %-15s = ", leaf->GetName());
      for (Int_t l = 0; l < len; l++) {
         leaf->PrintValue(l);
         if (l == (len - 1)) {
            printf("\n");
            continue;
         }
         printf(", ");
         if ((l % ltype) == 0) {
            printf("\n                  ");
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Start the TTreeViewer on this tree.
///
/// - ww is the width of the canvas in pixels
/// - wh is the height of the canvas in pixels

void TTree::StartViewer()
{
   GetPlayer();
   if (fPlayer) {
      fPlayer->StartViewer(600, 400);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stop the cache learning phase
///
/// Returns:
/// - 0 learning phase stopped or not active
/// - -1 on error

Int_t TTree::StopCacheLearningPhase()
{
   if (!GetTree()) {
      if (LoadTree(0)<0) {
         Error("StopCacheLearningPhase","Could not load a tree");
         return -1;
      }
   }
   if (GetTree()) {
      if (GetTree() != this) {
         return GetTree()->StopCacheLearningPhase();
      }
   } else {
      Error("StopCacheLearningPhase", "No tree is available. Could not stop cache learning phase");
      return -1;
   }

   TFile *f = GetCurrentFile();
   if (!f) {
      Error("StopCacheLearningPhase", "No file is available. Could not stop cache learning phase");
      return -1;
   }
   TTreeCache *tc = GetReadCache(f,kTRUE);
   if (!tc) {
      Error("StopCacheLearningPhase", "No cache is available. Could not stop learning phase");
      return -1;
   }
   tc->StopLearningPhase();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the fTree member for all branches and sub branches.

static void TBranch__SetTree(TTree *tree, TObjArray &branches)
{
   Int_t nb = branches.GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i) {
      TBranch* br = (TBranch*) branches.UncheckedAt(i);
      br->SetTree(tree);

      Int_t nBaskets = br->GetListOfBaskets()->GetEntries();
      Int_t writeBasket = br->GetWriteBasket();
      for (Int_t j=writeBasket,n=0;j>=0 && n<nBaskets;--j) {
         TBasket *bk = (TBasket*)br->GetListOfBaskets()->UncheckedAt(j);
         if (bk) {
            tree->IncrementTotalBuffers(bk->GetBufferSize());
            ++n;
         }
      }

      TBranch__SetTree(tree,*br->GetListOfBranches());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the fTree member for all friend elements.

void TFriendElement__SetTree(TTree *tree, TList *frlist)
{
   if (frlist) {
      TObjLink *lnk = frlist->FirstLink();
      while (lnk) {
         TFriendElement *elem = (TFriendElement*)lnk->GetObject();
         elem->fParentTree = tree;
         lnk = lnk->Next();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TTree::Streamer(TBuffer& b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      if (fDirectory) {
         fDirectory->Remove(this);
         //delete the file cache if it points to this Tree
         TFile *file = fDirectory->GetFile();
         MoveReadCache(file,0);
      }
      fDirectory = 0;
      fCacheDoAutoInit = kTRUE;
      fCacheUserSet = kFALSE;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 4) {
         b.ReadClassBuffer(TTree::Class(), this, R__v, R__s, R__c);

         fBranches.SetOwner(kTRUE); // True needed only for R__v < 19 and most R__v == 19

         if (fBranchRef) fBranchRef->SetTree(this);
         TBranch__SetTree(this,fBranches);
         TFriendElement__SetTree(this,fFriends);

         if (fTreeIndex) {
            fTreeIndex->SetTree(this);
         }
         if (fIndex.fN) {
            Warning("Streamer", "Old style index in this tree is deleted. Rebuild the index via TTree::BuildIndex");
            fIndex.Set(0);
            fIndexValues.Set(0);
         }
         if (fEstimate <= 10000) {
            fEstimate = 1000000;
         }

         if (fNClusterRange) {
            // The I/O allocated just enough memory to hold the
            // current set of ranges.
            fMaxClusterRange = fNClusterRange;
         }
         if (GetCacheAutoSize() != 0) {
            // a cache will be automatically created.
            // No need for TTreePlayer::Process to enable the cache
            fCacheSize = 0;
         } else if (fAutoFlush < 0) {
            // If there is no autoflush set, let's keep the cache completely
            // disable by default for now.
            fCacheSize = fAutoFlush;
         } else if (fAutoFlush != 0) {
            // Estimate the cluster size.
            // This will allow TTree::Process to enable the cache.
            Long64_t zipBytes = GetZipBytes();
            Long64_t totBytes = GetTotBytes();
            if (zipBytes != 0) {
               fCacheSize =  fAutoFlush*(zipBytes/fEntries);
            } else if (totBytes != 0) {
               fCacheSize =  fAutoFlush*(totBytes/fEntries);
            } else {
               fCacheSize = 30000000;
            }
            if (fCacheSize >= (INT_MAX / 4)) {
               fCacheSize = INT_MAX / 4;
            } else if (fCacheSize == 0) {
               fCacheSize = 30000000;
            }
         } else {
            fCacheSize = 0;
         }
         ResetBit(kMustCleanup);
         return;
      }
      //====process old versions before automatic schema evolution
      Stat_t djunk;
      Int_t ijunk;
      TNamed::Streamer(b);
      TAttLine::Streamer(b);
      TAttFill::Streamer(b);
      TAttMarker::Streamer(b);
      b >> fScanField;
      b >> ijunk; fMaxEntryLoop   = (Long64_t)ijunk;
      b >> ijunk; fMaxVirtualSize = (Long64_t)ijunk;
      b >> djunk; fEntries  = (Long64_t)djunk;
      b >> djunk; fTotBytes = (Long64_t)djunk;
      b >> djunk; fZipBytes = (Long64_t)djunk;
      b >> ijunk; fAutoSave = (Long64_t)ijunk;
      b >> ijunk; fEstimate = (Long64_t)ijunk;
      if (fEstimate <= 10000) fEstimate = 1000000;
      fBranches.Streamer(b);
      if (fBranchRef) fBranchRef->SetTree(this);
      TBranch__SetTree(this,fBranches);
      fLeaves.Streamer(b);
      fSavedBytes = fTotBytes;
      if (R__v > 1) fIndexValues.Streamer(b);
      if (R__v > 2) fIndex.Streamer(b);
      if (R__v > 3) {
         TList OldInfoList;
         OldInfoList.Streamer(b);
         OldInfoList.Delete();
      }
      fNClusterRange = 0;
      fDefaultEntryOffsetLen = 1000;
      ResetBit(kMustCleanup);
      b.CheckByteCount(R__s, R__c, TTree::IsA());
      //====end of old versions
   } else {
      if (fBranchRef) {
         fBranchRef->Clear();
      }
      TRefTable *table  = TRefTable::GetRefTable();
      if (table) TRefTable::SetRefTable(0);

      b.WriteClassBuffer(TTree::Class(), this);

      if (table) TRefTable::SetRefTable(table);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Unbinned fit of one or more variable(s) from a tree.
///
/// funcname is a TF1 function.
///
/// See TTree::Draw for explanations of the other parameters.
///
/// Fit the variable varexp using the function funcname using the
/// selection cuts given by selection.
///
/// The list of fit options is given in parameter option.
///
/// - option = "Q" Quiet mode (minimum printing)
/// - option = "V" Verbose mode (default is between Q and V)
/// - option = "E" Perform better Errors estimation using Minos technique
/// - option = "M" More. Improve fit results
///
/// You can specify boundary limits for some or all parameters via
/// ~~~ {.cpp}
///     func->SetParLimits(p_number, parmin, parmax);
/// ~~~
/// if parmin>=parmax, the parameter is fixed
///
/// Note that you are not forced to fix the limits for all parameters.
/// For example, if you fit a function with 6 parameters, you can do:
/// ~~~ {.cpp}
///     func->SetParameters(0,3.1,1.e-6,0.1,-8,100);
///     func->SetParLimits(4,-10,-4);
///     func->SetParLimits(5, 1,1);
/// ~~~
/// With this setup:
///
/// - Parameters 0->3 can vary freely
/// - Parameter 4 has boundaries [-10,-4] with initial value -8
/// - Parameter 5 is fixed to 100.
///
/// For the fit to be meaningful, the function must be self-normalized.
///
/// i.e. It must have the same integral regardless of the parameter
/// settings.  Otherwise the fit will effectively just maximize the
/// area.
///
/// It is mandatory to have a normalization variable
/// which is fixed for the fit.  e.g.
/// ~~~ {.cpp}
///     TF1* f1 = new TF1("f1", "gaus(0)/sqrt(2*3.14159)/[2]", 0, 5);
///     f1->SetParameters(1, 3.1, 0.01);
///     f1->SetParLimits(0, 1, 1); // fix the normalization parameter to 1
///     data->UnbinnedFit("f1", "jpsimass", "jpsipt>3.0");
/// ~~~
/// 1, 2 and 3 Dimensional fits are supported. See also TTree::Fit
///
/// Return status:
///
/// - The function return the status of the fit in the following form
///   fitResult = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult
/// - The fitResult is 0 is the fit is OK.
/// - The fitResult is negative in case of an error not connected with the fit.
/// - The number of entries used in the fit can be obtained via mytree.GetSelectedRows();
/// - If the number of selected entries is null the function returns -1

Int_t TTree::UnbinnedFit(const char* funcname, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   GetPlayer();
   if (fPlayer) {
      return fPlayer->UnbinnedFit(funcname, varexp, selection, option, nentries, firstentry);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace current attributes by current style.

void TTree::UseCurrentStyle()
{
   if (gStyle->IsReading()) {
      SetFillColor(gStyle->GetHistFillColor());
      SetFillStyle(gStyle->GetHistFillStyle());
      SetLineColor(gStyle->GetHistLineColor());
      SetLineStyle(gStyle->GetHistLineStyle());
      SetLineWidth(gStyle->GetHistLineWidth());
      SetMarkerColor(gStyle->GetMarkerColor());
      SetMarkerStyle(gStyle->GetMarkerStyle());
      SetMarkerSize(gStyle->GetMarkerSize());
   } else {
      gStyle->SetHistFillColor(GetFillColor());
      gStyle->SetHistFillStyle(GetFillStyle());
      gStyle->SetHistLineColor(GetLineColor());
      gStyle->SetHistLineStyle(GetLineStyle());
      gStyle->SetHistLineWidth(GetLineWidth());
      gStyle->SetMarkerColor(GetMarkerColor());
      gStyle->SetMarkerStyle(GetMarkerStyle());
      gStyle->SetMarkerSize(GetMarkerSize());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write this object to the current directory. For more see TObject::Write
/// If option & kFlushBasket, call FlushBasket before writing the tree.

Int_t TTree::Write(const char *name, Int_t option, Int_t bufsize) const
{
   FlushBasketsImpl();
   return TObject::Write(name, option, bufsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Write this object to the current directory. For more see TObject::Write
/// If option & kFlushBasket, call FlushBasket before writing the tree.

Int_t TTree::Write(const char *name, Int_t option, Int_t bufsize)
{
   return ((const TTree*)this)->Write(name, option, bufsize);
}

////////////////////////////////////////////////////////////////////////////////
/// \class TTreeFriendLeafIter
///
/// Iterator on all the leaves in a TTree and its friend

ClassImp(TTreeFriendLeafIter);

////////////////////////////////////////////////////////////////////////////////
/// Create a new iterator. By default the iteration direction
/// is kIterForward. To go backward use kIterBackward.

TTreeFriendLeafIter::TTreeFriendLeafIter(const TTree* tree, Bool_t dir)
: fTree(const_cast<TTree*>(tree))
, fLeafIter(0)
, fTreeIter(0)
, fDirection(dir)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.  Does NOT copy the 'cursor' location!

TTreeFriendLeafIter::TTreeFriendLeafIter(const TTreeFriendLeafIter& iter)
: TIterator(iter)
, fTree(iter.fTree)
, fLeafIter(0)
, fTreeIter(0)
, fDirection(iter.fDirection)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Overridden assignment operator. Does NOT copy the 'cursor' location!

TIterator& TTreeFriendLeafIter::operator=(const TIterator& rhs)
{
   if (this != &rhs && rhs.IsA() == TTreeFriendLeafIter::Class()) {
      const TTreeFriendLeafIter &rhs1 = (const TTreeFriendLeafIter &)rhs;
      fDirection = rhs1.fDirection;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Overridden assignment operator.  Does NOT copy the 'cursor' location!

TTreeFriendLeafIter& TTreeFriendLeafIter::operator=(const TTreeFriendLeafIter& rhs)
{
   if (this != &rhs) {
      fDirection = rhs.fDirection;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Go the next friend element

TObject* TTreeFriendLeafIter::Next()
{
   if (!fTree) return 0;

   TObject * next;
   TTree * nextTree;

   if (!fLeafIter) {
      TObjArray *list = fTree->GetListOfLeaves();
      if (!list) return 0; // Can happen with an empty chain.
      fLeafIter =  list->MakeIterator(fDirection);
      if (!fLeafIter) return 0;
   }

   next = fLeafIter->Next();
   if (!next) {
      if (!fTreeIter) {
         TCollection * list = fTree->GetListOfFriends();
         if (!list) return next;
         fTreeIter = list->MakeIterator(fDirection);
         if (!fTreeIter) return 0;
      }
      TFriendElement * nextFriend = (TFriendElement*) fTreeIter->Next();
      ///nextTree = (TTree*)fTreeIter->Next();
      if (nextFriend) {
         nextTree = const_cast<TTree*>(nextFriend->GetTree());
         if (!nextTree) return Next();
         SafeDelete(fLeafIter);
         fLeafIter = nextTree->GetListOfLeaves()->MakeIterator(fDirection);
         if (!fLeafIter) return 0;
         next = fLeafIter->Next();
      }
   }
   return next;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the object option stored in the list.

Option_t* TTreeFriendLeafIter::GetOption() const
{
   if (fLeafIter) return fLeafIter->GetOption();
   return "";
}
