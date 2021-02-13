// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTree
#define ROOT_TTree

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTree                                                                //
//                                                                      //
// A TTree object is a list of TBranch.                                 //
//   To Create a TTree object one must:                                 //
//    - Create the TTree header via the TTree constructor               //
//    - Call the TBranch constructor for every branch.                  //
//                                                                      //
//   To Fill this object, use member function Fill with no parameters.  //
//     The Fill function loops on all defined TBranch.                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Compression.h"
#include "ROOT/TIOFeatures.hxx"
#include "TArrayD.h"
#include "TArrayI.h"
#include "TAttFill.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TClass.h"
#include "TDataType.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TVirtualTreePlayer.h"

#ifdef R__LESS_INCLUDES
class TBranch;
class TList;
#else
#include "TBranch.h"
// #include "TBuffer.h"
#include "TList.h"
#endif

#include <array>
#include <atomic>
#include <vector>
#include <utility>

class TBuffer;
class TBrowser;
class TFile;
class TLeaf;
class TH1;
class TTreeFormula;
class TPolyMarker;
class TEventList;
class TEntryList;
class TSQLResult;
class TSelector;
class TPrincipal;
class TFriendElement;
class TCut;
class TVirtualIndex;
class TBranchRef;
class TBasket;
class TStreamerInfo;
class TTreeCache;
class TTreeCloner;
class TFileMergeInfo;
class TVirtualPerfStats;

class TTree : public TNamed, public TAttLine, public TAttFill, public TAttMarker {

   using TIOFeatures = ROOT::TIOFeatures;

protected:
   Long64_t       fEntries;               ///<  Number of entries
// NOTE: cannot use std::atomic for these counters as it cannot be serialized.
   Long64_t       fTotBytes;              ///<  Total number of bytes in all branches before compression
   Long64_t       fZipBytes;              ///<  Total number of bytes in all branches after compression
   Long64_t       fSavedBytes;            ///<  Number of autosaved bytes
   Long64_t       fFlushedBytes;          ///<  Number of auto-flushed bytes
   Double_t       fWeight;                ///<  Tree weight (see TTree::SetWeight)
   Int_t          fTimerInterval;         ///<  Timer interval in milliseconds
   Int_t          fScanField;             ///<  Number of runs before prompting in Scan
   Int_t          fUpdate;                ///<  Update frequency for EntryLoop
   Int_t          fDefaultEntryOffsetLen; ///<  Initial Length of fEntryOffset table in the basket buffers
   Int_t          fNClusterRange;         ///<  Number of Cluster range in addition to the one defined by 'AutoFlush'
   Int_t          fMaxClusterRange;       ///<! Memory allocated for the cluster range.
   Long64_t       fMaxEntries;            ///<  Maximum number of entries in case of circular buffers
   Long64_t       fMaxEntryLoop;          ///<  Maximum number of entries to process
   Long64_t       fMaxVirtualSize;        ///<  Maximum total size of buffers kept in memory
   Long64_t       fAutoSave;              ///<  Autosave tree when fAutoSave entries written or -fAutoSave (compressed) bytes produced
   Long64_t       fAutoFlush;             ///<  Auto-flush tree when fAutoFlush entries written or -fAutoFlush (compressed) bytes produced
   Long64_t       fEstimate;              ///<  Number of entries to estimate histogram limits
   Long64_t      *fClusterRangeEnd;       ///<[fNClusterRange] Last entry of a cluster range.
   Long64_t      *fClusterSize;           ///<[fNClusterRange] Number of entries in each cluster for a given range.
   Long64_t       fCacheSize;             ///<! Maximum size of file buffers
   Long64_t       fChainOffset;           ///<! Offset of 1st entry of this Tree in a TChain
   Long64_t       fReadEntry;             ///<! Number of the entry being processed
   std::atomic<Long64_t> fTotalBuffers;   ///<! Total number of bytes in branch buffers
   Int_t          fPacketSize;            ///<! Number of entries in one packet for parallel root
   Int_t          fNfill;                 ///<! Local for EntryLoop
   Int_t          fDebug;                 ///<! Debug level
   Long64_t       fDebugMin;              ///<! First entry number to debug
   Long64_t       fDebugMax;              ///<! Last entry number to debug
   TIOFeatures    fIOFeatures{0};         ///<  IO features to define for newly-written baskets and branches.
   Int_t          fMakeClass;             ///<! not zero when processing code generated by MakeClass
   Int_t          fFileNumber;            ///<! current file number (if file extensions)
   TObject       *fNotify;                ///<! Object to be notified when loading a Tree
   TDirectory    *fDirectory;             ///<! Pointer to directory holding this tree
   TObjArray      fBranches;              ///<  List of Branches
   TObjArray      fLeaves;                ///<  Direct pointers to individual branch leaves
   TList         *fAliases;               ///<  List of aliases for expressions based on the tree branches.
   TEventList    *fEventList;             ///<! Pointer to event selection list (if one)
   TEntryList    *fEntryList;             ///<! Pointer to event selection list (if one)
   TArrayD        fIndexValues;           ///<  Sorted index values
   TArrayI        fIndex;                 ///<  Index of sorted values
   TVirtualIndex *fTreeIndex;             ///<  Pointer to the tree Index (if any)
   TList         *fFriends;               ///<  pointer to list of friend elements
   TList         *fExternalFriends;       ///<! List of TFriendsElement pointing to us and need to be notified of LoadTree.  Content not owned.
   TVirtualPerfStats *fPerfStats;         ///<! pointer to the current perf stats object
   TList         *fUserInfo;              ///<  pointer to a list of user objects associated to this Tree
   TVirtualTreePlayer *fPlayer;           ///<! Pointer to current Tree player
   TList         *fClones;                ///<! List of cloned trees which share our addresses
   TBranchRef    *fBranchRef;             ///<  Branch supporting the TRefTable (if any)
   UInt_t         fFriendLockStatus;      ///<! Record which method is locking the friend recursion
   TBuffer       *fTransientBuffer;       ///<! Pointer to the current transient buffer.
   Bool_t         fCacheDoAutoInit;       ///<! true if cache auto creation or resize check is needed
   Bool_t         fCacheDoClusterPrefetch;///<! true if cache is prefetching whole clusters
   Bool_t         fCacheUserSet;          ///<! true if the cache setting was explicitly given by user
   Bool_t         fIMTEnabled;            ///<! true if implicit multi-threading is enabled for this tree
   UInt_t         fNEntriesSinceSorting;  ///<! Number of entries processed since the last re-sorting of branches
   std::vector<std::pair<Long64_t,TBranch*>> fSortedBranches; ///<! Branches to be processed in parallel when IMT is on, sorted by average task time
   std::vector<TBranch*> fSeqBranches;    ///<! Branches to be processed sequentially when IMT is on
   Float_t fTargetMemoryRatio{1.1f};      ///<! Ratio for memory usage in uncompressed buffers versus actual occupancy.  1.0
                                           /// indicates basket should be resized to exact memory usage, but causes significant
/// memory churn.
#ifdef R__TRACK_BASKET_ALLOC_TIME
   mutable std::atomic<ULong64_t> fAllocationTime{0}; ///<! Time spent reallocating basket memory buffers, in microseconds.
#endif
   mutable std::atomic<UInt_t> fAllocationCount{0};   ///<! Number of reallocations basket memory buffers.

   static Int_t     fgBranchStyle;        ///<  Old/New branch style
   static Long64_t  fgMaxTreeSize;        ///<  Maximum size of a file containing a Tree

private:
   // For simplicity, although fIMTFlush is always disabled in non-IMT builds, we don't #ifdef it out.
   mutable Bool_t fIMTFlush{false};               ///<! True if we are doing a multithreaded flush.
   mutable std::atomic<Long64_t> fIMTTotBytes;    ///<! Total bytes for the IMT flush baskets
   mutable std::atomic<Long64_t> fIMTZipBytes;    ///<! Zip bytes for the IMT flush baskets.

   void             InitializeBranchLists(bool checkLeafCount);
   void             SortBranchesByTime();
   Int_t            FlushBasketsImpl() const;
   void             MarkEventCluster();

protected:
   virtual void     KeepCircular();
   virtual TBranch *BranchImp(const char* branchname, const char* classname, TClass* ptrClass, void* addobj, Int_t bufsize, Int_t splitlevel);
   virtual TBranch *BranchImp(const char* branchname, TClass* ptrClass, void* addobj, Int_t bufsize, Int_t splitlevel);
   virtual TBranch *BranchImpRef(const char* branchname, const char* classname, TClass* ptrClass, void* addobj, Int_t bufsize, Int_t splitlevel);
   virtual TBranch *BranchImpRef(const char* branchname, TClass* ptrClass, EDataType datatype, void* addobj, Int_t bufsize, Int_t splitlevel);
   virtual TBranch *BranchImpArr(const char* branchname, EDataType datatype, std::size_t N, void* addobj, Int_t bufsize, Int_t splitlevel);
   virtual Int_t    CheckBranchAddressType(TBranch* branch, TClass* ptrClass, EDataType datatype, Bool_t ptr);
   virtual TBranch *BronchExec(const char* name, const char* classname, void* addobj, Bool_t isptrptr, Int_t bufsize, Int_t splitlevel);
   friend  TBranch *TTreeBranchImpRef(TTree *tree, const char* branchname, TClass* ptrClass, EDataType datatype, void* addobj, Int_t bufsize, Int_t splitlevel);
   Int_t    SetBranchAddressImp(TBranch *branch, void* addr, TBranch** ptr);
   virtual TLeaf   *GetLeafImpl(const char* branchname, const char* leafname);

   Long64_t         GetCacheAutoSize(Bool_t withDefault = kFALSE) const;
   char             GetNewlineValue(std::istream &inputStream);
   void             ImportClusterRanges(TTree *fromtree);
   void             MoveReadCache(TFile *src, TDirectory *dir);
   Int_t            SetCacheSizeAux(Bool_t autocache = kTRUE, Long64_t cacheSize = 0);

   class TFriendLock {
      // Helper class to prevent infinite recursion in the
      // usage of TTree Friends. Implemented in TTree.cxx.
      TTree  *fTree;      // Pointer to the locked tree
      UInt_t  fMethodBit; // BIT for the locked method
      Bool_t  fPrevious;  // Previous value of the BIT.

   protected:
      TFriendLock(const TFriendLock&);
      TFriendLock& operator=(const TFriendLock&);

   public:
      TFriendLock(TTree* tree, UInt_t methodbit);
      ~TFriendLock();
   };
   friend class TFriendLock;
   // So that the index class can use TFriendLock:
   friend class TTreeIndex;
   friend class TChainIndex;
   // So that the TTreeCloner can access the protected interfaces
   friend class TTreeCloner;

   // use to update fFriendLockStatus
   enum ELockStatusBits {
      kFindBranch        = BIT(0),
      kFindLeaf          = BIT(1),
      kGetAlias          = BIT(2),
      kGetBranch         = BIT(3),
      kGetEntry          = BIT(4),
      kGetEntryWithIndex = BIT(5),
      kGetFriend         = BIT(6),
      kGetFriendAlias    = BIT(7),
      kGetLeaf           = BIT(8),
      kLoadTree          = BIT(9),
      kPrint             = BIT(10),
      kRemoveFriend      = BIT(11),
      kSetBranchStatus   = BIT(12)
   };

public:
   // Used as the max value for any TTree range operation.
   static constexpr Long64_t kMaxEntries = TVirtualTreePlayer::kMaxEntries;

   // SetBranchAddress return values
   enum ESetBranchAddressStatus {
      kMissingBranch = -5,
      kInternalError = -4,
      kMissingCompiledCollectionProxy = -3,
      kMismatch = -2,
      kClassMismatch = -1,
      kMatch = 0,
      kMatchConversion = 1,
      kMatchConversionCollection = 2,
      kMakeClass = 3,
      kVoidPtr = 4,
      kNoCheck = 5,
      kNeedEnableDecomposedObj = BIT(29),   // DecomposedObj is the newer name of MakeClass mode
      kNeedDisableDecomposedObj = BIT(30),
      kDecomposedObjMask = kNeedEnableDecomposedObj | kNeedDisableDecomposedObj
   };

   // TTree status bits
   enum EStatusBits {
      kForceRead = BIT(11),
      kCircular = BIT(12),
      /// If set, the branch's buffers will grow until an event cluster boundary is hit,
      /// guaranteeing a basket per cluster.  This mode does not provide any guarantee on the
      /// memory bounds in the case of extremely large events.
      kOnlyFlushAtCluster = BIT(14),
      /// If set, signals that this TTree is the output of the processing of another TTree, and
      /// the entries are reshuffled w.r.t. to the original TTree. As a safety measure, a TTree
      /// with this bit set cannot add friends nor can be added as a friend. If you know what
      /// you are doing, you can manually unset this bit with `ResetBit(EStatusBits::kEntriesReshuffled)`.
      kEntriesReshuffled = BIT(19) // bits 15-18 are used by TChain
   };

   // Split level modifier
   enum {
      kSplitCollectionOfPointers = 100
   };

   class TClusterIterator
   {
   private:
      TTree    *fTree;         // TTree upon which we are iterating.
      Int_t    fClusterRange;  // Which cluster range are we looking at.
      Long64_t fStartEntry;    // Where does the cluster start.
      Long64_t fNextEntry;     // Where does the cluster end (exclusive).
      Long64_t fEstimatedSize; // If positive, the calculated estimated tree size.

      Long64_t GetEstimatedClusterSize();

   protected:
      friend class TTree;
      TClusterIterator(TTree *tree, Long64_t firstEntry);

   public:
      // Intentionally used the default copy constructor and default destructor
      // as the TClusterIterator does not own the TTree.
      //  TClusterIterator(const TClusterIterator&);
      // ~TClusterIterator();

      // No public constructors, the iterator must be
      // created via TTree::GetClusterIterator

      // Move on to the next cluster and return the starting entry
      // of this next cluster
      Long64_t Next();

      // Move on to the previous cluster and return the starting entry
      // of this previous cluster
      Long64_t Previous();

      // Return the start entry of the current cluster.
      Long64_t GetStartEntry() {
         return fStartEntry;
      }

      // Return the first entry of the next cluster.
      Long64_t GetNextEntry() {
         return fNextEntry;
      }

      Long64_t operator()() { return Next(); }
   };

   TTree();
   TTree(const char* name, const char* title, Int_t splitlevel = 99, TDirectory* dir = gDirectory);
   virtual ~TTree();

   TTree(const TTree& tt) = delete;
   TTree& operator=(const TTree& tt) = delete;

   virtual Int_t           AddBranchToCache(const char *bname, Bool_t subbranches = kFALSE);
   virtual Int_t           AddBranchToCache(TBranch *branch,   Bool_t subbranches = kFALSE);
   virtual Int_t           DropBranchFromCache(const char *bname, Bool_t subbranches = kFALSE);
   virtual Int_t           DropBranchFromCache(TBranch *branch,   Bool_t subbranches = kFALSE);
   void                    AddClone(TTree*);
   virtual TFriendElement *AddFriend(const char* treename, const char* filename = "");
   virtual TFriendElement *AddFriend(const char* treename, TFile* file);
   virtual TFriendElement *AddFriend(TTree* tree, const char* alias = "", Bool_t warn = kFALSE);
   // As the TBasket invokes Add{Tot,Zip}Bytes on its parent tree, we must do these updates in a thread-safe
   // manner only when we are flushing multiple baskets in parallel.
   virtual void            AddTotBytes(Int_t tot) { if (fIMTFlush) { fIMTTotBytes += tot; } else { fTotBytes += tot; } }
   virtual void            AddZipBytes(Int_t zip) { if (fIMTFlush) { fIMTZipBytes += zip; } else { fZipBytes += zip; } }
// NOTE: these counters aren't thread safe like the ones above.
#ifdef R__TRACK_BASKET_ALLOC_TIME
   void AddAllocationTime(ULong64_t time) { fAllocationTime += time; }
#endif
   void AddAllocationCount(UInt_t count) { fAllocationCount += count; }
   virtual Long64_t        AutoSave(Option_t* option = "");

   /// Add a new branch, and infer the data type from the type of `obj` being passed.
   ///
   /// \note This and the next overload should cover most cases for creating a branch. Try to use these two whenever
   /// possible, unless e.g. type conversions are needed.
   ///
   /// \param[in] name Name of the branch to be created.
   /// \param[in] obj Address of the object to be added. Make sure to pass a pointer to the actual type/class that
   /// should be stored in the tree (no pointers to base classes). When calling Fill(), the current value of the type/object will be saved.
   /// \param[in] bufsize The buffer size in bytes for this branch. When the buffer is full, it is compressed and written to disc.
   /// The default value of 32000 bytes and should be ok for most simple types. Larger buffers (e.g. 256000) if your Tree is not split and each entry is large (Megabytes).
   /// A small value for bufsize is beneficial if entries in the Tree are accessed randomly and the Tree is in split mode.
   /// \param[in] splitlevel If T is a class or struct and splitlevel > 0, the members of the object are serialised as separate branches.
   /// \return Pointer to the TBranch that was created. The branch is owned by the tree.
   template <class T> TBranch *Branch(const char* name, T* obj, Int_t bufsize = 32000, Int_t splitlevel = 99)
   {
      return BranchImpRef(name, TClass::GetClass<T>(), TDataType::GetType(typeid(T)), obj, bufsize, splitlevel);
   }

   /// Add a new branch, and infer the data type from the array `addobj` being passed.
   ///
   /// \note This and the previous overload should cover most cases for creating a branch. Try to use these two whenever
   /// possible, unless e.g. type conversions are needed.
   ///
   /// \param[in] name Name of the branch to be created.
   /// \param[in] addobj Array of the objects to be added. When calling Fill(), the current value of the type/object will be saved.
   /// \param[in] bufsize he buffer size in bytes for this branch. When the buffer is full, it is compressed and written to disc.
   /// The default value of 32000 bytes and should be ok for most simple types. Larger buffers (e.g. 256000) if your Tree is not split and each entry is large (Megabytes).
   /// A small value for bufsize is beneficial if entries in the Tree are accessed randomly and the Tree is in split mode.
   /// \param[in] splitlevel If T is a class or struct and splitlevel > 0, the members of the object are serialised as separate branches.
   /// \return Pointer to the TBranch that was created. The branch is owned by the tree.
   template <class T> TBranch *Branch(const char* name, T** addobj, Int_t bufsize = 32000, Int_t splitlevel = 99)
   {
      return BranchImp(name, TClass::GetClass<T>(), addobj, bufsize, splitlevel);
   }

   virtual Int_t           Branch(TCollection* list, Int_t bufsize = 32000, Int_t splitlevel = 99, const char* name = "");
   virtual Int_t           Branch(TList* list, Int_t bufsize = 32000, Int_t splitlevel = 99);
   virtual Int_t           Branch(const char* folder, Int_t bufsize = 32000, Int_t splitlevel = 99);
   virtual TBranch        *Branch(const char* name, void* address, const char* leaflist, Int_t bufsize = 32000);
           TBranch        *Branch(const char* name, char* address, const char* leaflist, Int_t bufsize = 32000)
   {
      // Overload to avoid confusion between this signature and the template instance.
      return Branch(name,(void*)address,leaflist,bufsize);
   }
   TBranch        *Branch(const char* name, Long_t address, const char* leaflist, Int_t bufsize = 32000)
   {
      // Overload to avoid confusion between this signature and the template instance.
      return Branch(name,(void*)address,leaflist,bufsize);
   }
   TBranch        *Branch(const char* name, int address, const char* leaflist, Int_t bufsize = 32000)
   {
      // Overload to avoid confusion between this signature and the template instance.
      return Branch(name,(void*)(Long_t)address,leaflist,bufsize);
   }
   virtual TBranch        *Branch(const char* name, const char* classname, void* addobj, Int_t bufsize = 32000, Int_t splitlevel = 99);
   template <class T> TBranch *Branch(const char* name, const char* classname, T* obj, Int_t bufsize = 32000, Int_t splitlevel = 99)
   {
      // See BranchImpRed for details. Here we __ignore
      return BranchImpRef(name, classname, TClass::GetClass<T>(), obj, bufsize, splitlevel);
   }
   template <class T> TBranch *Branch(const char* name, const char* classname, T** addobj, Int_t bufsize = 32000, Int_t splitlevel = 99)
   {
      // See BranchImp for details
      return BranchImp(name, classname, TClass::GetClass<T>(), addobj, bufsize, splitlevel);
   }
   template <typename T, std::size_t N> TBranch *Branch(const char* name, std::array<T, N> *obj, Int_t bufsize = 32000, Int_t splitlevel = 99)
   {
      TClass *cl = TClass::GetClass<T>();
      if (cl) {
         TClass *arrCl = TClass::GetClass<std::array<T, N>>();
         Error("Branch","std::array of objects not yet supported as top level branch object (the class is %s)",
               arrCl ? arrCl->GetName() : cl->GetName());
         return nullptr;
      }
      return BranchImpArr(name, TDataType::GetType(typeid(T)), N, obj, bufsize, splitlevel);
   }
   virtual TBranch        *Bronch(const char* name, const char* classname, void* addobj, Int_t bufsize = 32000, Int_t splitlevel = 99);
   virtual TBranch        *BranchOld(const char* name, const char* classname, void* addobj, Int_t bufsize = 32000, Int_t splitlevel = 1);
   virtual TBranch        *BranchRef();
   virtual void            Browse(TBrowser*);
   virtual Int_t           BuildIndex(const char* majorname, const char* minorname = "0");
   TStreamerInfo          *BuildStreamerInfo(TClass* cl, void* pointer = 0, Bool_t canOptimize = kTRUE);
   virtual TFile          *ChangeFile(TFile* file);
   virtual TTree          *CloneTree(Long64_t nentries = -1, Option_t* option = "");
   virtual void            CopyAddresses(TTree*,Bool_t undo = kFALSE);
   virtual Long64_t        CopyEntries(TTree* tree, Long64_t nentries = -1, Option_t *option = "");
   virtual TTree          *CopyTree(const char* selection, Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   virtual TBasket        *CreateBasket(TBranch*);
   virtual void            DirectoryAutoAdd(TDirectory *);
   Int_t                   Debug() const { return fDebug; }
   virtual void            Delete(Option_t* option = ""); // *MENU*
   virtual void            Draw(Option_t* opt) { Draw(opt, "", "", kMaxEntries, 0); }
   virtual Long64_t        Draw(const char* varexp, const TCut& selection, Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   virtual Long64_t        Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0); // *MENU*
   virtual void            DropBaskets();
   virtual void            DropBuffers(Int_t nbytes);
   virtual Int_t           Fill();
   virtual TBranch        *FindBranch(const char* name);
   virtual TLeaf          *FindLeaf(const char* name);
   virtual Int_t           Fit(const char* funcname, const char* varexp, const char* selection = "", Option_t* option = "", Option_t* goption = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0); // *MENU*
   virtual Int_t           FlushBaskets(Bool_t create_cluster = true) const;
   virtual const char     *GetAlias(const char* aliasName) const;
   UInt_t                  GetAllocationCount() const { return fAllocationCount; }
#ifdef R__TRACK_BASKET_ALLOC_TIME
   ULong64_t               GetAllocationTime() const { return fAllocationTime; }
#endif
   virtual Long64_t        GetAutoFlush() const {return fAutoFlush;}
   virtual Long64_t        GetAutoSave()  const {return fAutoSave;}
   virtual TBranch        *GetBranch(const char* name);
   virtual TBranchRef     *GetBranchRef() const { return fBranchRef; };
   virtual Bool_t          GetBranchStatus(const char* branchname) const;
   static  Int_t           GetBranchStyle();
   virtual Long64_t        GetCacheSize() const { return fCacheSize; }
   virtual TClusterIterator GetClusterIterator(Long64_t firstentry);
   virtual Long64_t        GetChainEntryNumber(Long64_t entry) const { return entry; }
   virtual Long64_t        GetChainOffset() const { return fChainOffset; }
   virtual Bool_t          GetClusterPrefetch() const { return fCacheDoClusterPrefetch; }
   TFile                  *GetCurrentFile() const;
           Int_t           GetDefaultEntryOffsetLen() const {return fDefaultEntryOffsetLen;}
           Long64_t        GetDebugMax()  const { return fDebugMax; }
           Long64_t        GetDebugMin()  const { return fDebugMin; }
   TDirectory             *GetDirectory() const { return fDirectory; }
   virtual Long64_t        GetEntries() const   { return fEntries; }
   virtual Long64_t        GetEntries(const char *selection);
   virtual Long64_t        GetEntriesFast() const   { return fEntries; }
   virtual Long64_t        GetEntriesFriend() const;
   virtual Long64_t        GetEstimate() const { return fEstimate; }
   virtual Int_t           GetEntry(Long64_t entry = 0, Int_t getall = 0);
           Int_t           GetEvent(Long64_t entry = 0, Int_t getall = 0) { return GetEntry(entry, getall); }
   virtual Int_t           GetEntryWithIndex(Int_t major, Int_t minor = 0);
   virtual Long64_t        GetEntryNumberWithBestIndex(Long64_t major, Long64_t minor = 0) const;
   virtual Long64_t        GetEntryNumberWithIndex(Long64_t major, Long64_t minor = 0) const;
   TEventList             *GetEventList() const { return fEventList; }
   virtual TEntryList     *GetEntryList();
   virtual Long64_t        GetEntryNumber(Long64_t entry) const;
   virtual Int_t           GetFileNumber() const { return fFileNumber; }
   virtual TTree          *GetFriend(const char*) const;
   virtual const char     *GetFriendAlias(TTree*) const;
   TH1                    *GetHistogram() { return GetPlayer()->GetHistogram(); }
   virtual Bool_t          GetImplicitMT() { return fIMTEnabled; }
   virtual Int_t          *GetIndex() { return &fIndex.fArray[0]; }
   virtual Double_t       *GetIndexValues() { return &fIndexValues.fArray[0]; }
           ROOT::TIOFeatures GetIOFeatures() const;
   virtual TIterator      *GetIteratorOnAllLeaves(Bool_t dir = kIterForward);
   virtual TLeaf          *GetLeaf(const char* branchname, const char* leafname);
   virtual TLeaf          *GetLeaf(const char* name);
   virtual TList          *GetListOfClones() { return fClones; }
   virtual TObjArray      *GetListOfBranches() { return &fBranches; }
   virtual TObjArray      *GetListOfLeaves() { return &fLeaves; }
   virtual TList          *GetListOfFriends() const { return fFriends; }
   virtual TList          *GetListOfAliases() const { return fAliases; }

   // GetMakeClass is left non-virtual for efficiency reason.
   // Making it virtual affects the performance of the I/O
           Int_t           GetMakeClass() const { return fMakeClass; }

   virtual Long64_t        GetMaxEntryLoop() const { return fMaxEntryLoop; }
   virtual Double_t        GetMaximum(const char* columname);
   static  Long64_t        GetMaxTreeSize();
   virtual Long64_t        GetMaxVirtualSize() const { return fMaxVirtualSize; }
   virtual Double_t        GetMinimum(const char* columname);
   virtual Int_t           GetNbranches() { return fBranches.GetEntriesFast(); }
   TObject                *GetNotify() const { return fNotify; }
   TVirtualTreePlayer     *GetPlayer();
   virtual Int_t           GetPacketSize() const { return fPacketSize; }
   virtual TVirtualPerfStats *GetPerfStats() const { return fPerfStats; }
           TTreeCache     *GetReadCache(TFile *file) const;
           TTreeCache     *GetReadCache(TFile *file, Bool_t create);
   virtual Long64_t        GetReadEntry()  const { return fReadEntry; }
   virtual Long64_t        GetReadEvent()  const { return fReadEntry; }
   virtual Int_t           GetScanField()  const { return fScanField; }
   TTreeFormula           *GetSelect()    { return GetPlayer()->GetSelect(); }
   virtual Long64_t        GetSelectedRows() { return GetPlayer()->GetSelectedRows(); }
   virtual Int_t           GetTimerInterval() const { return fTimerInterval; }
           TBuffer*        GetTransientBuffer(Int_t size);
   virtual Long64_t        GetTotBytes() const { return fTotBytes; }
   virtual TTree          *GetTree() const { return const_cast<TTree*>(this); }
   virtual TVirtualIndex  *GetTreeIndex() const { return fTreeIndex; }
   virtual Int_t           GetTreeNumber() const { return 0; }
   Float_t GetTargetMemoryRatio() const { return fTargetMemoryRatio; }
   virtual Int_t           GetUpdate() const { return fUpdate; }
   virtual TList          *GetUserInfo();
   // See TSelectorDraw::GetVar
   TTreeFormula           *GetVar(Int_t i)  { return GetPlayer()->GetVar(i); }
   // See TSelectorDraw::GetVar
   TTreeFormula           *GetVar1() { return GetPlayer()->GetVar1(); }
   // See TSelectorDraw::GetVar
   TTreeFormula           *GetVar2() { return GetPlayer()->GetVar2(); }
   // See TSelectorDraw::GetVar
   TTreeFormula           *GetVar3() { return GetPlayer()->GetVar3(); }
   // See TSelectorDraw::GetVar
   TTreeFormula           *GetVar4() { return GetPlayer()->GetVar4(); }
   // See TSelectorDraw::GetVal
   virtual Double_t       *GetVal(Int_t i)   { return GetPlayer()->GetVal(i); }
   // See TSelectorDraw::GetVal
   virtual Double_t       *GetV1()   { return GetPlayer()->GetV1(); }
   // See TSelectorDraw::GetVal
   virtual Double_t       *GetV2()   { return GetPlayer()->GetV2(); }
   // See TSelectorDraw::GetVal
   virtual Double_t       *GetV3()   { return GetPlayer()->GetV3(); }
   // See TSelectorDraw::GetVal
   virtual Double_t       *GetV4()   { return GetPlayer()->GetV4(); }
   virtual Double_t       *GetW()    { return GetPlayer()->GetW(); }
   virtual Double_t        GetWeight() const   { return fWeight; }
   virtual Long64_t        GetZipBytes() const { return fZipBytes; }
   virtual void            IncrementTotalBuffers(Int_t nbytes) { fTotalBuffers += nbytes; }
   Bool_t                  IsFolder() const { return kTRUE; }
   virtual Int_t           LoadBaskets(Long64_t maxmemory = 2000000000);
   virtual Long64_t        LoadTree(Long64_t entry);
   virtual Long64_t        LoadTreeFriend(Long64_t entry, TTree* T);
   virtual Int_t           MakeClass(const char* classname = 0, Option_t* option = "");
   virtual Int_t           MakeCode(const char* filename = 0);
   virtual Int_t           MakeProxy(const char* classname, const char* macrofilename = 0, const char* cutfilename = 0, const char* option = 0, Int_t maxUnrolling = 3);
   virtual Int_t           MakeSelector(const char* selector = 0, Option_t* option = "");
   Bool_t                  MemoryFull(Int_t nbytes);
   virtual Long64_t        Merge(TCollection* list, Option_t* option = "");
   virtual Long64_t        Merge(TCollection* list, TFileMergeInfo *info);
   static  TTree          *MergeTrees(TList* list, Option_t* option = "");
   virtual Bool_t          Notify();
   virtual void            OptimizeBaskets(ULong64_t maxMemory=10000000, Float_t minComp=1.1, Option_t *option="");
   TPrincipal             *Principal(const char* varexp = "", const char* selection = "", Option_t* option = "np", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   virtual void            Print(Option_t* option = "") const; // *MENU*
   virtual void            PrintCacheStats(Option_t* option = "") const;
   virtual Long64_t        Process(const char* filename, Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0); // *MENU*
   virtual Long64_t        Process(TSelector* selector, Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   virtual Long64_t        Project(const char* hname, const char* varexp, const char* selection = "", Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   virtual TSQLResult     *Query(const char* varexp = "", const char* selection = "", Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   virtual Long64_t        ReadFile(const char* filename, const char* branchDescriptor = "", char delimiter = ' ');
   virtual Long64_t        ReadStream(std::istream& inputStream, const char* branchDescriptor = "", char delimiter = ' ');
   virtual void            Refresh();
   virtual void            RegisterExternalFriend(TFriendElement *);
   virtual void            RemoveExternalFriend(TFriendElement *);
   virtual void            RemoveFriend(TTree*);
   virtual void            RecursiveRemove(TObject *obj);
   virtual void            Reset(Option_t* option = "");
   virtual void            ResetAfterMerge(TFileMergeInfo *);
   virtual void            ResetBranchAddress(TBranch *);
   virtual void            ResetBranchAddresses();
   virtual Long64_t        Scan(const char* varexp = "", const char* selection = "", Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0); // *MENU*
   virtual Bool_t          SetAlias(const char* aliasName, const char* aliasFormula);
   virtual void            SetAutoSave(Long64_t autos = -300000000);
   virtual void            SetAutoFlush(Long64_t autof = -30000000);
   virtual void            SetBasketSize(const char* bname, Int_t buffsize = 16000);
   virtual Int_t           SetBranchAddress(const char *bname,void *add, TBranch **ptr = 0);
   virtual Int_t           SetBranchAddress(const char *bname,void *add, TClass *realClass, EDataType datatype, Bool_t isptr);
   virtual Int_t           SetBranchAddress(const char *bname,void *add, TBranch **ptr, TClass *realClass, EDataType datatype, Bool_t isptr);
   template <class T> Int_t SetBranchAddress(const char *bname, T **add, TBranch **ptr = 0) {
      TClass *cl = TClass::GetClass<T>();
      EDataType type = kOther_t;
      if (cl==0) type = TDataType::GetType(typeid(T));
      return SetBranchAddress(bname,add,ptr,cl,type,true);
   }
#ifndef R__NO_CLASS_TEMPLATE_SPECIALIZATION
   // This can only be used when the template overload resolution can distinguish between
   // T* and T**
   template <class T> Int_t SetBranchAddress(const char *bname, T *add, TBranch **ptr = 0) {
      TClass *cl = TClass::GetClass<T>();
      EDataType type = kOther_t;
      if (cl==0) type = TDataType::GetType(typeid(T));
      return SetBranchAddress(bname,add,ptr,cl,type,false);
   }
#endif
   virtual void            SetBranchStatus(const char* bname, Bool_t status = 1, UInt_t* found = 0);
   static  void            SetBranchStyle(Int_t style = 1);  //style=0 for old branch, =1 for new branch style
   virtual Int_t           SetCacheSize(Long64_t cachesize = -1);
   virtual Int_t           SetCacheEntryRange(Long64_t first, Long64_t last);
   virtual void            SetCacheLearnEntries(Int_t n=10);
   virtual void            SetChainOffset(Long64_t offset = 0) { fChainOffset=offset; }
   virtual void            SetCircular(Long64_t maxEntries);
   virtual void            SetClusterPrefetch(Bool_t enabled) { fCacheDoClusterPrefetch = enabled; }
   virtual void            SetDebug(Int_t level = 1, Long64_t min = 0, Long64_t max = 9999999); // *MENU*
   virtual void            SetDefaultEntryOffsetLen(Int_t newdefault, Bool_t updateExisting = kFALSE);
   virtual void            SetDirectory(TDirectory* dir);
   virtual Long64_t        SetEntries(Long64_t n = -1);
   virtual void            SetEstimate(Long64_t nentries = 1000000);
           ROOT::TIOFeatures SetIOFeatures(const ROOT::TIOFeatures &);
   virtual void            SetFileNumber(Int_t number = 0);
   virtual void            SetEventList(TEventList* list);
   virtual void            SetEntryList(TEntryList* list, Option_t *opt="");
   virtual void            SetImplicitMT(Bool_t enabled) { fIMTEnabled = enabled; }
   virtual void            SetMakeClass(Int_t make);
   virtual void            SetMaxEntryLoop(Long64_t maxev = kMaxEntries) { fMaxEntryLoop = maxev; } // *MENU*
   static  void            SetMaxTreeSize(Long64_t maxsize = 100000000000LL);
   virtual void            SetMaxVirtualSize(Long64_t size = 0) { fMaxVirtualSize = size; } // *MENU*
   virtual void            SetName(const char* name); // *MENU*

   /**
    * @brief Sets the address of the object to be notified when the tree is loaded.
    *
    * The method TObject::Notify is called for the given object when the tree
    * is loaded. Specifically this occurs in the TTree::LoadTree method. To
    * remove the notification call this method with nullptr:
    * @code tree->SetNotify(nullptr); @endcode
    *
    * @param[in] obj Pointer to a TObject to be notified.
    */
   virtual void            SetNotify(TObject* obj) { fNotify = obj; }

   virtual void            SetObject(const char* name, const char* title);
   virtual void            SetParallelUnzip(Bool_t opt=kTRUE, Float_t RelSize=-1);
   virtual void            SetPerfStats(TVirtualPerfStats* perf);
   virtual void            SetScanField(Int_t n = 50) { fScanField = n; } // *MENU*
   void SetTargetMemoryRatio(Float_t ratio) { fTargetMemoryRatio = ratio; }
   virtual void            SetTimerInterval(Int_t msec = 333) { fTimerInterval=msec; }
   virtual void            SetTreeIndex(TVirtualIndex* index);
   virtual void            SetWeight(Double_t w = 1, Option_t* option = "");
   virtual void            SetUpdate(Int_t freq = 0) { fUpdate = freq; }
   virtual void            Show(Long64_t entry = -1, Int_t lenmax = 20);
   virtual void            StartViewer(); // *MENU*
   virtual Int_t           StopCacheLearningPhase();
   virtual Int_t           UnbinnedFit(const char* funcname, const char* varexp, const char* selection = "", Option_t* option = "", Long64_t nentries = kMaxEntries, Long64_t firstentry = 0);
   void                    UseCurrentStyle();
   virtual Int_t           Write(const char *name=0, Int_t option=0, Int_t bufsize=0);
   virtual Int_t           Write(const char *name=0, Int_t option=0, Int_t bufsize=0) const;

   ClassDef(TTree, 20) // Tree descriptor (the main ROOT I/O class)
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTreeFriendLeafIter                                                  //
//                                                                      //
// Iterator on all the leaves in a TTree and its friend                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TTreeFriendLeafIter : public TIterator {

protected:
   TTree             *fTree;         ///< tree being iterated
   TIterator         *fLeafIter;     ///< current leaf sub-iterator.
   TIterator         *fTreeIter;     ///< current tree sub-iterator.
   Bool_t             fDirection;    ///< iteration direction

   TTreeFriendLeafIter() : fTree(0), fLeafIter(0), fTreeIter(0),
       fDirection(0) { }

public:
   TTreeFriendLeafIter(const TTree* t, Bool_t dir = kIterForward);
   TTreeFriendLeafIter(const TTreeFriendLeafIter &iter);
   ~TTreeFriendLeafIter() { SafeDelete(fLeafIter); SafeDelete(fTreeIter); }
   TIterator &operator=(const TIterator &rhs);
   TTreeFriendLeafIter &operator=(const TTreeFriendLeafIter &rhs);

   const TCollection *GetCollection() const { return 0; }
   Option_t          *GetOption() const;
   TObject           *Next();
   void               Reset() { SafeDelete(fLeafIter); SafeDelete(fTreeIter); }
   Bool_t operator !=(const TIterator&) const {
      // TODO: Implement me
      return false;
   }
   Bool_t operator !=(const TTreeFriendLeafIter&) const {
      // TODO: Implement me
      return false;
   }
   TObject *operator*() const {
      // TODO: Implement me
      return nullptr;
   }
   ClassDef(TTreeFriendLeafIter,0)  //Linked list iterator
 };


#endif
