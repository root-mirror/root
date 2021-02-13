// @(#)root/geom:$Id$
// Author: Andrei Gheata   18/03/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
// Author : Andrei Gheata -           date Wed 12 Dec 2001 09:45:08 AM CET

#ifndef ROOT_TGeoCache
#define ROOT_TGeoCache

#include "TGeoNode.h"

#ifdef R__LESS_INCLUDES
struct TGeoStateInfo;
#else
#include "TGeoStateInfo.h"
#endif

// forward declarations
class TGeoManager;
class TGeoHMatrix;

class TGeoCacheState : public TObject
{
protected:
   Int_t                fCapacity;      // maximum level stored
   Int_t                fLevel;         // level in the current branch
   Int_t                fNmany;         // number of overlapping nodes on current branch
   Int_t                fStart;         // start level
   Int_t                fIdBranch[30];  // ID branch
   Double_t             fPoint[3];      // last point in master frame
   Bool_t               fOverlapping;   // overlap flag

   TGeoNode           **fNodeBranch;    // last node branch stored
   TGeoHMatrix        **fMatrixBranch;  // global matrices for last branch
   TGeoHMatrix        **fMatPtr;        // array of matrix pointers

   TGeoCacheState(const TGeoCacheState&);
   TGeoCacheState& operator=(const TGeoCacheState&);

public:
   TGeoCacheState();
   TGeoCacheState(Int_t capacity);
   virtual ~TGeoCacheState();

   void                 SetState(Int_t level, Int_t startlevel, Int_t nmany, Bool_t ovlp, Double_t *point=0);
   Bool_t               GetState(Int_t &level, Int_t &nmany, Double_t *point) const;

   ClassDef(TGeoCacheState, 0)       // class storing the cache state
};

class TGeoNodeCache : public TObject
{
private:
   Int_t                 fGeoCacheMaxLevels;// maximum supported number of levels
   Int_t                 fGeoCacheStackSize;// maximum size of the stack
   Int_t                 fGeoInfoStackSize; // maximum size of the stack of info states
   Int_t                 fLevel;            // level in the current branch
   Int_t                 fStackLevel;       // current level in the stack
   Int_t                 fInfoLevel;        // current level in the stack
   Int_t                 fCurrentID;        // unique ID of current node
   Int_t                 fIndex;            // index in array of ID's
   Int_t                 fIdBranch[100];    // current branch of indices
   TString               fPath;             // path for current branch
   TGeoNode             *fTop;              // top node
   TGeoNode             *fNode;             //! current node
   TGeoHMatrix          *fMatrix;           //! current matrix
   TObjArray            *fStack;            // stack of cache states
   TGeoHMatrix         **fMatrixBranch;     // current branch of global matrices
   TGeoHMatrix         **fMPB;              // pre-built matrices
   TGeoNode            **fNodeBranch;       // current branch of nodes
   TGeoStateInfo       **fInfoBranch;       // current branch of nodes
   TGeoStateInfo        *fPWInfo;           //! State info for the parallel world
   Int_t                *fNodeIdArray;      //! array of node id's

   TGeoNodeCache(const TGeoNodeCache&) = delete;
   TGeoNodeCache& operator=(const TGeoNodeCache&) = delete;

public:
   TGeoNodeCache();
   TGeoNodeCache(TGeoNode *top, Bool_t nodeid=kFALSE, Int_t capacity=30);
   virtual ~TGeoNodeCache();

   void                 BuildIdArray();
   void                 BuildInfoBranch();
   void                 CdNode(Int_t nodeid);
   Bool_t               CdDown(Int_t index);
   Bool_t               CdDown(TGeoNode *node);
   void                 CdTop() {fLevel=1; CdUp();}
   void                 CdUp();
   void                 FillIdBranch(const Int_t *br, Int_t startlevel=0) {memcpy(fIdBranch+startlevel,br,(fLevel+1-startlevel)*sizeof(Int_t)); fIndex=fIdBranch[fLevel];}
   const Int_t         *GetIdBranch() const {return fIdBranch;}
   void                *GetBranch() const   {return fNodeBranch;}
   void                 GetBranchNames(Int_t *names) const;
   void                 GetBranchNumbers(Int_t *copyNumbers, Int_t *volumeNumbers) const;
   void                 GetBranchOnlys(Int_t *isonly) const;
   void                *GetMatrices() const {return fMatrixBranch;}
   TGeoHMatrix         *GetCurrentMatrix() const {return fMatrix;}
   Int_t                GetCurrentNodeId() const;
   TGeoNode            *GetMother(Int_t up=1) const {return ((fLevel-up)>=0)?fNodeBranch[fLevel-up]:0;}
   TGeoHMatrix         *GetMotherMatrix(Int_t up=1) const {return ((fLevel-up)>=0)?fMatrixBranch[fLevel-up]:0;}
   TGeoNode            *GetNode() const        {return fNode;}
   TGeoNode            *GetTopNode() const     {return fTop;}
   TGeoStateInfo       *GetInfo();
   TGeoStateInfo       *GetMakePWInfo(Int_t nd);
   void                 ReleaseInfo();
   Int_t                GetLevel() const       {return fLevel;}
   const char          *GetPath();
   Int_t                GetStackLevel() const  {return fStackLevel;}
   Int_t                GetNodeId() const;
   Bool_t               HasIdArray() const { return fNodeIdArray ? kTRUE : kFALSE; }
   Bool_t               IsDummy() const {return kTRUE;}

   void                 LocalToMaster(const Double_t *local, Double_t *master) const;
   void                 MasterToLocal(const Double_t *master, Double_t *local) const;
   void                 LocalToMasterVect(const Double_t *local, Double_t *master) const;
   void                 MasterToLocalVect(const Double_t *master, Double_t *local) const;
   void                 LocalToMasterBomb(const Double_t *local, Double_t *master) const;
   void                 MasterToLocalBomb(const Double_t *master, Double_t *local) const;
   Int_t                PushState(Bool_t ovlp, Int_t ntmany=0, Int_t startlevel=0, Double_t *point=0);
   Bool_t               PopState(Int_t &nmany, Double_t *point=0);
   Bool_t               PopState(Int_t &nmany, Int_t level, Double_t *point=0);
   void                 PopDummy(Int_t ipop=9999) {fStackLevel=(ipop>fStackLevel)?(fStackLevel-1):(ipop-1);}
   void                 Refresh() {fNode=fNodeBranch[fLevel]; fMatrix=fMatrixBranch[fLevel];}
   Bool_t               RestoreState(Int_t &nmany, TGeoCacheState *state, Double_t *point=0);

   ClassDef(TGeoNodeCache, 0)        // cache of reusable physical nodes
};

#endif
