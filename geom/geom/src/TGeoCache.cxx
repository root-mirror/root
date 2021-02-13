// @(#)root/geom:$Id$
// Author: Andrei Gheata   18/03/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoCache.h"

#include "TGeoManager.h"
#include "TGeoStateInfo.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TObject.h"

//const Int_t kN3 = 3*sizeof(Double_t);

ClassImp(TGeoNodeCache);

/** \class TGeoNodeCache
\ingroup Geometry_classes

Special pool of reusable nodes

*/

////////////////////////////////////////////////////////////////////////////////
/// Dummy constructor

TGeoNodeCache::TGeoNodeCache()
{
   fGeoCacheMaxLevels    = 100;
   fGeoCacheStackSize    = 10;
   fGeoInfoStackSize     = 100;
   fLevel       = 0;
   fStackLevel  = 0;
   fInfoLevel   = 0;
   fCurrentID   = 0;
   fIndex       = 0;
   fPath        = "";
   fTop         = 0;
   fNode        = 0;
   fMatrix      = 0;
   fStack       = 0;
   fMatrixBranch = 0;
   fMPB         = 0;
   fNodeBranch  = 0;
   fInfoBranch  = 0;
   fPWInfo      = 0;
   fNodeIdArray = 0;
   for (Int_t i=0; i<100; i++) fIdBranch[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoNodeCache::TGeoNodeCache(TGeoNode *top, Bool_t nodeid, Int_t capacity)
{
   fGeoCacheMaxLevels    = capacity;
   fGeoCacheStackSize    = 10;
   fGeoInfoStackSize     = 100;
   fLevel       = 0;
   fStackLevel  = 0;
   fInfoLevel   = 0;
   fCurrentID   = 0;
   fIndex       = 0;
   fPath        = "";
   fTop         = top;
   fNode        = top;
   fStack = new TObjArray(fGeoCacheStackSize);
   for (Int_t ist=0; ist<fGeoCacheStackSize; ist++)
      fStack->Add(new TGeoCacheState(fGeoCacheMaxLevels)); // !obsolete 100
   fMatrixBranch = new TGeoHMatrix *[fGeoCacheMaxLevels];
   fMPB = new TGeoHMatrix *[fGeoCacheMaxLevels];
   fNodeBranch  = new TGeoNode*[fGeoCacheMaxLevels];
   fInfoBranch  = new TGeoStateInfo*[fGeoInfoStackSize];
   for (Int_t i=0; i<fGeoCacheMaxLevels; i++) {
      fMPB[i] = new TGeoHMatrix(TString::Format("global_%d",i));
      fMatrixBranch[i] = 0;
      fNodeBranch[i] = 0;
   }
   for (Int_t i=0; i<fGeoInfoStackSize; i++) {
      fInfoBranch[i] = 0;
   }
   fPWInfo      = 0;
   fMatrix = fMatrixBranch[0] = fMPB[0];
   fNodeBranch[0] = top;
   fNodeIdArray = 0;
   for (Int_t i=0; i<100; i++) fIdBranch[i] = 0;
   if (nodeid) BuildIdArray();
   CdTop();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoNodeCache::~TGeoNodeCache()
{
   if (fStack) {
      fStack->Delete();
      delete fStack;
   }
   if (fMatrixBranch) delete [] fMatrixBranch;
   if (fMPB) {
      for (Int_t i=0; i<fGeoCacheMaxLevels; i++) delete fMPB[i];
      delete [] fMPB;
   }
   delete [] fNodeBranch;
   if (fInfoBranch) {
      for (Int_t i=0; i<fGeoInfoStackSize; i++) delete fInfoBranch[i];
   }
   delete [] fInfoBranch;
   if (fNodeIdArray)  delete [] fNodeIdArray;
   delete fPWInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Builds node id array.

void TGeoNodeCache::BuildIdArray()
{
   Int_t nnodes = gGeoManager->GetNNodes();
   //if (nnodes>3E7) return;
   if (fNodeIdArray) delete [] fNodeIdArray;
   Info("BuildIDArray","--- node ID tracking enabled, size=%lu Bytes\n", ULong_t((2*nnodes+1)*sizeof(Int_t)));
   fNodeIdArray = new Int_t[2*nnodes+1];
   fNodeIdArray[0] = 0;
   Int_t ifree  = 1;
   Int_t nodeid = 0;
   gGeoManager->GetTopNode()->FillIdArray(ifree, nodeid, fNodeIdArray);
   gGeoManager->CdTop();
   fIdBranch[0] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Builds info branch. Navigation is possible only after this step.

void TGeoNodeCache::BuildInfoBranch()
{
   if (!fInfoBranch) fInfoBranch  = new TGeoStateInfo*[fGeoInfoStackSize];
   else if (fInfoBranch[0]) return;
   for (Int_t i=0; i<fGeoInfoStackSize; i++) {
      fInfoBranch[i] = new TGeoStateInfo();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the PW info, if none create one

TGeoStateInfo *TGeoNodeCache::GetMakePWInfo(Int_t nd)
{
   if (fPWInfo) return fPWInfo;
   fPWInfo = new TGeoStateInfo(nd);
   return fPWInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current path to point to the node having this id.
/// Node id has to be in range : 0 to fNNodes-1 (no check for performance reasons)

void TGeoNodeCache::CdNode(Int_t nodeid) {
   if (!fNodeIdArray) {
      Error("CdNode", "Navigation based on physical node unique id disabled.\n   To enable, use: gGeoManager->GetCache()->BuildIdArray()");
      return;
   }
   Int_t *arr = fNodeIdArray;
   if (nodeid == arr[fIndex]) return;
   while (fLevel>0) {
      gGeoManager->CdUp();
      if (nodeid == arr[fIndex]) return;
   }
   gGeoManager->CdTop();
   Int_t currentID = 0;
   Int_t nd = GetNode()->GetNdaughters();
   Int_t nabove, nbelow, middle;
   while (nodeid!=currentID && nd) {
      nabove = nd+1;
      nbelow = 0;
      while (nabove-nbelow > 1) {
         middle = (nabove+nbelow)>>1;
         currentID = arr[arr[fIndex+middle]];
         if (nodeid == currentID) {
            gGeoManager->CdDown(middle-1);
            return;
         }
         if (nodeid < currentID) nabove = middle;
         else                    nbelow = middle;
      }
      gGeoManager->CdDown(nbelow-1);
      currentID = arr[fIndex];
      nd = GetNode()->GetNdaughters();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make daughter INDEX of current node the active state. Compute global matrix.

Bool_t TGeoNodeCache::CdDown(Int_t index)
{
   TGeoNode *newnode = fNode->GetDaughter(index);
   if (!newnode) return kFALSE;
   fLevel++;
   if (fNodeIdArray) {
      fIndex = fNodeIdArray[fIndex+index+1];
      fIdBranch[fLevel] = fIndex;
   }
   fNode = newnode;
   fNodeBranch[fLevel] = fNode;
   TGeoMatrix  *local = newnode->GetMatrix();
   TGeoHMatrix *newmat = fMPB[fLevel];
   if (!local->IsIdentity()) {
      newmat->CopyFrom(fMatrix);
      newmat->Multiply(local);
      fMatrix = newmat;
   }
   fMatrixBranch[fLevel] = fMatrix;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make daughter INDEX of current node the active state. Compute global matrix.

Bool_t TGeoNodeCache::CdDown(TGeoNode *newnode)
{
   if (!newnode) return kFALSE;
   fLevel++;
   if (fNodeIdArray) {
      Int_t index = fNode->GetVolume()->GetIndex(newnode);
      fIndex = fNodeIdArray[fIndex+index+1];
      fIdBranch[fLevel] = fIndex;
   }
   fNode = newnode;
   fNodeBranch[fLevel] = fNode;
   TGeoMatrix  *local = newnode->GetMatrix();
   TGeoHMatrix *newmat = fMPB[fLevel];
   if (!local->IsIdentity()) {
      newmat->CopyFrom(fMatrix);
      newmat->Multiply(local);
      fMatrix = newmat;
   }
   fMatrixBranch[fLevel] = fMatrix;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make mother of current node the active state.

void TGeoNodeCache::CdUp()
{
   if (!fLevel) return;
   fLevel--;
   if (fNodeIdArray) fIndex = fIdBranch[fLevel];
   fNode = fNodeBranch[fLevel];
   fMatrix = fMatrixBranch[fLevel];
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a fixed ID for current physical node

Int_t TGeoNodeCache::GetCurrentNodeId() const
{
   if (fNodeIdArray) return fNodeIdArray[fIndex];
   return GetNodeId();
}

////////////////////////////////////////////////////////////////////////////////
/// Get unique node id.

Int_t TGeoNodeCache::GetNodeId() const
{
   Long_t id=0;
   for (Int_t level=0;level<fLevel+1; level++)
      id += (Long_t)fNodeBranch[level];
   return (Int_t)id;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill names with current branch volume names (4 char - used by GEANT3 interface).

void TGeoNodeCache::GetBranchNames(Int_t *names) const
{
   const char *name;
   for (Int_t i=0; i<fLevel+1; i++) {
      name = fNodeBranch[i]->GetVolume()->GetName();
      memcpy(&names[i], name, sizeof(Int_t));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill copy numbers of current branch nodes.

void TGeoNodeCache::GetBranchNumbers(Int_t *copyNumbers, Int_t *volumeNumbers) const
{
   for (Int_t i=0; i<fLevel+1; i++) {
      copyNumbers[i]   = fNodeBranch[i]->GetNumber();
      volumeNumbers[i] = fNodeBranch[i]->GetVolume()->GetNumber();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill copy numbers of current branch nodes.

void TGeoNodeCache::GetBranchOnlys(Int_t *isonly) const
{
   Bool_t ismany = kFALSE;
   for (Int_t i=0; i<fLevel+1; i++) {
      if (!fNodeBranch[i]->IsOffset()) ismany=fNodeBranch[i]->IsOverlapping();
      isonly[i] = (ismany)?0:1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get next state info pointer.

TGeoStateInfo *TGeoNodeCache::GetInfo()
{
   if (fInfoLevel==fGeoInfoStackSize-1) {
      TGeoStateInfo **infoBranch = new TGeoStateInfo*[2*fGeoInfoStackSize];
      memcpy(infoBranch, fInfoBranch, fGeoInfoStackSize*sizeof(TGeoStateInfo*));
      for (Int_t i=fGeoInfoStackSize; i<2*fGeoInfoStackSize; i++)
         infoBranch[i] = new TGeoStateInfo();
      delete [] fInfoBranch;
      fInfoBranch = infoBranch;
      fGeoInfoStackSize *= 2;
   }
   return fInfoBranch[fInfoLevel++];
}

////////////////////////////////////////////////////////////////////////////////
/// Release last used state info pointer.

void TGeoNodeCache::ReleaseInfo()
{
   fInfoLevel--;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current geometry path.

const char *TGeoNodeCache::GetPath()
{
   fPath = "";
   for (Int_t level=0;level<fLevel+1; level++) {
      fPath += "/";
      fPath += fNodeBranch[level]->GetName();
   }
   return fPath.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Push current state into heap.

Int_t TGeoNodeCache::PushState(Bool_t ovlp, Int_t startlevel, Int_t nmany, Double_t *point)
{
   if (fStackLevel>=fGeoCacheStackSize) {
      for (Int_t ist=0; ist<fGeoCacheStackSize; ist++)
         fStack->Add(new TGeoCacheState(fGeoCacheMaxLevels));
   }
   ((TGeoCacheState*)fStack->At(fStackLevel))->SetState(fLevel,startlevel,nmany,ovlp,point);
   return ++fStackLevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Pop next state/point from heap.

Bool_t TGeoNodeCache::PopState(Int_t &nmany, Double_t *point)
{
   if (!fStackLevel) return 0;
   Bool_t ovlp = ((TGeoCacheState*)fStack->At(--fStackLevel))->GetState(fLevel,nmany,point);
   Refresh();
//   return (fStackLevel+1);
   return ovlp;
}

////////////////////////////////////////////////////////////////////////////////
/// Pop next state/point from heap and restore matrices starting from LEVEL.

Bool_t TGeoNodeCache::PopState(Int_t &nmany, Int_t level, Double_t *point)
{
   if (level<=0) return 0;
   Bool_t ovlp = ((TGeoCacheState*)fStack->At(level-1))->GetState(fLevel,nmany,point);
   Refresh();
   return ovlp;
}

////////////////////////////////////////////////////////////////////////////////
/// Pop next state/point from a backed-up state.

Bool_t TGeoNodeCache::RestoreState(Int_t &nmany, TGeoCacheState *state, Double_t *point)
{
   Bool_t ovlp = state->GetState(fLevel,nmany,point);
   Refresh();
   return ovlp;
}

////////////////////////////////////////////////////////////////////////////////
/// Local point converted to master frame defined by current matrix.

void TGeoNodeCache::LocalToMaster(const Double_t *local, Double_t *master) const
{
   fMatrix->LocalToMaster(local, master);
}

////////////////////////////////////////////////////////////////////////////////
/// Point in master frame defined by current matrix converted to local one.

void TGeoNodeCache::MasterToLocal(const Double_t *master, Double_t *local) const
{
   fMatrix->MasterToLocal(master, local);
}

////////////////////////////////////////////////////////////////////////////////
/// Local vector converted to master frame defined by current matrix.

void TGeoNodeCache::LocalToMasterVect(const Double_t *local, Double_t *master) const
{
   fMatrix->LocalToMasterVect(local, master);
}

////////////////////////////////////////////////////////////////////////////////
/// Vector in master frame defined by current matrix converted to local one.

void TGeoNodeCache::MasterToLocalVect(const Double_t *master, Double_t *local) const
{
   fMatrix->MasterToLocalVect(master,local);
}

////////////////////////////////////////////////////////////////////////////////
/// Local point converted to master frame defined by current matrix and rescaled with bomb factor.

void TGeoNodeCache::LocalToMasterBomb(const Double_t *local, Double_t *master) const
{
   fMatrix->LocalToMasterBomb(local, master);
}

////////////////////////////////////////////////////////////////////////////////
/// Point in master frame defined by current matrix converted to local one and rescaled with bomb factor.

void TGeoNodeCache::MasterToLocalBomb(const Double_t *master, Double_t *local) const
{
   fMatrix->MasterToLocalBomb(master, local);
}

ClassImp(TGeoCacheState);

/** \class TGeoCacheState
\ingroup Geometry_classes

Class storing the state of the cache at a given moment

*/

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TGeoCacheState::TGeoCacheState()
{
   fCapacity = 0;
   fLevel = 0;
   fNmany = 0;
   fStart = 0;
   memset(fIdBranch, 0, 30*sizeof(Int_t));
   memset(fPoint, 0, 3*sizeof(Int_t));
   fOverlapping = kFALSE;
   fNodeBranch = 0;
   fMatrixBranch = 0;
   fMatPtr = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Ctor.

TGeoCacheState::TGeoCacheState(Int_t capacity)
{
   fCapacity = capacity;
   fLevel = 0;
   fNmany = 0;
   fStart = 0;
   memset(fIdBranch, 0, 30*sizeof(Int_t));
   memset(fPoint, 0, 3*sizeof(Int_t));
   fOverlapping = kFALSE;
   fNodeBranch = new TGeoNode *[capacity];
   fMatrixBranch = new TGeoHMatrix *[capacity];
   fMatPtr = new TGeoHMatrix *[capacity];
   for (Int_t i=0; i<capacity; i++) {
      fMatrixBranch[i] = new TGeoHMatrix("global");
      fNodeBranch[i] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoCacheState::TGeoCacheState(const TGeoCacheState& gcs) :
  TObject(gcs),
  fCapacity(gcs.fCapacity),
  fLevel(gcs.fLevel),
  fNmany(gcs.fNmany),
  fStart(gcs.fStart),
  fOverlapping(gcs.fOverlapping)
{
   Int_t i;
   for (i=0; i<3; i++) fPoint[i]=gcs.fPoint[i];
   for(i=0; i<30; i++) fIdBranch[i]=gcs.fIdBranch[i];
   fNodeBranch = new TGeoNode *[fCapacity];
   fMatrixBranch = new TGeoHMatrix *[fCapacity];
   fMatPtr = new TGeoHMatrix *[fCapacity];
   for (i=0; i<fCapacity; i++) {
      fNodeBranch[i] = gcs.fNodeBranch[i];
      fMatrixBranch[i] = new TGeoHMatrix(*gcs.fMatrixBranch[i]);
      fMatPtr[i] = gcs.fMatPtr[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoCacheState& TGeoCacheState::operator=(const TGeoCacheState& gcs)
{
   Int_t i;
   if(this!=&gcs) {
      TObject::operator=(gcs);
      fCapacity=gcs.fCapacity;
      fLevel=gcs.fLevel;
      fNmany=gcs.fNmany;
      fStart=gcs.fStart;
      for(i=0; i<30; i++) fIdBranch[i]=gcs.fIdBranch[i];
      for(i=0; i<3; i++) fPoint[i]=gcs.fPoint[i];
      fOverlapping=gcs.fOverlapping;
      fNodeBranch = new TGeoNode *[fCapacity];
      fMatrixBranch = new TGeoHMatrix *[fCapacity];
      fMatPtr = new TGeoHMatrix *[fCapacity];
      for (i=0; i<fCapacity; i++) {
         fNodeBranch[i] = gcs.fNodeBranch[i];
         fMatrixBranch[i] = new TGeoHMatrix(*gcs.fMatrixBranch[i]);
         fMatPtr[i] = gcs.fMatPtr[i];
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Dtor.

TGeoCacheState::~TGeoCacheState()
{
   if (fNodeBranch) {
      for (Int_t i=0; i<fCapacity; i++) {
         delete fMatrixBranch[i];
      }
      delete [] fNodeBranch;
      delete [] fMatrixBranch;
      delete [] fMatPtr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill current modeller state.

void TGeoCacheState::SetState(Int_t level, Int_t startlevel, Int_t nmany, Bool_t ovlp, Double_t *point)
{
   fLevel = level;
   fStart = startlevel;
   fNmany = nmany;
   TGeoNodeCache *cache = gGeoManager->GetCache();
   if (cache->HasIdArray()) memcpy(fIdBranch, cache->GetIdBranch()+fStart, (level+1-fStart)*sizeof(Int_t));
   TGeoNode **node_branch = (TGeoNode **) cache->GetBranch();
   TGeoHMatrix **mat_branch  = (TGeoHMatrix **) cache->GetMatrices();
   Int_t nelem = level+1-fStart;
   memcpy(fNodeBranch, node_branch+fStart, nelem*sizeof(TGeoNode *));
   memcpy(fMatPtr, mat_branch+fStart, nelem*sizeof(TGeoHMatrix *));
   TGeoHMatrix *last = 0;
   TGeoHMatrix *current;
   for (Int_t i=0; i<nelem; i++) {
      current = mat_branch[i+fStart];
      if (current == last) continue;
      *fMatrixBranch[i] = current;
      last = current;
   }
   fOverlapping = ovlp;
   if (point) memcpy(fPoint, point, 3*sizeof(Double_t));
}

////////////////////////////////////////////////////////////////////////////////
/// Restore a modeler state.

Bool_t TGeoCacheState::GetState(Int_t &level, Int_t &nmany, Double_t *point) const
{
   level = fLevel;
   nmany = fNmany;
   TGeoNodeCache *cache = gGeoManager->GetCache();
   if (cache->HasIdArray()) cache->FillIdBranch(fIdBranch, fStart);
   TGeoNode **node_branch = (TGeoNode **) cache->GetBranch();
   TGeoHMatrix **mat_branch  = (TGeoHMatrix **) cache->GetMatrices();
   Int_t nelem = level+1-fStart;
   memcpy(node_branch+fStart, fNodeBranch, nelem*sizeof(TGeoNode *));
   memcpy(mat_branch+fStart, fMatPtr, (level+1-fStart)*sizeof(TGeoHMatrix *));
   TGeoHMatrix *last = 0;
   TGeoHMatrix *current;
   for (Int_t i=0; i<nelem; i++) {
      current = mat_branch[i+fStart];
      if (current == last) continue;
      *current = fMatrixBranch[i];
      last = current;
   }
   if (point) memcpy(point, fPoint, 3*sizeof(Double_t));
   return fOverlapping;
}
