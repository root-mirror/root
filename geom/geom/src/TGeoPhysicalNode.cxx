// @(#)root/geom:$Id$
// Author: Andrei Gheata   17/02/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoPhysicalNode
\ingroup Geometry_classes

Physical nodes are the actual 'touchable' objects in the geometry, representing
a path of positioned volumes starting with the top node:
     path=/TOP/A_1/B_4/C_3 , where A, B, C represent names of volumes.

The number of physical nodes is given by the total number of possible of
branches in the geometry hierarchy. In case of detector geometries and
specially for calorimeters this number can be of the order 1e6-1e9, therefore
it is impossible to create all physical nodes as objects in memory. In TGeo,
physical nodes are represented by the class TGeoPhysicalNode and can be created
on demand for alignment purposes:

~~~ {.cpp}
   TGeoPhysicalNode *pn = new TGeoPhysicalNode("path_to_object")
~~~

Once created, a physical node can be misaligned, meaning that its position
or even shape can be changed:

~~~ {.cpp}
   pn->Align(TGeoMatrix* newmat, TGeoShape* newshape, Bool_t check=kFALSE)
~~~
*/

/** \class TGeoPNEntry
\ingroup Geometry_classes

The knowledge of the path to the objects that need to be misaligned is
essential since there is no other way of identifying them. One can however
create 'symbolic links' to any complex path to make it more representable
for the object it designates:

~~~ {.cpp}
   TGeoPNEntry *pne = new TGeoPNEntry("TPC_SECTOR_2", "path_to_tpc_sect2");
   pne->SetPhysicalNode(pn)
~~~

Such a symbolic link hides the complexity of the path to the align object and
replaces it with a more meaningful name. In addition, TGeoPNEntry objects are
faster to search by name and they may optionally store an additional user
matrix.

For more details please read the misalignment section in the Users Guide.
*/

#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoCache.h"
#include "TGeoMatrix.h"
#include "TGeoShapeAssembly.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"
#include "TGeoVolume.h"

#include "TGeoPhysicalNode.h"

// statics and globals

ClassImp(TGeoPhysicalNode);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPhysicalNode::TGeoPhysicalNode() : TNamed()
{
   fLevel        = 0;
   fMatrices     = 0;
   fNodes        = 0;
   fMatrixOrig   = 0;
   SetVisibility(kTRUE);
   SetVisibleFull(kFALSE);
   SetIsVolAtt(kTRUE);
   SetAligned(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoPhysicalNode::TGeoPhysicalNode(const char *path) : TNamed(path,"")
{
   if (!path[0]) {
      Error("ctor", "path not valid");
      return;
   }
   fLevel  = 0;
   fMatrices = new TObjArray(30);
   fNodes    = new TObjArray(30);
   fMatrixOrig   = 0;
   SetPath(path);
   SetVisibility(kTRUE);
   SetVisibleFull(kFALSE);
   SetIsVolAtt(kTRUE);
   SetAligned(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPhysicalNode::~TGeoPhysicalNode()
{
   if (fMatrices) {
      fMatrices->Delete();
      delete fMatrices;
   }
   if (fNodes) delete fNodes;
   if (fMatrixOrig) delete fMatrixOrig;
}

////////////////////////////////////////////////////////////////////////////////
/// Align a physical node with a new relative matrix/shape.
/// Example: /TOP_1/A_1/B_1/C_1
///    node->Align(transl_1, box) will perform:
///    - change RELATIVE translation of C_1 node (with respect to its
///      container volume B) to transl_1
///    - change the shape of the C volume
/// *NOTE* The operations will affect ONLY the LAST node in the branch. All
///   volumes/nodes in the branch represented by this physical node are
///   CLONED so the operation does not affect other possible replicas.

Bool_t TGeoPhysicalNode::Align(TGeoMatrix *newmat, TGeoShape *newshape, Bool_t check, Double_t ovlp)
{
   if (!newmat && !newshape) return kFALSE;
   if (TGeoManager::IsLocked()) {
      Error("Align", "Not performed. Geometry in LOCKED mode !");
      return kFALSE;
   }
   if (newmat == gGeoIdentity) {
      Error("Align", "Cannot align using gGeoIdentity. Use some default matrix constructor to represent identities.");
      return kFALSE;
   }
   TGeoNode *node = GetNode();
   if (node->IsOffset()) {
      Error("Align", "Cannot align division nodes: %s\n",node->GetName());
      return kFALSE;
   }
   // Refresh the node since other Align calls may have altered the stored nodes
   Refresh();
   TGeoNode *nnode = 0;
   TGeoVolume *vm = GetVolume(0);
   TGeoVolume *vd = 0;
   Int_t i;
   if (!IsAligned()) {
      Int_t *id = new Int_t[fLevel];
      for (i=0; i<fLevel; i++) {
         // Store daughter indexes
         vd = GetVolume(i);
         node = GetNode(i+1);
         id[i] = vd->GetIndex(node);
         if (id[i]<0) {
            Error("Align","%s cannot align node %s",GetName(), node->GetName());
            delete [] id;
            return kFALSE;
         }
      }
      for (i=0; i<fLevel; i++) {
         // Get daughter node and its id inside vm
         node = GetNode(i+1);
         // Clone daughter volume and node if not done yet
         if (node->IsCloned()) {
            vd = node->GetVolume();
            nnode = node;
         } else {
            vd = node->GetVolume()->CloneVolume();
            if (!vd) {
               delete [] id;
               Fatal("Align", "Cannot clone volume %s", node->GetVolume()->GetName());
               return kFALSE;
            }
            nnode = node->MakeCopyNode();
            if (!nnode) {
               delete [] id;
               Fatal("Align", "Cannot make copy node for %s", node->GetName());
               return kFALSE;
            }
            // Correct pointers to mother and volume
            nnode->SetVolume(vd);
            nnode->SetMotherVolume(vm);
            // Decouple old node from mother volume and connect new one
            if (vm->TestBit(TGeoVolume::kVolumeImportNodes)) {
               gGeoManager->GetListOfGShapes()->Add(nnode);
            }
            vm->GetNodes()->RemoveAt(id[i]);
            vm->GetNodes()->AddAt(nnode,id[i]);
            fNodes->RemoveAt(i+1);
            fNodes->AddAt(nnode,i+1);
     //       node->GetVolume()->Release();
         }
         // Consider new cloned volume as mother and continue
         vm = vd;
      }
      delete [] id;
   } else {
      nnode = GetNode();
   }
   // Now nnode is a cloned node of the one that need to be aligned
   TGeoNodeMatrix *aligned = (TGeoNodeMatrix*)nnode;
   vm = nnode->GetMotherVolume();
   vd = nnode->GetVolume();
   if (newmat) {
      // Check if the old matrix for this node was shared
      Bool_t shared = kFALSE;
      Int_t nd = vm->GetNdaughters();
      TGeoCompositeShape *cs;
      if (nnode->GetMatrix()->IsShared()) {
         // Now find the node having a composite shape using this shared matrix
         for (i=0; i<nd; i++) {
            node = vm->GetNode(i);
            if (node==nnode) continue;
            if (node->IsOffset()) continue;
            if (!node->GetVolume()->GetShape()->IsComposite()) continue;
            // We found a node having a composite shape, scan for the shared matrix
            cs = (TGeoCompositeShape*)node->GetVolume()->GetShape();
            if (cs->GetBoolNode()->GetRightMatrix() != nnode->GetMatrix()) continue;
            // The composite uses the matrix -> replace it
            TGeoCompositeShape *ncs = new TGeoCompositeShape(cs->GetName(), cs->GetBoolNode()->MakeClone());
            ncs->GetBoolNode()->ReplaceMatrix(nnode->GetMatrix(), newmat);
            // We have to clone the node/volume having the composite shape
            TGeoVolume *newvol = node->GetVolume()->CloneVolume();
            if (!newvol) {
               Error("Align", "Cannot clone volume %s", node->GetVolume()->GetName());
               return kFALSE;
            }
            newvol->SetShape(ncs);
            TGeoNode *newnode = node->MakeCopyNode();
            if (!newnode) {
               Error("Align", "Cannot clone node %s", node->GetName());
               return kFALSE;
            }
            newnode->SetVolume(newvol);
            newnode->SetMotherVolume(vm);
            if (vm->TestBit(TGeoVolume::kVolumeImportNodes)) {
               gGeoManager->GetListOfGShapes()->Add(newnode);
            }
            vm->GetNodes()->RemoveAt(i);
            vm->GetNodes()->AddAt(newnode,i);
            shared = kTRUE;
         }
         if (!shared) Error("Align", "The matrix replaced for %s is not actually shared", GetName());
      } else {
         // The aligned node may have a composite shape containing a shared matrix
         if (vd->GetShape()->IsComposite()) {
            cs = (TGeoCompositeShape*)vd->GetShape();
            if (cs->GetBoolNode()->GetRightMatrix()->IsShared()) {
               if (!nnode->GetMatrix()->IsIdentity()) {
                  Error("Align", "The composite shape having a shared matrix on the subtracted branch must be positioned using identity matrix.");
                  return kFALSE;
               }
               // We have to put the alignment matrix on top of the left branch
               // of the composite shape. The node is already decoupled from logical tree.
               TGeoCompositeShape *ncs = new TGeoCompositeShape(cs->GetName(), cs->GetBoolNode()->MakeClone());
               TGeoMatrix *oldmat = ncs->GetBoolNode()->GetLeftMatrix();
               TGeoHMatrix *newmat1 = new TGeoHMatrix(*newmat);
               newmat1->Multiply(oldmat);
               ncs->GetBoolNode()->ReplaceMatrix(oldmat, newmat1);
               vd->SetShape(ncs);
               // The right-side matrix pointer is preserved, so no need to update nodes.
               aligned = 0; // to prevent updating its matrix
            }
         }
      }
      // Register matrix and make it the active one
      if (!newmat->IsRegistered()) newmat->RegisterYourself();
      if (aligned) {
         aligned->SetMatrix(newmat);
         // Update the global matrix for the aligned node
         TGeoHMatrix *global = GetMatrix();
         TGeoHMatrix *up = GetMatrix(fLevel-1);
         *global = up;
         global->Multiply(newmat);
      }
   }
   // Change the shape for the aligned node
   if (newshape) vd->SetShape(newshape);

   // Re-compute bounding box of mother(s) if needed
   for (i=fLevel-1; i>0; i--) {
      Bool_t dassm = vd->IsAssembly(); // is daughter assembly ?
      vd = GetVolume(i);
      if (!vd) break;
      Bool_t cassm = vd->IsAssembly(); // is current assembly ?
      if (cassm) ((TGeoShapeAssembly*)vd->GetShape())->NeedsBBoxRecompute();
      if ((cassm || dassm) && vd->GetVoxels()) vd->GetVoxels()->SetNeedRebuild();
      if (!cassm) break;
   }

   // Now we have to re-voxelize the mother volume
   TGeoVoxelFinder *voxels = vm->GetVoxels();
   if (voxels) voxels->SetNeedRebuild();
   // Eventually check for overlaps
   if (check) {
      if (voxels) {
         voxels->Voxelize();
         vm->FindOverlaps();
      }
      // Set aligned node to be checked
      i = fLevel;
      node = GetNode(i);
      if (!node) return kTRUE;
      if (node->IsOverlapping()) {
         Info("Align", "The check for overlaps for node: \n%s\n cannot be performed since the node is declared possibly overlapping",
              GetName());
      } else {
         gGeoManager->SetCheckedNode(node);
         // Check overlaps for the first non-assembly parent node
         while ((node=GetNode(--i))) {
            if (!node->GetVolume()->IsAssembly()) break;
         }
         if (node && node->IsOverlapping()) {
            Info("Align", "The check for overlaps for assembly node: \n%s\n cannot be performed since the parent %s is declared possibly overlapping",
                 GetName(), node->GetName());
            node = 0;
         }
         if (node) node->CheckOverlaps(ovlp);
         gGeoManager->SetCheckedNode(0);
      }
   }
   // Clean current matrices from cache
   gGeoManager->CdTop();
   SetAligned(kTRUE);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

void TGeoPhysicalNode::cd() const
{
   if (GetNode(0) != gGeoManager->GetTopNode()) return;
   gGeoManager->cd(fName.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this node.

void TGeoPhysicalNode::Draw(Option_t * /*option*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Return parent at LEVUP generation

TGeoNode *TGeoPhysicalNode::GetMother(Int_t levup) const
{
   Int_t ind = fLevel-levup;
   if (ind<0) return 0;
   return (TGeoNode*)fNodes->UncheckedAt(ind);
}

////////////////////////////////////////////////////////////////////////////////
/// Return global matrix for node at LEVEL.

TGeoHMatrix *TGeoPhysicalNode::GetMatrix(Int_t level) const
{
   if (level<0) return (TGeoHMatrix*)fMatrices->UncheckedAt(fLevel);
   if (level>fLevel) return 0;
   return (TGeoHMatrix*)fMatrices->UncheckedAt(level);
}

////////////////////////////////////////////////////////////////////////////////
/// Return node in branch at LEVEL. If not specified, return last leaf.

TGeoNode *TGeoPhysicalNode::GetNode(Int_t level) const
{
   if (level<0) return (TGeoNode*)fNodes->UncheckedAt(fLevel);
   if (level>fLevel) return 0;
   return (TGeoNode*)fNodes->UncheckedAt(level);
}

////////////////////////////////////////////////////////////////////////////////
/// Return volume associated with node at LEVEL in the branch

TGeoVolume *TGeoPhysicalNode::GetVolume(Int_t level) const
{
   TGeoNode *node = GetNode(level);
   if (node) return node->GetVolume();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return shape associated with volume.

TGeoShape *TGeoPhysicalNode::GetShape(Int_t level) const
{
   TGeoVolume *vol = GetVolume(level);
   if (vol) return vol->GetShape();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this node and its content according to visualization settings.

void TGeoPhysicalNode::Paint(Option_t * /*option*/)
{
   TVirtualGeoPainter *painter = gGeoManager->GetGeomPainter();
   if (!painter) return;
//   painter->PaintNode(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about this node.

void TGeoPhysicalNode::Print(Option_t * /*option*/) const
{
   printf("TGeoPhysicalNode: %s level=%d aligned=%d\n", fName.Data(), fLevel, IsAligned());
   for (Int_t i=0; i<=fLevel; i++) {
      printf(" level %d: node %s\n", i, GetNode(i)->GetName());
      printf(" local matrix:\n");
      if (GetNode(i)->GetMatrix()->IsIdentity()) printf("   IDENTITY\n");
      else GetNode(i)->GetMatrix()->Print();
      printf(" global matrix:\n");
      if (GetMatrix(i)->IsIdentity()) printf("   IDENTITY\n");
      else GetMatrix(i)->Print();
   }
   if (IsAligned() && fMatrixOrig) {
      printf(" original local matrix:\n");
      fMatrixOrig->Print();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh this physical node. Called for all registered physical nodes
/// after an Align() call.

void TGeoPhysicalNode::Refresh()
{
   SetPath(fName.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Set node branch according to current state

void TGeoPhysicalNode::SetBranchAsState()
{
   TGeoNodeCache *cache = gGeoManager->GetCache();
   if (!cache) {
      Error("SetBranchAsState","no state available");
      return;
   }
   if (!cache->IsDummy()) {
      Error("SetBranchAsState", "not implemented for full cache");
      return;
   }
   if (!fNodes)    fNodes = new TObjArray(30);
   if (!fMatrices) fMatrices = new TObjArray(30);
   TGeoHMatrix **matrices = (TGeoHMatrix **) cache->GetMatrices();
   TGeoNode **branch = (TGeoNode **) cache->GetBranch();

   Bool_t refresh = (fLevel>0)?kTRUE:kFALSE;
   if (refresh) {
      TGeoHMatrix *current;
      for (Int_t i=0; i<=fLevel; i++) {
         fNodes->AddAtAndExpand(branch[i],i);
         current = (TGeoHMatrix*)fMatrices->UncheckedAt(i);
         *current = *matrices[i];
      }
      return;
   }
   fLevel = gGeoManager->GetLevel();
   for (Int_t i=0; i<=fLevel; i++) {
      fNodes->AddAtAndExpand(branch[i],i);
      fMatrices->AddAtAndExpand(new TGeoHMatrix(*matrices[i]),i);
   }
   TGeoNode *node = (TGeoNode*)fNodes->UncheckedAt(fLevel);
   if (!fMatrixOrig) fMatrixOrig = new TGeoHMatrix();
   *fMatrixOrig = node->GetMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Allows PN entries (or users) to preset the local original matrix for the
/// last node pointed by the path.

void TGeoPhysicalNode::SetMatrixOrig(const TGeoMatrix *local)
{
   if (!fMatrixOrig) fMatrixOrig = new TGeoHMatrix();
   if (!local) {
      fMatrixOrig->Clear();
      return;
   }
   *fMatrixOrig = local;
}

////////////////////////////////////////////////////////////////////////////////
/// Specify the path for this node.

Bool_t TGeoPhysicalNode::SetPath(const char *path)
{
   if (!gGeoManager->cd(path)) {
      Error("SetPath","wrong path -> maybe RestoreMasterVolume");
      return kFALSE;
   }
   SetBranchAsState();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if a given navigator state matches this physical node

Bool_t TGeoPhysicalNode::IsMatchingState(TGeoNavigator *nav) const
{
   TGeoNodeCache *cache = nav->GetCache();
   if (!cache) {
      Fatal("SetBranchAsState","no state available");
      return kFALSE;
   }
   TGeoNode **branch = (TGeoNode **) cache->GetBranch();
   for (Int_t i=1; i<=fLevel; i++)
      if (fNodes->At(i) != branch[i]) return kFALSE;
   return kTRUE;
}

ClassImp(TGeoPNEntry);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPNEntry::TGeoPNEntry()
{
   fNode = 0;
   fMatrix = 0;
   fGlobalOrig = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPNEntry::TGeoPNEntry(const char *name, const char *path)
            :TNamed(name, path)
{
   if (!gGeoManager || !gGeoManager->IsClosed() || !gGeoManager->CheckPath(path)) {
      TString errmsg("Cannot define a physical node link without a closed geometry and a valid path !");
      Error("ctor", "%s", errmsg.Data());
      throw errmsg;
      return;
   }
   gGeoManager->PushPath();
   gGeoManager->cd(path);
   fGlobalOrig = new TGeoHMatrix();
   *fGlobalOrig = gGeoManager->GetCurrentMatrix();
   gGeoManager->PopPath();
   fNode = 0;
   fMatrix = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPNEntry::~TGeoPNEntry()
{
   if (fMatrix && !fMatrix->IsRegistered()) delete fMatrix;
   delete fGlobalOrig;
}

////////////////////////////////////////////////////////////////////////////////
/// Setter for the corresponding physical node.

void TGeoPNEntry::SetPhysicalNode(TGeoPhysicalNode *node)
{
   if (fNode && node) {
      Warning("SetPhysicalNode", "Physical node changed for entry %s", GetName());
      Warning("SetPhysicalNode", "=== New path: %s", node->GetName());
   }
   fNode = node;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the additional matrix for this node entry. The matrix will be deleted
/// by this class unless registered by the user to gGeoManager

void TGeoPNEntry::SetMatrix(const TGeoHMatrix *mat)
{
   fMatrix = mat;
}
