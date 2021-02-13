// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoNode
\ingroup Geometry_classes

  A node represent a volume positioned inside another.They store links to both
volumes and to the TGeoMatrix representing the relative positioning. Node are
never instantiated directly by users, but created as a result of volume operations.
Adding a volume named A with a given user ID inside a volume B will create a node
node named A_ID. This will be added to the list of nodes stored by B. Also,
when applying a division operation in N slices to a volume A, a list of nodes
B_1, B_2, ..., B_N is also created. A node B_i does not represent a unique
object in the geometry because its container A might be at its turn positioned
as node inside several other volumes. Only when a complete branch of nodes
is fully defined up to the top node in the geometry, a given path like:

   /TOP_1/.../A_3/B_7 will represent an unique object.

Its global transformation matrix can be computed as the pile-up of all local
transformations in its branch. We will therefore call "logical graph" the
hierarchy defined by nodes and volumes. The expansion of the logical graph by
all possible paths defines a tree structure where all nodes are unique
"touchable" objects. We will call this the "physical tree". Unlike the logical
graph, the physical tree can become a huge structure with several milions of nodes
in case of complex geometries, therefore it is not always a good idea to keep it
transient in memory. Since a the logical and physical structures are correlated, the
modeller rather keeps track only of the current branch, updating the current
global matrix at each change of the level in geometry. The current physical node
is not an object that can be asked for at a given moment, but rather represented
by the combination: current node + current global matrix. However, physical nodes
have unique ID's that can be retrieved for a given modeler state. These can be
fed back to the modeler in order to force a physical node to become current.
The advantage of this comes from the fact that all navigation queries check
first the current node, therefore knowing the location of a point in the
geometry can be saved as a starting state for later use.

  Nodes can be declared as "overlapping" in case they do overlap with other
nodes inside the same container or extrude this container. Non-overlapping
nodes can be created with:

~~~ {.cpp}
     TGeoVolume::AddNode(TGeoVolume *daughter, Int_t copy_No, TGeoMatrix *matr);
~~~

The creation of overlapping nodes can be done with a similar prototype:

~~~ {.cpp}
     TGeoVolume::AddNodeOverlap(same arguments);
~~~

When closing the geometry, overlapping nodes perform a check of possible
overlaps with their neighbours. These are stored and checked all the time
during navigation, therefore navigation is slower when embedding such nodes
into geometry.

  Node have visualization attributes as volume have. When undefined by users,
painting a node on a pad will take the corresponding volume attributes.

\image html geom_t_node.png
*/

#include <iostream>

#include "TBrowser.h"
#include "TObjArray.h"
#include "TStyle.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TGeoExtension.h"

// statics and globals

ClassImp(TGeoNode);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoNode::TGeoNode()
{
   fVolume       = 0;
   fMother       = 0;
   fNumber       = 0;
   fNovlp        = 0;
   fOverlaps     = 0;
   fUserExtension = 0;
   fFWExtension = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoNode::TGeoNode(const TGeoVolume *vol)
{
   if (!vol) {
      Error("ctor", "volume not specified");
      return;
   }
   fVolume       = (TGeoVolume*)vol;
   if (fVolume->IsAdded()) fVolume->SetReplicated();
   fVolume->SetAdded();
   fMother       = 0;
   fNumber       = 0;
   fNovlp        = 0;
   fOverlaps     = 0;
   fUserExtension = 0;
   fFWExtension = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoNode::~TGeoNode()
{
   if (fOverlaps) delete [] fOverlaps;
   if (fUserExtension) {fUserExtension->Release(); fUserExtension=0;}
   if (fFWExtension) {fFWExtension->Release(); fFWExtension=0;}
}

////////////////////////////////////////////////////////////////////////////////
/// How-to-browse for a node.

void TGeoNode::Browse(TBrowser *b)
{
   if (!b) return;
   if (!GetNdaughters()) return;
   TGeoNode *daughter;
   TString title;
   for (Int_t i=0; i<GetNdaughters(); i++) {
      daughter = GetDaughter(i);
      b->Add(daughter, daughter->GetName(), daughter->IsVisible());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of daughters. Nodes pointing to same volume counted
/// once if unique_volumes is set.

Int_t TGeoNode::CountDaughters(Bool_t unique_volumes)
{
   static Int_t icall = 0;
   Int_t counter = 0;
   // Count this node
   if (unique_volumes) {
      if (!fVolume->IsSelected()) {
         counter++;
         fVolume->SelectVolume(kFALSE);
      }
   } else counter++;
   icall++;
   Int_t nd = fVolume->GetNdaughters();
   // Count daughters recursively
   for (Int_t i=0; i<nd; i++) counter += GetDaughter(i)->CountDaughters(unique_volumes);
   icall--;
   // Un-mark volumes
   if (icall == 0) fVolume->SelectVolume(kTRUE);
   return counter;
}

////////////////////////////////////////////////////////////////////////////////
/// Check overlaps bigger than OVLP hierarchically, starting with this node.

void TGeoNode::CheckOverlaps(Double_t ovlp, Option_t *option)
{
   Int_t icheck = 0;
   Int_t ncheck = 0;
   TStopwatch *timer;
   Int_t i;
   Bool_t sampling = kFALSE;
   TString opt(option);
   opt.ToLower();
   if (opt.Contains("s")) sampling = kTRUE;

   TGeoManager *geom = fVolume->GetGeoManager();
   ncheck = CountDaughters(kFALSE);
   timer = new TStopwatch();
   geom->ClearOverlaps();
   geom->SetCheckingOverlaps(kTRUE);
   Info("CheckOverlaps", "Checking overlaps for %s and daughters within %g", fVolume->GetName(),ovlp);
   if (sampling) {
      Info("CheckOverlaps", "Checking overlaps by sampling <%s> for %s and daughters", option, fVolume->GetName());
      Info("CheckOverlaps", "=== NOTE: Extrusions NOT checked with sampling option ! ===");
   }
   timer->Start();
   geom->GetGeomPainter()->OpProgress(fVolume->GetName(),icheck,ncheck,timer,kFALSE);
   fVolume->CheckOverlaps(ovlp,option);
   icheck++;
   TGeoIterator next(fVolume);
   TGeoNode *node;
   TString path;
   TObjArray *overlaps = geom->GetListOfOverlaps();
   Int_t novlps;
   TString msg;
   while ((node=next())) {
      next.GetPath(path);
      icheck++;
      if (!node->GetVolume()->IsSelected()) {
         msg = TString::Format("found %d overlaps", overlaps->GetEntriesFast());
         geom->GetGeomPainter()->OpProgress(node->GetVolume()->GetName(),icheck,ncheck,timer,kFALSE, msg);
         node->GetVolume()->SelectVolume(kFALSE);
         node->GetVolume()->CheckOverlaps(ovlp,option);
      }
   }
   fVolume->SelectVolume(kTRUE);
   geom->SetCheckingOverlaps(kFALSE);
   geom->SortOverlaps();
   novlps = overlaps->GetEntriesFast();
   TNamed *obj;
   for (i=0; i<novlps; i++) {
      obj = (TNamed*)overlaps->At(i);
      obj->SetName(TString::Format("ov%05d",i));
   }
   geom->GetGeomPainter()->OpProgress("Check overlaps:",icheck,ncheck,timer,kTRUE);
   Info("CheckOverlaps", "Number of illegal overlaps/extrusions : %d\n", novlps);
   delete timer;
}

////////////////////////////////////////////////////////////////////////////////
/// compute the closest distance of approach from point px,py to this node

Int_t TGeoNode::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t dist = 9999;
   if (!fVolume) return dist;
   if (gGeoManager != fVolume->GetGeoManager()) gGeoManager = fVolume->GetGeoManager();
   TVirtualGeoPainter *painter = gGeoManager->GetPainter();
   if (!painter) return dist;
   dist = painter->DistanceToPrimitiveVol(fVolume, px, py);
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse actions on this volume.

void TGeoNode::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!fVolume) return;
   TVirtualGeoPainter *painter = fVolume->GetGeoManager()->GetPainter();
   if (!painter) return;
   painter->ExecuteVolumeEvent(fVolume, event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Get node info for the browser.

char *TGeoNode::GetObjectInfo(Int_t px, Int_t py) const
{
   if (!fVolume) return 0;
   TVirtualGeoPainter *painter = fVolume->GetGeoManager()->GetPainter();
   if (!painter) return 0;
   return (char*)painter->GetVolumeInfo(fVolume, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// check if this node is drawn. Assumes that this node is current

Bool_t TGeoNode::IsOnScreen() const
{
   if (fVolume->TestAttBit(TGeoAtt::kVisOnScreen)) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Inspect this node.

void TGeoNode::InspectNode() const
{
   printf("== Inspecting node %s ", GetName());
   if (fMother) printf("mother volume %s. ", fMother->GetName());
   if (IsOverlapping()) printf("(Node is MANY)\n");
   else printf("\n");
   if (fOverlaps && fMother) {
      printf("   possibly overlapping with : ");
      for (Int_t i=0; i<fNovlp; i++)
         printf(" %s ", fMother->GetNode(fOverlaps[i])->GetName());
      printf("\n");
   }
   printf("Transformation matrix:\n");
   TGeoMatrix *matrix = GetMatrix();
   if (GetMatrix()) matrix->Print();
   fVolume->Print();
}

////////////////////////////////////////////////////////////////////////////////
/// check for wrong parameters in shapes

void TGeoNode::CheckShapes()
{
   fVolume->CheckShapes();
   Int_t nd = GetNdaughters();
   if (!nd) return;
   for (Int_t i=0; i<nd; i++) fVolume->GetNode(i)->CheckShapes();
}

////////////////////////////////////////////////////////////////////////////////
/// draw only this node independently of its vis options

void TGeoNode::DrawOnly(Option_t *option)
{
   fVolume->DrawOnly(option);
}

////////////////////////////////////////////////////////////////////////////////
/// draw current node according to option

void TGeoNode::Draw(Option_t *option)
{
   gGeoManager->FindNode();
   gGeoManager->CdUp();
   Double_t point[3];
   gGeoManager->MasterToLocal(gGeoManager->GetCurrentPoint(), &point[0]);
   gGeoManager->SetCurrentPoint(&point[0]);
   gGeoManager->GetCurrentVolume()->Draw(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Method drawing the overlap candidates with this node.

void TGeoNode::DrawOverlaps()
{
   if (!fNovlp) {printf("node %s is ONLY\n", GetName()); return;}
   if (!fOverlaps) {printf("node %s no overlaps\n", GetName()); return;}
   TGeoNode *node;
   Int_t i;
   Int_t nd = fMother->GetNdaughters();
   for (i=0; i<nd; i++) {
      node = fMother->GetNode(i);
      node->GetVolume()->SetVisibility(kFALSE);
   }
   fVolume->SetVisibility(kTRUE);
   for (i=0; i<fNovlp; i++) {
      node = fMother->GetNode(fOverlaps[i]);
      node->GetVolume()->SetVisibility(kTRUE);
   }
   gGeoManager->SetVisLevel(1);
   fMother->Draw();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill array with node id. Recursive on node branch.

void TGeoNode::FillIdArray(Int_t &ifree, Int_t &nodeid, Int_t *array) const
{
   Int_t nd = GetNdaughters();
   if (!nd) return;
   TGeoNode *daughter;
   Int_t istart = ifree; // start index for daughters
   ifree += nd;
   for (Int_t id=0; id<nd; id++) {
      daughter = GetDaughter(id);
      array[istart+id] = ifree;
      array[ifree++] = ++nodeid;
      daughter->FillIdArray(ifree, nodeid, array);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Search for a node within the branch of this one.

Int_t TGeoNode::FindNode(const TGeoNode *node, Int_t level)
{
   Int_t nd = GetNdaughters();
   if (!nd) return -1;
   TIter next(fVolume->GetNodes());
   TGeoNode *daughter;
   while ((daughter=(TGeoNode*)next())) {
      if (daughter==node) {
         gGeoManager->GetListOfNodes()->AddAt(daughter,level+1);
         return (level+1);
      }
   }
   next.Reset();
   Int_t new_level;
   while ((daughter=(TGeoNode*)next())) {
      new_level = daughter->FindNode(node, level+1);
      if (new_level>=0) {
         gGeoManager->GetListOfNodes()->AddAt(daughter, level+1);
         return new_level;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// save attributes for this node

void TGeoNode::SaveAttributes(std::ostream &out)
{
   if (IsVisStreamed()) return;
   SetVisStreamed(kTRUE);
   char quote='"';
   Bool_t voldef = kFALSE;
   if ((fVolume->IsVisTouched()) && (!fVolume->IsVisStreamed())) {
      fVolume->SetVisStreamed(kTRUE);
      out << "   vol = gGeoManager->GetVolume("<<quote<<fVolume->GetName()<<quote<<");"<<std::endl;
      voldef = kTRUE;
      if (!fVolume->IsVisDaughters())
         out << "   vol->SetVisDaughters(kFALSE);"<<std::endl;
      if (fVolume->IsVisible()) {
/*
         if (fVolume->GetLineColor() != gStyle->GetLineColor())
            out<<"   vol->SetLineColor("<<fVolume->GetLineColor()<<");"<<std::endl;
         if (fVolume->GetLineStyle() != gStyle->GetLineStyle())
            out<<"   vol->SetLineStyle("<<fVolume->GetLineStyle()<<");"<<std::endl;
         if (fVolume->GetLineWidth() != gStyle->GetLineWidth())
            out<<"   vol->SetLineWidth("<<fVolume->GetLineWidth()<<");"<<std::endl;
*/
      } else {
         out <<"   vol->SetVisibility(kFALSE);"<<std::endl;
      }
   }
   if (!IsVisDaughters()) return;
   Int_t nd = GetNdaughters();
   if (!nd) return;
   TGeoNode *node;
   for (Int_t i=0; i<nd; i++) {
      node = GetDaughter(i);
      if (node->IsVisStreamed()) continue;
      if (node->IsVisTouched()) {
         if (!voldef)
            out << "   vol = gGeoManager->GetVolume("<<quote<<fVolume->GetName()<<quote<<");"<<std::endl;
         out<<"   node = vol->GetNode("<<i<<");"<<std::endl;
         if (!node->IsVisDaughters()) {
            out<<"   node->VisibleDaughters(kFALSE);"<<std::endl;
            node->SetVisStreamed(kTRUE);
            continue;
         }
         if (!node->IsVisible())
            out<<"   node->SetVisibility(kFALSE);"<<std::endl;
      }
      node->SaveAttributes(out);
      node->SetVisStreamed(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Connect user-defined extension to the node. The node "grabs" a copy, so
/// the original object can be released by the producer. Release the previously
/// connected extension if any.
///
/// NOTE: This interface is intended for user extensions and is guaranteed not
/// to be used by TGeo

void TGeoNode::SetUserExtension(TGeoExtension *ext)
{
   if (fUserExtension) fUserExtension->Release();
   fUserExtension = 0;
   if (ext) fUserExtension = ext->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect framework defined extension to the node. The node "grabs" a copy,
/// so the original object can be released by the producer. Release the previously
/// connected extension if any.
///
/// NOTE: This interface is intended for the use by TGeo and the users should
///       NOT connect extensions using this method

void TGeoNode::SetFWExtension(TGeoExtension *ext)
{
   if (fFWExtension) fFWExtension->Release();
   fFWExtension = 0;
   if (ext) fFWExtension = ext->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Get a copy of the user extension pointer. The user must call Release() on
/// the copy pointer once this pointer is not needed anymore (equivalent to
/// delete() after calling new())

TGeoExtension *TGeoNode::GrabUserExtension() const
{
   if (fUserExtension) return fUserExtension->Grab();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a copy of the framework extension pointer. The user must call Release() on
/// the copy pointer once this pointer is not needed anymore (equivalent to
/// delete() after calling new())

TGeoExtension *TGeoNode::GrabFWExtension() const
{
   if (fFWExtension) return fFWExtension->Grab();
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
/// Check the overlab between the bounding box of the node overlaps with the one
/// the brother with index IOTHER.

Bool_t TGeoNode::MayOverlap(Int_t iother) const
{
   if (!fOverlaps) return kFALSE;
   for (Int_t i=0; i<fNovlp; i++) if (fOverlaps[i]==iother) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert the point coordinates from mother reference to local reference system

void TGeoNode::MasterToLocal(const Double_t *master, Double_t *local) const
{
   GetMatrix()->MasterToLocal(master, local);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a vector from mother reference to local reference system

void TGeoNode::MasterToLocalVect(const Double_t *master, Double_t *local) const
{
   GetMatrix()->MasterToLocalVect(master, local);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert the point coordinates from local reference system to mother reference

void TGeoNode::LocalToMaster(const Double_t *local, Double_t *master) const
{
   GetMatrix()->LocalToMaster(local, master);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a vector from local reference system to mother reference

void TGeoNode::LocalToMasterVect(const Double_t *local, Double_t *master) const
{
   GetMatrix()->LocalToMasterVect(local, master);
}

////////////////////////////////////////////////////////////////////////////////
/// Print the path (A/B/C/...) to this node on stdout

void TGeoNode::ls(Option_t * /*option*/) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this node and its content according to visualization settings.

void TGeoNode::Paint(Option_t *option)
{
   TVirtualGeoPainter *painter = gGeoManager->GetGeomPainter();
   if (!painter) return;
   painter->PaintNode(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// print daughters candidates for containing current point

void TGeoNode::PrintCandidates() const
{
   Double_t point[3];
   gGeoManager->MasterToLocal(gGeoManager->GetCurrentPoint(), &point[0]);
   printf("   Local : %g, %g, %g\n", point[0], point[1], point[2]);
   if (!fVolume->Contains(&point[0])) {
      printf("current point not inside this\n");
      return;
   }
   TGeoPatternFinder *finder = fVolume->GetFinder();
   TGeoNode *node;
   if (finder) {
      printf("current node divided\n");
      node = finder->FindNode(&point[0]);
      if (!node) {
         printf("point not inside division element\n");
         return;
      }
      printf("inside division element %s\n", node->GetName());
      return;
   }
   TGeoVoxelFinder *voxels = fVolume->GetVoxels();
   if (!voxels) {
      printf("volume not voxelized\n");
      return;
   }
   Int_t ncheck = 0;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   TGeoStateInfo &info = *nav->GetCache()->GetInfo();
   Int_t *check_list = voxels->GetCheckList(&point[0], ncheck, info);
   nav->GetCache()->ReleaseInfo();
   voxels->PrintVoxelLimits(&point[0]);
   if (!check_list) {
      printf("no candidates for current point\n");
      return;
   }
   TString overlap = "ONLY";
   for (Int_t id=0; id<ncheck; id++) {
      node = fVolume->GetNode(check_list[id]);
      if (node->IsOverlapping()) overlap = "MANY";
      else overlap = "ONLY";
      printf("%i %s %s\n", check_list[id], node->GetName(), overlap.Data());
   }
   PrintOverlaps();
}

////////////////////////////////////////////////////////////////////////////////
/// print possible overlapping nodes

void TGeoNode::PrintOverlaps() const
{
   if (!fOverlaps) {printf("node %s no overlaps\n", GetName()); return;}
   printf("Overlaps for node %s :\n", GetName());
   TGeoNode *node;
   for (Int_t i=0; i<fNovlp; i++) {
      node = fMother->GetNode(fOverlaps[i]);
      printf("   %s\n", node->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape

Double_t TGeoNode::Safety(const Double_t *point, Bool_t in) const
{
   Double_t local[3];
   GetMatrix()->MasterToLocal(point,local);
   return fVolume->GetShape()->Safety(local,in);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy content of lst of overlaps from source array

void TGeoNode::CopyOverlaps(Int_t *src, Int_t novlp)
{
   Int_t *ovlps = nullptr;
   if (src && (novlp > 0)) {
      ovlps = new Int_t[novlp];
      memcpy(ovlps, src, novlp*sizeof(Int_t));
   }
   SetOverlaps(ovlps, novlp);
}

////////////////////////////////////////////////////////////////////////////////
/// set the list of overlaps for this node (ovlp must be created with operator new)

void TGeoNode::SetOverlaps(Int_t *ovlp, Int_t novlp)
{
   if (fOverlaps) delete [] fOverlaps;
   fOverlaps = ovlp;
   fNovlp = novlp;
}

////////////////////////////////////////////////////////////////////////////////
/// Set visibility of the node (obsolete).

void TGeoNode::SetVisibility(Bool_t vis)
{
   if (gGeoManager->IsClosed()) SetVisTouched(kTRUE);
   TGeoAtt::SetVisibility(vis);
   if (vis && !fVolume->IsVisible()) fVolume->SetVisibility(vis);
   gGeoManager->ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Set visibility of the daughters (obsolete).

void TGeoNode::VisibleDaughters(Bool_t vis)
{
   if (gGeoManager->IsClosed()) SetVisTouched(kTRUE);
   SetVisDaughters(vis);
   gGeoManager->ModifiedPad();
}

/** \class TGeoNodeMatrix
\ingroup Geometry_classes
A node containing local transformation.
*/

ClassImp(TGeoNodeMatrix);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoNodeMatrix::TGeoNodeMatrix()
{
   fMatrix       = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoNodeMatrix::TGeoNodeMatrix(const TGeoVolume *vol, const TGeoMatrix *matrix) :
             TGeoNode(vol)
{
   fMatrix = (TGeoMatrix*)matrix;
   if (!fMatrix) fMatrix = gGeoIdentity;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoNodeMatrix::~TGeoNodeMatrix()
{
}

////////////////////////////////////////////////////////////////////////////////
/// return the total size in bytes of this node

Int_t TGeoNodeMatrix::GetByteCount() const
{
   Int_t count = 40 + 4; // TGeoNode + fMatrix
//   if (fMatrix) count += fMatrix->GetByteCount();
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns type of optimal voxelization for this node.
///  - type = 0 -> cartesian
///  - type = 1 -> cylindrical

Int_t TGeoNodeMatrix::GetOptimalVoxels() const
{
   Bool_t type = fVolume->GetShape()->IsCylType();
   if (!type) return 0;
   if (!fMatrix->IsRotAboutZ()) return 0;
   const Double_t *transl = fMatrix->GetTranslation();
   if (TMath::Abs(transl[0])>1E-10) return 0;
   if (TMath::Abs(transl[1])>1E-10) return 0;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this node.

TGeoNode *TGeoNodeMatrix::MakeCopyNode() const
{
   TGeoNodeMatrix *node = new TGeoNodeMatrix(fVolume, fMatrix);
   node->SetName(GetName());
   node->SetTitle(GetTitle());
   // set the mother
   node->SetMotherVolume(fMother);
   // set the copy number
   node->SetNumber(fNumber);
   // copy overlaps
   node->CopyOverlaps(fOverlaps, fNovlp);

   // copy VC
   if (IsVirtual()) node->SetVirtual();
   if (IsOverlapping()) node->SetOverlapping(); // <--- ADDED
   // Copy extensions
   node->SetUserExtension(fUserExtension);
   node->SetFWExtension(fFWExtension);
   node->SetCloned();
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Matrix setter.

void TGeoNodeMatrix::SetMatrix(const TGeoMatrix *matrix)
{
   fMatrix = (TGeoMatrix*)matrix;
   if (!fMatrix) fMatrix = gGeoIdentity;
}

/** \class TGeoNodeOffset
\ingroup Geometry_classes
Node containing an offset.
*/

ClassImp(TGeoNodeOffset);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoNodeOffset::TGeoNodeOffset()
{
   TObject::SetBit(kGeoNodeOffset);
   fOffset = 0;
   fIndex = 0;
   fFinder = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Null pointer to matrix means identity transformation

TGeoNodeOffset::TGeoNodeOffset(const TGeoVolume *vol, Int_t index, Double_t offset) :
           TGeoNode(vol)
{
   TObject::SetBit(kGeoNodeOffset);
   fOffset = offset;
   fIndex = index;
   fFinder = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoNodeOffset::~TGeoNodeOffset()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Get the index of this offset.

Int_t TGeoNodeOffset::GetIndex() const
{
   return (fIndex+fFinder->GetDivIndex());
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this node

TGeoNode *TGeoNodeOffset::MakeCopyNode() const
{
   TGeoNodeOffset *node = new TGeoNodeOffset(fVolume, GetIndex(), fOffset);
   node->SetName(GetName());
   node->SetTitle(GetTitle());
   // set the mother
   node->SetMotherVolume(fMother);
   // set the copy number
   node->SetNumber(fNumber);
   if (IsVirtual()) node->SetVirtual();
   // set the finder
   node->SetFinder(GetFinder());
   // set extensions
   node->SetUserExtension(fUserExtension);
   node->SetFWExtension(fFWExtension);
   return node;
}

/** \class TGeoIterator
\ingroup Geometry_classes
A geometry iterator.

A geometry iterator that sequentially follows all nodes of the geometrical
hierarchy of a volume. The iterator has to be initiated with a top volume
pointer:

~~~ {.cpp}
   TGeoIterator next(myVolume);
~~~

One can use the iterator as any other in ROOT:

~~~ {.cpp}
   TGeoNode *node;
   while ((node=next())) {
      ...
   }
~~~

The iterator can perform 2 types of iterations that can be selected via:

~~~ {.cpp}
   next.SetType(Int_t type);
~~~

Here TYPE can be:
  - 0 (default) - 'first daughter next' behavior
  - 1           - iteration at the current level only

Supposing the tree structure looks like:

~~~ {.cpp}
TOP ___ A_1 ___ A1_1 ___ A11_1
   |       |        |___ A12_1
   |      |_____A2_1 ___ A21_1
   |                |___ A21_2
   |___ B_1 ...
~~~

The order of iteration for TYPE=0 is: A_1, A1_1, A11_1, A12_1, A2_1, A21_1,
A21_2, B_1, ...

The order of iteration for TYPE=1 is: A_1, B_1, ...
At any moment during iteration, TYPE can be changed. If the last iterated node
is for instance A1_1 and the iteration type was 0, one can do:

~~~ {.cpp}
   next.SetType(1);
~~~

The next iterated nodes will be the rest of A daughters: A2,A3,... The iterator
will return 0 after finishing all daughters of A.

During iteration, the following can be retrieved:
  - Top volume where iteration started:    TGeoIterator::GetTopVolume()
  - Node at level I in the current branch: TGeoIterator::GetNode(Int_t i)
  - Iteration type:                        TGeoIterator::GetType()
  - Global matrix of the current node with respect to the top volume:
                                         TGeoIterator::GetCurrentMatrix()

The iterator can be reset by changing (or not) the top volume:

~~~ {.cpp}
   TGeoIterator::Reset(TGeoVolume *top);
~~~

### Example:

We want to find out a volume named "MyVol" in the hierarchy of TOP volume.

~~~ {.cpp}
   TIter next(TOP);
   TGeoNode *node;
   TString name("MyVol");
   while ((node=next()))
      if (name == node->GetVolume()->GetName()) return node->GetVolume();
~~~
*/

/** \class TGeoIteratorPlugin
\ingroup Geometry_classes
*/

ClassImp(TGeoIteratorPlugin);
ClassImp(TGeoIterator);

////////////////////////////////////////////////////////////////////////////////
/// Geometry iterator for a branch starting with a TOP node.

TGeoIterator::TGeoIterator(TGeoVolume *top)
{
   fTop = top;
   fLevel = 0;
   fMustResume = kFALSE;
   fMustStop = kFALSE;
   fType = 0;
   fArray = new Int_t[30];
   fMatrix = new TGeoHMatrix();
   fTopName = fTop->GetName();
   fPlugin = 0;
   fPluginAutoexec = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TGeoIterator::TGeoIterator(const TGeoIterator &iter)
{
   fTop = iter.GetTopVolume();
   fLevel = iter.GetLevel();
   fMustResume = kFALSE;
   fMustStop = kFALSE;
   fType = iter.GetType();
   fArray = new Int_t[30+ 30*Int_t(fLevel/30)];
   for (Int_t i=0; i<fLevel+1; i++) fArray[i] = iter.GetIndex(i);
   fMatrix = new TGeoHMatrix(*iter.GetCurrentMatrix());
   fTopName = fTop->GetName();
   fPlugin = iter.fPlugin;
   fPluginAutoexec = iter.fPluginAutoexec;;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoIterator::~TGeoIterator()
{
   if (fArray) delete [] fArray;
   delete fMatrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment.

TGeoIterator &TGeoIterator::operator=(const TGeoIterator &iter)
{
   if (&iter == this) return *this;
   fTop = iter.GetTopVolume();
   fLevel = iter.GetLevel();
   fMustResume = kFALSE;
   fMustStop = kFALSE;
   fType = iter.GetType();
   if (fArray) delete [] fArray;
   fArray = new Int_t[30+ 30*Int_t(fLevel/30)];
   for (Int_t i=0; i<fLevel+1; i++) fArray[i] = iter.GetIndex(i);
   if (!fMatrix) fMatrix = new TGeoHMatrix();
   *fMatrix = *iter.GetCurrentMatrix();
   fTopName = fTop->GetName();
   fPlugin = iter.fPlugin;
   fPluginAutoexec = iter.fPluginAutoexec;;
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns next node.

TGeoNode *TGeoIterator::Next()
{
   if (fMustStop) return 0;
   TGeoNode *mother = 0;
   TGeoNode *next = 0;
   Int_t i;
   Int_t nd = fTop->GetNdaughters();
   if (!nd) {
      fMustStop = kTRUE;
      return 0;
   }
   if (!fLevel) {
      fArray[++fLevel] = 0;
      next = fTop->GetNode(0);
      if (fPlugin && fPluginAutoexec) fPlugin->ProcessNode();
      return next;
   }
   next = fTop->GetNode(fArray[1]);
   // Move to current node
   for (i=2; i<fLevel+1; i++) {
      mother = next;
      next = mother->GetDaughter(fArray[i]);
   }
   if (fMustResume) {
      fMustResume = kFALSE;
      if (fPlugin && fPluginAutoexec) fPlugin->ProcessNode();
      return next;
   }

   switch (fType) {
      case 0:  // default next daughter behavior
         nd = next->GetNdaughters();
         if (nd) {
            // First daughter next
            fLevel++;
            if ((fLevel%30)==0) IncreaseArray();
            fArray[fLevel] = 0;
            if (fPlugin && fPluginAutoexec) fPlugin->ProcessNode();
            return next->GetDaughter(0);
         }
         // cd up and pick next
         while (next) {
            next = GetNode(fLevel-1);
            if (!next) {
               nd = fTop->GetNdaughters();
               if (fArray[fLevel]<nd-1) {
                  fArray[fLevel]++;
                  if (fPlugin && fPluginAutoexec) fPlugin->ProcessNode();
                  return fTop->GetNode(fArray[fLevel]);
               }
               fMustStop = kTRUE;
               return 0;
            } else {
               nd = next->GetNdaughters();
               if (fArray[fLevel]<nd-1) {
                  fArray[fLevel]++;
                  if (fPlugin && fPluginAutoexec) fPlugin->ProcessNode();
                  return next->GetDaughter(fArray[fLevel]);
               }
            }
            fLevel--;
         }
         break;
      case 1:  // one level search
         if (mother) nd = mother->GetNdaughters();
         if (fArray[fLevel]<nd-1) {
            fArray[fLevel]++;
            if (fPlugin && fPluginAutoexec) fPlugin->ProcessNode();
            if (!mother) return fTop->GetNode(fArray[fLevel]);
            else return mother->GetDaughter(fArray[fLevel]);
         }
   }
   fMustStop = kTRUE;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns next node.

TGeoNode *TGeoIterator::operator()()
{
   return Next();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns global matrix for current node.

const TGeoMatrix *TGeoIterator::GetCurrentMatrix() const
{
   fMatrix->Clear();
   if (!fLevel) return fMatrix;
   TGeoNode *node = fTop->GetNode(fArray[1]);
   fMatrix->Multiply(node->GetMatrix());
   for (Int_t i=2; i<fLevel+1; i++) {
      node = node->GetDaughter(fArray[i]);
      fMatrix->Multiply(node->GetMatrix());
   }
   return fMatrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current node at a given level.

TGeoNode *TGeoIterator::GetNode(Int_t level) const
{
   if (!level || level>fLevel) return 0;
   TGeoNode *node = fTop->GetNode(fArray[1]);
   for (Int_t i=2; i<level+1; i++) node = node->GetDaughter(fArray[i]);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the path for the current node.

void TGeoIterator::GetPath(TString &path) const
{
   path = fTopName;
   if (!fLevel) return;
   TGeoNode *node = fTop->GetNode(fArray[1]);
   path += "/";
   path += node->GetName();
   for (Int_t i=2; i<fLevel+1; i++) {
      node = node->GetDaughter(fArray[i]);
      path += "/";
      path += node->GetName();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Increase by 30 the size of the array.

void TGeoIterator::IncreaseArray()
{
   Int_t *array = new Int_t[fLevel+30];
   memcpy(array, fArray, fLevel*sizeof(Int_t));
   delete [] fArray;
   fArray = array;
}

////////////////////////////////////////////////////////////////////////////////
/// Resets the iterator for volume TOP.

void TGeoIterator::Reset(TGeoVolume *top)
{
   if (top) fTop = top;
   fLevel = 0;
   fMustResume = kFALSE;
   fMustStop = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the top name for path

void TGeoIterator::SetTopName(const char *name)
{
   fTopName = name;
}

////////////////////////////////////////////////////////////////////////////////
/// Stop iterating the current branch. The iteration of the next node will
/// behave as if the branch starting from the current node (included) is not existing.

void TGeoIterator::Skip()
{
   fMustResume = kTRUE;
   TGeoNode *next = GetNode(fLevel);
   if (!next) return;
   Int_t nd;
   switch (fType) {
      case 0:  // default next daughter behavior
         // cd up and pick next
         while (next) {
            next = GetNode(fLevel-1);
            nd = (next==0)?fTop->GetNdaughters():next->GetNdaughters();
            if (fArray[fLevel]<nd-1) {
               ++fArray[fLevel];
               return;
            }
            fLevel--;
            if (!fLevel) {
               fMustStop = kTRUE;
               return;
            }
         }
         break;
      case 1:  // one level search
         next = GetNode(fLevel-1);
         nd = (next==0)?fTop->GetNdaughters():next->GetNdaughters();
         if (fArray[fLevel]<nd-1) {
            ++fArray[fLevel];
            return;
         }
         fMustStop = kTRUE;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set a plugin.

void TGeoIterator::SetUserPlugin(TGeoIteratorPlugin *plugin)
{
   fPlugin = plugin;
   if (plugin) plugin->SetIterator(this);
}
