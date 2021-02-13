// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveGeoNode.h"
#include "TEveTrans.h"
#include "TEveManager.h"
#include "TEvePolygonSetProjected.h"

#include "TEveGeoShape.h"
#include "TEveGeoShapeExtract.h"
#include "TEvePad.h"
#include "TEveGeoPolyShape.h"
#include "TGLScenePad.h"
#include "TGLFaceSet.h"

#include "TROOT.h"
#include "TBuffer3D.h"
#include "TVirtualViewer3D.h"
#include "TColor.h"
#include "TFile.h"

#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoShapeAssembly.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TVirtualGeoPainter.h"

/** \class TEveGeoNode
\ingroup TEve
Wrapper for TGeoNode that allows it to be shown in GUI and controlled as a TEveElement.
*/

ClassImp(TEveGeoNode);

Int_t                 TEveGeoNode::fgCSGExportNSeg = 64;
std::list<TGeoShape*> TEveGeoNode::fgTemporaryStore;

////////////////////////////////////////////////////////////////////////////////
/// Returns number of segments used for CSG export.

Int_t TEveGeoNode::GetCSGExportNSeg()
{
   return fgCSGExportNSeg;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets number of segments used for CSG export.

void TEveGeoNode::SetCSGExportNSeg(Int_t nseg)
{
   fgCSGExportNSeg = nseg;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGeoNode::TEveGeoNode(TGeoNode* node) :
   TEveElement(),
   TObject(),
   fNode(node)
{
   // Hack!! Should use cint to retrieve TAttLine::fLineColor offset.
   char* l = (char*) dynamic_cast<TAttLine*>(node->GetVolume());
   SetMainColorPtr((Color_t*)(l + sizeof(void*)));
   SetMainTransparency(fNode->GetVolume()->GetTransparency());

   SetRnrSelfChildren(fNode->IsVisible(), fNode->IsVisDaughters());
}

////////////////////////////////////////////////////////////////////////////////
/// Return name, taken from geo-node. Used via TObject.

const char* TEveGeoNode::GetName()  const
{
   return fNode->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return title, taken from geo-node. Used via TObject.

const char* TEveGeoNode::GetTitle() const
{
   return fNode->GetTitle();
}

////////////////////////////////////////////////////////////////////////////////
/// Return name, taken from geo-node. Used via TEveElement.

const char* TEveGeoNode::GetElementName()  const
{
   return fNode->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return title, taken from geo-node. Used via TEveElement.

const char* TEveGeoNode::GetElementTitle() const
{
   return fNode->GetTitle();
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if child-nodes have been imported ... imports them if not.
/// Then calls TEveElement::ExpandIntoListTree.

void TEveGeoNode::ExpandIntoListTree(TGListTree* ltree,
                                     TGListTreeItem* parent)
{
   if ( ! HasChildren() && fNode->GetVolume()->GetNdaughters() > 0)
   {
      TIter next(fNode->GetVolume()->GetNodes());
      TGeoNode* dnode;
      while ((dnode = (TGeoNode*) next()) != 0)
      {
         TEveGeoNode* node_re = new TEveGeoNode(dnode);
         AddElement(node_re);
      }
   }
   TEveElement::ExpandIntoListTree(ltree, parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Expand children into all list-trees.

void TEveGeoNode::ExpandIntoListTrees()
{
   for (sLTI_i i = fItems.begin(); i != fItems.end(); ++i)
   {
      ExpandIntoListTree(i->fTree, i->fItem);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand children into all list-trees recursively.
/// This is useful if one wants to export extracted shapes.

void TEveGeoNode::ExpandIntoListTreesRecursively()
{
   ExpandIntoListTrees();
   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
   {
      TEveGeoNode *egn = dynamic_cast<TEveGeoNode*>(*i);
      if (egn)
         egn->ExpandIntoListTreesRecursively();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Override from TEveElement.
/// Process visibility changes and forward them to fNode.

void TEveGeoNode::AddStamp(UChar_t bits)
{
   TEveElement::AddStamp(bits);
   if (bits & kCBVisibility)
   {
      fNode->SetVisibility(fRnrSelf);
      fNode->VisibleDaughters(fRnrChildren);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Can edit main-color -- not available for assemblies.

Bool_t TEveGeoNode::CanEditMainColor() const
{
   return ! fNode->GetVolume()->IsAssembly();
}

////////////////////////////////////////////////////////////////////////////////
/// Set color, propagate to volume's line color.

void TEveGeoNode::SetMainColor(Color_t color)
{
   TEveElement::SetMainColor(color);
   fNode->GetVolume()->SetLineColor(color);
}

////////////////////////////////////////////////////////////////////////////////
/// Can edit main transparency -- not available for assemblies.

Bool_t TEveGeoNode::CanEditMainTransparency() const
{
   return ! fNode->GetVolume()->IsAssembly();
}

////////////////////////////////////////////////////////////////////////////////
/// Get transparency -- it is taken from the geo node.

Char_t TEveGeoNode::GetMainTransparency() const
{
   return fNode->GetVolume()->GetTransparency();
}

////////////////////////////////////////////////////////////////////////////////
/// Set transparency, propagate to volume's transparency.

void TEveGeoNode::SetMainTransparency(Char_t t)
{
   TEveElement::SetMainTransparency(t);
   fNode->GetVolume()->SetTransparency(t);
}

////////////////////////////////////////////////////////////////////////////////
/// Updates all reve-browsers having the node in their contents.
/// All 3D-pads updated if any change found.
///
/// Should (could?) be optimized with some assumptions about
/// volume/node structure (search for parent, know the same node can not
/// reoccur on lower level once found).

void TEveGeoNode::UpdateNode(TGeoNode* node)
{
   static const TEveException eH("TEveGeoNode::UpdateNode ");

   // printf("%s node %s %p\n", eH.Data(), node->GetName(), node);

   if (fNode == node)
      StampColorSelection();

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      ((TEveGeoNode*)(*i))->UpdateNode(node);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Updates all reve-browsers having the volume in their contents.
/// All 3D-pads updated if any change found.
///
/// Should (could?) be optimized with some assumptions about
/// volume/node structure (search for parent, know the same node can not
/// reoccur on lower level once found).

void TEveGeoNode::UpdateVolume(TGeoVolume* volume)
{
   static const TEveException eH("TEveGeoNode::UpdateVolume ");

   // printf("%s volume %s %p\n", eH.Data(), volume->GetName(), volume);

   if(fNode->GetVolume() == volume)
      StampColorSelection();

   for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      ((TEveGeoNode*)(*i))->UpdateVolume(volume);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the object.

void TEveGeoNode::Draw(Option_t* option)
{
   TString opt("SAME");
   opt += option;
   fNode->GetVolume()->Draw(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Save TEveGeoShapeExtract tree starting at this node.
/// This function is obsolete, use SaveExtract() instead.

void TEveGeoNode::Save(const char* file, const char* name, Bool_t leafs_only)
{
   Warning("Save()", "This function is deprecated, use SaveExtract() instead.");
   SaveExtract(file, name, leafs_only);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the shape tree as TEveGeoShapeExtract.
/// File is always recreated.

void TEveGeoNode::SaveExtract(const char* file, const char* name, Bool_t leafs_only)
{
   TEveGeoShapeExtract* gse = DumpShapeTree(this, 0, leafs_only);
   if (gse)
   {
      TFile f(file, "RECREATE");
      gse->Write(name);
      f.Close();
   }

   for (std::list<TGeoShape*>::iterator i = fgTemporaryStore.begin(); i != fgTemporaryStore.end(); ++i)
      delete *i;
   fgTemporaryStore.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Write the shape tree as TEveGeoShapeExtract to current directory.

void TEveGeoNode::WriteExtract(const char* name, Bool_t leafs_only)
{
   TEveGeoShapeExtract* gse = DumpShapeTree(this, 0, leafs_only);
   if (gse)
   {
      gse->Write(name);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Export the node hierarchy into tree of TEveGeoShapeExtract objects.

TEveGeoShapeExtract* TEveGeoNode::DumpShapeTree(TEveGeoNode*         geon,
                                                TEveGeoShapeExtract* parent,
                                                Bool_t               leafs_only)
{
   static const TEveException eh("TEveGeoNode::DumpShapeTree ");

   TGeoNode*   tnode   = 0;
   TGeoVolume* tvolume = 0;
   TGeoShape*  tshape  = 0;

   tnode = geon->GetNode();
   if (tnode == 0)
   {
      Info(eh, "Null TGeoNode for TEveGeoNode '%s': assuming it's a holder and descending.", geon->GetName());
   }
   else
   {
      tvolume = tnode->GetVolume();
      if (tvolume == 0) {
         Warning(eh, "Null TGeoVolume for TEveGeoNode '%s'; skipping its sub-tree.\n", geon->GetName());
         return 0;
      }
      tshape  = tvolume->GetShape();
      if (tshape->IsComposite())
      {
         TEvePad pad;
         TEvePadHolder gpad(kFALSE, &pad);
         pad.GetListOfPrimitives()->Add(tshape);
         TGLScenePad scene_pad(&pad);
         pad.SetViewer3D(&scene_pad);

         {
            TEveGeoManagerHolder gmgr(tvolume->GetGeoManager(), fgCSGExportNSeg);
            gGeoManager->SetPaintVolume(tvolume);

            TGeoMatrix *gst = TGeoShape::GetTransform();
            TGeoShape::SetTransform(TEveGeoShape::GetGeoHMatrixIdentity());

            scene_pad.BeginScene();
            dynamic_cast<TGeoCompositeShape*>(tshape)->PaintComposite();
            scene_pad.EndScene();

            TGeoShape::SetTransform(gst);
         }

         pad.SetViewer3D(0);

         TGLFaceSet* fs = dynamic_cast<TGLFaceSet*>(scene_pad.FindLogical(tvolume));
         if (!fs) {
            Warning(eh, "Failed extracting CSG tesselation TEveGeoNode '%s'; skipping its sub-tree.\n", geon->GetName());
            return 0;
         }

         TEveGeoPolyShape* egps = new TEveGeoPolyShape();
         egps->SetFromFaceSet(fs);
         tshape = egps;
         fgTemporaryStore.push_back(egps);
      }
   }

   // transformation
   TEveTrans trans;
   if (parent)
      trans.SetFromArray(parent->GetTrans());
   if (tnode)
   {
      TGeoMatrix     *gm = tnode->GetMatrix();
      const Double_t *rm = gm->GetRotationMatrix();
      const Double_t *tv = gm->GetTranslation();
      TEveTrans t;
      t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
      t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
      t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
      t(1,4) = tv[0]; t(2,4) = tv[1]; t(3,4) = tv[2];
      trans *= t;
   }

   TEveGeoShapeExtract* gse = new TEveGeoShapeExtract(geon->GetName(), geon->GetTitle());
   gse->SetTrans(trans.Array());
   Int_t  ci = 0;
   Char_t transp = 0;
   if (tvolume) {
      ci = tvolume->GetLineColor();
      transp = tvolume->GetTransparency();
   }
   TColor* c = gROOT->GetColor(ci);
   Float_t rgba[4] = {1, 0, 0, 1.0f - transp/100.0f};
   if (c) {
      rgba[0] = c->GetRed();
      rgba[1] = c->GetGreen();
      rgba[2] = c->GetBlue();
   }
   gse->SetRGBA(rgba);
   rgba[3] = 1;
   c = gROOT->GetColor(TColor::GetColorDark(ci));
   if (c) {
      rgba[0] = c->GetRed();
      rgba[1] = c->GetGreen();
      rgba[2] = c->GetBlue();
   }
   gse->SetRGBALine(rgba);

   // Keep default extract line color --> black.
   Bool_t rnr     = tnode ? tnode->IsVisible()      : geon->GetRnrSelf();
   Bool_t rnr_els = tnode ? tnode->IsVisDaughters() : geon->GetRnrChildren();
   if (tvolume) {
      rnr     = rnr     && tvolume->IsVisible();
      rnr_els = rnr_els && tvolume->IsVisDaughters();
   }
   gse->SetRnrSelf    (rnr);
   gse->SetRnrElements(rnr_els);
   gse->SetRnrFrame   (kTRUE);
   gse->SetMiniFrame  (kTRUE);

   gse->SetShape((leafs_only && geon->HasChildren()) ? 0 : tshape);

   if (geon->HasChildren())
   {
      TList* ele = new TList();
      gse->SetElements(ele);
      gse->GetElements()->SetOwner(true);

      TEveElement::List_i i = geon->BeginChildren();
      while (i != geon->EndChildren())
      {
         TEveGeoNode* l = dynamic_cast<TEveGeoNode*>(*i);
         DumpShapeTree(l, gse, leafs_only);
         ++i;
      }
   }

   if (parent)
      parent->GetElements()->Add(gse);

   return gse;
}



/** \class TEveGeoTopNode
\ingroup TEve
A wrapper over a TGeoNode, possibly displaced with a global
trasformation stored in TEveElement.

It holds a pointer to TGeoManager and controls for steering of
TGeoPainter, fVisOption, fVisLevel and fMaxVisNodes. They have the
same meaning as in TGeoManager/TGeoPainter.
*/

ClassImp(TEveGeoTopNode);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGeoTopNode::TEveGeoTopNode(TGeoManager* manager, TGeoNode* node,
                               Int_t visopt, Int_t vislvl, Int_t maxvisnds) :
   TEveGeoNode  (node),
   fManager     (manager),
   fVisOption   (visopt),
   fVisLevel    (vislvl),
   fMaxVisNodes (maxvisnds)
{
   InitMainTrans();
   fRnrSelf = kTRUE; // Override back from TEveGeoNode.
}

////////////////////////////////////////////////////////////////////////////////
/// Use transformation matrix from the TGeoNode.
/// Warning: this is local transformation of the node!

void TEveGeoTopNode::UseNodeTrans()
{
   RefMainTrans().SetFrom(*fNode->GetMatrix());
}

////////////////////////////////////////////////////////////////////////////////
/// Revert from TEveGeoNode back to standard behaviour, that is,
/// do not pass visibility changes to fNode as they are honoured
/// in Paint() method.

void TEveGeoTopNode::AddStamp(UChar_t bits)
{
   TEveElement::AddStamp(bits);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the top-node.

void TEveGeoTopNode::Draw(Option_t* option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the enclosed TGeo hierarchy with visibility level and
/// option given in data-members.
/// Uses TGeoPainter internally.

void TEveGeoTopNode::Paint(Option_t* option)
{
   if (fRnrSelf)
   {
      TEveGeoManagerHolder geo_holder(fManager);
      TVirtualPad *pad = gPad;
      gPad = 0;
      TGeoVolume* top_volume = fManager->GetTopVolume();
      if (fVisLevel > 0)
         fManager->SetVisLevel(fVisLevel);
      else
         fManager->SetMaxVisNodes(fMaxVisNodes);
      TVirtualGeoPainter* vgp = fManager->GetGeomPainter();
      fManager->SetTopVolume(fNode->GetVolume());
      switch (fVisOption)
      {
         case 0:
            fNode->GetVolume()->SetVisContainers(kTRUE);
            fManager->SetTopVisible(kTRUE);
            break;
         case 1:
            fNode->GetVolume()->SetVisLeaves(kTRUE);
            fManager->SetTopVisible(kFALSE);
            break;
         case 2:
            fNode->GetVolume()->SetVisOnly(kTRUE);
            break;
      }
      gPad = pad;
      if(vgp != 0) {
         vgp->SetVisOption(fVisOption);
         TGeoHMatrix geomat;
         if (HasMainTrans()) RefMainTrans().SetGeoHMatrix(geomat);
         vgp->PaintNode(fNode, option, &geomat);
      }
      fManager->SetTopVolume(top_volume);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Callback for propagating volume visibility changes.

void TEveGeoTopNode::VolumeVisChanged(TGeoVolume* volume)
{
   static const TEveException eh("TEveGeoTopNode::VolumeVisChanged ");
   printf("%s volume %s %p\n", eh.Data(), volume->GetName(), (void*)volume);
   UpdateVolume(volume);
}

////////////////////////////////////////////////////////////////////////////////
/// Callback for propagating volume parameter changes.

void TEveGeoTopNode::VolumeColChanged(TGeoVolume* volume)
{
   static const TEveException eh("TEveGeoTopNode::VolumeColChanged ");
   printf("%s volume %s %p\n", eh.Data(), volume->GetName(), (void*)volume);
   UpdateVolume(volume);
}

////////////////////////////////////////////////////////////////////////////////
/// Callback for propagating node visibility changes.

void TEveGeoTopNode::NodeVisChanged(TGeoNode* node)
{
   static const TEveException eh("TEveGeoTopNode::NodeVisChanged ");
   printf("%s node %s %p\n", eh.Data(), node->GetName(), (void*)node);
   UpdateNode(node);
}
