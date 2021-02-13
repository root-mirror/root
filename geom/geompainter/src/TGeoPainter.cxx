// @(#)root/geompainter:$Id: 58726ead32989b65bb2cbff2af4235fe9c6b12ae $
// Author: Andrei Gheata   05/03/02
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoPainter
\ingroup Geometry_classes

Class implementing all draw interfaces for a generic 3D viewer
using TBuffer3D mechanism.
*/

#include <map>
#include "TROOT.h"
#include "TClass.h"
#include "TColor.h"
#include "TPoint.h"
#include "TView.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TCanvasImp.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPluginManager.h"
#include "TVirtualPadEditor.h"
#include "TStopwatch.h"

#include "TPolyMarker3D.h"

#include "TGeoAtt.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TGeoTrack.h"
#include "TGeoOverlap.h"
#include "TGeoChecker.h"
#include "TGeoPhysicalNode.h"
#include "TGeoPolygon.h"
#include "TGeoCompositeShape.h"
#include "TGeoShapeAssembly.h"
#include "TGeoPainter.h"
#include "TMath.h"

#include "X3DBuffer.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualViewer3D.h"
#include "TVirtualX.h"

ClassImp(TGeoPainter);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoPainter::TGeoPainter(TGeoManager *manager) : TVirtualGeoPainter(manager)
{
   TVirtualGeoPainter::SetPainter(this);
   if (manager) fGeoManager = manager;
   else {
      Error("ctor", "No geometry loaded");
      return;
   }
   fNsegments = fGeoManager->GetNsegments();
   fNVisNodes = 0;
   fBombX = 1.3;
   fBombY = 1.3;
   fBombZ = 1.3;
   fBombR = 1.3;
   fVisLevel = fGeoManager->GetVisLevel();
   fVisOption = fGeoManager->GetVisOption();
   fExplodedView = fGeoManager->GetBombMode();
   fVisBranch = "";
   fVolInfo = "";
   fVisLock = kFALSE;
   fIsRaytracing = kFALSE;
   fTopVisible = kFALSE;
   fPaintingOverlaps = kFALSE;
   fPlugin = nullptr;
   fVisVolumes = new TObjArray();
   fOverlap = nullptr;
   fGlobal = new TGeoHMatrix();
   fBuffer = new TBuffer3D(TBuffer3DTypes::kGeneric,20,3*20,0,0,0,0);
   fClippingShape = nullptr;
   fLastVolume = nullptr;
   fTopVolume = nullptr;
   fIsPaintingShape = kFALSE;
   memset(&fCheckedBox[0], 0, 6*sizeof(Double_t));

   fCheckedNode = fGeoManager->GetTopNode();
   fChecker = new TGeoChecker(fGeoManager);
   fIsEditable = kFALSE;
   DefineColors();
}
////////////////////////////////////////////////////////////////////////////////
/// Default destructor.

TGeoPainter::~TGeoPainter()
{
   if (fChecker) delete fChecker;
   delete fVisVolumes;
   delete fGlobal;
   delete fBuffer;
   if (fPlugin) delete fPlugin;
}
////////////////////////////////////////////////////////////////////////////////
/// Add numpoints, numsegs, numpolys to the global 3D size.

void TGeoPainter::AddSize3D(Int_t numpoints, Int_t numsegs, Int_t numpolys)
{
   gSize3D.numPoints += numpoints;
   gSize3D.numSegs   += numsegs;
   gSize3D.numPolys  += numpolys;
}
////////////////////////////////////////////////////////////////////////////////
/// Create a primary TGeoTrack.

TVirtualGeoTrack *TGeoPainter::AddTrack(Int_t id, Int_t pdgcode, TObject *particle)
{
   return (TVirtualGeoTrack*)(new TGeoTrack(id,pdgcode,0,particle));
}

////////////////////////////////////////////////////////////////////////////////
/// Average center of view of all painted tracklets and compute view box.

void TGeoPainter::AddTrackPoint(Double_t *point, Double_t *box, Bool_t reset)
{
   static Int_t npoints = 0;
   static Double_t xmin[3] = {0,0,0};
   static Double_t xmax[3] = {0,0,0};
   Int_t i;
   if (reset) {
      memset(box, 0, 6*sizeof(Double_t));
      memset(xmin, 0, 3*sizeof(Double_t));
      memset(xmax, 0, 3*sizeof(Double_t));
      npoints = 0;
      return;
   }
   if (npoints==0) {
      for (i=0; i<3; i++) xmin[i]=xmax[i]=0;
      npoints++;
   }
   npoints++;
   Double_t  ninv = 1./Double_t(npoints);
   for (i=0; i<3; i++) {
      box[i] += ninv*(point[i]-box[i]);
      if (point[i]<xmin[i]) xmin[i]=point[i];
      if (point[i]>xmax[i]) xmax[i]=point[i];
      box[i+3] = 0.5*(xmax[i]-xmin[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the new 'bombed' translation vector according current exploded view mode.

void TGeoPainter::BombTranslation(const Double_t *tr, Double_t *bombtr)
{
   memcpy(bombtr, tr, 3*sizeof(Double_t));
   switch (fExplodedView) {
      case kGeoNoBomb:
         return;
      case kGeoBombXYZ:
         bombtr[0] *= fBombX;
         bombtr[1] *= fBombY;
         bombtr[2] *= fBombZ;
         return;
      case kGeoBombCyl:
         bombtr[0] *= fBombR;
         bombtr[1] *= fBombR;
         bombtr[2] *= fBombZ;
         return;
      case kGeoBombSph:
         bombtr[0] *= fBombR;
         bombtr[1] *= fBombR;
         bombtr[2] *= fBombR;
         return;
      default:
         return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check pushes and pulls needed to cross the next boundary with respect to the
/// position given by FindNextBoundary. If radius is not mentioned the full bounding
/// box will be sampled.

void TGeoPainter::CheckBoundaryErrors(Int_t ntracks, Double_t radius)
{
   fChecker->CheckBoundaryErrors(ntracks, radius);
}

////////////////////////////////////////////////////////////////////////////////
/// Check the boundary errors reference file created by CheckBoundaryErrors method.
/// The shape for which the crossing failed is drawn with the starting point in red
/// and the extrapolated point to boundary (+/- failing push/pull) in yellow.

void TGeoPainter::CheckBoundaryReference(Int_t icheck)
{
   fChecker->CheckBoundaryReference(icheck);
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry checking method (see: TGeoManager::CheckGeometry())

void TGeoPainter::CheckGeometryFull(Bool_t checkoverlaps, Bool_t checkcrossings, Int_t ntracks, const Double_t *vertex)
{
   fChecker->CheckGeometryFull(checkoverlaps,checkcrossings,ntracks,vertex);
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry checking method (see TGeoChecker).

void TGeoPainter::CheckGeometry(Int_t nrays, Double_t startx, Double_t starty, Double_t startz) const
{
   fChecker->CheckGeometry(nrays, startx, starty, startz);
}

////////////////////////////////////////////////////////////////////////////////
/// Check overlaps for the top volume of the geometry, within a limit OVLP.

void TGeoPainter::CheckOverlaps(const TGeoVolume *vol, Double_t ovlp, Option_t *option) const
{
   fChecker->CheckOverlaps(vol, ovlp, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Check current point in the geometry.

void TGeoPainter::CheckPoint(Double_t x, Double_t y, Double_t z, Option_t *option)
{
   fChecker->CheckPoint(x,y,z,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Test for shape navigation methods. Summary for test numbers:
///  - 1: DistFromInside/Outside. Sample points inside the shape. Generate
///    directions randomly in cos(theta). Compute DistFromInside and move the
///    point with bigger distance. Compute DistFromOutside back from new point.
///    Plot d-(d1+d2)

void TGeoPainter::CheckShape(TGeoShape *shape, Int_t testNo, Int_t nsamples, Option_t *option)
{
   fChecker->CheckShape(shape, testNo, nsamples, option);
}

////////////////////////////////////////////////////////////////////////////////
///Clear the list of visible volumes
///reset the kVisOnScreen bit for volumes previously in the list

void TGeoPainter::ClearVisibleVolumes()
{
   if (!fVisVolumes) return;
   TIter next(fVisVolumes);
   TGeoVolume *vol;
   while ((vol = (TGeoVolume*)next())) {
      vol->ResetAttBit(TGeoAtt::kVisOnScreen);
   }
   fVisVolumes->Clear();
}


////////////////////////////////////////////////////////////////////////////////
/// Define 100 colors with increasing light intensities for each basic color (1-7)
/// Register these colors at indexes starting with 1000.

void TGeoPainter::DefineColors() const
{
   static Int_t color = 0;
   if (!color) {
      TColor::InitializeColors();
      for (auto icol=1; icol<10; ++icol)
         color = GetColor(icol, 0.5);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get index of a base color with given light intensity (0,1)

Int_t TGeoPainter::GetColor(Int_t base, Float_t light) const
{
   using IntMap_t = std::map<Int_t, Int_t>;
   constexpr Int_t ncolors = 100;
   constexpr Float_t lmin = 0.25;
   constexpr Float_t lmax = 0.75;
   static IntMap_t colmap;
   Int_t color = base;
   // Search color in the map
   auto it = colmap.find(base);
   if (it != colmap.end()) return (it->second + light*(ncolors-1));
   // Get color pointer if stored
   TColor* col_base = gROOT->GetColor(base);
   if (!col_base) {
      // If color not defined, use gray palette
      it = colmap.find(kBlack);
      if (it != colmap.end()) return (it->second + light*(ncolors-1));
      col_base = gROOT->GetColor(kBlack);
      color = 1;
   }
   // Create a color palette for col_base
   Float_t r=0., g=0., b=0., h=0., l=0., s=0.;
   Double_t red[2], green[2], blue[2];
   Double_t stop[] = {0., 1.0};

   if (col_base) col_base->GetRGB(r,g,b);
   TColor::RGB2HLS(r,g,b,h,l,s);
   TColor::HLS2RGB(h,lmin,s,r,g,b);
   red[0] = r;
   green[0] = g;
   blue[0] = b;
   TColor::HLS2RGB(h,lmax,s,r,g,b);
   red[1] = r;
   green[1] = g;
   blue[1] = b;
   Int_t color_map_idx = TColor::CreateGradientColorTable(2, stop, red, green, blue, ncolors);
   colmap[color] = color_map_idx;
   return (color_map_idx + light*(ncolors-1));
}

////////////////////////////////////////////////////////////////////////////////
/// Get currently drawn volume.

TGeoVolume *TGeoPainter::GetDrawnVolume() const
{
   if (!gPad) return nullptr;
   return fTopVolume;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the closest distance of approach from point px,py to a volume.

Int_t TGeoPainter::DistanceToPrimitiveVol(TGeoVolume *volume, Int_t px, Int_t py)
{
   const Int_t big = 9999;
   const Int_t inaxis = 7;
   const Int_t maxdist = 5;

   if (fTopVolume != volume) fTopVolume = volume;
   TView *view = gPad->GetView();
   if (!view) return big;
   TGeoBBox *box;
   fGlobal->Clear();
   TGeoShape::SetTransform(fGlobal);

   Int_t puxmin = gPad->XtoAbsPixel(gPad->GetUxmin());
   Int_t puymin = gPad->YtoAbsPixel(gPad->GetUymin());
   Int_t puxmax = gPad->XtoAbsPixel(gPad->GetUxmax());
   Int_t puymax = gPad->YtoAbsPixel(gPad->GetUymax());
   // return if point not in user area
   if (px < puxmin - inaxis) return big;
   if (py > puymin + inaxis) return big;
   if (px > puxmax + inaxis) return big;
   if (py < puymax - inaxis) return big;

   fCheckedNode = fGeoManager->GetTopNode();
   gPad->SetSelected(view);
   Int_t dist = big;
//   Int_t id;

   if (fPaintingOverlaps) {
      TGeoVolume *crt;
      crt = fOverlap->GetFirstVolume();
      *fGlobal = fOverlap->GetFirstMatrix();
      dist = crt->GetShape()->DistancetoPrimitive(px,py);
      if (dist<maxdist) {
         gPad->SetSelected(crt);
         box = (TGeoBBox*)crt->GetShape();
         fGlobal->LocalToMaster(box->GetOrigin(), &fCheckedBox[0]);
         fCheckedBox[3] = box->GetDX();
         fCheckedBox[4] = box->GetDY();
         fCheckedBox[5] = box->GetDZ();
         return 0;
      }
      crt = fOverlap->GetSecondVolume();
      *fGlobal = fOverlap->GetSecondMatrix();
      dist = crt->GetShape()->DistancetoPrimitive(px,py);
      if (dist<maxdist) {
         gPad->SetSelected(crt);
         box = (TGeoBBox*)crt->GetShape();
         fGlobal->LocalToMaster(box->GetOrigin(), &fCheckedBox[0]);
         fCheckedBox[3] = box->GetDX();
         fCheckedBox[4] = box->GetDY();
         fCheckedBox[5] = box->GetDZ();
         return 0;
      }
      return big;
   }
      // Compute distance to the right edge
   if ((puxmax+inaxis-px) < 40) {
      if ((py-puymax+inaxis) < 40) {
         // when the mouse points to the (40x40) right corner of the pad, the manager class is selected
         gPad->SetSelected(fGeoManager);
         fVolInfo = fGeoManager->GetName();
         box = (TGeoBBox*)volume->GetShape();
         memcpy(fCheckedBox, box->GetOrigin(), 3*sizeof(Double_t));
         fCheckedBox[3] = box->GetDX();
         fCheckedBox[4] = box->GetDY();
         fCheckedBox[5] = box->GetDZ();
         return 0;
      }
      // when the mouse points to the (40 pix) right edge of the pad, the top volume is selected
      gPad->SetSelected(volume);
      fVolInfo = volume->GetName();
      box = (TGeoBBox*)volume->GetShape();
      memcpy(fCheckedBox, box->GetOrigin(), 3*sizeof(Double_t));
      fCheckedBox[3] = box->GetDX();
      fCheckedBox[4] = box->GetDY();
      fCheckedBox[5] = box->GetDZ();
      return 0;
   }

   TGeoVolume *vol = volume;
   Bool_t vis = vol->IsVisible();
//   Bool_t drawDaughters = kTRUE;
   // Do we need to check a branch only?
   if (volume->IsVisBranch()) {
      if (!fGeoManager->IsClosed()) return big;
      fGeoManager->PushPath();
      fGeoManager->cd(fVisBranch.Data());
      while (fGeoManager->GetLevel()) {
         vol = fGeoManager->GetCurrentVolume();
         *fGlobal = gGeoManager->GetCurrentMatrix();
         dist = vol->GetShape()->DistancetoPrimitive(px,py);
         if (dist<maxdist) {
            fVolInfo = fVisBranch;
            box = (TGeoBBox*)vol->GetShape();
            fGeoManager->LocalToMaster(box->GetOrigin(), &fCheckedBox[0]);
            fCheckedNode = gGeoManager->GetCurrentNode();
            if (fGeoManager->IsNodeSelectable()) gPad->SetSelected(fCheckedNode);
            else gPad->SetSelected(vol);
            fCheckedBox[3] = box->GetDX();
            fCheckedBox[4] = box->GetDY();
            fCheckedBox[5] = box->GetDZ();
            fGeoManager->PopPath();
            return 0;
         }
         fGeoManager->CdUp();
      }
      fGeoManager->PopPath();
      return dist;
   }

   // Do I need to look for the top volume ?
   if ((fTopVisible && vis) || !vol->GetNdaughters() || !vol->IsVisDaughters() || vol->IsVisOnly()) {
      dist = vol->GetShape()->DistancetoPrimitive(px,py);
      if (dist<maxdist) {
         fVolInfo = vol->GetName();
         gPad->SetSelected(vol);
         box = (TGeoBBox*)vol->GetShape();
         memcpy(fCheckedBox, box->GetOrigin(), 3*sizeof(Double_t));
         fCheckedBox[3] = box->GetDX();
         fCheckedBox[4] = box->GetDY();
         fCheckedBox[5] = box->GetDZ();
         return 0;
      }
      if (vol->IsVisOnly() || !vol->GetNdaughters() || !vol->IsVisDaughters())
         return dist;
   }

   // Iterate the volume content
   TGeoIterator next(vol);
   next.SetTopName(TString::Format("%s_1",vol->GetName()));
   TGeoNode *daughter;

   Int_t level, nd;
   Bool_t last;

   while ((daughter=next())) {
      vol = daughter->GetVolume();
      level = next.GetLevel();
      nd = daughter->GetNdaughters();
      vis = daughter->IsVisible();
      if (volume->IsVisContainers()) {
         if (vis && level<=fVisLevel) {
            *fGlobal = next.GetCurrentMatrix();
            dist = vol->GetShape()->DistancetoPrimitive(px,py);
            if (dist<maxdist) {
               next.GetPath(fVolInfo);
               box = (TGeoBBox*)vol->GetShape();
               fGlobal->LocalToMaster(box->GetOrigin(), &fCheckedBox[0]);
               fCheckedNode = daughter;
               if (fGeoManager->IsNodeSelectable()) gPad->SetSelected(fCheckedNode);
               else gPad->SetSelected(vol);
               fCheckedBox[3] = box->GetDX();
               fCheckedBox[4] = box->GetDY();
               fCheckedBox[5] = box->GetDZ();
               fGeoManager->PopPath();
               return 0;
            }
         }
         // Check if we have to skip this branch
         if (level==fVisLevel || !daughter->IsVisDaughters()) {
            next.Skip();
            continue;
         }
      } else if (volume->IsVisLeaves()) {
         last = ((nd==0) || (level==fVisLevel) || (!daughter->IsVisDaughters()))?kTRUE:kFALSE;
         if (vis && last) {
            *fGlobal = next.GetCurrentMatrix();
            dist = vol->GetShape()->DistancetoPrimitive(px,py);
            if (dist<maxdist) {
               next.GetPath(fVolInfo);
               box = (TGeoBBox*)vol->GetShape();
               fGlobal->LocalToMaster(box->GetOrigin(), &fCheckedBox[0]);
               fCheckedNode = daughter;
               if (fGeoManager->IsNodeSelectable()) gPad->SetSelected(fCheckedNode);
               else gPad->SetSelected(vol);
               fCheckedBox[3] = box->GetDX();
               fCheckedBox[4] = box->GetDY();
               fCheckedBox[5] = box->GetDZ();
               fGeoManager->PopPath();
               return 0;
            }
         }
         // Check if we have to skip the branch
         if (last || !daughter->IsVisDaughters()) next.Skip();
      }
   }
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default angles for the current view.

void TGeoPainter::DefaultAngles()
{
   if (gPad) {
      Int_t irep;
      TView *view = gPad->GetView();
      if (!view) return;
      view->SetView(-206,126,75,irep);
      ModifiedPad();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set default volume colors according to tracking media

void TGeoPainter::DefaultColors()
{
   TIter next(fGeoManager->GetListOfVolumes());
   TGeoVolume *vol;
   while ((vol=(TGeoVolume*)next()))
      vol->SetLineColor(vol->GetMaterial()->GetDefaultColor());
   ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Count number of visible nodes down to a given level.

Int_t TGeoPainter::CountNodes(TGeoVolume *volume, Int_t rlevel) const
{
   TGeoVolume *vol = volume;
   Int_t count = 0;
   Bool_t vis = vol->IsVisible();
   // Do I need to look for the top volume ?
   if ((fTopVisible && vis) || !vol->GetNdaughters() || !vol->IsVisDaughters() || vol->IsVisOnly())
      count++;
   // Is this the only volume?
   if (volume->IsVisOnly()) return count;

   // Do we need to check a branch only?
   if (volume->IsVisBranch()) {
      fGeoManager->PushPath();
      fGeoManager->cd(fVisBranch.Data());
      count = fGeoManager->GetLevel() + 1;
      fGeoManager->PopPath();
      return count;
   }
   // Iterate the volume content
   TGeoIterator next(vol);
   TGeoNode *daughter;
   Int_t level, nd;
   Bool_t last;

   while ((daughter=next())) {
//      vol = daughter->GetVolume();
      level = next.GetLevel();
      nd = daughter->GetNdaughters();
      vis = daughter->IsVisible();
      if (volume->IsVisContainers()) {
         if (vis && level<=rlevel) count++;
         // Check if we have to skip this branch
         if (level==rlevel || !daughter->IsVisDaughters()) {
            next.Skip();
            continue;
         }
      } else if (volume->IsVisLeaves()) {
         last = ((nd==0) || (level==rlevel) || (!daughter->IsVisDaughters()))?kTRUE:kFALSE;
         if (vis && last) count++;
         // Check if we have to skip the branch
         if (last) next.Skip();
      }
   }
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Count total number of visible nodes.

Int_t TGeoPainter::CountVisibleNodes()
{
   Int_t maxnodes = fGeoManager->GetMaxVisNodes();
   Int_t vislevel = fGeoManager->GetVisLevel();
//   TGeoVolume *top = fGeoManager->GetTopVolume();
   TGeoVolume *top = fTopVolume;
   if (maxnodes <= 0  && top) {
      fNVisNodes = CountNodes(top, vislevel);
      SetVisLevel(vislevel);
      return fNVisNodes;
   }
   //if (the total number of nodes of the top volume is less than maxnodes
   // we can visualize everything.
   //recompute the best visibility level
   if (!top) {
      SetVisLevel(vislevel);
      return 0;
   }
   fNVisNodes = -1;
   Bool_t again = kFALSE;
   for (Int_t level = 1;level<20;level++) {
      vislevel = level;
      Int_t nnodes = CountNodes(top, level);
      if (top->IsVisOnly() || top->IsVisBranch()) {
         vislevel = fVisLevel;
         fNVisNodes = nnodes;
         break;
      }
      if (nnodes > maxnodes) {
         vislevel--;
         break;
      }
      if (nnodes == fNVisNodes) {
         if (again) break;
         again = kTRUE;
      }
      fNVisNodes = nnodes;
   }
   SetVisLevel(vislevel);
   return fNVisNodes;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if Ged library is loaded and load geometry editor classe.

void TGeoPainter::CheckEdit()
{
   if (fIsEditable) return;
   if (!TClass::GetClass("TGedEditor")) return;
   TPluginHandler *h;
   if ((h = gROOT->GetPluginManager()->FindHandler("TGeoManagerEditor"))) {
      if (h->LoadPlugin() == -1) return;
      h->ExecPlugin(0);
   }
   fIsEditable = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Start the geometry editor.

void TGeoPainter::EditGeometry(Option_t *option)
{
   if (!gPad) return;
   if (!fIsEditable) {
      if (!option[0]) gPad->GetCanvas()->GetCanvasImp()->ShowEditor();
      else TVirtualPadEditor::ShowEditor();
      CheckEdit();
   }
   gPad->SetSelected(fGeoManager);
   gPad->GetCanvas()->Selected(gPad,fGeoManager,kButton1Down);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw method.

void TGeoPainter::Draw(Option_t *option)
{
   DrawVolume(fGeoManager->GetTopVolume(), option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the time evolution of a radionuclide.

void TGeoPainter::DrawBatemanSol(TGeoBatemanSol *sol, Option_t *option)
{
   Int_t ncoeff = sol->GetNcoeff();
   if (!ncoeff) return;
   Double_t tlo=0., thi=0.;
   Double_t cn=0., lambda=0.;
   Int_t i;
   sol->GetRange(tlo, thi);
   Bool_t autorange = (thi==0.)?kTRUE:kFALSE;

   // Try to find the optimum range in time.
   if (autorange) tlo = 0.;
   sol->GetCoeff(0, cn, lambda);
   Double_t lambdamin = lambda;
   TString formula = "";
   for (i=0; i<ncoeff; i++) {
      sol->GetCoeff(i, cn, lambda);
      formula += TString::Format("%g*exp(-%g*x)",cn, lambda);
      if (i < ncoeff-1) formula += "+";
      if (lambda < lambdamin &&
          lambda > 0.) lambdamin = lambda;
   }
   if (autorange) thi = 10./lambdamin;
   // Create a function
   TF1 *func = new TF1(TString::Format("conc%s",sol->GetElement()->GetName()), formula.Data(), tlo,thi);
   func->SetTitle(formula + ";time[s]" + TString::Format(";Concentration_of_%s",sol->GetElement()->GetName()));
   func->SetMinimum(1.e-3);
   func->SetMaximum(1.25*TMath::Max(sol->Concentration(tlo), sol->Concentration(thi)));
   func->SetLineColor(sol->GetLineColor());
   func->SetLineStyle(sol->GetLineStyle());
   func->SetLineWidth(sol->GetLineWidth());
   func->SetMarkerColor(sol->GetMarkerColor());
   func->SetMarkerStyle(sol->GetMarkerStyle());
   func->SetMarkerSize(sol->GetMarkerSize());
   func->Draw(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a polygon in 3D.

void TGeoPainter::DrawPolygon(const TGeoPolygon *poly)
{
   Int_t nvert = poly->GetNvert();
   if (!nvert) {
      Error("DrawPolygon", "No vertices defined");
      return;
   }
   Int_t nconv = poly->GetNconvex();
   Double_t *x = new Double_t[nvert+1];
   Double_t *y = new Double_t[nvert+1];
   poly->GetVertices(x,y);
   x[nvert] = x[0];
   y[nvert] = y[0];
   TGraph *g1 = new TGraph(nvert+1, x,y);
   g1->SetTitle(Form("Polygon with %d vertices (outscribed %d)",nvert, nconv));
   g1->SetLineColor(kRed);
   g1->SetMarkerColor(kRed);
   g1->SetMarkerStyle(4);
   g1->SetMarkerSize(0.8);
   delete [] x;
   delete [] y;
   Double_t *xc = 0;
   Double_t *yc = 0;
   TGraph *g2 = 0;
   if (nconv && !poly->IsConvex()) {
      xc = new Double_t[nconv+1];
      yc = new Double_t[nconv+1];
      poly->GetConvexVertices(xc,yc);
      xc[nconv] = xc[0];
      yc[nconv] = yc[0];
      g2 = new TGraph(nconv+1, xc,yc);
      g2->SetLineColor(kBlue);
      g2->SetLineColor(kBlue);
      g2->SetMarkerColor(kBlue);
      g2->SetMarkerStyle(21);
      g2->SetMarkerSize(0.4);
      delete [] xc;
      delete [] yc;
   }
   if (!gPad) {
      gROOT->MakeDefCanvas();
   }
   g1->Draw("ALP");
   if (g2) g2->Draw("LP");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw method.

void TGeoPainter::DrawVolume(TGeoVolume *vol, Option_t *option)
{
   fTopVolume = vol;
   fLastVolume = nullptr;
   fIsPaintingShape = kFALSE;
//   if (fVisOption==kGeoVisOnly ||
//       fVisOption==kGeoVisBranch) fGeoManager->SetVisOption(kGeoVisLeaves);
   CountVisibleNodes();
   TString opt = option;
   opt.ToLower();
   fPaintingOverlaps = kFALSE;
   fOverlap = nullptr;

   if (fVisLock) {
      ClearVisibleVolumes();
      fVisLock = kFALSE;
   }
   Bool_t has_pad = (gPad==0)?kFALSE:kTRUE;
   // Clear pad if option "same" not given
   if (!gPad) {
      gROOT->MakeDefCanvas();
   }
   if (!opt.Contains("same")) gPad->Clear();
   // append this volume to pad
   fTopVolume->AppendPad(option);

   // Create a 3-D view
   TView *view = gPad->GetView();
   if (!view) {
      view = TView::CreateView(11,0,0);
      // Set the view to perform a first autorange (frame) draw.
      // TViewer3DPad will revert view to normal painting after this
      view->SetAutoRange(kTRUE);
      if (has_pad) gPad->Update();
   }
   if (!opt.Contains("same")) Paint("range");
   else Paint(opt);
   view->SetAutoRange(kFALSE);
   // If we are drawing into the pad, then the view needs to be
   // set to perspective
//   if (!view->IsPerspective()) view->SetPerspective();

   fLastVolume = fTopVolume;

         // Create a 3D viewer to paint us
   gPad->GetViewer3D(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a shape.

void TGeoPainter::DrawShape(TGeoShape *shape, Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   fPaintingOverlaps = kFALSE;
   fOverlap = nullptr;
   fIsPaintingShape = kTRUE;

   Bool_t has_pad = (gPad==0)?kFALSE:kTRUE;
   // Clear pad if option "same" not given
   if (!gPad) {
      gROOT->MakeDefCanvas();
   }
   if (!opt.Contains("same")) gPad->Clear();
   // append this shape to pad
   shape->AppendPad(option);

   // Create a 3-D view
   TView *view = gPad->GetView();
   if (!view) {
      view = TView::CreateView(11,0,0);
      // Set the view to perform a first autorange (frame) draw.
      // TViewer3DPad will revert view to normal painting after this
      view->SetAutoRange(kTRUE);
      if (has_pad) gPad->Update();
   }
   PaintShape(shape,"range");
   view->SetAutoRange(kFALSE);
   view->SetPerspective();
   // Create a 3D viewer to paint us
   gPad->GetViewer3D(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw an overlap.

void TGeoPainter::DrawOverlap(void *ovlp, Option_t *option)
{
   TString opt = option;
   fIsPaintingShape = kFALSE;
   TGeoOverlap *overlap = (TGeoOverlap*)ovlp;
   if (!overlap) return;

   fPaintingOverlaps = kTRUE;
   fOverlap = overlap;
   opt.ToLower();
   if (fVisLock) {
      ClearVisibleVolumes();
      fVisLock = kFALSE;
   }
   Bool_t has_pad = (gPad==0)?kFALSE:kTRUE;
   // Clear pad if option "same" not given
   if (!gPad) {
      gROOT->MakeDefCanvas();
   }
   if (!opt.Contains("same")) gPad->Clear();
   // append this volume to pad
   overlap->AppendPad(option);

   // Create a 3-D view
         // Create a 3D viewer to paint us
   gPad->GetViewer3D(option);
   TView *view = gPad->GetView();
   if (!view) {
      view = TView::CreateView(11,0,0);
      // Set the view to perform a first autorange (frame) draw.
      // TViewer3DPad will revert view to normal painting after this
      view->SetAutoRange(kTRUE);
      PaintOverlap(ovlp, "range");
      overlap->GetPolyMarker()->Draw("SAME");
      if (has_pad) gPad->Update();
   }

   // If we are drawing into the pad, then the view needs to be
   // set to perspective
//   if (!view->IsPerspective()) view->SetPerspective();
   fVisLock = kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw only one volume.

void TGeoPainter::DrawOnly(Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   if (fVisLock) {
      ClearVisibleVolumes();
      fVisLock = kFALSE;
   }
   fPaintingOverlaps = kFALSE;
   fIsPaintingShape = kFALSE;
   Bool_t has_pad = (gPad==0)?kFALSE:kTRUE;
   // Clear pad if option "same" not given
   if (!gPad) {
      gROOT->MakeDefCanvas();
   }
   if (!opt.Contains("same")) gPad->Clear();
   // append this volume to pad
   fTopVolume = fGeoManager->GetCurrentVolume();
   fTopVolume->AppendPad(option);

   // Create a 3-D view
   TView *view = gPad->GetView();
   if (!view) {
      view = TView::CreateView(11,0,0);
      // Set the view to perform a first autorange (frame) draw.
      // TViewer3DPad will revert view to normal painting after this
      view->SetAutoRange(kTRUE);
      fVisOption = kGeoVisOnly;
      if (has_pad) gPad->Update();
   }

   // If we are drawing into the pad, then the view needs to be
   // set to perspective
//   if (!view->IsPerspective()) view->SetPerspective();
   fVisLock = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw current point in the same view.

void TGeoPainter::DrawCurrentPoint(Int_t color)
{
   if (!gPad) return;
   if (!gPad->GetView()) return;
   TPolyMarker3D *pm = new TPolyMarker3D();
   pm->SetMarkerColor(color);
   const Double_t *point = fGeoManager->GetCurrentPoint();
   pm->SetNextPoint(point[0], point[1], point[2]);
   pm->SetMarkerStyle(8);
   pm->SetMarkerSize(0.5);
   pm->Draw("SAME");
}

////////////////////////////////////////////////////////////////////////////////

void TGeoPainter::DrawPanel()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draw all volumes for a given path.

void TGeoPainter::DrawPath(const char *path, Option_t *option)
{
   fVisOption=kGeoVisBranch;
   fVisBranch=path;
   fIsPaintingShape = kFALSE;
   fTopVolume = fGeoManager->GetTopVolume();
   fTopVolume->SetVisRaytrace(kFALSE);
   DrawVolume(fTopVolume,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Estimate camera movement between tmin and tmax for best track display

void TGeoPainter::EstimateCameraMove(Double_t tmin, Double_t tmax, Double_t *start, Double_t *end)
{
   if (!gPad) return;
   TIter next(gPad->GetListOfPrimitives());
   TVirtualGeoTrack *track;
   TObject *obj;
   Int_t ntracks = 0;
   Double_t *point = 0;
   AddTrackPoint(point, start, kTRUE);
   while ((obj=next())) {
      if (strcmp(obj->ClassName(), "TGeoTrack")) continue;
      track = (TVirtualGeoTrack*)obj;
      ntracks++;
      track->PaintCollect(tmin, start);
   }

   if (!ntracks) return;
   next.Reset();
   AddTrackPoint(point, end, kTRUE);
   while ((obj=next())) {
      if (strcmp(obj->ClassName(), "TGeoTrack")) continue;
      track = (TVirtualGeoTrack*)obj;
      if (!track) continue;
      track->PaintCollect(tmax, end);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse actions on a given volume.

void TGeoPainter::ExecuteManagerEvent(TGeoManager * /*geom*/, Int_t event, Int_t /*px*/, Int_t /*py*/)
{
   if (!gPad) return;
   gPad->SetCursor(kPointer);
   switch (event) {
      case kButton1Down:
         if (!fIsEditable) CheckEdit();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse actions on a given shape.

void TGeoPainter::ExecuteShapeEvent(TGeoShape * /*shape*/, Int_t event, Int_t /*px*/, Int_t /*py*/)
{
   if (!gPad) return;
   gPad->SetCursor(kHand);
   switch (event) {
      case kButton1Down:
         if (!fIsEditable) CheckEdit();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse actions on a given volume.

void TGeoPainter::ExecuteVolumeEvent(TGeoVolume * /*volume*/, Int_t event, Int_t /*px*/, Int_t /*py*/)
{
   if (!gPad) return;
   if (!fIsEditable) CheckEdit();
//   if (fIsRaytracing) return;
//   Bool_t istop = (volume==fTopVolume)?kTRUE:kFALSE;
//   if (istop) gPad->SetCursor(kHand);
//   else gPad->SetCursor(kPointer);
   gPad->SetCursor(kHand);
//   static Int_t width, color;
   switch (event) {
   case kMouseEnter:
//      width = volume->GetLineWidth();
//      color = volume->GetLineColor();
      break;

   case kMouseLeave:
//      volume->SetLineWidth(width);
//      volume->SetLineColor(color);
      break;

   case kButton1Down:
//      volume->SetLineWidth(3);
//      volume->SetLineColor(2);
//      gPad->Modified();
//      gPad->Update();
      break;

   case kButton1Up:
//      volume->SetLineWidth(width);
//      volume->SetLineColor(color);
//      gPad->Modified();
//      gPad->Update();
      break;

   case kButton1Double:
      gPad->SetCursor(kWatch);
      GrabFocus();
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get some info about the current selected volume.

const char *TGeoPainter::GetVolumeInfo(const TGeoVolume *volume, Int_t /*px*/, Int_t /*py*/) const
{
   static TString info;
   info = "";
   if (!gPad) return info;
   if (fPaintingOverlaps) {
      if (!fOverlap) {
         info =  "wrong overlapping flag";
         return info;
      }
      TString ovtype, name;
      if (fOverlap->IsExtrusion()) ovtype="EXTRUSION";
      else ovtype = "OVERLAP";
      if (volume==fOverlap->GetFirstVolume()) name=volume->GetName();
      else name=fOverlap->GetSecondVolume()->GetName();
      info = TString::Format("%s: %s of %g", name.Data(), ovtype.Data(), fOverlap->GetOverlap());
      return info;
   }
   else info = TString::Format("%s, shape=%s", fVolInfo.Data(), volume->GetShape()->ClassName());
   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// Create/return geometry checker.

TGeoChecker *TGeoPainter::GetChecker()
{
   if (!fChecker) fChecker = new TGeoChecker(fGeoManager);
   return fChecker;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the current view angles.

void TGeoPainter::GetViewAngles(Double_t &longitude, Double_t &latitude, Double_t &psi)
{
   if (!gPad) return;
   TView *view = gPad->GetView();
   if (!view) return;
   longitude = view->GetLongitude();
   latitude = view->GetLatitude();
   psi = view->GetPsi();
}

////////////////////////////////////////////////////////////////////////////////
/// Move focus to current volume

void TGeoPainter::GrabFocus(Int_t nfr, Double_t dlong, Double_t dlat, Double_t dpsi)
{
   if (!gPad) return;
   TView *view = gPad->GetView();
   if (!view) return;
   if (!fCheckedNode && !fPaintingOverlaps) {
      printf("Woops!!!\n");
      TGeoBBox *box = (TGeoBBox*)fGeoManager->GetTopVolume()->GetShape();
      memcpy(&fCheckedBox[0], box->GetOrigin(), 3*sizeof(Double_t));
      fCheckedBox[3] = box->GetDX();
      fCheckedBox[4] = box->GetDY();
      fCheckedBox[5] = box->GetDZ();
   }
   view->SetPerspective();
   Int_t nvols = fVisVolumes->GetEntriesFast();
   Int_t nframes = nfr;
   if (nfr==0) {
      nframes = 1;
      if (nvols<1500) nframes=10;
      if (nvols<1000) nframes=20;
      if (nvols<200) nframes = 50;
      if (nvols<100) nframes = 100;
   }
   view->MoveFocus(&fCheckedBox[0], fCheckedBox[3], fCheckedBox[4], fCheckedBox[5], nframes, dlong, dlat, dpsi);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a lego plot fot the top volume, according to option.

TH2F *TGeoPainter::LegoPlot(Int_t ntheta, Double_t themin, Double_t themax,
                            Int_t nphi,   Double_t phimin, Double_t phimax,
                            Double_t rmin, Double_t rmax, Option_t *option)
{
   return fChecker->LegoPlot(ntheta, themin, themax, nphi, phimin, phimax, rmin, rmax, option);
}
////////////////////////////////////////////////////////////////////////////////
/// Convert a local vector according view rotation matrix

void TGeoPainter::LocalToMasterVect(const Double_t *local, Double_t *master) const
{
   for (Int_t i=0; i<3; i++)
      master[i] = -local[0]*fMat[i]-local[1]*fMat[i+3]-local[2]*fMat[i+6];
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a pad and view are present and send signal "Modified" to pad.

void TGeoPainter::ModifiedPad(Bool_t update) const
{
   if (!gPad) return;
   if (update) {
      gPad->Update();
      return;
   }
   TView *view = gPad->GetView();
   if (!view) return;
   view->SetViewChanged();
   gPad->Modified();
   if (gROOT->FromPopUp()) gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Paint current geometry according to option.

void TGeoPainter::Paint(Option_t *option)
{
   if (!fGeoManager || !fTopVolume) return;
   Bool_t is_padviewer = kTRUE;
   if (gPad) is_padviewer = (!strcmp(gPad->GetViewer3D()->ClassName(),"TViewer3DPad"))?kTRUE:kFALSE;

   fIsRaytracing = fTopVolume->IsRaytracing();
   if (fTopVolume->IsVisContainers()) fVisOption = kGeoVisDefault;
   else if (fTopVolume->IsVisLeaves()) fVisOption = kGeoVisLeaves;
   else if (fTopVolume->IsVisOnly()) fVisOption = kGeoVisOnly;
   else if (fTopVolume->IsVisBranch()) fVisOption = kGeoVisBranch;


   if (!fIsRaytracing || !is_padviewer) {
      if (fGeoManager->IsDrawingExtra()) {
         // loop the list of physical volumes
         fGeoManager->CdTop();
         TObjArray *nodeList = fGeoManager->GetListOfPhysicalNodes();
         Int_t nnodes = nodeList->GetEntriesFast();
         Int_t inode;
         TGeoPhysicalNode *node;
         for (inode=0; inode<nnodes; inode++) {
            node = (TGeoPhysicalNode*)nodeList->UncheckedAt(inode);
            PaintPhysicalNode(node, option);
         }
      } else {
         PaintVolume(fTopVolume,option);
      }
      fVisLock = kTRUE;
   }
   // Check if we have to raytrace (only in pad)
   if (fIsRaytracing && is_padviewer) Raytrace();
}

////////////////////////////////////////////////////////////////////////////////
/// Paint an overlap.

void TGeoPainter::PaintOverlap(void *ovlp, Option_t *option)
{
   if (!fGeoManager) return;
   TGeoOverlap *overlap = (TGeoOverlap *)ovlp;
   if (!overlap) return;
   Int_t color, transparency;
   if (fOverlap != overlap) fOverlap = overlap;
   TGeoShape::SetTransform(fGlobal);
   TGeoHMatrix *hmat = fGlobal;
   TGeoVolume *vol;
   TGeoVolume *vol1 = overlap->GetFirstVolume();
   TGeoVolume *vol2 = overlap->GetSecondVolume();
   TGeoHMatrix *matrix1 = overlap->GetFirstMatrix();
   TGeoHMatrix *matrix2 = overlap->GetSecondMatrix();
   //
   vol = vol1;
   *hmat = matrix1;
   fGeoManager->SetMatrixReflection(matrix1->IsReflection());
   if (!fVisLock) fVisVolumes->Add(vol);
   fGeoManager->SetPaintVolume(vol);
   color = vol->GetLineColor();
   transparency = vol->GetTransparency();
   vol->SetLineColor(kGreen);
   vol->SetTransparency(40);
   if (!strstr(option,"range")) ((TAttLine*)vol)->Modify();
   PaintShape(*(vol->GetShape()),option);
   vol->SetLineColor(color);
   vol->SetTransparency(transparency);
   vol = vol2;
   *hmat = matrix2;
   fGeoManager->SetMatrixReflection(matrix2->IsReflection());
   if (!fVisLock) fVisVolumes->Add(vol);
   fGeoManager->SetPaintVolume(vol);
   color = vol->GetLineColor();
   transparency = vol->GetTransparency();
   vol->SetLineColor(kBlue);
   vol->SetTransparency(40);
   if (!strstr(option,"range")) ((TAttLine*)vol)->Modify();
   PaintShape(*(vol->GetShape()),option);
   vol->SetLineColor(color);
   vol->SetTransparency(transparency);
   fGeoManager->SetMatrixReflection(kFALSE);
   fVisLock = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint recursively a node and its content according to visualization options.

void TGeoPainter::PaintNode(TGeoNode *node, Option_t *option, TGeoMatrix* global)
{
   PaintVolume(node->GetVolume(), option, global);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint recursively a node and its content according to visualization options.

void TGeoPainter::PaintVolume(TGeoVolume *top, Option_t *option, TGeoMatrix* global)
{
   if (fTopVolume != top) {
      ClearVisibleVolumes();
      fVisLock = kFALSE;
   }
   fTopVolume = top;
   if (!fVisLevel) return;
   TGeoVolume *vol = top;
   if(global)
      *fGlobal = *global;
   else
      fGlobal->Clear();
   TGeoShape::SetTransform(fGlobal);
   Bool_t drawDaughters = kTRUE;
   Bool_t vis = (top->IsVisible() && !top->IsAssembly());
   Int_t transparency = 0;

   // Update pad attributes in case we need to paint VOL
   if (!strstr(option,"range")) ((TAttLine*)vol)->Modify();

   // Do we need to draw a branch ?
   if (top->IsVisBranch()) {
      fGeoManager->PushPath();
      fGeoManager->cd(fVisBranch.Data());
//      while (fGeoManager->GetLevel()) {
         vol = fGeoManager->GetCurrentVolume();
         if (!fVisLock) {
            fVisVolumes->Add(vol);
            vol->SetAttBit(TGeoAtt::kVisOnScreen);
         }
         fGeoManager->SetPaintVolume(vol);
         transparency = vol->GetTransparency();
         vol->SetTransparency(40);
         if (!strstr(option,"range")) ((TAttLine*)vol)->Modify();
         if (global) {
            *fGlobal  = *global;
            *fGlobal *= *fGeoManager->GetCurrentMatrix();
         } else {
            *fGlobal = fGeoManager->GetCurrentMatrix();
         }
         fGeoManager->SetMatrixReflection(fGlobal->IsReflection());
         PaintShape(*(vol->GetShape()),option);
         vol->SetTransparency(transparency);
         fGeoManager->CdUp();
//      }
      fVisLock = kTRUE;
      fGeoManager->PopPath();
      fGeoManager->SetMatrixReflection(kFALSE);
      return;
   }

   // Do I need to draw the top volume ?
   if ((fTopVisible && vis) || !top->GetNdaughters() || !top->IsVisDaughters() || top->IsVisOnly()) {
      fGeoManager->SetPaintVolume(vol);
      fGeoManager->SetMatrixReflection(fGlobal->IsReflection());
      PaintShape(*(vol->GetShape()),option);
      if (!fVisLock && !vol->TestAttBit(TGeoAtt::kVisOnScreen)) {
         fVisVolumes->Add(vol);
         vol->SetAttBit(TGeoAtt::kVisOnScreen);
      }
      if (top->IsVisOnly() || !top->GetNdaughters() || !top->IsVisDaughters()) {
         fVisLock = kTRUE;
         return;
      }
   }

   // Iterate the volume content
   TGeoIterator next(vol);
   if (fPlugin) next.SetUserPlugin(fPlugin);
   TGeoNode *daughter;
//   TGeoMatrix *glmat;
   Int_t level, nd;
   Bool_t last;
   Int_t line_color=0, line_width=0, line_style=0;
   while ((daughter=next())) {
      vol = daughter->GetVolume();
      fGeoManager->SetPaintVolume(vol);
      level = next.GetLevel();
      nd = daughter->GetNdaughters();
      vis = daughter->IsVisible();
      drawDaughters = kTRUE;
      if (top->IsVisContainers()) {
         if (vis && level<=fVisLevel) {
            if (fPlugin) {
               line_color = vol->GetLineColor();
               line_width = vol->GetLineWidth();
               line_style = vol->GetLineStyle();
               transparency = vol->GetTransparency();
               fPlugin->ProcessNode();
            }
            if (!strstr(option,"range")) ((TAttLine*)vol)->Modify();
            if (global) {
               *fGlobal  = *global;
               *fGlobal *= *next.GetCurrentMatrix();
            } else {
               *fGlobal = next.GetCurrentMatrix();
            }
            fGeoManager->SetMatrixReflection(fGlobal->IsReflection());
            drawDaughters = PaintShape(*(vol->GetShape()),option);
            if (fPlugin) {
               vol->SetLineColor(line_color);
               vol->SetLineWidth(line_width);
               vol->SetLineStyle(line_style);
               vol->SetTransparency(transparency);
            }
            if (!fVisLock && !daughter->IsOnScreen()) {
               fVisVolumes->Add(vol);
               vol->SetAttBit(TGeoAtt::kVisOnScreen);
            }
         }
         // Check if we have to skip this branch
         if (!drawDaughters || level==fVisLevel || !daughter->IsVisDaughters()) {
            next.Skip();
            continue;
         }
      } else if (top->IsVisLeaves()) {
         last = ((nd==0) || (level==fVisLevel) || (!daughter->IsVisDaughters()))?kTRUE:kFALSE;
         if (vis && last) {
            if (fPlugin) {
               line_color = vol->GetLineColor();
               line_width = vol->GetLineWidth();
               line_style = vol->GetLineStyle();
               transparency = vol->GetTransparency();
               fPlugin->ProcessNode();
            }
            if (!strstr(option,"range")) ((TAttLine*)vol)->Modify();
            if (global) {
               *fGlobal  = *global;
               *fGlobal *= *next.GetCurrentMatrix();
            } else {
               *fGlobal = next.GetCurrentMatrix();
            }
            fGeoManager->SetMatrixReflection(fGlobal->IsReflection());
            drawDaughters = PaintShape(*(vol->GetShape()),option);
            if (fPlugin) {
               vol->SetLineColor(line_color);
               vol->SetLineWidth(line_width);
               vol->SetLineStyle(line_style);
               vol->SetTransparency(transparency);
            }
            if (!fVisLock && !daughter->IsOnScreen()) {
               fVisVolumes->Add(vol);
               vol->SetAttBit(TGeoAtt::kVisOnScreen);
            }
         }
         // Check if we have to skip the branch
         if (!drawDaughters || last || !daughter->IsVisDaughters()) next.Skip();
      }
   }
   if (fPlugin) fPlugin->SetIterator(0);
   fGeoManager->SetMatrixReflection(kFALSE);
   fVisLock = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the supplied shape into the current 3D viewer

Bool_t TGeoPainter::PaintShape(const TGeoShape &shape, Option_t *option) const
{
   Bool_t addDaughters = kTRUE;

   TVirtualViewer3D *viewer = gPad->GetViewer3D();

   if (!viewer || shape.IsA()==TGeoShapeAssembly::Class()) {
      return addDaughters;
   }

   // For non-composite shapes we are the main paint method & perform the negotiation
   // with the viewer here
   if (!shape.IsComposite()) {
      // Does viewer prefer local frame positions?
      Bool_t localFrame = viewer->PreferLocalFrame();
      // Perform first fetch of buffer from the shape and try adding it
      // to the viewer
      const TBuffer3D & buffer =
         shape.GetBuffer3D(TBuffer3D::kCore|TBuffer3D::kBoundingBox|TBuffer3D::kShapeSpecific, localFrame);
      Int_t reqSections = viewer->AddObject(buffer, &addDaughters);

      // If the viewer requires additional sections fetch from the shape (if possible)
      // and add again
      if (reqSections != TBuffer3D::kNone) {
         shape.GetBuffer3D(reqSections, localFrame);
         viewer->AddObject(buffer, &addDaughters);
      }
   }
   // Composite shapes have their own internal hierarchy of shapes, each
   // of which generate a filled TBuffer3D. Therefore we can't pass up a
   // single buffer to here. So as a special case the TGeoCompositeShape
   // performs it's own painting & negotiation with the viewer.
   else {
      const TGeoCompositeShape * composite = static_cast<const TGeoCompositeShape *>(&shape);

      // We need the addDaughters flag returned from the viewer from paint
      // so can't use the normal TObject::Paint()
//      TGeoHMatrix *matrix = (TGeoHMatrix*)TGeoShape::GetTransform();
//      if (viewer->PreferLocalFrame()) matrix->Clear();
      addDaughters = composite->PaintComposite(option);
   }

   return addDaughters;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint an overlap.

void TGeoPainter::PaintShape(TGeoShape *shape, Option_t *option)
{
   TGeoShape::SetTransform(fGlobal);
   fGlobal->Clear();
   fGeoManager->SetPaintVolume(0);
   PaintShape(*shape,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Paints a physical node associated with a path.

void TGeoPainter::PaintPhysicalNode(TGeoPhysicalNode *node, Option_t *option)
{
   if (!node->IsVisible()) return;
   Int_t level = node->GetLevel();
   Int_t i, col, wid, sty;
   TGeoShape *shape;
   TGeoShape::SetTransform(fGlobal);
   TGeoHMatrix *matrix = fGlobal;
   TGeoVolume *vcrt;
   if (!node->IsVisibleFull()) {
      // Paint only last node in the branch
      vcrt  = node->GetVolume();
      if (!strstr(option,"range")) ((TAttLine*)vcrt)->Modify();
      shape = vcrt->GetShape();
      *matrix = node->GetMatrix();
      fGeoManager->SetMatrixReflection(matrix->IsReflection());
      fGeoManager->SetPaintVolume(vcrt);
      if (!node->IsVolAttributes() && !strstr(option,"range")) {
         col = vcrt->GetLineColor();
         wid = vcrt->GetLineWidth();
         sty = vcrt->GetLineStyle();
         vcrt->SetLineColor(node->GetLineColor());
         vcrt->SetLineWidth(node->GetLineWidth());
         vcrt->SetLineStyle(node->GetLineStyle());
         ((TAttLine*)vcrt)->Modify();
         PaintShape(*shape,option);
         vcrt->SetLineColor(col);
         vcrt->SetLineWidth(wid);
         vcrt->SetLineStyle(sty);
      } else {
         PaintShape(*shape,option);
      }
   } else {
      // Paint full branch, except top node
      for (i=1;i<=level; i++) {
         vcrt  = node->GetVolume(i);
         if (!strstr(option,"range")) ((TAttLine*)vcrt)->Modify();
         shape = vcrt->GetShape();
         *matrix = node->GetMatrix(i);
         fGeoManager->SetMatrixReflection(matrix->IsReflection());
         fGeoManager->SetPaintVolume(vcrt);
         if (!node->IsVolAttributes() && !strstr(option,"range")) {
            col = vcrt->GetLineColor();
            wid = vcrt->GetLineWidth();
            sty = vcrt->GetLineStyle();
            vcrt->SetLineColor(node->GetLineColor());
            vcrt->SetLineWidth(node->GetLineWidth());
            vcrt->SetLineStyle(node->GetLineStyle());
            ((TAttLine*)vcrt)->Modify();
            PaintShape(*shape,option);
            vcrt->SetLineColor(col);
            vcrt->SetLineWidth(wid);
            vcrt->SetLineStyle(sty);
         } else {
            PaintShape(*shape,option);
         }
      }
   }
   fGeoManager->SetMatrixReflection(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Print overlaps (see TGeoChecker::PrintOverlaps())

void TGeoPainter::PrintOverlaps() const
{
   fChecker->PrintOverlaps();
}

////////////////////////////////////////////////////////////////////////////////
/// Text progress bar.

void TGeoPainter::OpProgress(const char *opname, Long64_t current, Long64_t size, TStopwatch *watch, Bool_t last, Bool_t refresh, const char *msg)
{
   fChecker->OpProgress(opname,current,size,watch,last,refresh, msg);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw random points in the bounding box of a volume.

void TGeoPainter::RandomPoints(const TGeoVolume *vol, Int_t npoints, Option_t *option)
{
   fChecker->RandomPoints((TGeoVolume*)vol, npoints, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Shoot nrays in the current drawn geometry

void TGeoPainter::RandomRays(Int_t nrays, Double_t startx, Double_t starty, Double_t startz, const char *target_vol, Bool_t check_norm)
{
   fChecker->RandomRays(nrays, startx, starty, startz, target_vol, check_norm);
}

////////////////////////////////////////////////////////////////////////////////
/// Raytrace current drawn geometry

void TGeoPainter::Raytrace(Option_t *)
{
   if (!gPad || gPad->IsBatch()) return;
   TView *view = gPad->GetView();
   if (!view) return;
   Int_t rtMode = fGeoManager->GetRTmode();
   TGeoVolume *top = fGeoManager->GetTopVolume();
   if (top != fTopVolume) fGeoManager->SetTopVolume(fTopVolume);
   if (!view->IsPerspective()) view->SetPerspective();
   gVirtualX->SetMarkerSize(1);
   gVirtualX->SetMarkerStyle(1);
   Bool_t inclipst=kFALSE, inclip=kFALSE;
   Double_t krad = TMath::DegToRad();
   Double_t lat = view->GetLatitude();
   Double_t longit = view->GetLongitude();
   Double_t psi = view->GetPsi();
   Double_t c1 = TMath::Cos(psi*krad);
   Double_t s1 = TMath::Sin(psi*krad);
   Double_t c2 = TMath::Cos(lat*krad);
   Double_t s2 = TMath::Sin(lat*krad);
   Double_t s3 = TMath::Cos(longit*krad);
   Double_t c3 = -TMath::Sin(longit*krad);
   fMat[0] =  c1*c3 - s1*c2*s3;
   fMat[1] =  c1*s3 + s1*c2*c3;
   fMat[2] =  s1*s2;

   fMat[3] =  -s1*c3 - c1*c2*s3;
   fMat[4] = -s1*s3 + c1*c2*c3;
   fMat[5] =  c1*s2;

   fMat[6] =  s2*s3;
   fMat[7] =  -s2*c3;
   fMat[8] = c2;
   Double_t u0, v0, du, dv;
   view->GetWindow(u0,v0,du,dv);
   Double_t dview = view->GetDview();
   Double_t dproj = view->GetDproj();
   Double_t local[3] = {0,0,1};
   Double_t dir[3], normal[3];
   LocalToMasterVect(local,dir);
   Double_t min[3], max[3];
   view->GetRange(min, max);
   Double_t cov[3];
   for (Int_t i=0; i<3; i++) cov[i] = 0.5*(min[i]+max[i]);
   Double_t cop[3];
   for (Int_t i=0; i<3; i++) cop[i] = cov[i] - dir[i]*dview;
   fGeoManager->InitTrack(cop, dir);
   Bool_t outside = fGeoManager->IsOutside();
   fGeoManager->DoBackupState();
   if (fClippingShape) inclipst = inclip = fClippingShape->Contains(cop);
   Int_t px, py;
   Double_t xloc, yloc, modloc;
   Int_t pxmin,pxmax, pymin,pymax;
   pxmin = gPad->UtoAbsPixel(0);
   pxmax = gPad->UtoAbsPixel(1);
   pymin = gPad->VtoAbsPixel(1);
   pymax = gPad->VtoAbsPixel(0);
   TGeoNode *next = nullptr;
   TGeoNode *nextnode = nullptr;
   Double_t step,steptot;
   Double_t *norm;
   const Double_t *point = fGeoManager->GetCurrentPoint();
   Double_t *ppoint = (Double_t*)point;
   Double_t tosource[3];
   Double_t calf;
   Double_t phi = 45.*krad;
   tosource[0] = -dir[0]*TMath::Cos(phi)+dir[1]*TMath::Sin(phi);
   tosource[1] = -dir[0]*TMath::Sin(phi)-dir[1]*TMath::Cos(phi);
   tosource[2] = -dir[2];

   Bool_t done;
//   Int_t istep;
   Int_t base_color, color;
   Double_t light;
   Double_t stemin=0, stemax=TGeoShape::Big();
   TPoint *pxy = new TPoint[1];
   TGeoVolume *nextvol;
   Int_t up;
   Int_t ntotal = pxmax*pymax;
   Int_t nrays = 0;
   TStopwatch *timer = new TStopwatch();
   timer->Start();
   for (px=pxmin; px<pxmax; px++) {
      for (py=pymin; py<pymax; py++) {
         if ((nrays%100)==0) OpProgress("Raytracing",nrays,ntotal,timer,kFALSE);
         nrays++;
         base_color = 1;
         steptot = 0;
         inclip = inclipst;
         xloc = gPad->AbsPixeltoX(pxmin+pxmax-px);
         xloc = xloc*du-u0;
         yloc = gPad->AbsPixeltoY(pymin+pymax-py);
         yloc = yloc*dv-v0;
         modloc = TMath::Sqrt(xloc*xloc+yloc*yloc+dproj*dproj);
         local[0] = xloc/modloc;
         local[1] = yloc/modloc;
         local[2] = dproj/modloc;
         LocalToMasterVect(local,dir);
         fGeoManager->DoRestoreState();
         fGeoManager->SetOutside(outside);
         fGeoManager->SetCurrentPoint(cop);
         fGeoManager->SetCurrentDirection(dir);
//         fGeoManager->InitTrack(cop,dir);
         // current ray pointing to pixel (px,py)
         done = kFALSE;
         norm = 0;
         // propagate to the clipping shape if any
         if (fClippingShape) {
            if (inclip) {
               stemin = fClippingShape->DistFromInside(cop,dir,3);
               stemax = TGeoShape::Big();
            } else {
               stemax = fClippingShape->DistFromOutside(cop,dir,3);
               stemin = 0;
            }
         }

         while (!done) {
            if (fClippingShape) {
               if (stemin>1E10) break;
               if (stemin>0) {
                  // we are inside clipping shape
                  fGeoManager->SetStep(stemin);
                  next = fGeoManager->Step();
                  steptot = 0;
                  stemin = 0;
                  up = 0;
                  while (next) {
                     // we found something after clipping region
                     nextvol = next->GetVolume();
                     if (nextvol->TestAttBit(TGeoAtt::kVisOnScreen)) {
                        done = kTRUE;
                        base_color = nextvol->GetLineColor();
                        fClippingShape->ComputeNormal(ppoint, dir, normal);
                        norm = normal;
                        break;
                     }
                     up++;
                     next = fGeoManager->GetMother(up);
                  }
                  if (done) break;
                  inclip = fClippingShape->Contains(ppoint);
                  fGeoManager->SetStep(1E-3);
                  while (inclip) {
                     fGeoManager->Step();
                     inclip = fClippingShape->Contains(ppoint);
                  }
                  stemax = fClippingShape->DistFromOutside(ppoint,dir,3);
               }
            }
            nextnode = fGeoManager->FindNextBoundaryAndStep();
            step = fGeoManager->GetStep();
            if (step>1E10) break;
            steptot += step;
            next = nextnode;
            // Check the step
            if (fClippingShape) {
               if (steptot>stemax) {
                  steptot = 0;
                  inclip = fClippingShape->Contains(ppoint);
                  if (inclip) {
                     stemin = fClippingShape->DistFromInside(ppoint,dir,3);
                     stemax = TGeoShape::Big();
                     continue;
                  } else {
                     stemin = 0;
                     stemax = fClippingShape->DistFromOutside(ppoint,dir,3);
                  }
               }
            }
            // Check if next node is visible
            if (!nextnode) continue;
            nextvol = nextnode->GetVolume();
            if (nextvol->TestAttBit(TGeoAtt::kVisOnScreen)) {
               done = kTRUE;
               base_color = nextvol->GetLineColor();
               next = nextnode;
               break;
            }
         }
         if (!done) continue;
         // current ray intersect a visible volume having color=base_color
         if (rtMode > 0) {
            fGeoManager->MasterToLocal(gGeoManager->GetCurrentPoint(), local);
            fGeoManager->MasterToLocalVect(gGeoManager->GetCurrentDirection(), dir);
            for (Int_t i=0; i<3; ++i) local[i] += 1.E-8*dir[i];
            step = next->GetVolume()->GetShape()->DistFromInside(local,dir,3);
            for (Int_t i=0; i<3; ++i) local[i] += step*dir[i];
            next->GetVolume()->GetShape()->ComputeNormal(local, dir, normal);
            norm = normal;
         } else {
            if (!norm) norm = fGeoManager->FindNormalFast();
            if (!norm) continue;
         }
         calf = norm[0]*tosource[0]+norm[1]*tosource[1]+norm[2]*tosource[2];
         light = TMath::Abs(calf);
         color = GetColor(base_color, light);
         // Now we know the color of the pixel, just draw it
         gVirtualX->SetMarkerColor(color);
         pxy[0].fX = px;
         pxy[0].fY = py;
         gVirtualX->DrawPolyMarker(1,pxy);
      }
   }
   delete [] pxy;
   timer->Stop();
   fChecker->OpProgress("Raytracing",nrays,ntotal,timer,kTRUE);
   delete timer;
}

////////////////////////////////////////////////////////////////////////////////
/// Shoot npoints randomly in a box of 1E-5 around current point.
/// Return minimum distance to points outside.

TGeoNode *TGeoPainter::SamplePoints(Int_t npoints, Double_t &dist, Double_t epsil,
                                    const char* g3path)
{
   return fChecker->SamplePoints(npoints, dist, epsil, g3path);
}

////////////////////////////////////////////////////////////////////////////////
/// Set cartesian and radial bomb factors for translations.

void TGeoPainter::SetBombFactors(Double_t bombx, Double_t bomby, Double_t bombz, Double_t bombr)
{
   fBombX = bombx;
   fBombY = bomby;
   fBombZ = bombz;
   fBombR = bombr;
   if (IsExplodedView()) ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Set type of exploding view.

void TGeoPainter::SetExplodedView(Int_t ibomb)
{
   if ((ibomb<0) || (ibomb>3)) {
      Warning("SetExplodedView", "exploded view can be 0-3");
      return;
   }
   if ((Int_t)ibomb==fExplodedView) return;
   Bool_t change = (gPad==0)?kFALSE:kTRUE;

   if (ibomb==kGeoNoBomb) {
      change &= ((fExplodedView==kGeoNoBomb)?kFALSE:kTRUE);
   }
   if (ibomb==kGeoBombXYZ) {
      change &= ((fExplodedView==kGeoBombXYZ)?kFALSE:kTRUE);
   }
   if (ibomb==kGeoBombCyl) {
      change &= ((fExplodedView==kGeoBombCyl)?kFALSE:kTRUE);
   }
   if (ibomb==kGeoBombSph) {
      change &= ((fExplodedView==kGeoBombSph)?kFALSE:kTRUE);
   }
   fExplodedView = ibomb;
   if (change) ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of segments to approximate circles.

void TGeoPainter::SetNsegments(Int_t nseg)
{
   if (nseg<3) {
      Warning("SetNsegments", "number of segments should be > 2");
      return;
   }
   if (fNsegments==nseg) return;
   fNsegments = nseg;
   ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of points to be generated on the shape outline when checking for overlaps.

void TGeoPainter::SetNmeshPoints(Int_t npoints) {
   fChecker->SetNmeshPoints(npoints);
}

////////////////////////////////////////////////////////////////////////////////
/// Select a node to be checked for overlaps. All overlaps not involving it will
/// be ignored.

void TGeoPainter::SetCheckedNode(TGeoNode *node) {
   fChecker->SetSelectedNode(node);
}

////////////////////////////////////////////////////////////////////////////////
/// Set default level down to which visualization is performed

void TGeoPainter::SetVisLevel(Int_t level) {
   if (level==fVisLevel && fLastVolume==fTopVolume) return;
   fVisLevel=level;
   if (!fTopVolume) return;
   if (fVisLock) {
      ClearVisibleVolumes();
      fVisLock = kFALSE;
   }
   if (!fLastVolume) {
//      printf("--- Drawing   %6d nodes with %d visible levels\n",fNVisNodes,fVisLevel);
      return;
   }
   if (!gPad) return;
   if (gPad->GetView()) {
//      printf("--- Drawing   %6d nodes with %d visible levels\n",fNVisNodes,fVisLevel);
      ModifiedPad();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set top geometry volume as visible.

void TGeoPainter::SetTopVisible(Bool_t vis)
{
   if (fTopVisible==vis) return;
   fTopVisible = vis;
   ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Set drawing mode :
///  - option=0 (default) all nodes drawn down to vislevel
///  - option=1           leaves and nodes at vislevel drawn
///  - option=2           path is drawn

void TGeoPainter::SetVisOption(Int_t option) {
   if ((fVisOption<0) || (fVisOption>4)) {
      Warning("SetVisOption", "wrong visualization option");
      return;
   }

   if (option == kGeoVisChanged) {
      if (fVisLock) {
         ClearVisibleVolumes();
         fVisLock = kFALSE;
      }
      ModifiedPad();
      return;
   }

   if (fTopVolume) {
      TGeoAtt *att = (TGeoAtt*)fTopVolume;
      att->SetAttBit(TGeoAtt::kVisBranch,kFALSE);
      att->SetAttBit(TGeoAtt::kVisContainers,kFALSE);
      att->SetAttBit(TGeoAtt::kVisOnly,kFALSE);
      switch (option) {
         case kGeoVisDefault:
            att->SetAttBit(TGeoAtt::kVisContainers,kTRUE);
            break;
         case kGeoVisLeaves:
            break;
         case kGeoVisOnly:
            att->SetAttBit(TGeoAtt::kVisOnly,kTRUE);
            break;
      }
   }

   if (fVisOption==option) return;
   fVisOption=option;
   if (fVisLock) {
      ClearVisibleVolumes();
      fVisLock = kFALSE;
   }
   ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
///  Returns distance between point px,py on the pad an a shape.

Int_t TGeoPainter::ShapeDistancetoPrimitive(const TGeoShape *shape, Int_t numpoints, Int_t px, Int_t py) const
{
   const Int_t inaxis = 7;
   const Int_t maxdist = 5;
   const Int_t big = 9999;
   Int_t dist = big;
   if (!gPad) return dist;
   TView *view = gPad->GetView();
   if (!(numpoints && view)) return dist;
   if (shape->IsA()==TGeoShapeAssembly::Class()) return dist;

   if (fIsPaintingShape) {
      Int_t puxmin = gPad->XtoAbsPixel(gPad->GetUxmin());
      Int_t puymin = gPad->YtoAbsPixel(gPad->GetUymin());
      Int_t puxmax = gPad->XtoAbsPixel(gPad->GetUxmax());
      Int_t puymax = gPad->YtoAbsPixel(gPad->GetUymax());
      // return if point not in user area
      if (px < puxmin - inaxis) return big;
      if (py > puymin + inaxis) return big;
      if (px > puxmax + inaxis) return big;
      if (py < puymax - inaxis) return big;
      if ((puxmax+inaxis-px) < 40) {
         // when the mouse points to the (40 pix) right edge of the pad, the manager class is selected
         gPad->SetSelected(fGeoManager);
         return 0;
      }
   }

   fBuffer->SetRawSizes(numpoints, 3*numpoints, 0, 0, 0, 0);
   Double_t *points = fBuffer->fPnts;
   shape->SetPoints(points);
   Double_t dpoint2, x1, y1, xndc[3];
   Double_t dmaster[3];
   Int_t j;
   for (Int_t i=0; i<numpoints; i++) {
      j = 3*i;
      TGeoShape::GetTransform()->LocalToMaster(&points[j], dmaster);
      points[j]=dmaster[0]; points[j+1]=dmaster[1]; points[j+2]=dmaster[2];
      view->WCtoNDC(&points[j], xndc);
      x1 = gPad->XtoAbsPixel(xndc[0]);
      y1 = gPad->YtoAbsPixel(xndc[1]);
      dpoint2 = (px-x1)*(px-x1) + (py-y1)*(py-y1);
      if (dpoint2 < dist) dist=(Int_t)dpoint2;
   }
   if (dist > 100) return dist;
   dist = Int_t(TMath::Sqrt(Double_t(dist)));
   if (dist<maxdist && fIsPaintingShape) gPad->SetSelected((TObject*)shape);
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Check time of finding "Where am I" for n points.

void TGeoPainter::Test(Int_t npoints, Option_t *option)
{
   fChecker->Test(npoints, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry overlap checker based on sampling.

void TGeoPainter::TestOverlaps(const char* path)
{
   fChecker->TestOverlaps(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Check voxels efficiency per volume.

Bool_t TGeoPainter::TestVoxels(TGeoVolume *vol)
{
   return fChecker->TestVoxels(vol);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the new 'unbombed' translation vector according current exploded view mode.

void TGeoPainter::UnbombTranslation(const Double_t *tr, Double_t *bombtr)
{
   memcpy(bombtr, tr, 3*sizeof(Double_t));
   switch (fExplodedView) {
      case kGeoNoBomb:
         return;
      case kGeoBombXYZ:
         bombtr[0] /= fBombX;
         bombtr[1] /= fBombY;
         bombtr[2] /= fBombZ;
         return;
      case kGeoBombCyl:
         bombtr[0] /= fBombR;
         bombtr[1] /= fBombR;
         bombtr[2] /= fBombZ;
         return;
      case kGeoBombSph:
         bombtr[0] /= fBombR;
         bombtr[1] /= fBombR;
         bombtr[2] /= fBombR;
         return;
      default:
         return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute weight [kg] of the current volume.

Double_t TGeoPainter::Weight(Double_t precision, Option_t *option)
{
   return fChecker->Weight(precision, option);
}


