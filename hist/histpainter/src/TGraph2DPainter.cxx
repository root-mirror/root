// @(#)root/histpainter:$Id: TGraph2DPainter.cxx,v 1.00
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGraph2DPainter.h"
#include "TMath.h"
#include "TList.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TGraphDelaunay2D.h"
#include "TPolyLine.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "THLimitsFinder.h"
#include "TStyle.h"
#include "Hoption.h"
#include "TH1.h"

R__EXTERN TH1  *gCurrentHist;
R__EXTERN Hoption_t Hoption;

ClassImp(TGraph2DPainter);


////////////////////////////////////////////////////////////////////////////////
/*! \class TGraph2DPainter
    \ingroup Histpainter
    \brief The TGraphDelaunay painting class.

TGraph2DPainter paints a TGraphDelaunay using triangles or clouds of points.

*/


////////////////////////////////////////////////////////////////////////////////
/// TGraph2DPainter default constructor

TGraph2DPainter::TGraph2DPainter()
{
   fX        = 0;
   fY        = 0;
   fZ        = 0;
   fEX       = 0;
   fEY       = 0;
   fEZ       = 0;
   fXN       = 0;
   fYN       = 0;
   fPTried   = 0;
   fNTried   = 0;
   fMTried   = 0;
   fGraph2D  = 0;
   fDelaunay = 0;
   fDelaunay2D = 0;
   fXmin     = 0.;
   fXmax     = 0.;
   fYmin     = 0.;
   fYmax     = 0.;
   fZmin     = 0.;
   fZmax     = 0.;
   fXNmin    = 0;
   fXNmax    = 0;
   fYNmin    = 0;
   fYNmax    = 0;
   fNdt      = 0;
   fNpoints  = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraph2DPainter constructor using the old Delaunay algorithm

TGraph2DPainter::TGraph2DPainter(TGraphDelaunay *gd)
{
   fDelaunay = gd;
   fDelaunay2D = 0;
   fGraph2D  = fDelaunay->GetGraph2D();
   fNpoints  = fGraph2D->GetN();
   fX        = fGraph2D->GetX();
   fY        = fGraph2D->GetY();
   fZ        = fGraph2D->GetZ();
   fEX       = fGraph2D->GetEX();
   fEY       = fGraph2D->GetEY();
   fEZ       = fGraph2D->GetEZ();
   fNdt      = 0;
   fXN       = 0;
   fYN       = 0;
   fXNmin    = 0;
   fXNmax    = 0;
   fYNmin    = 0;
   fYNmax    = 0;
   fPTried   = 0;
   fNTried   = 0;
   fMTried   = 0;
   fXmin     = 0.;
   fXmax     = 0.;
   fYmin     = 0.;
   fYmax     = 0.;
   fZmin     = 0.;
   fZmax     = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// TGraph2DPainter constructor using the new Delaunay algorithm

TGraph2DPainter::TGraph2DPainter(TGraphDelaunay2D *gd)
{
   fDelaunay = 0;
   fDelaunay2D = gd;
   fGraph2D  = fDelaunay2D->GetGraph2D();
   fNpoints  = fGraph2D->GetN();
   fX        = fGraph2D->GetX();
   fY        = fGraph2D->GetY();
   fZ        = fGraph2D->GetZ();
   fEX       = fGraph2D->GetEX();
   fEY       = fGraph2D->GetEY();
   fEZ       = fGraph2D->GetEZ();
   fNdt      = 0;
   fXN       = 0;
   fYN       = 0;
   fXNmin    = 0;
   fXNmax    = 0;
   fYNmin    = 0;
   fYNmax    = 0;
   fPTried   = 0;
   fNTried   = 0;
   fMTried   = 0;
   fXmin     = 0.;
   fXmax     = 0.;
   fYmin     = 0.;
   fYmax     = 0.;
   fZmin     = 0.;
   fZmax     = 0.;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraph2DPainter destructor.

TGraph2DPainter::~TGraph2DPainter()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Find triangles in fDelaunay and initialise the TGraph2DPainter values
/// needed to paint triangles or find contours.

void TGraph2DPainter::FindTriangles()
{
   if (fDelaunay) {
      fDelaunay->FindAllTriangles();
      fNdt    = fDelaunay->GetNdt();
      fXN     = fDelaunay->GetXN();
      fYN     = fDelaunay->GetYN();
      fXNmin  = fDelaunay->GetXNmin();
      fXNmax  = fDelaunay->GetXNmax();
      fYNmin  = fDelaunay->GetYNmin();
      fYNmax  = fDelaunay->GetYNmax();
      fPTried = fDelaunay->GetPTried();
      fNTried = fDelaunay->GetNTried();
      fMTried = fDelaunay->GetMTried();
   }
   else if (fDelaunay2D) {
      fDelaunay2D->FindAllTriangles();
      fNdt    = fDelaunay2D->GetNdt();
      fXNmin  = fDelaunay2D->GetXNmin();
      fXNmax  = fDelaunay2D->GetXNmax();
      fYNmin  = fDelaunay2D->GetYNmin();
      fYNmax  = fDelaunay2D->GetYNmax();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the X and Y graphs building a contour. A contour level may
/// consist in several parts not connected to each other. This function
/// finds them and returns them in a graphs' list.

TList *TGraph2DPainter::GetContourList(Double_t contour)
{
   // Exit if the contour is outside the Z range.
   Double_t zmin = gCurrentHist->GetMinimum();
   Double_t zmax = gCurrentHist->GetMaximum();
   if (Hoption.Logz) {
      if (zmin > 0) {
         zmin = TMath::Log10(zmin);
         zmax = TMath::Log10(zmax);
      } else {
         return 0;
      }
   }
   if(contour<zmin || contour>zmax) {
      Error("GetContourList", "Contour level (%g) outside the Z scope [%g,%g]",
      contour,zmin,zmax);
      return 0;
   }

   if (!fNdt) FindTriangles();

   TGraph *graph = nullptr;     // current graph
   Int_t npg     = 0;           // number of points in the current graph

   // Find all the segments making the contour

   Double_t r21, r20, r10;
   Int_t p0, p1, p2;
   Double_t x0, y0, z0;
   Double_t x1, y1, z1;
   Double_t x2, y2, z2;
   Int_t t[3],i,it,i0,i1,i2;

   // Allocate space to store the segments. They cannot be more than the
   // number of triangles.
   Double_t xs0c, ys0c, xs1c, ys1c;
   Double_t *xs0 = new Double_t[fNdt];
   Double_t *ys0 = new Double_t[fNdt];
   Double_t *xs1 = new Double_t[fNdt];
   Double_t *ys1 = new Double_t[fNdt];
   for (i=0;i<fNdt;i++) {
      xs0[i] = 0.;
      ys0[i] = 0.;
      xs1[i] = 0.;
      ys1[i] = 0.;
   }
   Int_t nbSeg   = 0;

   // Loop over all the triangles in order to find all the line segments
   // making the contour.

   // old implementation
   if (fDelaunay) {
      for(it=0; it<fNdt; it++) {
         t[0] = fPTried[it];
         t[1] = fNTried[it];
         t[2] = fMTried[it];
         p0   = t[0]-1;
         p1   = t[1]-1;
         p2   = t[2]-1;
         x0   = fX[p0]; x2 = fX[p0];
         y0   = fY[p0]; y2 = fY[p0];
         z0   = fZ[p0]; z2 = fZ[p0];

         // Order along Z axis the points (xi,yi,zi) where "i" belongs to {0,1,2}
         // After this z0 < z1 < z2
         i0 = i1 = i2 = 0;
         if (fZ[p1]<=z0) {z0=fZ[p1]; x0=fX[p1]; y0=fY[p1]; i0=1;}
         if (fZ[p1]>z2)  {z2=fZ[p1]; x2=fX[p1]; y2=fY[p1]; i2=1;}
         if (fZ[p2]<=z0) {z0=fZ[p2]; x0=fX[p2]; y0=fY[p2]; i0=2;}
         if (fZ[p2]>z2)  {z2=fZ[p2]; x2=fX[p2]; y2=fY[p2]; i2=2;}
         if (i0==0 && i2==0) {
            Error("GetContourList", "wrong vertices ordering");
            delete [] xs0;
            delete [] ys0;
            delete [] xs1;
            delete [] ys1;
            return nullptr;
         } else {
            i1 = 3-i2-i0;
         }
         x1 = fX[t[i1]-1];
         y1 = fY[t[i1]-1];
         z1 = fZ[t[i1]-1];

         if (Hoption.Logz) {
            z0 = TMath::Log10(z0);
            z1 = TMath::Log10(z1);
            z2 = TMath::Log10(z2);
         }

         if(contour >= z0 && contour <=z2) {
            r20 = (contour-z0)/(z2-z0);
            xs0c = r20*(x2-x0)+x0;
            ys0c = r20*(y2-y0)+y0;
            if(contour >= z1 && contour <=z2) {
               r21 = (contour-z1)/(z2-z1);
               xs1c = r21*(x2-x1)+x1;
               ys1c = r21*(y2-y1)+y1;
            } else {
               r10 = (contour-z0)/(z1-z0);
               xs1c = r10*(x1-x0)+x0;
               ys1c = r10*(y1-y0)+y0;
            }
            // do not take the segments equal to a point
            if(xs0c != xs1c || ys0c != ys1c) {
               nbSeg++;
               xs0[nbSeg-1] = xs0c;
               ys0[nbSeg-1] = ys0c;
               xs1[nbSeg-1] = xs1c;
               ys1[nbSeg-1] = ys1c;
            }
         }
      }
   }
   else if (fDelaunay2D) {

      Int_t p[3];
      for(const auto & face : *fDelaunay2D) {
         p[0]   = face.idx[0];
         p[1]   = face.idx[1];
         p[2]   = face.idx[2];
         x0   = fX[p[0]]; x2 = fX[p[0]];
         y0   = fY[p[0]]; y2 = fY[p[0]];
         z0   = fZ[p[0]]; z2 = fZ[p[0]];

         // Order along Z axis the points (xi,yi,zi) where "i" belongs to {0,1,2}
         // After this z0 < z1 < z2
         i0 = i1 = i2 = 0;
         if (fZ[p[1]]<=z0) {z0=fZ[p[1]]; x0=fX[p[1]]; y0=fY[p[1]]; i0=1;}
         if (fZ[p[1]]>z2)  {z2=fZ[p[1]]; x2=fX[p[1]]; y2=fY[p[1]]; i2=1;}
         if (fZ[p[2]]<=z0) {z0=fZ[p[2]]; x0=fX[p[2]]; y0=fY[p[2]]; i0=2;}
         if (fZ[p[2]]>z2)  {z2=fZ[p[2]]; x2=fX[p[2]]; y2=fY[p[2]]; i2=2;}
         if (i0==0 && i2==0) {
            Error("GetContourList", "wrong vertices ordering");
            delete [] xs0;
            delete [] ys0;
            delete [] xs1;
            delete [] ys1;
            return nullptr;
         } else {
            i1 = 3-i2-i0;
         }
         x1 = fX[p[i1]];
         y1 = fY[p[i1]];
         z1 = fZ[p[i1]];

         if (Hoption.Logz) {
            z0 = TMath::Log10(z0);
            z1 = TMath::Log10(z1);
            z2 = TMath::Log10(z2);
         }

         if(contour >= z0 && contour <=z2) {
            r20 = (contour-z0)/(z2-z0);
            xs0c = r20*(x2-x0)+x0;
            ys0c = r20*(y2-y0)+y0;
            if(contour >= z1 && contour <=z2) {
               r21 = (contour-z1)/(z2-z1);
               xs1c = r21*(x2-x1)+x1;
               ys1c = r21*(y2-y1)+y1;
            } else {
               r10 = (contour-z0)/(z1-z0);
               xs1c = r10*(x1-x0)+x0;
               ys1c = r10*(y1-y0)+y0;
            }
            // do not take the segments equal to a point
            if(xs0c != xs1c || ys0c != ys1c) {
               nbSeg++;
               xs0[nbSeg-1] = xs0c;
               ys0[nbSeg-1] = ys0c;
               xs1[nbSeg-1] = xs1c;
               ys1[nbSeg-1] = ys1c;
            }
         }
      }
   }

   TList *list   = new TList(); // list holding all the graphs

   Int_t *segUsed = new Int_t[fNdt];
   for(i=0; i<fNdt; i++) segUsed[i]=kFALSE;

   // Find all the graphs making the contour. There is two kind of graphs,
   // either they are "opened" or they are "closed"

   // Find the opened graphs
   Double_t xc=0, yc=0, xnc=0, ync=0;
   Bool_t findNew;
   Bool_t s0, s1;
   Int_t is, js;
   for (is=0; is<nbSeg; is++) {
      if (segUsed[is]) continue;
      s0 = s1 = kFALSE;

      // Find to which segment is is connected. It can be connected
      // via 0, 1 or 2 vertices.
      for (js=0; js<nbSeg; js++) {
         if (is==js) continue;
         if (xs0[is]==xs0[js] && ys0[is]==ys0[js]) s0 = kTRUE;
         if (xs0[is]==xs1[js] && ys0[is]==ys1[js]) s0 = kTRUE;
         if (xs1[is]==xs0[js] && ys1[is]==ys0[js]) s1 = kTRUE;
         if (xs1[is]==xs1[js] && ys1[is]==ys1[js]) s1 = kTRUE;
      }

      // Segment is is alone, not connected. It is stored in the
      // list and the next segment is examined.
      if (!s0 && !s1) {
         graph = new TGraph();
         graph->SetPoint(npg,xs0[is],ys0[is]); npg++;
         graph->SetPoint(npg,xs1[is],ys1[is]); npg++;
         segUsed[is] = kTRUE;
         list->Add(graph); npg = 0;
         continue;
      }

      // Segment is is connected via 1 vertex only and can be considered
      // as the starting point of an opened contour.
      if (!s0 || !s1) {
         // Find all the segments connected to segment is
         graph = new TGraph();
         if (s0) {xc = xs0[is]; yc = ys0[is]; xnc = xs1[is]; ync = ys1[is];}
         if (s1) {xc = xs1[is]; yc = ys1[is]; xnc = xs0[is]; ync = ys0[is];}
         graph->SetPoint(npg,xnc,ync); npg++;
         segUsed[is] = kTRUE;
         js = 0;
L01:
         findNew = kFALSE;
         if (js < nbSeg && segUsed[js]) {
            js++;
            goto L01;
         } else if (xc==xs0[js] && yc==ys0[js]) {
            xc      = xs1[js];
            yc      = ys1[js];
            findNew = kTRUE;
         } else if (xc==xs1[js] && yc==ys1[js]) {
            xc      = xs0[js];
            yc      = ys0[js];
            findNew = kTRUE;
         }
         if (findNew) {
            segUsed[js] = kTRUE;
            graph->SetPoint(npg,xc,yc); npg++;
            js = 0;
            goto L01;
         }
         js++;
         if (js<nbSeg) goto L01;
         list->Add(graph); npg = 0;
      }
   }


   // Find the closed graphs. At this point all the remaining graphs
   // are closed. Any segment can be used to start the search.
   for (is=0; is<nbSeg; is++) {
      if (segUsed[is]) continue;

      // Find all the segments connected to segment is
      graph = new TGraph();
      segUsed[is] = kTRUE;
      xc = xs0[is];
      yc = ys0[is];
      js = 0;
      graph->SetPoint(npg,xc,yc); npg++;
L02:
      findNew = kFALSE;
      if (js < nbSeg && segUsed[js]) {
         js++;
         goto L02;
      } else if (xc==xs0[js] && yc==ys0[js]) {
         xc      = xs1[js];
         yc      = ys1[js];
         findNew = kTRUE;
      } else if (xc==xs1[js] && yc==ys1[js]) {
         xc      = xs0[js];
         yc      = ys0[js];
         findNew = kTRUE;
      }
      if (findNew) {
         segUsed[js] = kTRUE;
         graph->SetPoint(npg,xc,yc); npg++;
         js = 0;
         goto L02;
      }
      js++;
      if (js<nbSeg) goto L02;
      // Close the contour
      graph->SetPoint(npg,xs0[is],ys0[is]); npg++;
      list->Add(graph); npg = 0;
   }

   delete [] xs0;
   delete [] ys0;
   delete [] xs1;
   delete [] ys1;
   delete [] segUsed;
   return list;
}


////////////////////////////////////////////////////////////////////////////////
/// Paint a TGraphDelaunay according to the value of "option":
///
/// | Option   | Description                                                   |
/// |----------|---------------------------------------------------------------|
/// | "TRI"    | The Delaunay triangles are drawn using filled area. An hidden surface drawing technique is used. The surface is painted with the current fill area color. The edges of each triangles are painted with the current line color. |
/// | "TRIW"   | The Delaunay triangles are drawn as wire frame. |
/// | "TRI1"   | The Delaunay triangles are painted with color levels. The edges of each triangles are painted with the current line color. |
/// | "TRI2"   | the Delaunay triangles are painted with color levels. |
/// | "P"      | Draw a marker at each vertex. |
/// | "P0"     | Draw a circle at each vertex. Each circle background is white. |
/// | "PCOL"   | Draw a marker at each vertex. The color of each marker is defined according to its Z position. |
/// | "CONT"   | Draw contours. |
/// | "LINE"   | Draw a 3D polyline. |

void TGraph2DPainter::Paint(Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   Bool_t triangles = opt.Contains("tri")  ||
                      opt.Contains("tri1") ||
                      opt.Contains("tri2");
   if (opt.Contains("tri0")) triangles = kFALSE;

   Bool_t markers   = opt.Contains("p") && !triangles;
   Bool_t contour   = opt.Contains("cont");
   Bool_t line      = opt.Contains("line");
   Bool_t err       = opt.Contains("err");

   fGraph2D->TAttLine::Modify();
   fGraph2D->TAttFill::Modify();
   fGraph2D->TAttMarker::Modify();

   // Compute minimums and maximums
   TAxis *xaxis = gCurrentHist->GetXaxis();
   Int_t first = xaxis->GetFirst();
   fXmin = xaxis->GetBinLowEdge(first);
   if (Hoption.Logx && fXmin <= 0) fXmin = xaxis->GetBinUpEdge(xaxis->FindFixBin(0.01*xaxis->GetBinWidth(first)));
   fXmax = xaxis->GetBinUpEdge(xaxis->GetLast());
   TAxis *yaxis = gCurrentHist->GetYaxis();
   first = yaxis->GetFirst();
   fYmin = yaxis->GetBinLowEdge(first);
   if (Hoption.Logy && fYmin <= 0) fYmin = yaxis->GetBinUpEdge(yaxis->FindFixBin(0.01*yaxis->GetBinWidth(first)));
   fYmax = yaxis->GetBinUpEdge(yaxis->GetLast());
   fZmax = fGraph2D->GetZmaxE();
   fZmin = fGraph2D->GetZminE();
   if (Hoption.Logz && fZmin <= 0) fZmin = TMath::Min((Double_t)1, (Double_t)0.001*fGraph2D->GetZmax());

   if (triangles) PaintTriangles(option);
   if (contour)   PaintContour(option);
   if (line)      PaintPolyLine(option);
   if (err)       PaintErrors(option);
   if (markers)   PaintPolyMarker(option);
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as a contour plot. Delaunay triangles are used
/// to compute the contours.

void TGraph2DPainter::PaintContour(Option_t * /*option*/)
{
   // Initialize the levels on the Z axis
   Int_t ncolors  = gStyle->GetNumberOfColors();
   Int_t ndiv   = gCurrentHist->GetContour();
   if (ndiv == 0 ) {
      ndiv = gStyle->GetNumberContours();
      gCurrentHist->SetContour(ndiv);
   }
   Int_t ndivz  = TMath::Abs(ndiv);
   if (gCurrentHist->TestBit(TH1::kUserContour) == 0) gCurrentHist->SetContour(ndiv);

   Int_t theColor;
   TList *l;
   TGraph *g;
   TObject *obj;
   Double_t c;

   if (!fNdt) FindTriangles();

   for (Int_t k=0; k<ndiv; k++) {
      c = gCurrentHist->GetContourLevelPad(k);
      l = GetContourList(c);
      TIter next(l);
      while ((obj = next())) {
         if(obj->InheritsFrom(TGraph::Class()) ) {
            g=(TGraph*)obj;
            g->SetLineWidth(fGraph2D->GetLineWidth());
            g->SetLineStyle(fGraph2D->GetLineStyle());
            theColor = Int_t((k+0.99)*Float_t(ncolors)/Float_t(ndivz));
            g->SetLineColor(gStyle->GetColorPalette(theColor));
            g->Paint("l");
         }
      }
      if (l) { l->Delete(); delete l; }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as error bars

void TGraph2DPainter::PaintErrors(Option_t * /* option */)
{
   Double_t temp1[3],temp2[3];

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintErrors", "No TView in current pad");
      return;
   }

   Int_t  it;

   Double_t *xm = new Double_t[2];
   Double_t *ym = new Double_t[2];

   fGraph2D->SetLineStyle(fGraph2D->GetLineStyle());
   fGraph2D->SetLineWidth(fGraph2D->GetLineWidth());
   fGraph2D->SetLineColor(fGraph2D->GetLineColor());
   fGraph2D->TAttLine::Modify();

   for (it=0; it<fNpoints; it++) {
      if(fX[it] < fXmin || fX[it] > fXmax) continue;
      if(fY[it] < fYmin || fY[it] > fYmax) continue;
      if (fEX) {
         temp1[0] = fX[it]-fEX[it];
         temp1[1] = fY[it];
         temp1[2] = fZ[it];
         temp1[0] = TMath::Max(temp1[0],fXmin);
         temp1[1] = TMath::Max(temp1[1],fYmin);
         temp1[2] = TMath::Max(temp1[2],fZmin);
         temp1[2] = TMath::Min(temp1[2],fZmax);
         if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
         if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
         if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
         view->WCtoNDC(temp1, &temp2[0]);
         xm[0] = temp2[0];
         ym[0] = temp2[1];

         temp1[0] = fX[it]+fEX[it];
         temp1[0] = TMath::Max(temp1[0],fXmin);
         if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
         view->WCtoNDC(temp1, &temp2[0]);
         xm[1] = temp2[0];
         ym[1] = temp2[1];
         gPad->PaintPolyLine(2,xm,ym);
      }
      if (fEY) {
         temp1[0] = fX[it];
         temp1[1] = fY[it]-fEY[it];
         temp1[2] = fZ[it];
         temp1[0] = TMath::Max(temp1[0],fXmin);
         temp1[1] = TMath::Max(temp1[1],fYmin);
         temp1[2] = TMath::Max(temp1[2],fZmin);
         temp1[2] = TMath::Min(temp1[2],fZmax);
         if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
         if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
         if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
         view->WCtoNDC(temp1, &temp2[0]);
         xm[0] = temp2[0];
         ym[0] = temp2[1];

         temp1[1] = fY[it]+fEY[it];
         temp1[1] = TMath::Max(temp1[1],fYmin);
         if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
         view->WCtoNDC(temp1, &temp2[0]);
         xm[1] = temp2[0];
         ym[1] = temp2[1];
         gPad->PaintPolyLine(2,xm,ym);
      }
      if (fEZ) {
         temp1[0] = fX[it];
         temp1[1] = fY[it];
         temp1[2] = fZ[it]-fEZ[it];
         temp1[0] = TMath::Max(temp1[0],fXmin);
         temp1[1] = TMath::Max(temp1[1],fYmin);
         temp1[2] = TMath::Max(temp1[2],fZmin);
         temp1[2] = TMath::Min(temp1[2],fZmax);
         if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
         if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
         if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
         view->WCtoNDC(temp1, &temp2[0]);
         xm[0] = temp2[0];
         ym[0] = temp2[1];

         temp1[2] = fZ[it]+fEZ[it];
         temp1[2] = TMath::Max(temp1[2],fZmin);
         temp1[2] = TMath::Min(temp1[2],fZmax);
         if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
         view->WCtoNDC(temp1, &temp2[0]);
         xm[1] = temp2[0];
         ym[1] = temp2[1];
         gPad->PaintPolyLine(2,xm,ym);
      }
   }
   delete [] xm;
   delete [] ym;
}


////////////////////////////////////////////////////////////////////////////////
/// Paints one triangle.
///
/// - nblev  = 0 : paint the color levels
/// - nblev != 0 : paint the grid

void TGraph2DPainter::PaintLevels(Int_t *t,Double_t *x, Double_t *y,
                           Int_t nblev, Double_t *glev)
{
   Int_t i, fillColor, ncolors, theColor0, theColor2;

   Int_t p[3];
   if (fDelaunay) {
      p[0]=t[0]-1;
      p[1]=t[1]-1;
      p[2]=t[2]-1;
   }
   else  {
      p[0]=t[0];
      p[1]=t[1];
      p[2]=t[2];
   }
   Double_t xl[2],yl[2];
   Double_t zl, r21, r20, r10;
   Double_t x0 = x[0]  , x2 = x[0];
   Double_t y0 = y[0]  , y2 = y[0];
   Double_t z0 = fZ[p[0]], z2 = fZ[p[0]];
   Double_t zmin = fGraph2D->GetMinimum();
   Double_t zmax = fGraph2D->GetMaximum();
   if (zmin==-1111 && zmax==-1111) {
      zmin = TMath::Min(fZmin, 0.);
      if (Hoption.Logz && zmin <= 0) zmin = TMath::Min((Double_t)1, (Double_t)0.001*fGraph2D->GetZmax());
      zmax = fZmax;
   }

   // Order along Z axis the points (xi,yi,zi) where "i" belongs to {0,1,2}
   // After this z0 < z1 < z2
   Int_t i0=0, i1=0, i2=0;
   if (fZ[p[1]]<=z0) {z0=fZ[p[1]]; x0=x[1]; y0=y[1]; i0=1;}
   if (fZ[p[1]]>z2)  {z2=fZ[p[1]]; x2=x[1]; y2=y[1]; i2=1;}
   if (fZ[p[2]]<=z0) {z0=fZ[p[2]]; x0=x[2]; y0=y[2]; i0=2;}
   if (fZ[p[2]]>z2)  {z2=fZ[p[2]]; x2=x[2]; y2=y[2]; i2=2;}
   i1 = 3-i2-i0;
   Double_t x1 = x[i1];
   Double_t y1 = y[i1];
   Double_t z1 = fZ[p[i1]];

   if (z0>zmax) z0 = zmax;
   if (z2>zmax) z2 = zmax;
   if (z0<zmin) z0 = zmin;
   if (z2<zmin) z2 = zmin;
   if (z1>zmax) z1 = zmax;
   if (z1<zmin) z1 = zmin;

   if (Hoption.Logz) {
      z0   = TMath::Log10(z0);
      z1   = TMath::Log10(z1);
      z2   = TMath::Log10(z2);
      zmin = TMath::Log10(zmin);
      zmax = TMath::Log10(zmax);
   }

   // zi  = Z values of the stripe number i
   // zip = Previous zi
   Double_t zi=0, zip=0;

   if (nblev <= 0) {
      // Paint the colors levels

      // Compute the color associated to z0 (theColor0) and z2 (theColor2)
      ncolors   = gStyle->GetNumberOfColors();
      theColor0 = (Int_t)( ((z0-zmin)/(zmax-zmin))*(ncolors-1) );
      theColor2 = (Int_t)( ((z2-zmin)/(zmax-zmin))*(ncolors-1) );

      // The stripes drawn to fill the triangles may have up to 5 points
      Double_t xp[5], yp[5];

      // rl = Ratio between z0 and z2 (long)
      // rs = Ratio between z0 and z1 or z1 and z2 (short)
      Double_t rl,rs;

      // ci = Color of the stripe number i
      // npf = number of point needed to draw the current stripe
      Int_t ci,npf;

      fillColor = fGraph2D->GetFillColor();

      // If the z0's color and z2's colors are the same, the whole triangle
      // can be painted in one go.
      if(theColor0 == theColor2) {
         fGraph2D->SetFillColor(gStyle->GetColorPalette(theColor0));
         fGraph2D->TAttFill::Modify();
         gPad->PaintFillArea(3,x,y);

      // The triangle must be painted with several colors
      } else {
         for(ci=theColor0; ci<=theColor2; ci++) {
            fGraph2D->SetFillColor(gStyle->GetColorPalette(ci));
            fGraph2D->TAttFill::Modify();
            if (ci==theColor0) {
               zi    = (((ci+1)*(zmax-zmin))/(ncolors-1))+zmin;
               xp[0] = x0;
               yp[0] = y0;
               rl    = (zi-z0)/(z2-z0);
               xp[1] = rl*(x2-x0)+x0;
               yp[1] = rl*(y2-y0)+y0;
               if (zi>=z1 || z0==z1) {
                  rs    = (zi-z1)/(z2-z1);
                  xp[2] = rs*(x2-x1)+x1;
                  yp[2] = rs*(y2-y1)+y1;
                  xp[3] = x1;
                  yp[3] = y1;
                  npf   = 4;
               } else {
                  rs    = (zi-z0)/(z1-z0);
                  xp[2] = rs*(x1-x0)+x0;
                  yp[2] = rs*(y1-y0)+y0;
                  npf   = 3;
               }
            } else if (ci==theColor2) {
               xp[0] = xp[1];
               yp[0] = yp[1];
               xp[1] = x2;
               yp[1] = y2;
               if (zi<z1 || z2==z1) {
                  xp[3] = xp[2];
                  yp[3] = yp[2];
                  xp[2] = x1;
                  yp[2] = y1;
                  npf   = 4;
               } else {
                  npf   = 3;
               }
            } else {
               zi    = (((ci+1)*(zmax-zmin))/(ncolors-1))+zmin;
               xp[0] = xp[1];
               yp[0] = yp[1];
               rl    = (zi-z0)/(z2-z0);
               xp[1] = rl*(x2-x0)+x0;
               yp[1] = rl*(y2-y0)+y0;
               if ( zi>=z1 && zip<=z1) {
                  xp[3] = x1;
                  yp[3] = y1;
                  xp[4] = xp[2];
                  yp[4] = yp[2];
                  npf   = 5;
               } else {
                  xp[3] = xp[2];
                  yp[3] = yp[2];
                  npf   = 4;
               }
               if (zi<z1) {
                  rs    = (zi-z0)/(z1-z0);
                  xp[2] = rs*(x1-x0)+x0;
                  yp[2] = rs*(y1-y0)+y0;
               } else {
                  rs    = (zi-z1)/(z2-z1);
                  xp[2] = rs*(x2-x1)+x1;
                  yp[2] = rs*(y2-y1)+y1;
               }
            }
            zip = zi;
            // Paint a stripe
            gPad->PaintFillArea(npf,xp,yp);
         }
      }
      fGraph2D->SetFillColor(fillColor);
      fGraph2D->TAttFill::Modify();

   } else {
      // Paint the grid levels
      fGraph2D->SetLineStyle(3);
      fGraph2D->TAttLine::Modify();
      for(i=0; i<nblev; i++){
         zl=glev[i];
         if(zl >= z0 && zl <=z2) {
            r21=(zl-z1)/(z2-z1);
            r20=(zl-z0)/(z2-z0);
            r10=(zl-z0)/(z1-z0);
            xl[0]=r20*(x2-x0)+x0;
            yl[0]=r20*(y2-y0)+y0;
            if(zl >= z1 && zl <=z2) {
               xl[1]=r21*(x2-x1)+x1;
               yl[1]=r21*(y2-y1)+y1;
            } else {
               xl[1]=r10*(x1-x0)+x0;
               yl[1]=r10*(y1-y0)+y0;
            }
            gPad->PaintPolyLine(2,xl,yl);
         }
      }
      fGraph2D->SetLineStyle(1);
      fGraph2D->TAttLine::Modify();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as PaintPolyMarker

void TGraph2DPainter::PaintPolyMarker(Option_t *option)
{
   Double_t temp1[3],temp2[3];

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintPolyMarker", "No TView in current pad");
      return;
   }

   TString opt = option;
   opt.ToLower();
   Bool_t markers0 = opt.Contains("p0");
   Bool_t colors   = opt.Contains("pcol");
   Int_t  ncolors  = gStyle->GetNumberOfColors();
   Int_t  it, theColor;

   // Initialize the levels on the Z axis
   if (colors) {
      Int_t ndiv   = gCurrentHist->GetContour();
      if (ndiv == 0 ) {
         ndiv = gStyle->GetNumberContours();
         gCurrentHist->SetContour(ndiv);
      }
      if (gCurrentHist->TestBit(TH1::kUserContour) == 0) gCurrentHist->SetContour(ndiv);
   }

   Double_t *xm = new Double_t[fNpoints];
   Double_t *ym = new Double_t[fNpoints];
   Double_t *zm = new Double_t[fNpoints];
   Double_t hzmin = gCurrentHist->GetMinimum();
   Double_t hzmax = gCurrentHist->GetMaximum();

   // min and max for colors
   Double_t hzmincol = hzmin;
   Double_t hzmaxcol = hzmax;
   if (hzmincol==-1111 && hzmaxcol==-1111) {
      hzmincol = TMath::Min(hzmincol, 0.);
      if (Hoption.Logz && hzmincol <= 0) hzmincol = TMath::Min((Double_t)1, (Double_t)0.001*fGraph2D->GetZmax());
      hzmaxcol = fZmax;
   }
   if (Hoption.Logz) {
      hzmincol = TMath::Log10(hzmincol);
      hzmaxcol = TMath::Log10(hzmaxcol);
   }

   Double_t Xeps = (fXmax-fXmin)*0.0001;
   Double_t Yeps = (fYmax-fYmin)*0.0001;
   Double_t Zeps = (hzmax-hzmin)*0.0001;

   Int_t    npd = 0;
   for (it=0; it<fNpoints; it++) {
      xm[it] = 0;
      ym[it] = 0;
      if(fXmin - fX[it] > Xeps || fX[it] - fXmax > Xeps) continue;
      if(fYmin - fY[it] > Yeps || fY[it] - fYmax > Yeps) continue;
      if(hzmin - fZ[it] > Zeps || fZ[it] - hzmax > Zeps) continue;
      temp1[0] = fX[it];
      temp1[1] = fY[it];
      temp1[2] = fZ[it];
      temp1[0] = TMath::Max(temp1[0],fXmin);
      temp1[1] = TMath::Max(temp1[1],fYmin);
      temp1[2] = TMath::Max(temp1[2],hzmin);
      temp1[2] = TMath::Min(temp1[2],hzmax);
      if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
      if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
      if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
      view->WCtoNDC(temp1, &temp2[0]);
      xm[npd] = temp2[0];
      ym[npd] = temp2[1];
      zm[npd] = temp1[2];
      npd++;
   }
   if (markers0) {
      PaintPolyMarker0(npd,xm,ym);
   } else if (colors) {
      Int_t cols = fGraph2D->GetMarkerColor();
      for (it=0; it<npd; it++) {
         theColor = (Int_t)( ((zm[it]-hzmincol)/(hzmaxcol-hzmincol))*(ncolors-1) );
         fGraph2D->SetMarkerColor(gStyle->GetColorPalette(theColor));
         fGraph2D->TAttMarker::Modify();
         gPad->PaintPolyMarker(1,&xm[it],&ym[it]);
      }
      fGraph2D->SetMarkerColor(cols);
   } else {
      fGraph2D->SetMarkerStyle(fGraph2D->GetMarkerStyle());
      fGraph2D->SetMarkerSize(fGraph2D->GetMarkerSize());
      fGraph2D->SetMarkerColor(fGraph2D->GetMarkerColor());
      fGraph2D->TAttMarker::Modify();
      gPad->PaintPolyMarker(npd,xm,ym);
   }
   delete [] xm;
   delete [] ym;
   delete [] zm;
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as PaintPolyLine

void TGraph2DPainter::PaintPolyLine(Option_t * /* option */)
{
   Double_t temp1[3],temp2[3];

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintPolyLine", "No TView in current pad");
      return;
   }

   Int_t  it;
   Double_t Xeps = (fXmax - fXmin) * 0.0001;
   Double_t Yeps = (fYmax - fYmin) * 0.0001;

   Double_t *xm = new Double_t[fNpoints];
   Double_t *ym = new Double_t[fNpoints];
   Int_t    npd = 0;

   for (it=0; it<fNpoints; it++) {
      if (fXmin - fX[it] > Xeps || fX[it] - fXmax > Xeps)
         continue;
      if (fYmin - fY[it] > Yeps || fY[it] - fYmax > Yeps)
         continue;
      npd++;
      temp1[0] = fX[it];
      temp1[1] = fY[it];
      temp1[2] = fZ[it];
      temp1[0] = TMath::Max(temp1[0],fXmin);
      temp1[1] = TMath::Max(temp1[1],fYmin);
      temp1[2] = TMath::Max(temp1[2],fZmin);
      temp1[2] = TMath::Min(temp1[2],fZmax);
      if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
      if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
      if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
      view->WCtoNDC(temp1, &temp2[0]);
      xm[npd - 1] = temp2[0];
      ym[npd - 1] = temp2[1];
   }
   fGraph2D->SetLineStyle(fGraph2D->GetLineStyle());
   fGraph2D->SetLineWidth(fGraph2D->GetLineWidth());
   fGraph2D->SetLineColor(fGraph2D->GetLineColor());
   fGraph2D->TAttLine::Modify();
   gPad->PaintPolyLine(npd,xm,ym);
   delete [] xm;
   delete [] ym;
}


////////////////////////////////////////////////////////////////////////////////
/// Paints a circle at each vertex. Each circle background is white.

void TGraph2DPainter::PaintPolyMarker0(Int_t n, Double_t *x, Double_t *y)
{
   fGraph2D->SetMarkerSize(fGraph2D->GetMarkerSize());
   Int_t mc = fGraph2D->GetMarkerColor();
   Int_t ms = fGraph2D->GetMarkerStyle();
   for (Int_t i=0; i<n; i++) {
      fGraph2D->SetMarkerStyle(20);
      fGraph2D->SetMarkerColor(0);
      fGraph2D->TAttMarker::Modify();
      gPad->PaintPolyMarker(1,&x[i],&y[i]);
      fGraph2D->SetMarkerStyle(24);
      fGraph2D->SetMarkerColor(mc);
      fGraph2D->TAttMarker::Modify();
      gPad->PaintPolyMarker(1,&x[i],&y[i]);
   }
   fGraph2D->SetMarkerStyle(ms);
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as triangles

void TGraph2DPainter::PaintTriangles(Option_t *option)
{
   if (fDelaunay)
      PaintTriangles_old(option);
   else if (fDelaunay2D)
      PaintTriangles_new(option);
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as triangles (old implementation)

void TGraph2DPainter::PaintTriangles_old(Option_t *option)
{
   Double_t x[4], y[4], temp1[3],temp2[3];
   Int_t it,t[3];
   Int_t *order = 0;
   Double_t *dist = 0;

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintTriangles", "No TView in current pad");
      return;
   }

   TString opt = option;
   opt.ToLower();
   Bool_t tri1      = opt.Contains("tri1");
   Bool_t tri2      = opt.Contains("tri2");
   Bool_t markers   = opt.Contains("p");
   Bool_t markers0  = opt.Contains("p0");
   Bool_t wire      = opt.Contains("w");

   // Define the grid levels drawn on the triangles.
   // The grid levels are aligned on the Z axis' main tick marks.
   Int_t nblev=0;
   Double_t *glev=0;
   if (!tri1 && !tri2 && !wire) {
      Int_t ndivz = gCurrentHist->GetZaxis()->GetNdivisions()%100;
      Int_t nbins;
      Double_t binLow = 0, binHigh = 0, binWidth = 0;

      // Find the main tick marks positions.
      Double_t *r0 = view->GetRmin();
      Double_t *r1 = view->GetRmax();
      if (!r0 || !r1) return;

      if (ndivz > 0) {
         THLimitsFinder::Optimize(r0[2], r1[2], ndivz,
                                  binLow, binHigh, nbins, binWidth, " ");
      } else {
         nbins = TMath::Abs(ndivz);
         binLow = r0[2];
         binHigh = r1[2];
         binWidth = (binHigh-binLow)/nbins;
      }
      // Define the grid levels
      nblev = nbins+1;
      glev = new Double_t[nblev];
      for (Int_t i = 0; i < nblev; ++i) glev[i] = binLow+i*binWidth;
   }

   // Initialize the levels on the Z axis
   if (tri1 || tri2) {
      Int_t ndiv   = gCurrentHist->GetContour();
      if (ndiv == 0 ) {
         ndiv = gStyle->GetNumberContours();
         gCurrentHist->SetContour(ndiv);
      }
      if (gCurrentHist->TestBit(TH1::kUserContour) == 0) gCurrentHist->SetContour(ndiv);
   }

   // For each triangle, compute the distance between the triangle centre
   // and the back planes. Then these distances are sorted in order to draw
   // the triangles from back to front.
   if (!fNdt) FindTriangles();
   Double_t cp = TMath::Cos(view->GetLongitude()*TMath::Pi()/180.);
   Double_t sp = TMath::Sin(view->GetLongitude()*TMath::Pi()/180.);
   order = new Int_t[fNdt];
   dist  = new Double_t[fNdt];
   Double_t xd,yd;
   Int_t p, n, m;
   Bool_t o = kFALSE;
   for (it=0; it<fNdt; it++) {
      p = fPTried[it];
      n = fNTried[it];
      m = fMTried[it];
      xd = (fXN[p]+fXN[n]+fXN[m])/3;
      yd = (fYN[p]+fYN[n]+fYN[m])/3;
      if ((cp >= 0) && (sp >= 0.)) {
         dist[it] = -(fXNmax-xd+fYNmax-yd);
      } else if ((cp <= 0) && (sp >= 0.)) {
         dist[it] = -(fXNmax-xd+yd-fYNmin);
         o = kTRUE;
      } else if ((cp <= 0) && (sp <= 0.)) {
         dist[it] = -(xd-fXNmin+yd-fYNmin);
      } else {
         dist[it] = -(xd-fXNmin+fYNmax-yd);
         o = kTRUE;
      }
   }
   TMath::Sort(fNdt, dist, order, o);

   // Draw the triangles and markers if requested
   fGraph2D->SetFillColor(fGraph2D->GetFillColor());
   Int_t fs = fGraph2D->GetFillStyle();
   fGraph2D->SetFillStyle(1001);
   fGraph2D->TAttFill::Modify();
   fGraph2D->SetLineColor(fGraph2D->GetLineColor());
   fGraph2D->TAttLine::Modify();
   int lst = fGraph2D->GetLineStyle();
   for (it=0; it<fNdt; it++) {
      t[0] = fPTried[order[it]];
      t[1] = fNTried[order[it]];
      t[2] = fMTried[order[it]];
      for (Int_t k=0; k<3; k++) {
         if(fX[t[k]-1] < fXmin || fX[t[k]-1] > fXmax) goto endloop;
         if(fY[t[k]-1] < fYmin || fY[t[k]-1] > fYmax) goto endloop;
         temp1[0] = fX[t[k]-1];
         temp1[1] = fY[t[k]-1];
         temp1[2] = fZ[t[k]-1];
         temp1[0] = TMath::Max(temp1[0],fXmin);
         temp1[1] = TMath::Max(temp1[1],fYmin);
         temp1[2] = TMath::Max(temp1[2],fZmin);
         temp1[2] = TMath::Min(temp1[2],fZmax);
         if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
         if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
         if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
         view->WCtoNDC(temp1, &temp2[0]);
         x[k] = temp2[0];
         y[k] = temp2[1];
      }
      x[3] = x[0];
      y[3] = y[0];
      if (tri1 || tri2) PaintLevels(t,x,y);
      if (!tri1 && !tri2 && !wire) {
         gPad->PaintFillArea(3,x,y);
         PaintLevels(t,x,y,nblev,glev);
      }
      if (!tri2) gPad->PaintPolyLine(4,x,y);
      if (markers) {
         if (markers0) {
            PaintPolyMarker0(3,x,y);
         } else {
            fGraph2D->SetMarkerStyle(fGraph2D->GetMarkerStyle());
            fGraph2D->SetMarkerSize(fGraph2D->GetMarkerSize());
            fGraph2D->SetMarkerColor(fGraph2D->GetMarkerColor());
            fGraph2D->TAttMarker::Modify();
            gPad->PaintPolyMarker(3,x,y);
         }
      }
endloop:
      continue;
   }
   fGraph2D->SetFillStyle(fs);
   fGraph2D->SetLineStyle(lst);
   fGraph2D->TAttLine::Modify();
   fGraph2D->TAttFill::Modify();
   delete [] order;
   delete [] dist;
   if (glev) delete [] glev;
}


////////////////////////////////////////////////////////////////////////////////
/// Paints the 2D graph as triangles (new implementation)

void TGraph2DPainter::PaintTriangles_new(Option_t *option)
{

   Double_t x[4], y[4], temp1[3],temp2[3];
   Int_t p[3];

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintTriangles", "No TView in current pad");
      return;
   }

   TString opt = option;
   opt.ToLower();
   Bool_t tri1      = opt.Contains("tri1");
   Bool_t tri2      = opt.Contains("tri2");
   Bool_t markers   = opt.Contains("p");
   Bool_t markers0  = opt.Contains("p0");
   Bool_t wire      = opt.Contains("w");

   // Define the grid levels drawn on the triangles.
   // The grid levels are aligned on the Z axis' main tick marks.
   Int_t nblev=0;
   Double_t *glev=0;
   if (!tri1 && !tri2 && !wire) {
      Int_t ndivz = gCurrentHist->GetZaxis()->GetNdivisions()%100;
      Int_t nbins;
      Double_t binLow = 0, binHigh = 0, binWidth = 0;

      // Find the main tick marks positions.
      Double_t *r0 = view->GetRmin();
      Double_t *r1 = view->GetRmax();
      if (!r0 || !r1) return;

      if (ndivz > 0) {
         THLimitsFinder::Optimize(r0[2], r1[2], ndivz,
                                  binLow, binHigh, nbins, binWidth, " ");
      } else {
         nbins = TMath::Abs(ndivz);
         binLow = r0[2];
         binHigh = r1[2];
         binWidth = (binHigh-binLow)/nbins;
      }
      // Define the grid levels
      nblev = nbins+1;
      glev = new Double_t[nblev];
      for (Int_t i = 0; i < nblev; ++i) glev[i] = binLow+i*binWidth;
   }

   // Initialize the levels on the Z axis
   if (tri1 || tri2) {
      Int_t ndiv   = gCurrentHist->GetContour();
      if (ndiv == 0 ) {
         ndiv = gStyle->GetNumberContours();
         gCurrentHist->SetContour(ndiv);
      }
      if (gCurrentHist->TestBit(TH1::kUserContour) == 0) gCurrentHist->SetContour(ndiv);
   }

   // For each triangle, compute the distance between the triangle centre
   // and the back planes. Then these distances are sorted in order to draw
   // the triangles from back to front.
   if (!fNdt) FindTriangles();
   Double_t cp = TMath::Cos(view->GetLongitude()*TMath::Pi()/180.);
   Double_t sp = TMath::Sin(view->GetLongitude()*TMath::Pi()/180.);

   Bool_t reverse = kFALSE;
   std::function<Double_t(Double_t, Double_t)> fDist;

   if ((cp >= 0) && (sp >= 0.)) {
      fDist = [&](Double_t xd, Double_t yd) -> Double_t { return -(fXNmax-xd+fYNmax-yd);};
   } else if ((cp <= 0) && (sp >= 0.)) {
      fDist = [&](Double_t xd, Double_t yd) -> Double_t { return -(fXNmax-xd+yd-fYNmin);};
      reverse = kTRUE;
   } else if ((cp <= 0) && (sp <= 0.)) {
      fDist = [&](Double_t xd, Double_t yd) -> Double_t { return -(xd-fXNmin+yd-fYNmin);};
   } else {
      fDist = [&](Double_t xd, Double_t yd) -> Double_t { return -(xd-fXNmin+fYNmax-yd);};
      reverse = kTRUE;
   }

   typedef std::pair<Double_t, TGraphDelaunay2D::Triangles::const_iterator> DistEntry;
   std::vector<DistEntry> dist;
   for(auto it = fDelaunay2D->begin(); it != fDelaunay2D->end(); ++it){
      auto face = *it;
      Double_t xd = (face.x[0] + face.x[1] + face.x[2]) / 3;
      Double_t yd = (face.y[0] + face.y[1] + face.y[2]) / 3;

      dist.emplace_back(fDist(xd, yd), it);
   }

   std::sort(dist.begin(), dist.end(),
             [&](const DistEntry & a, const DistEntry & b){ return !reverse ? (a.first < b.first) : (b.first < a .first); });

   // Draw the triangles and markers if requested
   fGraph2D->SetFillColor(fGraph2D->GetFillColor());
   Int_t fs = fGraph2D->GetFillStyle();
   fGraph2D->SetFillStyle(1001);
   fGraph2D->TAttFill::Modify();
   fGraph2D->SetLineColor(fGraph2D->GetLineColor());
   fGraph2D->TAttLine::Modify();
   int lst = fGraph2D->GetLineStyle();
   for (const auto & it : dist) {
      p[0] = it.second->idx[0];
      p[1] = it.second->idx[1];
      p[2] = it.second->idx[2];
      for (Int_t k=0; k<3; k++) {
         if(fX[p[k]] < fXmin || fX[p[k]] > fXmax) goto endloop;
         if(fY[p[k]] < fYmin || fY[p[k]] > fYmax) goto endloop;
         temp1[0] = fX[p[k]];
         temp1[1] = fY[p[k]];
         temp1[2] = fZ[p[k]];
         temp1[0] = TMath::Max(temp1[0],fXmin);
         temp1[1] = TMath::Max(temp1[1],fYmin);
         temp1[2] = TMath::Max(temp1[2],fZmin);
         temp1[2] = TMath::Min(temp1[2],fZmax);
         if (Hoption.Logx) temp1[0] = TMath::Log10(temp1[0]);
         if (Hoption.Logy) temp1[1] = TMath::Log10(temp1[1]);
         if (Hoption.Logz) temp1[2] = TMath::Log10(temp1[2]);
         view->WCtoNDC(temp1, &temp2[0]);
         x[k] = temp2[0];
         y[k] = temp2[1];
      }
      x[3] = x[0];
      y[3] = y[0];
      if (tri1 || tri2) PaintLevels(p,x,y);
      if (!tri1 && !tri2 && !wire) {
         gPad->PaintFillArea(3,x,y);
         PaintLevels(p,x,y,nblev,glev);
      }
      if (!tri2) gPad->PaintPolyLine(4,x,y);
      if (markers) {
         if (markers0) {
            PaintPolyMarker0(3,x,y);
         } else {
            fGraph2D->SetMarkerStyle(fGraph2D->GetMarkerStyle());
            fGraph2D->SetMarkerSize(fGraph2D->GetMarkerSize());
            fGraph2D->SetMarkerColor(fGraph2D->GetMarkerColor());
            fGraph2D->TAttMarker::Modify();
            gPad->PaintPolyMarker(3,x,y);
         }
      }
endloop:
      continue;
   }
   fGraph2D->SetFillStyle(fs);
   fGraph2D->SetLineStyle(lst);
   fGraph2D->TAttLine::Modify();
   fGraph2D->TAttFill::Modify();

   if (glev) delete [] glev;
}
