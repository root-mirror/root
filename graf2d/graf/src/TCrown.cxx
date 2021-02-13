// @(#)root/graf:$Id$
// Author: Rene Brun   108/08/2002

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "TCrown.h"
#include "TVirtualPad.h"

ClassImp(TCrown);

/** \class TCrown
\ingroup BasicGraphics

To draw a Crown.

A crown is specified with the position of its centre, its inner/outer radius
a minimum and maximum angle. The attributes of the outline line are given via
TAttLine. The attributes of the fill area are given via TAttFill.

Example:

Begin_Macro(source)
{
   TCanvas *c1 = new TCanvas("c1","c1",400,400);
   TCrown cr1(.5,.5,.3,.4);
   cr1.SetLineStyle(2);
   cr1.SetLineWidth(4);
   cr1.Draw();
   TCrown cr2(.5,.5,.2,.3,45,315);
   cr2.SetFillColor(38);
   cr2.SetFillStyle(3010);
   cr2.Draw();
   TCrown cr3(.5,.5,.2,.3,-45,45);
   cr3.SetFillColor(50);
   cr3.SetFillStyle(3025);
   cr3.Draw();
   TCrown cr4(.5,.5,.0,.2);
   cr4.SetFillColor(4);
   cr4.SetFillStyle(3008);
   cr4.Draw();
   return c1;
}
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// Crown default constructor.

TCrown::TCrown(): TEllipse()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Crown normal constructor.
///
/// \param[in] x1,y1    coordinates of centre of crown
/// \param[in] radin    inner crown radius
/// \param[in] radout   outer crown radius
/// \param[in] phimin   min angle in degrees (default is 0)
/// \param[in] phimax   max angle in degrees (default is 360)
///
/// When a crown sector only is drawn, the lines connecting the center
/// of the crown to the edges are drawn by default. One can specify
/// the drawing option "only" to not draw these lines.

TCrown::TCrown(Double_t x1, Double_t y1,Double_t radin, Double_t radout,Double_t phimin,Double_t phimax)
      :TEllipse(x1,y1,radin,radout,phimin,phimax,0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Crown copy constructor.

TCrown::TCrown(const TCrown &crown) : TEllipse(crown)
{

   ((TCrown&)crown).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Crown default destructor.

TCrown::~TCrown()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this crown to crown.

void TCrown::Copy(TObject &crown) const
{

   TEllipse::Copy(crown);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a crown.
///
/// If crown is filled, return OK if we are inside
/// otherwise, crown is found if near the crown edges.

Int_t TCrown::DistancetoPrimitive(Int_t px, Int_t py)
{

   const Double_t kPI = TMath::Pi();
   Double_t x = gPad->PadtoX(gPad->AbsPixeltoX(px)) - fX1;
   Double_t y = gPad->PadtoY(gPad->AbsPixeltoY(py)) - fY1;

   Double_t r1 = fR1;
   Double_t r2 = fR2;
   Double_t r  = TMath::Sqrt(x*x+y*y);

   if (r1>r2) {
      r1 = fR2;
      r2 = fR1;
   }

   Int_t dist = 9999;
   if (r > r2) return dist;
   if (r < r1) return dist;
   if (fPhimax-fPhimin < 360) {
      Double_t phi = 180*TMath::ACos(x/r)/kPI;
      if (y<0) phi = 360-phi;
      Double_t phi1 = fPhimin;
      Double_t phi2 = fPhimax;
      if (phi1<0) phi1=phi1+360;
      if (phi2<0) phi2=phi2+360;
      if (phi2<phi1) {
         if (phi < phi1 && phi > phi2) return dist;
      } else {
         if (phi < phi1) return dist;
         if (phi > phi2) return dist;
      }
   }

   if (GetFillColor() && GetFillStyle()) {
      return 0;
   } else {
      if (TMath::Abs(r2-r)/r2 < 0.02) return 0;
      if (TMath::Abs(r1-r)/r1 < 0.02) return 0;
   }
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this crown with new coordinates.

void TCrown::DrawCrown(Double_t x1, Double_t y1,Double_t radin,Double_t radout,Double_t phimin,Double_t phimax,Option_t *option)
{

   TCrown *newcrown = new TCrown(x1, y1, radin, radout, phimin, phimax);
   TAttLine::Copy(*newcrown);
   TAttFill::Copy(*newcrown);
   newcrown->SetBit(kCanDelete);
   newcrown->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event
///
/// For the time being TEllipse::ExecuteEvent is used.

void TCrown::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{

   TEllipse::ExecuteEvent(event,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this crown with its current attributes.

void TCrown::Paint(Option_t *)
{

   const Double_t kPI = TMath::Pi();
   const Int_t np = 40;
   static Double_t x[2*np+3], y[2*np+3];
   TAttLine::Modify();
   TAttFill::Modify();

   Double_t angle,dx,dy;
   Double_t dphi = (fPhimax-fPhimin)*kPI/(180*np);
   Double_t ct   = TMath::Cos(kPI*fTheta/180);
   Double_t st   = TMath::Sin(kPI*fTheta/180);
   Int_t i;
   //compute outer points
   for (i=0;i<=np;i++) {
      angle = fPhimin*kPI/180 + Double_t(i)*dphi;
      dx    = fR2*TMath::Cos(angle);
      dy    = fR2*TMath::Sin(angle);
      x[i]  = fX1 + dx*ct - dy*st;
      y[i]  = fY1 + dx*st + dy*ct;
   }
   //compute inner points
   for (i=0;i<=np;i++) {
      angle = fPhimin*kPI/180 + Double_t(i)*dphi;
      dx    = fR1*TMath::Cos(angle);
      dy    = fR1*TMath::Sin(angle);
      x[2*np-i+1]  = fX1 + dx*ct - dy*st;
      y[2*np-i+1]  = fY1 + dx*st + dy*ct;
   }
   x[2*np+2]  = x[0];
   y[2*np+2]  = y[0];
   if (fPhimax-fPhimin >= 360 ) {
      // a complete filled crown
      if (GetFillColor()  && GetFillStyle()) {
         gPad->PaintFillArea(2*np+2,x,y);
      }
      // a complete empty crown
      if (GetLineStyle()) {
         gPad->PaintPolyLine(np+1,x,y);
         gPad->PaintPolyLine(np+1,&x[np+1],&y[np+1]);
      }
   } else {
      //crown segment
      if (GetFillColor()  && GetFillStyle()) gPad->PaintFillArea(2*np+2,x,y);
      if (GetLineStyle()) gPad->PaintPolyLine(2*np+3,x,y);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TCrown::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TCrown::Class())) {
      out<<"   ";
   } else {
      out<<"   TCrown *";
   }
   out<<"crown = new TCrown("<<fX1<<","<<fY1<<","<<fR1<<","<<fR2
      <<","<<fPhimin<<","<<fPhimax<<");"<<std::endl;

   SaveFillAttributes(out,"crown",0,1001);
   SaveLineAttributes(out,"crown",1,1,1);

   if (GetNoEdges()) out<<"   crown->SetNoEdges();"<<std::endl;

   out<<"   crown->Draw();"<<std::endl;
}
