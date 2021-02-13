// @(#)root/graf:$Id$
// Author: Rene Brun   16/10/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TArc.h"

#include <iostream>

ClassImp(TArc);

/** \class TArc
\ingroup BasicGraphics

Create an Arc.

An arc is specified with the position of its centre, its radius a minimum and
maximum angle. The attributes of the outline line are given via TAttLine. The
attributes of the fill area are given via TAttFill
*/

////////////////////////////////////////////////////////////////////////////////
/// Arc  default constructor.

TArc::TArc(): TEllipse()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Arc  normal constructor.
///
/// \param[in] x1,y1    coordinates of centre of arc
/// \param[in] r1       arc radius
/// \param[in] phimin   min angle in degrees (default is 0-->360)
/// \param[in] phimax   max angle in degrees (default is 0-->360)
///
/// When a circle sector only is drawn, the lines connecting the center
/// of the arc to the edges are drawn by default. One can specify
/// the drawing option "only" to not draw these lines.

TArc::TArc(Double_t x1, Double_t y1,Double_t r1,Double_t phimin,Double_t phimax)
      :TEllipse(x1,y1,r1,r1,phimin,phimax,0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TArc::TArc(const TArc &arc) : TEllipse(arc)
{
   ((TArc&)arc).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Arc default destructor.

TArc::~TArc()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this arc to arc.

void TArc::Copy(TObject &arc) const
{
   TEllipse::Copy(arc);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this arc with new coordinates.

void TArc::DrawArc(Double_t x1, Double_t y1,Double_t r1,Double_t phimin,Double_t phimax,Option_t *option)
{
   TArc *newarc = new TArc(x1, y1, r1, phimin, phimax);
   TAttLine::Copy(*newarc);
   TAttFill::Copy(*newarc);
   newarc->SetBit(kCanDelete);
   newarc->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TArc::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TArc::Class())) {
      out<<"   ";
   } else {
      out<<"   TArc *";
   }
   out<<"arc = new TArc("<<fX1<<","<<fY1<<","<<fR1
      <<","<<fPhimin<<","<<fPhimax<<");"<<std::endl;

   SaveFillAttributes(out,"arc",0,1001);
   SaveLineAttributes(out,"arc",1,1,1);

   if (GetNoEdges()) out<<"   arc->SetNoEdges();"<<std::endl;

   out<<"   arc->Draw();"<<std::endl;
}
