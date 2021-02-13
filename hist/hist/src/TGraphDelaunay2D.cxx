// @(#)root/hist:$Id: TGraphDelaunay2D.cxx,v 1.00
// Author: Olivier Couet, Luke Jones (Royal Holloway, University of London)

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGraph2D.h"
#include "TGraphDelaunay2D.h"

ClassImp(TGraphDelaunay2D);


/** \class TGraphDelaunay2D
    \ingroup Hist
TGraphDelaunay2D generates a Delaunay triangulation of a TGraph2D. This
triangulation code derives from an implementation done by Luke Jones
(Royal Holloway, University of London) in April 2002 in the PAW context.

This software cannot be guaranteed to work under all circumstances. They
were originally written to work with a few hundred points in an XY space
with similar X and Y ranges.

Definition of Delaunay triangulation (After B. Delaunay):
For a set S of points in the Euclidean plane, the unique triangulation DT(S)
of S such that no point in S is inside the circumcircle of any triangle in
DT(S). DT(S) is the dual of the Voronoi diagram of S. If n is the number of
points in S, the Voronoi diagram of S is the partitioning of the plane
containing S points into n convex polygons such that each polygon contains
exactly one point and every point in a given polygon is closer to its
central point than to any other. A Voronoi diagram is sometimes also known
as a Dirichlet tessellation.

\image html tgraph2d_delaunay.png

[This applet](http://www.cs.cornell.edu/Info/People/chew/Delaunay.html)
gives a nice practical view of Delaunay triangulation and Voronoi diagram.
*/

/// TGraphDelaunay2D normal constructor
TGraphDelaunay2D::TGraphDelaunay2D(TGraph2D *g ) :
   TNamed("TGraphDelaunay2D","TGraphDelaunay2D"),
   fGraph2D(g),
   fDelaunay((g) ? g->GetN() : 0, (g) ? g->GetX() : nullptr, (g) ? g->GetY() : nullptr, (g) ? g->GetZ() : nullptr ,
             (g) ? g->GetXmin() : 0, (g) ? g->GetXmax() : 0,
             (g) ? g->GetYmin() : 0, (g) ? g->GetYmax() : 0 )

{}

