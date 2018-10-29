// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT7_REveGluTess
#define ROOT7_REveGluTess

#include "Rtypes.h"

#include <vector>

struct GLUtesselator;

namespace ROOT {
namespace Experimental {
namespace EveGlu {

//==============================================================================
// TriangleCollector
//==============================================================================

class TestTriangleHandler;

class TriangleCollector {
   friend class TestTriangleHandler;
protected:
   GLUtesselator *fTess{nullptr};
   Int_t fNTriangles{0};
   Int_t fNVertices{0};
   Int_t fV0{-1}, fV1{-1};
   Int_t fType{0};
   std::vector<Int_t> fPolyDesc;

   void add_triangle(Int_t v0, Int_t v1, Int_t v2);
   void process_vertex(Int_t vi);

public:
   TriangleCollector();
   ~TriangleCollector();

   // Process polygons
   void ProcessData(const std::vector<Double_t> &verts, const std::vector<Int_t> &polys, const Int_t n_polys);

   // Get output
   Int_t GetNTrianlges() { return fNTriangles; }
   std::vector<Int_t> &RefPolyDesc() { return fPolyDesc; }
};

} // namespace EveGlu
} // namespace Experimental
} // namespace ROOT

#endif
