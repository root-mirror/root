// @(#)root/geom:$Id$
// Author: Andrei Gheata   20/01/2020

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoTypedefs
#define ROOT_TGeoTypedefs

/// Typedefs used by the geometry group
#include <vector>

namespace Tessellated {

  using Vertex_t    = ROOT::Geom::Vertex_t;
  using VertexVec_t = std::vector<Vertex_t>;

} // namespace Tessellated

#endif
