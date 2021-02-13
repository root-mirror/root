// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


// This is a minimal-change import of GLU libtess from:
//   git://anongit.freedesktop.org/git/mesa/glu
//
// Changes to make it build in the horrible one-file-hack way:
// - remove include gl.h from glu.h (and replace with the mini-gl defs below);
// - comment out clashing typedef EdgePair in tess.c;
// - use -Wno-unused-parameter for this cxx file.

// Sergey: first include gl code before any normal ROOT includes,
// try to avoid clash with other system includes through Rtypes.h

#include "glu/GL_glu.h"

#include <ROOT/REveGluTess.hxx>

#include <math.h>
#include <stdlib.h>
#include <stdexcept>

namespace ROOT {
namespace Experimental {
namespace EveGlu {

//////////////////////////////////////////////////////////////////////////////////////
/// TestTriangleHandler is just helper class to get access to protected members of TriangleCollector
/// Hide static declarations, let use "native" GL types
//////////////////////////////////////////////////////////////////////////////////////

class TestTriangleHandler {
public:
   static void tess_begin(GLenum type, TriangleCollector *tc);
   static void tess_vertex(Int_t *vi, TriangleCollector *tc);
   static void
   tess_combine(GLdouble coords[3], void *vertex_data[4], GLfloat weight[4], void **outData, TriangleCollector *tc);
   static void tess_end(TriangleCollector *tc);
};


//////////////////////////////////////////////////////////////////////////////////////
/// tess_begin

void TestTriangleHandler::tess_begin(GLenum type, TriangleCollector *tc)
{
   tc->fNVertices = 0;
   tc->fV0 = tc->fV1 = -1;
   tc->fType = type;
}

//////////////////////////////////////////////////////////////////////////////////////
/// tess_vertex

void TestTriangleHandler::tess_vertex(Int_t *vi, TriangleCollector *tc)
{
   tc->process_vertex(*vi);
}

//////////////////////////////////////////////////////////////////////////////////////
/// tess_combine

void TestTriangleHandler::tess_combine(GLdouble /*coords*/[3], void * /*vertex_data*/[4],
                                     GLfloat  /*weight*/[4], void ** /*outData*/,
                                     TriangleCollector * /*tc*/)
{
   throw std::runtime_error("GLU::TriangleCollector tesselator requested vertex combining -- not supported yet.");
}

//////////////////////////////////////////////////////////////////////////////////////
/// tess_end

void TestTriangleHandler::tess_end(TriangleCollector *tc)
{
   tc->fType = GL_NONE;
}


//==============================================================================
// TriangleCollector
//==============================================================================

//////////////////////////////////////////////////////////////////////////////////////
/// constructor

TriangleCollector::TriangleCollector()
{
   fType = GL_NONE;
   fTess = gluNewTess();
   if (!fTess) throw std::bad_alloc();

   gluTessCallback(fTess, (GLenum)GLU_TESS_BEGIN_DATA,   (_GLUfuncptr) TestTriangleHandler::tess_begin);
   gluTessCallback(fTess, (GLenum)GLU_TESS_VERTEX_DATA,  (_GLUfuncptr) TestTriangleHandler::tess_vertex);
   gluTessCallback(fTess, (GLenum)GLU_TESS_COMBINE_DATA, (_GLUfuncptr) TestTriangleHandler::tess_combine);
   gluTessCallback(fTess, (GLenum)GLU_TESS_END_DATA,     (_GLUfuncptr) TestTriangleHandler::tess_end);
}

//////////////////////////////////////////////////////////////////////////////////////
/// destructor

TriangleCollector::~TriangleCollector()
{
   gluDeleteTess(fTess);
}

//////////////////////////////////////////////////////////////////////////////////////
/// add triangle

void TriangleCollector::add_triangle(UInt_t v0, UInt_t v1, UInt_t v2)
{
   fPolyDesc.emplace_back(3);
   fPolyDesc.emplace_back(v0);
   fPolyDesc.emplace_back(v1);
   fPolyDesc.emplace_back(v2);
   ++fNTriangles;
}

//////////////////////////////////////////////////////////////////////////////////////
/// process_vertex

void ROOT::Experimental::EveGlu::TriangleCollector::process_vertex(UInt_t vi)
{
   ++fNVertices;

   if (fV0 == -1) {
      fV0 = vi;
      return;
   }
   if (fV1 == -1) {
      fV1 = vi;
      return;
   }

   switch (fType)
   {
      case GL_TRIANGLES:
      {
         add_triangle(fV0, fV1, vi);
         fV0 = fV1 = -1;
         break;
      }
      case GL_TRIANGLE_STRIP:
      {
         if (fNVertices % 2 == 0)
            add_triangle(fV1, fV0, vi);
         else
            add_triangle(fV0, fV1, vi);
         fV0 = fV1;
         fV1 = vi;
         break;
      }
      case GL_TRIANGLE_FAN:
      {
         add_triangle(fV0, fV1, vi);
         fV1 = vi;
         break;
      }
      default:
      {
         throw std::runtime_error("GLU::TriangleCollector unexpected type in tess_vertex callback.");
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////
/// ProcessData

void ROOT::Experimental::EveGlu::TriangleCollector::ProcessData(const std::vector<Double_t>& verts,
                                    const std::vector<UInt_t>  & polys,
                                    const Int_t                  n_polys)
{
   const Double_t *pnts = &verts[0];
   const UInt_t   *pols = &polys[0];

   for (Int_t i = 0, j = 0; i < n_polys; ++i)
   {
      Int_t n_points = pols[j++];

      gluTessBeginPolygon(fTess, this);
      gluTessBeginContour(fTess);

      for (Int_t k = 0; k < n_points; ++k, ++j)
      {
         gluTessVertex(fTess, (Double_t*) pnts + pols[j] * 3, (GLvoid*) &pols[j]);
      }

      gluTessEndContour(fTess);

      static int except_cnt = 0;

      try {
         gluTessEndPolygon(fTess);
      } catch(...) {
         if (except_cnt++ < 100) printf("Catch exception gluTessEndPolygon!\n");
      }
   }
}

} // namespace EveGlu
} // namespace Experimental
} // namespace ROOT
