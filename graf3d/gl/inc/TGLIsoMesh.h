// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  06/01/2009

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLIsoMesh
#define ROOT_TGLIsoMesh

#include <vector>

#include "RtypesCore.h"
#include "TAxis.h"

class TGLBoxCut;

namespace Rgl {
namespace Mc {

/*
TIsoMesh - set of vertices, per-vertex normals, "triangles".
Each "triangle" is a triplet of indices, pointing into vertices
and normals arrays. For example, triangle t = {1, 4, 6}
has vertices &fVerts[1 * 3], &fVerts[4 * 3], &fVerts[6 * 3];
and normals &fNorms[1 * 3], &fNorms[4 * 3], &fNorms[6 * 3]
"V" parameter should be Float_t or Double_t (or some
integral type?).

Prefix "T" in a class name only for code-style checker.
*/

template<class V>
class TIsoMesh {
public:
   UInt_t AddVertex(const V *v)
   {
      const UInt_t index = UInt_t(fVerts.size() / 3);
      fVerts.push_back(v[0]);
      fVerts.push_back(v[1]);
      fVerts.push_back(v[2]);

      return index;
   }

   void AddNormal(const V *n)
   {
      fNorms.push_back(n[0]);
      fNorms.push_back(n[1]);
      fNorms.push_back(n[2]);
   }

   UInt_t AddTriangle(const UInt_t *t)
   {
      const UInt_t index = UInt_t(fTris.size() / 3);
      fTris.push_back(t[0]);
      fTris.push_back(t[1]);
      fTris.push_back(t[2]);

      return index;
   }

   void Swap(TIsoMesh &rhs)
   {
      std::swap(fVerts, rhs.fVerts);
      std::swap(fNorms, rhs.fNorms);
      std::swap(fTris, rhs.fTris);
   }

   void ClearMesh()
   {
      fVerts.clear();
      fNorms.clear();
      fTris.clear();
   }

   std::vector<V>      fVerts;
   std::vector<V>      fNorms;
   std::vector<UInt_t> fTris;
};

/*
TGridGeometry describes ranges and cell steps (scales are
already in steps and ranges).
*/
template<class V>
class TGridGeometry {
public:
   enum EVertexPosition{
      kBinCenter,
      kBinEdge
   };

   TGridGeometry() : fMinX(0),  fStepX(0),
                     fMinY(0),  fStepY(0),
                     fMinZ(0),  fStepZ(0),
                     fXScaleInverted(1.),
                     fYScaleInverted(1.),
                     fZScaleInverted(1.)
   {
      //Default constructor.
   }

   TGridGeometry(const TAxis *x, const TAxis *y, const TAxis *z,
                 Double_t xs = 1., Double_t ys = 1., Double_t zs = 1.,
                 EVertexPosition pos = kBinCenter)
         : fMinX(0),  fStepX(0),
           fMinY(0),  fStepY(0),
           fMinZ(0),  fStepZ(0),
           fXScaleInverted(1.),
           fYScaleInverted(1.),
           fZScaleInverted(1.)
   {
      //Define geometry using TAxis.
      if (pos == kBinCenter) {
         fMinX  = V(x->GetBinCenter(x->GetFirst()));
         fStepX = V((x->GetBinCenter(x->GetLast()) - fMinX) / (x->GetNbins() - 1));
         fMinY  = V(y->GetBinCenter(y->GetFirst()));
         fStepY = V((y->GetBinCenter(y->GetLast()) - fMinY) / (y->GetNbins() - 1));
         fMinZ  = V(z->GetBinCenter(z->GetFirst()));
         fStepZ = V((z->GetBinCenter(z->GetLast()) - fMinZ) / (z->GetNbins() - 1));

         fMinX *= xs, fStepX *= xs;
         fMinY *= ys, fStepY *= ys;
         fMinZ *= zs, fStepZ *= zs;
      } else if (pos == kBinEdge) {
         fMinX  = V(x->GetBinLowEdge(x->GetFirst()));
         fStepX = V((x->GetBinUpEdge(x->GetLast()) - fMinX) / (x->GetNbins()));
         fMinY  = V(y->GetBinLowEdge(y->GetFirst()));
         fStepY = V((y->GetBinUpEdge(y->GetLast()) - fMinY) / (y->GetNbins()));
         fMinZ  = V(z->GetBinLowEdge(z->GetFirst()));
         fStepZ = V((z->GetBinUpEdge(z->GetLast()) - fMinZ) / (z->GetNbins()));

         fMinX *= xs, fStepX *= xs;
         fMinY *= ys, fStepY *= ys;
         fMinZ *= zs, fStepZ *= zs;
      }

      fXScaleInverted = 1. / xs;
      fYScaleInverted = 1. / ys;
      fZScaleInverted = 1. / zs;
   }

   V fMinX;
   V fStepX;

   V fMinY;
   V fStepY;

   V fMinZ;
   V fStepZ;

   V fXScaleInverted;
   V fYScaleInverted;
   V fZScaleInverted;
};

}//namespace Mc

//Auxilary functions to draw an iso mesh in different modes.
void DrawMesh(const std::vector<Float_t> &vs, const std::vector<Float_t> &ns,
              const std::vector<UInt_t> &ts);
void DrawMesh(const std::vector<Double_t> &vs, const std::vector<Double_t> &ns,
              const std::vector<UInt_t> &ts);

void DrawMesh(const std::vector<Float_t> &vs, const std::vector<UInt_t> &fTS);
void DrawMesh(const std::vector<Double_t> &vs, const std::vector<UInt_t> &fTS);

void DrawMesh(const std::vector<Float_t> &vs, const std::vector<Float_t> &ns,
              const std::vector<UInt_t> &ts, const TGLBoxCut &box);
void DrawMesh(const std::vector<Double_t> &vs, const std::vector<Double_t> &ns,
              const std::vector<UInt_t> &ts, const TGLBoxCut &box);

void DrawMesh(const std::vector<Float_t> &vs, const std::vector<UInt_t> &ts,
              const TGLBoxCut &box);
void DrawMesh(const std::vector<Double_t> &vs, const std::vector<UInt_t> &ts,
              const TGLBoxCut &box);

void DrawMesh(const std::vector<Double_t> &vs, const std::vector<UInt_t> &ts,
              const TGLBoxCut &box);
void DrawMesh(const std::vector<Float_t> &vs, const std::vector<UInt_t> &ts,
              const TGLBoxCut &box);

void DrawMapleMesh(const std::vector<Double_t> &vs, const std::vector<Double_t> &ns,
                   const std::vector<UInt_t> &ts);
void DrawMapleMesh(const std::vector<Double_t> &vs, const std::vector<Double_t> &ns,
                   const std::vector<UInt_t> &ts, const TGLBoxCut & box);

}//namespace Rgl

#endif
