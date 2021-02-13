// @(#)root/geom:$Id$
// Author: Andrei Gheata   20/12/19

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoTessellated
#define ROOT_TGeoTessellated

#include "TGeoVector3.h"
#include "TGeoTypedefs.h"
#include "TGeoBBox.h"

class TGeoFacet {
   using Vertex_t = Tessellated::Vertex_t;
   using VertexVec_t = Tessellated::VertexVec_t;

private:
   int fIvert[4] = {0, 0, 0, 0};     // Vertex indices in the array
   VertexVec_t *fVertices = nullptr; //! array of vertices
   int fNvert = 0;                   // number of vertices (can be 3 or 4)
   bool fShared = false;             // Vector of vertices shared flag

public:
   TGeoFacet() {}
   TGeoFacet(const TGeoFacet &other);

   ~TGeoFacet()
   {
      if (!fShared)
         delete fVertices;
   }

   const TGeoFacet &operator=(const TGeoFacet &other);

   // Triangular facet
   TGeoFacet(const Vertex_t &pt0, const Vertex_t &pt1, const Vertex_t &pt2) : fIvert{0, 1, 2}
   {
      fVertices = new VertexVec_t;
      fVertices->push_back(pt0);
      fVertices->push_back(pt1);
      fVertices->push_back(pt2);
      fNvert = 3;
   }

   // Quadrilateral facet
   TGeoFacet(const Vertex_t &pt0, const Vertex_t &pt1, const Vertex_t &pt2, const Vertex_t &pt3) : fIvert{0, 1, 2, 3}
   {
      fVertices = new VertexVec_t;
      fVertices->push_back(pt0);
      fVertices->push_back(pt1);
      fVertices->push_back(pt2);
      fVertices->push_back(pt3);
      fNvert = 4;
   }

   TGeoFacet(VertexVec_t *vertices, int nvert, int i0 = -1, int i1 = -1, int i2 = -1, int i3 = -1)
   {
      fShared = true;
      SetVertices(vertices, nvert, i0, i1, i2, i3);
   }

   static int CompactFacet(Vertex_t *vert, int nvertices);

   void SetVertices(VertexVec_t *vertices, int nvert = 0, int i0 = -1, int i1 = -1, int i2 = -1, int i3 = -1)
   {
      if (!fShared)
         delete fVertices;
      fVertices = vertices;
      if (nvert > 0) {
         fIvert[0] = i0;
         fIvert[1] = i1;
         fIvert[2] = i2;
         fIvert[3] = i3;
      }
      fNvert = nvert;
      fShared = true;
   }

   Vertex_t ComputeNormal(bool &degenerated) const;
   int GetNvert() const { return fNvert; }

   Vertex_t &GetVertex(int ivert) { return fVertices->operator[](fIvert[ivert]); }
   const Vertex_t &GetVertex(int ivert) const { return fVertices->operator[](fIvert[ivert]); }

   int GetVertexIndex(int ivert) const { return fIvert[ivert]; }

   bool Check() const;
   void Flip()
   {
      int iv = fIvert[0];
      fIvert[0] = fIvert[2];
      fIvert[2] = iv;
   }
   bool IsNeighbour(const TGeoFacet &other, bool &flip) const;
};

std::ostream &operator<<(std::ostream &os, TGeoFacet const &facet);

class TGeoTessellated : public TGeoBBox {

public:
   using Vertex_t = Tessellated::Vertex_t;

private:
   int fNfacets = 0;                // Number of facets
   int fNvert = 0;                  // Number of vertices
   int fNseg = 0;                   // Number of segments
   bool fDefined = false;           //! Shape fully defined
   bool fClosedBody = false;        // The faces are making a closed body
   std::vector<Vertex_t> fVertices; // List of vertices
   std::vector<TGeoFacet> fFacets;  // List of facets

public:
   // constructors
   TGeoTessellated() {}
   TGeoTessellated(const char *name, int nfacets = 0);
   TGeoTessellated(const char *name, const std::vector<Vertex_t> &vertices);
   // destructor
   virtual ~TGeoTessellated() {}

   void ComputeBBox();
   void CloseShape(bool check = true, bool fixFlipped = true, bool verbose = true);

   bool AddFacet(const Vertex_t &pt0, const Vertex_t &pt1, const Vertex_t &pt2);
   bool AddFacet(const Vertex_t &pt0, const Vertex_t &pt1, const Vertex_t &pt2, const Vertex_t &pt3);
   bool AddFacet(int i1, int i2, int i3);
   bool AddFacet(int i1, int i2, int i3, int i4);

   int GetNfacets() const { return fFacets.size(); }
   int GetNsegments() const { return fNseg; }
   int GetNvertices() const { return fNvert; }
   bool IsClosedBody() const { return fClosedBody; }
   bool IsDefined() const { return fDefined; }

   const TGeoFacet &GetFacet(int i) { return fFacets[i]; }
   const Vertex_t &GetVertex(int i) { return fVertices[i]; }

   virtual void AfterStreamer();
   virtual int DistancetoPrimitive(int, int) { return 99999; }
   virtual const TBuffer3D &GetBuffer3D(int reqSections, Bool_t localFrame) const;
   virtual void GetMeshNumbers(int &nvert, int &nsegs, int &npols) const;
   virtual int GetNmeshVertices() const { return fNvert; }
   virtual void InspectShape() const {}
   virtual TBuffer3D *MakeBuffer3D() const;
   virtual void Print(Option_t *option = "") const;
   virtual void SavePrimitive(std::ostream &, Option_t *) {}
   virtual void SetPoints(double *points) const;
   virtual void SetPoints(float *points) const;
   virtual void SetSegsAndPols(TBuffer3D &buff) const;
   virtual void Sizeof3D() const {}

   /// Resize and center the shape in a box of size maxsize
   void ResizeCenter(double maxsize);

   /// Flip all facets
   void FlipFacets()
   {
      for (auto facet : fFacets)
         facet.Flip();
   }

   bool CheckClosure(bool fixFlipped = true, bool verbose = true);

   /// Reader from .obj format
   static TGeoTessellated *ImportFromObjFormat(const char *objfile, bool check = false, bool verbose = false);

   ClassDef(TGeoTessellated, 1) // tessellated shape class
};

#endif
