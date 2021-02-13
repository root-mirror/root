// @(#)root/gl:$Id$
// Author:  Richard Maunder  25/05/2005

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLBoundingBox
#define ROOT_TGLBoundingBox

#include "TGLUtil.h"

#include <vector>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLBoundingBox                                                       //
//                                                                      //
// Concrete class describing an orientated (free) or axis aligned box   //
// of 8 verticies. Supports methods for setting aligned or orientated   //
// boxes, find volume, axes, extents, centers, face planes etc.         //
// Also tests for overlap testing of planes and other bounding boxes,   //
// with fast sphere approximation.                                      //
//////////////////////////////////////////////////////////////////////////

// TODO: Create more compact version + axis aligned version, both with lazy
// sphere testing.
class TGLBoundingBox
{
private:
   // Fields

   // Box vertices are indexed thus (OpenGL is left handed by default)
   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
   //
   // 0123 'far' face
   // 4567 'near' face
   //
   // This could be more compact:
   // For orientated box 3 vertices which form plane cutting box
   // diagonally (e.g. 0,5,6 or 1,3,6 etc) would fix in space.
   // For axis aligned 2 verticies would suffice.
   // Rest could be calculated on demand - however speed more important
   // than memory considerations
   TGLVertex3              fVertex[8];  //! the 8 bounding box vertices
   Double_t                fVolume;     //! box volume - cached for speed
   Double_t                fDiagonal;   //! max box diagonal - cached for speed
   TGLVector3              fAxes[3];    //! box axes in global frame - cached for speed
   TGLVector3              fAxesNorm[3];//! normalised box axes in global frame - cached for speed

   // Methods
   void     UpdateCache();
   Bool_t   ValidIndex(UInt_t index) const { return (index < 8); }
   Double_t Min(UInt_t index) const;
   Double_t Max(UInt_t index) const;

public:
   TGLBoundingBox();
   TGLBoundingBox(const TGLVertex3 vertex[8]);
   TGLBoundingBox(const Double_t vertex[8][3]);
   TGLBoundingBox(const TGLVertex3 & lowVertex, const TGLVertex3 & highVertex);
   TGLBoundingBox(const TGLBoundingBox & other);
   virtual ~TGLBoundingBox(); // ClassDef introduces virtual fns

   // Set orientated box
   TGLBoundingBox & operator =(const TGLBoundingBox & other);
   void Set(const TGLVertex3 vertex[8]);
   void Set(const Double_t vertex[8][3]);
   void Set(const TGLBoundingBox & other);
   void SetEmpty();

   // Set axis aligned box
   void SetAligned(const TGLVertex3 & lowVertex, const TGLVertex3 & highVertex); // axis aligned
   void SetAligned(UInt_t nbPnts, const Double_t * pnts); // axis aligned
   void MergeAligned(const TGLBoundingBox & other);
   void ExpandAligned(const TGLVertex3 & point);

   // Manipulation
   void Transform(const TGLMatrix & matrix);
   void Scale(Double_t factor);
   void Scale(Double_t xFactor, Double_t yFactor, Double_t zFactor);
   void Translate(const TGLVector3 & offset);

   // Single vertex accessors
   const TGLVertex3 & operator [] (UInt_t index) const;
   const TGLVertex3 & Vertex(UInt_t index) const;
   Double_t XMin() const { return Min(0); }
   Double_t XMax() const { return Max(0); }
   Double_t YMin() const { return Min(1); }
   Double_t YMax() const { return Max(1); }
   Double_t ZMin() const { return Min(2); }
   Double_t ZMax() const { return Max(2); }
   TGLVertex3 MinAAVertex() const;
   TGLVertex3 MaxAAVertex() const;

   // Multiple vertices accessors
   const TGLVertex3* Vertices() const;           // All 8 box vertices
   Int_t             NumVertices() const { return 8; }

   enum EFace { kFaceLowX, kFaceHighX, kFaceLowY, kFaceHighY, kFaceLowZ, kFaceHighZ, kFaceCount };
   const std::vector<UInt_t> & FaceVertices(EFace face) const; // 4 box face vertices

   // Other properties
   TGLVertex3   Center() const;
   TGLVector3   Extents() const;
   const  TGLVector3 & Axis(UInt_t i, Bool_t normalised = kTRUE) const;
   Bool_t       IsEmpty()  const;
   Double_t     Volume()   const { return fVolume;   }
   Double_t     Diagonal() const { return fDiagonal; }
   void         PlaneSet(TGLPlaneSet_t & planeSet) const;
   TGLPlane     GetNearPlane() const;

   // Overlap testing
   Rgl::EOverlap Overlap(const TGLPlane & plane) const;
   Rgl::EOverlap Overlap(const TGLBoundingBox & box) const;

   void Draw(Bool_t solid = kFALSE) const;
   void Dump() const;

   ClassDef(TGLBoundingBox,0); // a 3D orientated bounding box
};

//______________________________________________________________________________
inline TGLBoundingBox & TGLBoundingBox::operator =(const TGLBoundingBox & other)
{
   // Check for self-assignment
   if (this != &other) {
      Set(other);
   }
   return *this;
}

//______________________________________________________________________________
inline const TGLVertex3 & TGLBoundingBox::operator [] (UInt_t index) const
{
   return fVertex[index];
}

//______________________________________________________________________________
inline const TGLVertex3 & TGLBoundingBox::Vertex(UInt_t index) const
{
   return fVertex[index];
}

//______________________________________________________________________________
inline const TGLVertex3* TGLBoundingBox::Vertices() const
{
   return fVertex;
}

//______________________________________________________________________________
inline TGLVector3 TGLBoundingBox::Extents() const
{
   // Return the local axis entents of the box
   return TGLVector3(Axis(0,kFALSE).Mag(),
                     Axis(1,kFALSE).Mag(),
                     Axis(2,kFALSE).Mag());
}

//______________________________________________________________________________
inline TGLVertex3 TGLBoundingBox::Center() const
{
   // Return the center vertex of the box
   return TGLVertex3((fVertex[0].X() + fVertex[6].X())/2.0,
                     (fVertex[0].Y() + fVertex[6].Y())/2.0,
                     (fVertex[0].Z() + fVertex[6].Z())/2.0);
}

//______________________________________________________________________________
inline const TGLVector3 & TGLBoundingBox::Axis(UInt_t i, Bool_t normalised) const
{
   // Return a vector representing axis of index i (0:X, 1:Y, 2:Z).
   // Vector can be as-is (edge, magnitude == extent) or normalised (default)
   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
   //

   if (normalised) {
      return fAxesNorm[i];
   } else {
      return fAxes[i];
   }
}

//______________________________________________________________________________
inline Bool_t TGLBoundingBox::IsEmpty() const
{
   // Return kTRUE if box has zero diagonal - kFALSE otherwise

   // TODO: Round errors - should have epsilon test
   return (Diagonal() == 0.0);
}

#endif // ROOT_TGLBoundingBox
