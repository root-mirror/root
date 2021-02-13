// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  03/08/2004
// NOTE: This code moved from obsoleted TGLSceneObject.h / .cxx - see these
// attic files for previous CVS history

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLFaceSet
#define ROOT_TGLFaceSet

#include "TGLLogicalShape.h"
#include "CsgOps.h"
#include <vector>

///////////////////////////////////////////////////////////////////////
class TGLFaceSet : public TGLLogicalShape
{
private:
   std::vector<Double_t> fVertices;
   std::vector<Double_t> fNormals;
   std::vector<Int_t>    fPolyDesc;
   UInt_t                fNbPols;

   static Bool_t fgEnforceTriangles;

public:
   TGLFaceSet(const TBuffer3D & buffer);

   virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

   void SetFromMesh(const RootCsg::TBaseMesh *m);
   void CalculateNormals();
   void EnforceTriangles();

   std::vector<Double_t>& GetVertices() { return fVertices; }
   std::vector<Double_t>& GetNormals()  { return fNormals;  }
   std::vector<Int_t>&    GetPolyDesc() { return fPolyDesc; }
   UInt_t                 GetNbPols()   { return fNbPols;   }

   static Bool_t GetEnforceTriangles();
   static void   SetEnforceTriangles(Bool_t e);

private:
   void  GLDrawPolys() const;
   Int_t CheckPoints(const Int_t *source, Int_t *dest) const;

   static Bool_t Eq(const Double_t *p1, const Double_t *p2);

   ClassDef(TGLFaceSet,0) // a faceset logical shape
};

#endif
