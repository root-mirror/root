// @(#)root/g3d:$Id$
// Author: Rene Brun   22/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TBuffer.h"
#include "THashList.h"
#include "TObjArray.h"
#include "TGeometry.h"
#include "TNode.h"
#include "TMaterial.h"
#include "TBrowser.h"

TGeometry *gGeometry = nullptr;

ClassImp(TGeometry);

/** \class TGeometry
\ingroup g3d
TGeometry description.

The Geometry class describes the geometry of a detector.
The current implementation supports the GEANT3 style description.
A special program provided in the ROOT utilities (toroot) can be used
to automatically translate a GEANT detector geometry into a ROOT geometry.

a Geometry object is entered into the list of geometries into the
ROOT main object (see TROOT description) when the TGeometry
constructor is invoked.
Several geometries may coexist in memory.
/
A Geometry object consist of the following linked lists:

  - the TMaterial list (material definition only).
  - the TRotmatrix list (Rotation matrices definition only).
  - the TShape list (volume definition only).
  - the TNode list assembling all detector elements.

Only the Build and Draw functions for a geometry are currently supported.

The conversion program from Geant to Root has been added in the list
of utilities in utils directory.(see g2root)
The executable module of g2root can be found in $ROOTSYS/bin/g2root.

To use this conversion program, type the shell command:

~~~ {.cpp}
      g2root  geant_rzfile macro_name
~~~

for example

~~~ {.cpp}
      g2root na49.geom na49.C
~~~

will convert the GEANT RZ file na49.geom into a ROOT macro na49.C

To generate the Geometry structure within Root, do:

~~~ {.cpp}
  Root > .x na49.C
  Root > na49.Draw()
  Root > wh.x3d()    (this invokes the 3-d Root viewer)
  Root > TFile gna49("na49.root","NEW")  //open a new root file
  Root > na49.Write()                    //Write the na49 geometry structure
  Root > gna49.Write()                   //Write all keys (in this case only one)
~~~

Note: all keys are also written on closing of the file, gna49.Close or
when the program exits, Root closes all open files correctly.
Once this file has been written, in a subsequent session, simply do:

~~~ {.cpp}
  Root > TFile gna49("na49.root")
  Root > na49.Draw()
~~~

The figure below shows the geometry above using the x3d viewer.
This x3d viewer is invoked by selecting "View x3d" in the View menu
of a canvas (See example of this tool bar in TCanvas).

\image html g3d_na49.png
*/

////////////////////////////////////////////////////////////////////////////////
/// Geometry default constructor.

TGeometry::TGeometry()
{
   fMaterials       = new THashList(100,3);
   fMatrices        = new THashList(100,3);
   fShapes          = new THashList(500,3);
   fNodes           = new TList;
   fCurrentNode     = 0;
   fMaterialPointer = 0;
   fMatrixPointer   = 0;
   fShapePointer    = 0;
   gGeometry = this;
   fBomb            = 1;
   fMatrix          = 0;
   fX=fY=fZ         =0.0;
   fGeomLevel       =0;
   fIsReflection[fGeomLevel] = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry normal constructor.

TGeometry::TGeometry(const char *name,const char *title ) : TNamed (name, title)
{
   fMaterials       = new THashList(1000,3);
   fMatrices        = new THashList(1000,3);
   fShapes          = new THashList(5000,3);
   fNodes           = new TList;
   fCurrentNode     = 0;
   fMaterialPointer = 0;
   fMatrixPointer   = 0;
   fShapePointer    = 0;
   gGeometry = this;
   fBomb            = 1;
   fMatrix          = 0;
   fX=fY=fZ         =0.0;
   gROOT->GetListOfGeometries()->Add(this);
   fGeomLevel       =0;
   fIsReflection[fGeomLevel] = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TGeometry::TGeometry(const TGeometry& geo) :
  TNamed(geo),
  fMaterials(geo.fMaterials),
  fMatrices(geo.fMatrices),
  fShapes(geo.fShapes),
  fNodes(geo.fNodes),
  fMatrix(geo.fMatrix),
  fCurrentNode(geo.fCurrentNode),
  fMaterialPointer(geo.fMaterialPointer),
  fMatrixPointer(geo.fMatrixPointer),
  fShapePointer(geo.fShapePointer),
  fBomb(geo.fBomb),
  fGeomLevel(geo.fGeomLevel),
  fX(geo.fX),
  fY(geo.fY),
  fZ(geo.fZ)
{
   for(Int_t i=0; i<kMAXLEVELS; i++) {
      for(Int_t j=0; j<kVectorSize; j++)
         fTranslation[i][j]=geo.fTranslation[i][j];
      for(Int_t j=0; j<kMatrixSize; j++)
         fRotMatrix[i][j]=geo.fRotMatrix[i][j];
      fIsReflection[i]=geo.fIsReflection[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TGeometry& TGeometry::operator=(const TGeometry& geo)
{
   if(this!=&geo) {
      TNamed::operator=(geo);
      fMaterials=geo.fMaterials;
      fMatrices=geo.fMatrices;
      fShapes=geo.fShapes;
      fNodes=geo.fNodes;
      fMatrix=geo.fMatrix;
      fCurrentNode=geo.fCurrentNode;
      fMaterialPointer=geo.fMaterialPointer;
      fMatrixPointer=geo.fMatrixPointer;
      fShapePointer=geo.fShapePointer;
      fBomb=geo.fBomb;
      fGeomLevel=geo.fGeomLevel;
      fX=geo.fX;
      fY=geo.fY;
      fZ=geo.fZ;
      for(Int_t i=0; i<kMAXLEVELS; i++) {
         for(Int_t j=0; j<kVectorSize; j++)
            fTranslation[i][j]=geo.fTranslation[i][j];
         for(Int_t j=0; j<kMatrixSize; j++)
            fRotMatrix[i][j]=geo.fRotMatrix[i][j];
         fIsReflection[i]=geo.fIsReflection[i];
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry default destructor.

TGeometry::~TGeometry()
{
   if (!fMaterials) return;
   fMaterials->Delete();
   fMatrices->Delete();
   fShapes->Delete();
   fNodes->Delete();
   delete fMaterials;
   delete fMatrices;
   delete fShapes;
   delete fNodes;
   delete [] fMaterialPointer;
   delete [] fMatrixPointer;
   delete [] fShapePointer;
   fMaterials       = 0;
   fMatrices        = 0;
   fShapes          = 0;
   fNodes           = 0;
   fMaterialPointer = 0;
   fMatrixPointer   = 0;
   fShapePointer    = 0;

   if (gGeometry == this) {
      gGeometry = (TGeometry*) gROOT->GetListOfGeometries()->First();
      if (gGeometry == this)
         gGeometry = (TGeometry*) gROOT->GetListOfGeometries()->After(gGeometry);
   }
   gROOT->GetListOfGeometries()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Browse.

void TGeometry::Browse(TBrowser *b)
{
   if( b ) {
      b->Add( fMaterials, "Materials" );
      b->Add( fMatrices, "Rotation Matrices" );
      b->Add( fShapes, "Shapes" );
      b->Add( fNodes, "Nodes" );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change Current Geometry to this.

void TGeometry::cd(const char *)
{
   gGeometry = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this Geometry.

void TGeometry::Draw(Option_t *option)
{
   TNode *node1 = (TNode*)fNodes->First();
   if (node1) node1->Draw(option);

}

////////////////////////////////////////////////////////////////////////////////
/// Find object in a geometry node, material, etc

TObject *TGeometry::FindObject(const TObject *) const
{
   Error("FindObject","Not yet implemented");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Search object identified by name in the geometry tree

TObject *TGeometry::FindObject(const char *name) const
{
   TObjArray *loc = TGeometry::Get(name);
   if (loc) return loc->At(0);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function called by TROOT to search name in the geometry.
/// Returns a TObjArray containing a pointer to the found object
/// and a pointer to the container where the object was found.

TObjArray *TGeometry::Get(const char *name)
{
   static TObjArray *locs = 0;
   if (!locs) locs = new TObjArray(2);
   TObjArray &loc = *locs;
   loc[0] = 0;
   loc[1] = 0;

   if (!gGeometry) return &loc;

   TObject *temp;
   TObject *where;

   temp  = gGeometry->GetListOfMaterials()->FindObject(name);
   where = gGeometry->GetListOfMaterials();

   if (!temp) {
      temp  = gGeometry->GetListOfShapes()->FindObject(name);
      where = gGeometry->GetListOfShapes();
   }
   if (!temp) {
      temp  = gGeometry->GetListOfMatrices()->FindObject(name);
      where = gGeometry->GetListOfMatrices();
   }
   if (!temp) {
      temp  = gGeometry->GetNode(name);
      where = gGeometry;
   }
   loc[0] = temp;
   loc[1] = where;

   return &loc;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to Material with name.

TMaterial *TGeometry::GetMaterial(const char *name) const
{
   return (TMaterial*)fMaterials->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to Material with number.

TMaterial *TGeometry::GetMaterialByNumber(Int_t number) const
{
   TMaterial *mat;
   if (number < 0 || number >= fMaterials->GetSize()) return 0;
   if (fMaterialPointer)  return fMaterialPointer[number];
   TIter next(fMaterials);
   while ((mat = (TMaterial*) next())) {
      if (mat->GetNumber() == number) return mat;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to node with name in the geometry tree.

TNode *TGeometry::GetNode(const char *name) const
{
   TNode *node= (TNode*)GetListOfNodes()->First();
   if (!node) return 0;
   if (node->TestBit(kNotDeleted))  return node->GetNode(name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to RotMatrix with name.

TRotMatrix *TGeometry::GetRotMatrix(const char *name) const
{
   return (TRotMatrix*)fMatrices->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to RotMatrix with number.

TRotMatrix *TGeometry::GetRotMatrixByNumber(Int_t number) const
{
   TRotMatrix *matrix;
   if (number < 0 || number >= fMatrices->GetSize()) return 0;
   if (fMatrixPointer)  return fMatrixPointer[number];
   TIter next(fMatrices);
   while ((matrix = (TRotMatrix*) next())) {
      if (matrix->GetNumber() == number) return matrix;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to Shape with name.

TShape *TGeometry::GetShape(const char *name) const
{
   return (TShape*)fShapes->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to Shape with number.

TShape *TGeometry::GetShapeByNumber(Int_t number) const
{
   TShape *shape;
   if (number < 0 || number >= fShapes->GetSize()) return 0;
   if (fShapePointer)  return fShapePointer[number];
   TIter next(fShapes);
   while ((shape = (TShape*) next())) {
      if (shape->GetNumber() == number) return shape;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert one point from local system to master reference system.
///
///  Note that before invoking this function, the global rotation matrix
///  and translation vector for this node must have been computed.
///  This is automatically done by the Paint functions.
///  Otherwise TNode::UpdateMatrix should be called before.

void TGeometry::Local2Master(Double_t *local, Double_t *master)
{
   if (GeomLevel()) {
      Double_t x,y,z;
      Double_t bomb = GetBomb();
      Double_t *matrix = &fRotMatrix[GeomLevel()][0];
      x = bomb*fX
        + local[0]*matrix[0]
        + local[1]*matrix[3]
        + local[2]*matrix[6];

      y = bomb*fY
        + local[0]*matrix[1]
        + local[1]*matrix[4]
        + local[2]*matrix[7];

      z = bomb*fZ
        + local[0]*matrix[2]
        + local[1]*matrix[5]
        + local[2]*matrix[8];
      master[0] = x; master[1] = y; master[2] = z;
   }
   else
      for (Int_t i=0;i<3;i++) master[i] = local[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Convert one point from local system to master reference system.
///
///  Note that before invoking this function, the global rotation matrix
///  and translation vector for this node must have been computed.
///  This is automatically done by the Paint functions.
///  Otherwise TNode::UpdateMatrix should be called before.

void TGeometry::Local2Master(Float_t *local, Float_t *master)
{
   if (GeomLevel()) {
      Float_t x,y,z;
      Float_t bomb = GetBomb();

      Double_t *matrix = &fRotMatrix[GeomLevel()][0];

      x = bomb*fX
        + local[0]*matrix[0]
        + local[1]*matrix[3]
        + local[2]*matrix[6];

      y = bomb*fY
        + local[0]*matrix[1]
        + local[1]*matrix[4]
        + local[2]*matrix[7];

      z = bomb*fZ
        + local[0]*matrix[2]
        + local[1]*matrix[5]
        + local[2]*matrix[8];

      master[0] = x; master[1] = y; master[2] = z;
   }
   else
      for (Int_t i=0;i<3;i++) master[i] = local[i];
}

////////////////////////////////////////////////////////////////////////////////
/// List this geometry.

void TGeometry::ls(Option_t *option) const
{
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("m")) {
      Printf("=================List of Materials================");
      fMaterials->ls(option);
   }
   if (opt.Contains("r")) {
      Printf("=================List of RotationMatrices================");
      fMatrices->ls(option);
   }
   if (opt.Contains("s")) {
      Printf("=================List of Shapes==========================");
      fShapes->ls(option);
   }
   if (opt.Contains("n")) {
      Printf("=================List of Nodes===========================");
      fNodes->ls(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert one point from master system to local reference system.
///
///  Note that before invoking this function, the global rotation matrix
///  and translation vector for this node must have been computed.
///  This is automatically done by the Paint functions.
///  Otherwise TNode::UpdateMatrix should be called before.

void TGeometry::Master2Local(Double_t *master, Double_t *local)
{
   if (GeomLevel()) {
      Double_t x,y,z;
      Double_t bomb = GetBomb();
      Double_t *matrix = &fRotMatrix[GeomLevel()][0];

      Double_t xms = master[0] - bomb*fX;
      Double_t yms = master[1] - bomb*fY;
      Double_t zms = master[2] - bomb*fZ;

      x = xms*matrix[0] + yms*matrix[1] + zms*matrix[2];
      y = xms*matrix[3] + yms*matrix[4] + zms*matrix[5];
      z = xms*matrix[6] + yms*matrix[7] + zms*matrix[8];

      local[0] = x; local[1] = y; local[2] = z;
   }
   else
      memcpy(local,master,sizeof(Double_t)* kVectorSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert one point from master system to local reference system.
///
///  Note that before invoking this function, the global rotation matrix
///  and translation vector for this node must have been computed.
///  This is automatically done by the Paint functions.
///  Otherwise TNode::UpdateMatrix should be called before.

void TGeometry::Master2Local(Float_t *master, Float_t *local)
{
   if (GeomLevel()) {
      Float_t x,y,z;
      Float_t bomb = GetBomb();

      Double_t *matrix = &fRotMatrix[GeomLevel()][0];

      Double_t xms = master[0] - bomb*fX;
      Double_t yms = master[1] - bomb*fY;
      Double_t zms = master[2] - bomb*fZ;

      x = xms*matrix[0] + yms*matrix[1] + zms*matrix[2];
      y = xms*matrix[3] + yms*matrix[4] + zms*matrix[5];
      z = xms*matrix[6] + yms*matrix[7] + zms*matrix[8];

      local[0] = x; local[1] = y; local[2] = z;
   }
   else
      memcpy(local,master,sizeof(Float_t)* kVectorSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a node to the current node in this geometry.

void TGeometry::Node(const char *name, const char *title, const char *shapename, Double_t x, Double_t y, Double_t z, const char *matrixname, Option_t *option)
{
   new TNode(name,title,shapename,x,y,z,matrixname,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from a Geometry list.

void TGeometry::RecursiveRemove(TObject *obj)
{
   if (fNodes) fNodes->RecursiveRemove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TGeometry::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         b.ReadClassBuffer(TGeometry::Class(), this, R__v, R__s, R__c);
      } else {
         //====process old versions before automatic schema evolution
         TNamed::Streamer(b);
         fMaterials->Streamer(b);
         fMatrices->Streamer(b);
         fShapes->Streamer(b);
         fNodes->Streamer(b);
         b >> fBomb;
         b.CheckByteCount(R__s, R__c, TGeometry::IsA());
         //====end of old versions
      }
      // Build direct access pointers to individual materials,matrices and shapes
      Int_t i;
      TMaterial *onemat;
      TRotMatrix *onematrix;
      TShape *oneshape;
      Int_t nmat = fMaterials->GetSize();
      if (nmat) fMaterialPointer = new TMaterial* [nmat];
      TIter nextmat(fMaterials);
      i = 0;
      while ((onemat = (TMaterial*) nextmat())) {
         fMaterialPointer[i] = onemat;
         i++;
      }

      Int_t nrot = fMatrices->GetSize();
      if (nrot) fMatrixPointer = new TRotMatrix* [nrot];
      TIter nextmatrix(fMatrices);
      i = 0;
      while ((onematrix = (TRotMatrix*) nextmatrix())) {
         fMatrixPointer[i] = onematrix;
         i++;
      }

      Int_t nsha = fShapes->GetSize();
      if (nsha) fShapePointer = new TShape* [nsha];
      TIter nextshape(fShapes);
      i = 0;
      while ((oneshape = (TShape*) nextshape())) {
         fShapePointer[i] = oneshape;
         i++;
      }

      gROOT->GetListOfGeometries()->Add(this);

      fCurrentNode = (TNode*)GetListOfNodes()->First();
   } else {
      b.WriteClassBuffer(TGeometry::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update global rotation matrix/translation vector for this node
/// this function must be called before invoking Local2Master

void TGeometry::UpdateMatrix(TNode *node)
{
   TNode *nodes[kMAXLEVELS];
   for (Int_t i=0;i<kVectorSize;i++) fTranslation[0][i] = 0;
   for (Int_t i=0;i<kMatrixSize;i++) fRotMatrix[0][i] = 0;
   fRotMatrix[0][0] = 1;   fRotMatrix[0][4] = 1;   fRotMatrix[0][8] = 1;

   fGeomLevel  = 0;
   //build array of parent nodes
   while (node) {
      nodes[fGeomLevel] = node;
      node = node->GetParent();
      fGeomLevel++;
   }
   fGeomLevel--;
   Int_t saveGeomLevel = fGeomLevel;
   //Update matrices in the hierarchy
   for (fGeomLevel=1;fGeomLevel<=saveGeomLevel;fGeomLevel++) {
      node = nodes[fGeomLevel-1];
      UpdateTempMatrix(node->GetX(),node->GetY(),node->GetZ(),node->GetMatrix());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update temp matrix.

void TGeometry::UpdateTempMatrix(Double_t x, Double_t y, Double_t z, TRotMatrix *rotMatrix)
{
   Double_t *matrix = 0;
   Bool_t isReflection = kFALSE;
   if (rotMatrix && rotMatrix->GetType()) {
      matrix = rotMatrix->GetMatrix();
      isReflection = rotMatrix->IsReflection();
   }
   UpdateTempMatrix( x,y,z, matrix,isReflection);
}

////////////////////////////////////////////////////////////////////////////////
/// Update temp matrix.

void TGeometry::UpdateTempMatrix(Double_t x, Double_t y, Double_t z, Double_t *matrix,Bool_t isReflection)
{
   Int_t i=GeomLevel();
   if (i) {
      if(matrix) {
         UpdateTempMatrix(&(fTranslation[i-1][0]),&fRotMatrix[i-1][0]
                          ,x,y,z,matrix
                          ,&fTranslation[i][0],&fRotMatrix[i][0]);
         fX = fTranslation[i][0];
         fY = fTranslation[i][1];
         fZ = fTranslation[i][2];
         fIsReflection[i] = fIsReflection[i-1] ^ isReflection;
      } else {
         fX = fTranslation[i][0] = fTranslation[i-1][0] + x;
         fY = fTranslation[i][1] = fTranslation[i-1][1] + y;
         fZ = fTranslation[i][2] = fTranslation[i-1][2] + z;
      }
   } else {
      fX=fY=fZ=0;
      fIsReflection[0] = kFALSE;
      for (i=0;i<kVectorSize;i++) fTranslation[0][i] = 0;
      for (i=0;i<kMatrixSize;i++) fRotMatrix[0][i] = 0;
      fRotMatrix[0][0] = 1;   fRotMatrix[0][4] = 1;   fRotMatrix[0][8] = 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute new translation vector and global matrix.
///
///  - dx      old translation vector
///  - rmat    old global matrix
///  - x,y,z   offset of new local system with respect to mother
///  - dxnew   new translation vector
///  - rmatnew new global rotation matrix

void TGeometry::UpdateTempMatrix(Double_t *dx,Double_t *rmat
                         , Double_t x, Double_t y, Double_t z, Double_t *matrix
                         , Double_t *dxnew, Double_t *rmatnew)
{
   dxnew[0] = dx[0] + x*rmat[0] + y*rmat[3] + z*rmat[6];
   dxnew[1] = dx[1] + x*rmat[1] + y*rmat[4] + z*rmat[7];
   dxnew[2] = dx[2] + x*rmat[2] + y*rmat[5] + z*rmat[8];

   rmatnew[0] = rmat[0]*matrix[0] + rmat[3]*matrix[1] + rmat[6]*matrix[2];
   rmatnew[1] = rmat[1]*matrix[0] + rmat[4]*matrix[1] + rmat[7]*matrix[2];
   rmatnew[2] = rmat[2]*matrix[0] + rmat[5]*matrix[1] + rmat[8]*matrix[2];
   rmatnew[3] = rmat[0]*matrix[3] + rmat[3]*matrix[4] + rmat[6]*matrix[5];
   rmatnew[4] = rmat[1]*matrix[3] + rmat[4]*matrix[4] + rmat[7]*matrix[5];
   rmatnew[5] = rmat[2]*matrix[3] + rmat[5]*matrix[4] + rmat[8]*matrix[5];
   rmatnew[6] = rmat[0]*matrix[6] + rmat[3]*matrix[7] + rmat[6]*matrix[8];
   rmatnew[7] = rmat[1]*matrix[6] + rmat[4]*matrix[7] + rmat[7]*matrix[8];
   rmatnew[8] = rmat[2]*matrix[6] + rmat[5]*matrix[7] + rmat[8]*matrix[8];
}
