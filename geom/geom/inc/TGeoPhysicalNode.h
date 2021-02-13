// @(#)root/geom:$Id$
// Author: Andrei Gheata   17/02/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoPhysicalNode
#define ROOT_TGeoPhysicalNode

#include "TNamed.h"

#include "TAttLine.h"

// forward declarations
class TGeoHMatrix;
class TGeoMatrix;
class TGeoVolume;
class TGeoNode;
class TGeoShape;
class TGeoNavigator;
class TObjArray;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TGeoPhysicalNode - class representing an unique object associated with a //
//   path.                                                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TGeoPhysicalNode : public TNamed,
                         public TAttLine
{
protected:
   Int_t             fLevel;          // depth in the geometry tree
   TObjArray        *fMatrices;       // global transformation matrices
   TObjArray        *fNodes;          // branch of nodes
   TGeoHMatrix      *fMatrixOrig;     // original local matrix of the last node in the path

   TGeoPhysicalNode(const TGeoPhysicalNode&) = delete;
   TGeoPhysicalNode& operator=(const TGeoPhysicalNode&) = delete;

   void              SetAligned(Bool_t flag=kTRUE) {TObject::SetBit(kGeoPNodeAligned,flag);}
   Bool_t            SetPath(const char *path);
   void              SetBranchAsState();

public:
   enum {
      kGeoPNodeFull    = BIT(10),     // full branch is visible (default only last node)
      kGeoPNodeVisible = BIT(11),     // this node is visible (default)
      kGeoPNodeVolAtt  = BIT(12),     // preserve volume attributes (default)
      kGeoPNodeAligned = BIT(13)      // alignment bit
   };

   // constructors
   TGeoPhysicalNode();
   TGeoPhysicalNode(const char *path);
   // destructor
   virtual ~TGeoPhysicalNode();

   Bool_t            Align(TGeoMatrix *newmat=0, TGeoShape *newshape=0, Bool_t check=kFALSE, Double_t ovlp=0.001);
   void              cd() const;
   void              Draw(Option_t *option="");
   Int_t             GetLevel() const {return fLevel;}
   TGeoHMatrix      *GetMatrix(Int_t level=-1) const;
   TGeoHMatrix      *GetOriginalMatrix() const {return fMatrixOrig;}
   TGeoNode         *GetMother(Int_t levup=1) const;
   TGeoNode         *GetNode(Int_t level=-1) const;
   TGeoShape        *GetShape(Int_t level=-1) const;
   TGeoVolume       *GetVolume(Int_t level=-1) const;


   Bool_t            IsAligned() const {return TObject::TestBit(kGeoPNodeAligned);}
   Bool_t            IsMatchingState(TGeoNavigator *nav) const;
   Bool_t            IsVolAttributes() const {return TObject::TestBit(kGeoPNodeVolAtt);}
   Bool_t            IsVisible() const {return TObject::TestBit(kGeoPNodeVisible);}
   Bool_t            IsVisibleFull() const {return TObject::TestBit(kGeoPNodeFull);}

   virtual void      Print(Option_t *option="") const;
   void              Refresh();

   void              SetMatrixOrig(const TGeoMatrix *local);
   void              SetIsVolAtt(Bool_t flag=kTRUE) {TObject::SetBit(kGeoPNodeVolAtt,flag);}
   void              SetVisibility(Bool_t flag=kTRUE)  {TObject::SetBit(kGeoPNodeVisible,flag);}
   void              SetVisibleFull(Bool_t flag=kTRUE) {TObject::SetBit(kGeoPNodeFull,flag);}
   virtual void      Paint(Option_t *option = "");


   ClassDef(TGeoPhysicalNode, 1)               // base class for physical nodes
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TGeoPNEntry - class representing physical node entry having a unique name //
//   associated to a path.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TGeoPNEntry : public TNamed
{
private:
   enum EPNEntryFlags {
      kPNEntryOwnMatrix = BIT(14)
   };
   TGeoPhysicalNode   *fNode;        // Physical node to which this applies
   const TGeoHMatrix  *fMatrix;      // Additional matrix
   TGeoHMatrix        *fGlobalOrig;  // Original global matrix for the linked physical node

protected:
   TGeoPNEntry(const TGeoPNEntry& pne)
     : TNamed(pne), fNode(pne.fNode), fMatrix(NULL), fGlobalOrig(NULL) { }
   TGeoPNEntry& operator=(const TGeoPNEntry& pne)
     {if(this!=&pne) {TNamed::operator=(pne); fNode=pne.fNode; fMatrix=pne.fMatrix;}
     return *this;}

public:
   TGeoPNEntry();
   TGeoPNEntry(const char *unique_name, const char *path);
   virtual ~TGeoPNEntry();

   inline const char   *GetPath() const {return GetTitle();}
   const TGeoHMatrix   *GetMatrix() const {return fMatrix;}
   TGeoHMatrix      *GetMatrixOrig() const {if (fNode) return fNode->GetOriginalMatrix(); else return NULL;};
   TGeoHMatrix      *GetGlobalOrig() const {return fGlobalOrig;}
   TGeoPhysicalNode *GetPhysicalNode() const {return fNode;}
   void              SetMatrix(const TGeoHMatrix *matrix);
   void              SetPhysicalNode(TGeoPhysicalNode *node);

   ClassDef(TGeoPNEntry, 4)                  // a physical node entry with unique name
};

#endif

