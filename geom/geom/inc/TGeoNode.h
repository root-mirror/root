// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoNode
#define ROOT_TGeoNode

#include "TGeoAtt.h"

#include "TNamed.h"

#include "TGeoVolume.h"

#include "TGeoPatternFinder.h"

// forward declarations
class TString;
class TGeoVolume;
class TGeoShape;
class TGeoMedium;
class TGeoMatrix;
class TGeoHMatrix;
class TGeoExtension;

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoNode - base class for logical nodes. They represent volumes        //
//   positioned inside a mother volume                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoNode : public TNamed,
                 public TGeoAtt
{
protected:
   TGeoVolume       *fVolume = nullptr;        // volume associated with this
   TGeoVolume       *fMother = nullptr;        // mother volume
   Int_t             fNumber = 0;              // copy number
   Int_t             fNovlp = 0;               // number of overlaps
   Int_t            *fOverlaps = nullptr;      //[fNovlp] list of indices for overlapping brothers
   TGeoExtension    *fUserExtension = nullptr; //! Transient user-defined extension to volumes
   TGeoExtension    *fFWExtension = nullptr;   //! Transient framework-defined extension to volumes

   void              CopyOverlaps(Int_t *ovlp, Int_t novlp);

private:
   TGeoNode(const TGeoNode&) = delete;
   TGeoNode& operator=(const TGeoNode&) = delete;

public:
   enum {
      kGeoNodeMatrix = BIT(14),
      kGeoNodeOffset = BIT(15),
      kGeoNodeVC     = BIT(16),
      kGeoNodeOverlap = BIT(17),
      kGeoNodeCloned = BIT(18)
   };

   // constructors
   TGeoNode();
   TGeoNode(const TGeoVolume *vol);
   // destructor
   virtual ~TGeoNode();

   void              Browse(TBrowser *b);
   virtual void      cd() const {;}
   void              CheckOverlaps(Double_t ovlp=0.1, Option_t *option=""); // *MENU*
   void              CheckShapes();
   Int_t             CountDaughters(Bool_t unique_volumes=kFALSE);
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   void              Draw(Option_t *option="");
   void              DrawOnly(Option_t *option="");
   void              DrawOverlaps();
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   void              FillIdArray(Int_t &ifree, Int_t &nodeid, Int_t *array) const;
   Int_t             FindNode(const TGeoNode *node, Int_t level);
   virtual Int_t     GetByteCount() const {return 44;}
   TGeoNode         *GetDaughter(Int_t ind) const {return fVolume->GetNode(ind);}
   virtual TGeoMatrix *GetMatrix() const = 0;

   Int_t             GetColour() const {return fVolume->GetLineColor();}
   virtual Int_t     GetIndex() const                    {return 0;}
   virtual TGeoPatternFinder *GetFinder() const          {return 0;}
   TGeoMedium       *GetMedium() const                   {return fVolume->GetMedium();}
   TGeoVolume       *GetMotherVolume() const             {return fMother;}
   Int_t             GetNdaughters() const {return fVolume->GetNdaughters();}
   TObjArray        *GetNodes() const {return fVolume->GetNodes();}
   Int_t             GetNumber() const {return fNumber;}
   Int_t            *GetOverlaps(Int_t &novlp) const {novlp=fNovlp; return fOverlaps;}
   TGeoVolume       *GetVolume() const                   {return fVolume;}
   virtual char     *GetObjectInfo(Int_t px, Int_t py) const;
   virtual Int_t     GetOptimalVoxels() const {return 0;}
   void              InspectNode() const; // *MENU*
   Bool_t            IsCloned() const {return TObject::TestBit(kGeoNodeCloned);}
   virtual Bool_t    IsFolder() const {return (GetNdaughters()?kTRUE:kFALSE);}
   Bool_t            IsOffset() const {return TObject::TestBit(kGeoNodeOffset);}
   Bool_t            IsOnScreen() const; // *MENU*
   Bool_t            IsOverlapping() const {return TObject::TestBit(kGeoNodeOverlap);}
   Bool_t            IsVirtual() const {return TObject::TestBit(kGeoNodeVC);}
   Bool_t            IsVisible() const {return (TGeoAtt::IsVisible() && fVolume->IsVisible());}
   Bool_t            IsVisDaughters() const {return (TGeoAtt::IsVisDaughters() && fVolume->IsVisDaughters());}
   Bool_t            MayOverlap(Int_t iother) const;

   virtual TGeoNode *MakeCopyNode() const {return 0;}
   Double_t          Safety(const Double_t *point, Bool_t in=kTRUE) const;
   void              SaveAttributes(std::ostream &out);
   void              SetCurrentPoint(Double_t x, Double_t y, Double_t z) {fVolume->SetCurrentPoint(x,y,z);}// *MENU*
   void              SetVolume(TGeoVolume *volume)       {fVolume = volume;}
   void              SetNumber(Int_t number)             {fNumber=number;}
   void              SetCloned(Bool_t flag=kTRUE)        {TObject::SetBit(kGeoNodeCloned, flag);}
   void              SetOverlapping(Bool_t flag=kTRUE)   {TObject::SetBit(kGeoNodeOverlap, flag);}
   void              SetVirtual()                        {TObject::SetBit(kGeoNodeVC, kTRUE);}
   void              SetVisibility(Bool_t vis=kTRUE); // *MENU*
   void              SetInvisible()                      {SetVisibility(kFALSE);} // *MENU*
   void              SetAllInvisible()                   {VisibleDaughters(kFALSE);} // *MENU*
   void              SetMotherVolume(TGeoVolume *mother) {fMother = mother;}
   void              SetOverlaps(Int_t *ovlp, Int_t novlp);
   void              SetUserExtension(TGeoExtension *ext);
   void              SetFWExtension(TGeoExtension *ext);
   TGeoExtension    *GetUserExtension() const {return fUserExtension;}
   TGeoExtension    *GetFWExtension() const   {return fFWExtension;}
   TGeoExtension    *GrabUserExtension() const;
   TGeoExtension    *GrabFWExtension() const;

   virtual void      MasterToLocal(const Double_t *master, Double_t *local) const;
   virtual void      MasterToLocalVect(const Double_t *master, Double_t *local) const;
   virtual void      LocalToMaster(const Double_t *local, Double_t *master) const;
   virtual void      LocalToMasterVect(const Double_t *local, Double_t *master) const;

   virtual void      ls(Option_t *option = "") const;
   virtual void      Paint(Option_t *option = "");
   void              PrintCandidates() const; // *MENU*
   void              PrintOverlaps() const; // *MENU*
   void              VisibleDaughters(Bool_t vis=kTRUE); // *MENU*

   ClassDef(TGeoNode, 2)               // base class for all geometry nodes
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoNodeMatrix - node containing a general transformation              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoNodeMatrix : public TGeoNode
{
private:
   TGeoMatrix       *fMatrix = nullptr; // transf. matrix of fNode in fMother system

   TGeoNodeMatrix(const TGeoNodeMatrix& gnm) = delete;
   TGeoNodeMatrix& operator=(const TGeoNodeMatrix& gnm) = delete;

public:
   // constructors
   TGeoNodeMatrix();
   TGeoNodeMatrix(const TGeoVolume *vol, const TGeoMatrix *matrix);
   // destructor
   virtual ~TGeoNodeMatrix();

   virtual Int_t     GetByteCount() const;
   virtual Int_t     GetOptimalVoxels() const;
   virtual Bool_t    IsFolder() const {return kTRUE;}
   virtual TGeoMatrix *GetMatrix() const   {return fMatrix;}
   virtual TGeoNode *MakeCopyNode() const;
   void              SetMatrix(const TGeoMatrix *matrix);

   ClassDef(TGeoNodeMatrix, 1)               // a geometry node in the general case
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoNodeOffset - node containing only an translation offset            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoNodeOffset : public TGeoNode
{
private:
   Double_t          fOffset = 0.;       // X offset for this node with respect to its mother
   Int_t             fIndex = 0;         // index of this node in the division
   TGeoPatternFinder *fFinder = nullptr; // finder for this node

   TGeoNodeOffset(const TGeoNodeOffset&) = delete;
   TGeoNodeOffset& operator=(const TGeoNodeOffset&) = delete;

public:
   // constructors
   TGeoNodeOffset();
   TGeoNodeOffset(const TGeoVolume *vol, Int_t index, Double_t offset);
   // destructor
   virtual ~TGeoNodeOffset();

   virtual void      cd() const           {fFinder->cd(fIndex);}
   Double_t          GetOffset() const {return fOffset;}
   virtual Int_t     GetIndex() const;
   virtual TGeoPatternFinder *GetFinder() const {return fFinder;}
   virtual TGeoMatrix *GetMatrix() const {cd(); return fFinder->GetMatrix();}
   virtual TGeoNode *MakeCopyNode() const;
   void              SetFinder(TGeoPatternFinder *finder) {fFinder = finder;}

   ClassDef(TGeoNodeOffset, 1)      // a geometry node with just an offset
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoIteratorPlugin - Plugin for a TGeoIterator providing the method    //
//                      ProcessNode each time Next is called.             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoIterator;

class TGeoIteratorPlugin : public TObject
{
protected:
   const TGeoIterator *fIterator = nullptr; // Caller iterator
private:
   // No copy
   TGeoIteratorPlugin(const TGeoIteratorPlugin &);
   TGeoIteratorPlugin &operator=(const TGeoIteratorPlugin &);
public:
   TGeoIteratorPlugin() : TObject(),fIterator(0) {}
   virtual ~TGeoIteratorPlugin() {}

   virtual void      ProcessNode() = 0;
   void              SetIterator(const TGeoIterator *iter) {fIterator = iter;}

   ClassDef(TGeoIteratorPlugin, 0)  // ABC for user plugins connecter to a geometry iterator.
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoIterator - iterator for the node tree                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoIterator
{
private:
   TGeoVolume       *fTop = nullptr;        // Top volume of the iterated branch
   Bool_t            fMustResume = kFALSE;  // Private flag to resume from current node.
   Bool_t            fMustStop = kFALSE;    // Private flag to signal that the iterator has finished.
   Int_t             fLevel = 0;            // Current level in the tree
   Int_t             fType = 0;             // Type of iteration
   Int_t            *fArray = nullptr;      // Array of node indices for the current path
   TGeoHMatrix      *fMatrix = nullptr;     // Current global matrix
   TString           fTopName;              // User name for top
   TGeoIteratorPlugin
                    *fPlugin = nullptr;     // User iterator plugin
   Bool_t            fPluginAutoexec = kFALSE; // Plugin automatically executed during next()

   void            IncreaseArray();
protected:
   TGeoIterator() : fTop(0), fMustResume(0), fMustStop(0), fLevel(0), fType(0),
                    fArray(0), fMatrix(0), fTopName(), fPlugin(0), fPluginAutoexec(kFALSE) { }

public:
   TGeoIterator(TGeoVolume *top);
   TGeoIterator(const TGeoIterator &iter);
   virtual          ~TGeoIterator();

   TGeoIterator   &operator=(const TGeoIterator &iter);
   TGeoNode       *operator()();
   TGeoNode       *Next();
   void            Up() { if (fLevel > 0) fLevel--; }

   const TGeoMatrix *GetCurrentMatrix() const;
   Int_t           GetIndex(Int_t i) const {return ((i<=fLevel)?fArray[i]:-1);}
   Int_t           GetLevel() const {return fLevel;}
   TGeoNode       *GetNode(Int_t level) const;
   void            GetPath(TString &path) const;
   TGeoIteratorPlugin
                  *GetUserPlugin() const {return fPlugin;}

   TGeoVolume     *GetTopVolume() const {return fTop;}
   Int_t           GetType() const {return fType;}
   void            Reset(TGeoVolume *top=0);
   void            SetUserPlugin(TGeoIteratorPlugin *plugin);
   void            SetPluginAutoexec(Bool_t mode) {fPluginAutoexec = mode;}
   void            SetType(Int_t type) {fType = type;}
   void            SetTopName(const char* name);
   void            Skip();

   ClassDef(TGeoIterator,0)  //Iterator for geometry.
};

#endif
