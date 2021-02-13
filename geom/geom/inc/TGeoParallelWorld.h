/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Author: Andrei Gheata   30/06/14

#ifndef ROOT_TGeoParallelWorld
#define ROOT_TGeoParallelWorld

#include "TNamed.h"

// forward declarations
class TGeoManager;
class TGeoPhysicalNode;
class TGeoVolume;

class TGeoParallelWorld : public TNamed
{
protected :
   TGeoManager       *fGeoManager;     // base geometry
   TObjArray         *fPaths;          // array of paths
   Bool_t             fUseOverlaps;    // Activated if user defined overlapping candidates
   Bool_t             fIsClosed;       //! Closed flag
   TGeoVolume        *fVolume;         //! helper volume
   TGeoPhysicalNode  *fLastState;      //! Last PN touched
   TObjArray         *fPhysical;       //! array of physical nodes

   TGeoParallelWorld(const TGeoParallelWorld&) = delete;
   TGeoParallelWorld& operator=(const TGeoParallelWorld&) = delete;

public:
   // constructors
   TGeoParallelWorld() : TNamed(),fGeoManager(0),fPaths(0),fUseOverlaps(kFALSE),fIsClosed(kFALSE),fVolume(0),fLastState(0),fPhysical(0) {}
   TGeoParallelWorld(const char *name, TGeoManager *mgr);

   // destructor
   virtual ~TGeoParallelWorld();
   // API for adding components nodes
   void              AddNode(const char *path);
   // Activate/deactivate  overlap usage
   void              SetUseOverlaps(Bool_t flag) {fUseOverlaps = flag;}
   Bool_t            IsUsingOverlaps() const {return fUseOverlaps;}
   void              ResetOverlaps() const;
   // Adding overlap candidates can highly improve performance.
   void              AddOverlap(TGeoVolume *vol, Bool_t activate=kTRUE);
   void              AddOverlap(const char *volname, Bool_t activate=kTRUE);
   // The normal PW mode (without declaring overlaps) does detect them
   Int_t             PrintDetectedOverlaps() const;

   // Closing a parallel geometry is mandatory
   Bool_t            CloseGeometry();
   // Refresh structures in case of re-alignment
   void              RefreshPhysicalNodes();

   // Navigation interface
   TGeoPhysicalNode *FindNode(Double_t point[3]);
   TGeoPhysicalNode *FindNextBoundary(Double_t point[3], Double_t dir[3], Double_t &step, Double_t stepmax=1.E30);
   Double_t          Safety(Double_t point[3], Double_t safmax=1.E30);

   // Getters
   TGeoManager      *GetGeometry() const {return fGeoManager;}
   Bool_t            IsClosed() const    {return fIsClosed;}
   TGeoVolume       *GetVolume() const   {return fVolume;}

   // Utilities
   void              CheckOverlaps(Double_t ovlp=0.001); // default 10 microns
   void              Draw(Option_t *option);

   ClassDef(TGeoParallelWorld, 3)     // parallel world base class
};

#endif

