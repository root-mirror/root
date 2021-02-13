// @(#)root/geom:$Id$
// Author: Andrei Gheata   25/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoManager
\ingroup Geometry_classes

The manager class for any TGeo geometry. Provides user
interface for geometry creation, navigation, state querying,
visualization, IO, geometry checking and other utilities.

## General architecture

  The ROOT geometry package is a tool designed for building, browsing,
tracking and visualizing a detector geometry. The code is independent from
other external MC for simulation, therefore it does not contain any
constraints related to physics. However, the package defines a number of
hooks for tracking, such as media, materials, magnetic field or track state flags,
in order to allow interfacing to tracking MC's. The final goal is to be
able to use the same geometry for several purposes, such as tracking,
reconstruction or visualization, taking advantage of the ROOT features
related to bookkeeping, I/O, histogramming, browsing and GUI's.

  The geometrical modeler is the most important component of the package and
it provides answers to the basic questions like "Where am I ?" or "How far
from the next boundary ?", but also to more complex ones like "How far from
the closest surface ?" or "Which is the next crossing along a helix ?".

  The architecture of the modeler is a combination between a GEANT-like
containment scheme and a normal CSG binary tree at the level of shapes. An
important common feature of all detector geometry descriptions is the
mother-daughter concept. This is the most natural approach when tracking
is concerned and imposes a set of constraints to the way geometry is defined.
Constructive solid geometry composition is used only in order to create more
complex shapes from an existing set of primitives through boolean operations.
This feature is not implemented yet but in future full definition of boolean
expressions will be supported.

  Practically every geometry defined in GEANT style can be mapped by the modeler.
The basic components used for building the logical hierarchy of the geometry
are called "volumes" and "nodes". Volumes (sometimes called "solids") are fully
defined geometrical objects having a given shape and medium and possibly
containing a list of nodes. Nodes represent just positioned instances of volumes
inside a container volume and they are not directly defined by user. They are
automatically created as a result of adding one volume inside other or dividing
a volume. The geometrical transformation hold by nodes is always defined with
respect to their mother (relative positioning). Reflection matrices are allowed.
All volumes have to be fully aware of their containees when the geometry is
closed. They will build additional structures (voxels) in order to fasten-up
the search algorithms. Finally, nodes can be regarded as bidirectional links
between containers and containees objects.

  The structure defined in this way is a graph structure since volumes are
replicable (same volume can become daughter node of several other volumes),
every volume becoming a branch in this graph. Any volume in the logical graph
can become the actual top volume at run time (see TGeoManager::SetTopVolume()).
All functionalities of the modeler will behave in this case as if only the
corresponding branch starting from this volume is the registered geometry.

\image html geom_graf.jpg

  A given volume can be positioned several times in the geometry. A volume
can be divided according default or user-defined patterns, creating automatically
the list of division nodes inside. The elementary volumes created during the
dividing process follow the same scheme as usual volumes, therefore it is possible
to position further geometrical structures inside or to divide them further more
(see TGeoVolume::Divide()).

  The primitive shapes supported by the package are basically the GEANT3
shapes (see class TGeoShape), arbitrary wedges with eight vertices on two parallel
planes. All basic primitives inherits from class TGeoBBox since the bounding box
of a solid is essential for the tracking algorithms. They also implement the
virtual methods defined in the virtual class TGeoShape (point and segment
classification). User-defined primitives can be directly plugged into the modeler
provided that they override these methods. Composite shapes will be soon supported
by the modeler. In order to build a TGeoCompositeShape, one will have to define
first the primitive components. The object that handle boolean
operations among components is called TGeoBoolCombinator and it has to be
constructed providing a string boolean expression between the components names.


## Example for building a simple geometry

Begin_Macro(source)
../../../tutorials/geom/rootgeom.C
End_Macro

## TGeoManager - the manager class for the geometry package.

  TGeoManager class is embedding all the API needed for building and tracking
a geometry. It defines a global pointer (gGeoManager) in order to be fully
accessible from external code. The mechanism of handling multiple geometries
at the same time will be soon implemented.

  TGeoManager is the owner of all geometry objects defined in a session,
therefore users must not try to control their deletion. It contains lists of
media, materials, transformations, shapes and volumes. Logical nodes (positioned
volumes) are created and destroyed by the TGeoVolume class. Physical
nodes and their global transformations are subjected to a caching mechanism
due to the sometimes very large memory requirements of logical graph expansion.
The caching mechanism is triggered by the total number of physical instances
of volumes and the cache manager is a client of TGeoManager. The manager class
also controls the painter client. This is linked with ROOT graphical libraries
loaded on demand in order to control visualization actions.

## Rules for building a valid geometry

  A given geometry can be built in various ways, but there are mandatory steps
that have to be followed in order to be validated by the modeler. There are
general rules : volumes needs media and shapes in order to be created,
both container and containee volumes must be created before linking them together,
and the relative transformation matrix must be provided. All branches must
have an upper link point otherwise they will not be considered as part of the
geometry. Visibility or tracking properties of volumes can be provided both
at build time or after geometry is closed, but global visualization settings
(see TGeoPainter class) should not be provided at build time, otherwise the
drawing package will be loaded. There is also a list of specific rules :
positioned daughters should not extrude their mother or intersect with sisters
unless this is specified (see TGeoVolume::AddNodeOverlap()), the top volume
(containing all geometry tree) must be specified before closing the geometry
and must not be positioned - it represents the global reference frame. After
building the full geometry tree, the geometry must be closed
(see TGeoManager::CloseGeometry()). Voxelization can be redone per volume after
this process.


  Below is the general scheme of the manager class.

\image html geom_mgr.jpg

## An interactive session

  Provided that a geometry was successfully built and closed (for instance the
previous example $ROOTSYS/tutorials/geom/rootgeom.C ), the manager class will register
itself to ROOT and the logical/physical structures will become immediately browsable.
The ROOT browser will display starting from the geometry folder : the list of
transformations and media, the top volume and the top logical node. These last
two can be fully expanded, any intermediate volume/node in the browser being subject
of direct access context menu operations (right mouse button click). All user
utilities of classes TGeoManager, TGeoVolume and TGeoNode can be called via the
context menu.

\image html geom_browser.jpg

### Drawing the geometry

  Any logical volume can be drawn via TGeoVolume::Draw() member function.
This can be directly accessed from the context menu of the volume object
directly from the browser.
  There are several drawing options that can be set with
TGeoManager::SetVisOption(Int_t opt) method :

#### opt=0
   only the content of the volume is drawn, N levels down (default N=3).
   This is the default behavior. The number of levels to be drawn can be changed
   via TGeoManager::SetVisLevel(Int_t level) method.

\image html geom_frame0.jpg

#### opt=1
   the final leaves (e.g. daughters with no containment) of the branch
   starting from volume are drawn down to the current number of levels.
                                    WARNING : This mode is memory consuming
   depending of the size of geometry, so drawing from top level within this mode
   should be handled with care for expensive geometries. In future there will be
   a limitation on the maximum number of nodes to be visualized.

\image html geom_frame1.jpg

#### opt=2
   only the clicked volume is visualized. This is automatically set by
   TGeoVolume::DrawOnly() method

#### opt=3 - only a given path is visualized. This is automatically set by
   TGeoVolume::DrawPath(const char *path) method

   The current view can be exploded in cartesian, cylindrical or spherical
coordinates :
  TGeoManager::SetExplodedView(Int_t opt). Options may be :
- 0  - default (no bombing)
- 1  - cartesian coordinates. The bomb factor on each axis can be set with
       TGeoManager::SetBombX(Double_t bomb) and corresponding Y and Z.
- 2  - bomb in cylindrical coordinates. Only the bomb factors on Z and R
       are considered
      \image html geom_frameexp.jpg

- 3  - bomb in radial spherical coordinate : TGeoManager::SetBombR()

Volumes themselves support different visualization settings :
   - TGeoVolume::SetVisibility() : set volume visibility.
   - TGeoVolume::VisibleDaughters() : set daughters visibility.
All these actions automatically updates the current view if any.

### Checking the geometry

 Several checking methods are accessible from the volume context menu. They
generally apply only to the visible parts of the drawn geometry in order to
ease geometry checking, and their implementation is in the TGeoChecker class
from the painting package.

#### Checking a given point.
  Can be called from TGeoManager::CheckPoint(Double_t x, Double_t y, Double_t z).
This method is drawing the daughters of the volume containing the point one
level down, printing the path to the deepest physical node holding this point.
It also computes the closest distance to any boundary. The point will be drawn
in red.

\image html geom_checkpoint.jpg

#### Shooting random points.
  Can be called from TGeoVolume::RandomPoints() (context menu function) and
it will draw this volume with current visualization settings. Random points
are generated in the bounding box of the top drawn volume. The points are
classified and drawn with the color of their deepest container. Only points
in visible nodes will be drawn.

\image html geom_random1.jpg


#### Raytracing.
  Can be called from TGeoVolume::RandomRays() (context menu of volumes) and
will shoot rays from a given point in the local reference frame with random
directions. The intersections with displayed nodes will appear as segments
having the color of the touched node. Drawn geometry will be then made invisible
in order to enhance rays.

\image html geom_random2.jpg
*/

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TGeoManager.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TBrowser.h"
#include "TFile.h"
#include "TKey.h"
#include "THashList.h"
#include "TClass.h"
#include "ThreadLocalStorage.h"
#include "TBufferText.h"

#include "TGeoVoxelFinder.h"
#include "TGeoElement.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoPhysicalNode.h"
#include "TGeoPara.h"
#include "TGeoParaboloid.h"
#include "TGeoTube.h"
#include "TGeoEltu.h"
#include "TGeoHype.h"
#include "TGeoCone.h"
#include "TGeoSphere.h"
#include "TGeoArb8.h"
#include "TGeoPgon.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoTorus.h"
#include "TGeoXtru.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"
#include "TGeoBuilder.h"
#include "TVirtualGeoPainter.h"
#include "TPluginManager.h"
#include "TVirtualGeoTrack.h"
#include "TQObject.h"
#include "TMath.h"
#include "TEnv.h"
#include "TGeoParallelWorld.h"
#include "TGeoRegion.h"
#include "TGDMLMatrix.h"
#include "TGeoOpticalSurface.h"

// statics and globals

TGeoManager *gGeoManager = nullptr;

ClassImp(TGeoManager);

std::mutex TGeoManager::fgMutex;
Bool_t TGeoManager::fgLock            = kFALSE;
Bool_t TGeoManager::fgLockNavigators  = kFALSE;
Int_t  TGeoManager::fgVerboseLevel    = 1;
Int_t  TGeoManager::fgMaxLevel        = 1;
Int_t  TGeoManager::fgMaxDaughters    = 1;
Int_t  TGeoManager::fgMaxXtruVert     = 1;
Int_t  TGeoManager::fgNumThreads      = 0;
UInt_t TGeoManager::fgExportPrecision = 17;
TGeoManager::EDefaultUnits TGeoManager::fgDefaultUnits = TGeoManager::kG4Units;
TGeoManager::ThreadsMap_t *TGeoManager::fgThreadId = 0;
static Bool_t gGeometryLocked = kTRUE;

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoManager::TGeoManager()
{
   if (!fgThreadId) fgThreadId = new TGeoManager::ThreadsMap_t;
   if (TClass::IsCallingNew() == TClass::kDummyNew) {
      fTimeCut = kFALSE;
      fTmin = 0.;
      fTmax = 999.;
      fPhiCut = kFALSE;
      fPhimin = 0;
      fPhimax = 360;
      fDrawExtra = kFALSE;
      fStreamVoxels = kFALSE;
      fIsGeomReading = kFALSE;
      fIsGeomCleaning = kFALSE;
      fClosed = kFALSE;
      fLoopVolumes = kFALSE;
      fBits = 0;
      fCurrentNavigator = 0;
      fMaterials = 0;
      fHashPNE = 0;
      fArrayPNE = 0;
      fMatrices = 0;
      fNodes = 0;
      fOverlaps = 0;
      fRegions = 0;
      fNNodes = 0;
      fMaxVisNodes = 10000;
      fVolumes = 0;
      fPhysicalNodes = 0;
      fShapes = 0;
      fGVolumes = 0;
      fGShapes = 0;
      fTracks = 0;
      fMedia = 0;
      fNtracks = 0;
      fNpdg = 0;
      fPdgNames = 0;
      fGDMLMatrices = 0;
      fOpticalSurfaces = 0;
      fSkinSurfaces = 0;
      fBorderSurfaces = 0;
      memset(fPdgId, 0, 1024*sizeof(Int_t));
//   TObjArray            *fNavigators;       //! list of navigators
      fCurrentTrack = 0;
      fCurrentVolume = 0;
      fTopVolume = 0;
      fTopNode = 0;
      fMasterVolume = 0;
      fPainter = 0;
      fActivity = kFALSE;
      fIsNodeSelectable = kFALSE;
      fVisDensity = 0.;
      fVisLevel = 3;
      fVisOption = 1;
      fExplodedView = 0;
      fNsegments = 20;
      fNLevel = 0;
      fUniqueVolumes = 0;
      fClippingShape = 0;
      fMatrixTransform = kFALSE;
      fMatrixReflection = kFALSE;
      fGLMatrix = 0;
      fPaintVolume = 0;
      fUserPaintVolume = 0;
      fElementTable = 0;
      fHashVolumes = 0;
      fHashGVolumes = 0;
      fSizePNEId = 0;
      fNPNEId = 0;
      fKeyPNEId = 0;
      fValuePNEId = 0;
      fMultiThread = kFALSE;
      fRaytraceMode = 0;
      fMaxThreads = 0;
      fUsePWNav = kFALSE;
      fParallelWorld = 0;
      ClearThreadsMap();
   } else {
      Init();
      if (!gGeoIdentity && TClass::IsCallingNew() == TClass::kRealNew) gGeoIdentity = new TGeoIdentity("Identity");
      BuildDefaultMaterials();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoManager::TGeoManager(const char *name, const char *title)
            :TNamed(name, title)
{
   if (!gROOT->GetListOfGeometries()->FindObject(this)) gROOT->GetListOfGeometries()->Add(this);
   if (!gROOT->GetListOfBrowsables()->FindObject(this)) gROOT->GetListOfBrowsables()->Add(this);
   Init();
   gGeoIdentity = new TGeoIdentity("Identity");
   BuildDefaultMaterials();
   if (fgVerboseLevel>0) Info("TGeoManager","Geometry %s, %s created", GetName(), GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize manager class.

void TGeoManager::Init()
{
   if (gGeoManager) {
      Warning("Init","Deleting previous geometry: %s/%s",gGeoManager->GetName(),gGeoManager->GetTitle());
      delete gGeoManager;
      if (fgLock) Fatal("Init", "New geometry created while the old one locked !!!");
   }

   gGeoManager = this;
   if (!fgThreadId) fgThreadId = new TGeoManager::ThreadsMap_t;
   fTimeCut = kFALSE;
   fTmin = 0.;
   fTmax = 999.;
   fPhiCut = kFALSE;
   fPhimin = 0;
   fPhimax = 360;
   fDrawExtra = kFALSE;
   fStreamVoxels = kFALSE;
   fIsGeomReading = kFALSE;
   fIsGeomCleaning = kFALSE;
   fClosed = kFALSE;
   fLoopVolumes = kFALSE;
   fBits = new UChar_t[50000]; // max 25000 nodes per volume
   fCurrentNavigator = 0;
   fHashPNE = new THashList(256,3);
   fArrayPNE = 0;
   fMaterials = new THashList(200,3);
   fMatrices = new TObjArray(256);
   fNodes = new TObjArray(30);
   fOverlaps = new TObjArray(256);
   fRegions = new TObjArray(256);
   fNNodes = 0;
   fMaxVisNodes = 10000;
   fVolumes = new TObjArray(256);
   fPhysicalNodes = new TObjArray(256);
   fShapes = new TObjArray(256);
   fGVolumes = new TObjArray(256);
   fGShapes = new TObjArray(256);
   fTracks = new TObjArray(256);
   fMedia = new THashList(200,3);
   fNtracks = 0;
   fNpdg = 0;
   fPdgNames = 0;
   fGDMLMatrices = new TObjArray();
   fOpticalSurfaces = new TObjArray();
   fSkinSurfaces = new TObjArray();
   fBorderSurfaces = new TObjArray();
   memset(fPdgId, 0, 1024*sizeof(Int_t));
   fCurrentTrack = 0;
   fCurrentVolume = 0;
   fTopVolume = 0;
   fTopNode = 0;
   fMasterVolume = 0;
   fPainter = 0;
   fActivity = kFALSE;
   fIsNodeSelectable = kFALSE;
   fVisDensity = 0.;
   fVisLevel = 3;
   fVisOption = 1;
   fExplodedView = 0;
   fNsegments = 20;
   fNLevel = 0;
   fUniqueVolumes = new TObjArray(256);
   fClippingShape = 0;
   fMatrixTransform = kFALSE;
   fMatrixReflection = kFALSE;
   fGLMatrix = new TGeoHMatrix();
   fPaintVolume = 0;
   fUserPaintVolume = 0;
   fElementTable = 0;
   fHashVolumes = 0;
   fHashGVolumes = 0;
   fSizePNEId = 0;
   fNPNEId = 0;
   fKeyPNEId = 0;
   fValuePNEId = 0;
   fMultiThread = kFALSE;
   fRaytraceMode = 0;
   fMaxThreads = 0;
   fUsePWNav = kFALSE;
   fParallelWorld = 0;
   ClearThreadsMap();
}

////////////////////////////////////////////////////////////////////////////////
///   Destructor

TGeoManager::~TGeoManager()
{
   if (gGeoManager != this) gGeoManager = this;
   fIsGeomCleaning = kTRUE;

   if (gROOT->GetListOfFiles()) { //in case this function is called from TROOT destructor
      gROOT->GetListOfGeometries()->Remove(this);
      gROOT->GetListOfBrowsables()->Remove(this);
   }
//   TSeqCollection *brlist = gROOT->GetListOfBrowsers();
//   TIter next(brlist);
//   TBrowser *browser = 0;
//   while ((browser=(TBrowser*)next())) browser->RecursiveRemove(this);
   ClearThreadsMap();
   ClearThreadData();
   delete TGeoBuilder::Instance(this);
   if (fBits)  delete [] fBits;
   SafeDelete(fNodes);
   SafeDelete(fTopNode);
   if (fOverlaps) {fOverlaps->Delete(); SafeDelete(fOverlaps);}
   if (fRegions) {fRegions->Delete(); SafeDelete(fRegions);}
   if (fMaterials) {fMaterials->Delete(); SafeDelete(fMaterials);}
   SafeDelete(fElementTable);
   if (fMedia) {fMedia->Delete(); SafeDelete(fMedia);}
   if (fHashVolumes) fHashVolumes->Clear("nodelete"); SafeDelete(fHashVolumes);
   if (fHashGVolumes) fHashGVolumes->Clear("nodelete"); SafeDelete(fHashGVolumes);
   if (fHashPNE) {fHashPNE->Delete(); SafeDelete(fHashPNE);}
   if (fArrayPNE) {delete fArrayPNE;}
   if (fVolumes) {fVolumes->Delete(); SafeDelete(fVolumes);}
   if (fShapes) {fShapes->Delete(); SafeDelete( fShapes );}
   if (fPhysicalNodes) {fPhysicalNodes->Delete(); SafeDelete( fPhysicalNodes );}
   if (fMatrices) {fMatrices->Delete(); SafeDelete( fMatrices );}
   if (fTracks) {fTracks->Delete(); SafeDelete( fTracks );}
   SafeDelete( fUniqueVolumes );
   if (fPdgNames) {fPdgNames->Delete(); SafeDelete( fPdgNames );}
   if (fGDMLMatrices) {fGDMLMatrices->Delete(); SafeDelete( fGDMLMatrices );}
   if (fOpticalSurfaces) {fOpticalSurfaces->Delete(); SafeDelete( fOpticalSurfaces );}
   if (fSkinSurfaces) {fSkinSurfaces->Delete(); SafeDelete( fSkinSurfaces );}
   if (fBorderSurfaces) {fBorderSurfaces->Delete(); SafeDelete( fBorderSurfaces );}
   ClearNavigators();
   CleanGarbage();
   SafeDelete( fPainter );
   SafeDelete( fGLMatrix );
   if (fSizePNEId) {
      delete [] fKeyPNEId;
      delete [] fValuePNEId;
   }
   delete fParallelWorld;
   fIsGeomCleaning = kFALSE;
   gGeoIdentity = 0;
   gGeoManager = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a material to the list. Returns index of the material in list.

Int_t TGeoManager::AddMaterial(const TGeoMaterial *material)
{
   return TGeoBuilder::Instance(this)->AddMaterial((TGeoMaterial*)material);
}

////////////////////////////////////////////////////////////////////////////////
/// Add an illegal overlap/extrusion to the list.

Int_t TGeoManager::AddOverlap(const TNamed *ovlp)
{
   Int_t size = fOverlaps->GetEntriesFast();
   fOverlaps->Add((TObject*)ovlp);
   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new region of volumes.
Int_t TGeoManager::AddRegion(TGeoRegion *region)
{
  Int_t size = fRegions->GetEntriesFast();
  fRegions->Add(region);
  return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a user-defined property. Returns true if added, false if existing.

Bool_t TGeoManager::AddProperty(const char* property, Double_t value)
{
   auto pos = fProperties.insert(ConstPropMap_t::value_type(property, value));
   if (!pos.second) {
      Warning("AddProperty", "Property \"%s\" already exists with value %g", property, (pos.first)->second);
      return false;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a user-defined property

Double_t TGeoManager::GetProperty(const char *property, Bool_t *error) const
{
   auto pos = fProperties.find(property);
   if (pos == fProperties.end()) {
      if (error) *error = kTRUE;
      return 0.;
   }
   if (error) *error = kFALSE;
   return pos->second;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a user-defined property from a given index

Double_t TGeoManager::GetProperty(size_t i, TString &name, Bool_t *error) const
{
   // This is a quite inefficient way to access map elements, but needed for the GDML writer to
   if (i >= fProperties.size()) {
      if (error) *error = kTRUE;
      return 0.;
   }
   size_t pos = 0;
   auto it = fProperties.begin();
   while (pos < i) { ++it; ++pos; }
   if (error) *error = kFALSE;
   name = (*it).first;
   return (*it).second;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a matrix to the list. Returns index of the matrix in list.

Int_t TGeoManager::AddTransformation(const TGeoMatrix *matrix)
{
   return TGeoBuilder::Instance(this)->AddTransformation((TGeoMatrix*)matrix);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a shape to the list. Returns index of the shape in list.

Int_t TGeoManager::AddShape(const TGeoShape *shape)
{
   return TGeoBuilder::Instance(this)->AddShape((TGeoShape*)shape);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a track to the list of tracks. Use this for primaries only. For secondaries,
/// add them to the parent track. The method create objects that are registered
/// to the analysis manager but have to be cleaned-up by the user via ClearTracks().

Int_t TGeoManager::AddTrack(Int_t id, Int_t pdgcode, TObject *particle)
{
   Int_t index = fNtracks;
   fTracks->AddAtAndExpand(GetGeomPainter()->AddTrack(id,pdgcode,particle),fNtracks++);
   return index;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a track to the list of tracks

Int_t TGeoManager::AddTrack(TVirtualGeoTrack *track)
{
   Int_t index = fNtracks;
   fTracks->AddAtAndExpand(track,fNtracks++);
   return index;
}

////////////////////////////////////////////////////////////////////////////////
/// Makes a primary track but do not attach it to the list of tracks. The track
/// can be attached as daughter to another one with TVirtualGeoTrack::AddTrack

TVirtualGeoTrack *TGeoManager::MakeTrack(Int_t id, Int_t pdgcode, TObject *particle)
{
   TVirtualGeoTrack *track = GetGeomPainter()->AddTrack(id,pdgcode,particle);
   return track;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a volume to the list. Returns index of the volume in list.

Int_t TGeoManager::AddVolume(TGeoVolume *volume)
{
   if (!volume) {
      Error("AddVolume", "invalid volume");
      return -1;
   }
   Int_t uid = fUniqueVolumes->GetEntriesFast();
   if (!uid) uid++;
   if (!fCurrentVolume) {
      fCurrentVolume = volume;
      fUniqueVolumes->AddAtAndExpand(volume,uid);
   } else {
      if (!strcmp(volume->GetName(), fCurrentVolume->GetName())) {
         uid = fCurrentVolume->GetNumber();
      } else {
         fCurrentVolume = volume;
         Int_t olduid = GetUID(volume->GetName());
         if (olduid<0) {
            fUniqueVolumes->AddAtAndExpand(volume,uid);
         } else {
            uid = olduid;
         }
      }
   }
   volume->SetNumber(uid);
   if (!fHashVolumes) {
      fHashVolumes = new THashList(256);
      fHashGVolumes = new THashList(256);
   }
   TObjArray *list = fVolumes;
   if (!volume->GetShape() || volume->IsRunTime() || volume->IsVolumeMulti()) {
      list = fGVolumes;
      fHashGVolumes->Add(volume);
   } else {
      fHashVolumes->Add(volume);
   }
   Int_t index = list->GetEntriesFast();
   list->AddAtAndExpand(volume,index);
   return uid;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a navigator in the list of navigators. If it is the first one make it
/// current navigator.

TGeoNavigator *TGeoManager::AddNavigator()
{
   if (fMultiThread) { TGeoManager::ThreadId(); fgMutex.lock(); }
   std::thread::id threadId = std::this_thread::get_id();
   NavigatorsMap_t::const_iterator it = fNavigators.find(threadId);
   TGeoNavigatorArray *array = 0;
   if (it != fNavigators.end()) array = it->second;
   else {
      array = new TGeoNavigatorArray(this);
      fNavigators.insert(NavigatorsMap_t::value_type(threadId, array));
   }
   TGeoNavigator *nav = array->AddNavigator();
   if (fClosed) nav->GetCache()->BuildInfoBranch();
   if (fMultiThread) fgMutex.unlock();
   return nav;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current navigator for the calling thread.

TGeoNavigator *TGeoManager::GetCurrentNavigator() const
{
   TTHREAD_TLS(TGeoNavigator*) tnav = 0;
   if (!fMultiThread) return fCurrentNavigator;
   TGeoNavigator *nav = tnav; // TTHREAD_TLS_GET(TGeoNavigator*,tnav);
   if (nav) return nav;
   std::thread::id threadId = std::this_thread::get_id();
   NavigatorsMap_t::const_iterator it = fNavigators.find(threadId);
   if (it == fNavigators.end()) return 0;
   TGeoNavigatorArray *array = it->second;
   nav = array->GetCurrentNavigator();
   tnav = nav; // TTHREAD_TLS_SET(TGeoNavigator*,tnav,nav);
   return nav;
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of navigators for the calling thread.

TGeoNavigatorArray *TGeoManager::GetListOfNavigators() const
{
   std::thread::id threadId = std::this_thread::get_id();
   NavigatorsMap_t::const_iterator it = fNavigators.find(threadId);
   if (it == fNavigators.end()) return 0;
   TGeoNavigatorArray *array = it->second;
   return array;
}

////////////////////////////////////////////////////////////////////////////////
/// Switch to another existing navigator for the calling thread.

Bool_t TGeoManager::SetCurrentNavigator(Int_t index)
{
   std::thread::id threadId = std::this_thread::get_id();
   NavigatorsMap_t::const_iterator it = fNavigators.find(threadId);
   if (it == fNavigators.end()) {
      Error("SetCurrentNavigator", "No navigator defined for this thread\n");
      std::cout << "  thread id: " << threadId << std::endl;
      return kFALSE;
   }
   TGeoNavigatorArray *array = it->second;
   TGeoNavigator *nav = array->SetCurrentNavigator(index);
   if (!nav) {
      Error("SetCurrentNavigator", "Navigator %d not existing for this thread\n", index);
      std::cout << "  thread id: " << threadId << std::endl;
      return kFALSE;
   }
   if (!fMultiThread) fCurrentNavigator = nav;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the lock for navigators.

void TGeoManager::SetNavigatorsLock(Bool_t flag)
{
   fgLockNavigators = flag;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear all navigators.

void TGeoManager::ClearNavigators()
{
   if (fMultiThread) fgMutex.lock();
   TGeoNavigatorArray *arr = 0;
   for (NavigatorsMap_t::iterator it = fNavigators.begin();
        it != fNavigators.end(); ++it) {
      arr = (*it).second;
      if (arr) delete arr;
   }
   fNavigators.clear();
   if (fMultiThread) fgMutex.unlock();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear a single navigator.

void TGeoManager::RemoveNavigator(const TGeoNavigator *nav)
{
   if (fMultiThread) fgMutex.lock();
   for (NavigatorsMap_t::iterator it = fNavigators.begin(); it != fNavigators.end(); ++it) {
      TGeoNavigatorArray *arr = (*it).second;
      if (arr) {
         if ((TGeoNavigator*)arr->Remove((TObject*)nav)) {
            delete nav;
            if (!arr->GetEntries()) fNavigators.erase(it);
            if (fMultiThread) fgMutex.unlock();
            return;
         }
      }
   }
   Error("Remove navigator", "Navigator %p not found", nav);
   if (fMultiThread) fgMutex.unlock();
}

////////////////////////////////////////////////////////////////////////////////
/// Set maximum number of threads for navigation.

void TGeoManager::SetMaxThreads(Int_t nthreads)
{
   if (!fClosed) {
      Error("SetMaxThreads", "Cannot set maximum number of threads before closing the geometry");
      return;
   }
   if (!fMultiThread) {
      ROOT::EnableThreadSafety();
      std::thread::id threadId = std::this_thread::get_id();
      NavigatorsMap_t::const_iterator it = fNavigators.find(threadId);
      if (it != fNavigators.end()) {
         TGeoNavigatorArray *array = it->second;
         fNavigators.erase(it);
         fNavigators.insert(NavigatorsMap_t::value_type(threadId, array));
      }
   }
   if (fMaxThreads) {
      ClearThreadsMap();
      ClearThreadData();
   }
   fMaxThreads = nthreads+1;
   if (fMaxThreads>0) {
      fMultiThread = kTRUE;
      CreateThreadData();
   }
}

////////////////////////////////////////////////////////////////////////////////

void TGeoManager::ClearThreadData() const
{
   if (!fMaxThreads) return;
   fgMutex.lock();
   TIter next(fVolumes);
   TGeoVolume *vol;
   while ((vol=(TGeoVolume*)next())) vol->ClearThreadData();
   fgMutex.unlock();
}

////////////////////////////////////////////////////////////////////////////////
/// Create thread private data for all geometry objects.

void TGeoManager::CreateThreadData() const
{
   if (!fMaxThreads) return;
   fgMutex.lock();
   TIter next(fVolumes);
   TGeoVolume *vol;
   while ((vol=(TGeoVolume*)next())) vol->CreateThreadData(fMaxThreads);
   fgMutex.unlock();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the current map of threads. This will be filled again by the calling
/// threads via ThreadId calls.

void TGeoManager::ClearThreadsMap()
{
   if (gGeoManager && !gGeoManager->IsMultiThread()) return;
   fgMutex.lock();
   if (!fgThreadId->empty()) fgThreadId->clear();
   fgNumThreads = 0;
   fgMutex.unlock();
}

////////////////////////////////////////////////////////////////////////////////
/// Translates the current thread id to an ordinal number. This can be used to
/// manage data which is specific for a given thread.

Int_t TGeoManager::ThreadId()
{
   TTHREAD_TLS(Int_t) tid = -1;
   Int_t ttid = tid; // TTHREAD_TLS_GET(Int_t,tid);
   if (ttid > -1) return ttid;
   if (gGeoManager && !gGeoManager->IsMultiThread()) return 0;
   std::thread::id threadId = std::this_thread::get_id();
   TGeoManager::ThreadsMapIt_t it = fgThreadId->find(threadId);
   if (it != fgThreadId->end()) return it->second;
   // Map needs to be updated.
   fgMutex.lock();
   (*fgThreadId)[threadId] = fgNumThreads;
   tid = fgNumThreads; // TTHREAD_TLS_SET(Int_t,tid,fgNumThreads);
   ttid = fgNumThreads++;
   fgMutex.unlock();
   return ttid;
}

////////////////////////////////////////////////////////////////////////////////
/// Describe how to browse this object.

void TGeoManager::Browse(TBrowser *b)
{
   if (!b) return;
   if (fMaterials) b->Add(fMaterials, "Materials");
   if (fMedia)     b->Add(fMedia,     "Media");
   if (fMatrices)  b->Add(fMatrices, "Local transformations");
   if (fOverlaps)  b->Add(fOverlaps, "Illegal overlaps");
   if (fTracks)    b->Add(fTracks,   "Tracks");
   if (fMasterVolume) b->Add(fMasterVolume, "Master Volume", fMasterVolume->IsVisible());
   if (fTopVolume) b->Add(fTopVolume, "Top Volume", fTopVolume->IsVisible());
   if (fTopNode)   b->Add(fTopNode);
   TString browserImp(gEnv->GetValue("Browser.Name", "TRootBrowserLite"));
   TQObject::Connect(browserImp.Data(), "Checked(TObject*,Bool_t)",
                     "TGeoManager", this, "SetVisibility(TObject*,Bool_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Append a pad for this geometry.

void TGeoManager::Edit(Option_t *option) {
   AppendPad("");
   GetGeomPainter()->EditGeometry(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set visibility for a volume.

void TGeoManager::SetVisibility(TObject *obj, Bool_t vis)
{
   if(obj->IsA() == TGeoVolume::Class()) {
      TGeoVolume *vol = (TGeoVolume *) obj;
      vol->SetVisibility(vis);
   } else {
      if (obj->InheritsFrom(TGeoNode::Class())) {
         TGeoNode *node = (TGeoNode *) obj;
         node->SetVisibility(vis);
      } else return;
   }
   GetGeomPainter()->ModifiedPad(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the new 'bombed' translation vector according current exploded view mode.

void TGeoManager::BombTranslation(const Double_t *tr, Double_t *bombtr)
{
   if (fPainter) fPainter->BombTranslation(tr, bombtr);
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the new 'unbombed' translation vector according current exploded view mode.

void TGeoManager::UnbombTranslation(const Double_t *tr, Double_t *bombtr)
{
   if (fPainter) fPainter->UnbombTranslation(tr, bombtr);
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Backup the current state without affecting the cache stack.

void TGeoManager::DoBackupState()
{
   GetCurrentNavigator()->DoBackupState();
}

////////////////////////////////////////////////////////////////////////////////
/// Restore a backed-up state without affecting the cache stack.

void TGeoManager::DoRestoreState()
{
   GetCurrentNavigator()->DoRestoreState();
}

////////////////////////////////////////////////////////////////////////////////
/// Register a matrix to the list of matrices. It will be cleaned-up at the
/// destruction TGeoManager.

void TGeoManager::RegisterMatrix(const TGeoMatrix *matrix)
{
   return TGeoBuilder::Instance(this)->RegisterMatrix((TGeoMatrix*)matrix);
}

////////////////////////////////////////////////////////////////////////////////
/// Replaces all occurrences of VORIG with VNEW in the geometry tree. The volume VORIG
/// is not replaced from the list of volumes, but all node referencing it will reference
/// VNEW instead. Returns number of occurrences changed.

Int_t TGeoManager::ReplaceVolume(TGeoVolume *vorig, TGeoVolume *vnew)
{
   Int_t nref = 0;
   if (!vorig || !vnew) return nref;
   TGeoMedium *morig = vorig->GetMedium();
   Bool_t checkmed = kFALSE;
   if (morig) checkmed = kTRUE;
   TGeoMedium *mnew = vnew->GetMedium();
   // Try to limit the damage produced by incorrect usage.
   if (!mnew && !vnew->IsAssembly()) {
      Error("ReplaceVolume","Replacement volume %s has no medium and it is not an assembly",
             vnew->GetName());
      return nref;
   }
   if (mnew && checkmed) {
      if (mnew->GetId() != morig->GetId())
         Warning("ReplaceVolume","Replacement volume %s has different medium than original volume %s",
                 vnew->GetName(), vorig->GetName());
      checkmed = kFALSE;
   }

   // Medium checking now performed only if replacement is an assembly and old volume a real one.
   // Check result is dependent on positioning.
   Int_t nvol = fVolumes->GetEntriesFast();
   Int_t i,j,nd;
   Int_t ierr = 0;
   TGeoVolume *vol;
   TGeoNode *node;
   TGeoVoxelFinder *voxels;
   for (i=0; i<nvol; i++) {
      vol = (TGeoVolume*)fVolumes->At(i);
      if (!vol) continue;
      if (vol==vorig || vol==vnew) continue;
      nd = vol->GetNdaughters();
      for (j=0; j<nd; j++) {
         node = vol->GetNode(j);
         if (node->GetVolume() == vorig) {
            if (checkmed) {
               mnew = node->GetMotherVolume()->GetMedium();
               if (mnew && mnew->GetId()!=morig->GetId()) ierr++;
            }
            nref++;
            if (node->IsOverlapping()) {
               node->SetOverlapping(kFALSE);
               Info("ReplaceVolume","%s replaced with assembly and declared NON-OVERLAPPING!",node->GetName());
            }
            node->SetVolume(vnew);
            voxels = node->GetMotherVolume()->GetVoxels();
            if (voxels) voxels->SetNeedRebuild();
         } else {
            if (node->GetMotherVolume() == vorig) {
               nref++;
               node->SetMotherVolume(vnew);
               if (node->IsOverlapping()) {
                  node->SetOverlapping(kFALSE);
                  Info("ReplaceVolume","%s inside substitute assembly %s declared NON-OVERLAPPING!",node->GetName(),vnew->GetName());
               }
            }
         }
      }
   }
   if (ierr) Warning("ReplaceVolume", "Volumes should not be replaced with assemblies if they are positioned in containers having a different medium ID.\n %i occurrences for assembly replacing volume %s",
                     ierr, vorig->GetName());
   return nref;
}

////////////////////////////////////////////////////////////////////////////////
/// Transform all volumes named VNAME to assemblies. The volumes must be virtual.

Int_t TGeoManager::TransformVolumeToAssembly(const char *vname)
{
   TGeoVolume *toTransform = FindVolumeFast(vname);
   if (!toTransform) {
      Warning("TransformVolumeToAssembly", "Volume %s not found", vname);
      return 0;
   }
   Int_t index = fVolumes->IndexOf(toTransform);
   Int_t count = 0;
   Int_t indmax = fVolumes->GetEntries();
   Bool_t replace = kTRUE;
   TGeoVolume *transformed;
   while (index<indmax) {
      if (replace) {
         replace = kFALSE;
         transformed = TGeoVolumeAssembly::MakeAssemblyFromVolume(toTransform);
         if (transformed) {
            ReplaceVolume(toTransform, transformed);
            count++;
         } else {
            if (toTransform->IsAssembly())
               Warning("TransformVolumeToAssembly", "Volume %s already assembly", toTransform->GetName());
            if (!toTransform->GetNdaughters())
               Warning("TransformVolumeToAssembly", "Volume %s has no daughters, cannot transform", toTransform->GetName());
            if (toTransform->IsVolumeMulti())
               Warning("TransformVolumeToAssembly", "Volume %s divided, cannot transform", toTransform->GetName());
         }
      }
      index++;
      if (index >= indmax) return count;
      toTransform = (TGeoVolume*)fVolumes->At(index);
      if (!strcmp(toTransform->GetName(),vname)) replace = kTRUE;
   }
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new volume by dividing an existing one (GEANT3 like)
///
/// Divides MOTHER into NDIV divisions called NAME
/// along axis IAXIS starting at coordinate value START
/// and having size STEP. The created volumes will have tracking
/// media ID=NUMED (if NUMED=0 -> same media as MOTHER)
///    The behavior of the division operation can be triggered using OPTION :
///
/// OPTION (case insensitive) :
///  - N  - divide all range in NDIV cells (same effect as STEP<=0) (GSDVN in G3)
///  - NX - divide range starting with START in NDIV cells          (GSDVN2 in G3)
///  - S  - divide all range with given STEP. NDIV is computed and divisions will be centered
///           in full range (same effect as NDIV<=0)                (GSDVS, GSDVT in G3)
///  - SX - same as DVS, but from START position.                   (GSDVS2, GSDVT2 in G3)

TGeoVolume *TGeoManager::Division(const char *name, const char *mother, Int_t iaxis,
                                  Int_t ndiv, Double_t start, Double_t step, Int_t numed, Option_t *option)
{
   return TGeoBuilder::Instance(this)->Division(name, mother, iaxis, ndiv, start, step, numed, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Create rotation matrix named 'mat<index>'.
///
///  - index    rotation matrix number
///  - theta1   polar angle for axis X
///  - phi1     azimuthal angle for axis X
///  - theta2   polar angle for axis Y
///  - phi2     azimuthal angle for axis Y
///  - theta3   polar angle for axis Z
///  - phi3     azimuthal angle for axis Z
///

void TGeoManager::Matrix(Int_t index, Double_t theta1, Double_t phi1,
                         Double_t theta2, Double_t phi2,
                         Double_t theta3, Double_t phi3)
{
   TGeoBuilder::Instance(this)->Matrix(index, theta1, phi1, theta2, phi2, theta3, phi3);
}

////////////////////////////////////////////////////////////////////////////////
/// Create material with given A, Z and density, having an unique id.

TGeoMaterial *TGeoManager::Material(const char *name, Double_t a, Double_t z, Double_t dens, Int_t uid,Double_t radlen, Double_t intlen)
{
   return TGeoBuilder::Instance(this)->Material(name, a, z, dens, uid, radlen, intlen);

}

////////////////////////////////////////////////////////////////////////////////
/// Create mixture OR COMPOUND IMAT as composed by THE BASIC nelem
/// materials defined by arrays A,Z and WMAT, having an unique id.

TGeoMaterial *TGeoManager::Mixture(const char *name, Float_t *a, Float_t *z, Double_t dens,
                                   Int_t nelem, Float_t *wmat, Int_t uid)
{
   return TGeoBuilder::Instance(this)->Mixture(name, a, z, dens, nelem, wmat, uid);
}

////////////////////////////////////////////////////////////////////////////////
/// Create mixture OR COMPOUND IMAT as composed by THE BASIC nelem
/// materials defined by arrays A,Z and WMAT, having an unique id.

TGeoMaterial *TGeoManager::Mixture(const char *name, Double_t *a, Double_t *z, Double_t dens,
                                   Int_t nelem, Double_t *wmat, Int_t uid)
{
   return TGeoBuilder::Instance(this)->Mixture(name, a, z, dens, nelem, wmat, uid);
}

////////////////////////////////////////////////////////////////////////////////
/// Create tracking medium
///
///  - numed      tracking medium number assigned
///  - name      tracking medium name
///  - nmat      material number
///  - isvol     sensitive volume flag
///  - ifield    magnetic field
///  - fieldm    max. field value (kilogauss)
///  - tmaxfd    max. angle due to field (deg/step)
///  - stemax    max. step allowed
///  - deemax    max. fraction of energy lost in a step
///  - epsil     tracking precision (cm)
///  - stmin     min. step due to continuous processes (cm)
///
///  - ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
///  - ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
///     performed with g3helix; ifield = 3 if tracking performed with g3helx3.
///

TGeoMedium *TGeoManager::Medium(const char *name, Int_t numed, Int_t nmat, Int_t isvol,
                                Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                                Double_t stemax, Double_t deemax, Double_t epsil,
                                Double_t stmin)
{
   return TGeoBuilder::Instance(this)->Medium(name, numed, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a node called <name_nr> pointing to the volume called <name>
/// as daughter of the volume called <mother> (gspos). The relative matrix is
/// made of : a translation (x,y,z) and a rotation matrix named <matIROT>.
/// In case npar>0, create the volume to be positioned in mother, according
/// its actual parameters (gsposp).
///  - NAME   Volume name
///  - NUMBER Copy number of the volume
///  - MOTHER Mother volume name
///  - X      X coord. of the volume in mother ref. sys.
///  - Y      Y coord. of the volume in mother ref. sys.
///  - Z      Z coord. of the volume in mother ref. sys.
///  - IROT   Rotation matrix number w.r.t. mother ref. sys.
///  - ISONLY ONLY/MANY flag

void TGeoManager::Node(const char *name, Int_t nr, const char *mother,
                       Double_t x, Double_t y, Double_t z, Int_t irot,
                       Bool_t isOnly, Float_t *upar, Int_t npar)
{
   TGeoBuilder::Instance(this)->Node(name, nr, mother, x, y, z, irot, isOnly, upar, npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a node called <name_nr> pointing to the volume called <name>
/// as daughter of the volume called <mother> (gspos). The relative matrix is
/// made of : a translation (x,y,z) and a rotation matrix named <matIROT>.
/// In case npar>0, create the volume to be positioned in mother, according
/// its actual parameters (gsposp).
///  - NAME   Volume name
///  - NUMBER Copy number of the volume
///  - MOTHER Mother volume name
///  - X      X coord. of the volume in mother ref. sys.
///  - Y      Y coord. of the volume in mother ref. sys.
///  - Z      Z coord. of the volume in mother ref. sys.
///  - IROT   Rotation matrix number w.r.t. mother ref. sys.
///  - ISONLY ONLY/MANY flag

void TGeoManager::Node(const char *name, Int_t nr, const char *mother,
                       Double_t x, Double_t y, Double_t z, Int_t irot,
                       Bool_t isOnly, Double_t *upar, Int_t npar)
{
   TGeoBuilder::Instance(this)->Node(name, nr, mother, x, y, z, irot, isOnly, upar, npar);

}

////////////////////////////////////////////////////////////////////////////////
/// Create a volume in GEANT3 style.
///  - NAME   Volume name
///  - SHAPE  Volume type
///  - NMED   Tracking medium number
///  - NPAR   Number of shape parameters
///  - UPAR   Vector containing shape parameters

TGeoVolume *TGeoManager::Volume(const char *name, const char *shape, Int_t nmed,
                                Float_t *upar, Int_t npar)
{
   return TGeoBuilder::Instance(this)->Volume(name, shape, nmed, upar, npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a volume in GEANT3 style.
///  - NAME   Volume name
///  - SHAPE  Volume type
///  - NMED   Tracking medium number
///  - NPAR   Number of shape parameters
///  - UPAR   Vector containing shape parameters

TGeoVolume *TGeoManager::Volume(const char *name, const char *shape, Int_t nmed,
                                Double_t *upar, Int_t npar)
{
   return TGeoBuilder::Instance(this)->Volume(name, shape, nmed, upar, npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Assigns uid's for all materials,media and matrices.

void TGeoManager::SetAllIndex()
{
   Int_t index = 1;
   TIter next(fMaterials);
   TGeoMaterial *mater;
   while ((mater=(TGeoMaterial*)next())) {
      mater->SetUniqueID(index++);
      mater->ResetBit(TGeoMaterial::kMatSavePrimitive);
   }
   index = 1;
   TIter next1(fMedia);
   TGeoMedium *med;
   while ((med=(TGeoMedium*)next1())) {
      med->SetUniqueID(index++);
      med->ResetBit(TGeoMedium::kMedSavePrimitive);
   }
   index = 1;
   TIter next2(fShapes);
   TGeoShape *shape;
   while ((shape=(TGeoShape*)next2())) {
      shape->SetUniqueID(index++);
      if (shape->IsComposite()) ((TGeoCompositeShape*)shape)->GetBoolNode()->RegisterMatrices();
   }

   TIter next3(fMatrices);
   TGeoMatrix *matrix;
   while ((matrix=(TGeoMatrix*)next3())) {
      matrix->RegisterYourself();
   }
   TIter next4(fMatrices);
   index = 1;
   while ((matrix=(TGeoMatrix*)next4())) {
      matrix->SetUniqueID(index++);
      matrix->ResetBit(TGeoMatrix::kGeoSavePrimitive);
   }
   TIter next5(fVolumes);
   TGeoVolume *vol;
   while ((vol=(TGeoVolume*)next5())) vol->UnmarkSaved();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset all attributes to default ones. Default attributes for visualization
/// are those defined before closing the geometry.

void TGeoManager::ClearAttributes()
{
   if (gPad) delete gPad;
   gPad = 0;
   SetVisOption(0);
   SetVisLevel(3);
   SetExplodedView(0);
   SetBombFactors();
   if (!gStyle) return;
   TIter next(fVolumes);
   TGeoVolume *vol = 0;
   while ((vol=(TGeoVolume*)next())) {
      if (!vol->IsVisTouched()) continue;
      vol->SetVisTouched(kFALSE);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Closing geometry implies checking the geometry validity, fixing shapes
/// with negative parameters (run-time shapes)building the cache manager,
/// voxelizing all volumes, counting the total number of physical nodes and
/// registering the manager class to the browser.

void TGeoManager::CloseGeometry(Option_t *option)
{
   if (fClosed) {
      Warning("CloseGeometry", "geometry already closed");
      return;
   }
   if (!fMasterVolume) {
      Error("CloseGeometry","you MUST call SetTopVolume() first !");
      return;
   }
   if (!gROOT->GetListOfGeometries()->FindObject(this)) gROOT->GetListOfGeometries()->Add(this);
   if (!gROOT->GetListOfBrowsables()->FindObject(this)) gROOT->GetListOfBrowsables()->Add(this);
//   TSeqCollection *brlist = gROOT->GetListOfBrowsers();
//   TIter next(brlist);
//   TBrowser *browser = 0;
//   while ((browser=(TBrowser*)next())) browser->Refresh();
   TString opt(option);
   opt.ToLower();
//   Bool_t dummy = opt.Contains("d");
   Bool_t nodeid = opt.Contains("i");
   // Create a geometry navigator if not present
   TGeoNavigator *nav = 0;
   Int_t nnavigators = 0;
   // Check if the geometry is streamed from file
   if (fIsGeomReading) {
      if (fgVerboseLevel>0) Info("CloseGeometry","Geometry loaded from file...");
      gGeoIdentity=(TGeoIdentity *)fMatrices->At(0);
      if (!fElementTable) fElementTable = new TGeoElementTable(200);
      if (!fTopNode) {
         if (!fMasterVolume) {
            Error("CloseGeometry", "Master volume not streamed");
            return;
         }
         SetTopVolume(fMasterVolume);
         if (fStreamVoxels && fgVerboseLevel>0) Info("CloseGeometry","Voxelization retrieved from file");
      }
      // Create a geometry navigator if not present
      if (!GetCurrentNavigator()) fCurrentNavigator = AddNavigator();
      nnavigators = GetListOfNavigators()->GetEntriesFast();
      Voxelize("ALL");
      CountLevels();
      for (Int_t i=0; i<nnavigators; i++) {
         nav = (TGeoNavigator*)GetListOfNavigators()->At(i);
         nav->GetCache()->BuildInfoBranch();
         if (nodeid) nav->GetCache()->BuildIdArray();
      }
      if (!fHashVolumes) {
         Int_t nvol = fVolumes->GetEntriesFast();
         Int_t ngvol = fGVolumes->GetEntriesFast();
         fHashVolumes = new THashList(nvol+1);
         fHashGVolumes = new THashList(ngvol+1);
         Int_t i;
         for (i=0; i<ngvol; i++) fHashGVolumes->AddLast(fGVolumes->At(i));
         for (i=0; i<nvol; i++) fHashVolumes->AddLast(fVolumes->At(i));
      }
      fClosed = kTRUE;
      if (fParallelWorld) {
         if (fgVerboseLevel>0) Info("CloseGeometry","Recreating parallel world %s ...",fParallelWorld->GetName());
         fParallelWorld->CloseGeometry();
      }

      if (fgVerboseLevel>0) Info("CloseGeometry","%i nodes/ %i volume UID's in %s", fNNodes, fUniqueVolumes->GetEntriesFast()-1, GetTitle());
      if (fgVerboseLevel>0) Info("CloseGeometry","----------------modeler ready----------------");
      return;
   }

   // Create a geometry navigator if not present
   if (!GetCurrentNavigator()) fCurrentNavigator = AddNavigator();
   nnavigators = GetListOfNavigators()->GetEntriesFast();
   SelectTrackingMedia();
   CheckGeometry();
   if (fgVerboseLevel>0) Info("CloseGeometry","Counting nodes...");
   fNNodes = CountNodes();
   fNLevel = fMasterVolume->CountNodes(1,3)+1;
   if (fNLevel<30) fNLevel = 100;

//   BuildIdArray();
   Voxelize("ALL");
   if (fgVerboseLevel>0) Info("CloseGeometry","Building cache...");
   CountLevels();
   for (Int_t i=0; i<nnavigators; i++) {
      nav = (TGeoNavigator*)GetListOfNavigators()->At(i);
      nav->GetCache()->BuildInfoBranch();
      if (nodeid) nav->GetCache()->BuildIdArray();
   }
   fClosed = kTRUE;
   if (fgVerboseLevel>0) {
      Info("CloseGeometry","%i nodes/ %i volume UID's in %s", fNNodes, fUniqueVolumes->GetEntriesFast()-1, GetTitle());
      Info("CloseGeometry","----------------modeler ready----------------");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the list of overlaps.

void TGeoManager::ClearOverlaps()
{
   if (fOverlaps) {
      fOverlaps->Delete();
      delete fOverlaps;
   }
   fOverlaps = new TObjArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a shape from the list of shapes.

void TGeoManager::ClearShape(const TGeoShape *shape)
{
   if (fShapes->FindObject(shape)) fShapes->Remove((TGeoShape*)shape);
   delete shape;
}

////////////////////////////////////////////////////////////////////////////////
/// Clean temporary volumes and shapes from garbage collection.

void TGeoManager::CleanGarbage()
{
   if (!fGVolumes && !fGShapes) return;
   Int_t i,nentries;
   if (fGVolumes) {
      nentries = fGVolumes->GetEntries();
      TGeoVolume *vol = 0;
      for (i=0; i<nentries; i++) {
         vol=(TGeoVolume*)fGVolumes->At(i);
         if (vol) vol->SetFinder(0);
      }
      fGVolumes->Delete();
      delete fGVolumes;
      fGVolumes = 0;
   }
   if (fGShapes) {
      fGShapes->Delete();
      delete fGShapes;
      fGShapes = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change current path to point to the node having this id.
/// Node id has to be in range : 0 to fNNodes-1 (no check for performance reasons)

void TGeoManager::CdNode(Int_t nodeid)
{
   GetCurrentNavigator()->CdNode(nodeid);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the unique ID of the current node.

Int_t TGeoManager::GetCurrentNodeId() const
{
   return GetCurrentNavigator()->GetCurrentNodeId();
}

////////////////////////////////////////////////////////////////////////////////
/// Make top level node the current node. Updates the cache accordingly.
/// Determine the overlapping state of current node.

void TGeoManager::CdTop()
{
   GetCurrentNavigator()->CdTop();
}

////////////////////////////////////////////////////////////////////////////////
/// Go one level up in geometry. Updates cache accordingly.
/// Determine the overlapping state of current node.

void TGeoManager::CdUp()
{
   GetCurrentNavigator()->CdUp();
}

////////////////////////////////////////////////////////////////////////////////
/// Make a daughter of current node current. Can be called only with a valid
/// daughter index (no check). Updates cache accordingly.

void TGeoManager::CdDown(Int_t index)
{
   GetCurrentNavigator()->CdDown(index);
}

////////////////////////////////////////////////////////////////////////////////
/// Do a cd to the node found next by FindNextBoundary

void TGeoManager::CdNext()
{
   GetCurrentNavigator()->CdNext();
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the tree of nodes starting from fTopNode according to pathname.
/// Changes the path accordingly.

Bool_t TGeoManager::cd(const char *path)
{
   return GetCurrentNavigator()->cd(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a geometry path is valid without changing the state of the current navigator.

Bool_t TGeoManager::CheckPath(const char *path) const
{
   return GetCurrentNavigator()->CheckPath(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert all reflections in geometry to normal rotations + reflected shapes.

void TGeoManager::ConvertReflections()
{
   if (!fTopNode) return;
   if (fgVerboseLevel>0) Info("ConvertReflections", "Converting reflections in: %s - %s ...", GetName(), GetTitle());
   TGeoIterator next(fTopVolume);
   TGeoNode *node;
   TGeoNodeMatrix *nodematrix;
   TGeoMatrix *matrix, *mclone;
   TGeoVolume *reflected;
   while ((node=next())) {
      matrix = node->GetMatrix();
      if (matrix->IsReflection()) {
//         printf("%s before\n", node->GetName());
//         matrix->Print();
         mclone = new TGeoCombiTrans(*matrix);
         mclone->RegisterYourself();
         // Reflect just the rotation component
         mclone->ReflectZ(kFALSE, kTRUE);
         nodematrix = (TGeoNodeMatrix*)node;
         nodematrix->SetMatrix(mclone);
//         printf("%s after\n", node->GetName());
//         node->GetMatrix()->Print();
         reflected = node->GetVolume()->MakeReflectedVolume();
         node->SetVolume(reflected);
      }
   }
   if (fgVerboseLevel>0) Info("ConvertReflections", "Done");
}

////////////////////////////////////////////////////////////////////////////////
/// Count maximum number of nodes per volume, maximum depth and maximum
/// number of xtru vertices.

void TGeoManager::CountLevels()
{
   if (!fTopNode) {
      Error("CountLevels", "Top node not defined.");
      return;
   }
   TGeoIterator next(fTopVolume);
   Bool_t fixrefs = fIsGeomReading && (fMasterVolume->GetRefCount()==1);
   if (fMasterVolume->GetRefCount()>1) fMasterVolume->Release();
   if (fgVerboseLevel>1 && fixrefs) Info("CountLevels", "Fixing volume reference counts");
   TGeoNode *node;
   Int_t maxlevel = 1;
   Int_t maxnodes = fTopVolume->GetNdaughters();
   Int_t maxvertices = 1;
   while ((node=next())) {
      if (fixrefs) {
         node->GetVolume()->Grab();
         for (Int_t ibit=10; ibit<14; ibit++) {
            node->SetBit(BIT(ibit+4), node->TestBit(BIT(ibit)));
//            node->ResetBit(BIT(ibit)); // cannot overwrite old crap for reproducibility
         }
      }
      if (node->GetVolume()->GetVoxels()) {
         if (node->GetNdaughters()>maxnodes) maxnodes = node->GetNdaughters();
      }
      if (next.GetLevel()>maxlevel) maxlevel = next.GetLevel();
      if (node->GetVolume()->GetShape()->IsA()==TGeoXtru::Class()) {
         TGeoXtru *xtru = (TGeoXtru*)node->GetVolume()->GetShape();
         if (xtru->GetNvert()>maxvertices) maxvertices = xtru->GetNvert();
      }
   }
   fgMaxLevel = maxlevel;
   fgMaxDaughters = maxnodes;
   fgMaxXtruVert = maxvertices;
   if (fgVerboseLevel>0) Info("CountLevels", "max level = %d, max placements = %d", fgMaxLevel, fgMaxDaughters);
}

////////////////////////////////////////////////////////////////////////////////
/// Count the total number of nodes starting from a volume, nlevels down.

Int_t TGeoManager::CountNodes(const TGeoVolume *vol, Int_t nlevels, Int_t option)
{
   TGeoVolume *top;
   if (!vol) {
      top = fTopVolume;
   } else {
      top = (TGeoVolume*)vol;
   }
   Int_t count = top->CountNodes(nlevels, option);
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default angles for a given view.

void TGeoManager::DefaultAngles()
{
   if (fPainter) fPainter->DefaultAngles();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw current point in the same view.

void TGeoManager::DrawCurrentPoint(Int_t color)
{
   if (fPainter) fPainter->DrawCurrentPoint(color);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw animation of tracks

void TGeoManager::AnimateTracks(Double_t tmin, Double_t tmax, Int_t nframes, Option_t *option)
{
   SetAnimateTracks();
   GetGeomPainter();
   if (tmin<0 || tmin>=tmax || nframes<1) return;
   Double_t *box = fPainter->GetViewBox();
   box[0] = box[1] = box[2] = 0;
   box[3] = box[4] = box[5] = 100;
   Double_t dt = (tmax-tmin)/Double_t(nframes);
   Double_t delt = 2E-9;
   Double_t t = tmin;
   Int_t i, j;
   TString opt(option);
   Bool_t save = kFALSE, geomanim=kFALSE;
   TString fname;
   if (opt.Contains("/S")) save = kTRUE;

   if (opt.Contains("/G")) geomanim = kTRUE;
   SetTminTmax(0,0);
   DrawTracks(opt.Data());
   Double_t start[6], end[6];
   Double_t dd[6] = {0,0,0,0,0,0};
   Double_t dlat=0, dlong=0, dpsi=0;
   if (geomanim) {
      fPainter->EstimateCameraMove(tmin+5*dt, tmin+15*dt, start, end);
      for (i=0; i<3; i++) {
         start[i+3] = 20 + 1.3*start[i+3];
         end[i+3] = 20 + 0.9*end[i+3];
      }
      for (i=0; i<6; i++) {
         dd[i] = (end[i]-start[i])/10.;
      }
      memcpy(box, start, 6*sizeof(Double_t));
      fPainter->GetViewAngles(dlong,dlat,dpsi);
      dlong = (-206-dlong)/Double_t(nframes);
      dlat  = (126-dlat)/Double_t(nframes);
      dpsi  = (75-dpsi)/Double_t(nframes);
      fPainter->GrabFocus();
   }

   for (i=0; i<nframes; i++) {
      if (t-delt<0) SetTminTmax(t-delt,t);
      else gGeoManager->SetTminTmax(t-delt,t);
      if (geomanim) {
         for (j=0; j<6; j++) box[j]+=dd[j];
         fPainter->GrabFocus(1,dlong,dlat,dpsi);
      } else {
         ModifiedPad();
      }
      if (save) {
         fname = TString::Format("anim%04d.gif", i);
         gPad->Print(fname);
      }
      t += dt;
   }
   SetAnimateTracks(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw tracks over the geometry, according to option. By default, only
/// primaries are drawn. See TGeoTrack::Draw() for additional options.

void TGeoManager::DrawTracks(Option_t *option)
{
   TVirtualGeoTrack *track;
   //SetVisLevel(1);
   //SetVisOption(1);
   SetAnimateTracks();
   for (Int_t i=0; i<fNtracks; i++) {
      track = GetTrack(i);
      if (track) track->Draw(option);
   }
   SetAnimateTracks(kFALSE);
   ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw current path

void TGeoManager::DrawPath(const char *path, Option_t *option)
{
   if (!fTopVolume) return;
   fTopVolume->SetVisBranch();
   GetGeomPainter()->DrawPath(path, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw random points in the bounding box of a volume.

void TGeoManager::RandomPoints(const TGeoVolume *vol, Int_t npoints, Option_t *option)
{
   GetGeomPainter()->RandomPoints((TGeoVolume*)vol, npoints, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Check time of finding "Where am I" for n points.

void TGeoManager::Test(Int_t npoints, Option_t *option)
{
   GetGeomPainter()->Test(npoints, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry overlap checker based on sampling.

void TGeoManager::TestOverlaps(const char* path)
{
   GetGeomPainter()->TestOverlaps(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill volume names of current branch into an array.

void TGeoManager::GetBranchNames(Int_t *names) const
{
   GetCurrentNavigator()->GetBranchNames(names);
}

////////////////////////////////////////////////////////////////////////////////
/// Get name for given pdg code;

const char *TGeoManager::GetPdgName(Int_t pdg) const
{
   static char defaultname[5] = { "XXX" };
   if (!fPdgNames || !pdg) return defaultname;
   for (Int_t i=0; i<fNpdg; i++) {
      if (fPdgId[i]==pdg) return fPdgNames->At(i)->GetName();
   }
   return defaultname;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a name for a particle having a given pdg.

void TGeoManager::SetPdgName(Int_t pdg, const char *name)
{
   if (!pdg) return;
   if (!fPdgNames) {
      fPdgNames = new TObjArray(1024);
   }
   if (!strcmp(name, GetPdgName(pdg))) return;
   // store pdg name
   if (fNpdg>1023) {
      Warning("SetPdgName", "No more than 256 different pdg codes allowed");
      return;
   }
   fPdgId[fNpdg] = pdg;
   TNamed *pdgname = new TNamed(name, "");
   fPdgNames->AddAtAndExpand(pdgname, fNpdg++);
}

////////////////////////////////////////////////////////////////////////////////
/// Get GDML matrix with a given name;

TGDMLMatrix *TGeoManager::GetGDMLMatrix(const char *name) const
{
   return (TGDMLMatrix*)fGDMLMatrices->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Add GDML matrix;
void TGeoManager::AddGDMLMatrix(TGDMLMatrix *mat)
{
   if (GetGDMLMatrix(mat->GetName())) {
      Error("AddGDMLMatrix", "Matrix %s already added to manager", mat->GetName());
      return;
   }
   fGDMLMatrices->Add(mat);
}

////////////////////////////////////////////////////////////////////////////////
/// Get optical surface with a given name;

TGeoOpticalSurface *TGeoManager::GetOpticalSurface(const char *name) const
{
   return (TGeoOpticalSurface*)fOpticalSurfaces->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Add optical surface;
void TGeoManager::AddOpticalSurface(TGeoOpticalSurface *optsurf)
{
   if (GetOpticalSurface(optsurf->GetName())) {
      Error("AddOpticalSurface", "Surface %s already added to manager", optsurf->GetName());
      return;
   }
   fOpticalSurfaces->Add(optsurf);
}

////////////////////////////////////////////////////////////////////////////////
/// Get skin surface with a given name;

TGeoSkinSurface *TGeoManager::GetSkinSurface(const char *name) const
{
   return (TGeoSkinSurface*)fSkinSurfaces->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Add skin surface;
void TGeoManager::AddSkinSurface(TGeoSkinSurface *surf)
{
   if (GetSkinSurface(surf->GetName())) {
      Error("AddSkinSurface", "Surface %s already added to manager", surf->GetName());
      return;
   }
   fSkinSurfaces->Add(surf);
}

////////////////////////////////////////////////////////////////////////////////
/// Get border surface with a given name;

TGeoBorderSurface *TGeoManager::GetBorderSurface(const char *name) const
{
   return (TGeoBorderSurface*)fBorderSurfaces->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Add border surface;
void TGeoManager::AddBorderSurface(TGeoBorderSurface *surf)
{
   if (GetBorderSurface(surf->GetName())) {
      Error("AddBorderSurface", "Surface %s already added to manager", surf->GetName());
      return;
   }
   fBorderSurfaces->Add(surf);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill node copy numbers of current branch into an array.

void TGeoManager::GetBranchNumbers(Int_t *copyNumbers, Int_t *volumeNumbers) const
{
   GetCurrentNavigator()->GetBranchNumbers(copyNumbers, volumeNumbers);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill node copy numbers of current branch into an array.

void TGeoManager::GetBranchOnlys(Int_t *isonly) const
{
   GetCurrentNavigator()->GetBranchOnlys(isonly);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve cartesian and radial bomb factors.

void TGeoManager::GetBombFactors(Double_t &bombx, Double_t &bomby, Double_t &bombz, Double_t &bombr) const
{
   if (fPainter) {
      fPainter->GetBombFactors(bombx, bomby, bombz, bombr);
      return;
   }
   bombx = bomby = bombz = bombr = 1.3;
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum number of daughters of a volume used in the geometry.

Int_t TGeoManager::GetMaxDaughters()
{
   return fgMaxDaughters;
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum number of levels used in the geometry.

Int_t TGeoManager::GetMaxLevels()
{
   return fgMaxLevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum number of vertices for an xtru shape used.

Int_t TGeoManager::GetMaxXtruVert()
{
   return fgMaxXtruVert;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of threads that were set to use geometry.

Int_t TGeoManager::GetNumThreads()
{
   return fgNumThreads;
}

////////////////////////////////////////////////////////////////////////////////
/// Return stored current matrix (global matrix of the next touched node).

TGeoHMatrix *TGeoManager::GetHMatrix()
{
   if (!GetCurrentNavigator()) return NULL;
   return GetCurrentNavigator()->GetHMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current depth to which geometry is drawn.

Int_t TGeoManager::GetVisLevel() const
{
   return fVisLevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current depth to which geometry is drawn.

Int_t TGeoManager::GetVisOption() const
{
   return fVisOption;
}

////////////////////////////////////////////////////////////////////////////////
/// Find level of virtuality of current overlapping node (number of levels
/// up having the same tracking media.

Int_t TGeoManager::GetVirtualLevel()
{
   return GetCurrentNavigator()->GetVirtualLevel();
}

////////////////////////////////////////////////////////////////////////////////
/// Search the track hierarchy to find the track with the
/// given id
///
/// if 'primsFirst' is true, then:
/// first tries TGeoManager::GetTrackOfId, then does a
/// recursive search if that fails. this would be faster
/// if the track is somehow known to be a primary

TVirtualGeoTrack *TGeoManager::FindTrackWithId(Int_t id) const
{
   TVirtualGeoTrack* trk = 0;
   trk = GetTrackOfId(id);
   if (trk) return trk;
   // need recursive search
   TIter next(fTracks);
   TVirtualGeoTrack* prim;
   while ((prim = (TVirtualGeoTrack*)next())) {
      trk = prim->FindTrackWithId(id);
      if (trk) return trk;
   }
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Get track with a given ID.

TVirtualGeoTrack *TGeoManager::GetTrackOfId(Int_t id) const
{
   TVirtualGeoTrack *track;
   for (Int_t i=0; i<fNtracks; i++) {
      if ((track = (TVirtualGeoTrack *)fTracks->UncheckedAt(i))) {
         if (track->GetId() == id) return track;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get parent track with a given ID.

TVirtualGeoTrack *TGeoManager::GetParentTrackOfId(Int_t id) const
{
   TVirtualGeoTrack *track = fCurrentTrack;
   while ((track=track->GetMother())) {
      if (track->GetId()==id) return track;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get index for track id, -1 if not found.

Int_t TGeoManager::GetTrackIndex(Int_t id) const
{
   TVirtualGeoTrack *track;
   for (Int_t i=0; i<fNtracks; i++) {
      if ((track = (TVirtualGeoTrack *)fTracks->UncheckedAt(i))) {
         if (track->GetId() == id) return i;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Go upwards the tree until a non-overlapping node

Bool_t TGeoManager::GotoSafeLevel()
{
   return GetCurrentNavigator()->GotoSafeLevel();
}

////////////////////////////////////////////////////////////////////////////////
/// Go upwards the tree until a non-overlapping node

Int_t TGeoManager::GetSafeLevel() const
{
   return GetCurrentNavigator()->GetSafeLevel();
}

////////////////////////////////////////////////////////////////////////////////
/// Set default volume colors according to A of material

void TGeoManager::DefaultColors()
{
   const Int_t nmax = 110;
   Int_t col[nmax];
   for (Int_t i=0;i<nmax;i++) col[i] = kGray;

   //here we should create a new TColor with the same rgb as in the default
   //ROOT colors used below
   col[ 3] = kYellow-10;
   col[ 4] = col[ 5] = kGreen-10;
   col[ 6] = col[ 7] = kBlue-7;
   col[ 8] = col[ 9] = kMagenta-3;
   col[10] = col[11] = kRed-10;
   col[12] = kGray+1;
   col[13] = kBlue-10;
   col[14] = kOrange+7;
   col[16] = kYellow+1;
   col[20] = kYellow-10;
   col[24] = col[25] = col[26] = kBlue-8;
   col[29] = kOrange+9;
   col[79] = kOrange-2;

   TGeoVolume *vol;
   TIter next(fVolumes);
   while ((vol=(TGeoVolume*)next())) {
      TGeoMedium *med = vol->GetMedium();
      if (!med) continue;
      TGeoMaterial *mat = med->GetMaterial();
      Int_t matZ = (Int_t)mat->GetZ();
      vol->SetLineColor(col[matZ]);
      if (mat->GetDensity()<0.1) vol->SetTransparency(60);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from the current point. This represent the distance
/// from POINT to the closest boundary.

Double_t TGeoManager::Safety(Bool_t inside)
{
   return GetCurrentNavigator()->Safety(inside);
}

////////////////////////////////////////////////////////////////////////////////
/// Set volume attributes in G3 style.

void TGeoManager::SetVolumeAttribute(const char *name, const char *att, Int_t val)
{
   TGeoVolume *volume;
   Bool_t all = kFALSE;
   if (strstr(name,"*")) all=kTRUE;
   Int_t ivo=0;
   TIter next(fVolumes);
   TString chatt = att;
   chatt.ToLower();
   while ((volume=(TGeoVolume*)next())) {
      if (strcmp(volume->GetName(), name) && !all) continue;
      ivo++;
      if (chatt.Contains("colo")) volume->SetLineColor(val);
      if (chatt.Contains("lsty")) volume->SetLineStyle(val);
      if (chatt.Contains("lwid")) volume->SetLineWidth(val);
      if (chatt.Contains("fill")) volume->SetFillColor(val);
      if (chatt.Contains("seen")) volume->SetVisibility(val);
   }
   TIter next1(fGVolumes);
   while ((volume=(TGeoVolume*)next1())) {
      if (strcmp(volume->GetName(), name) && !all) continue;
      ivo++;
      if (chatt.Contains("colo")) volume->SetLineColor(val);
      if (chatt.Contains("lsty")) volume->SetLineStyle(val);
      if (chatt.Contains("lwid")) volume->SetLineWidth(val);
      if (chatt.Contains("fill")) volume->SetFillColor(val);
      if (chatt.Contains("seen")) volume->SetVisibility(val);
   }
   if (!ivo) {
      Warning("SetVolumeAttribute","volume: %s does not exist",name);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set factors that will "bomb" all translations in cartesian and cylindrical coordinates.

void TGeoManager::SetBombFactors(Double_t bombx, Double_t bomby, Double_t bombz, Double_t bombr)
{
   if (fPainter) fPainter->SetBombFactors(bombx, bomby, bombz, bombr);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a user-defined shape as clipping for ray tracing.

void TGeoManager::SetClippingShape(TGeoShape *shape)
{
   TVirtualGeoPainter *painter = GetGeomPainter();
   if (shape) {
      if (fClippingShape && (fClippingShape!=shape)) ClearShape(fClippingShape);
      fClippingShape = shape;
   }
   painter->SetClippingShape(shape);
}

////////////////////////////////////////////////////////////////////////////////
/// set the maximum number of visible nodes.

void TGeoManager::SetMaxVisNodes(Int_t maxnodes) {
   fMaxVisNodes = maxnodes;
   if (maxnodes>0 && fgVerboseLevel>0)
      Info("SetMaxVisNodes","Automatic visible depth for %d visible nodes", maxnodes);
   if (!fPainter) return;
   fPainter->CountVisibleNodes();
   Int_t level = fPainter->GetVisLevel();
   if (level != fVisLevel) fVisLevel = level;
}

////////////////////////////////////////////////////////////////////////////////
/// make top volume visible on screen

void TGeoManager::SetTopVisible(Bool_t vis) {
   GetGeomPainter();
   fPainter->SetTopVisible(vis);
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a given node to be checked for overlaps. Any other overlaps will be ignored.

void TGeoManager::SetCheckedNode(TGeoNode *node) {
   GetGeomPainter()->SetCheckedNode(node);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of points to be generated on the shape outline when checking
/// for overlaps.

void TGeoManager::SetNmeshPoints(Int_t npoints)
{
   GetGeomPainter()->SetNmeshPoints(npoints);
}

////////////////////////////////////////////////////////////////////////////////
/// set drawing mode :
///  - option=0 (default) all nodes drawn down to vislevel
///  - option=1           leaves and nodes at vislevel drawn
///  - option=2           path is drawn
///  - option=4           visibility changed

void TGeoManager::SetVisOption(Int_t option) {
   if ((option>=0) && (option<3)) fVisOption=option;
   if (fPainter) fPainter->SetVisOption(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set visualization option (leaves only OR all volumes)

void TGeoManager::ViewLeaves(Bool_t flag)
{
   if (flag) SetVisOption(1);
   else      SetVisOption(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set density threshold. Volumes with densities lower than this become
/// transparent.

void TGeoManager::SetVisDensity(Double_t density)
{
   fVisDensity = density;
   if (fPainter) fPainter->ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// set default level down to which visualization is performed

void TGeoManager::SetVisLevel(Int_t level) {
   if (level>0) {
      fVisLevel = level;
      fMaxVisNodes = 0;
      if (fgVerboseLevel>0)
         Info("SetVisLevel","Automatic visible depth disabled");
      if (fPainter) fPainter->CountVisibleNodes();
   } else {
      SetMaxVisNodes();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sort overlaps by decreasing overlap distance. Extrusions comes first.

void TGeoManager::SortOverlaps()
{
   fOverlaps->Sort();
}

////////////////////////////////////////////////////////////////////////////////
/// Optimize voxelization type for all volumes. Save best choice in a macro.

void TGeoManager::OptimizeVoxels(const char *filename)
{
   if (!fTopNode) {
      Error("OptimizeVoxels","Geometry must be closed first");
      return;
   }
   std::ofstream out;
   TString fname = filename;
   if (fname.IsNull()) fname = "tgeovox.C";
   out.open(fname, std::ios::out);
   if (!out.good()) {
      Error("OptimizeVoxels", "cannot open file");
      return;
   }
   // write header
   TDatime t;
   TString sname(fname);
   sname.ReplaceAll(".C", "");
   out << sname.Data()<<"()"<<std::endl;
   out << "{" << std::endl;
   out << "//=== Macro generated by ROOT version "<< gROOT->GetVersion()<<" : "<<t.AsString()<<std::endl;
   out << "//=== Voxel optimization for " << GetTitle() << " geometry"<<std::endl;
   out << "//===== <run this macro JUST BEFORE closing the geometry>"<<std::endl;
   out << "   TGeoVolume *vol = 0;"<<std::endl;
   out << "   // parse all voxelized volumes"<<std::endl;
   TGeoVolume *vol = 0;
   Bool_t cyltype;
   TIter next(fVolumes);
   while ((vol=(TGeoVolume*)next())) {
      if (!vol->GetVoxels()) continue;
      out<<"   vol = gGeoManager->GetVolume(\""<<vol->GetName()<<"\");"<<std::endl;
      cyltype = vol->OptimizeVoxels();
      if (cyltype) {
         out<<"   vol->SetCylVoxels();"<<std::endl;
      } else {
         out<<"   vol->SetCylVoxels(kFALSE);"<<std::endl;
      }
   }
   out << "}" << std::endl;
   out.close();
}
////////////////////////////////////////////////////////////////////////////////
/// Parse a string boolean expression and do a syntax check. Find top
/// level boolean operator and returns its type. Fill the two
/// substrings to which this operator applies. The returned integer is :
///  - -1 : parse error
///  - 0 : no boolean operator
///  - 1 : union - represented as '+' in expression
///  - 2 : difference (subtraction) - represented as '-' in expression
///  - 3 : intersection - represented as '*' in expression.
/// Parentheses should be used to avoid ambiguities. For instance :
///  - A+B-C will be interpreted as (A+B)-C which is not the same as A+(B-C)
/// eliminate not needed parentheses

Int_t TGeoManager::Parse(const char *expr, TString &expr1, TString &expr2, TString &expr3)
{
   TString startstr(expr);
   Int_t len = startstr.Length();
   Int_t i;
   TString e0 = "";
   expr3 = "";
   // eliminate blanks
   for (i=0; i< len; i++) {
      if (startstr(i)==' ') continue;
      e0 += startstr(i, 1);
   }
   Int_t level = 0;
   Int_t levmin = 999;
   Int_t boolop = 0;
   Int_t indop = 0;
   Int_t iloop = 1;
   Int_t lastop = 0;
   Int_t lastdp = 0;
   Int_t lastpp = 0;
   Bool_t foundmat = kFALSE;
   // check/eliminate parentheses
   while (iloop==1) {
      iloop = 0;
      lastop = 0;
      lastdp = 0;
      lastpp = 0;
      len = e0.Length();
      for (i=0; i<len; i++) {
         if (e0(i)=='(') {
            if (!level) iloop++;
            level++;
            continue;
         }
         if  (e0(i)==')') {
            level--;
            if (level==0) lastpp=i;
            continue;
         }
         if ((e0(i)=='+') || (e0(i)=='-') || (e0(i)=='*')) {
            lastop = i;
            if (level<levmin) {
               levmin = level;
               indop = i;
            }
            continue;
         }
         if  ((e0(i)==':') && (level==0)) {
            lastdp = i;
            continue;
         }
      }
      if (level!=0) {
         if (gGeoManager) gGeoManager->Error("Parse","parentheses does not match");
         return -1;
      }
      if (iloop==1 && (e0(0)=='(') && (e0(len-1)==')')) {
         // eliminate extra parentheses
         e0=e0(1, len-2);
         continue;
      }
      if (foundmat) break;
      if (((lastop==0) && (lastdp>0)) || ((lastpp>0) && (lastdp>lastpp) && (indop<lastpp))) {
         expr3 = e0(lastdp+1, len-lastdp);
         e0=e0(0, lastdp);
         foundmat = kTRUE;
         iloop = 1;
         continue;
      } else break;
   }
   // loop expression and search parentheses/operators
   levmin = 999;
   for (i=0; i<len; i++) {
      if (e0(i)=='(') {
         level++;
         continue;
      }
      if  (e0(i)==')') {
         level--;
         continue;
      }
      // Take LAST operator at lowest level (revision 28/07/08)
      if (level<=levmin) {
         if (e0(i)=='+') {
            boolop = 1; // union
            levmin = level;
            indop = i;
         }
         if (e0(i)=='-') {
            boolop = 2; // difference
            levmin = level;
            indop = i;
         }
         if (e0(i)=='*') {
            boolop = 3; // intersection
            levmin = level;
            indop = i;
         }
      }
   }
   if (indop==0) {
      expr1=e0;
      return indop;
   }
   expr1 = e0(0, indop);
   expr2 = e0(indop+1, len-indop);
   return boolop;
}


////////////////////////////////////////////////////////////////////////////////
/// Save current attributes in a macro

void TGeoManager::SaveAttributes(const char *filename)
{
   if (!fTopNode) {
      Error("SaveAttributes","geometry must be closed first\n");
      return;
   }
   std::ofstream out;
   TString fname(filename);
   if (fname.IsNull()) fname = "tgeoatt.C";
   out.open(fname, std::ios::out);
   if (!out.good()) {
      Error("SaveAttributes", "cannot open file");
      return;
   }
   // write header
   TDatime t;
   TString sname(fname);
   sname.ReplaceAll(".C", "");
   out << sname.Data()<<"()"<<std::endl;
   out << "{" << std::endl;
   out << "//=== Macro generated by ROOT version "<< gROOT->GetVersion()<<" : "<<t.AsString()<<std::endl;
   out << "//=== Attributes for " << GetTitle() << " geometry"<<std::endl;
   out << "//===== <run this macro AFTER loading the geometry in memory>"<<std::endl;
   // save current top volume
   out << "   TGeoVolume *top = gGeoManager->GetVolume(\""<<fTopVolume->GetName()<<"\");"<<std::endl;
   out << "   TGeoVolume *vol = 0;"<<std::endl;
   out << "   TGeoNode *node = 0;"<<std::endl;
   out << "   // clear all volume attributes and get painter"<<std::endl;
   out << "   gGeoManager->ClearAttributes();"<<std::endl;
   out << "   gGeoManager->GetGeomPainter();"<<std::endl;
   out << "   // set visualization modes and bomb factors"<<std::endl;
   out << "   gGeoManager->SetVisOption("<<GetVisOption()<<");"<<std::endl;
   out << "   gGeoManager->SetVisLevel("<<GetVisLevel()<<");"<<std::endl;
   out << "   gGeoManager->SetExplodedView("<<GetBombMode()<<");"<<std::endl;
   Double_t bombx, bomby, bombz, bombr;
   GetBombFactors(bombx, bomby, bombz, bombr);
   out << "   gGeoManager->SetBombFactors("<<bombx<<","<<bomby<<","<<bombz<<","<<bombr<<");"<<std::endl;
   out << "   // iterate volumes container and set new attributes"<<std::endl;
//   out << "   TIter next(gGeoManager->GetListOfVolumes());"<<std::endl;
   TGeoVolume *vol = 0;
   fTopNode->SaveAttributes(out);

   TIter next(fVolumes);
   while ((vol=(TGeoVolume*)next())) {
      vol->SetVisStreamed(kFALSE);
   }
   out << "   // draw top volume with new settings"<<std::endl;
   out << "   top->Draw();"<<std::endl;
   out << "   gPad->x3d();"<<std::endl;
   out << "}" << std::endl;
   out.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the deepest node containing fPoint, which must be set a priori.

TGeoNode *TGeoManager::SearchNode(Bool_t downwards, const TGeoNode *skipnode)
{
   return GetCurrentNavigator()->SearchNode(downwards, skipnode);
}

////////////////////////////////////////////////////////////////////////////////
/// Cross next boundary and locate within current node
/// The current point must be on the boundary of fCurrentNode.

TGeoNode *TGeoManager::CrossBoundaryAndLocate(Bool_t downwards, TGeoNode *skipnode)
{
   return GetCurrentNavigator()->CrossBoundaryAndLocate(downwards, skipnode);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance to next boundary within STEPMAX. If no boundary is found,
/// propagate current point along current direction with fStep=STEPMAX. Otherwise
/// propagate with fStep=SNEXT (distance to boundary) and locate/return the next
/// node.

TGeoNode *TGeoManager::FindNextBoundaryAndStep(Double_t stepmax, Bool_t compsafe)
{
   return GetCurrentNavigator()->FindNextBoundaryAndStep(stepmax, compsafe);
}

////////////////////////////////////////////////////////////////////////////////
/// Find distance to next boundary and store it in fStep. Returns node to which this
/// boundary belongs. If PATH is specified, compute only distance to the node to which
/// PATH points. If STEPMAX is specified, compute distance only in case fSafety is smaller
/// than this value. STEPMAX represent the step to be made imposed by other reasons than
/// geometry (usually physics processes). Therefore in this case this method provides the
/// answer to the question : "Is STEPMAX a safe step ?" returning a NULL node and filling
/// fStep with a big number.
/// In case frombdr=kTRUE, the isotropic safety is set to zero.
///
/// Note : safety distance for the current point is computed ONLY in case STEPMAX is
///        specified, otherwise users have to call explicitly TGeoManager::Safety() if
///        they want this computed for the current point.

TGeoNode *TGeoManager::FindNextBoundary(Double_t stepmax, const char *path, Bool_t frombdr)
{
   // convert current point and direction to local reference
   return GetCurrentNavigator()->FindNextBoundary(stepmax,path, frombdr);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes as fStep the distance to next daughter of the current volume.
/// The point and direction must be converted in the coordinate system of the current volume.
/// The proposed step limit is fStep.

TGeoNode *TGeoManager::FindNextDaughterBoundary(Double_t *point, Double_t *dir, Int_t &idaughter, Bool_t compmatrix)
{
   return GetCurrentNavigator()->FindNextDaughterBoundary(point, dir, idaughter, compmatrix);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset current state flags.

void TGeoManager::ResetState()
{
   GetCurrentNavigator()->ResetState();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns deepest node containing current point.

TGeoNode *TGeoManager::FindNode(Bool_t safe_start)
{
   return GetCurrentNavigator()->FindNode(safe_start);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns deepest node containing current point.

TGeoNode *TGeoManager::FindNode(Double_t x, Double_t y, Double_t z)
{
   return GetCurrentNavigator()->FindNode(x, y, z);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes fast normal to next crossed boundary, assuming that the current point
/// is close enough to the boundary. Works only after calling FindNextBoundary.

Double_t *TGeoManager::FindNormalFast()
{
   return GetCurrentNavigator()->FindNormalFast();
}

////////////////////////////////////////////////////////////////////////////////
/// Computes normal vector to the next surface that will be or was already
/// crossed when propagating on a straight line from a given point/direction.
/// Returns the normal vector cosines in the MASTER coordinate system. The dot
/// product of the normal and the current direction is positive defined.

Double_t *TGeoManager::FindNormal(Bool_t forward)
{
   return GetCurrentNavigator()->FindNormal(forward);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if point (x,y,z) is still in the current node.

Bool_t TGeoManager::IsSameLocation(Double_t x, Double_t y, Double_t z, Bool_t change)
{
   return GetCurrentNavigator()->IsSameLocation(x,y,z,change);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a new point with given coordinates is the same as the last located one.

Bool_t TGeoManager::IsSamePoint(Double_t x, Double_t y, Double_t z) const
{
   return GetCurrentNavigator()->IsSamePoint(x,y,z);
}

////////////////////////////////////////////////////////////////////////////////
/// True if current node is in phi range

Bool_t TGeoManager::IsInPhiRange() const
{
   if (!fPhiCut) return kTRUE;
   const Double_t *origin;
   if (!GetCurrentNavigator() || !GetCurrentNavigator()->GetCurrentNode()) return kFALSE;
   origin = ((TGeoBBox*)GetCurrentNavigator()->GetCurrentVolume()->GetShape())->GetOrigin();
   Double_t point[3];
   LocalToMaster(origin, &point[0]);
   Double_t phi = TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
   if (phi<0) phi+=360.;
   if ((phi>=fPhimin) && (phi<=fPhimax)) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize current point and current direction vector (normalized)
/// in MARS. Return corresponding node.

TGeoNode *TGeoManager::InitTrack(const Double_t *point, const Double_t *dir)
{
   return GetCurrentNavigator()->InitTrack(point, dir);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize current point and current direction vector (normalized)
/// in MARS. Return corresponding node.

TGeoNode *TGeoManager::InitTrack(Double_t x, Double_t y, Double_t z, Double_t nx, Double_t ny, Double_t nz)
{
   return GetCurrentNavigator()->InitTrack(x,y,z,nx,ny,nz);
}

////////////////////////////////////////////////////////////////////////////////
/// Inspects path and all flags for the current state.

void TGeoManager::InspectState() const
{
   GetCurrentNavigator()->InspectState();
}

////////////////////////////////////////////////////////////////////////////////
/// Get path to the current node in the form /node0/node1/...

const char *TGeoManager::GetPath() const
{
   return GetCurrentNavigator()->GetPath();
}

////////////////////////////////////////////////////////////////////////////////
/// Get total size of geometry in bytes.

Int_t TGeoManager::GetByteCount(Option_t * /*option*/)
{
   Int_t count = 0;
   TIter next(fVolumes);
   TGeoVolume *vol;
   while ((vol=(TGeoVolume*)next())) count += vol->GetByteCount();
   TIter next1(fMatrices);
   TGeoMatrix *matrix;
   while ((matrix=(TGeoMatrix*)next1())) count += matrix->GetByteCount();
   TIter next2(fMaterials);
   TGeoMaterial *mat;
   while ((mat=(TGeoMaterial*)next2())) count += mat->GetByteCount();
   TIter next3(fMedia);
   TGeoMedium *med;
   while ((med=(TGeoMedium*)next3())) count += med->GetByteCount();
   if (fgVerboseLevel>0) Info("GetByteCount","Total size of logical tree : %i bytes", count);
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a default painter if none present. Returns pointer to it.

TVirtualGeoPainter *TGeoManager::GetGeomPainter()
{
   if (!fPainter) {
      const char *kind = "root";
      if (gROOT->IsWebDisplay() && !gROOT->IsWebDisplayBatch()) kind = "web";
      if (auto h = gROOT->GetPluginManager()->FindHandler("TVirtualGeoPainter", kind)) {
         if (h->LoadPlugin() == -1) {
            Error("GetGeomPainter", "could not load plugin for %s geo_painter", kind);
            return nullptr;
         }
         fPainter = (TVirtualGeoPainter*)h->ExecPlugin(1,this);
         if (!fPainter) {
            Error("GetGeomPainter", "could not create %s geo_painter", kind);
            return nullptr;
         }
      } else {
         Error("GetGeomPainter", "not found plugin %s for geo_painter", kind);
      }
   }
   return fPainter;
}

////////////////////////////////////////////////////////////////////////////////
/// Search for a named volume. All trailing blanks stripped.

TGeoVolume *TGeoManager::GetVolume(const char *name) const
{
   TString sname = name;
   sname = sname.Strip();
   TGeoVolume *vol = (TGeoVolume*)fVolumes->FindObject(sname.Data());
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Fast search for a named volume. All trailing blanks stripped.

TGeoVolume *TGeoManager::FindVolumeFast(const char *name, Bool_t multi)
{
   if (!fHashVolumes) {
      Int_t nvol = fVolumes->GetEntriesFast();
      Int_t ngvol = fGVolumes->GetEntriesFast();
      fHashVolumes = new THashList(nvol+1);
      fHashGVolumes = new THashList(ngvol+1);
      Int_t i;
      for (i=0; i<ngvol; i++) fHashGVolumes->AddLast(fGVolumes->At(i));
      for (i=0; i<nvol; i++) fHashVolumes->AddLast(fVolumes->At(i));
   }
   TString sname = name;
   sname = sname.Strip();
   THashList *list = fHashVolumes;
   if (multi) list = fHashGVolumes;
   TGeoVolume *vol = (TGeoVolume*)list->FindObject(sname.Data());
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve unique id for a volume name. Return -1 if name not found.

Int_t TGeoManager::GetUID(const char *volname) const
{
   TGeoManager *geom = (TGeoManager*)this;
   TGeoVolume *vol = geom->FindVolumeFast(volname, kFALSE);
   if (!vol) vol = geom->FindVolumeFast(volname, kTRUE);
   if (!vol) return -1;
   return vol->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Find if a given material duplicates an existing one.

TGeoMaterial *TGeoManager::FindDuplicateMaterial(const TGeoMaterial *mat) const
{
   Int_t index = fMaterials->IndexOf(mat);
   if (index <= 0) return 0;
   TGeoMaterial *other;
   for (Int_t i=0; i<index; i++) {
      other = (TGeoMaterial*)fMaterials->At(i);
      if (other == mat) continue;
      if (other->IsEq(mat)) return other;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Search for a named material. All trailing blanks stripped.

TGeoMaterial *TGeoManager::GetMaterial(const char *matname) const
{
   TString sname = matname;
   sname = sname.Strip();
   TGeoMaterial *mat = (TGeoMaterial*)fMaterials->FindObject(sname.Data());
   return mat;
}

////////////////////////////////////////////////////////////////////////////////
/// Search for a named tracking medium. All trailing blanks stripped.

TGeoMedium *TGeoManager::GetMedium(const char *medium) const
{
   TString sname = medium;
   sname = sname.Strip();
   TGeoMedium *med = (TGeoMedium*)fMedia->FindObject(sname.Data());
   return med;
}

////////////////////////////////////////////////////////////////////////////////
/// Search for a tracking medium with a given ID.

TGeoMedium *TGeoManager::GetMedium(Int_t numed) const
{
   TIter next(fMedia);
   TGeoMedium *med;
   while ((med=(TGeoMedium*)next())) {
      if (med->GetId()==numed) return med;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return material at position id.

TGeoMaterial *TGeoManager::GetMaterial(Int_t id) const
{
   if (id<0 || id >= fMaterials->GetSize()) return 0;
   TGeoMaterial *mat = (TGeoMaterial*)fMaterials->At(id);
   return mat;
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of named material.

Int_t TGeoManager::GetMaterialIndex(const char *matname) const
{
   TIter next(fMaterials);
   TGeoMaterial *mat;
   Int_t id = 0;
   TString sname = matname;
   sname = sname.Strip();
   while ((mat = (TGeoMaterial*)next())) {
      if (!strcmp(mat->GetName(),sname.Data()))
         return id;
      id++;
   }
   return -1;  // fail
}

////////////////////////////////////////////////////////////////////////////////
/// Randomly shoot nrays and plot intersections with surfaces for current
/// top node.

void TGeoManager::RandomRays(Int_t nrays, Double_t startx, Double_t starty, Double_t startz, const char *target_vol, Bool_t check_norm)
{
   GetGeomPainter()->RandomRays(nrays, startx, starty, startz, target_vol, check_norm);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove material at given index.

void TGeoManager::RemoveMaterial(Int_t index)
{
   TObject *obj = fMaterials->At(index);
   if (obj) fMaterials->Remove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets all pointers TGeoVolume::fField to NULL. User data becomes decoupled
/// from geometry. Deletion has to be managed by users.

void TGeoManager::ResetUserData()
{
   TIter next(fVolumes);
   TGeoVolume *vol;
   while ((vol=(TGeoVolume*)next())) vol->SetField(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Change raytracing mode.

void TGeoManager::SetRTmode(Int_t mode)
{
   fRaytraceMode = mode;
   if (fPainter && fPainter->IsRaytracing()) ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Restore the master volume of the geometry.

void TGeoManager::RestoreMasterVolume()
{
   if (fTopVolume == fMasterVolume) return;
   if (fMasterVolume) SetTopVolume(fMasterVolume);
}

////////////////////////////////////////////////////////////////////////////////
/// Voxelize all non-divided volumes.

void TGeoManager::Voxelize(Option_t *option)
{
   TGeoVolume *vol;
//   TGeoVoxelFinder *vox = 0;
   if (!fStreamVoxels && fgVerboseLevel>0) Info("Voxelize","Voxelizing...");
//   Int_t nentries = fVolumes->GetSize();
   TIter next(fVolumes);
   while ((vol = (TGeoVolume*)next())) {
      if (!fIsGeomReading) vol->SortNodes();
      if (!fStreamVoxels) {
         vol->Voxelize(option);
      }
      if (!fIsGeomReading) vol->FindOverlaps();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Send "Modified" signal to painter.

void TGeoManager::ModifiedPad() const
{
   if (!fPainter) return;
   fPainter->ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Make an TGeoArb8 volume.

TGeoVolume *TGeoManager::MakeArb8(const char *name, TGeoMedium *medium,
                                  Double_t dz, Double_t *vertices)
{
   return TGeoBuilder::Instance(this)->MakeArb8(name, medium, dz, vertices);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a box shape with given medium.

TGeoVolume *TGeoManager::MakeBox(const char *name, TGeoMedium *medium,
                                    Double_t dx, Double_t dy, Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeBox(name, medium, dx, dy, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a parallelepiped shape with given medium.

TGeoVolume *TGeoManager::MakePara(const char *name, TGeoMedium *medium,
                                    Double_t dx, Double_t dy, Double_t dz,
                                    Double_t alpha, Double_t theta, Double_t phi)
{
   return TGeoBuilder::Instance(this)->MakePara(name, medium, dx, dy, dz, alpha, theta, phi);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a sphere shape with given medium

TGeoVolume *TGeoManager::MakeSphere(const char *name, TGeoMedium *medium,
                                    Double_t rmin, Double_t rmax, Double_t themin, Double_t themax,
                                    Double_t phimin, Double_t phimax)
{
   return TGeoBuilder::Instance(this)->MakeSphere(name, medium, rmin, rmax, themin, themax, phimin, phimax);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a torus shape with given medium.

TGeoVolume *TGeoManager::MakeTorus(const char *name, TGeoMedium *medium, Double_t r,
                                   Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
{
   return TGeoBuilder::Instance(this)->MakeTorus(name, medium, r, rmin, rmax, phi1, dphi);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium.

TGeoVolume *TGeoManager::MakeTube(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax, Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeTube(name, medium, rmin, rmax, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube segment shape with given medium.
/// The segment will be from phiStart to phiEnd, the angles are expressed in degree

TGeoVolume *TGeoManager::MakeTubs(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax, Double_t dz,
                                     Double_t phiStart, Double_t phiEnd)
{
   return TGeoBuilder::Instance(this)->MakeTubs(name, medium, rmin, rmax, dz, phiStart, phiEnd);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium

TGeoVolume *TGeoManager::MakeEltu(const char *name, TGeoMedium *medium,
                                     Double_t a, Double_t b, Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeEltu(name, medium, a, b, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium

TGeoVolume *TGeoManager::MakeHype(const char *name, TGeoMedium *medium,
                                        Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeHype(name, medium, rin, stin, rout, stout, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium

TGeoVolume *TGeoManager::MakeParaboloid(const char *name, TGeoMedium *medium,
                                        Double_t rlo, Double_t rhi, Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeParaboloid(name, medium, rlo, rhi, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube segment shape with given medium

TGeoVolume *TGeoManager::MakeCtub(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                                     Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
{
   return TGeoBuilder::Instance(this)->MakeCtub(name, medium, rmin, rmax, dz, phi1, phi2, lx, ly, lz, tx, ty, tz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a cone shape with given medium.

TGeoVolume *TGeoManager::MakeCone(const char *name, TGeoMedium *medium,
                                     Double_t dz, Double_t rmin1, Double_t rmax1,
                                     Double_t rmin2, Double_t rmax2)
{
   return TGeoBuilder::Instance(this)->MakeCone(name, medium, dz, rmin1, rmax1, rmin2, rmax2);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a cone segment shape with given medium

TGeoVolume *TGeoManager::MakeCons(const char *name, TGeoMedium *medium,
                                     Double_t dz, Double_t rmin1, Double_t rmax1,
                                     Double_t rmin2, Double_t rmax2,
                                     Double_t phi1, Double_t phi2)
{
   return TGeoBuilder::Instance(this)->MakeCons(name, medium, dz, rmin1, rmax1, rmin2, rmax2, phi1, phi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a polycone shape with given medium.

TGeoVolume *TGeoManager::MakePcon(const char *name, TGeoMedium *medium,
                                     Double_t phi, Double_t dphi, Int_t nz)
{
   return TGeoBuilder::Instance(this)->MakePcon(name, medium, phi, dphi, nz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a polygone shape with given medium.

TGeoVolume *TGeoManager::MakePgon(const char *name, TGeoMedium *medium,
                                     Double_t phi, Double_t dphi, Int_t nedges, Int_t nz)
{
   return TGeoBuilder::Instance(this)->MakePgon(name, medium, phi, dphi, nedges, nz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a TGeoTrd1 shape with given medium.

TGeoVolume *TGeoManager::MakeTrd1(const char *name, TGeoMedium *medium,
                                  Double_t dx1, Double_t dx2, Double_t dy, Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeTrd1(name, medium, dx1, dx2, dy, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a TGeoTrd2 shape with given medium.

TGeoVolume *TGeoManager::MakeTrd2(const char *name, TGeoMedium *medium,
                                  Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2,
                                  Double_t dz)
{
   return TGeoBuilder::Instance(this)->MakeTrd2(name, medium, dx1, dx2, dy1, dy2, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a trapezoid shape with given medium.

TGeoVolume *TGeoManager::MakeTrap(const char *name, TGeoMedium *medium,
                                  Double_t dz, Double_t theta, Double_t phi, Double_t h1,
                                  Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                  Double_t tl2, Double_t alpha2)
{
   return TGeoBuilder::Instance(this)->MakeTrap(name, medium, dz, theta, phi, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2);
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a twisted trapezoid shape with given medium.

TGeoVolume *TGeoManager::MakeGtra(const char *name, TGeoMedium *medium,
                                  Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
                                  Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                  Double_t tl2, Double_t alpha2)
{
   return TGeoBuilder::Instance(this)->MakeGtra(name, medium, dz, theta, phi, twist, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a TGeoXtru-shaped volume with nz planes

TGeoVolume *TGeoManager::MakeXtru(const char *name, TGeoMedium *medium, Int_t nz)
{
   return TGeoBuilder::Instance(this)->MakeXtru(name, medium, nz);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates an alignable object with unique name corresponding to a path
/// and adds it to the list of alignables. An optional unique ID can be
/// provided, in which case PN entries can be searched fast by uid.

TGeoPNEntry *TGeoManager::SetAlignableEntry(const char *unique_name, const char *path,
                                            Int_t uid)
{
   if (!CheckPath(path)) return NULL;
   if (!fHashPNE) fHashPNE = new THashList(256,3);
   if (!fArrayPNE) fArrayPNE = new TObjArray(256);
   TGeoPNEntry *entry = GetAlignableEntry(unique_name);
   if (entry) {
      Error("SetAlignableEntry", "An alignable object with name %s already existing. NOT ADDED !", unique_name);
      return 0;
   }
   entry = new TGeoPNEntry(unique_name, path);
   Int_t ientry = fHashPNE->GetSize();
   fHashPNE->Add(entry);
   fArrayPNE->AddAtAndExpand(entry, ientry);
   if (uid>=0) {
      Bool_t added = InsertPNEId(uid, ientry);
      if (!added) Error("SetAlignableEntry", "A PN entry: has already uid=%i", uid);
   }
   return entry;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieves an existing alignable object.

TGeoPNEntry *TGeoManager::GetAlignableEntry(const char *name) const
{
   if (!fHashPNE) return 0;
   return (TGeoPNEntry*)fHashPNE->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieves an existing alignable object at a given index.

TGeoPNEntry *TGeoManager::GetAlignableEntry(Int_t index) const
{
   if (!fArrayPNE && !InitArrayPNE()) return 0;
   return (TGeoPNEntry*)fArrayPNE->At(index);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieves an existing alignable object having a preset UID.

TGeoPNEntry *TGeoManager::GetAlignableEntryByUID(Int_t uid) const
{
   if (!fNPNEId || (!fArrayPNE && !InitArrayPNE())) return NULL;
   Int_t index = TMath::BinarySearch(fNPNEId, fKeyPNEId, uid);
   if (index<0 || fKeyPNEId[index]!=uid) return NULL;
   return (TGeoPNEntry*)fArrayPNE->At(fValuePNEId[index]);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieves number of PN entries with or without UID.

Int_t TGeoManager::GetNAlignable(Bool_t with_uid) const
{
   if (!fHashPNE) return 0;
   if (with_uid) return fNPNEId;
   return fHashPNE->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Insert a PN entry in the sorted array of indexes.

Bool_t TGeoManager::InsertPNEId(Int_t uid, Int_t ientry)
{
   if (!fSizePNEId) {
      // Create the arrays.
      fSizePNEId = 128;
      fKeyPNEId = new Int_t[fSizePNEId];
      memset(fKeyPNEId, 0, fSizePNEId*sizeof(Int_t));
      fValuePNEId = new Int_t[fSizePNEId];
      memset(fValuePNEId, 0, fSizePNEId*sizeof(Int_t));
      fKeyPNEId[fNPNEId] = uid;
      fValuePNEId[fNPNEId++] = ientry;
      return kTRUE;
   }
   // Search id in the existing array and return false if it already exists.
   Int_t index = TMath::BinarySearch(fNPNEId, fKeyPNEId, uid);
   if (index>0 && fKeyPNEId[index]==uid) return kFALSE;
   // Resize the arrays and insert the value
   Bool_t resize = (fNPNEId==fSizePNEId)?kTRUE:kFALSE;
   if (resize) {
      // Double the size of the array
      fSizePNEId *= 2;
      // Create new arrays of keys and values
      Int_t *keys = new Int_t[fSizePNEId];
      memset(keys, 0, fSizePNEId*sizeof(Int_t));
      Int_t *values = new Int_t[fSizePNEId];
      memset(values, 0, fSizePNEId*sizeof(Int_t));
      // Copy all keys<uid in the new keys array (0 to index)
      memcpy(keys,   fKeyPNEId,   (index+1)*sizeof(Int_t));
      memcpy(values, fValuePNEId, (index+1)*sizeof(Int_t));
      // Insert current key at index+1
      keys[index+1]   = uid;
      values[index+1] = ientry;
      // Copy all remaining keys from the old to new array
      memcpy(&keys[index+2],   &fKeyPNEId[index+1],   (fNPNEId-index-1)*sizeof(Int_t));
      memcpy(&values[index+2], &fValuePNEId[index+1], (fNPNEId-index-1)*sizeof(Int_t));
      delete [] fKeyPNEId;
      fKeyPNEId = keys;
      delete [] fValuePNEId;
      fValuePNEId = values;
      fNPNEId++;
      return kTRUE;
   }
   // Insert the value in the existing arrays
   Int_t i;
   for (i=fNPNEId-1; i>index; i--) {
      fKeyPNEId[i+1] = fKeyPNEId[i];
      fValuePNEId[i+1] = fValuePNEId[i];
   }
   fKeyPNEId[index+1] = uid;
   fValuePNEId[index+1] = ientry;
   fNPNEId++;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a physical node from the path pointed by an alignable object with a given name.

TGeoPhysicalNode *TGeoManager::MakeAlignablePN(const char *name)
{
   TGeoPNEntry *entry = GetAlignableEntry(name);
   if (!entry) {
      Error("MakeAlignablePN","No alignable object named %s found !", name);
      return 0;
   }
   return MakeAlignablePN(entry);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a physical node from the path pointed by a given alignable object.

TGeoPhysicalNode *TGeoManager::MakeAlignablePN(TGeoPNEntry *entry)
{
   if (!entry) {
      Error("MakeAlignablePN","No alignable object specified !");
      return 0;
   }
   const char *path = entry->GetTitle();
   if (!cd(path)) {
      Error("MakeAlignablePN", "Alignable object %s poins to invalid path: %s",
            entry->GetName(), path);
      return 0;
   }
   TGeoPhysicalNode *node = MakePhysicalNode(path);
   entry->SetPhysicalNode(node);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Makes a physical node corresponding to a path. If PATH is not specified,
/// makes physical node matching current modeller state.

TGeoPhysicalNode *TGeoManager::MakePhysicalNode(const char *path)
{
   TGeoPhysicalNode *node;
   if (path) {
      if (!CheckPath(path)) {
         Error("MakePhysicalNode", "path: %s not valid", path);
         return NULL;
      }
      node = new TGeoPhysicalNode(path);
   } else {
      node = new TGeoPhysicalNode(GetPath());
   }
   fPhysicalNodes->Add(node);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh physical nodes to reflect the actual geometry paths after alignment
/// was applied. Optionally locks physical nodes (default).

void TGeoManager::RefreshPhysicalNodes(Bool_t lock)
{
   TIter next(gGeoManager->GetListOfPhysicalNodes());
   TGeoPhysicalNode *pn;
   while ((pn=(TGeoPhysicalNode*)next())) pn->Refresh();
   if (fParallelWorld && fParallelWorld->IsClosed()) fParallelWorld->RefreshPhysicalNodes();
   if (lock) LockGeometry();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the current list of physical nodes, so that we can start over with a new list.
/// If MUSTDELETE is true, delete previous nodes.

void TGeoManager::ClearPhysicalNodes(Bool_t mustdelete)
{
   if (mustdelete) fPhysicalNodes->Delete();
   else fPhysicalNodes->Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Make an assembly of volumes.

TGeoVolumeAssembly *TGeoManager::MakeVolumeAssembly(const char *name)
{
   return TGeoBuilder::Instance(this)->MakeVolumeAssembly(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a TGeoVolumeMulti handling a list of volumes.

TGeoVolumeMulti *TGeoManager::MakeVolumeMulti(const char *name, TGeoMedium *medium)
{
   return TGeoBuilder::Instance(this)->MakeVolumeMulti(name, medium);
}

////////////////////////////////////////////////////////////////////////////////
/// Set type of exploding view (see TGeoPainter::SetExplodedView())

void TGeoManager::SetExplodedView(Int_t ibomb)
{
   if ((ibomb>=0) && (ibomb<4)) fExplodedView = ibomb;
   if (fPainter) fPainter->SetExplodedView(ibomb);
}

////////////////////////////////////////////////////////////////////////////////
/// Set cut phi range

void TGeoManager::SetPhiRange(Double_t phimin, Double_t phimax)
{
   if ((phimin==0) && (phimax==360)) {
      fPhiCut = kFALSE;
      return;
   }
   fPhiCut = kTRUE;
   fPhimin = phimin;
   fPhimax = phimax;
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of segments for approximating circles in drawing.

void TGeoManager::SetNsegments(Int_t nseg)
{
   if (fNsegments==nseg) return;
   if (nseg>2) fNsegments = nseg;
   if (fPainter) fPainter->SetNsegments(nseg);
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of segments approximating circles

Int_t TGeoManager::GetNsegments() const
{
   return fNsegments;
}

////////////////////////////////////////////////////////////////////////////////
/// Now just a shortcut for GetElementTable.

void TGeoManager::BuildDefaultMaterials()
{
   GetElementTable();
   TGeoVolume::CreateDummyMedium();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns material table. Creates it if not existing.

TGeoElementTable *TGeoManager::GetElementTable()
{
   if (!fElementTable) fElementTable = new TGeoElementTable(200);
   return fElementTable;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a rectilinear step of length fStep from current point (fPoint) on current
/// direction (fDirection). If the step is imposed by geometry, is_geom flag
/// must be true (default). The cross flag specifies if the boundary should be
/// crossed in case of a geometry step (default true). Returns new node after step.
/// Set also on boundary condition.

TGeoNode *TGeoManager::Step(Bool_t is_geom, Bool_t cross)
{
   return GetCurrentNavigator()->Step(is_geom, cross);
}

////////////////////////////////////////////////////////////////////////////////
/// shoot npoints randomly in a box of 1E-5 around current point.
/// return minimum distance to points outside

TGeoNode *TGeoManager::SamplePoints(Int_t npoints, Double_t &dist, Double_t epsil,
                                    const char* g3path)
{
   return GetGeomPainter()->SamplePoints(npoints, dist, epsil, g3path);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the top volume and corresponding node as starting point of the geometry.

void TGeoManager::SetTopVolume(TGeoVolume *vol)
{
   if (fTopVolume==vol) return;

   TSeqCollection *brlist = gROOT->GetListOfBrowsers();
   TIter next(brlist);
   TBrowser *browser = 0;

   if (fTopVolume) fTopVolume->SetTitle("");
   fTopVolume = vol;
   vol->SetTitle("Top volume");
   if (fTopNode) {
      TGeoNode *topn = fTopNode;
      fTopNode = 0;
      while ((browser=(TBrowser*)next())) browser->RecursiveRemove(topn);
      delete topn;
   } else {
      fMasterVolume = vol;
      fMasterVolume->Grab();
      fUniqueVolumes->AddAtAndExpand(vol,0);
      if (fgVerboseLevel>0) Info("SetTopVolume","Top volume is %s. Master volume is %s", fTopVolume->GetName(),
           fMasterVolume->GetName());
   }
//   fMasterVolume->FindMatrixOfDaughterVolume(vol);
//   fCurrentMatrix->Print();
   fTopNode = new TGeoNodeMatrix(vol, gGeoIdentity);
   fTopNode->SetName(TString::Format("%s_1",vol->GetName()));
   fTopNode->SetNumber(1);
   fTopNode->SetTitle("Top logical node");
   fNodes->AddAt(fTopNode, 0);
   if (!GetCurrentNavigator()) {
      fCurrentNavigator = AddNavigator();
      return;
   }
   Int_t nnavigators = 0;
   TGeoNavigatorArray *arr = GetListOfNavigators();
   if (!arr) return;
   nnavigators = arr->GetEntriesFast();
   for (Int_t i=0; i<nnavigators; i++) {
      TGeoNavigator *nav = (TGeoNavigator*)arr->At(i);
      nav->ResetAll();
      if (fClosed) nav->GetCache()->BuildInfoBranch();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Define different tracking media.

void TGeoManager::SelectTrackingMedia()
{
/*
   Int_t nmat = fMaterials->GetSize();
   if (!nmat) {printf(" No materials !\n"); return;}
   Int_t *media = new Int_t[nmat];
   memset(media, 0, nmat*sizeof(Int_t));
   Int_t imedia = 1;
   TGeoMaterial *mat, *matref;
   mat = (TGeoMaterial*)fMaterials->At(0);
   if (mat->GetMedia()) {
      for (Int_t i=0; i<nmat; i++) {
         mat = (TGeoMaterial*)fMaterials->At(i);
         mat->Print();
      }
      return;
   }
   mat->SetMedia(imedia);
   media[0] = imedia++;
   mat->Print();
   for (Int_t i=0; i<nmat; i++) {
      mat = (TGeoMaterial*)fMaterials->At(i);
      for (Int_t j=0; j<i; j++) {
         matref = (TGeoMaterial*)fMaterials->At(j);
         if (mat->IsEq(matref)) {
            mat->SetMedia(media[j]);
            break;
         }
         if (j==(i-1)) {
         // different material
            mat->SetMedia(imedia);
            media[i] = imedia++;
            mat->Print();
         }
      }
   }
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Check pushes and pulls needed to cross the next boundary with respect to the
/// position given by FindNextBoundary. If radius is not mentioned the full bounding
/// box will be sampled.

void TGeoManager::CheckBoundaryErrors(Int_t ntracks, Double_t radius)
{
   GetGeomPainter()->CheckBoundaryErrors(ntracks, radius);
}

////////////////////////////////////////////////////////////////////////////////
/// Check the boundary errors reference file created by CheckBoundaryErrors method.
/// The shape for which the crossing failed is drawn with the starting point in red
/// and the extrapolated point to boundary (+/- failing push/pull) in yellow.

void TGeoManager::CheckBoundaryReference(Int_t icheck)
{
   GetGeomPainter()->CheckBoundaryReference(icheck);
}

////////////////////////////////////////////////////////////////////////////////
/// Classify a given point. See TGeoChecker::CheckPoint().

void TGeoManager::CheckPoint(Double_t x, Double_t y, Double_t z, Option_t *option)
{
   GetGeomPainter()->CheckPoint(x,y,z,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Test for shape navigation methods. Summary for test numbers:
///  - 1: DistFromInside/Outside. Sample points inside the shape. Generate
///    directions randomly in cos(theta). Compute DistFromInside and move the
///    point with bigger distance. Compute DistFromOutside back from new point.
///    Plot d-(d1+d2)
///

void TGeoManager::CheckShape(TGeoShape *shape, Int_t testNo, Int_t nsamples, Option_t *option)
{
   GetGeomPainter()->CheckShape(shape, testNo, nsamples, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Geometry checking.
/// - if option contains 'o': Optional overlap checkings (by sampling and by mesh).
/// - if option contains 'b': Optional boundary crossing check + timing per volume.
///
/// STAGE 1: extensive overlap checking by sampling per volume. Stdout need to be
///  checked by user to get report, then TGeoVolume::CheckOverlaps(0.01, "s") can
///  be called for the suspicious volumes.
///
/// STAGE 2: normal overlap checking using the shapes mesh - fills the list of
///  overlaps.
///
/// STAGE 3: shooting NRAYS rays from VERTEX and counting the total number of
///  crossings per volume (rays propagated from boundary to boundary until
///  geometry exit). Timing computed and results stored in a histo.
///
/// STAGE 4: shooting 1 mil. random rays inside EACH volume and calling
///  FindNextBoundary() + Safety() for each call. The timing is normalized by the
///  number of crossings computed at stage 2 and presented as percentage.
///  One can get a picture on which are the most "burned" volumes during
///  transportation from geometry point of view. Another plot of the timing per
///  volume vs. number of daughters is produced.

void TGeoManager::CheckGeometryFull(Int_t ntracks, Double_t vx, Double_t vy, Double_t vz, Option_t *option)
{
   TString opt(option);
   opt.ToLower();
   if (!opt.Length()) {
      Error("CheckGeometryFull","The option string must contain a letter. See method documentation.");
      return;
   }
   Bool_t checkoverlaps  = opt.Contains("o");
   Bool_t checkcrossings = opt.Contains("b");
   Double_t vertex[3];
   vertex[0] = vx;
   vertex[1] = vy;
   vertex[2] = vz;
   GetGeomPainter()->CheckGeometryFull(checkoverlaps,checkcrossings,ntracks,vertex);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform last checks on the geometry

void TGeoManager::CheckGeometry(Option_t * /*option*/)
{
   if (fgVerboseLevel>0) Info("CheckGeometry","Fixing runtime shapes...");
   TIter next(fShapes);
   TIter nextv(fVolumes);
   TGeoShape *shape;
   TGeoVolume *vol;
   Bool_t has_runtime = kFALSE;
   while ((shape = (TGeoShape*)next())) {
      if (shape->IsRunTimeShape()) {
         has_runtime = kTRUE;
      }
      if (fIsGeomReading) shape->AfterStreamer();
      if (shape->TestShapeBit(TGeoShape::kGeoPcon) || shape->TestShapeBit(TGeoShape::kGeoArb8))
         if (!shape->TestShapeBit(TGeoShape::kGeoClosedShape)) shape->ComputeBBox();
   }
   if (has_runtime) fTopNode->CheckShapes();
   else if (fgVerboseLevel>0) Info("CheckGeometry","...Nothing to fix");
   // Compute bounding  box for assemblies
   TGeoMedium *dummy = TGeoVolume::DummyMedium();
   while ((vol = (TGeoVolume*)nextv())) {
      if (vol->IsAssembly()) vol->GetShape()->ComputeBBox();
      else if (vol->GetMedium() == dummy) {
         Warning("CheckGeometry", "Volume \"%s\" has no medium: assigned dummy medium and material", vol->GetName());
         vol->SetMedium(dummy);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check all geometry for illegal overlaps within a limit OVLP.

void TGeoManager::CheckOverlaps(Double_t ovlp, Option_t * option)
{
   if (!fTopNode) {
      Error("CheckOverlaps","Top node not set");
      return;
   }
   fTopNode->CheckOverlaps(ovlp,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Prints the current list of overlaps.

void TGeoManager::PrintOverlaps() const
{
   if (!fOverlaps) return;
   Int_t novlp = fOverlaps->GetEntriesFast();
   if (!novlp) return;
   TGeoManager *geom = (TGeoManager*)this;
   geom->GetGeomPainter()->PrintOverlaps();
}

////////////////////////////////////////////////////////////////////////////////
/// Estimate weight of volume VOL with a precision SIGMA(W)/W better than PRECISION.
/// Option can be "v" - verbose (default)

Double_t TGeoManager::Weight(Double_t precision, Option_t *option)
{
   GetGeomPainter();
   TString opt(option);
   opt.ToLower();
   Double_t weight;
   TGeoVolume *volume = fTopVolume;
   if (opt.Contains("v")) {
      if (opt.Contains("a")) {
         if (fgVerboseLevel>0) Info("Weight", "Computing analytically weight of %s", volume->GetName());
         weight = volume->WeightA();
         if (fgVerboseLevel>0) Info("Weight", "Computed weight: %f [kg]\n", weight);
         return weight;
      }
      if (fgVerboseLevel>0) {
         Info("Weight", "Estimating weight of %s with %g %% precision", fTopVolume->GetName(), 100.*precision);
         printf("    event         weight         err\n");
         printf("========================================\n");
      }
   }
   weight = fPainter->Weight(precision, option);
   return weight;
}

////////////////////////////////////////////////////////////////////////////////
/// computes the total size in bytes of the branch starting with node.
/// The option can specify if all the branch has to be parsed or only the node

ULong_t TGeoManager::SizeOf(const TGeoNode * /*node*/, Option_t * /*option*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGeoManager.

void TGeoManager::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TGeoManager::Class(), this);
      fIsGeomReading = kTRUE;
      CloseGeometry();
      fStreamVoxels = kFALSE;
      fIsGeomReading = kFALSE;
   } else {
      R__b.WriteClassBuffer(TGeoManager::Class(), this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse actions on this manager.

void TGeoManager::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!fPainter) return;
   fPainter->ExecuteManagerEvent(this, event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Export this geometry to a file
///
///  - Case 1: root file or root/xml file
///    if filename end with ".root". The key will be named name
///    By default the geometry is saved without the voxelisation info.
///    Use option 'v" to save the voxelisation info.
///    if filename end with ".xml" a root/xml file is produced.
///
///  - Case 2: C++ script
///    if filename end with ".C"
///
///  - Case 3: gdml file
///    if filename end with ".gdml"
///    NOTE that to use this option, the PYTHONPATH must be defined like
///      export PYTHONPATH=$ROOTSYS/lib:$ROOTSYS/geom/gdml
///

Int_t TGeoManager::Export(const char *filename, const char *name, Option_t *option)
{
   TString sfile(filename);
   if (sfile.Contains(".C")) {
      //Save geometry as a C++ script
      if (fgVerboseLevel>0) Info("Export","Exporting %s %s as C++ code", GetName(), GetTitle());
      fTopVolume->SaveAs(filename);
      return 1;
   }
   if (sfile.Contains(".gdml")) {
      //Save geometry as a gdml file
      if (fgVerboseLevel>0) Info("Export","Exporting %s %s as gdml code", GetName(), GetTitle());
      //C++ version
      TString cmd ;
      cmd = TString::Format("TGDMLWrite::StartGDMLWriting(gGeoManager,\"%s\",\"%s\")", filename, option);
      gROOT->ProcessLineFast(cmd);
      return 1;
   }
   if (sfile.Contains(".root") || sfile.Contains(".xml")) {
      //Save geometry as a root file
      TFile *f = TFile::Open(filename,"recreate");
      if (!f || f->IsZombie()) {
         Error("Export","Cannot open file");
         return 0;
      }
      TString keyname = name;
      if (keyname.IsNull()) keyname = GetName();
      TString opt = option;
      opt.ToLower();
      if (opt.Contains("v")) {
         fStreamVoxels = kTRUE;
         if (fgVerboseLevel>0) Info("Export","Exporting %s %s as root file. Optimizations streamed.", GetName(), GetTitle());
      } else {
         fStreamVoxels = kFALSE;
         if (fgVerboseLevel>0) Info("Export","Exporting %s %s as root file. Optimizations not streamed.", GetName(), GetTitle());
      }

      const char *precision_dbl = TBufferText::GetDoubleFormat();
      const char *precision_flt = TBufferText::GetFloatFormat();
      TString new_format_dbl = TString::Format("%%.%dg", TGeoManager::GetExportPrecision());
      if (sfile.Contains(".xml")) {
        TBufferText::SetDoubleFormat(new_format_dbl.Data());
        TBufferText::SetFloatFormat(new_format_dbl.Data());
      }
      Int_t nbytes = Write(keyname);
      if (sfile.Contains(".xml")) {
        TBufferText::SetFloatFormat(precision_dbl);
        TBufferText::SetDoubleFormat(precision_flt);
      }

      fStreamVoxels = kFALSE;
      delete f;
      return nbytes;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Lock current geometry so that no other geometry can be imported.

void TGeoManager::LockGeometry()
{
   fgLock = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Unlock current geometry.

void TGeoManager::UnlockGeometry()
{
   fgLock = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check lock state.

Bool_t TGeoManager::IsLocked()
{
   return fgLock;
}

////////////////////////////////////////////////////////////////////////////////
/// Set verbosity level (static function).
///  - 0 - suppress messages related to geom-painter visibility level
///  - 1 - default value

Int_t TGeoManager::GetVerboseLevel()
{
   return fgVerboseLevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current verbosity level (static function).

void TGeoManager::SetVerboseLevel(Int_t vl)
{
   fgVerboseLevel = vl;
}

////////////////////////////////////////////////////////////////////////////////
///static function
///Import a geometry from a gdml or ROOT file
///
///  - Case 1: gdml
///    if filename ends with ".gdml" the foreign geometry described with gdml
///    is imported executing some python scripts in $ROOTSYS/gdml.
///    NOTE that to use this option, the PYTHONPATH must be defined like
///      export PYTHONPATH=$ROOTSYS/lib:$ROOTSYS/gdml
///
///  - Case 2: root file (.root) or root/xml file (.xml)
///    Import in memory from filename the geometry with key=name.
///    if name="" (default), the first TGeoManager object in the file is returned.
///
/// Note that this function deletes the current gGeoManager (if one)
/// before importing the new object.

TGeoManager *TGeoManager::Import(const char *filename, const char *name, Option_t * /*option*/)
{
   if (fgLock) {
      ::Warning("TGeoManager::Import", "TGeoMananager in lock mode. NOT IMPORTING new geometry");
      return NULL;
   }
   if (!filename) return 0;
   if (fgVerboseLevel>0) ::Info("TGeoManager::Import","Reading geometry from file: %s",filename);

   if (gGeoManager) delete gGeoManager;
   gGeoManager = 0;

   if (strstr(filename,".gdml")) {
      // import from a gdml file
      new TGeoManager("GDMLImport", "Geometry imported from GDML");
      TString cmd = TString::Format("TGDMLParse::StartGDML(\"%s\")", filename);
      TGeoVolume* world = (TGeoVolume*)gROOT->ProcessLineFast(cmd);

      if(world == 0) {
         ::Error("TGeoManager::Import", "Cannot open file");
      }
      else {
         gGeoManager->SetTopVolume(world);
         gGeoManager->CloseGeometry();
         gGeoManager->DefaultColors();
      }
   } else {
      // import from a root file
      TDirectory::TContext ctxt;
      // in case a web file is specified, use the cacheread option to cache
      // this file in the cache directory
      TFile *f = 0;
      if (strstr(filename,"http")) f = TFile::Open(filename,"CACHEREAD");
      else                         f = TFile::Open(filename);
      if (!f || f->IsZombie()) {
         ::Error("TGeoManager::Import", "Cannot open file");
         return 0;
      }
      if (name && strlen(name) > 0) {
         gGeoManager = (TGeoManager*)f->Get(name);
      } else {
         TIter next(f->GetListOfKeys());
         TKey *key;
         while ((key = (TKey*)next())) {
            if (strcmp(key->GetClassName(),"TGeoManager") != 0) continue;
            gGeoManager = (TGeoManager*)key->ReadObj();
            break;
         }
      }
      delete f;
   }
   if (!gGeoManager) return 0;
   if (!gROOT->GetListOfGeometries()->FindObject(gGeoManager)) gROOT->GetListOfGeometries()->Add(gGeoManager);
   if (!gROOT->GetListOfBrowsables()->FindObject(gGeoManager)) gROOT->GetListOfBrowsables()->Add(gGeoManager);
   gGeoManager->UpdateElements();
   return gGeoManager;
}

////////////////////////////////////////////////////////////////////////////////
/// Update element flags when geometry is loaded from a file.

void TGeoManager::UpdateElements()
{
   if (!fElementTable) return;
   TIter next(fMaterials);
   TGeoMaterial *mat;
   TGeoMixture *mix;
   TGeoElement *elem, *elem_table;
   Int_t i, nelem;
   while ((mat=(TGeoMaterial*)next())) {
      if (mat->IsMixture()) {
         mix = (TGeoMixture*)mat;
         nelem = mix->GetNelements();
         for (i=0; i<nelem; i++) {
            elem = mix->GetElement(i);
            if (!elem) continue;
            elem_table = fElementTable->GetElement(elem->Z());
            if (!elem_table) continue;
            if (elem != elem_table) {
               elem_table->SetDefined(elem->IsDefined());
               elem_table->SetUsed(elem->IsUsed());
            } else {
               elem_table->SetDefined();
            }
         }
      } else {
         elem = mat->GetElement();
         if (!elem) continue;
         elem_table = fElementTable->GetElement(elem->Z());
         if (!elem_table) continue;
         if (elem != elem_table) {
            elem_table->SetDefined(elem->IsDefined());
            elem_table->SetUsed(elem->IsUsed());
         } else {
            elem_table->SetUsed();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize PNE array for fast access via index and unique-id.

Bool_t TGeoManager::InitArrayPNE() const
{
   if (fHashPNE) {
     fArrayPNE = new TObjArray(fHashPNE->GetSize());
     TIter next(fHashPNE);
     TObject *obj;
     while ((obj = next())) {
       fArrayPNE->Add(obj);
     }
     return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get time cut for drawing tracks.

Bool_t TGeoManager::GetTminTmax(Double_t &tmin, Double_t &tmax) const
{
   tmin = fTmin;
   tmax = fTmax;
   return fTimeCut;
}

////////////////////////////////////////////////////////////////////////////////
/// Set time cut interval for drawing tracks. If called with no arguments, time
/// cut will be disabled.

void TGeoManager::SetTminTmax(Double_t tmin, Double_t tmax)
{
   fTmin = tmin;
   fTmax = tmax;
   if (tmin==0 && tmax==999) fTimeCut = kFALSE;
   else fTimeCut = kTRUE;
   if (fTracks && !IsAnimatingTracks()) ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Convert coordinates from master volume frame to top.

void TGeoManager::MasterToTop(const Double_t *master, Double_t *top) const
{
   GetCurrentNavigator()->MasterToLocal(master, top);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert coordinates from top volume frame to master.

void TGeoManager::TopToMaster(const Double_t *top, Double_t *master) const
{
   GetCurrentNavigator()->LocalToMaster(top, master);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a parallel world for prioritised navigation. This can be populated
/// with physical nodes and can be navigated independently using its API.
/// In case the flag SetUseParallelWorldNav is set, any navigation query in the
/// main geometry is checked against the parallel geometry, which gets priority
/// in case of overlaps with the main geometry volumes.

TGeoParallelWorld *TGeoManager::CreateParallelWorld(const char *name)
{
   fParallelWorld = new TGeoParallelWorld(name, this);
   return fParallelWorld;
}

////////////////////////////////////////////////////////////////////////////////
/// Activate/deactivate usage of parallel world navigation. Can only be done if
/// there is a parallel world. Activating navigation will automatically close
/// the parallel geometry.

void TGeoManager::SetUseParallelWorldNav(Bool_t flag)
{
   if (!fParallelWorld) {
      Error("SetUseParallelWorldNav", "No parallel world geometry defined. Use CreateParallelWorld.");
      return;
   }
   if (!flag) {
      fUsePWNav = flag;
      return;
   }
   if (!fClosed) {
      Error("SetUseParallelWorldNav", "The geometry must be closed first");
      return;
   }
   // Closing the parallel world geometry is mandatory
   if (fParallelWorld->CloseGeometry()) fUsePWNav=kTRUE;
}

Bool_t TGeoManager::LockDefaultUnits(Bool_t new_value)    {
  Bool_t val = gGeometryLocked;
  gGeometryLocked = new_value;
  return val;
}

TGeoManager::EDefaultUnits TGeoManager::GetDefaultUnits()
{
  return fgDefaultUnits;
}

void TGeoManager::SetDefaultUnits(EDefaultUnits new_value)
{
   if ( fgDefaultUnits == new_value )   {
      return;
   }
   else if ( gGeometryLocked )    {
      ::Fatal("TGeoManager","The system of units may only be changed once, \n"
	      "BEFORE any elements and materials are created! \n"
	      "Alternatively unlock the default units at own risk.");
   }
   else if ( new_value == kG4Units )   {
      ::Warning("TGeoManager","Changing system of units to Geant4 units (mm, ns, MeV).");
   }
   else if ( new_value == kRootUnits )   {
      ::Warning("TGeoManager","Changing system of units to ROOT units (cm, s, GeV).");
   }
   fgDefaultUnits = new_value;
}
