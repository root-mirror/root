// @(#)root/geom:$Id$
// Author: Andrei Gheata   25/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoManager
#define ROOT_TGeoManager

#include <mutex>
#include <thread>
#include <map>

#include "TNamed.h"
#include "TObjArray.h"
#include "TGeoNavigator.h"

// forward declarations
class TVirtualGeoTrack;
class TGeoNode;
class TGeoPhysicalNode;
class TGeoPNEntry;
class TGeoVolume;
class TGeoVolumeMulti;
class TGeoMatrix;
class TGeoHMatrix;
class TGeoMaterial;
class TGeoMedium;
class TGeoShape;
class TVirtualGeoPainter;
class THashList;
class TGeoParallelWorld;
class TGeoRegion;
class TGDMLMatrix;
class TGeoOpticalSurface;
class TGeoSkinSurface;
class TGeoBorderSurface;

class TGeoManager : public TNamed
{
public:
  enum EDefaultUnits {
     kG4Units   = 0,
     kRootUnits = 1
  };

protected:
   static std::mutex     fgMutex;           //! mutex for navigator booking in MT mode
   static Bool_t         fgLock;            //! Lock preventing a second geometry to be loaded
   static Int_t          fgVerboseLevel;    //! Verbosity level for Info messages (no IO).
   static Int_t          fgMaxLevel;        //! Maximum level in geometry
   static Int_t          fgMaxDaughters;    //! Maximum number of daughters
   static Int_t          fgMaxXtruVert;     //! Maximum number of Xtru vertices
   static UInt_t         fgExportPrecision; //! Precision to be used in ASCII exports
   static EDefaultUnits  fgDefaultUnits;    //! Default units in GDML if not explicit in some tags

   TGeoManager(const TGeoManager&) = delete;
   TGeoManager& operator=(const TGeoManager&) = delete;

private :
   Double_t              fPhimin;           //! lowest range for phi cut
   Double_t              fPhimax;           //! highest range for phi cut
   Double_t              fTmin;             //! lower time limit for tracks drawing
   Double_t              fTmax;             //! upper time limit for tracks drawing
   Int_t                 fNNodes;           // total number of physical nodes
   TString               fPath;             //! path to current node
   TString               fParticleName;     //! particles to be drawn
   Double_t              fVisDensity;       // transparency threshold by density
   Int_t                 fExplodedView;     // exploded view mode
   Int_t                 fVisOption;        // global visualization option
   Int_t                 fVisLevel;         // maximum visualization depth
   Int_t                 fNsegments;        // number of segments to approximate circles
   Int_t                 fNtracks;          // number of tracks
   Int_t                 fMaxVisNodes;      // maximum number of visible nodes
   TVirtualGeoTrack     *fCurrentTrack;     //! current track
   Int_t                 fNpdg;             // number of different pdg's stored
   Int_t                 fPdgId[1024];      // pdg conversion table
   Bool_t                fClosed;           //! flag that geometry is closed
   Bool_t                fLoopVolumes;      //! flag volume lists loop
   Bool_t                fStreamVoxels;     // flag to allow voxelization I/O
   Bool_t                fIsGeomReading;    //! flag set when reading geometry
   Bool_t                fIsGeomCleaning;   //! flag to notify that the manager is being destructed
   Bool_t                fPhiCut;           // flag for phi cuts
   Bool_t                fTimeCut;          // time cut for tracks
   Bool_t                fDrawExtra;        //! flag that the list of physical nodes has to be drawn
   Bool_t                fMatrixTransform;  //! flag for using GL matrix
   Bool_t                fMatrixReflection; //! flag for GL reflections
   Bool_t                fActivity;         //! switch ON/OFF volume activity (default OFF - all volumes active))
   Bool_t                fIsNodeSelectable; //! flag that nodes are the selected objects in pad rather than volumes
   TVirtualGeoPainter   *fPainter;          //! current painter

   TObjArray            *fMatrices;         //-> list of local transformations
   TObjArray            *fShapes;           //-> list of shapes
   TObjArray            *fVolumes;          //-> list of volumes
   TObjArray            *fPhysicalNodes;    //-> list of physical nodes
   TObjArray            *fGShapes;          //! list of runtime shapes
   TObjArray            *fGVolumes;         //! list of runtime volumes
   TObjArray            *fTracks;           //-> list of tracks attached to geometry
   TObjArray            *fPdgNames;         //-> list of pdg names for tracks
   TObjArray            *fGDMLMatrices;     //-> list of matrices read from GDML
   TObjArray            *fOpticalSurfaces;  //-> list of optical surfaces read from GDML
   TObjArray            *fSkinSurfaces;     //-> list of skin surfaces read from GDML
   TObjArray            *fBorderSurfaces;   //-> list of border surfaces read from GDML
   TList                *fMaterials;        //-> list of materials
   TList                *fMedia;            //-> list of tracking media
   TObjArray            *fNodes;            //-> current branch of nodes
   TObjArray            *fOverlaps;         //-> list of geometrical overlaps
   TObjArray            *fRegions;          //-> list of regions
   UChar_t              *fBits;             //! bits used for voxelization
   // Map of navigator arrays per thread
   typedef std::map<std::thread::id, TGeoNavigatorArray *>   NavigatorsMap_t;
   typedef NavigatorsMap_t::iterator                         NavigatorsMapIt_t;
   typedef std::map<std::thread::id, Int_t>                  ThreadsMap_t;
   typedef ThreadsMap_t::const_iterator                      ThreadsMapIt_t;
   // Map of constant properties
   typedef std::map<std::string, Double_t>                   ConstPropMap_t;

   NavigatorsMap_t       fNavigators;       //! Map between thread id's and navigator arrays
   static ThreadsMap_t  *fgThreadId;        //! Thread id's map
   static Int_t          fgNumThreads;      //! Number of registered threads
   static Bool_t         fgLockNavigators;   //! Lock existing navigators
   TGeoNavigator        *fCurrentNavigator; //! current navigator
   TGeoVolume           *fCurrentVolume;    //! current volume
   TGeoVolume           *fTopVolume;        //! top level volume in geometry
   TGeoNode             *fTopNode;          //! top physical node
   TGeoVolume           *fMasterVolume;     // master volume
   TGeoHMatrix          *fGLMatrix;         // matrix to be used for view transformations
   TObjArray            *fUniqueVolumes;    //-> list of unique volumes
   TGeoShape            *fClippingShape;    //! clipping shape for raytracing
   TGeoElementTable     *fElementTable;     //! table of elements

   Int_t                 fNLevel;           // maximum accepted level in geometry
   TGeoVolume           *fPaintVolume;      //! volume currently painted
   TGeoVolume           *fUserPaintVolume;  //!
   THashList            *fHashVolumes;      //! hash list of volumes providing fast search
   THashList            *fHashGVolumes;     //! hash list of group volumes providing fast search
   THashList            *fHashPNE;          //-> hash list of physical node entries
   mutable TObjArray    *fArrayPNE;         //! array of physical node entries
   Int_t                 fSizePNEId;        // size of the array of unique ID's for PN entries
   Int_t                 fNPNEId;           // number of PN entries having a unique ID
   Int_t                *fKeyPNEId;         //[fSizePNEId] array of uid values for PN entries
   Int_t                *fValuePNEId;       //[fSizePNEId] array of pointers to PN entries with ID's
   Int_t                 fMaxThreads;       //! Max number of threads
   Bool_t                fMultiThread;      //! Flag for multi-threading
   Int_t                 fRaytraceMode;     //! Raytrace mode: 0=normal, 1=pass through, 2=transparent
   Bool_t                fUsePWNav;         // Activate usage of parallel world in navigation
   TGeoParallelWorld    *fParallelWorld;    // Parallel world
   ConstPropMap_t        fProperties;       // Map of user-defined constant properties
//--- private methods
   Bool_t                IsLoopingVolumes() const     {return fLoopVolumes;}
   void                  Init();
   Bool_t                InitArrayPNE() const;
   Bool_t                InsertPNEId(Int_t uid, Int_t ientry);
   void                  SetLoopVolumes(Bool_t flag=kTRUE) {fLoopVolumes=flag;}
   void                  UpdateElements();
   void                  Voxelize(Option_t *option = 0);

public:
   // constructors
   TGeoManager();
   TGeoManager(const char *name, const char *title);
   // destructor
   virtual ~TGeoManager();
   //--- adding geometrical objects
   Int_t                  AddMaterial(const TGeoMaterial *material);
   Int_t                  AddOverlap(const TNamed *ovlp);
   Int_t                  AddRegion(TGeoRegion *region);
   Int_t                  AddTransformation(const TGeoMatrix *matrix);
   Int_t                  AddShape(const TGeoShape *shape);
   Int_t                  AddTrack(Int_t id, Int_t pdgcode, TObject *particle=0);
   Int_t                  AddTrack(TVirtualGeoTrack *track);
   Int_t                  AddVolume(TGeoVolume *volume);
   TGeoNavigator         *AddNavigator();
   Bool_t                 AddProperty(const char *property, Double_t value);
   Double_t               GetProperty(const char *name, Bool_t *error = nullptr) const;
   Double_t               GetProperty(size_t i, TString &name, Bool_t *error = nullptr) const;
   Int_t                  GetNproperties() const { return fProperties.size(); }
   void                   ClearOverlaps();
   void                   RegisterMatrix(const TGeoMatrix *matrix);
   void                   SortOverlaps();
   //--- browsing and tree navigation
   void                   Browse(TBrowser *b);
   void                   SetVisibility(TObject *obj, Bool_t vis);
   virtual Bool_t         cd(const char *path=""); // *MENU*
   Bool_t                 CheckPath(const char *path) const;
   void                   CdNode(Int_t nodeid);
   void                   CdDown(Int_t index);
   void                   CdUp();
   void                   CdTop();
   void                   CdNext();
   void                   GetBranchNames(Int_t *names) const;
   void                   GetBranchNumbers(Int_t *copyNumbers, Int_t *volumeNumbers) const;
   void                   GetBranchOnlys(Int_t *isonly) const;
   Int_t                  GetNmany() const {return GetCurrentNavigator()->GetNmany();}
   const char            *GetPdgName(Int_t pdg) const;
   void                   SetPdgName(Int_t pdg, const char *name);
   Bool_t                 IsFolder() const { return kTRUE; }
   //--- visualization settings
   virtual void           Edit(Option_t *option="");  // *MENU*
   void                   BombTranslation(const Double_t *tr, Double_t *bombtr);
   void                   UnbombTranslation(const Double_t *tr, Double_t *bombtr);
   void                   ClearAttributes(); // *MENU*
   void                   DefaultAngles();   // *MENU*
   void                   DefaultColors();   // *MENU*
   TGeoShape             *GetClippingShape() const {return fClippingShape;}
   Int_t                  GetNsegments() const;
   TVirtualGeoPainter    *GetGeomPainter();
   TVirtualGeoPainter    *GetPainter() const {return fPainter;}
   Int_t                  GetBombMode() const  {return fExplodedView;}
   void                   GetBombFactors(Double_t &bombx, Double_t &bomby, Double_t &bombz, Double_t &bombr) const;
   Int_t                  GetMaxVisNodes() const {return fMaxVisNodes;}
   Bool_t                 GetTminTmax(Double_t &tmin, Double_t &tmax) const;
   Double_t               GetTmax() const {return fTmax;}
   TGeoVolume            *GetPaintVolume() const {return fPaintVolume;}
   TGeoVolume            *GetUserPaintVolume() const {return fUserPaintVolume;}
   Double_t               GetVisDensity() const  {return fVisDensity;}
   Int_t                  GetVisLevel() const;
   Int_t                  GetVisOption() const;
   Bool_t                 IsInPhiRange() const;
   Bool_t                 IsDrawingExtra() const {return fDrawExtra;}
   Bool_t                 IsNodeSelectable() const {return fIsNodeSelectable;}
   Bool_t                 IsVisLeaves() const {return (fVisOption==1)?kTRUE:kFALSE;}
   void                   ModifiedPad() const;
   void                   OptimizeVoxels(const char *filename="tgeovox.C"); // *MENU*
   void                   SetClipping(Bool_t flag=kTRUE) {SetClippingShape(((flag)?fClippingShape:0));} // *MENU*
   void                   SetClippingShape(TGeoShape *clip);
   void                   SetExplodedView(Int_t iopt=0); // *MENU*
   void                   SetPhiRange(Double_t phimin=0., Double_t phimax=360.);
   void                   SetNsegments(Int_t nseg); // *MENU*
   Bool_t                 SetCurrentNavigator(Int_t index);
   void                   SetBombFactors(Double_t bombx=1.3, Double_t bomby=1.3, Double_t bombz=1.3, Double_t bombr=1.3); // *MENU*
   void                   SetPaintVolume(TGeoVolume *vol) {fPaintVolume = vol;}
   void                   SetUserPaintVolume(TGeoVolume *vol) {fUserPaintVolume = vol;}
   void                   SetTopVisible(Bool_t vis=kTRUE);
   void                   SetTminTmax(Double_t tmin=0, Double_t tmax=999);
   void                   SetDrawExtraPaths(Bool_t flag=kTRUE) {fDrawExtra=flag;}
   void                   SetNodeSelectable(Bool_t flag=kTRUE) {fIsNodeSelectable=flag;}
   void                   SetVisDensity(Double_t dens=0.01); // *MENU*
   void                   SetVisLevel(Int_t level=3);   // *MENU*
   void                   SetVisOption(Int_t option=0);
   void                   ViewLeaves(Bool_t flag=kTRUE); // *TOGGLE* *GETTER=IsVisLeaves
   void                   SaveAttributes(const char *filename="tgeoatt.C"); // *MENU*
   void                   RestoreMasterVolume(); // *MENU*
   void                   SetMaxVisNodes(Int_t maxnodes=10000); // *MENU*
   //--- geometry checking
   void                   AnimateTracks(Double_t tmin=0, Double_t tmax=5E-8, Int_t nframes=200, Option_t *option="/*"); // *MENU*
   void                   CheckBoundaryErrors(Int_t ntracks=1000000, Double_t radius=-1.); // *MENU*
   void                   CheckBoundaryReference(Int_t icheck=-1);
   void                   CheckGeometryFull(Int_t ntracks=1000000, Double_t vx=0., Double_t vy=0., Double_t vz=0., Option_t *option="ob"); // *MENU*
   void                   CheckGeometry(Option_t *option="");
   void                   CheckOverlaps(Double_t ovlp=0.1, Option_t *option=""); // *MENU*
   void                   CheckPoint(Double_t x=0,Double_t y=0, Double_t z=0, Option_t *option=""); // *MENU*
   void                   CheckShape(TGeoShape *shape, Int_t testNo, Int_t nsamples, Option_t *option);
   void                   ConvertReflections();
   void                   DrawCurrentPoint(Int_t color=2); // *MENU*
   void                   DrawTracks(Option_t *option=""); // *MENU*
   void                   SetParticleName(const char *pname) {fParticleName=pname;}
   const char            *GetParticleName() const {return fParticleName.Data();}
   void                   DrawPath(const char *path, Option_t *option="");
   void                   PrintOverlaps() const; // *MENU*
   void                   RandomPoints(const TGeoVolume *vol, Int_t npoints=10000, Option_t *option="");
   void                   RandomRays(Int_t nrays=1000, Double_t startx=0, Double_t starty=0, Double_t startz=0, const char *target_vol=0, Bool_t check_norm=kFALSE);
   TGeoNode              *SamplePoints(Int_t npoints, Double_t &dist, Double_t epsil=1E-5,
                                       const char *g3path="");
   void                   SetNmeshPoints(Int_t npoints=1000);
   void                   SetCheckedNode(TGeoNode *node);
   void                   Test(Int_t npoints=1000000, Option_t *option=""); // *MENU*
   void                   TestOverlaps(const char* path=""); // *MENU*
   Double_t               Weight(Double_t precision=0.01, Option_t *option="va"); // *MENU*

   //--- GEANT3-like geometry creation
   TGeoVolume            *Division(const char *name, const char *mother, Int_t iaxis, Int_t ndiv,
                                         Double_t start, Double_t step, Int_t numed=0, Option_t *option="");
   void                   Matrix(Int_t index, Double_t theta1, Double_t phi1,
                                       Double_t theta2, Double_t phi2,
                                       Double_t theta3, Double_t phi3);
   TGeoMaterial          *Material(const char *name, Double_t a, Double_t z, Double_t dens, Int_t uid, Double_t radlen=0, Double_t intlen=0);
   TGeoMaterial          *Mixture(const char *name, Float_t *a, Float_t *z, Double_t dens,
                                        Int_t nelem, Float_t *wmat, Int_t uid);
   TGeoMaterial          *Mixture(const char *name, Double_t *a, Double_t *z, Double_t dens,
                                        Int_t nelem, Double_t *wmat, Int_t uid);
   TGeoMedium            *Medium(const char *name, Int_t numed, Int_t nmat, Int_t isvol,
                                       Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                                       Double_t stemax, Double_t deemax, Double_t epsil,
                                       Double_t stmin);
   void                   Node(const char *name, Int_t nr, const char *mother,
                                     Double_t x, Double_t y, Double_t z, Int_t irot,
                                     Bool_t isOnly, Float_t *upar, Int_t npar=0);
   void                   Node(const char *name, Int_t nr, const char *mother,
                                     Double_t x, Double_t y, Double_t z, Int_t irot,
                                     Bool_t isOnly, Double_t *upar, Int_t npar=0);
   TGeoVolume            *Volume(const char *name, const char *shape, Int_t nmed,
                                       Float_t *upar, Int_t npar=0);
   TGeoVolume            *Volume(const char *name, const char *shape, Int_t nmed,
                                       Double_t *upar, Int_t npar=0);
   void                   SetVolumeAttribute(const char *name, const char *att, Int_t val);
   //--- geometry building
   void                   BuildDefaultMaterials();
   void                   CloseGeometry(Option_t *option="d");
   Bool_t                 IsClosed() const {return fClosed;}
   TGeoVolume            *MakeArb8(const char *name, TGeoMedium *medium,
                                     Double_t dz, Double_t *vertices=0);
   TGeoVolume            *MakeBox(const char *name, TGeoMedium *medium,
                                     Double_t dx, Double_t dy, Double_t dz);
   TGeoVolume            *MakeCone(const char *name, TGeoMedium *medium,
                                      Double_t dz, Double_t rmin1, Double_t rmax1,
                                      Double_t rmin2, Double_t rmax2);
   TGeoVolume            *MakeCons(const char *name, TGeoMedium *medium,
                                      Double_t dz, Double_t rmin1, Double_t rmax1,
                                      Double_t rmin2, Double_t rmax2,
                                      Double_t phi1, Double_t phi2);
   TGeoVolume            *MakeCtub(const char *name, TGeoMedium *medium,
                                      Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                                      Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz);
   TGeoVolume            *MakeEltu(const char *name, TGeoMedium *medium,
                                      Double_t a, Double_t b, Double_t dz);
   TGeoVolume            *MakeGtra(const char *name, TGeoMedium *medium,
                                   Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
                                   Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                   Double_t tl2, Double_t alpha2);
   TGeoVolume            *MakePara(const char *name, TGeoMedium *medium,
                                     Double_t dx, Double_t dy, Double_t dz,
                                     Double_t alpha, Double_t theta, Double_t phi);
   TGeoVolume            *MakePcon(const char *name, TGeoMedium *medium,
                                      Double_t phi, Double_t dphi, Int_t nz);
   TGeoVolume            *MakeParaboloid(const char *name, TGeoMedium *medium,
                                      Double_t rlo, Double_t rhi, Double_t dz);
   TGeoVolume            *MakeHype(const char *name, TGeoMedium *medium,
                                      Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz);
   TGeoVolume            *MakePgon(const char *name, TGeoMedium *medium,
                                      Double_t phi, Double_t dphi, Int_t nedges, Int_t nz);
   TGeoVolume            *MakeSphere(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax,
                                     Double_t themin=0, Double_t themax=180,
                                     Double_t phimin=0, Double_t phimax=360);
   TGeoVolume            *MakeTorus(const char *name, TGeoMedium *medium, Double_t r,
                                    Double_t rmin, Double_t rmax, Double_t phi1=0, Double_t dphi=360);
   TGeoVolume            *MakeTrap(const char *name, TGeoMedium *medium,
                                   Double_t dz, Double_t theta, Double_t phi, Double_t h1,
                                   Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                   Double_t tl2, Double_t alpha2);
   TGeoVolume            *MakeTrd1(const char *name, TGeoMedium *medium,
                                      Double_t dx1, Double_t dx2, Double_t dy, Double_t dz);
   TGeoVolume            *MakeTrd2(const char *name, TGeoMedium *medium,
                                      Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2,
                                      Double_t dz);
   TGeoVolume            *MakeTube(const char *name, TGeoMedium *medium,
                                      Double_t rmin, Double_t rmax, Double_t dz);
   TGeoVolume            *MakeTubs(const char *name, TGeoMedium *medium,
                                      Double_t rmin, Double_t rmax, Double_t dz,
                                      Double_t phi1, Double_t phi2);
   TGeoVolume            *MakeXtru(const char *name, TGeoMedium *medium,
                                   Int_t nz);

   TGeoPNEntry           *SetAlignableEntry(const char *unique_name, const char *path, Int_t uid=-1);
   TGeoPNEntry           *GetAlignableEntry(const char *name) const;
   TGeoPNEntry           *GetAlignableEntry(Int_t index) const;
   TGeoPNEntry           *GetAlignableEntryByUID(Int_t uid) const;
   Int_t                  GetNAlignable(Bool_t with_uid=kFALSE) const;
   TGeoPhysicalNode      *MakeAlignablePN(const char *name);
   TGeoPhysicalNode      *MakeAlignablePN(TGeoPNEntry *entry);
   TGeoPhysicalNode      *MakePhysicalNode(const char *path=0);
   void                   ClearPhysicalNodes(Bool_t mustdelete=kFALSE);
   void                   RefreshPhysicalNodes(Bool_t lock=kTRUE);
   TVirtualGeoTrack      *MakeTrack(Int_t id, Int_t pdgcode, TObject *particle);
   TGeoVolumeAssembly    *MakeVolumeAssembly(const char *name);
   TGeoVolumeMulti       *MakeVolumeMulti(const char *name, TGeoMedium *medium);
   void                   SetTopVolume(TGeoVolume *vol);

   //--- geometry queries
   TGeoNode              *CrossBoundaryAndLocate(Bool_t downwards, TGeoNode *skipnode);
   TGeoNode              *FindNextBoundary(Double_t stepmax=TGeoShape::Big(),const char *path="", Bool_t frombdr=kFALSE);
   TGeoNode              *FindNextDaughterBoundary(Double_t *point, Double_t *dir, Int_t &idaughter, Bool_t compmatrix=kFALSE);
   TGeoNode              *FindNextBoundaryAndStep(Double_t stepmax=TGeoShape::Big(), Bool_t compsafe=kFALSE);
   TGeoNode              *FindNode(Bool_t safe_start=kTRUE);
   TGeoNode              *FindNode(Double_t x, Double_t y, Double_t z);
   Double_t              *FindNormal(Bool_t forward=kTRUE);
   Double_t              *FindNormalFast();
   TGeoNode              *InitTrack(const Double_t *point, const Double_t *dir);
   TGeoNode              *InitTrack(Double_t x, Double_t y, Double_t z, Double_t nx, Double_t ny, Double_t nz);
   void                   ResetState();
   Double_t               Safety(Bool_t inside=kFALSE);
   TGeoNode              *SearchNode(Bool_t downwards=kFALSE, const TGeoNode *skipnode=0);
   TGeoNode              *Step(Bool_t is_geom=kTRUE, Bool_t cross=kTRUE);
   void                   DisableInactiveVolumes() {fActivity=kTRUE;}
   void                   EnableInactiveVolumes()  {fActivity=kFALSE;}
   void                   SetCurrentTrack(Int_t i) {fCurrentTrack = (TVirtualGeoTrack*)fTracks->At(i);}
   void                   SetCurrentTrack(TVirtualGeoTrack *track) {fCurrentTrack=track;}
   Int_t                  GetNtracks() const {return fNtracks;}
   TVirtualGeoTrack      *GetCurrentTrack() {return fCurrentTrack;}
   TVirtualGeoTrack      *GetLastTrack() {return (TVirtualGeoTrack*)((fNtracks>0)?fTracks->At(fNtracks-1):NULL);}
   const Double_t        *GetLastPoint() const {return GetCurrentNavigator()->GetLastPoint();}
   TVirtualGeoTrack      *GetTrack(Int_t index)         {return (index<fNtracks)?(TVirtualGeoTrack*)fTracks->At(index):0;}
   Int_t                  GetTrackIndex(Int_t id) const;
   TVirtualGeoTrack      *GetTrackOfId(Int_t id) const;
   TVirtualGeoTrack      *FindTrackWithId(Int_t id) const;
   TVirtualGeoTrack      *GetParentTrackOfId(Int_t id) const;
   Int_t                  GetVirtualLevel();
   Bool_t                 GotoSafeLevel();
   Int_t                  GetSafeLevel() const;
   Double_t               GetSafeDistance() const      {return GetCurrentNavigator()->GetSafeDistance();}
   Double_t               GetLastSafety() const        {return GetCurrentNavigator()->GetLastSafety();}
   Double_t               GetStep() const              {return GetCurrentNavigator()->GetStep();}
   void                   InspectState() const;
   Bool_t                 IsAnimatingTracks() const    {return fIsGeomReading;}
   Bool_t                 IsCheckingOverlaps() const   {return GetCurrentNavigator()->IsCheckingOverlaps();}
   Bool_t                 IsMatrixTransform() const    {return fMatrixTransform;}
   Bool_t                 IsMatrixReflection() const   {return fMatrixReflection;}
   Bool_t                 IsSameLocation(Double_t x, Double_t y, Double_t z, Bool_t change=kFALSE);
   Bool_t                 IsSameLocation() const {return GetCurrentNavigator()->IsSameLocation();}
   Bool_t                 IsSamePoint(Double_t x, Double_t y, Double_t z) const;
   Bool_t                 IsStartSafe() const {return GetCurrentNavigator()->IsStartSafe();}
   void                   SetCheckingOverlaps(Bool_t flag=kTRUE) {GetCurrentNavigator()->SetCheckingOverlaps(flag);}
   void                   SetStartSafe(Bool_t flag=kTRUE)   {GetCurrentNavigator()->SetStartSafe(flag);}
   void                   SetMatrixTransform(Bool_t on=kTRUE) {fMatrixTransform = on;}
   void                   SetMatrixReflection(Bool_t flag=kTRUE) {fMatrixReflection = flag;}
   void                   SetStep(Double_t step) {GetCurrentNavigator()->SetStep(step);}
   Bool_t                 IsCurrentOverlapping() const {return GetCurrentNavigator()->IsCurrentOverlapping();}
   Bool_t                 IsEntering() const           {return GetCurrentNavigator()->IsEntering();}
   Bool_t                 IsExiting() const            {return GetCurrentNavigator()->IsExiting();}
   Bool_t                 IsStepEntering() const       {return GetCurrentNavigator()->IsStepEntering();}
   Bool_t                 IsStepExiting() const        {return GetCurrentNavigator()->IsStepExiting();}
   Bool_t                 IsOutside() const            {return GetCurrentNavigator()->IsOutside();}
   Bool_t                 IsOnBoundary() const         {return GetCurrentNavigator()->IsOnBoundary();}
   Bool_t                 IsNullStep() const           {return GetCurrentNavigator()->IsNullStep();}
   Bool_t                 IsActivityEnabled() const    {return fActivity;}
   void                   SetOutside(Bool_t flag=kTRUE) {GetCurrentNavigator()->SetOutside(flag);}


   //--- cleaning
   void                   CleanGarbage();
   void                   ClearShape(const TGeoShape *shape);
   void                   ClearTracks() {fTracks->Delete(); fNtracks=0;}
   void                   ClearNavigators();
   void                   RemoveMaterial(Int_t index);
   void                   RemoveNavigator(const TGeoNavigator *nav);
   void                   ResetUserData();


   //--- utilities
   Int_t                  CountNodes(const TGeoVolume *vol=0, Int_t nlevels=10000, Int_t option=0);
   void                   CountLevels();
   virtual void           ExecuteEvent(Int_t event, Int_t px, Int_t py);
   static Int_t           Parse(const char* expr, TString &expr1, TString &expr2, TString &expr3);
   Int_t                  ReplaceVolume(TGeoVolume *vorig, TGeoVolume *vnew);
   Int_t                  TransformVolumeToAssembly(const char *vname);
   UChar_t               *GetBits() {return fBits;}
   virtual Int_t          GetByteCount(Option_t *option=0);
   void                   SetAllIndex();
   static Int_t           GetMaxDaughters();
   static Int_t           GetMaxLevels();
   static Int_t           GetMaxXtruVert();
   Int_t                  GetMaxThreads() const {return fMaxThreads-1;}
   void                   SetMaxThreads(Int_t nthreads);
   Int_t                  GetRTmode() const {return fRaytraceMode;}
   void                   SetRTmode(Int_t mode); // *MENU*
   Bool_t                 IsMultiThread() const {return fMultiThread;}
   static void            SetNavigatorsLock(Bool_t flag);
   static Int_t           ThreadId();
   static Int_t           GetNumThreads();
   static void            ClearThreadsMap();
   void                   ClearThreadData() const;
   void                   CreateThreadData() const;

   //--- I/O
   virtual Int_t          Export(const char *filename, const char *name="", Option_t *option="vg");
   static  void           LockGeometry();
   static  void           UnlockGeometry();
   static  Int_t          GetVerboseLevel();
   static  void           SetVerboseLevel(Int_t vl);
   static TGeoManager    *Import(const char *filename, const char *name="", Option_t *option="");
   static Bool_t          IsLocked();
   static void            SetExportPrecision(UInt_t prec) {fgExportPrecision = prec;}
   static UInt_t          GetExportPrecision() {return fgExportPrecision;}
   static void            SetDefaultUnits(EDefaultUnits new_value);
   static EDefaultUnits   GetDefaultUnits();
   static Bool_t          LockDefaultUnits(Bool_t new_value);
   Bool_t                 IsStreamingVoxels() const {return fStreamVoxels;}
   Bool_t                 IsCleaning() const {return fIsGeomCleaning;}

   //--- list getters
   TObjArray             *GetListOfNodes()              {return fNodes;}
   TObjArray             *GetListOfPhysicalNodes()      {return fPhysicalNodes;}
   TObjArray             *GetListOfOverlaps()           {return fOverlaps;}
   TObjArray             *GetListOfMatrices() const     {return fMatrices;}
   TList                 *GetListOfMaterials() const    {return fMaterials;}
   TList                 *GetListOfMedia() const        {return fMedia;}
   TObjArray             *GetListOfVolumes() const      {return fVolumes;}
   TObjArray             *GetListOfGVolumes() const     {return fGVolumes;}
   TObjArray             *GetListOfShapes() const       {return fShapes;}
   TObjArray             *GetListOfGShapes() const      {return fGShapes;}
   TObjArray             *GetListOfUVolumes() const     {return fUniqueVolumes;}
   TObjArray             *GetListOfTracks() const       {return fTracks;}
   TObjArray             *GetListOfGDMLMatrices() const {return fGDMLMatrices;}
   TObjArray             *GetListOfOpticalSurfaces() const {return fOpticalSurfaces;}
   TObjArray             *GetListOfSkinSurfaces() const    {return fSkinSurfaces;}
   TObjArray             *GetListOfBorderSurfaces() const  {return fBorderSurfaces;}
   TObjArray             *GetListOfRegions() const      {return fRegions;}
   TGeoNavigatorArray    *GetListOfNavigators() const;
   TGeoElementTable      *GetElementTable();

   //--- modeler state getters/setters
   void                   DoBackupState();
   void                   DoRestoreState();
   TGeoNode              *GetNode(Int_t level) const  {return (TGeoNode*)fNodes->UncheckedAt(level);}
   Int_t                  GetNodeId() const           {return GetCurrentNavigator()->GetNodeId();}
   TGeoNode              *GetNextNode() const         {return GetCurrentNavigator()->GetNextNode();}
   TGeoNode              *GetMother(Int_t up=1) const {return GetCurrentNavigator()->GetMother(up);}
   TGeoHMatrix           *GetMotherMatrix(Int_t up=1) const {return GetCurrentNavigator()->GetMotherMatrix(up);}
   TGeoHMatrix           *GetHMatrix();
   TGeoHMatrix           *GetCurrentMatrix() const    {return GetCurrentNavigator()->GetCurrentMatrix();}
   TGeoHMatrix           *GetGLMatrix() const         {return fGLMatrix;}
   TGeoNavigator         *GetCurrentNavigator() const;
   TGeoNode              *GetCurrentNode() const      {return GetCurrentNavigator()->GetCurrentNode();}
   Int_t                  GetCurrentNodeId() const;
   const Double_t        *GetCurrentPoint() const     {return GetCurrentNavigator()->GetCurrentPoint();}
   const Double_t        *GetCurrentDirection() const {return GetCurrentNavigator()->GetCurrentDirection();}
   TGeoVolume            *GetCurrentVolume() const {return GetCurrentNavigator()->GetCurrentVolume();}
   const Double_t        *GetCldirChecked() const  {return GetCurrentNavigator()->GetCldirChecked();}
   const Double_t        *GetCldir() const         {return GetCurrentNavigator()->GetCldir();}
   const Double_t        *GetNormal() const        {return GetCurrentNavigator()->GetNormal();}
   Int_t                  GetLevel() const         {return GetCurrentNavigator()->GetLevel();}
   Int_t                  GetMaxLevel() const      {return fNLevel;}
   const char            *GetPath() const;
   Int_t                  GetStackLevel() const    {return GetCurrentNavigator()->GetStackLevel();}
   TGeoVolume            *GetMasterVolume() const  {return fMasterVolume;}
   TGeoVolume            *GetTopVolume() const     {return fTopVolume;}
   TGeoNode              *GetTopNode() const       {return fTopNode;}
   TGeoPhysicalNode      *GetPhysicalNode(Int_t i) const {return (TGeoPhysicalNode*)fPhysicalNodes->UncheckedAt(i);}
   void                   SetCurrentPoint(Double_t *point) {GetCurrentNavigator()->SetCurrentPoint(point);}
   void                   SetCurrentPoint(Double_t x, Double_t y, Double_t z) {GetCurrentNavigator()->SetCurrentPoint(x,y,z);}
   void                   SetLastPoint(Double_t x, Double_t y, Double_t z) {GetCurrentNavigator()->SetLastPoint(x,y,z);}
   void                   SetCurrentDirection(Double_t *dir) {GetCurrentNavigator()->SetCurrentDirection(dir);}
   void                   SetCurrentDirection(Double_t nx, Double_t ny, Double_t nz) {GetCurrentNavigator()->SetCurrentDirection(nx,ny,nz);}
   void                   SetCldirChecked(Double_t *dir) {GetCurrentNavigator()->SetCldirChecked(dir);}

   //--- point/vector reference frame conversion
   void                   LocalToMaster(const Double_t *local, Double_t *master) const {GetCurrentNavigator()->LocalToMaster(local, master);}
   void                   LocalToMasterVect(const Double_t *local, Double_t *master) const {GetCurrentNavigator()->LocalToMasterVect(local, master);}
   void                   LocalToMasterBomb(const Double_t *local, Double_t *master) const {GetCurrentNavigator()->LocalToMasterBomb(local, master);}
   void                   MasterToLocal(const Double_t *master, Double_t *local) const {GetCurrentNavigator()->MasterToLocal(master, local);}
   void                   MasterToLocalVect(const Double_t *master, Double_t *local) const {GetCurrentNavigator()->MasterToLocalVect(master, local);}
   void                   MasterToLocalBomb(const Double_t *master, Double_t *local) const {GetCurrentNavigator()->MasterToLocalBomb(master, local);}
   void                   MasterToTop(const Double_t *master, Double_t *top) const;
   void                   TopToMaster(const Double_t *top, Double_t *master) const;

   //--- general use getters/setters
   TGeoMaterial          *FindDuplicateMaterial(const TGeoMaterial *mat) const;
   TGeoVolume            *FindVolumeFast(const char*name, Bool_t multi=kFALSE);
   TGeoMaterial          *GetMaterial(const char *matname) const;
   TGeoMaterial          *GetMaterial(Int_t id) const;
   TGeoMedium            *GetMedium(const char *medium) const;
   TGeoMedium            *GetMedium(Int_t numed) const;
   Int_t                  GetMaterialIndex(const char *matname) const;

   //--- GDML object accessors
   TGDMLMatrix           *GetGDMLMatrix(const char *name) const;
   void                   AddGDMLMatrix(TGDMLMatrix *mat);
   TGeoOpticalSurface    *GetOpticalSurface(const char *name) const;
   void                   AddOpticalSurface(TGeoOpticalSurface *optsurf);
   TGeoSkinSurface       *GetSkinSurface(const char *name) const;
   void                   AddSkinSurface(TGeoSkinSurface *surf);
   TGeoBorderSurface     *GetBorderSurface(const char *name) const;
   void                   AddBorderSurface(TGeoBorderSurface *surf);

//   TGeoShape             *GetShape(const char *name) const;
   TGeoVolume            *GetVolume(const char*name) const;
   TGeoVolume            *GetVolume(Int_t uid) const {return (TGeoVolume*)fUniqueVolumes->At(uid);}
   int                    GetNregions() const {return fRegions->GetEntriesFast();}
   TGeoRegion            *GetRegion(int i) {return (TGeoRegion*)fRegions->At(i);}
   Int_t                  GetUID(const char *volname) const;
   Int_t                  GetNNodes() {if (!fNNodes) CountNodes(); return fNNodes;}
   TGeoNodeCache         *GetCache() const         {return GetCurrentNavigator()->GetCache();}
//   void                   SetCache(const TGeoNodeCache *cache) {fCache = (TGeoNodeCache*)cache;}
   void                   SetAnimateTracks(Bool_t flag=kTRUE) {fIsGeomReading=flag;}
   virtual ULong_t        SizeOf(const TGeoNode *node, Option_t *option); // size of the geometry in memory
   void                   SelectTrackingMedia();

   //--- stack manipulation
   Int_t                  PushPath(Int_t startlevel=0) {return GetCurrentNavigator()->PushPath(startlevel);}
   Bool_t                 PopPath() {return GetCurrentNavigator()->PopPath();}
   Bool_t                 PopPath(Int_t index) {return GetCurrentNavigator()->PopPath(index);}
   Int_t                  PushPoint(Int_t startlevel=0) {return GetCurrentNavigator()->PushPoint(startlevel);}
   Bool_t                 PopPoint() {return GetCurrentNavigator()->PopPoint();}
   Bool_t                 PopPoint(Int_t index) {return GetCurrentNavigator()->PopPoint(index);}
   void                   PopDummy(Int_t ipop=9999) {return GetCurrentNavigator()->PopDummy(ipop);}

   //--- parallel world navigation
   TGeoParallelWorld    *CreateParallelWorld(const char *name);
   TGeoParallelWorld    *GetParallelWorld() const {return fParallelWorld;}
   void                  SetUseParallelWorldNav(Bool_t flag);
   Bool_t                IsParallelWorldNav() const {return fUsePWNav;}

   ClassDef(TGeoManager, 17)          // geometry manager
};

R__EXTERN TGeoManager *gGeoManager;

#endif

