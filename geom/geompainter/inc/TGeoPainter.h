// Author: Andrei Gheata   05/03/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TGeoPainter
#define ROOT_TGeoPainter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGeoPainter                                                          //
//                                                                      //
// Painter class utility TGeo geometries. Interfaces visualization      //
// queries with the viewers.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TVirtualGeoPainter.h"

#include "TGeoManager.h"

class TString;
class TGeoHMatrix;
class TGeoNode;
class TGeoVolume;
class TGeoShape;
class TVirtualGeoTrack;
class TGeoPhysicalNode;
class TGeoChecker;
class TGeoOverlap;
class TH2F;
class TGeoBatemanSol;
class TGeoPolygon;

class TGeoPainter : public TVirtualGeoPainter {
private:
   Double_t           fBombX;            // bomb factor on X
   Double_t           fBombY;            // bomb factor on Y
   Double_t           fBombZ;            // bomb factor on Z
   Double_t           fBombR;            // bomb factor on radius (cyl or sph)
   Double_t           fCheckedBox[6];    // bounding box of checked node
   Double_t           fMat[9];           // view rotation matrix
   Int_t              fNsegments;        // number of segments approximating circles
   Int_t              fNVisNodes;        // number of visible nodes
   Int_t              fVisLevel;         // depth for drawing
   Int_t              fVisOption;        // global visualization option
   Int_t              fExplodedView;     // type of exploding current view
   Bool_t             fVisLock;          // lock for adding visible volumes
   Bool_t             fTopVisible;       // set top volume visible
   Bool_t             fPaintingOverlaps; // lock overlaps painting
   Bool_t             fIsRaytracing;     // raytracing flag
   Bool_t             fIsPaintingShape;  // flag for shape painting
   TString            fVisBranch;        // drawn branch
   TString            fVolInfo;          // volume info
   TGeoNode          *fCheckedNode;      // checked node
   TGeoOverlap       *fOverlap;          // current overlap
   TGeoHMatrix       *fGlobal;           // current global matrix
   TBuffer3D         *fBuffer;           // buffer used for painting
   TGeoManager       *fGeoManager;       // geometry to which applies
   TGeoChecker       *fChecker;          // geometry checker
   TGeoShape         *fClippingShape;    // clipping shape
   TGeoVolume        *fTopVolume;        // top drawn volume
   TGeoVolume        *fLastVolume;       // last drawn volume
   TGeoIteratorPlugin *fPlugin;          // User iterator plugin for changing pain volume properties
   TObjArray         *fVisVolumes;       // list of visible volumes
   Bool_t             fIsEditable;       // flag that geometry is editable

   void       DefineColors() const;
   void       LocalToMasterVect(const Double_t *local, Double_t *master) const;

protected:
   void       ClearVisibleVolumes();

public:
   TGeoPainter(TGeoManager *manager);
   virtual ~TGeoPainter();

   void       AddSize3D(Int_t numpoints, Int_t numsegs, Int_t numpolys) override;
   TVirtualGeoTrack *AddTrack(Int_t id, Int_t pdgcode, TObject *part) override;
   void       AddTrackPoint(Double_t *point, Double_t *box, Bool_t reset=kFALSE) override;
   void       BombTranslation(const Double_t *tr, Double_t *bombtr) override;
   void       CheckBoundaryErrors(Int_t ntracks=1000000, Double_t radius=-1.) override;
   void       CheckBoundaryReference(Int_t icheck=-1) override;
   void       CheckGeometryFull(Bool_t checkoverlaps=kTRUE, Bool_t checkcrossings=kTRUE, Int_t nrays=10000, const Double_t *vertex=nullptr) override;
   void       CheckGeometry(Int_t nrays, Double_t startx, Double_t starty, Double_t startz) const override;
   void       CheckEdit();
   void       CheckPoint(Double_t x=0, Double_t y=0, Double_t z=0, Option_t *option="") override;
   void       CheckShape(TGeoShape *shape, Int_t testNo, Int_t nsamples, Option_t *option) override;
   void       CheckOverlaps(const TGeoVolume *vol, Double_t ovlp=0.1, Option_t *option="") const override;
   Int_t      CountNodes(TGeoVolume *vol, Int_t level) const;
   Int_t      CountVisibleNodes() override;
   void       DefaultAngles() override;
   void       DefaultColors() override;
   Int_t      DistanceToPrimitiveVol(TGeoVolume *vol, Int_t px, Int_t py) override;
   void       Draw(Option_t *option="") override;
   void       DrawBatemanSol(TGeoBatemanSol *sol, Option_t *option="") override;
   void       DrawOverlap(void *ovlp, Option_t *option="") override;
   void       DrawCurrentPoint(Int_t color) override;
   void       DrawOnly(Option_t *option="") override;
   void       DrawPanel() override;
   void       DrawPath(const char *path, Option_t *option="") override;
   void       DrawPolygon(const TGeoPolygon *poly) override;
   void       DrawShape(TGeoShape *shape, Option_t *option="") override;
   void       DrawVolume(TGeoVolume *vol, Option_t *option="") override;
   void       EditGeometry(Option_t *option="") override;
   void       EstimateCameraMove(Double_t tmin, Double_t tmax, Double_t *start, Double_t *end) override;
   void       ExecuteManagerEvent(TGeoManager *geom, Int_t event, Int_t px, Int_t py) override;
   void       ExecuteShapeEvent(TGeoShape *shape, Int_t event, Int_t px, Int_t py) override;
   void       ExecuteVolumeEvent(TGeoVolume *volume, Int_t event, Int_t px, Int_t py) override;
   const char *GetVolumeInfo(const TGeoVolume *volume, Int_t px, Int_t py) const override;
   void       GetBombFactors(Double_t &bombx, Double_t &bomby, Double_t &bombz, Double_t &bombr) const override
                                    {bombx=fBombX; bomby=fBombY; bombz=fBombZ; bombr=fBombR;}
   Int_t      GetBombMode() const  override {return fExplodedView;}
   TGeoNode  *GetCheckedNode() {return fCheckedNode;}
   TGeoChecker       *GetChecker();
   Int_t      GetColor(Int_t base, Float_t light) const override;
   const char *GetDrawPath() const  override {return fVisBranch.Data();}
   TGeoVolume *GetDrawnVolume() const override;
   TGeoVolume *GetTopVolume() const override {return fTopVolume;}
   Int_t      GetVisLevel() const  override     {return fVisLevel;}
   Int_t      GetVisOption() const  override    {return fVisOption;}
   Int_t      GetNsegments() const   override   {return fNsegments;}
   void       GrabFocus(Int_t nfr=0, Double_t dlong=0, Double_t dlat=0, Double_t dpsi=0) override;
   Double_t  *GetViewBox() override {return &fCheckedBox[0];}
   void       GetViewAngles(Double_t &longitude, Double_t &latitude, Double_t &psi) override;
   Bool_t     IsExplodedView() const  override {return (fExplodedView==kGeoVisDefault)?kFALSE:kTRUE;}
   Bool_t     IsRaytracing() const override {return fIsRaytracing;}
   Bool_t     IsPaintingShape() const  override {return fIsPaintingShape;}
   TH2F      *LegoPlot(Int_t ntheta=60, Double_t themin=0., Double_t themax=180.,
                       Int_t nphi=90, Double_t phimin=0., Double_t phimax=360.,
                       Double_t rmin=0., Double_t rmax=9999999, Option_t *option="") override;
   void       Lock(Bool_t flag = kTRUE) { fVisLock = flag; }
   void       ModifiedPad(Bool_t update=kFALSE) const override;
   void       OpProgress(const char *opname, Long64_t current, Long64_t size, TStopwatch *watch=nullptr, Bool_t last=kFALSE, Bool_t refresh=kFALSE, const char *msg="") override;
   void       Paint(Option_t *option="") override;
   void       PaintNode(TGeoNode *node, Option_t *option="", TGeoMatrix *global=nullptr) override;
   Bool_t     PaintShape(const TGeoShape & shape, Option_t * option) const;
   void       PaintShape(TGeoShape *shape, Option_t *option="") override;
   void       PaintOverlap(void *ovlp, Option_t *option="") override;
   void       PaintVolume(TGeoVolume *vol, Option_t *option="", TGeoMatrix *global=nullptr) override;
   void       PrintOverlaps() const override;
   void       PaintPhysicalNode(TGeoPhysicalNode *node, Option_t *option="");
   void       RandomPoints(const TGeoVolume *vol, Int_t npoints, Option_t *option="") override;
   void       RandomRays(Int_t nrays, Double_t startx, Double_t starty, Double_t startz, const char *target_vol = nullptr, Bool_t check_norm = kFALSE) override;
   void       Raytrace(Option_t *option="") override;
   TGeoNode  *SamplePoints(Int_t npoints, Double_t &dist, Double_t epsil, const char* g3path) override;
   void       SetBombFactors(Double_t bombx=1.3, Double_t bomby=1.3, Double_t bombz=1.3, Double_t bombr=1.3) override;
   void       SetClippingShape(TGeoShape *shape) override {fClippingShape = shape;}
   void       SetExplodedView(Int_t iopt=0) override;
   void       SetNsegments(Int_t nseg=20) override;
   void       SetNmeshPoints(Int_t npoints) override;
   void       SetGeoManager(TGeoManager *geom) override {fGeoManager=geom;}
   void       SetIteratorPlugin(TGeoIteratorPlugin *plugin) override {fPlugin = plugin; ModifiedPad();}
   void       SetCheckedNode(TGeoNode *node) override;
   void       SetRaytracing(Bool_t flag=kTRUE) override {fIsRaytracing = flag;}
   void       SetTopVisible(Bool_t vis=kTRUE) override;
   void       SetTopVolume(TGeoVolume *vol) override {fTopVolume = vol;}
   void       SetVisLevel(Int_t level=3) override;
   void       SetVisOption(Int_t option=0) override;
   Int_t      ShapeDistancetoPrimitive(const TGeoShape *shape, Int_t numpoints, Int_t px, Int_t py) const override;
   void       Test(Int_t npoints, Option_t *option) override;
   void       TestOverlaps(const char *path) override;
   Bool_t     TestVoxels(TGeoVolume *vol) override;
   void       UnbombTranslation(const Double_t *tr, Double_t *bombtr) override;
   Double_t   Weight(Double_t precision, Option_t *option="v") override;

   ClassDefOverride(TGeoPainter,0)  //geometry painter
};

#endif
