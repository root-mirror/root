// @(#)root/gl:$Id$
// Author:  Richard Maunder  25/05/2005

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLCamera
#define ROOT_TGLCamera

#include "TGLUtil.h"
#include "TGLBoundingBox.h"
#include "TPoint.h"
#include "TObject.h"

#include <cmath>
#include <utility>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLCamera                                                            //
//                                                                      //
// Abstract base camera class - concrete classes for orthographic and   //
// persepctive cameras derive from it. This class maintains values for  //
// the current:                                                         //
// i)   Viewport                                                        //
// ii)  Projection, modelview and clip matricies - extracted from GL    //
// iii) The 6 frustum planes                                            //
// iv)  Expanded frustum interest box                                   //
//                                                                      //
// It provides methods for various projection, overlap and intersection //
// tests for viewport and world locations, against the true frustum and //
// expanded interest box, and for extracting eye position and direction.//
//                                                                      //
// It also defines the pure virtual manipulation interface methods the  //
// concrete ortho and prespective classes must implement.               //
//////////////////////////////////////////////////////////////////////////

class TGLCamera : public TObject
{
public:
   enum EFrustumPlane
   {
      kNear             = 0,
      kLeft             = 1,
      kRight            = 2,
      kTop              = 3,
      kBottom           = 4,
      kFar              = 5,
      kPlanesPerFrustum = 6
   };

private:
   // Fields

   // Debuging visual aids
   TGLBoundingBox   fPreviousInterestBox;  //! previous interest box (DEBUG)
   TGLBoundingBox   fInterestFrustum;      //! frustum basis of current interest box - NOT a true BB! (DEBUG)
   TGLBoundingBox   fInterestFrustumAsBox; //! frustum basis (as box) of current interest box (DEBUG)

   static const Double_t fgInterestBoxExpansion; //! expansion c.f. aligned current frustum box

   // Methods
   TGLBoundingBox Frustum(Bool_t asBox = kTRUE) const; // current frustum

   // Non-copyable class
   TGLCamera(const TGLCamera &);
   TGLCamera & operator=(const TGLCamera &);

protected:
    // Fields
    TGLMatrix   fCamBase;         //  tranformation to center and rotation from up to x vector
    TGLMatrix   fCamTrans;        //  transformation relative to fCamTrans
    Bool_t      fExternalCenter;  //  use external center insead of scene center
    Bool_t      fFixDefCenter;    //  use fixed default center
    Bool_t      fWasArcBalled;    //  set when arc-ball rotation is used
    TGLVector3  fExtCenter;       //  external camera center
    TGLVector3  fDefCenter;       //  default camera center
    TGLVector3  fFDCenter;        //  fixed default camera center
    TGLVector3 *fCenter;          //! current camera center

    mutable Double_t fNearClip;   //! last applied near-clip
    mutable Double_t fFarClip;    //! last applied far-clip

    // Set in Setup()
    Double_t    fDollyDefault;    // default distnce from viewing centre
    Double_t    fDollyDistance;   // unit distance for camera movement in fwd/bck direction
    Float_t     fVAxisMinAngle;   // minimal allowed angle between up and fCamTrans Z vector

   // Internal cached matrices and frustum planes
   mutable Bool_t    fCacheDirty;                      //! cached items dirty?
   mutable UInt_t    fTimeStamp;                       //! timestamp
   mutable TGLMatrix fLastNoPickProjM;                 //! no-pick projection matrix (cached)
   mutable TGLMatrix fProjM;                           //! projection matrix        (cached)
   mutable TGLMatrix fModVM;                           //! modelView matrix         (cached)
   mutable TGLMatrix fClipM;                           //! object space clip matrix (cached)
   mutable TGLPlane fFrustumPlanes[kPlanesPerFrustum]; //! frustum planes           (cached)

   TGLRect   fViewport;    //! viewport (GL coords - origin bottom left)

   TGLBoundingBox   fInterestBox;          //! the interest box - created in UpdateInterest()
   mutable Double_t fLargestSeen;          //! largest box diagonal seen in OfInterest() - used when
                                           //! bootstrapping interest box

   // Internal cache update - const as the actual camera configuration is unaffected
   void       UpdateCache() const;

   static     UInt_t   fgDollyDeltaSens;
public:
   TGLCamera();
   TGLCamera(const TGLVector3 & hAxis, const TGLVector3 & vAxis);
   virtual ~TGLCamera();

   virtual Bool_t IsOrthographic() const {  return kFALSE; }
   virtual Bool_t IsPerspective() const { return kFALSE; }

   const TGLMatrix& RefModelViewMatrix() const { return fModVM; }

   Bool_t IsCacheDirty() const { return fCacheDirty; }
   void   IncTimeStamp()       { fCacheDirty = kTRUE; ++fTimeStamp; }
   UInt_t TimeStamp()    const { return fTimeStamp; }

   void           SetViewport(const TGLRect & viewport);
   TGLRect&       RefViewport()       { return fViewport; }
   const TGLRect& RefViewport() const { return fViewport; }

   // Camera manipulation interface (GL coord - origin bottom left)
   virtual void   Setup(const TGLBoundingBox & box, Bool_t reset=kTRUE) = 0;
   virtual void   Reset() = 0;

   virtual Bool_t Dolly(Int_t delta, Bool_t mod1, Bool_t mod2);
   virtual Bool_t Zoom (Int_t delta, Bool_t mod1, Bool_t mod2) = 0;
   virtual Bool_t Truck(Double_t xDelta, Double_t yDelta);
   virtual Bool_t Truck(Int_t xDelta, Int_t yDelta, Bool_t mod1, Bool_t mod2) = 0;
   virtual Bool_t Rotate(Int_t xDelta, Int_t yDelta, Bool_t mod1, Bool_t mod2);
   virtual Bool_t RotateRad(Double_t hRotate, Double_t vRotate);
   virtual Bool_t RotateArcBall(Int_t xDelta, Int_t yDelta, Bool_t mod1, Bool_t mod2);
   virtual Bool_t RotateArcBallRad(Double_t hRotate, Double_t vRotate);

   virtual void   Apply(const TGLBoundingBox & sceneBox, const TGLRect * pickRect = 0) const = 0;

   Bool_t     AdjustAndClampVal(Double_t & val, Double_t min, Double_t max,
                                Int_t screenShift, Int_t screenShiftRange,
                                Bool_t mod1, Bool_t mod2) const;
   Double_t   AdjustDelta(Double_t screenShift, Double_t deltaFactor,
                          Bool_t mod1, Bool_t mod2) const;

   void    SetExternalCenter(Bool_t x);
   Bool_t  GetExternalCenter(){ return fExternalCenter; }

   void    SetCenterVec(Double_t x, Double_t y, Double_t z);
   void    SetCenterVecWarp(Double_t x, Double_t y, Double_t z);
   Double_t* GetCenterVec() { return fCenter->Arr(); }

   void    SetFixDefCenter(Bool_t x) { fFixDefCenter = x; }
   void    SetFixDefCenterVec(Double_t x, Double_t y, Double_t z) { fFDCenter.Set(x, y, z); }
   Double_t* GetFixDefCenterVec() { return fFDCenter.Arr(); }

   Double_t GetNearClip() const { return fNearClip; }
   Double_t GetFarClip()  const { return fFarClip;  }

   const TGLMatrix& GetCamBase()  const { return fCamBase;  }
   const TGLMatrix& GetCamTrans() const { return fCamTrans; }
   // If you manipulate camera ... also call IncTimeStamp() before redraw.
   TGLMatrix& RefCamBase()  { return fCamBase;  }
   TGLMatrix& RefCamTrans() { return fCamTrans; }

   Double_t GetTheta() const;

   TGLMatrix& RefLastNoPickProjM() const { return fLastNoPickProjM; }

   // Current orientation and frustum
   TGLVertex3 EyePoint() const;
   TGLVector3 EyeDirection() const;
   TGLVertex3 FrustumCenter() const;
   const TGLPlane & FrustumPlane(EFrustumPlane plane) const;

   // Overlap / projection / intersection tests
   // Viewport is GL coorinate system - origin bottom/left
   Rgl::EOverlap FrustumOverlap (const TGLBoundingBox & box) const; // box/frustum overlap test
   Rgl::EOverlap ViewportOverlap(const TGLBoundingBox & box) const; // box/viewport overlap test
   TGLRect    ViewportRect   (const TGLBoundingBox & box, TGLBoundingBox::EFace face) const;
   TGLRect    ViewportRect   (const TGLBoundingBox & box, const TGLBoundingBox::EFace * face = 0) const;
   TGLVertex3 WorldToViewport(const TGLVertex3 & worldVertex, TGLMatrix* modviewMat=0) const;
   TGLVector3 WorldDeltaToViewport(const TGLVertex3 & worldRef, const TGLVector3 & worldDelta) const;
   TGLVertex3 ViewportToWorld(const TGLVertex3 & viewportVertex, TGLMatrix* modviewMat=0) const;
   TGLLine3   ViewportToWorld(Double_t viewportX, Double_t viewportY) const;
   TGLLine3   ViewportToWorld(const TPoint & viewport) const;
   TGLVector3 ViewportDeltaToWorld(const TGLVertex3 & worldRef, Double_t viewportXDelta, Double_t viewportYDelta, TGLMatrix* modviewMat=0) const;
   std::pair<Bool_t, TGLVertex3> ViewportPlaneIntersection(Double_t viewportX, Double_t viewportY, const TGLPlane & worldPlane) const;
   std::pair<Bool_t, TGLVertex3> ViewportPlaneIntersection(const TPoint & viewport, const TGLPlane & worldPlane) const;

   // Window to GL viewport conversion - invert Y
   void WindowToViewport(Int_t & /* x */, Int_t & y) const { y = fViewport.Height() - y; }
   void WindowToViewport(TPoint & point)             const { point.SetY(fViewport.Height() - point.GetY()); }
   void WindowToViewport(TGLRect & rect)             const { rect.Y() = fViewport.Height() - rect.Y(); }
   void WindowToViewport(TGLVertex3 & vertex)        const { vertex.Y() = fViewport.Height() - vertex.Y(); }

   Float_t GetVAxisMinAngle(){return fVAxisMinAngle;}
   void    SetVAxisMinAngle(Float_t x){fVAxisMinAngle = x;}

   virtual void Configure(Double_t zoom, Double_t dolly, Double_t center[3],
                          Double_t hRotate, Double_t vRotate) = 0;
   // Cameras expanded-frustum interest box
   Bool_t OfInterest(const TGLBoundingBox & box, Bool_t ignoreSize) const;
   Bool_t UpdateInterest(Bool_t force);
   void   ResetInterest();

   // Debuging - draw frustum and interest boxes
   void  DrawDebugAids() const;

   ClassDef(TGLCamera,1); // Camera abstract base class.
};

inline const TGLPlane & TGLCamera::FrustumPlane(EFrustumPlane plane) const
{
   // Return one of the planes forming the camera frustum
   if (fCacheDirty) {
      Error("TGLCamera::FrustumBox()", "cache dirty");
   }
   return fFrustumPlanes[plane];
}


#endif // ROOT_TGLCamera
