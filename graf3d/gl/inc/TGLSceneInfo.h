// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Feb 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLSceneInfo_H
#define ROOT_TGLSceneInfo_H

#include "Rtypes.h"

#include "TGLBoundingBox.h"
#include "TGLUtil.h"
#include <vector>

class TGLViewerBase;
class TGLSceneBase;
class TGLClip;
class TGLRenderContext;
class TGLCamera;

class TGLSceneInfo
{
   friend class TGLSceneBase;

public:
   enum EClipTest { kClipNone, kClipOutside, kClipInside };

private:
   TGLSceneInfo(const TGLSceneInfo&) = delete;
   TGLSceneInfo& operator=(const TGLSceneInfo&) = delete;

protected:
   TGLViewerBase  * fViewer;
   TGLSceneBase   * fScene;
   Bool_t           fActive;      // Show fScene in fViewer

   Short_t          fLOD;         // Optional override of scene lod
   Short_t          fStyle;       // Optional override of scene style
   Float_t          fWFLineW;     // Optional override of scene wire-frame line-width
   Float_t          fOLLineW;     // Optional override of scene outline line-width
   TGLClip        * fClip;        // Optional override of clipping-plane

   Short_t          fLastLOD;     // Last combined viewer/scene lod   (set in scene::lodify-scene-info).
   Short_t          fLastStyle;   // Last combined viewer/scene style (set in scene::pre-draw).
   Float_t          fLastWFLineW; // Last combined viewer/scene wire-frame line-width (set in scene::pre-draw).
   Float_t          fLastOLLineW; // Last combined viewer/scene outline line-width (set in scene::pre-draw).
   TGLClip        * fLastClip;    // Last combined viewer/scene clip  (set in scene::update)
   TGLCamera      * fLastCamera;  // Last camera used.

   UInt_t           fSceneStamp;  // Scene's time-stamp on last update.
   UInt_t           fClipStamp;   // Clip's time-stamp on last update.
   UInt_t           fCameraStamp; // Camera's time-stamp on last update.
   Bool_t           fUpdateTimeouted; // Set if update was interrupted.

   // Eventually we will allow additional per-scene transforamtion.
   // TGLMatrix  fSceneTrans;
   // Also needed:
   // *) transformed clipping planes of the camera
   // *) transformed bounding-box of the scene
   TGLBoundingBox   fTransformedBBox;

   Bool_t           fViewCheck;     // Viewer side check if render is necessary.
   Bool_t           fInFrustum;     // Is scene intersecting view-frustum.
   Bool_t           fInClip;        // Is scene contained within clipping-volume.
   Char_t           fClipMode;      // Clipping mode, can be disbled.
   TGLPlaneSet_t    fFrustumPlanes; // Clipping planes defined by frustum; only those intersecting the scene volume are kept.
   TGLPlaneSet_t    fClipPlanes;    // Clipping planes from clip-object; which planes are kept depends on inside/outside mode.

   // Additional stuff (scene-class specific) can be added by sub-classing.
   // For TGLScene these data include draw-lists after clipping.

public:
   TGLSceneInfo(TGLViewerBase* view=0, TGLSceneBase* scene=0);
   virtual ~TGLSceneInfo() {}

   TGLViewerBase * GetViewer() const { return  fViewer; }
   TGLViewerBase & RefViewer() const { return *fViewer; }
   TGLSceneBase  * GetScene()  const { return  fScene;  }
   TGLSceneBase  & RefScene()  const { return *fScene;  }

   Bool_t GetActive() const { return fActive; }
   void   SetActive(Bool_t a);

   void  SetupTransformsAndBBox();

   const TGLBoundingBox& GetTransformedBBox() { return fTransformedBBox; }

   virtual void SetSceneTrans(TGLMatrix&) { ResetSceneStamp(); }

   Bool_t   ViewCheck() const   { return fViewCheck; }
   void     ViewCheck(Bool_t c) { fViewCheck = c;    }
   Bool_t   IsInFrustum() const { return fInFrustum; }
   void     InFrustum(Bool_t f) { fInFrustum = f;    }
   Bool_t   IsInClip()    const { return fInClip;    }
   void     InClip(Bool_t c)    { fInClip = c;       }
   Char_t   ClipMode()    const { return fClipMode;  }
   void     ClipMode(Char_t m)  { fClipMode = m;     }

   Bool_t   ShouldClip()  const { return fClipMode != kClipNone; }
   Bool_t   IsVisible()   const { return fInFrustum && fInClip;  }

   std::vector<TGLPlane>& FrustumPlanes() { return fFrustumPlanes; }
   std::vector<TGLPlane>& ClipPlanes()    { return fClipPlanes;    }

   Short_t  LOD()          const { return fLOD; }
   void     SetLOD(Short_t lod)  { fLOD = lod;  }

   Short_t  Style()        const { return fStyle; }
   void     SetStyle(Short_t st) { fStyle = st;   }

   Float_t  WFLineW()       const { return fWFLineW; }
   void     SetWFLineW(Float_t w) { fWFLineW = w;    }
   Float_t  OLLineW()       const { return fOLLineW; }
   void     SetOLLineW(Float_t w) { fOLLineW = w;    }

   TGLClip* Clip()         const { return fClip; }
   void     SetClip(TGLClip *p)  { fClip = p;    }

   Short_t  LastLOD()        const { return fLastLOD; }
   void     SetLastLOD(Short_t ld) { fLastLOD = ld;   }

   Short_t  LastStyle()      const   { return fLastStyle; }
   void     SetLastStyle(Short_t st) { fLastStyle = st;   }

   Float_t  LastWFLineW()       const { return fLastWFLineW; }
   void     SetLastWFLineW(Float_t w) { fLastWFLineW = w;    }
   Float_t  LastOLLineW()       const { return fLastOLLineW; }
   void     SetLastOLLineW(Float_t w) { fLastOLLineW = w;    }

   TGLClip* LastClip()         const { return fLastClip; }
   void     SetLastClip(TGLClip *p)  { fLastClip = p;    }

   TGLCamera* LastCamera()        const { return fLastCamera; }
   void     SetLastCamera(TGLCamera *p) { fLastCamera = p;    }

   UInt_t   SceneStamp()       const { return fSceneStamp; }
   void     SetSceneStamp(UInt_t ts) { fSceneStamp = ts;   }
   void     ResetSceneStamp()        { fSceneStamp = 0;    }

   UInt_t   ClipStamp()       const { return fClipStamp; }
   void     SetClipStamp(UInt_t ts) { fClipStamp = ts;   }
   void     ResetClipStamp()        { fClipStamp = 0;    }

   UInt_t   CameraStamp()       const { return fCameraStamp; }
   void     SetCameraStamp(UInt_t ts) { fCameraStamp = ts;   }
   void     ResetCameraStamp()        { fCameraStamp = 0;    }

   Bool_t   HasUpdateTimeouted() const { return fUpdateTimeouted;   }
   void     UpdateTimeouted()          { fUpdateTimeouted = kTRUE;  }
   void     ResetUpdateTimeouted()     { fUpdateTimeouted = kFALSE; }

   ClassDef(TGLSceneInfo, 0) // Data about a scene within a viewer context.
}; // endclass TGLSceneInfo


#endif
