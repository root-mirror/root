// @(#)root/gl:$Id$
// Author:  Richard Maunder  25/05/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLViewer.h"
#include "TGLIncludes.h"
#include "TGLStopwatch.h"
#include "TGLRnrCtx.h"
#include "TGLSelectBuffer.h"
#include "TGLLightSet.h"
#include "TGLManipSet.h"
#include "TGLCameraOverlay.h"
#include "TGLAutoRotator.h"

#include "TGLScenePad.h"
#include "TGLLogicalShape.h"
#include "TGLPhysicalShape.h"
#include "TGLObject.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

#include "TGLOutput.h"

#include "TROOT.h"
#include "TVirtualMutex.h"

#include "TVirtualPad.h" // Remove when pad removed - use signal
#include "TVirtualX.h"

#include "TMath.h"
#include "TColor.h"
#include "TError.h"
#include "TEnv.h"

// For event type translation ExecuteEvent
#include "Buttons.h"
#include "GuiTypes.h"

#include "TVirtualGL.h"

#include "TGLWidget.h"
#include "TGLFBO.h"
#include "TGLViewerEditor.h"
#include "TGedEditor.h"
#include "TGLPShapeObj.h"

#include "KeySymbols.h"
#include "TContextMenu.h"
#include "TImage.h"

#include <stdexcept>

#ifndef GL_BGRA
#define GL_BGRA GL_BGRA_EXT
#endif

/** \class TGLViewer
\ingroup opengl
Base GL viewer object - used by both standalone and embedded (in pad)
GL. Contains core viewer objects :

GL scene - collection of main drawn objects - see TGLStdScene
Cameras (fXyzzCamera) - ortho and perspective cameras - see TGLCamera
Clipping (fClipXyzz) - collection of clip objects - see TGLClip
Manipulators (fXyzzManip) - collection of manipulators - see TGLManip

It maintains the current active draw styles, clipping object,
manipulator, camera etc.

TGLViewer is 'GUI free' in that it does not derive from any ROOT GUI
TGFrame etc - see TGLSAViewer for this. However it contains GUI
GUI style methods HandleButton() etc to which GUI events can be
directed from standalone frame or embedding pad to perform
interaction.

Also, the TGLWidget needs to be created externally. It is not owned
by the viewer.

For embedded (pad) GL this viewer is created directly by plugin
manager. For standalone the derived TGLSAViewer is.
*/

ClassImp(TGLViewer);

TGLColorSet TGLViewer::fgDefaultColorSet;
Bool_t      TGLViewer::fgUseDefaultColorSetForNewViewers = kFALSE;

////////////////////////////////////////////////////////////////////////////////

TGLViewer::TGLViewer(TVirtualPad * pad, Int_t x, Int_t y,
                     Int_t width, Int_t height) :
   fPad(pad),
   fContextMenu(0),
   fPerspectiveCameraXOZ(TGLVector3(-1.0, 0.0, 0.0), TGLVector3(0.0, 1.0, 0.0)), // XOZ floor
   fPerspectiveCameraYOZ(TGLVector3( 0.0,-1.0, 0.0), TGLVector3(1.0, 0.0, 0.0)), // YOZ floor
   fPerspectiveCameraXOY(TGLVector3(-1.0, 0.0, 0.0), TGLVector3(0.0, 0.0, 1.0)), // XOY floor
   fOrthoXOYCamera (TGLOrthoCamera::kXOY,  TGLVector3( 0.0, 0.0, 1.0), TGLVector3(0.0, 1.0, 0.0)), // Looking down  Z axis,  X horz, Y vert
   fOrthoXOZCamera (TGLOrthoCamera::kXOZ,  TGLVector3( 0.0,-1.0, 0.0), TGLVector3(0.0, 0.0, 1.0)), // Looking along Y axis,  X horz, Z vert
   fOrthoZOYCamera (TGLOrthoCamera::kZOY,  TGLVector3(-1.0, 0.0, 0.0), TGLVector3(0.0, 1.0, 0.0)), // Looking along X axis,  Z horz, Y vert
   fOrthoZOXCamera (TGLOrthoCamera::kZOX,  TGLVector3( 0.0,-1.0, 0.0), TGLVector3(1.0, 0.0, 0.0)), // Looking along Y axis,  Z horz, X vert
   fOrthoXnOYCamera(TGLOrthoCamera::kXnOY, TGLVector3( 0.0, 0.0,-1.0), TGLVector3(0.0, 1.0, 0.0)), // Looking along Z axis, -X horz, Y vert
   fOrthoXnOZCamera(TGLOrthoCamera::kXnOZ, TGLVector3( 0.0, 1.0, 0.0), TGLVector3(0.0, 0.0, 1.0)), // Looking down  Y axis, -X horz, Z vert
   fOrthoZnOYCamera(TGLOrthoCamera::kZnOY, TGLVector3( 1.0, 0.0, 0.0), TGLVector3(0.0, 1.0, 0.0)), // Looking down  X axis, -Z horz, Y vert
   fOrthoZnOXCamera(TGLOrthoCamera::kZnOX, TGLVector3( 0.0, 1.0, 0.0), TGLVector3(1.0, 0.0, 0.0)), // Looking down  Y axis, -Z horz, X vert
   fCurrentCamera(&fPerspectiveCameraXOZ),
   fAutoRotator(0),

   fStereo               (kFALSE),
   fStereoQuadBuf        (kFALSE),
   fStereoZeroParallax   (0.03f),
   fStereoEyeOffsetFac   (1.0f),
   fStereoFrustumAsymFac (1.0f),

   fLightSet          (0),
   fClipSet           (0),
   fSelectedPShapeRef (0),
   fCurrentOvlElm     (0),

   fEventHandler(0),
   fGedEditor(0),
   fPShapeWrap(0),
   fPushAction(kPushStd), fDragAction(kDragNone),
   fRedrawTimer(0),
   fMaxSceneDrawTimeHQ(5000),
   fMaxSceneDrawTimeLQ(100),
   fPointScale (1), fLineScale(1), fSmoothPoints(kFALSE), fSmoothLines(kFALSE),
   fAxesType(TGLUtil::kAxesNone),
   fAxesDepthTest(kTRUE),
   fReferenceOn(kFALSE),
   fReferencePos(0.0, 0.0, 0.0),
   fDrawCameraCenter(kFALSE),
   fCameraOverlay(0),
   fSmartRefresh(kFALSE),
   fDebugMode(kFALSE),
   fIsPrinting(kFALSE),
   fPictureFileName("viewer.jpg"),
   fFader(0),
   fGLWidget(0),
   fGLDevice(-1),
   fGLCtxId(0),
   fIgnoreSizesOnUpdate(kFALSE),
   fResetCamerasOnUpdate(kTRUE),
   fResetCamerasOnNextUpdate(kFALSE)
{
   // Construct the viewer object, with following arguments:
   //    'pad' - external pad viewer is bound to
   //    'x', 'y' - initial top left position
   //    'width', 'height' - initial width/height

   InitSecondaryObjects();

   SetViewport(x, y, width, height);
}

////////////////////////////////////////////////////////////////////////////////

TGLViewer::TGLViewer(TVirtualPad * pad) :
   fPad(pad),
   fContextMenu(0),
   fPerspectiveCameraXOZ(TGLVector3(-1.0, 0.0, 0.0), TGLVector3(0.0, 1.0, 0.0)), // XOZ floor
   fPerspectiveCameraYOZ(TGLVector3( 0.0,-1.0, 0.0), TGLVector3(1.0, 0.0, 0.0)), // YOZ floor
   fPerspectiveCameraXOY(TGLVector3(-1.0, 0.0, 0.0), TGLVector3(0.0, 0.0, 1.0)), // XOY floor
   fOrthoXOYCamera (TGLOrthoCamera::kXOY,  TGLVector3( 0.0, 0.0, 1.0), TGLVector3(0.0, 1.0, 0.0)), // Looking down  Z axis,  X horz, Y vert
   fOrthoXOZCamera (TGLOrthoCamera::kXOZ,  TGLVector3( 0.0,-1.0, 0.0), TGLVector3(0.0, 0.0, 1.0)), // Looking along Y axis,  X horz, Z vert
   fOrthoZOYCamera (TGLOrthoCamera::kZOY,  TGLVector3(-1.0, 0.0, 0.0), TGLVector3(0.0, 1.0, 0.0)), // Looking along X axis,  Z horz, Y vert
   fOrthoZOXCamera (TGLOrthoCamera::kZOX,  TGLVector3( 0.0,-1.0, 0.0), TGLVector3(1.0, 0.0, 0.0)), // Looking along Y axis,  Z horz, X vert
   fOrthoXnOYCamera(TGLOrthoCamera::kXnOY, TGLVector3( 0.0, 0.0,-1.0), TGLVector3(0.0, 1.0, 0.0)), // Looking along Z axis, -X horz, Y vert
   fOrthoXnOZCamera(TGLOrthoCamera::kXnOZ, TGLVector3( 0.0, 1.0, 0.0), TGLVector3(0.0, 0.0, 1.0)), // Looking down  Y axis, -X horz, Z vert
   fOrthoZnOYCamera(TGLOrthoCamera::kZnOY, TGLVector3( 1.0, 0.0, 0.0), TGLVector3(0.0, 1.0, 0.0)), // Looking down  X axis, -Z horz, Y vert
   fOrthoZnOXCamera(TGLOrthoCamera::kZnOX, TGLVector3( 0.0, 1.0, 0.0), TGLVector3(1.0, 0.0, 0.0)), // Looking down  Y axis, -Z horz, X vert
   fCurrentCamera(&fPerspectiveCameraXOZ),
   fAutoRotator(0),

   fStereo               (kFALSE),
   fStereoQuadBuf        (kFALSE),
   fStereoZeroParallax   (0.03f),
   fStereoEyeOffsetFac   (1.0f),
   fStereoFrustumAsymFac (1.0f),

   fLightSet          (0),
   fClipSet           (0),
   fSelectedPShapeRef (0),
   fCurrentOvlElm     (0),

   fEventHandler(0),
   fGedEditor(0),
   fPShapeWrap(0),
   fPushAction(kPushStd), fDragAction(kDragNone),
   fRedrawTimer(0),
   fMaxSceneDrawTimeHQ(5000),
   fMaxSceneDrawTimeLQ(100),
   fPointScale (1), fLineScale(1), fSmoothPoints(kFALSE), fSmoothLines(kFALSE),
   fAxesType(TGLUtil::kAxesNone),
   fAxesDepthTest(kTRUE),
   fReferenceOn(kFALSE),
   fReferencePos(0.0, 0.0, 0.0),
   fDrawCameraCenter(kFALSE),
   fCameraOverlay(0),
   fSmartRefresh(kFALSE),
   fDebugMode(kFALSE),
   fIsPrinting(kFALSE),
   fPictureFileName("viewer.jpg"),
   fFader(0),
   fGLWidget(0),
   fGLDevice(fPad->GetGLDevice()),
   fGLCtxId(0),
   fIgnoreSizesOnUpdate(kFALSE),
   fResetCamerasOnUpdate(kTRUE),
   fResetCamerasOnNextUpdate(kFALSE)
{
   //gl-embedded viewer's ctor
   // Construct the viewer object, with following arguments:
   //    'pad' - external pad viewer is bound to
   //    'x', 'y' - initial top left position
   //    'width', 'height' - initial width/height

   InitSecondaryObjects();

   if (fGLDevice != -1) {
      // For the moment instantiate a fake context identity.
      fGLCtxId = new TGLContextIdentity;
      fGLCtxId->AddRef(0);
      Int_t viewport[4] = {0};
      gGLManager->ExtractViewport(fGLDevice, viewport);
      SetViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Common initialization.

void TGLViewer::InitSecondaryObjects()
{
   fLightSet = new TGLLightSet;
   fClipSet  = new TGLClipSet;
   AddOverlayElement(fClipSet);

   fSelectedPShapeRef = new TGLManipSet;
   fSelectedPShapeRef->SetDrawBBox(kTRUE);
   AddOverlayElement(fSelectedPShapeRef);

   fPShapeWrap = new TGLPShapeObj(0, this);

   fLightColorSet.StdLightBackground();
   if (fgUseDefaultColorSetForNewViewers) {
      fRnrCtx->ChangeBaseColorSet(&fgDefaultColorSet);
   } else {
      if (fPad) {
         fRnrCtx->ChangeBaseColorSet(&fLightColorSet);
         fLightColorSet.Background().SetColor(fPad->GetFillColor());
         fLightColorSet.Foreground().SetColor(fPad->GetLineColor());
      } else {
         fRnrCtx->ChangeBaseColorSet(&fDarkColorSet);
      }
   }

   fCameraOverlay = new TGLCameraOverlay(kFALSE, kFALSE);
   AddOverlayElement(fCameraOverlay);

   fRedrawTimer = new TGLRedrawTimer(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy viewer object.

TGLViewer::~TGLViewer()
{
   delete fAutoRotator;

   delete fLightSet;
   // fClipSet, fSelectedPShapeRef and fCameraOverlay deleted via overlay.

   delete fContextMenu;
   delete fRedrawTimer;

   if (fEventHandler) {
      if (fGLWidget)
         fGLWidget->SetEventHandler(0);
      delete fEventHandler;
   }

   if (fPad)
      fPad->ReleaseViewer3D();
   if (fGLDevice != -1)
      fGLCtxId->Release(0);
}


////////////////////////////////////////////////////////////////////////////////
/// Entry point for updating viewer contents via VirtualViewer3D
/// interface.
/// We search and forward the request to appropriate TGLScenePad.
/// If it is not found we create a new TGLScenePad so this can
/// potentially also be used for registration of new pads.

void TGLViewer::PadPaint(TVirtualPad* pad)
{
   TGLScenePad* scenepad = 0;
   for (SceneInfoList_i si = fScenes.begin(); si != fScenes.end(); ++si)
   {
      scenepad = dynamic_cast<TGLScenePad*>((*si)->GetScene());
      if (scenepad && scenepad->GetPad() == pad)
         break;
      scenepad = 0;
   }
   if (scenepad == 0)
   {
      scenepad = new TGLScenePad(pad);
      AddScene(scenepad);
   }

   scenepad->PadPaintFromViewer(this);

   PostSceneBuildSetup(fResetCamerasOnNextUpdate || fResetCamerasOnUpdate);
   fResetCamerasOnNextUpdate = kFALSE;

   RequestDraw();
}


/**************************************************************************/
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Force update of pad-scenes. Eventually this could be generalized
/// to all scene-types via a virtual function in TGLSceneBase.

void TGLViewer::UpdateScene(Bool_t redraw)
{
   // Cancel any pending redraw timer.
   fRedrawTimer->Stop();

   for (SceneInfoList_i si = fScenes.begin(); si != fScenes.end(); ++si)
   {
      TGLScenePad* scenepad = dynamic_cast<TGLScenePad*>((*si)->GetScene());
      if (scenepad)
         scenepad->PadPaintFromViewer(this);
   }

   PostSceneBuildSetup(fResetCamerasOnNextUpdate || fResetCamerasOnUpdate);
   fResetCamerasOnNextUpdate = kFALSE;

   if (redraw)
      RequestDraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Resets position/rotation of current camera to default values.

void TGLViewer::ResetCurrentCamera()
{
   MergeSceneBBoxes(fOverallBoundingBox);
   CurrentCamera().Setup(fOverallBoundingBox, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Setup cameras for current bounding box.

void TGLViewer::SetupCameras(Bool_t reset)
{
   if (IsLocked()) {
      Error("TGLViewer::SetupCameras", "expected kUnlocked, found %s", LockName(CurrentLock()));
      return;
   }

   // Setup cameras if scene box is not empty
   const TGLBoundingBox & box = fOverallBoundingBox;
   if (!box.IsEmpty()) {
      fPerspectiveCameraYOZ.Setup(box, reset);
      fPerspectiveCameraXOZ.Setup(box, reset);
      fPerspectiveCameraXOY.Setup(box, reset);
      fOrthoXOYCamera.Setup(box, reset);
      fOrthoXOZCamera.Setup(box, reset);
      fOrthoZOYCamera.Setup(box, reset);
      fOrthoXnOYCamera.Setup(box, reset);
      fOrthoXnOZCamera.Setup(box, reset);
      fOrthoZnOYCamera.Setup(box, reset);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform post scene-build setup.

void TGLViewer::PostSceneBuildSetup(Bool_t resetCameras)
{
   MergeSceneBBoxes(fOverallBoundingBox);
   SetupCameras(resetCameras);

   // Set default reference to center
   fReferencePos.Set(fOverallBoundingBox.Center());
   RefreshPadEditor(this);
}


/**************************************************************************/
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Initialise GL state.

void TGLViewer::InitGL()
{
   glEnable(GL_LIGHTING);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);
   glClearColor(0.f, 0.f, 0.f, 0.f);
   glClearDepth(1.0);
   glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);
   glMaterialf(GL_BACK, GL_SHININESS, 0.0);
   glPolygonMode(GL_FRONT, GL_FILL);
   glDisable(GL_BLEND);

   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
   Float_t lmodelAmb[] = {0.5f, 0.5f, 1.f, 1.f};
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodelAmb);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

   glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

   TGLUtil::CheckError("TGLViewer::InitGL");
}

////////////////////////////////////////////////////////////////////////////////
/// Post request for redraw of viewer at level of detail 'LOD'
/// Request is directed via cross thread gVirtualGL object.

void TGLViewer::RequestDraw(Short_t LODInput)
{
   fRedrawTimer->Stop();
   // Ignore request if GL window or context not yet available or shown.
   if ((!fGLWidget && fGLDevice == -1) || (fGLWidget && !fGLWidget->IsMapped()))
   {
      return;
   }

   // Take scene draw lock - to be revisited
   if ( ! TakeLock(kDrawLock)) {
      // If taking drawlock fails the previous draw is still in progress
      // set timer to do this one later
      if (gDebug>3) {
         Info("TGLViewer::RequestDraw", "viewer locked - requesting another draw.");
      }
      fRedrawTimer->RequestDraw(100, LODInput);
      return;
   }
   fLOD = LODInput;

   if (!gVirtualX->IsCmdThread())
      gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoDraw()", (ULong_t)this));
   else
      DoDraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Setup clip-object. Protected virtual method.

void TGLViewer::SetupClipObject()
{
   if (GetClipAutoUpdate())
   {
      fClipSet->SetupCurrentClip(fOverallBoundingBox);
   }
   else
   {
      fClipSet->SetupCurrentClipIfInvalid(fOverallBoundingBox);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Initialize objects that influence rendering.
/// Called before every render.

void TGLViewer::PreRender()
{
   fCamera = fCurrentCamera;
   fClip   = fClipSet->GetCurrentClip();
   if (fGLDevice != -1)
   {
      fRnrCtx->SetGLCtxIdentity(fGLCtxId);
      fGLCtxId->DeleteGLResources();
   }

   TGLUtil::SetPointSizeScale(fPointScale * fRnrCtx->GetRenderScale());
   TGLUtil::SetLineWidthScale(fLineScale  * fRnrCtx->GetRenderScale());

   if (fSmoothPoints) glEnable(GL_POINT_SMOOTH); else glDisable(GL_POINT_SMOOTH);
   if (fSmoothLines)  glEnable(GL_LINE_SMOOTH);  else glDisable(GL_LINE_SMOOTH);
   if (fSmoothPoints || fSmoothLines)
   {
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
   }
   else
   {
      glDisable(GL_BLEND);
   }

   TGLViewerBase::PreRender();

   // Setup lighting
   fLightSet->StdSetupLights(fOverallBoundingBox, *fCamera, fDebugMode);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal rendering, used by mono and stereo rendering.

void TGLViewer::Render()
{
   TGLViewerBase::Render();

   DrawGuides();
   RenderOverlay(TGLOverlayElement::kAllVisible, kFALSE);

   if ( ! fRnrCtx->Selection())
   {
      RenderSelectedForHighlight();
   }

   glClear(GL_DEPTH_BUFFER_BIT);
   DrawDebugInfo();
}

////////////////////////////////////////////////////////////////////////////////
/// Restore state set in PreRender().
/// Called after every render.

void TGLViewer::PostRender()
{
   TGLViewerBase::PostRender();

   TGLUtil::SetPointSizeScale(1);
   TGLUtil::SetLineWidthScale(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw out the viewer.

void TGLViewer::DoDraw(Bool_t swap_buffers)
{
   // Locking mainly for Win32 multi thread safety - but no harm in all using it
   // During normal draws a draw lock is taken in other thread (Win32) in RequestDraw()
   // to ensure thread safety. For PrintObjects repeated Draw() calls are made.
   // If no draw lock taken get one now.

   R__LOCKGUARD(gROOTMutex);

   fRedrawTimer->Stop();

   if (CurrentLock() != kDrawLock) {
      if ( ! TakeLock(kDrawLock)) {
         Error("TGLViewer::DoDraw", "viewer is %s", LockName(CurrentLock()));
         return;
      }
   }

   TUnlocker ulck(this);

   if (fGLDevice == -1 && (fViewport.Width() <= 1 || fViewport.Height() <= 1)) {
      if (gDebug > 2) {
         Info("TGLViewer::DoDraw()", "zero surface area, draw skipped.");
      }
      return;
   }

   if (fGLDevice != -1) {
      Int_t viewport[4] = {};
      gGLManager->ExtractViewport(fGLDevice, viewport);
      SetViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
   }

   TGLStopwatch timer;
   if (gDebug>2) {
      timer.Start();
   }

   // Setup scene draw time
   fRnrCtx->SetRenderTimeOut(fLOD == TGLRnrCtx::kLODHigh ?
                             fMaxSceneDrawTimeHQ :
                             fMaxSceneDrawTimeLQ);

   if (fStereo && fCurrentCamera->IsPerspective() && !fRnrCtx->GetGrabImage() &&
       !fIsPrinting)
   {
      DoDrawStereo(swap_buffers);
   }
   else
   {
      DoDrawMono(swap_buffers);
   }

   ReleaseLock(kDrawLock);

   if (gDebug>2) {
      Info("TGLViewer::DoDraw()", "Took %f msec", timer.End());
   }

   // Check if further redraws are needed and schedule them.

   if (CurrentCamera().UpdateInterest(kFALSE)) {
      // Reset major view-dependant cache.
      ResetSceneInfos();
      fRedrawTimer->RequestDraw(0, fLOD);
   }

   if (fLOD != TGLRnrCtx::kLODHigh &&
       (fDragAction < kDragCameraRotate || fDragAction > kDragCameraDolly))
   {
      // Request final draw pass.
      fRedrawTimer->RequestDraw(100, TGLRnrCtx::kLODHigh);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw out in monoscopic mode.

void TGLViewer::DoDrawMono(Bool_t swap_buffers)
{
   MakeCurrent();

   if (!fIsPrinting) PreDraw();
   PreRender();

   fRnrCtx->StartStopwatch();
   if (fFader < 1)
   {
      Render();
   }
   fRnrCtx->StopStopwatch();

   PostRender();

   if (fFader > 0)
   {
      FadeView(fFader);
   }

   PostDraw();

   if (swap_buffers)
   {
      SwapBuffers();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw out in stereoscopic mode.

void TGLViewer::DoDrawStereo(Bool_t swap_buffers)
{
   TGLPerspectiveCamera &c = *dynamic_cast<TGLPerspectiveCamera*>(fCurrentCamera);

   Float_t gl_near, gl_far, zero_p_dist;
   Float_t h_half, w_half;
   Float_t x_len_at_zero_parallax;
   Float_t stereo_offset;
   Float_t frustum_asym;

   MakeCurrent();

   // Draw left
   if (fStereoQuadBuf)
   {
      glDrawBuffer(GL_BACK_LEFT);
   }
   else
   {
      glScissor(0, 0, fViewport.Width(), fViewport.Height());
      glEnable(GL_SCISSOR_TEST);
   }
   PreDraw();
   PreRender();

   gl_near = c.GetNearClip();
   gl_far  = c.GetFarClip();
   zero_p_dist = gl_near + fStereoZeroParallax*(gl_far-gl_near);

   h_half = TMath::Tan(0.5*TMath::DegToRad()*c.GetFOV()) * gl_near;
   w_half = h_half * fViewport.Aspect();

   x_len_at_zero_parallax = 2.0f * w_half * zero_p_dist / gl_near;
   stereo_offset = 0.035f * x_len_at_zero_parallax * fStereoEyeOffsetFac;

   frustum_asym = stereo_offset * gl_near / zero_p_dist * fStereoFrustumAsymFac;

   TGLMatrix  abs_trans(c.RefCamBase());
   abs_trans *= c.RefCamTrans();
   TGLVector3 left_vec = abs_trans.GetBaseVec(2);

   glTranslatef(stereo_offset*left_vec[0], stereo_offset*left_vec[1], stereo_offset*left_vec[2]);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-w_half + frustum_asym, w_half + frustum_asym,
             -h_half, h_half, gl_near, gl_far);
   glMatrixMode(GL_MODELVIEW);

   fRnrCtx->StartStopwatch();
   if (fFader < 1)
   {
      Render();
   }
   fRnrCtx->StopStopwatch();

   PostRender();

   if (fFader > 0)
   {
      FadeView(fFader);
   }
   PostDraw();

   // Draw right
   if (fStereoQuadBuf)
   {
      glDrawBuffer(GL_BACK_RIGHT);
   }
   else
   {
      glScissor(fViewport.Width(), 0, fViewport.Width(), fViewport.Height());
   }
   PreDraw();
   PreRender();
   if ( ! fStereoQuadBuf)
   {
      glViewport(fViewport.Width(), 0, fViewport.Width(), fViewport.Height());
   }

   glTranslatef(-stereo_offset*left_vec[0], -stereo_offset*left_vec[1], -stereo_offset*left_vec[2]);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-w_half - frustum_asym, w_half - frustum_asym,
             -h_half, h_half, gl_near, gl_far);
   glMatrixMode(GL_MODELVIEW);

   fRnrCtx->StartStopwatch();
   if (fFader < 1)
   {
      Render();
   }
   fRnrCtx->StopStopwatch();

   PostRender();

   if (fFader > 0)
   {
      FadeView(fFader);
   }
   PostDraw();

   // End
   if (swap_buffers)
   {
      SwapBuffers();
   }

   if (fStereoQuadBuf)
   {
      glDrawBuffer(GL_BACK);
   }
   else
   {
      glDisable(GL_SCISSOR_TEST);
      glViewport(0, 0, fViewport.Width(), fViewport.Height());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save current image using the default file name which can be set
/// via SetPictureFileName() and defaults to "viewer.jpg".
/// Really useful for the files ending with 'gif+'.

Bool_t TGLViewer::SavePicture()
{
   return SavePicture(fPictureFileName);
}

////////////////////////////////////////////////////////////////////////////////
/// Save current image in various formats (gif, gif+, jpg, png, eps, pdf).
/// 'gif+' will append image to an existing file (animated gif).
/// 'eps' and 'pdf' do not fully support transparency and texturing.
/// The viewer window most be fully contained within the desktop but
/// can be covered by other windows.
/// Returns false if something obvious goes wrong, true otherwise.
///
/// The mage is saved using a frame-buffer object if the GL implementation
/// claims to support it -- this claim is not always true, especially when
/// running over ssh with drastically different GL implementations on the
/// client and server sides. Set this in .rootrc to enforce creation of
/// pictures using the back-buffer:
///   OpenGL.SavePicturesViaFBO: off

Bool_t TGLViewer::SavePicture(const TString &fileName)
{
   if (fileName.EndsWith(".eps"))
   {
      return TGLOutput::Capture(*this, TGLOutput::kEPS_BSP, fileName.Data());
   }
   else if (fileName.EndsWith(".pdf"))
   {
      return TGLOutput::Capture(*this, TGLOutput::kPDF_BSP, fileName.Data());
   }
   else
   {
      if (GLEW_EXT_framebuffer_object && gEnv->GetValue("OpenGL.SavePicturesViaFBO", 1))
      {
         return SavePictureUsingFBO(fileName, fViewport.Width(), fViewport.Height(), kFALSE);
      }
      else
      {
         return SavePictureUsingBB(fileName);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save current image in various formats (gif, gif+, jpg, png).
/// 'gif+' will append image to an existing file (animated gif).
/// Back-Buffer is used for capturing of the image.
/// The viewer window most be fully contained within the desktop but
/// can be covered by other windows.
/// Returns false if something obvious goes wrong, true otherwise.

Bool_t TGLViewer::SavePictureUsingBB(const TString &fileName)
{
   static const TString eh("TGLViewer::SavePictureUsingBB");

   if (! fileName.EndsWith(".gif") && ! fileName.Contains(".gif+") &&
       ! fileName.EndsWith(".jpg") && ! fileName.EndsWith(".png"))
   {
      Warning(eh, "file %s cannot be saved with this extension.", fileName.Data());
      return kFALSE;
   }

   if ( ! TakeLock(kDrawLock)) {
      Error(eh, "viewer locked - try later.");
      return kFALSE;
   }

   TUnlocker ulck(this);

   fLOD = TGLRnrCtx::kLODHigh;
   fRnrCtx->SetGrabImage(kTRUE);

   if (!gVirtualX->IsCmdThread())
      gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoDraw(kFALSE)", (ULong_t)this));
   else
      DoDraw(kFALSE);

   fRnrCtx->SetGrabImage(kFALSE);

   glReadBuffer(GL_BACK);

   UChar_t* xx = new UChar_t[4 * fViewport.Width() * fViewport.Height()];
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   glReadPixels(0, 0, fViewport.Width(), fViewport.Height(),
                GL_BGRA, GL_UNSIGNED_BYTE, xx);

   std::unique_ptr<TImage> image(TImage::Create());
   image->FromGLBuffer(xx, fViewport.Width(), fViewport.Height());
   image->WriteImage(fileName);

   delete [] xx;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save current image in various formats (gif, gif+, jpg, png).
/// 'gif+' will append image to an existing file (animated gif).
/// Frame-Buffer-Object is used for capturing of the image - OpenGL
/// 1.5 is required.
/// The viewer window does not have to be visible at all.
/// Returns false if something obvious goes wrong, true otherwise.
///
/// pixel_object_scale is used to scale (as much as possible) the
/// objects whose representation size is pixel based (point-sizes,
/// line-widths, bitmap/pixmap font-sizes).
/// If set to 0 (default) no scaling is applied.

Bool_t TGLViewer::SavePictureUsingFBO(const TString &fileName, Int_t w, Int_t h,
                                      Float_t pixel_object_scale)
{
   static const TString eh("TGLViewer::SavePictureUsingFBO");

   if (! fileName.EndsWith(".gif") && ! fileName.Contains(".gif+") &&
       ! fileName.EndsWith(".jpg") && ! fileName.EndsWith(".png"))
   {
      Warning(eh, "file %s cannot be saved with this extension.", fileName.Data());
      return kFALSE;
   }

   if ( ! TakeLock(kDrawLock)) {
      Error(eh, "viewer locked - try later.");
      return kFALSE;
   }

   TUnlocker ulck(this);

   MakeCurrent();

   TGLFBO *fbo = new TGLFBO();
   try
   {
      fbo->Init(w, h, fGLWidget->GetPixelFormat()->GetSamples());
   }
   catch (std::runtime_error& exc)
   {
      Error(eh, "%s",exc.what());
      if (gEnv->GetValue("OpenGL.SavePictureFallbackToBB", 1)) {
         Info(eh, "Falling back to saving image via back-buffer. Window must be fully visible.");
         if (w != fViewport.Width() || h != fViewport.Height())
            Warning(eh, "Back-buffer does not support image scaling, window size will be used.");
         return SavePictureUsingBB(fileName);
      } else {
         return kFALSE;
      }
   }

   TGLRect old_vp(fViewport);
   SetViewport(0, 0, w, h);

   Float_t old_scale = 1;
   if (pixel_object_scale != 0)
   {
      old_scale = fRnrCtx->GetRenderScale();
      fRnrCtx->SetRenderScale(old_scale * pixel_object_scale);
   }

   fbo->Bind();

   fLOD = TGLRnrCtx::kLODHigh;
   fRnrCtx->SetGrabImage(kTRUE);

   if (!gVirtualX->IsCmdThread())
      gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoDraw(kFALSE)", (ULong_t)this));
   else
      DoDraw(kFALSE);

   fRnrCtx->SetGrabImage(kFALSE);

   fbo->Unbind();

   fbo->SetAsReadBuffer();

   UChar_t* xx = new UChar_t[4 * fViewport.Width() * fViewport.Height()];
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   glReadPixels(0, 0, fViewport.Width(), fViewport.Height(),
                GL_BGRA, GL_UNSIGNED_BYTE, xx);

   std::unique_ptr<TImage> image(TImage::Create());
   image->FromGLBuffer(xx, fViewport.Width(), fViewport.Height());
   image->WriteImage(fileName);

   delete [] xx;

   delete fbo;

   if (pixel_object_scale != 0)
   {
      fRnrCtx->SetRenderScale(old_scale);
   }

   SetViewport(old_vp);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current image.
/// Back-Buffer is used for capturing of the image.
/// The viewer window most be fully contained within the desktop but
/// can be covered by other windows.

TImage* TGLViewer::GetPictureUsingBB()
{
    static const TString eh("TGLViewer::GetPictureUsingBB");

    if ( ! TakeLock(kDrawLock)) {
        Error(eh, "viewer locked - try later.");
        return NULL;
    }

    TUnlocker ulck(this);

    fLOD = TGLRnrCtx::kLODHigh;
    fRnrCtx->SetGrabImage(kTRUE);

    if (!gVirtualX->IsCmdThread())
        gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoDraw(kFALSE)", (ULong_t)this));
    else
        DoDraw(kFALSE);

    fRnrCtx->SetGrabImage(kFALSE);

    glReadBuffer(GL_BACK);

    UChar_t* xx = new UChar_t[4 * fViewport.Width() * fViewport.Height()];
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, fViewport.Width(), fViewport.Height(),
                 GL_BGRA, GL_UNSIGNED_BYTE, xx);

    TImage *image(TImage::Create());
    image->FromGLBuffer(xx, fViewport.Width(), fViewport.Height());

    delete [] xx;

    return image;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current image.
/// Frame-Buffer-Object is used for capturing of the image - OpenGL
/// 1.5 is required.
/// The viewer window does not have to be visible at all.
///
/// pixel_object_scale is used to scale (as much as possible) the
/// objects whose representation size is pixel based (point-sizes,
/// line-widths, bitmap/pixmap font-sizes).
/// If set to 0 (default) no scaling is applied.

TImage* TGLViewer::GetPictureUsingFBO(Int_t w, Int_t h,Float_t pixel_object_scale)
{
    static const TString eh("TGLViewer::GetPictureUsingFBO");

    if ( ! TakeLock(kDrawLock)) {
        Error(eh, "viewer locked - try later.");
        return NULL;
    }

    TUnlocker ulck(this);

    MakeCurrent();

    TGLFBO *fbo = new TGLFBO();
    try
    {
        fbo->Init(w, h, fGLWidget->GetPixelFormat()->GetSamples());
    }
    catch (std::runtime_error& exc)
    {
        Error(eh, "%s",exc.what());
        if (gEnv->GetValue("OpenGL.GetPictureFallbackToBB", 1)) {
            Info(eh, "Falling back to saving image via back-buffer. Window must be fully visible.");
            if (w != fViewport.Width() || h != fViewport.Height())
                Warning(eh, "Back-buffer does not support image scaling, window size will be used.");
            return GetPictureUsingBB();
        } else {
            return NULL;
        }
    }

    TGLRect old_vp(fViewport);
    SetViewport(0, 0, w, h);

    Float_t old_scale = 1;
    if (pixel_object_scale != 0)
    {
        old_scale = fRnrCtx->GetRenderScale();
        fRnrCtx->SetRenderScale(old_scale * pixel_object_scale);
    }

    fbo->Bind();

    fLOD = TGLRnrCtx::kLODHigh;
    fRnrCtx->SetGrabImage(kTRUE);

    if (!gVirtualX->IsCmdThread())
        gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoDraw(kFALSE)", (ULong_t)this));
    else
        DoDraw(kFALSE);

    fRnrCtx->SetGrabImage(kFALSE);

    fbo->Unbind();

    fbo->SetAsReadBuffer();

    UChar_t* xx = new UChar_t[4 * fViewport.Width() * fViewport.Height()];
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, fViewport.Width(), fViewport.Height(),
                 GL_BGRA, GL_UNSIGNED_BYTE, xx);

    TImage *image(TImage::Create());
    image->FromGLBuffer(xx, fViewport.Width(), fViewport.Height());

    delete [] xx;
    delete fbo;

    if (pixel_object_scale != 0)
    {
        fRnrCtx->SetRenderScale(old_scale);
    }

    SetViewport(old_vp);

    return image;
}


////////////////////////////////////////////////////////////////////////////////
/// Save picture with given width (height scaled proportionally).
/// If pixel_object_scale is true (default), the corresponding
/// scaling gets calculated from the current window size.

Bool_t TGLViewer::SavePictureWidth(const TString &fileName, Int_t width,
                                   Bool_t pixel_object_scale)
{
   Float_t scale  = Float_t(width) / fViewport.Width();
   Int_t   height = TMath::Nint(scale*fViewport.Height());

   return SavePictureUsingFBO(fileName, width, height, pixel_object_scale ? scale : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Save picture with given height (width scaled proportionally).
/// If pixel_object_scale is true (default), the corresponding
/// scaling gets calculated from the current window size.

Bool_t TGLViewer::SavePictureHeight(const TString &fileName, Int_t height,
                                    Bool_t pixel_object_scale)
{
   Float_t scale = Float_t(height) / fViewport.Height();
   Int_t   width = TMath::Nint(scale*fViewport.Width());

   return SavePictureUsingFBO(fileName, width, height, pixel_object_scale ? scale : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Save picture with given scale to current window size.
/// If pixel_object_scale is true (default), the same scaling is
/// used.

Bool_t TGLViewer::SavePictureScale (const TString &fileName, Float_t scale,
                                    Bool_t pixel_object_scale)
{
   Int_t w = TMath::Nint(scale*fViewport.Width());
   Int_t h = TMath::Nint(scale*fViewport.Height());

   return SavePictureUsingFBO(fileName, w, h, pixel_object_scale ? scale : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw reference marker and coordinate axes.

void TGLViewer::DrawGuides()
{
   Bool_t disabled = kFALSE;
   if (fReferenceOn)
   {
      glDisable(GL_DEPTH_TEST);
      TGLUtil::DrawReferenceMarker(*fCamera, fReferencePos);
      disabled = kTRUE;
   }
   if (fDrawCameraCenter)
   {
      glDisable(GL_DEPTH_TEST);
      Float_t radius = fCamera->ViewportDeltaToWorld(TGLVertex3(fCamera->GetCenterVec()), 3, 3).Mag();
      const UChar_t rgba[4] = { 0, 255, 255, 255 };
      TGLUtil::DrawSphere(fCamera->GetCenterVec(), radius, rgba);
      disabled = kTRUE;
   }
   if (fAxesDepthTest && disabled)
   {
      glEnable(GL_DEPTH_TEST);
      disabled = kFALSE;
   }
   else if (fAxesDepthTest == kFALSE && disabled == kFALSE)
   {
      glDisable(GL_DEPTH_TEST);
      disabled = kTRUE;
   }
   TGLUtil::DrawSimpleAxes(*fCamera, fOverallBoundingBox, fAxesType);
   if (disabled)
      glEnable(GL_DEPTH_TEST);
}

////////////////////////////////////////////////////////////////////////////////
/// If in debug mode draw camera aids and overall bounding box.

void TGLViewer::DrawDebugInfo()
{
   if (fDebugMode)
   {
      glDisable(GL_LIGHTING);
      CurrentCamera().DrawDebugAids();

      // Green scene bounding box
      glColor3d(0.0, 1.0, 0.0);
      fOverallBoundingBox.Draw();

      // Scene bounding box center sphere (green) and
      glDisable(GL_DEPTH_TEST);
      Double_t size = fOverallBoundingBox.Extents().Mag() / 200.0;
      TGLUtil::DrawSphere(TGLVertex3(0.0, 0.0, 0.0), size, TGLUtil::fgWhite);
      const TGLVertex3 & center = fOverallBoundingBox.Center();
      TGLUtil::DrawSphere(center, size, TGLUtil::fgGreen);
      glEnable(GL_DEPTH_TEST);

      glEnable(GL_LIGHTING);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform GL work which must be done before each draw.

void TGLViewer::PreDraw()
{
   InitGL();

   // For embedded gl clear color must be pad's background color.
   {
      Color_t ci = (fGLDevice != -1) ? gPad->GetFillColor() : fRnrCtx->ColorSet().Background().GetColorIndex();
      TColor *color = gROOT->GetColor(ci);
      Float_t rgb[3];
      if (color)
         color->GetRGB(rgb[0], rgb[1], rgb[2]);
      else
         rgb[0] = rgb[1] = rgb[2] = 0.0f;

      glClearColor(rgb[0], rgb[1], rgb[2], 0.0f);
   }

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

   TGLUtil::CheckError("TGLViewer::PreDraw");
}

////////////////////////////////////////////////////////////////////////////////
/// Perform GL work which must be done after each draw.

void TGLViewer::PostDraw()
{
   glFlush();
   TGLUtil::CheckError("TGLViewer::PostDraw");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a rectangle (background color and given alpha) across the
/// whole viewport.

void TGLViewer::FadeView(Float_t alpha)
{
   static const Float_t z = -1.0f;

   glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
   glMatrixMode(GL_MODELVIEW);  glPushMatrix(); glLoadIdentity();

   {
      TGLCapabilitySwitch blend(GL_BLEND,    kTRUE);
      TGLCapabilitySwitch light(GL_LIGHTING, kFALSE);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      TGLUtil::ColorAlpha(fRnrCtx->ColorSet().Background(), alpha);
      glBegin(GL_QUADS);
      glVertex3f(-1, -1, z);  glVertex3f( 1, -1, z);
      glVertex3f( 1,  1, z);  glVertex3f(-1,  1, z);
      glEnd();
   }

   glMatrixMode(GL_PROJECTION); glPopMatrix();
   glMatrixMode(GL_MODELVIEW);  glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Make GL context current

void TGLViewer::MakeCurrent() const
{
   if (fGLDevice == -1)
      fGLWidget->MakeCurrent();
   else
      gGLManager->MakeCurrent(fGLDevice);
}

////////////////////////////////////////////////////////////////////////////////
/// Swap GL buffers

void TGLViewer::SwapBuffers() const
{
   if ( ! IsDrawOrSelectLock()) {
      Error("TGLViewer::SwapBuffers", "viewer is %s", LockName(CurrentLock()));
   }
   if (fGLDevice == -1)
      fGLWidget->SwapBuffers();
   else {
      gGLManager->ReadGLBuffer(fGLDevice);
      gGLManager->Flush(fGLDevice);
      gGLManager->MarkForDirectCopy(fGLDevice, kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Post request for selection render pass viewer, picking objects
/// around the window point (x,y).

Bool_t TGLViewer::RequestSelect(Int_t x, Int_t y)
{
   // Take select lock on scene immediately we enter here - it is released
   // in the other (drawing) thread - see TGLViewer::DoSelect()

   if ( ! TakeLock(kSelectLock)) {
      return kFALSE;
   }

   if (!gVirtualX->IsCmdThread())
      return Bool_t(gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoSelect(%d, %d)", (ULong_t)this, x, y)));
   else
      return DoSelect(x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform GL selection, picking objects overlapping WINDOW
/// area described by 'rect'. Return kTRUE if selection should be
/// changed, kFALSE otherwise.
/// Select lock should already been taken in other thread in
/// TGLViewer::ReqSelect().

Bool_t TGLViewer::DoSelect(Int_t x, Int_t y)
{
   R__LOCKGUARD(gROOTMutex);

   if (CurrentLock() != kSelectLock) {
      Error("TGLViewer::DoSelect", "expected kSelectLock, found %s", LockName(CurrentLock()));
      return kFALSE;
   }

   TGLUtil::PointToViewport(x, y);

   TUnlocker ulck(this);

   MakeCurrent();

   fRnrCtx->BeginSelection(x, y, TGLUtil::GetPickingRadius());
   glRenderMode(GL_SELECT);

   PreRender();
   TGLViewerBase::Render();
   PostRender();

   Int_t nHits = glRenderMode(GL_RENDER);
   fRnrCtx->EndSelection(nHits);

   // Process selection.
   if (gDebug > 0) Info("TGLViewer::DoSelect", "Primary select nHits=%d.", nHits);

   if (nHits > 0)
   {
      Int_t idx = 0;
      if (FindClosestRecord(fSelRec, idx))
      {
         if (fSelRec.GetTransparent() && fRnrCtx->SelectTransparents() != TGLRnrCtx::kIfClosest)
         {
            TGLSelectRecord opaque;
            if (FindClosestOpaqueRecord(opaque, ++idx))
               fSelRec = opaque;
            else if (fRnrCtx->SelectTransparents() == TGLRnrCtx::kNever)
               fSelRec.Reset();
         }
         if (gDebug > 1) fSelRec.Print();
      }
   } else {
      fSelRec.Reset();
   }

   ReleaseLock(kSelectLock);
   return ! TGLSelectRecord::AreSameSelectionWise(fSelRec, fCurrentSelRec);
}

////////////////////////////////////////////////////////////////////////////////
/// Request secondary select.

Bool_t TGLViewer::RequestSecondarySelect(Int_t x, Int_t y)
{
   if ( ! TakeLock(kSelectLock)) {
      return kFALSE;
   }

   if (!gVirtualX->IsCmdThread())
      return Bool_t(gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoSecondarySelect(%d, %d)", (ULong_t)this, x, y)));
   else
      return DoSecondarySelect(x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Secondary selection.

Bool_t TGLViewer::DoSecondarySelect(Int_t x, Int_t y)
{
   R__LOCKGUARD(gROOTMutex);

   if (CurrentLock() != kSelectLock) {
      Error("TGLViewer::DoSecondarySelect", "expected kSelectLock, found %s", LockName(CurrentLock()));
      return kFALSE;
   }

   TGLUtil::PointToViewport(x, y);

   TUnlocker ulck(this);

   if (! fSelRec.GetSceneInfo() || ! fSelRec.GetPhysShape() ||
       ! fSelRec.GetLogShape()->SupportsSecondarySelect())
   {
      if (gDebug > 0)
         Info("TGLViewer::SecondarySelect", "Skipping secondary selection "
              "(sinfo=0x%lx, pshape=0x%lx).\n",
              (Long_t)fSelRec.GetSceneInfo(), (Long_t)fSelRec.GetPhysShape());
      fSecSelRec.Reset();
      return kFALSE;
   }

   MakeCurrent();

   TGLSceneInfo*    sinfo = fSelRec.GetSceneInfo();
   TGLSceneBase*    scene = sinfo->GetScene();
   TGLPhysicalShape* pshp = fSelRec.GetPhysShape();

   SceneInfoList_t foo;
   foo.push_back(sinfo);
   fScenes.swap(foo);
   fRnrCtx->BeginSelection(x, y, TGLUtil::GetPickingRadius());
   fRnrCtx->SetSecSelection(kTRUE);
   glRenderMode(GL_SELECT);

   PreRender();
   fRnrCtx->SetSceneInfo(sinfo);
   scene->PreRender(*fRnrCtx);
   fRnrCtx->SetDrawPass(TGLRnrCtx::kPassFill);
   fRnrCtx->SetShapeLOD(TGLRnrCtx::kLODHigh);
   glPushName(pshp->ID());
   // !!! Hack: does not use clipping and proper draw-pass settings.
   pshp->Draw(*fRnrCtx);
   glPopName();
   scene->PostRender(*fRnrCtx);
   fRnrCtx->SetSceneInfo(0);
   PostRender();

   Int_t nSecHits = glRenderMode(GL_RENDER);
   fRnrCtx->EndSelection(nSecHits);
   fScenes.swap(foo);

   if (gDebug > 0) Info("TGLViewer::DoSelect", "Secondary select nSecHits=%d.", nSecHits);

   ReleaseLock(kSelectLock);

   if (nSecHits > 0)
   {
      fSecSelRec = fSelRec;
      fSecSelRec.SetRawOnly(fRnrCtx->GetSelectBuffer()->RawRecord(0));
      if (gDebug > 1) fSecSelRec.Print();
      return kTRUE;
   } else {
      fSecSelRec.Reset();
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process result from last selection (in fSelRec) and
/// extract a new current selection from it.
/// Here we only use physical shape.

void TGLViewer::ApplySelection()
{
   fCurrentSelRec = fSelRec;

   TGLPhysicalShape *selPhys = fSelRec.GetPhysShape();
   fSelectedPShapeRef->SetPShape(selPhys);

   // Inform external client selection has been modified.
   SelectionChanged();

   RequestDraw(TGLRnrCtx::kLODHigh);
}

////////////////////////////////////////////////////////////////////////////////
/// Post request for secondary selection rendering of selected object
/// around the window point (x,y).

Bool_t TGLViewer::RequestOverlaySelect(Int_t x, Int_t y)
{
   // Take select lock on viewer immediately - it is released
   // in the other (drawing) thread - see TGLViewer::DoSecondarySelect().

   if ( ! TakeLock(kSelectLock)) {
      return kFALSE;
   }

   if (!gVirtualX->IsCmdThread())
      return Bool_t(gROOT->ProcessLineFast(Form("((TGLViewer *)0x%lx)->DoOverlaySelect(%d, %d)", (ULong_t)this, x, y)));
   else
      return DoOverlaySelect(x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform GL selection, picking overlay objects only.
/// Return TRUE if the selected overlay-element has changed.

Bool_t TGLViewer::DoOverlaySelect(Int_t x, Int_t y)
{
   R__LOCKGUARD(gROOTMutex);

   if (CurrentLock() != kSelectLock) {
      Error("TGLViewer::DoOverlaySelect", "expected kSelectLock, found %s", LockName(CurrentLock()));
      return kFALSE;
   }

   TGLUtil::PointToViewport(x, y);

   TUnlocker ulck(this);

   MakeCurrent();

   fRnrCtx->BeginSelection(x, y, TGLUtil::GetPickingRadius());
   glRenderMode(GL_SELECT);

   PreRenderOverlaySelection();
   RenderOverlay(TGLOverlayElement::kActive, kTRUE);
   PostRenderOverlaySelection();

   Int_t nHits = glRenderMode(GL_RENDER);
   fRnrCtx->EndSelection(nHits);

   // Process overlay selection.
   TGLOverlayElement * selElm = 0;
   if (nHits > 0)
   {
      Int_t idx = 0;
      while (idx < nHits && FindClosestOverlayRecord(fOvlSelRec, idx))
      {
         TGLOverlayElement* el = fOvlSelRec.GetOvlElement();
         if (el == fCurrentOvlElm)
         {
            if (el->MouseStillInside(fOvlSelRec))
            {
               selElm = el;
               break;
            }
         }
         else if (el->MouseEnter(fOvlSelRec))
         {
            selElm = el;
            break;
         }
         ++idx;
      }
   }
   else
   {
      fOvlSelRec.Reset();
   }

   ReleaseLock(kSelectLock);

   if (fCurrentOvlElm != selElm)
   {
      if (fCurrentOvlElm) fCurrentOvlElm->MouseLeave();
      fCurrentOvlElm = selElm;
      return kTRUE;
   }
   else
   {
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make one fading step and request redraw.

void TGLFaderHelper::MakeFadeStep()
{
   Float_t fade = fViewer->GetFader();

   if (fade == fFadeTarget) {
      delete this; return;
   }
   if (TMath::Abs(fFadeTarget - fade) < 1e-3) {
      fViewer->SetFader(fFadeTarget);
      fViewer->RequestDraw(TGLRnrCtx::kLODHigh);
      delete this;
      return;
   }

   Float_t dt = fTime/fNSteps;
   Float_t df = (fFadeTarget - fade)/fNSteps;
   fViewer->SetFader(fade + df);
   fViewer->RequestDraw(TGLRnrCtx::kLODHigh);
   fTime -= dt; --fNSteps;
   TTimer::SingleShot(TMath::CeilNint(1000*dt),
                      "TGLFaderHelper", this, "MakeFadeStep()");
}

////////////////////////////////////////////////////////////////////////////////
/// Animate fading from current value to fade over given time (sec)
/// and number of steps.

void TGLViewer::AutoFade(Float_t fade, Float_t time, Int_t steps)
{
   TGLFaderHelper* fh = new TGLFaderHelper(this, fade, time, steps);
   fh->MakeFadeStep();
}

////////////////////////////////////////////////////////////////////////////////
/// Use the dark color-set.

void TGLViewer::UseDarkColorSet()
{
   fRnrCtx->ChangeBaseColorSet(&fDarkColorSet);
   RefreshPadEditor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Use the light color-set.

void TGLViewer::UseLightColorSet()
{
   fRnrCtx->ChangeBaseColorSet(&fLightColorSet);
   RefreshPadEditor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Switch between dark and light colorsets.

void TGLViewer::SwitchColorSet()
{
   if (IsUsingDefaultColorSet())
   {
      Info("SwitchColorSet()", "Global color-set is in use, switch not supported.");
      return;
   }

   if (fRnrCtx->GetBaseColorSet() == &fLightColorSet)
      UseDarkColorSet();
   else
      UseLightColorSet();
}

////////////////////////////////////////////////////////////////////////////////
/// Set usage of the default color set.

void TGLViewer::UseDefaultColorSet(Bool_t x)
{
   if (x)
      fRnrCtx->ChangeBaseColorSet(&fgDefaultColorSet);
   else
      fRnrCtx->ChangeBaseColorSet(&fDarkColorSet);
   RefreshPadEditor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the viewer is using the default color set.
/// If yes, some operations might be disabled.

Bool_t TGLViewer::IsUsingDefaultColorSet() const
{
   return fRnrCtx->GetBaseColorSet() == &fgDefaultColorSet;
}

////////////////////////////////////////////////////////////////////////////////
/// Set background method.
/// Deprecated method - set background color in the color-set.

void TGLViewer::SetClearColor(Color_t col)
{
   fRnrCtx->GetBaseColorSet()->Background().SetColor(col);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns reference to the default color-set.
/// Static function.

TGLColorSet& TGLViewer::GetDefaultColorSet()
{
   return fgDefaultColorSet;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets static flag that determines if new viewers should use the
/// default color-set.
/// This is false at startup.

void TGLViewer::UseDefaultColorSetForNewViewers(Bool_t x)
{
   fgUseDefaultColorSetForNewViewers = x;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the value of the static flag that determines if new
/// viewers should use the default color-set.
/// This is false at startup.

Bool_t TGLViewer::IsUsingDefaultColorSetForNewViewers()
{
   return fgUseDefaultColorSetForNewViewers;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if current color set is dark.

Bool_t TGLViewer::IsColorSetDark() const
{
   return fRnrCtx->GetBaseColorSet() == &fDarkColorSet;
}

/**************************************************************************/
// Viewport
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Set viewer viewport (window area) with bottom/left at (x,y), with
/// dimensions 'width'/'height'

void TGLViewer::SetViewport(Int_t x, Int_t y, Int_t width, Int_t height)
{
   if (fStereo && ! fStereoQuadBuf) width /= 2;

   // Only process if changed
   if (fViewport.X() == x && fViewport.Y() == y &&
       fViewport.Width() == width && fViewport.Height() == height) {
      return;
   }

   fViewport.Set(x, y, width, height);
   fCurrentCamera->SetViewport(fViewport);

   if (gDebug > 2) {
      Info("TGLViewer::SetViewport", "updated - corner %d,%d dimensions %d,%d", x, y, width, height);
   }
}

void TGLViewer::SetViewport(const TGLRect& vp)
{
   // Set viewer viewport from TGLRect.

   SetViewport(vp.X(), vp.Y(), vp.Width(), vp.Height());
}

/**************************************************************************/
// Camera methods
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Return camera reference by type.

TGLCamera& TGLViewer::RefCamera(ECameraType cameraType)
{
   // TODO: Move these into a vector!
   switch(cameraType) {
      case kCameraPerspXOZ:
         return fPerspectiveCameraXOZ;
      case kCameraPerspYOZ:
         return fPerspectiveCameraYOZ;
      case kCameraPerspXOY:
         return fPerspectiveCameraXOY;
      case kCameraOrthoXOY:
         return fOrthoXOYCamera;
      case kCameraOrthoXOZ:
         return fOrthoXOZCamera;
      case kCameraOrthoZOY:
         return fOrthoZOXCamera;
      case kCameraOrthoZOX:
         return fOrthoZOYCamera;
      case kCameraOrthoXnOY:
         return fOrthoXnOYCamera;
      case kCameraOrthoXnOZ:
         return fOrthoXnOZCamera;
      case kCameraOrthoZnOY:
         return fOrthoZnOYCamera;
      case kCameraOrthoZnOX:
         return fOrthoZnOXCamera;
      default:
         Error("TGLViewer::SetCurrentCamera", "invalid camera type");
         return *fCurrentCamera;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set current active camera - 'cameraType' one of:
///   kCameraPerspX,    kCameraPerspY,    kCameraPerspZ,
///   kCameraOrthoXOY,  kCameraOrthoXOZ,  kCameraOrthoZOY,
///   kCameraOrthoXnOY, kCameraOrthoXnOZ, kCameraOrthoZnOY

void TGLViewer::SetCurrentCamera(ECameraType cameraType)
{
   if (IsLocked()) {
      Error("TGLViewer::SetCurrentCamera", "expected kUnlocked, found %s", LockName(CurrentLock()));
      return;
   }

   // TODO: Move these into a vector!
   TGLCamera *prev = fCurrentCamera;
   switch (cameraType)
   {
      case kCameraPerspXOZ: {
         fCurrentCamera = &fPerspectiveCameraXOZ;
         break;
      }
      case kCameraPerspYOZ: {
         fCurrentCamera = &fPerspectiveCameraYOZ;
         break;
      }
      case kCameraPerspXOY: {
         fCurrentCamera = &fPerspectiveCameraXOY;
         break;
      }
      case kCameraOrthoXOY: {
         fCurrentCamera = &fOrthoXOYCamera;
         break;
      }
      case kCameraOrthoXOZ: {
         fCurrentCamera = &fOrthoXOZCamera;
         break;
      }
      case kCameraOrthoZOY: {
         fCurrentCamera = &fOrthoZOYCamera;
         break;
      }
      case kCameraOrthoZOX: {
         fCurrentCamera = &fOrthoZOXCamera;
         break;
      }
      case kCameraOrthoXnOY: {
         fCurrentCamera = &fOrthoXnOYCamera;
         break;
      }
      case kCameraOrthoXnOZ: {
         fCurrentCamera = &fOrthoXnOZCamera;
         break;
      }
      case kCameraOrthoZnOY: {
         fCurrentCamera = &fOrthoZnOYCamera;
         break;
      }
      case kCameraOrthoZnOX: {
         fCurrentCamera = &fOrthoZnOXCamera;
         break;
      }
      default: {
         Error("TGLViewer::SetCurrentCamera", "invalid camera type");
         break;
      }
   }

   if (fCurrentCamera != prev)
   {
      // Ensure any viewport has been propagated to the current camera
      fCurrentCamera->SetViewport(fViewport);
      RefreshPadEditor(this);

      if (fAutoRotator)
      {
         if (fAutoRotator->IsRunning())
         {
            fAutoRotator->Stop();
         }
         else
         {
            if (fAutoRotator->GetCamera() == fCurrentCamera)
            {
               fAutoRotator->Start();
            }
         }
      }

      RequestDraw(TGLRnrCtx::kLODHigh);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set an orthographic camera to supplied configuration - note this
/// does not need to be the current camera - though you will not see
/// the effect if it is not.
///
/// 'camera' defines the ortho camera - one of kCameraOrthoXOY / XOZ / ZOY
/// 'left' / 'right' / 'top' / 'bottom' define the WORLD coordinates which
/// correspond with the left/right/top/bottom positions on the GL viewer viewport
/// E.g. for kCameraOrthoXOY camera left/right are X world coords,
/// top/bottom are Y world coords
/// As this is an orthographic camera the other axis (in eye direction) is
/// no relevant. The near/far clip planes are set automatically based in scene
/// contents

void TGLViewer::SetOrthoCamera(ECameraType camera,
                               Double_t zoom, Double_t dolly,
                               Double_t center[3],
                               Double_t hRotate, Double_t vRotate)
{
   // TODO: Move these into a vector!
   switch(camera) {
      case kCameraOrthoXOY: {
         fOrthoXOYCamera.Configure(zoom, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fOrthoXOYCamera) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      case kCameraOrthoXOZ: {
         fOrthoXOZCamera.Configure(zoom, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fOrthoXOZCamera) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      case kCameraOrthoZOY: {
         fOrthoZOYCamera.Configure(zoom, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fOrthoZOYCamera) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      case kCameraOrthoZOX: {
         fOrthoZOXCamera.Configure(zoom, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fOrthoZOXCamera) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      default: {
         Error("TGLViewer::SetOrthoCamera", "invalid camera type");
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set a perspective camera to supplied configuration - note this
/// does not need to be the current camera - though you will not see
/// the effect if it is not.
///
///  - 'camera' defines the persp camera - one of kCameraPerspXOZ, kCameraPerspYOZ, kCameraPerspXOY
///  - 'fov' - field of view (lens angle) in degrees (clamped to 0.1 - 170.0)
///  - 'dolly' - distance from 'center'
///  - 'center' - world position from which dolly/hRotate/vRotate are measured
///                camera rotates round this, always facing in (in center of viewport)
///  - 'hRotate' - horizontal rotation from initial configuration in degrees
///  - 'hRotate' - vertical rotation from initial configuration in degrees

void TGLViewer::SetPerspectiveCamera(ECameraType camera,
                                     Double_t fov, Double_t dolly,
                                     Double_t center[3],
                                     Double_t hRotate, Double_t vRotate)
{
   // TODO: Move these into a vector!
   switch(camera) {
      case kCameraPerspXOZ: {
         fPerspectiveCameraXOZ.Configure(fov, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fPerspectiveCameraXOZ) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      case kCameraPerspYOZ: {
         fPerspectiveCameraYOZ.Configure(fov, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fPerspectiveCameraYOZ) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      case kCameraPerspXOY: {
         fPerspectiveCameraXOY.Configure(fov, dolly, center, hRotate, vRotate);
         if (fCurrentCamera == &fPerspectiveCameraXOY) {
            RequestDraw(TGLRnrCtx::kLODHigh);
         }
         break;
      }
      default: {
         Error("TGLViewer::SetPerspectiveCamera", "invalid camera type");
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change base-vectors defining the camera-base transformation of current
/// camera. hAxis and vAxis are the default directions for forward
/// (inverted) and upwards.

void TGLViewer::ReinitializeCurrentCamera(const TGLVector3& hAxis, const TGLVector3& vAxis, Bool_t redraw)
{
   TGLMatrix& cb = fCurrentCamera->RefCamBase();
   cb.Set(cb.GetTranslation(), vAxis, hAxis);
   fCurrentCamera->Setup(fOverallBoundingBox, kTRUE);
   if (redraw)
      RequestDraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the auto-rotator for this viewer.

TGLAutoRotator* TGLViewer::GetAutoRotator()
{
   if (fAutoRotator == 0)
      fAutoRotator = new TGLAutoRotator(this);
   return fAutoRotator;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the auto-rotator for this viewer. The old rotator is deleted.

void TGLViewer::SetAutoRotator(TGLAutoRotator* ar)
{
   delete fAutoRotator;
   fAutoRotator = ar;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable stereo rendering.
/// If quad_buf is true rendering is done into separate left and right GL
/// buffers. This requires hardware support. Otherwise left and right images
/// get rendered into left and right half of the window.
/// Note that mouse highlighting and selection will not work exactly right
/// as image for each eye gets slightly shifted and there are two different
/// directions through the mouse pointer, one for each eye.

void TGLViewer::SetStereo(Bool_t stereo, Bool_t quad_buf)
{
   if (stereo != fStereo)
   {
      fStereo = stereo;
      fStereoQuadBuf = quad_buf;
      if (fStereo)
         SetViewport(fViewport.X(), fViewport.Y(),   fViewport.Width(), fViewport.Height());
      else
         SetViewport(fViewport.X(), fViewport.Y(), 2*fViewport.Width(), fViewport.Height());
   }
   RequestDraw(TGLRnrCtx::kLODHigh);
}

/**************************************************************************/
// Guide methods
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Fetch the state of guides (axes & reference markers) into arguments

void TGLViewer::GetGuideState(Int_t & axesType, Bool_t & axesDepthTest, Bool_t & referenceOn, Double_t referencePos[3]) const
{
   axesType = fAxesType;
   axesDepthTest = fAxesDepthTest;

   referenceOn = fReferenceOn;
   referencePos[0] = fReferencePos.X();
   referencePos[1] = fReferencePos.Y();
   referencePos[2] = fReferencePos.Z();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the state of guides (axes & reference markers) from arguments.

void TGLViewer::SetGuideState(Int_t axesType, Bool_t axesDepthTest, Bool_t referenceOn, const Double_t referencePos[3])
{
   fAxesType    = axesType;
   fAxesDepthTest = axesDepthTest;
   fReferenceOn = referenceOn;
   if (referencePos)
      fReferencePos.Set(referencePos[0], referencePos[1], referencePos[2]);
   if (fGLDevice != -1)
      gGLManager->MarkForDirectCopy(fGLDevice, kTRUE);
   RequestDraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw camera look at and rotation point.

void TGLViewer::SetDrawCameraCenter(Bool_t x)
{
   fDrawCameraCenter = x;
   RequestDraw();
}

// Selected physical
////////////////////////////////////////////////////////////////////////////////
/// Return selected physical shape.

const TGLPhysicalShape * TGLViewer::GetSelected() const
{
   return fSelectedPShapeRef->GetPShape();
}

/**************************************************************************/
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Emit MouseOver signal.

void TGLViewer::MouseOver(TGLPhysicalShape *shape)
{
   Emit("MouseOver(TGLPhysicalShape*)", (Long_t)shape);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit MouseOver signal.

void TGLViewer::MouseOver(TGLPhysicalShape *shape, UInt_t state)
{
   Long_t args[2];
   args[0] = (Long_t)shape;
   args[1] = state;
   Emit("MouseOver(TGLPhysicalShape*,UInt_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit MouseOver signal.

void TGLViewer::MouseOver(TObject *obj, UInt_t state)
{
   Long_t args[2];
   args[0] = (Long_t)obj;
   args[1] = state;
   Emit("MouseOver(TObject*,UInt_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit MouseOver signal.

void TGLViewer::ReMouseOver(TObject *obj, UInt_t state)
{
   Long_t args[2];
   args[0] = (Long_t)obj;
   args[1] = state;
   Emit("ReMouseOver(TObject*,UInt_t)", args);
}


////////////////////////////////////////////////////////////////////////////////
/// Emit UnMouseOver signal.

void TGLViewer::UnMouseOver(TObject *obj, UInt_t state)
{
   Long_t args[2];
   args[0] = (Long_t)obj;
   args[1] = state;
   Emit("UnMouseOver(TObject*,UInt_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Clicked signal.

void TGLViewer::Clicked(TObject *obj)
{
   Emit("Clicked(TObject*)", (Long_t)obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Clicked signal with button id and modifier state.

void TGLViewer::Clicked(TObject *obj, UInt_t button, UInt_t state)
{
   Long_t args[3];
   args[0] = (Long_t)obj;
   args[1] = button;
   args[2] = state;
   Emit("Clicked(TObject*,UInt_t,UInt_t)", args);
}


////////////////////////////////////////////////////////////////////////////////
/// Emit ReClicked signal with button id and modifier state.

void TGLViewer::ReClicked(TObject *obj, UInt_t button, UInt_t state)
{
   Long_t args[3];
   args[0] = (Long_t)obj;
   args[1] = button;
   args[2] = state;
   Emit("ReClicked(TObject*,UInt_t,UInt_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit UnClicked signal with button id and modifier state.

void TGLViewer::UnClicked(TObject *obj, UInt_t button, UInt_t state)
{
   Long_t args[3];
   args[0] = (Long_t)obj;
   args[1] = button;
   args[2] = state;
   Emit("UnClicked(TObject*,UInt_t,UInt_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit MouseIdle signal.

void TGLViewer::MouseIdle(TGLPhysicalShape *shape, UInt_t posx, UInt_t posy)
{
   Long_t args[3];
   static UInt_t oldx = 0, oldy = 0;

   if (oldx != posx || oldy != posy) {
      args[0] = (Long_t)shape;
      args[1] = posx;
      args[2] = posy;
      Emit("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)", args);
      oldx = posx;
      oldy = posy;
   }
}

/**************************************************************************/
/**************************************************************************/
////////////////////////////////////////////////////////////////////////////////
/// Calculate and return pixel distance to nearest viewer object from
/// window location px, py
/// This is provided for use when embedding GL viewer into pad

Int_t TGLViewer::DistancetoPrimitive(Int_t /*px*/, Int_t /*py*/)
{
   // Can't track the indvidual objects in rollover. Just set the viewer as the
   // selected object, and return 0 (object identified) so we receive ExecuteEvent calls
   gPad->SetSelected(this);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Process event of type 'event' - one of EEventType types,
/// occurring at window location px, py
/// This is provided for use when embedding GL viewer into pad

void TGLViewer::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (fEventHandler)
      return fEventHandler->ExecuteEvent(event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Pass viewer for print capture by TGLOutput.

void TGLViewer::PrintObjects()
{
   TGLOutput::Capture(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Update GUI components for embedded viewer selection change.

void TGLViewer::SelectionChanged()
{
   if (!fGedEditor)
      return;

   TGLPhysicalShape *selected = const_cast<TGLPhysicalShape*>(GetSelected());

   if (selected) {
      fPShapeWrap->fPShape = selected;
      fGedEditor->SetModel(fPad, fPShapeWrap, kButton1Down);
   } else {
      fPShapeWrap->fPShape = 0;
      fGedEditor->SetModel(fPad, this, kButton1Down);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// An overlay operation can result in change to an object.
/// Refresh geditor.

void TGLViewer::OverlayDragFinished()
{
   if (fGedEditor)
   {
      fGedEditor->SetModel(fPad, fGedEditor->GetModel(), kButton1Down);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update GED editor if it is set.

void TGLViewer::RefreshPadEditor(TObject* obj)
{
   if (fGedEditor && (obj == 0 || fGedEditor->GetModel() == obj))
   {
      fGedEditor->SetModel(fPad, fGedEditor->GetModel(), kButton1Down);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the event-handler. The event-handler is owned by the viewer.
/// If GLWidget is set, the handler is propagated to it.
///
/// If called with handler=0, the current handler will be deleted
/// (also from TGLWidget).

void TGLViewer::SetEventHandler(TGEventHandler *handler)
{
   if (fEventHandler)
      delete fEventHandler;

   fEventHandler = handler;
   if (fGLWidget)
      fGLWidget->SetEventHandler(fEventHandler);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove overlay element.

void  TGLViewer::RemoveOverlayElement(TGLOverlayElement* el)
{
   if (el == fCurrentOvlElm)
   {
      fCurrentOvlElm = 0;
   }
   TGLViewerBase::RemoveOverlayElement(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset current overlay-element to zero, eventually notifying the
/// old one that the mouse has left.
/// Usually called when mouse leaves the window.

void TGLViewer::ClearCurrentOvlElm()
{
   if (fCurrentOvlElm)
   {
      fCurrentOvlElm->MouseLeave();
      fCurrentOvlElm = 0;
      RequestDraw();
   }
}
