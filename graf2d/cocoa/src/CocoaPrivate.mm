// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov   29/11/2011

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#define DEBUG_ROOT_COCOA

//#define NDEBUG

#ifdef DEBUG_ROOT_COCOA
#include <algorithm>
#endif

#include <stdexcept>
#include <cassert>

#include <OpenGL/OpenGL.h>
#include <Cocoa/Cocoa.h>

#include "ROOTApplicationDelegate.h"
#include "ROOTOpenGLView.h"
#include "CocoaPrivate.h"
#include "QuartzWindow.h"
#include "CocoaUtils.h"

namespace ROOT {
namespace MacOSX {
namespace Details {

//______________________________________________________________________________
CocoaPrivate::CocoaPrivate()
               : fCurrentDrawableID(GetRootWindowID() + 1), //Any real window has id > rootID.
                                                            //0 is also used by some X11 functions as None.
                 fFreeGLContextID(1),
                 fApplicationDelegate([[ROOTApplicationDelegate alloc] init])
{
   //Init NSApplication, if it was not done yet.
   Util::AutoreleasePool pool;
   [NSApplication sharedApplication];
}

//______________________________________________________________________________
CocoaPrivate::~CocoaPrivate()
{
}

//______________________________________________________________________________
Window_t CocoaPrivate::GetRootWindowID()const
{
   //First I had root ID == 0, but this is None in X11 and
   //it can be used by ROOT, for example, I had trouble with
   //gClient able to found TGWindow for None - crash!
   return 1;
}

//______________________________________________________________________________
bool CocoaPrivate::IsRootWindow(Window_t windowID)const
{
   return windowID == GetRootWindowID();
}

//______________________________________________________________________________
Drawable_t CocoaPrivate::RegisterDrawable(NSObject *nsObj)
{
   //Return integer identifier for a new "drawable" (like in X11)

   if (fCurrentDrawableID == 999)//I have to skip this, many thanks to ROOT who uses 999 as "all windows".
      ++fCurrentDrawableID;

   Drawable_t newID = fCurrentDrawableID;

   if (fFreeDrawableIDs.size()) {
      newID = fFreeDrawableIDs.back();
      fFreeDrawableIDs.pop_back();
   } else
      fCurrentDrawableID++;

   assert(fDrawables.find(newID) == fDrawables.end() && "RegisterDrawable, id for new drawable is still in use");
   fDrawables[newID] = nsObj;

   return newID;
}

//______________________________________________________________________________
NSObject<X11Drawable> *CocoaPrivate::GetDrawable(Drawable_t drawableID)const
{
   const_drawable_iterator drawableIter = fDrawables.find(drawableID);

#ifdef DEBUG_ROOT_COCOA
   if (drawableIter == fDrawables.end()) {
      NSLog(@"Fatal error: requested non-existing drawable %lu", drawableID);
      //We do not care about efficiency, ROOT's gonna die on assert :)
      std::vector<Drawable_t>::const_iterator deletedDrawable = std::find(fFreeDrawableIDs.begin(), fFreeDrawableIDs.end(), drawableID);
      if (deletedDrawable != fFreeDrawableIDs.end()) {
         NSLog(@"This drawable was deleted already");
      } else {
         NSLog(@"This drawable not found among allocated/deleted drawables");
      }
   }
#endif
   assert(drawableIter != fDrawables.end() && "GetDrawable, non-existing drawable requested");
   return drawableIter->second.Get();
}

//______________________________________________________________________________
NSObject<X11Window> *CocoaPrivate::GetWindow(Window_t windowID)const
{
   const_drawable_iterator winIter = fDrawables.find(windowID);
#ifdef DEBUG_ROOT_COCOA
   if (winIter == fDrawables.end()) {
      NSLog(@"Fatal error: requested non-existing drawable %lu", windowID);
      //We do not care about efficiency, ROOT's gonna die on assert :)
      std::vector<Drawable_t>::const_iterator deletedDrawable = std::find(fFreeDrawableIDs.begin(), fFreeDrawableIDs.end(), windowID);
      if (deletedDrawable != fFreeDrawableIDs.end()) {
         NSLog(@"This window was deleted already");
      } else {
         NSLog(@"This window not found among allocated/deleted drawables");
      }
      return 0;
   }
#endif
   assert(winIter != fDrawables.end() && "GetWindow, non-existing window requested");
   return (NSObject<X11Window> *)winIter->second.Get();
}

//______________________________________________________________________________
void CocoaPrivate::DeleteDrawable(Drawable_t drawableID)
{
   drawable_iterator drawableIter = fDrawables.find(drawableID);
   assert(drawableIter != fDrawables.end() && "DeleteDrawable, non existing drawableID");

   NSObject<X11Drawable> * const base = drawableIter->second.Get();
   if ([base isKindOfClass : [QuartzView class]]) {
      [(QuartzView *)base removeFromSuperview];
      ((QuartzView *)base).fParentView = nil;
   } else if ([base isKindOfClass : [QuartzWindow class]]) {
      QuartzWindow *qw = (QuartzWindow *)base;
      qw.fContentView.fParentView = nil;
      [qw.fContentView removeFromSuperview];
      qw.contentView = nil;
      qw.fIsDeleted = YES;

      if (qw.fMainWindow) {
         [qw.fMainWindow removeChildWindow : qw];
         qw.fMainWindow = nil;
      }

      [qw orderOut:nil];
   }

   fDrawables.erase(drawableIter);//StrongReference should do work here.
}

//______________________________________________________________________________
Handle_t CocoaPrivate::RegisterGLContext(NSOpenGLContext *glContext)
{
   assert(fGLContextToHandle.find(glContext) == fGLContextToHandle.end() && "RegisterGLContext, context was registered already");

   //Strong es-guarantee guarantee - if we have an exception, everything is rolled-back.

   bool contextInserted = false;
   try {
      fHandleToGLContext[fFreeGLContextID] = glContext;
      contextInserted = true;
      fGLContextToHandle[glContext] = fFreeGLContextID;
   } catch (const std::exception &) {//bad alloc in one of two insertions.
      if (contextInserted)
         fHandleToGLContext.erase(fHandleToGLContext.find(fFreeGLContextID));
      throw;
   }

   return fFreeGLContextID++;
}

//______________________________________________________________________________
void CocoaPrivate::DeleteGLContext(Handle_t contextID)
{
   assert(fHandleToGLContext.find(contextID) != fHandleToGLContext.end() && "DeleteGLContext, bad context id");

   handle2ctx_map::iterator h2cIt = fHandleToGLContext.find(contextID);

   ctx2handle_map::iterator c2hIt = fGLContextToHandle.find(h2cIt->second.Get());
   assert(c2hIt != fGLContextToHandle.end() && "DeleteGLContext, inconsistent context map");

   fGLContextToHandle.erase(c2hIt);
   fHandleToGLContext.erase(h2cIt);//RAII does work here.
}

//______________________________________________________________________________
NSOpenGLContext *CocoaPrivate::GetGLContextForHandle(Handle_t ctxID)
{
   if (fHandleToGLContext.find(ctxID) == fHandleToGLContext.end())
      return nil;

   return fHandleToGLContext[ctxID].Get();
}

//______________________________________________________________________________
Handle_t CocoaPrivate::GetHandleForGLContext(NSOpenGLContext *glContext)
{
   if (fGLContextToHandle.find(glContext) == fGLContextToHandle.end())
      return Handle_t();

   return fGLContextToHandle[glContext];
}

//______________________________________________________________________________
void CocoaPrivate::SetFakeGLWindow(QuartzWindow *fakeWin)
{
   fFakeGLWindow.Reset(fakeWin);
}

//______________________________________________________________________________
QuartzWindow *CocoaPrivate::GetFakeGLWindow()
{
   return fFakeGLWindow.Get();
}

//______________________________________________________________________________
void CocoaPrivate::ReplaceDrawable(Drawable_t drawableID, NSObject *nsObj)
{
   drawable_iterator drawableIter = fDrawables.find(drawableID);
   assert(drawableIter != fDrawables.end() && "ReplaceDrawable, can not replace non existing drawable");
   drawableIter->second.Reset(nsObj);
}

}//Details
}//MacOSX
}//ROOT
