// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveScene_hxx
#define ROOT_TEveScene_hxx

#include "ROOT/TEveElement.hxx"

class TExMap;

namespace ROOT { namespace Experimental
{

// class TEvePad;

/******************************************************************************/
// TEveScene
/******************************************************************************/

class TEveScene : public TEveElementList
{
private:
   TEveScene(const TEveScene&);            // Not implemented
   TEveScene& operator=(const TEveScene&); // Not implemented

protected:
   //TEvePad     *fPad;

   Bool_t       fChanged;
   Bool_t       fSmartRefresh;
   Bool_t       fHierarchical;

   // void RetransHierarchicallyRecurse(TEveElement* el, const TEveTrans& tp);

public:
   TEveScene(const char* n="TEveScene", const char* t="");
   virtual ~TEveScene();

   virtual void CollectSceneParents(List_t& scenes);

   virtual Bool_t SingleRnrState() const { return kTRUE; }

   void   Changed()         { fChanged = kTRUE; }
   Bool_t IsChanged() const { return fChanged;  }

   void   SetHierarchical(Bool_t h) { fHierarchical = h;    }
   Bool_t GetHierarchical()   const { return fHierarchical; }

   // void   Repaint(Bool_t dropLogicals=kFALSE);
   // void   RetransHierarchically();

   virtual void SetName(const char* n);

   // virtual void Paint(Option_t* option = "");

   // void DestroyElementRenderers(TEveElement* element);
   // void DestroyElementRenderers(TObject* rnrObj);

   ClassDef(ROOT::Experimental::TEveScene, 0); // Reve representation of TGLScene.
};


/******************************************************************************/
// TEveSceneList
/******************************************************************************/

class TEveSceneList : public TEveElementList
{
private:
   TEveSceneList(const TEveSceneList&);            // Not implemented
   TEveSceneList& operator=(const TEveSceneList&); // Not implemented

protected:

public:
   TEveSceneList(const char* n="TEveSceneList", const char* t="");
   virtual ~TEveSceneList() {}

   void DestroyScenes();

   // void RepaintChangedScenes(Bool_t dropLogicals);
   // void RepaintAllScenes(Bool_t dropLogicals);

   // void DestroyElementRenderers(TEveElement* element);

   // void ProcessSceneChanges(Bool_t dropLogicals, TExMap* stampMap);

   ClassDef(ROOT::Experimental::TEveSceneList, 0); // List of Scenes providing common operations on TEveScene collections.
};

}}

#endif
