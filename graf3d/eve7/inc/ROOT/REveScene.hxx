// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveScene
#define ROOT7_REveScene

#include <ROOT/REveElement.hxx>

#include "TClass.h"

#include <vector>
#include <memory>

namespace ROOT {
namespace Experimental {

class REveClient;
class REveManager;

/******************************************************************************/
// REveScene
// REve representation of TGLScene.
/******************************************************************************/

class REveScene : public REveElement
{
   friend class REveManager;

private:
   REveScene(const REveScene &) = delete;
   REveScene &operator=(const REveScene &) = delete;

protected:
   struct SceneCommand
   {
      std::string fName;
      std::string fIcon;
      std::string fElementClass;
      std::string fAction;
      ElementId_t fElementId{0};

      SceneCommand(const std::string& name, const std::string& icon,
                   const REveElement* element, const std::string& action) :
         fName(name),
         fIcon(icon),
         fElementClass(element->IsA()->GetName()),
         fAction(action),
         fElementId(element->GetElementId())
      {}
   };

   Bool_t fSmartRefresh{kTRUE};            ///<!
   Bool_t fHierarchical{kFALSE};           ///<!

   Bool_t fAcceptingChanges{kFALSE};       ///<!
   Bool_t fChanged{kFALSE};                ///<!
   // Changed or/and added.
   // XXXX can change to vector (element checks if already registered now).
   List_t  fChangedElements;               ///<!
   // For the following two have to re-think how the hierarchy will be handled.
   // If I remove a parent, i have to remove all the children.
   // So this has to be done right on both sides (on eve element and here).
   // I might need a set, so i can easily check if parent is in the removed / added list already.
   std::vector<ElementId_t> fRemovedElements; ///<!

   std::vector<std::unique_ptr<REveClient>> fSubscribers; ///<!

   List_t fElsWithBinaryData;
   std::string fOutputJson;               ///<!
   std::vector<char> fOutputBinary;       ///<!
   Int_t fTotalBinarySize;                ///<!

   std::vector<SceneCommand> fCommands;   ///<!

   // void RetransHierarchicallyRecurse(REveElement* el, const REveTrans& tp);

public:
   REveScene(const std::string &n = "REveScene", const std::string &t = "");
   virtual ~REveScene();

   Bool_t SingleRnrState() const override { return kTRUE; }

   void   SetHierarchical(Bool_t h) { fHierarchical = h; }
   Bool_t GetHierarchical() const { return fHierarchical; }

   void   Changed() { fChanged = kTRUE; } // AMT ??? depricated
   Bool_t IsChanged() const;

   Bool_t IsAcceptingChanges() const { return fAcceptingChanges; }
   void BeginAcceptingChanges();
   void SceneElementChanged(REveElement *element);
   void SceneElementRemoved(ElementId_t id);
   void EndAcceptingChanges();
   void ProcessChanges();

   void StreamElements();
   void StreamJsonRecurse(REveElement *el, nlohmann::json &jobj);

   // void   Repaint(Bool_t dropLogicals=kFALSE);
   // void   RetransHierarchically();

   // virtual void Paint(Option_t* option = "");

   // void DestroyElementRenderers(REveElement* element);
   // void DestroyElementRenderers(TObject* rnrObj);
   void StreamRepresentationChanges();
   void SendChangesToSubscribers();

   Bool_t HasSubscribers() const { return !fSubscribers.empty(); }
   void AddSubscriber(std::unique_ptr<REveClient> &&sub);
   void RemoveSubscriber(unsigned int);

   void AddCommand(const std::string &name, const std::string &icon, const REveElement *element, const std::string &action);
};

/******************************************************************************/
// REveSceneList
// List of Scenes providing common operations on REveScene collections.
/******************************************************************************/

class REveSceneList : public REveElement
{
private:
   REveSceneList(const REveSceneList &) = delete;
   REveSceneList &operator=(const REveSceneList &) = delete;

protected:
public:
   REveSceneList(const std::string &n = "REveSceneList", const std::string &t = "");
   virtual ~REveSceneList() {}

   void DestroyScenes();

   // void RepaintChangedScenes(Bool_t dropLogicals);
   // void RepaintAllScenes(Bool_t dropLogicals);

   // void DestroyElementRenderers(REveElement* element);
   void AcceptChanges(bool);

   void ProcessSceneChanges();
};

} // namespace Experimental
} // namespace ROOT

#endif
