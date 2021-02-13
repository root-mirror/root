// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveSceneInfo.hxx>
#include <ROOT/REveScene.hxx>

using namespace ROOT::Experimental;

/** \class REveSceneInfo
\ingroup REve
Representation of a REveScene in a REveViewer. This allows for
viewer specific settings to be applied to a scene, e.g., global position
in viewer coordinate system.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveSceneInfo::REveSceneInfo(REveViewer* viewer, REveScene* scene) :
   REveElement (Form("SI - %s", scene->GetCName()),
                Form("REveSceneInfo of scene '%s'", scene->GetCName())),
   fViewer     (viewer),
   fScene      (scene)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveSceneInfo::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);

   j["fSceneId"] = fScene->GetElementId();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveElement.
/// REveSceneInfo does not accept children.

Bool_t REveSceneInfo::AcceptElement(REveElement* /*el*/)
{
   static const REveException eH("REveSceneInfo::AcceptElement ");

   // gEve->SetStatusLine(eH + "this class does not accept children.");
   return kFALSE;
}
