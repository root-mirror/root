/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RProvider.hxx>

#include "TString.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TBranchBrowsable.h"
#include "TTree.h"
#include "TH1.h"
#include "TDirectory.h"

using namespace ROOT::Experimental::Browsable;

/** Provider for drawing of branches / leafs in the TTree */

template<typename T>
class TLeafProvider : public RProvider {
public:

   virtual ~TLeafProvider() = default;

   TH1 *DrawTree(TTree *ttree, const std::string &expr, const std::string &hname)
   {
      if (!ttree)
         return nullptr;

      std::string expr2 = expr + ">>htemp_tree_draw";

      ttree->Draw(expr2.c_str(),"","goff");

      if (!gDirectory)
         return nullptr;

      auto htemp = dynamic_cast<TH1*>(gDirectory->FindObject("htemp_tree_draw"));

      if (!htemp)
         return nullptr;

      htemp->SetDirectory(nullptr);
      htemp->SetName(hname.c_str());

      return htemp;
   }

   TH1 *DrawLeaf(std::unique_ptr<RHolder> &obj)
   {
      auto tleaf = obj->get_object<TLeaf>();
      if (!tleaf)
         return nullptr;

      return DrawTree(tleaf->GetBranch()->GetTree(), tleaf->GetName(), tleaf->GetName());
   }

   TH1 *DrawBranchElement(std::unique_ptr<RHolder> &obj)
   {
      auto tbranch = obj->get_object<TBranchElement>();
      if (!tbranch)
         return nullptr;

      // there are sub-branches, plain TTree::Draw does not work
      if (tbranch->GetListOfBranches()->GetEntriesFast() > 0)
         return nullptr;

      // just copy and paste code from TBranchElement::Browse
      TString slash("/");
      TString escapedSlash("\\/");
      TString name = tbranch->GetName();
      Int_t pos = name.First('[');
      if (pos != kNPOS)
         name.Remove(pos);
      if (tbranch->GetMother()) {
         TString mothername = tbranch->GetMother()->GetName();
         pos = mothername.First('[');
         if (pos != kNPOS) {
            mothername.Remove(pos);
         }
         Int_t len = mothername.Length();
         if (len) {
            if (mothername(len-1) != '.') {
               // We do not know for sure whether the mother's name is
               // already preprended.  So we need to check:
               //    a) it is prepended
               //    b) it is NOT the name of a daugher (i.e. mothername.mothername exist)
               TString doublename = mothername;
               doublename.Append(".");
               Int_t isthere = (name.Index(doublename) == 0);
               if (!isthere) {
                  name.Prepend(doublename);
               } else {
                  if (tbranch->GetMother()->FindBranch(mothername)) {
                     doublename.Append(mothername);
                     isthere = (name.Index(doublename) == 0);
                     if (!isthere) {
                        mothername.Append(".");
                        name.Prepend(mothername);
                     }
                  } else {
                     // Nothing to do because the mother's name is
                     // already in the name.
                  }
               }
            } else {
               // If the mother's name end with a dot then
               // the daughter probably already contains the mother's name
               if (name.Index(mothername) == kNPOS) {
                  name.Prepend(mothername);
               }
            }
         }
      }
      name.ReplaceAll(slash, escapedSlash);

      return DrawTree(tbranch->GetTree(), name.Data(), tbranch->GetName());
   }

   TH1 *DrawBranchBrowsable(std::unique_ptr<RHolder> &obj)
   {
      auto browsable = obj->get_object<TVirtualBranchBrowsable>();
      if (!browsable)
         return nullptr;

      auto cl = browsable->GetClassType();

      bool can_draw  = (!cl || (cl->GetCollectionProxy() && cl->GetCollectionProxy()->GetType() > 0));
      if (!can_draw)
         return nullptr;

      auto br = browsable->GetBranch();
      if (!br)
         return nullptr;

      TString name;
      browsable->GetScope(name);

      // If this is meant to be run on the collection
      // we need to "move" the "@" from branch.@member
      // to branch@.member
      name.ReplaceAll(".@","@.");
      name.ReplaceAll("->@","@->");

      return DrawTree(br->GetTree(), name.Data(), browsable->GetName());
   }

};
