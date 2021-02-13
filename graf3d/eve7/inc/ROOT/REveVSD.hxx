// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveVSD
#define ROOT7_REveVSD

#include <ROOT/REveUtil.hxx>
#include <ROOT/REveVSDStructs.hxx>

class TTree;
class TFile;
class TDirectory;

namespace ROOT {
namespace Experimental {

class REveVSD : public TObject {
   REveVSD(const REveVSD &) = delete;
   REveVSD &operator=(const REveVSD &) = delete;

protected:
   TFile *fFile{nullptr};           //!
   TDirectory *fDirectory{nullptr}; //!

   Int_t fBuffSize{0}; //!
   Int_t fVerbose{0};  //!

public:
   TTree *fTreeK{nullptr};  //! Kinematics.
   TTree *fTreeH{nullptr};  //! Hits.
   TTree *fTreeC{nullptr};  //! Clusters.
   TTree *fTreeR{nullptr};  //! Reconstructed tracks.
   TTree *fTreeKK{nullptr}; //! Kinks.
   TTree *fTreeV0{nullptr}; //! VO's.
   TTree *fTreeCC{nullptr}; //! Cascades.
   TTree *fTreeGI{nullptr}; //! Sim-Rec cross references.

   REveMCTrack fK, *fpK{nullptr};         //!
   REveHit fH, *fpH{nullptr};             //!
   REveCluster fC, *fpC{nullptr};         //!
   REveRecTrack fR, *fpR{nullptr};        //!
   REveRecKink fKK, *fpKK{nullptr};       //!
   REveRecV0 fV0, *fpV0{nullptr};         //!
   REveRecCascade fCC, *fpCC{nullptr};    //!
   REveMCRecCrossRef fGI, *fpGI{nullptr}; //!

public:
   REveVSD(const char *name = "REveVSD", const char *title = "");
   virtual ~REveVSD();

   virtual void SetDirectory(TDirectory *dir);

   virtual void CreateTrees();
   virtual void DeleteTrees();

   virtual void CreateBranches();
   virtual void WriteTrees();

   virtual void LoadTrees();
   virtual void SetBranchAddresses();

   static void DisableTObjectStreamersForVSDStruct();

   ClassDef(REveVSD, 1); // Visualization Summary Data - a collection of trees holding standard event data in experiment independent format.
};

} // namespace Experimental
} // namespace ROOT

#endif
