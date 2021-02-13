// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveVSD.hxx>

#include "TTree.h"

using namespace ROOT::Experimental;

/** \class REveVSD
\ingroup REve
Visualization Summary Data - a collection of trees holding standard
event data in experiment independent format.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveVSD::REveVSD(const char* , const char*) :
   TObject(),

   fFile      (nullptr),
   fDirectory (nullptr),
   fBuffSize  (32000),
   fVerbose   (0),

   fTreeK  (nullptr),
   fTreeH  (nullptr),
   fTreeC  (nullptr),
   fTreeR  (nullptr),
   fTreeKK (nullptr),
   fTreeV0 (nullptr),
   fTreeCC (nullptr),
   fTreeGI (nullptr),

   fK(),  fpK (&fK),
   fH(),  fpH (&fH),
   fC(),  fpC (&fC),
   fR(),  fpR (&fR),
   fKK(), fpKK(&fKK),
   fV0(), fpV0(&fV0),
   fCC(), fpCC(&fCC),
   fGI(), fpGI(&fGI)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveVSD::~REveVSD()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set directory in which the trees are (or will be) created.

void REveVSD::SetDirectory(TDirectory* dir)
{
   fDirectory = dir;
}

////////////////////////////////////////////////////////////////////////////////
/// Create internal trees.

void REveVSD::CreateTrees()
{
   fDirectory->cd();
   fTreeK  = new TTree("Kinematics", "Simulated tracks.");
   fTreeH  = new TTree("Hits",       "Combined detector hits.");
   fTreeC  = new TTree("Clusters",   "Reconstructed clusters.");
   fTreeR  = new TTree("RecTracks",  "Reconstructed tracks.");
   fTreeKK = new TTree("RecKinks",   "Reconstructed kinks.");
   fTreeV0 = new TTree("RecV0s",     "Reconstructed V0s.");
   fTreeCC = new TTree("RecCascades","Reconstructed cascades.");
   fTreeGI = new TTree("REveMCRecCrossRef", "Objects prepared for cross query.");
}

////////////////////////////////////////////////////////////////////////////////
/// Delete internal trees.

void REveVSD::DeleteTrees()
{
   delete fTreeK;  fTreeK  = nullptr;
   delete fTreeH;  fTreeH  = nullptr;
   delete fTreeC;  fTreeC  = nullptr;
   delete fTreeR;  fTreeR  = nullptr;
   delete fTreeV0; fTreeV0 = nullptr;
   delete fTreeKK; fTreeKK = nullptr;
   delete fTreeGI; fTreeGI = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Create internal VSD branches.

void REveVSD::CreateBranches()
{
   if (fTreeK)  fTreeK ->Branch("K",  "REveMCTrack",  &fpK);
   if (fTreeH)  fTreeH ->Branch("H",  "REveHit",      &fpH);
   if (fTreeC)  fTreeC ->Branch("C",  "REveCluster",  &fpC);
   if (fTreeR)  fTreeR ->Branch("R",  "REveRecTrack", &fpR);
   if (fTreeKK) fTreeKK->Branch("KK", "REveRecKink",  &fpKK);
   if (fTreeV0) fTreeV0->Branch("V0", "REveRecV0",    &fpV0);

   if (fTreeGI)
   {
      fTreeGI->Branch("GI", "REveMCRecCrossRef",  &fpGI);
      fTreeGI->Branch("K.", "REveMCTrack",  &fpK);
      fTreeGI->Branch("R.", "REveRecTrack", &fpR);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set branche addresses of internal trees.

void REveVSD::SetBranchAddresses()
{
   if (fTreeK)  fTreeK ->SetBranchAddress("K",  &fpK);
   if (fTreeH)  fTreeH ->SetBranchAddress("H",  &fpH);
   if (fTreeC)  fTreeC ->SetBranchAddress("C",  &fpC);
   if (fTreeR)  fTreeR ->SetBranchAddress("R",  &fpR);
   if (fTreeKK) fTreeKK->SetBranchAddress("KK", &fpKK);
   if (fTreeV0) fTreeV0->SetBranchAddress("V0", &fpV0);

   if (fTreeGI)
   {
      fTreeGI->SetBranchAddress("GI", &fpGI);
      fTreeGI->SetBranchAddress("K.", &fpK);
      fTreeGI->SetBranchAddress("R.", &fpR);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Does nothing here ... reimplemented in sub-classes.

void REveVSD::WriteTrees()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Load internal trees from directory.

void REveVSD::LoadTrees()
{
   static const REveException eH("REveVSD::LoadTrees");

   if (!fDirectory)
      throw eH + " directory not set.";

   fTreeK = (TTree*) fDirectory->Get("Kinematics");
   if (!fTreeK && fVerbose)
      Error("REveVSD::LoadTrees","Kinematics not available in fDirectory %s.", fDirectory->GetName());

   fTreeH = (TTree*) fDirectory->Get("Hits");
   if (!fTreeH && fVerbose)
      Error("REveVSD::LoadTrees", "Hits not available in fDirectory %s.", fDirectory->GetName());

   fTreeC = (TTree*) fDirectory->Get("Clusters");
   if (!fTreeC && fVerbose)
      Error("REveVSD::LoadTrees", "Clusters not available in fDirectory %s.", fDirectory->GetName());

   fTreeR = (TTree*) fDirectory->Get("RecTracks");
   if (!fTreeR && fVerbose)
      Error("REveVSD::LoadTrees", "RecTracks not available in fDirectory %s.", fDirectory->GetName());

   fTreeKK =  (TTree*) fDirectory->Get("RecKinks");
   if (!fTreeKK && fVerbose)
      Error("REveVSD::LoadTrees", "RecKinks not available in fDirectory %s.", fDirectory->GetName());

   fTreeV0 =  (TTree*) fDirectory->Get("RecV0s");
   if (!fTreeV0  && fVerbose)
      Error("REveVSD::LoadTrees", "RecV0 not available in fDirectory %s.", fDirectory->GetName());

   fTreeGI = (TTree*)fDirectory->Get("REveMCRecCrossRef");
   if(!fTreeGI && fVerbose)
      Error("REveVSD::LoadTrees", "REveMCRecCrossRef not available in fDirectory %s.", fDirectory->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Disable TObject streamers for those VSD structs that inherit from
/// TObject directly.

void REveVSD::DisableTObjectStreamersForVSDStruct()
{
   // REveVector is not TObject

   // REveMCTrack derives from TParticle
   TParticle::Class()->IgnoreTObjectStreamer(true);

   // REveHit::Class()->IgnoreTObjectStreamer(true);
   // REveCluster::Class()->IgnoreTObjectStreamer(true);

   // REveRecTrack::Class()->IgnoreTObjectStreamer(true);
   // REveRecKink derives from REveRecTrack

   // REveRecV0::Class()->IgnoreTObjectStreamer(true);

   // REveMCRecCrossRef::Class()->IgnoreTObjectStreamer(true);
}
