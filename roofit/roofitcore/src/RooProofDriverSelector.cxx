#define RooProofDriverSelector_cxx
// The class definition in RooProofDriverSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("RooProofDriverSelector.C")
// Root > T->Process("RooProofDriverSelector.C","some options")
// Root > T->Process("RooProofDriverSelector.C+")
//

#include "RooProofDriverSelector.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TRandom.h"
#include "RooRandom.h"
#include "RooAbsStudy.h"
#include "RooStudyPackage.h"
#include "RooGlobalFunc.h"

using namespace RooFit;
using namespace std;

void RooProofDriverSelector::SlaveBegin(TTree * /*tree*/) 
{  
  // Retrieve study pack 
  _pkg=0 ;
  if (fInput) { 
    TIterator* iter = fInput->MakeIterator() ;
    TObject* obj ;
    while((obj=iter->Next())) {
      RooStudyPackage* tmp = dynamic_cast<RooStudyPackage*>(obj) ;
      if (tmp) {
	_pkg = tmp ;
      }
     }
  }
  if (_pkg==0) {
    cout << "RooProofDriverSelector::SlaveBegin() no RooStudyPackage found, aborting process" << endl ;
    fStatus = kAbortProcess ;
  } else {
    cout << "workspace contents = " << endl ;
    _pkg->wspace().Print() ;
    
    // Initialize study pack
    seed = _pkg->initRandom() ;
    _pkg->initialize() ;
  }
    
}

Bool_t RooProofDriverSelector::Process(Long64_t entry)
{
  cout << "RooProofDriverSelector::Process(" << entry << ")" << endl ;
  _pkg->runOne() ;
  return kTRUE;
}


void RooProofDriverSelector::SlaveTerminate() 
{
  _pkg->finalize() ;
  _pkg->exportData(fOutput,seed) ;
}



void RooProofDriverSelector::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("i", &i, &b_i);
}

Bool_t RooProofDriverSelector::Notify()
{
   return kTRUE;
}

