//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 10 11:01:34 2009 by ROOT version 5.23/05
// from TTree t/t
// found on file: Memory Directory
//////////////////////////////////////////////////////////

#ifndef RooProofDriverSelector_h
#define RooProofDriverSelector_h

#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
class RooStudyPackage ;
class TIterator ;

class RooProofDriverSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           i;

   // List of branches
   TBranch        *b_i;   //!

   RooProofDriverSelector(TTree * /*tree*/ =0) { b_i = 0 ; _pkg = 0 ; fChain = 0 ; }
   virtual ~RooProofDriverSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree* tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual void    SlaveTerminate() ;
   virtual TList  *GetOutputList() const { return fOutput; }

   RooStudyPackage* _pkg ;
   Int_t      seed ;

   ClassDef(RooProofDriverSelector,0);
};

#endif

