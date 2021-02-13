// @(#)root/treeplayer:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreePlayer

Implement some of the functionality of the class TTree requiring access to
extra libraries (Histogram, display, etc).
*/

#include "TTreePlayer.h"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TFile.h"
#include "TEventList.h"
#include "TEntryList.h"
#include "TBranchObject.h"
#include "TBranchElement.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TLeafObject.h"
#include "TLeafF.h"
#include "TLeafD.h"
#include "TLeafC.h"
#include "TLeafB.h"
#include "TLeafI.h"
#include "TLeafS.h"
#include "TMath.h"
#include "TH1.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TText.h"
#include "TDirectory.h"
#include "TClonesArray.h"
#include "TClass.h"
#include "TVirtualPad.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "TStyle.h"
#include "Foption.h"
#include "TTreeResult.h"
#include "TTreeRow.h"
#include "TPrincipal.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "THLimitsFinder.h"
#include "TSelectorDraw.h"
#include "TSelectorEntries.h"
#include "TPluginManager.h"
#include "TObjString.h"
#include "TTreeProxyGenerator.h"
#include "TTreeReaderGenerator.h"
#include "TTreeIndex.h"
#include "TChainIndex.h"
#include "TRefProxy.h"
#include "TRefArrayProxy.h"
#include "TVirtualMonitoring.h"
#include "TTreeCache.h"
#include "TVirtualMutex.h"
#include "ThreadLocalStorage.h"
#include "strlcpy.h"
#include "snprintf.h"

#include "HFitInterface.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Math/MinimizerOptions.h"


R__EXTERN Foption_t Foption;

TVirtualFitter *tFitter = nullptr;

ClassImp(TTreePlayer);

////////////////////////////////////////////////////////////////////////////////
/// Default Tree constructor.

TTreePlayer::TTreePlayer()
{
   fTree           = 0;
   fScanFileName   = 0;
   fScanRedirect   = kFALSE;
   fSelectedRows   = 0;
   fDimension      = 0;
   fHistogram      = 0;
   fFormulaList    = new TList();
   fFormulaList->SetOwner(kTRUE);
   fSelector         = new TSelectorDraw();
   fSelectorFromFile = 0;
   fSelectorClass    = 0;
   fSelectorUpdate   = 0;
   fInput            = new TList();
   fInput->Add(new TNamed("varexp",""));
   fInput->Add(new TNamed("selection",""));
   fSelector->SetInputList(fInput);
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfCleanups()->Add(this);
   }
   TClass::GetClass("TRef")->AdoptReferenceProxy(new TRefProxy());
   TClass::GetClass("TRefArray")->AdoptReferenceProxy(new TRefArrayProxy());
}

////////////////////////////////////////////////////////////////////////////////
/// Tree destructor.

TTreePlayer::~TTreePlayer()
{
   delete fFormulaList;
   delete fSelector;
   DeleteSelectorFromFile();
   fInput->Delete();
   delete fInput;
   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfCleanups()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Build the index for the tree (see TTree::BuildIndex)

TVirtualIndex *TTreePlayer::BuildIndex(const TTree *T, const char *majorname, const char *minorname)
{
   TVirtualIndex *index;
   if (dynamic_cast<const TChain*>(T)) {
      index = new TChainIndex(T, majorname, minorname);
      if (index->IsZombie()) {
         delete index;
         Error("BuildIndex", "Creating a TChainIndex unsuccessful - switching to TTreeIndex");
      }
      else
         return index;
   }
   return new TTreeIndex(T,majorname,minorname);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a Tree with selection, make a clone of this Tree header, then copy the
/// selected entries.
///
/// -  selection is a standard selection expression (see TTreePlayer::Draw)
/// -  option is reserved for possible future use
/// -  nentries is the number of entries to process (default is all)
/// -  first is the first entry to process (default is 0)
///
/// IMPORTANT: The copied tree stays connected with this tree until this tree
/// is deleted.  In particular, any changes in branch addresses
/// in this tree are forwarded to the clone trees.  Any changes
/// made to the branch addresses of the copied trees are over-ridden
/// anytime this tree changes its branch addresses.
/// Once this tree is deleted, all the addresses of the copied tree
/// are reset to their default values.
///
/// The following example illustrates how to copy some events from the Tree
/// generated in $ROOTSYS/test/Event
/// ~~~{.cpp}
///   gSystem->Load("libEvent");
///   TFile f("Event.root");
///   TTree *T = (TTree*)f.Get("T");
///   Event *event = new Event();
///   T->SetBranchAddress("event",&event);
///   TFile f2("Event2.root","recreate");
///   TTree *T2 = T->CopyTree("fNtrack<595");
///   T2->Write();
/// ~~~

TTree *TTreePlayer::CopyTree(const char *selection, Option_t *, Long64_t nentries,
                             Long64_t firstentry)
{

   // we make a copy of the tree header
   TTree *tree = fTree->CloneTree(0);
   if (tree == 0) return 0;

   // The clone should not delete any shared i/o buffers.
   TObjArray* branches = tree->GetListOfBranches();
   Int_t nb = branches->GetEntriesFast();
   for (Int_t i = 0; i < nb; ++i) {
      TBranch* br = (TBranch*) branches->UncheckedAt(i);
      if (br->InheritsFrom(TBranchElement::Class())) {
         ((TBranchElement*) br)->ResetDeleteObject();
      }
   }

   Long64_t entry,entryNumber;
   nentries = GetEntriesToProcess(firstentry, nentries);

   // Compile selection expression if there is one
   TTreeFormula *select = 0; // no need to interfere with fSelect since we
                             // handle the loop explicitly below and can call
                             // UpdateFormulaLeaves ourselves.
   if (strlen(selection)) {
      select = new TTreeFormula("Selection",selection,fTree);
      if (!select || !select->GetNdim()) {
         delete select;
         delete tree;
         return 0;
      }
      fFormulaList->Add(select);
   }

   //loop on the specified entries
   Int_t tnumber = -1;
   for (entry=firstentry;entry<firstentry+nentries;entry++) {
      entryNumber = fTree->GetEntryNumber(entry);
      if (entryNumber < 0) break;
      Long64_t localEntry = fTree->LoadTree(entryNumber);
      if (localEntry < 0) break;
      if (tnumber != fTree->GetTreeNumber()) {
         tnumber = fTree->GetTreeNumber();
         if (select) select->UpdateFormulaLeaves();
      }
      if (select) {
         Int_t ndata = select->GetNdata();
         Bool_t keep = kFALSE;
         for(Int_t current = 0; current<ndata && !keep; current++) {
            keep |= (select->EvalInstance(current) != 0);
         }
         if (!keep) continue;
      }
      fTree->GetEntry(entryNumber);
      tree->Fill();
   }
   fFormulaList->Clear();
   return tree;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete any selector created by this object.
/// The selector has been created using TSelector::GetSelector(file)

void TTreePlayer::DeleteSelectorFromFile()
{
   if (fSelectorFromFile && fSelectorClass) {
      if (fSelectorClass->IsLoaded()) {
         delete fSelectorFromFile;
      }
   }
   fSelectorFromFile = 0;
   fSelectorClass = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the result of a C++ script.
///
/// The macrofilename and optionally cutfilename are assumed to contain
/// at least a method with the same name as the file.  The method
/// should return a value that can be automatically cast to
/// respectively a double and a boolean.
///
/// Both methods will be executed in a context such that the
/// branch names can be used as C++ variables. This is
/// accomplished by generating a TTreeProxy (see MakeProxy)
/// and including the files in the proper location.
///
/// If the branch name can not be used a proper C++ symbol name,
/// it will be modified as follow:
///    - white spaces are removed
///    - if the leading character is not a letter, an underscore is inserted
///    - < and > are replace by underscores
///    - * is replaced by st
///    - & is replaced by rf
///
/// If a cutfilename is specified, for each entry, we execute
/// ~~~{.cpp}
///   if (cutfilename()) htemp->Fill(macrofilename());
/// ~~~
/// If no cutfilename is specified, for each entry we execute
/// ~~~{.cpp}
///   htemp(macrofilename());
/// ~~~
/// The default for the histogram are the same as for
/// TTreePlayer::DrawSelect

Long64_t TTreePlayer::DrawScript(const char* wrapperPrefix,
                                 const char *macrofilename, const char *cutfilename,
                                 Option_t *option, Long64_t nentries, Long64_t firstentry)
{
   if (!macrofilename || strlen(macrofilename)==0) return 0;

   TString aclicMode;
   TString arguments;
   TString io;
   TString realcutname;
   if (cutfilename && strlen(cutfilename))
      realcutname =  gSystem->SplitAclicMode(cutfilename, aclicMode, arguments, io);

   // we ignore the aclicMode for the cutfilename!
   TString realname = gSystem->SplitAclicMode(macrofilename, aclicMode, arguments, io);

   TString selname = wrapperPrefix;

   ROOT::Internal::TTreeProxyGenerator gp(fTree,realname,realcutname,selname,option,3);

   selname = gp.GetFileName();
   if (aclicMode.Length()==0) {
      Warning("DrawScript","TTreeProxy does not work in interpreted mode yet. The script will be compiled.");
      aclicMode = "+";
   }
   selname.Append(aclicMode);

   Info("DrawScript","%s",Form("Will process tree/chain using %s",selname.Data()));
   Long64_t result = fTree->Process(selname,option,nentries,firstentry);
   fTree->SetNotify(0);

   // could delete the file selname+".h"
   // However this would remove the optimization of avoiding a useless
   // recompilation if the user ask for the same thing twice!

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expression varexp for specified entries that matches the selection.
/// Returns -1 in case of error or number of selected events in case of success.
///
/// See the documentation of TTree::Draw for the complete details.

Long64_t TTreePlayer::DrawSelect(const char *varexp0, const char *selection, Option_t *option,Long64_t nentries, Long64_t firstentry)
{
   if (fTree->GetEntriesFriend() == 0) return 0;

   // Let's see if we have a filename as arguments instead of
   // a TTreeFormula expression.

   TString possibleFilename = varexp0;
   Ssiz_t dot_pos = possibleFilename.Last('.');
   if ( dot_pos != kNPOS
       && possibleFilename.Index("Alt$")<0 && possibleFilename.Index("Entries$")<0
       && possibleFilename.Index("LocalEntries$")<0
       && possibleFilename.Index("Length$")<0  && possibleFilename.Index("Entry$")<0
       && possibleFilename.Index("LocalEntry$")<0
       && possibleFilename.Index("Min$")<0 && possibleFilename.Index("Max$")<0
       && possibleFilename.Index("MinIf$")<0 && possibleFilename.Index("MaxIf$")<0
       && possibleFilename.Index("Iteration$")<0 && possibleFilename.Index("Sum$")<0
       && possibleFilename.Index(">")<0 && possibleFilename.Index("<")<0
       && gSystem->IsFileInIncludePath(possibleFilename.Data())) {

      if (selection && strlen(selection) && !gSystem->IsFileInIncludePath(selection)) {
         Error("DrawSelect",
               "Drawing using a C++ file currently requires that both the expression and the selection are files\n\t\"%s\" is not a file",
               selection);
         return 0;
      }
      return DrawScript("generatedSel",varexp0,selection,option,nentries,firstentry);

   } else {
      possibleFilename = selection;
      if (possibleFilename.Index("Alt$")<0 && possibleFilename.Index("Entries$")<0
          && possibleFilename.Index("LocalEntries$")<0
          && possibleFilename.Index("Length$")<0  && possibleFilename.Index("Entry$")<0
          && possibleFilename.Index("LocalEntry$")<0
          && possibleFilename.Index("Min$")<0 && possibleFilename.Index("Max$")<0
          && possibleFilename.Index("MinIf$")<0 && possibleFilename.Index("MaxIf$")<0
          && possibleFilename.Index("Iteration$")<0 && possibleFilename.Index("Sum$")<0
          && possibleFilename.Index(">")<0 && possibleFilename.Index("<")<0
          && gSystem->IsFileInIncludePath(possibleFilename.Data())) {

         Error("DrawSelect",
               "Drawing using a C++ file currently requires that both the expression and the selection are files\n\t\"%s\" is not a file",
               varexp0);
         return 0;
      }
   }

   Long64_t oldEstimate  = fTree->GetEstimate();
   TEventList *evlist  = fTree->GetEventList();
   TEntryList *elist = fTree->GetEntryList();
   if (evlist && elist){
      elist->SetBit(kCanDelete, kTRUE);
   }
   TNamed *cvarexp    = (TNamed*)fInput->FindObject("varexp");
   TNamed *cselection = (TNamed*)fInput->FindObject("selection");
   if (cvarexp) cvarexp->SetTitle(varexp0);
   if (cselection) cselection->SetTitle(selection);

   TString opt = option;
   opt.ToLower();
   Bool_t optpara   = kFALSE;
   Bool_t optcandle = kFALSE;
   Bool_t optgl5d   = kFALSE;
   Bool_t optnorm   = kFALSE;
   if (opt.Contains("norm")) {optnorm = kTRUE; opt.ReplaceAll("norm",""); opt.ReplaceAll(" ","");}
   if (opt.Contains("para")) optpara = kTRUE;
   if (opt.Contains("candle")) optcandle = kTRUE;
   if (opt.Contains("gl5d")) optgl5d = kTRUE;
   Bool_t pgl = gStyle->GetCanvasPreferGL();
   if (optgl5d) {
      fTree->SetEstimate(fTree->GetEntries());
      if (!gPad) {
         if (pgl == kFALSE) gStyle->SetCanvasPreferGL(kTRUE);
         gROOT->ProcessLineFast("new TCanvas();");
      }
   }

   // Do not process more than fMaxEntryLoop entries
   if (nentries > fTree->GetMaxEntryLoop()) nentries = fTree->GetMaxEntryLoop();

   // invoke the selector
   Long64_t nrows = Process(fSelector,option,nentries,firstentry);
   fSelectedRows = nrows;
   fDimension = fSelector->GetDimension();

   //*-* an Event List
   if (fDimension <= 0) {
      fTree->SetEstimate(oldEstimate);
      if (fSelector->GetCleanElist()) {
         // We are in the case where the input list was reset!
         fTree->SetEntryList(elist);
         delete fSelector->GetObject();
      }
      return nrows;
   }

   // Draw generated histogram
   Long64_t drawflag = fSelector->GetDrawFlag();
   Int_t action   = fSelector->GetAction();
   Bool_t draw = kFALSE;
   if (!drawflag && !opt.Contains("goff")) draw = kTRUE;
   if (!optcandle && !optpara) fHistogram = (TH1*)fSelector->GetObject();
   if (optnorm) {
      Double_t sumh= fHistogram->GetSumOfWeights();
      if (sumh != 0) fHistogram->Scale(1./sumh);
   }

   if (drawflag) {
      if (gPad) {
         if (!opt.Contains("same") && !opt.Contains("goff")) {
            gPad->DrawFrame(-1.,-1.,1.,1.);
            TText *text_empty = new TText(0.,0.,"Empty");
            text_empty->SetTextAlign(22);
            text_empty->SetTextFont(42);
            text_empty->SetTextSize(0.1);
            text_empty->SetTextColor(1);
            text_empty->Draw();
         }
      } else {
         Warning("DrawSelect", "The selected TTree subset is empty.");
      }
   }

   //*-*- 1-D distribution
   if (fDimension == 1 && !(optpara||optcandle)) {
      if (fSelector->GetVar1()->IsInteger()) fHistogram->LabelsDeflate("X");
      if (draw) fHistogram->Draw(opt.Data());

   //*-*- 2-D distribution
   } else if (fDimension == 2 && !(optpara||optcandle)) {
      if (fSelector->GetVar1()->IsInteger()) fHistogram->LabelsDeflate("Y");
      if (fSelector->GetVar2()->IsInteger()) fHistogram->LabelsDeflate("X");
      if (action == 4) {
         if (draw) fHistogram->Draw(opt.Data());
      } else {
         Bool_t graph = kFALSE;
         Int_t l = opt.Length();
         if (l == 0 || opt == "same") graph = kTRUE;
         if (opt.Contains("p")     || opt.Contains("*")    || opt.Contains("l"))    graph = kTRUE;
         if (opt.Contains("surf")  || opt.Contains("lego") || opt.Contains("cont")) graph = kFALSE;
         if (opt.Contains("col")   || opt.Contains("hist") || opt.Contains("scat")) graph = kFALSE;
         if (!graph) {
            if (draw) fHistogram->Draw(opt.Data());
         } else {
            if (fSelector->GetOldHistogram() && draw) fHistogram->Draw(opt.Data());
         }
      }
   //*-*- 3-D distribution
   } else if (fDimension == 3 && !(optpara||optcandle)) {
      if (fSelector->GetVar1()->IsInteger()) fHistogram->LabelsDeflate("Z");
      if (fSelector->GetVar2()->IsInteger()) fHistogram->LabelsDeflate("Y");
      if (fSelector->GetVar3()->IsInteger()) fHistogram->LabelsDeflate("X");
      if (action == 23) {
         if (draw) fHistogram->Draw(opt.Data());
      } else if (action == 33) {
         if (draw) {
            if (opt.Contains("z")) fHistogram->Draw("func z");
            else                   fHistogram->Draw("func");
         }
      } else {
         Int_t noscat = opt.Length();
         if (opt.Contains("same")) noscat -= 4;
         if (noscat) {
            if (draw) fHistogram->Draw(opt.Data());
         } else {
            if (fSelector->GetOldHistogram() && draw) fHistogram->Draw(opt.Data());
         }
      }
   //*-*- 4-D distribution
   } else if (fDimension == 4 && !(optpara||optcandle)) {
      if (fSelector->GetVar1()->IsInteger()) fHistogram->LabelsDeflate("Z");
      if (fSelector->GetVar2()->IsInteger()) fHistogram->LabelsDeflate("Y");
      if (fSelector->GetVar3()->IsInteger()) fHistogram->LabelsDeflate("X");
      if (draw) fHistogram->Draw(opt.Data());
      Int_t ncolors  = gStyle->GetNumberOfColors();
      TObjArray *pms = (TObjArray*)fHistogram->GetListOfFunctions()->FindObject("polymarkers");
      for (Int_t col=0;col<ncolors;col++) {
         if (!pms) continue;
         TPolyMarker3D *pm3d = (TPolyMarker3D*)pms->UncheckedAt(col);
         if (draw) pm3d->Draw();
      }
   //*-*- Parallel Coordinates or Candle chart.
   } else if (fDimension > 1 && (optpara || optcandle)) {
      if (draw) {
         TObject* para = fSelector->GetObject();
         fTree->Draw(">>enlist",selection,"entrylist",nentries,firstentry);
         TObject *enlist = gDirectory->FindObject("enlist");
         gROOT->ProcessLine(Form("TParallelCoord::SetEntryList((TParallelCoord*)0x%lx,(TEntryList*)0x%lx)",
                                     (ULong_t)para, (ULong_t)enlist));
      }
   //*-*- 5d with gl
   } else if (fDimension == 5 && optgl5d) {
      gROOT->ProcessLineFast(Form("(new TGL5DDataSet((TTree *)0x%lx))->Draw(\"%s\");", (ULong_t)fTree, opt.Data()));
      gStyle->SetCanvasPreferGL(pgl);
   }

   if (fHistogram) fHistogram->SetCanExtend(TH1::kNoAxis);
   return fSelectedRows;
}

////////////////////////////////////////////////////////////////////////////////
/// Fit  a projected item(s) from a Tree.
/// Returns -1 in case of error or number of selected events in case of success.
///
///  The formula is a TF1 expression.
///
///  See TTree::Draw for explanations of the other parameters.
///
///  By default the temporary histogram created is called htemp.
///  If varexp contains >>hnew , the new histogram created is called hnew
///  and it is kept in the current directory.
///  Example:
/// ~~~{.cpp}
///    tree.Fit("pol4","sqrt(x)>>hsqrt","y>0")
///    will fit sqrt(x) and save the histogram as "hsqrt" in the current
///    directory.
/// ~~~
///
/// The function returns the status of the histogram fit (see TH1::Fit)
/// If no entries were selected, the function returns -1;
///   (i.e. fitResult is null if the fit is OK)

Int_t TTreePlayer::Fit(const char *formula ,const char *varexp, const char *selection,Option_t *option ,Option_t *goption,Long64_t nentries, Long64_t firstentry)
{
   Int_t nch = option ? strlen(option) + 10 : 10;
   char *opt = new char[nch];
   if (option) strlcpy(opt,option,nch-1);
   else        strlcpy(opt,"goff",5);

   Long64_t nsel = DrawSelect(varexp,selection,opt,nentries,firstentry);

   delete [] opt;
   Int_t fitResult = -1;

   if (fHistogram && nsel > 0) {
      fitResult = fHistogram->Fit(formula,option,goption);
   }
   return fitResult;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of entries matching the selection.
/// Return -1 in case of errors.
///
/// If the selection uses any arrays or containers, we return the number
/// of entries where at least one element match the selection.
/// GetEntries is implemented using the selector class TSelectorEntries,
/// which can be used directly (see code in TTreePlayer::GetEntries) for
/// additional option.
/// If SetEventList was used on the TTree or TChain, only that subset
/// of entries will be considered.

Long64_t TTreePlayer::GetEntries(const char *selection)
{
   TSelectorEntries s(selection);
   fTree->Process(&s);
   fTree->SetNotify(0);
   return s.GetSelectedRows();
}

////////////////////////////////////////////////////////////////////////////////
/// return the number of entries to be processed
/// this function checks that nentries is not bigger than the number
/// of entries in the Tree or in the associated TEventlist

Long64_t TTreePlayer::GetEntriesToProcess(Long64_t firstentry, Long64_t nentries) const
{
   Long64_t lastentry = firstentry + nentries - 1;
   if (lastentry > fTree->GetEntriesFriend()-1) {
      lastentry  = fTree->GetEntriesFriend() - 1;
      nentries   = lastentry - firstentry + 1;
   }
   //TEventList *elist = fTree->GetEventList();
   //if (elist && elist->GetN() < nentries) nentries = elist->GetN();
   TEntryList *elist = fTree->GetEntryList();
   if (elist && elist->GetN() < nentries) nentries = elist->GetN();
   return nentries;
}

////////////////////////////////////////////////////////////////////////////////
/// Return name corresponding to colindex in varexp.
///
/// -  varexp is a string of names separated by :
/// -  index is an array with pointers to the start of name[i] in varexp

const char *TTreePlayer::GetNameByIndex(TString &varexp, Int_t *index,Int_t colindex)
{
   TTHREAD_TLS_DECL(std::string,column);
   if (colindex<0 ) return "";
   Int_t i1,n;
   i1 = index[colindex] + 1;
   n  = index[colindex+1] - i1;
   column = varexp(i1,n).Data();
   //  return (const char*)Form((const char*)column);
   return column.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the name of the branch pointer needed by MakeClass/MakeSelector

static TString R__GetBranchPointerName(TLeaf *leaf, Bool_t replace = kTRUE)
{
   TLeaf *leafcount = leaf->GetLeafCount();
   TBranch *branch = leaf->GetBranch();

   TString branchname( branch->GetName() );

   if ( branch->GetNleaves() <= 1 ) {
       if (branch->IsA() != TBranchObject::Class()) {
         if (!leafcount) {
            TBranch *mother = branch->GetMother();
            const char* ltitle = leaf->GetTitle();
            if (mother && mother!=branch) {
               branchname = mother->GetName();
               if (branchname[branchname.Length()-1]!='.') {
                  branchname += ".";
               }
               if (strncmp(branchname.Data(),ltitle,branchname.Length())==0) {
                  branchname = "";
               }
            } else {
               branchname = "";
            }
            branchname += ltitle;
         }
      }
   }
   if (replace) {
      char *bname = (char*)branchname.Data();
      char *twodim = (char*)strstr(bname,"[");
      if (twodim) *twodim = 0;
      while (*bname) {
         if (*bname == '.') *bname='_';
         if (*bname == ',') *bname='_';
         if (*bname == ':') *bname='_';
         if (*bname == '<') *bname='_';
         if (*bname == '>') *bname='_';
         bname++;
      }
   }
   return branchname;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate skeleton analysis class for this Tree.
///
/// The following files are produced: classname.h and classname.C
/// If classname is 0, classname will be called "nameoftree.
///
/// The generated code in classname.h includes the following:
///    - Identification of the original Tree and Input file name
///    - Definition of analysis class (data and functions)
///    - the following class functions:
///       - constructor (connecting by default the Tree file)
///       - GetEntry(Long64_t entry)
///       - Init(TTree *tree) to initialize a new TTree
///       - Show(Long64_t entry) to read and Dump entry
///
/// The generated code in classname.C includes only the main
/// analysis function Loop.
///
/// To use this function:
///    - connect your Tree file (eg: TFile f("myfile.root");)
///    - T->MakeClass("MyClass");
///
/// where T is the name of the Tree in file myfile.root
/// and MyClass.h, MyClass.C the name of the files created by this function.
/// In a ROOT session, you can do:
/// ~~~{.cpp}
///    root> .L MyClass.C
///    root> MyClass t
///    root> t.GetEntry(12); // Fill t data members with entry number 12
///    root> t.Show();       // Show values of entry 12
///    root> t.Show(16);     // Read and show values of entry 16
///    root> t.Loop();       // Loop on all entries
/// ~~~
///  NOTE: Do not use the code generated for one Tree in case of a TChain.
///        Maximum dimensions calculated on the basis of one TTree only
///        might be too small when processing all the TTrees in one TChain.
///        Instead of myTree.MakeClass(..,  use myChain.MakeClass(..

Int_t TTreePlayer::MakeClass(const char *classname, const char *option)
{
   TString opt = option;
   opt.ToLower();

   // Connect output files
   if (!classname) classname = fTree->GetName();

   TString thead;
   thead.Form("%s.h", classname);
   FILE *fp = fopen(thead, "w");
   if (!fp) {
      Error("MakeClass","cannot open output file %s", thead.Data());
      return 3;
   }
   TString tcimp;
   tcimp.Form("%s.C", classname);
   FILE *fpc = fopen(tcimp, "w");
   if (!fpc) {
      Error("MakeClass","cannot open output file %s", tcimp.Data());
      fclose(fp);
      return 3;
   }
   TString treefile;
   if (fTree->GetDirectory() && fTree->GetDirectory()->GetFile()) {
      treefile = fTree->GetDirectory()->GetFile()->GetName();
   } else {
      treefile = "Memory Directory";
   }
   // In the case of a chain, the GetDirectory information usually does
   // pertain to the Chain itself but to the currently loaded tree.
   // So we can not rely on it.
   Bool_t ischain = fTree->InheritsFrom(TChain::Class());
   Bool_t isHbook = fTree->InheritsFrom("THbookTree");
   if (isHbook)
      treefile = fTree->GetTitle();

//======================Generate classname.h=====================
   // Print header
   TObjArray *leaves = fTree->GetListOfLeaves();
   Int_t nleaves = leaves ? leaves->GetEntriesFast() : 0;
   TDatime td;
   fprintf(fp,"//////////////////////////////////////////////////////////\n");
   fprintf(fp,"// This class has been automatically generated on\n");
   fprintf(fp,"// %s by ROOT version %s\n",td.AsString(),gROOT->GetVersion());
   if (!ischain) {
      fprintf(fp,"// from TTree %s/%s\n",fTree->GetName(),fTree->GetTitle());
      fprintf(fp,"// found on file: %s\n",treefile.Data());
   } else {
      fprintf(fp,"// from TChain %s/%s\n",fTree->GetName(),fTree->GetTitle());
   }
   fprintf(fp,"//////////////////////////////////////////////////////////\n");
   fprintf(fp,"\n");
   fprintf(fp,"#ifndef %s_h\n",classname);
   fprintf(fp,"#define %s_h\n",classname);
   fprintf(fp,"\n");
   fprintf(fp,"#include <TROOT.h>\n");
   fprintf(fp,"#include <TChain.h>\n");
   fprintf(fp,"#include <TFile.h>\n");
   if (isHbook) fprintf(fp,"#include <THbookFile.h>\n");
   if (opt.Contains("selector")) fprintf(fp,"#include <TSelector.h>\n");

   // See if we can add any #include about the user data.
   Int_t l;
   fprintf(fp,"\n// Header file for the classes stored in the TTree if any.\n");
   TList listOfHeaders;
   listOfHeaders.SetOwner();
   for (l=0;l<nleaves;l++) {
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      TBranch *branch = leaf->GetBranch();
      TClass *cl = TClass::GetClass(branch->GetClassName());
      if (cl && cl->IsLoaded() && !listOfHeaders.FindObject(cl->GetName())) {
         const char *declfile = cl->GetDeclFileName();
         if (declfile && declfile[0]) {
            static const char *precstl = "prec_stl/";
            static const unsigned int precstl_len = strlen(precstl);
            static const char *rootinclude = "include/";
            static const unsigned int rootinclude_len = strlen(rootinclude);
            if (strncmp(declfile,precstl,precstl_len) == 0) {
               fprintf(fp,"#include <%s>\n",declfile+precstl_len);
               listOfHeaders.Add(new TNamed(cl->GetName(),declfile+precstl_len));
            } else if (strncmp(declfile,"/usr/include/",13) == 0) {
               fprintf(fp,"#include <%s>\n",declfile+strlen("/include/c++/"));
               listOfHeaders.Add(new TNamed(cl->GetName(),declfile+strlen("/include/c++/")));
            } else if (strstr(declfile,"/include/c++/") != 0) {
               fprintf(fp,"#include <%s>\n",declfile+strlen("/include/c++/"));
               listOfHeaders.Add(new TNamed(cl->GetName(),declfile+strlen("/include/c++/")));
            } else if (strncmp(declfile,rootinclude,rootinclude_len) == 0) {
               fprintf(fp,"#include <%s>\n",declfile+rootinclude_len);
               listOfHeaders.Add(new TNamed(cl->GetName(),declfile+rootinclude_len));
            } else {
               fprintf(fp,"#include \"%s\"\n",declfile);
               listOfHeaders.Add(new TNamed(cl->GetName(),declfile));
            }
         }
      }
   }

   // First loop on all leaves to generate dimension declarations
   Int_t len, lenb;
   char blen[1024];
   char *bname;
   Int_t *leaflen = new Int_t[nleaves];
   TObjArray *leafs = new TObjArray(nleaves);
   for (l=0;l<nleaves;l++) {
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      leafs->AddAt(new TObjString(leaf->GetName()),l);
      leaflen[l] = leaf->GetMaximum();
   }
   if (ischain) {
      // In case of a chain, one must find the maximum dimension of each leaf
      // One must be careful and not assume that all Trees in the chain
      // have the same leaves and in the same order!
      TChain *chain = (TChain*)fTree;
      Int_t ntrees = chain->GetNtrees();
      for (Int_t file=0;file<ntrees;file++) {
         Long64_t first = chain->GetTreeOffset()[file];
         chain->LoadTree(first);
         for (l=0;l<nleaves;l++) {
            TObjString *obj = (TObjString*)leafs->At(l);
            TLeaf *leaf = chain->GetLeaf(obj->GetName());
            if (leaf) {
               leaflen[l] = TMath::Max(leaflen[l],leaf->GetMaximum());
            }
         }
      }
      chain->LoadTree(0);
   }

   fprintf(fp,"\n");
   if (opt.Contains("selector")) {
      fprintf(fp,"class %s : public TSelector {\n",classname);
      fprintf(fp,"public :\n");
      fprintf(fp,"   TTree          *fChain;   //!pointer to the analyzed TTree or TChain\n");
   } else {
      fprintf(fp,"class %s {\n",classname);
      fprintf(fp,"public :\n");
      fprintf(fp,"   TTree          *fChain;   //!pointer to the analyzed TTree or TChain\n");
      fprintf(fp,"   Int_t           fCurrent; //!current Tree number in a TChain\n");
   }

   fprintf(fp,"\n// Fixed size dimensions of array or collections stored in the TTree if any.\n");
   leaves = fTree->GetListOfLeaves();
   for (l=0;l<nleaves;l++) {
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      strlcpy(blen,leaf->GetName(),sizeof(blen));
      bname = &blen[0];
      while (*bname) {
         if (*bname == '.') *bname='_';
         if (*bname == ',') *bname='_';
         if (*bname == ':') *bname='_';
         if (*bname == '<') *bname='_';
         if (*bname == '>') *bname='_';
         bname++;
      }
      lenb = strlen(blen);
      if (blen[lenb-1] == '_') {
         blen[lenb-1] = 0;
         len = leaflen[l];
         if (len <= 0) len = 1;
         fprintf(fp,"   static constexpr Int_t kMax%s = %d;\n",blen,len);
      }
   }
   delete [] leaflen;
   leafs->Delete();
   delete leafs;

// second loop on all leaves to generate type declarations
   fprintf(fp,"\n   // Declaration of leaf types\n");
   TLeaf *leafcount;
   TLeafObject *leafobj;
   TBranchElement *bre=0;
   const char *headOK  = "   ";
   const char *headcom = " //";
   const char *head;
   char branchname[1024];
   char aprefix[1024];
   TObjArray branches(100);
   TObjArray mustInit(100);
   TObjArray mustInitArr(100);
   mustInitArr.SetOwner(kFALSE);
   Int_t *leafStatus = new Int_t[nleaves];
   for (l=0;l<nleaves;l++) {
      Int_t kmax = 0;
      head = headOK;
      leafStatus[l] = 0;
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      len = leaf->GetLen(); if (len<=0) len = 1;
      leafcount =leaf->GetLeafCount();
      TBranch *branch = leaf->GetBranch();
      branchname[0] = 0;
      strlcpy(branchname,branch->GetName(),sizeof(branchname));
      strlcpy(aprefix,branch->GetName(),sizeof(aprefix));
      if (!branches.FindObject(branch)) branches.Add(branch);
      else leafStatus[l] = 1;
      if ( branch->GetNleaves() > 1) {
         // More than one leaf for the branch we need to distinguish them
         strlcat(branchname,".",sizeof(branchname));
         strlcat(branchname,leaf->GetTitle(),sizeof(branchname));
         if (leafcount) {
            // remove any dimension in title
            char *dim =  (char*)strstr(branchname,"["); if (dim) dim[0] = 0;
         }
      } else {
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
      }
      char *twodim = (char*)strstr(leaf->GetTitle(),"][");
      bname = branchname;
      while (*bname) {
         if (*bname == '.') *bname='_';
         if (*bname == ',') *bname='_';
         if (*bname == ':') *bname='_';
         if (*bname == '<') *bname='_';
         if (*bname == '>') *bname='_';
         bname++;
      }
      if (branch->IsA() == TBranchObject::Class()) {
         if (branch->GetListOfBranches()->GetEntriesFast()) {leafStatus[l] = 1; continue;}
         leafobj = (TLeafObject*)leaf;
         if (!leafobj->GetClass()) {leafStatus[l] = 1; head = headcom;}
         fprintf(fp,"%s%-15s *%s;\n",head,leafobj->GetTypeName(), leafobj->GetName());
         if (leafStatus[l] == 0) mustInit.Add(leafobj);
         continue;
      }
      if (leafcount) {
         len = leafcount->GetMaximum();
         if (len<=0) len = 1;
         strlcpy(blen,leafcount->GetName(),sizeof(blen));
         bname = &blen[0];
         while (*bname) {
            if (*bname == '.') *bname='_';
            if (*bname == ',') *bname='_';
            if (*bname == ':') *bname='_';
            if (*bname == '<') *bname='_';
            if (*bname == '>') *bname='_';
            bname++;
         }
         lenb = strlen(blen);
         if (blen[lenb-1] == '_') {blen[lenb-1] = 0; kmax = 1;}
         else                     snprintf(blen,sizeof(blen),"%d",len);
      }
      if (branch->IsA() == TBranchElement::Class()) {
         bre = (TBranchElement*)branch;
         if (bre->GetType() != 3 && bre->GetType() != 4
             && bre->GetStreamerType() <= 0 && bre->GetListOfBranches()->GetEntriesFast()) {
            leafStatus[l] = 0;
         }
         if (bre->GetType() == 3 || bre->GetType() == 4) {
            fprintf(fp,"   %-15s %s_;\n","Int_t", branchname);
            continue;
         }
         if (bre->IsBranchFolder()) {
            fprintf(fp,"   %-15s *%s;\n",bre->GetClassName(), branchname);
            mustInit.Add(bre);
            continue;
         } else {
            if (branch->GetListOfBranches()->GetEntriesFast()) {leafStatus[l] = 1;}
         }
         if (bre->GetStreamerType() < 0) {
            if (branch->GetListOfBranches()->GetEntriesFast()) {
               fprintf(fp,"%s%-15s *%s;\n",headcom,bre->GetClassName(), branchname);
            } else {
               fprintf(fp,"%s%-15s *%s;\n",head,bre->GetClassName(), branchname);
               mustInit.Add(bre);
            }
            continue;
         }
         if (bre->GetStreamerType() == 0) {
            if (!TClass::GetClass(bre->GetClassName())->HasInterpreterInfo()) {leafStatus[l] = 1; head = headcom;}
            fprintf(fp,"%s%-15s *%s;\n",head,bre->GetClassName(), branchname);
            if (leafStatus[l] == 0) mustInit.Add(bre);
            continue;
         }
         if (bre->GetStreamerType() > 60) {
            TClass *cle = TClass::GetClass(bre->GetClassName());
            if (!cle) {leafStatus[l] = 1; continue;}
            if (bre->GetStreamerType() == 66) leafStatus[l] = 0;
            char brename[256];
            strlcpy(brename,bre->GetName(),255);
            char *bren = brename;
            char *adot = strrchr(bren,'.');
            if (adot) bren = adot+1;
            char *brack = strchr(bren,'[');
            if (brack) *brack = 0;
            TStreamerElement *elem = (TStreamerElement*)cle->GetStreamerInfo()->GetElements()->FindObject(bren);
            if (elem) {
               if (elem->IsA() == TStreamerBase::Class()) {leafStatus[l] = 1; continue;}
               if (!TClass::GetClass(elem->GetTypeName())) {leafStatus[l] = 1; continue;}
               if (!TClass::GetClass(elem->GetTypeName())->HasInterpreterInfo()) {leafStatus[l] = 1; head = headcom;}
               if (leafcount) fprintf(fp,"%s%-15s %s[kMax%s];\n",head,elem->GetTypeName(), branchname,blen);
               else           fprintf(fp,"%s%-15s %s;\n",head,elem->GetTypeName(), branchname);
            } else {
               if (!TClass::GetClass(bre->GetClassName())->HasInterpreterInfo()) {leafStatus[l] = 1; head = headcom;}
               fprintf(fp,"%s%-15s %s;\n",head,bre->GetClassName(), branchname);
            }
            continue;
         }
      }
      if (strlen(leaf->GetTypeName()) == 0) {leafStatus[l] = 1; continue;}
      if (leafcount) {
         //len = leafcount->GetMaximum();
         //strlcpy(blen,leafcount->GetName(),sizeof(blen));
         //bname = &blen[0];
         //while (*bname) {if (*bname == '.') *bname='_'; bname++;}
         //lenb = strlen(blen);
         //Int_t kmax = 0;
         //if (blen[lenb-1] == '_') {blen[lenb-1] = 0; kmax = 1;}
         //else                     sprintf(blen,"%d",len);

         const char *stars = " ";
         if (bre && bre->GetBranchCount2()) {
            stars = "*";
         }
         // Dimensions can be in the branchname for a split Object with a fix length C array.
         // Theses dimensions HAVE TO be placed after the dimension explicited by leafcount
         TString dimensions;
         char *dimInName = (char*) strstr(branchname,"[");
         if ( twodim || dimInName ) {
            if (dimInName) {
               dimensions = dimInName;
               dimInName[0] = 0; // terminate branchname before the array dimensions.
            }
            if (twodim) dimensions += (char*)(twodim+1);
         }
         const char* leafcountName = leafcount->GetName();
         char b2len[1024];
         if (bre && bre->GetBranchCount2()) {
            TLeaf * l2 = (TLeaf*)bre->GetBranchCount2()->GetListOfLeaves()->At(0);
            strlcpy(b2len,l2->GetName(),sizeof(b2len));
            bname = &b2len[0];
            while (*bname) {
               if (*bname == '.') *bname='_';
               if (*bname == ',') *bname='_';
               if (*bname == ':') *bname='_';
               if (*bname == '<') *bname='_';
               if (*bname == '>') *bname='_';
               bname++;
            }
            leafcountName = b2len;
         }
         if (dimensions.Length()) {
            if (kmax) fprintf(fp,"   %-14s %s%s[kMax%s]%s;   //[%s]\n",leaf->GetTypeName(), stars,
                              branchname,blen,dimensions.Data(),leafcountName);
            else      fprintf(fp,"   %-14s %s%s[%d]%s;   //[%s]\n",leaf->GetTypeName(), stars,
                              branchname,len,dimensions.Data(),leafcountName);
         } else {
            if (kmax) fprintf(fp,"   %-14s %s%s[kMax%s];   //[%s]\n",leaf->GetTypeName(), stars, branchname,blen,leafcountName);
            else      fprintf(fp,"   %-14s %s%s[%d];   //[%s]\n",leaf->GetTypeName(), stars, branchname,len,leafcountName);
         }
         if (stars[0]=='*') {
            TNamed *n;
            if (kmax) n = new TNamed(branchname, Form("kMax%s",blen));
            else n = new TNamed(branchname, Form("%d",len));
            mustInitArr.Add(n);
         }
      } else {
         if (strstr(branchname,"[")) len = 1;
         if (len < 2) fprintf(fp,"   %-15s %s;\n",leaf->GetTypeName(), branchname);
         else {
            if (twodim) fprintf(fp,"   %-15s %s%s;\n",leaf->GetTypeName(), branchname,(char*)strstr(leaf->GetTitle(),"["));
            else        fprintf(fp,"   %-15s %s[%d];\n",leaf->GetTypeName(), branchname,len);
         }
      }
   }

// generate list of branches
   fprintf(fp,"\n");
   fprintf(fp,"   // List of branches\n");
   for (l=0;l<nleaves;l++) {
      if (leafStatus[l]) continue;
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      fprintf(fp,"   TBranch        *b_%s;   //!\n",R__GetBranchPointerName(leaf).Data());
   }

// generate class member functions prototypes
   if (opt.Contains("selector")) {
      fprintf(fp,"\n");
      fprintf(fp,"   %s(TTree * /*tree*/ =0) : fChain(0) { }\n",classname) ;
      fprintf(fp,"   virtual ~%s() { }\n",classname);
      fprintf(fp,"   virtual Int_t   Version() const { return 2; }\n");
      fprintf(fp,"   virtual void    Begin(TTree *tree);\n");
      fprintf(fp,"   virtual void    SlaveBegin(TTree *tree);\n");
      fprintf(fp,"   virtual void    Init(TTree *tree);\n");
      fprintf(fp,"   virtual Bool_t  Notify();\n");
      fprintf(fp,"   virtual Bool_t  Process(Long64_t entry);\n");
      fprintf(fp,"   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }\n");
      fprintf(fp,"   virtual void    SetOption(const char *option) { fOption = option; }\n");
      fprintf(fp,"   virtual void    SetObject(TObject *obj) { fObject = obj; }\n");
      fprintf(fp,"   virtual void    SetInputList(TList *input) { fInput = input; }\n");
      fprintf(fp,"   virtual TList  *GetOutputList() const { return fOutput; }\n");
      fprintf(fp,"   virtual void    SlaveTerminate();\n");
      fprintf(fp,"   virtual void    Terminate();\n\n");
      fprintf(fp,"   ClassDef(%s,0);\n",classname);
      fprintf(fp,"};\n");
      fprintf(fp,"\n");
      fprintf(fp,"#endif\n");
      fprintf(fp,"\n");
   } else {
      fprintf(fp,"\n");
      fprintf(fp,"   %s(TTree *tree=0);\n",classname);
      fprintf(fp,"   virtual ~%s();\n",classname);
      fprintf(fp,"   virtual Int_t    Cut(Long64_t entry);\n");
      fprintf(fp,"   virtual Int_t    GetEntry(Long64_t entry);\n");
      fprintf(fp,"   virtual Long64_t LoadTree(Long64_t entry);\n");
      fprintf(fp,"   virtual void     Init(TTree *tree);\n");
      fprintf(fp,"   virtual void     Loop();\n");
      fprintf(fp,"   virtual Bool_t   Notify();\n");
      fprintf(fp,"   virtual void     Show(Long64_t entry = -1);\n");
      fprintf(fp,"};\n");
      fprintf(fp,"\n");
      fprintf(fp,"#endif\n");
      fprintf(fp,"\n");
   }
// generate code for class constructor
   fprintf(fp,"#ifdef %s_cxx\n",classname);
   if (!opt.Contains("selector")) {
      fprintf(fp,"%s::%s(TTree *tree) : fChain(0) \n",classname,classname);
      fprintf(fp,"{\n");
      fprintf(fp,"// if parameter tree is not specified (or zero), connect the file\n");
      fprintf(fp,"// used to generate this class and read the Tree.\n");
      fprintf(fp,"   if (tree == 0) {\n");
      if (ischain) {
         fprintf(fp,"\n#ifdef SINGLE_TREE\n");
         fprintf(fp,"      // The following code should be used if you want this class to access\n");
         fprintf(fp,"      // a single tree instead of a chain\n");
      }
      if (isHbook) {
         fprintf(fp,"      THbookFile *f = (THbookFile*)gROOT->GetListOfBrowsables()->FindObject(\"%s\");\n",
                    treefile.Data());
         fprintf(fp,"      if (!f) {\n");
         fprintf(fp,"         f = new THbookFile(\"%s\");\n",treefile.Data());
         fprintf(fp,"      }\n");
         Int_t hid;
         sscanf(fTree->GetName(),"h%d",&hid);
         fprintf(fp,"      tree = (TTree*)f->Get(%d);\n\n",hid);
      } else {
         fprintf(fp,"      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(\"%s\");\n",treefile.Data());
         fprintf(fp,"      if (!f || !f->IsOpen()) {\n");
         fprintf(fp,"         f = new TFile(\"%s\");\n",treefile.Data());
         fprintf(fp,"      }\n");
         if (fTree->GetDirectory() != fTree->GetCurrentFile()) {
            fprintf(fp,"      TDirectory * dir = (TDirectory*)f->Get(\"%s\");\n",fTree->GetDirectory()->GetPath());
            fprintf(fp,"      dir->GetObject(\"%s\",tree);\n\n",fTree->GetName());
         } else {
            fprintf(fp,"      f->GetObject(\"%s\",tree);\n\n",fTree->GetName());
         }
      }
      if (ischain) {
         fprintf(fp,"#else // SINGLE_TREE\n\n");
         fprintf(fp,"      // The following code should be used if you want this class to access a chain\n");
         fprintf(fp,"      // of trees.\n");
         fprintf(fp,"      TChain * chain = new TChain(\"%s\",\"%s\");\n",
                 fTree->GetName(),fTree->GetTitle());
         {
            R__LOCKGUARD(gROOTMutex);
            TIter next(((TChain*)fTree)->GetListOfFiles());
            TChainElement *element;
            while ((element = (TChainElement*)next())) {
               fprintf(fp,"      chain->Add(\"%s/%s\");\n",element->GetTitle(),element->GetName());
            }
         }
         fprintf(fp,"      tree = chain;\n");
         fprintf(fp,"#endif // SINGLE_TREE\n\n");
      }
      fprintf(fp,"   }\n");
      fprintf(fp,"   Init(tree);\n");
      fprintf(fp,"}\n");
      fprintf(fp,"\n");
   }

// generate code for class destructor()
   if (!opt.Contains("selector")) {
      fprintf(fp,"%s::~%s()\n",classname,classname);
      fprintf(fp,"{\n");
      fprintf(fp,"   if (!fChain) return;\n");
      if (isHbook) {
         //fprintf(fp,"   delete fChain->GetCurrentFile();\n");
      } else {
         fprintf(fp,"   delete fChain->GetCurrentFile();\n");
      }
      fprintf(fp,"}\n");
      fprintf(fp,"\n");
   }
// generate code for class member function GetEntry()
   if (!opt.Contains("selector")) {
      fprintf(fp,"Int_t %s::GetEntry(Long64_t entry)\n",classname);
      fprintf(fp,"{\n");
      fprintf(fp,"// Read contents of entry.\n");

      fprintf(fp,"   if (!fChain) return 0;\n");
      fprintf(fp,"   return fChain->GetEntry(entry);\n");
      fprintf(fp,"}\n");
   }
// generate code for class member function LoadTree()
   if (!opt.Contains("selector")) {
      fprintf(fp,"Long64_t %s::LoadTree(Long64_t entry)\n",classname);
      fprintf(fp,"{\n");
      fprintf(fp,"// Set the environment to read one entry\n");
      fprintf(fp,"   if (!fChain) return -5;\n");
      fprintf(fp,"   Long64_t centry = fChain->LoadTree(entry);\n");
      fprintf(fp,"   if (centry < 0) return centry;\n");
      fprintf(fp,"   if (fChain->GetTreeNumber() != fCurrent) {\n");
      fprintf(fp,"      fCurrent = fChain->GetTreeNumber();\n");
      fprintf(fp,"      Notify();\n");
      fprintf(fp,"   }\n");
      fprintf(fp,"   return centry;\n");
      fprintf(fp,"}\n");
      fprintf(fp,"\n");
   }

// generate code for class member function Init(), first pass = get branch pointer
   fprintf(fp,"void %s::Init(TTree *tree)\n",classname);
   fprintf(fp,"{\n");
   fprintf(fp,"   // The Init() function is called when the selector needs to initialize\n"
              "   // a new tree or chain. Typically here the branch addresses and branch\n"
              "   // pointers of the tree will be set.\n"
              "   // It is normally not necessary to make changes to the generated\n"
              "   // code, but the routine can be extended by the user if needed.\n"
              "   // Init() will be called many times when running on PROOF\n"
              "   // (once per file to be processed).\n\n");
   if (mustInit.Last()) {
      TIter next(&mustInit);
      TObject *obj;
      fprintf(fp,"   // Set object pointer\n");
      while( (obj = next()) ) {
         if (obj->InheritsFrom(TBranch::Class())) {
            strlcpy(branchname,((TBranch*)obj)->GetName(),sizeof(branchname));
         } else if (obj->InheritsFrom(TLeaf::Class())) {
            strlcpy(branchname,((TLeaf*)obj)->GetName(),sizeof(branchname));
         }
         branchname[1023]=0;
         bname = branchname;
         while (*bname) {
            if (*bname == '.') *bname='_';
            if (*bname == ',') *bname='_';
            if (*bname == ':') *bname='_';
            if (*bname == '<') *bname='_';
            if (*bname == '>') *bname='_';
            bname++;
         }
         fprintf(fp,"   %s = 0;\n",branchname );
      }
   }
   if (mustInitArr.Last()) {
      TIter next(&mustInitArr);
      TNamed *info;
      fprintf(fp,"   // Set array pointer\n");
      while( (info = (TNamed*)next()) ) {
         fprintf(fp,"   for(int i=0; i<%s; ++i) %s[i] = 0;\n",info->GetTitle(),info->GetName());
      }
      fprintf(fp,"\n");
   }
   fprintf(fp,"   // Set branch addresses and branch pointers\n");
   fprintf(fp,"   if (!tree) return;\n");
   fprintf(fp,"   fChain = tree;\n");
   if (!opt.Contains("selector")) fprintf(fp,"   fCurrent = -1;\n");
   fprintf(fp,"   fChain->SetMakeClass(1);\n");
   fprintf(fp,"\n");
   for (l=0;l<nleaves;l++) {
      if (leafStatus[l]) continue;
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      len = leaf->GetLen();
      leafcount =leaf->GetLeafCount();
      TBranch *branch = leaf->GetBranch();
      strlcpy(aprefix,branch->GetName(),sizeof(aprefix));

      if ( branch->GetNleaves() > 1) {
         // More than one leaf for the branch we need to distinguish them
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
         strlcat(branchname,".",sizeof(branchname));
         strlcat(branchname,leaf->GetTitle(),sizeof(branchname));
         if (leafcount) {
            // remove any dimension in title
            char *dim =  (char*)strstr(branchname,"["); if (dim) dim[0] = 0;
         }
      } else {
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
         if (branch->IsA() == TBranchElement::Class()) {
            bre = (TBranchElement*)branch;
            if (bre->GetType() == 3 || bre->GetType()==4) strlcat(branchname,"_",sizeof(branchname));
         }
      }
      bname = branchname;
      char *brak = strstr(branchname,"[");     if (brak) *brak = 0;
      char *twodim = (char*)strstr(bname,"["); if (twodim) *twodim = 0;
      while (*bname) {
         if (*bname == '.') *bname='_';
         if (*bname == ',') *bname='_';
         if (*bname == ':') *bname='_';
         if (*bname == '<') *bname='_';
         if (*bname == '>') *bname='_';
         bname++;
      }
      const char *maybedisable = "";
      if (branch != fTree->GetBranch(branch->GetName())) {
         Error("MakeClass","The branch named %s (full path name: %s) is hidden by another branch of the same name and its data will not be loaded.",branch->GetName(),R__GetBranchPointerName(leaf,kFALSE).Data());
         maybedisable = "// ";
      }
      if (branch->IsA() == TBranchObject::Class()) {
         if (branch->GetListOfBranches()->GetEntriesFast()) {
            fprintf(fp,"%s   fChain->SetBranchAddress(\"%s\",(void*)-1,&b_%s);\n",maybedisable,branch->GetName(),R__GetBranchPointerName(leaf).Data());
            continue;
         }
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
      }
      if (branch->IsA() == TBranchElement::Class()) {
         if (((TBranchElement*)branch)->GetType() == 3) len =1;
         if (((TBranchElement*)branch)->GetType() == 4) len =1;
      }
      if (leafcount) len = leafcount->GetMaximum()+1;
      if (len > 1) fprintf(fp,"%s   fChain->SetBranchAddress(\"%s\", %s, &b_%s);\n",
                           maybedisable,branch->GetName(), branchname, R__GetBranchPointerName(leaf).Data());
      else         fprintf(fp,"%s   fChain->SetBranchAddress(\"%s\", &%s, &b_%s);\n",
                           maybedisable,branch->GetName(), branchname, R__GetBranchPointerName(leaf).Data());
   }
   //must call Notify in case of MakeClass
   if (!opt.Contains("selector")) {
      fprintf(fp,"   Notify();\n");
   }

   fprintf(fp,"}\n");
   fprintf(fp,"\n");

// generate code for class member function Notify()
   fprintf(fp,"Bool_t %s::Notify()\n",classname);
   fprintf(fp,"{\n");
   fprintf(fp,"   // The Notify() function is called when a new file is opened. This\n"
              "   // can be either for a new TTree in a TChain or when when a new TTree\n"
              "   // is started when using PROOF. It is normally not necessary to make changes\n"
              "   // to the generated code, but the routine can be extended by the\n"
              "   // user if needed. The return value is currently not used.\n\n");
   fprintf(fp,"   return kTRUE;\n");
   fprintf(fp,"}\n");
   fprintf(fp,"\n");

// generate code for class member function Show()
   if (!opt.Contains("selector")) {
      fprintf(fp,"void %s::Show(Long64_t entry)\n",classname);
      fprintf(fp,"{\n");
      fprintf(fp,"// Print contents of entry.\n");
      fprintf(fp,"// If entry is not specified, print current entry\n");

      fprintf(fp,"   if (!fChain) return;\n");
      fprintf(fp,"   fChain->Show(entry);\n");
      fprintf(fp,"}\n");
   }
// generate code for class member function Cut()
   if (!opt.Contains("selector")) {
      fprintf(fp,"Int_t %s::Cut(Long64_t entry)\n",classname);
      fprintf(fp,"{\n");
      fprintf(fp,"// This function may be called from Loop.\n");
      fprintf(fp,"// returns  1 if entry is accepted.\n");
      fprintf(fp,"// returns -1 otherwise.\n");

      fprintf(fp,"   return 1;\n");
      fprintf(fp,"}\n");
   }
   fprintf(fp,"#endif // #ifdef %s_cxx\n",classname);

//======================Generate classname.C=====================
   if (!opt.Contains("selector")) {
      // generate code for class member function Loop()
      fprintf(fpc,"#define %s_cxx\n",classname);
      fprintf(fpc,"#include \"%s\"\n",thead.Data());
      fprintf(fpc,"#include <TH2.h>\n");
      fprintf(fpc,"#include <TStyle.h>\n");
      fprintf(fpc,"#include <TCanvas.h>\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"void %s::Loop()\n",classname);
      fprintf(fpc,"{\n");
      fprintf(fpc,"//   In a ROOT session, you can do:\n");
      fprintf(fpc,"//      root> .L %s.C\n",classname);
      fprintf(fpc,"//      root> %s t\n",classname);
      fprintf(fpc,"//      root> t.GetEntry(12); // Fill t data members with entry number 12\n");
      fprintf(fpc,"//      root> t.Show();       // Show values of entry 12\n");
      fprintf(fpc,"//      root> t.Show(16);     // Read and show values of entry 16\n");
      fprintf(fpc,"//      root> t.Loop();       // Loop on all entries\n");
      fprintf(fpc,"//\n");
      fprintf(fpc,"\n//     This is the loop skeleton where:\n");
      fprintf(fpc,"//    jentry is the global entry number in the chain\n");
      fprintf(fpc,"//    ientry is the entry number in the current Tree\n");
      fprintf(fpc,"//  Note that the argument to GetEntry must be:\n");
      fprintf(fpc,"//    jentry for TChain::GetEntry\n");
      fprintf(fpc,"//    ientry for TTree::GetEntry and TBranch::GetEntry\n");
      fprintf(fpc,"//\n");
      fprintf(fpc,"//       To read only selected branches, Insert statements like:\n");
      fprintf(fpc,"// METHOD1:\n");
      fprintf(fpc,"//    fChain->SetBranchStatus(\"*\",0);  // disable all branches\n");
      fprintf(fpc,"//    fChain->SetBranchStatus(\"branchname\",1);  // activate branchname\n");
      fprintf(fpc,"// METHOD2: replace line\n");
      fprintf(fpc,"//    fChain->GetEntry(jentry);       //read all branches\n");
      fprintf(fpc,"//by  b_branchname->GetEntry(ientry); //read only this branch\n");
      fprintf(fpc,"   if (fChain == 0) return;\n");
      fprintf(fpc,"\n   Long64_t nentries = fChain->GetEntriesFast();\n");
      fprintf(fpc,"\n   Long64_t nbytes = 0, nb = 0;\n");
      fprintf(fpc,"   for (Long64_t jentry=0; jentry<nentries;jentry++) {\n");
      fprintf(fpc,"      Long64_t ientry = LoadTree(jentry);\n");
      fprintf(fpc,"      if (ientry < 0) break;\n");
      fprintf(fpc,"      nb = fChain->GetEntry(jentry);   nbytes += nb;\n");
      fprintf(fpc,"      // if (Cut(ientry) < 0) continue;\n");
      fprintf(fpc,"   }\n");
      fprintf(fpc,"}\n");
   }
   if (opt.Contains("selector")) {
      // generate usage comments and list of includes
      fprintf(fpc,"#define %s_cxx\n",classname);
      fprintf(fpc,"// The class definition in %s.h has been generated automatically\n",classname);
      fprintf(fpc,"// by the ROOT utility TTree::MakeSelector(). This class is derived\n");
      fprintf(fpc,"// from the ROOT class TSelector. For more information on the TSelector\n"
                  "// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.\n\n");
      fprintf(fpc,"// The following methods are defined in this file:\n");
      fprintf(fpc,"//    Begin():        called every time a loop on the tree starts,\n");
      fprintf(fpc,"//                    a convenient place to create your histograms.\n");
      fprintf(fpc,"//    SlaveBegin():   called after Begin(), when on PROOF called only on the\n"
                  "//                    slave servers.\n");
      fprintf(fpc,"//    Process():      called for each event, in this function you decide what\n");
      fprintf(fpc,"//                    to read and fill your histograms.\n");
      fprintf(fpc,"//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF\n"
                  "//                    called only on the slave servers.\n");
      fprintf(fpc,"//    Terminate():    called at the end of the loop on the tree,\n");
      fprintf(fpc,"//                    a convenient place to draw/fit your histograms.\n");
      fprintf(fpc,"//\n");
      fprintf(fpc,"// To use this file, try the following session on your Tree T:\n");
      fprintf(fpc,"//\n");
      fprintf(fpc,"// root> T->Process(\"%s.C\")\n",classname);
      fprintf(fpc,"// root> T->Process(\"%s.C\",\"some options\")\n",classname);
      fprintf(fpc,"// root> T->Process(\"%s.C+\")\n",classname);
      fprintf(fpc,"//\n\n");
      fprintf(fpc,"#include \"%s\"\n",thead.Data());
      fprintf(fpc,"#include <TH2.h>\n");
      fprintf(fpc,"#include <TStyle.h>\n");
      fprintf(fpc,"\n");
      // generate code for class member function Begin
      fprintf(fpc,"\n");
      fprintf(fpc,"void %s::Begin(TTree * /*tree*/)\n",classname);
      fprintf(fpc,"{\n");
      fprintf(fpc,"   // The Begin() function is called at the start of the query.\n");
      fprintf(fpc,"   // When running with PROOF Begin() is only called on the client.\n");
      fprintf(fpc,"   // The tree argument is deprecated (on PROOF 0 is passed).\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"   TString option = GetOption();\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"}\n");
      // generate code for class member function SlaveBegin
      fprintf(fpc,"\n");
      fprintf(fpc,"void %s::SlaveBegin(TTree * /*tree*/)\n",classname);
      fprintf(fpc,"{\n");
      fprintf(fpc,"   // The SlaveBegin() function is called after the Begin() function.\n");
      fprintf(fpc,"   // When running with PROOF SlaveBegin() is called on each slave server.\n");
      fprintf(fpc,"   // The tree argument is deprecated (on PROOF 0 is passed).\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"   TString option = GetOption();\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"}\n");
      // generate code for class member function Process
      fprintf(fpc,"\n");
      fprintf(fpc,"Bool_t %s::Process(Long64_t entry)\n",classname);
      fprintf(fpc,"{\n");
      fprintf(fpc,"   // The Process() function is called for each entry in the tree (or possibly\n"
                  "   // keyed object in the case of PROOF) to be processed. The entry argument\n"
                  "   // specifies which entry in the currently loaded tree is to be processed.\n"
                  "   // It can be passed to either %s::GetEntry() or TBranch::GetEntry()\n"
                  "   // to read either all or the required parts of the data. When processing\n"
                  "   // keyed objects with PROOF, the object is already loaded and is available\n"
                  "   // via the fObject pointer.\n"
                  "   //\n"
                  "   // This function should contain the \"body\" of the analysis. It can contain\n"
                  "   // simple or elaborate selection criteria, run algorithms on the data\n"
                  "   // of the event and typically fill histograms.\n"
                  "   //\n"
                  "   // The processing can be stopped by calling Abort().\n"
                  "   //\n"
                  "   // Use fStatus to set the return value of TTree::Process().\n"
                  "   //\n"
                  "   // The return value is currently not used.\n\n", classname);
      fprintf(fpc,"\n");
      fprintf(fpc,"   return kTRUE;\n");
      fprintf(fpc,"}\n");
      // generate code for class member function SlaveTerminate
      fprintf(fpc,"\n");
      fprintf(fpc,"void %s::SlaveTerminate()\n",classname);
      fprintf(fpc,"{\n");
      fprintf(fpc,"   // The SlaveTerminate() function is called after all entries or objects\n"
                  "   // have been processed. When running with PROOF SlaveTerminate() is called\n"
                  "   // on each slave server.");
      fprintf(fpc,"\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"}\n");
      // generate code for class member function Terminate
      fprintf(fpc,"\n");
      fprintf(fpc,"void %s::Terminate()\n",classname);
      fprintf(fpc,"{\n");
      fprintf(fpc,"   // The Terminate() function is the last function to be called during\n"
                  "   // a query. It always runs on the client, it can be used to present\n"
                  "   // the results graphically or save the results to file.");
      fprintf(fpc,"\n");
      fprintf(fpc,"\n");
      fprintf(fpc,"}\n");
   }
   Info("MakeClass","Files: %s and %s generated from TTree: %s",thead.Data(),tcimp.Data(),fTree->GetName());
   delete [] leafStatus;
   fclose(fp);
   fclose(fpc);

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate skeleton function for this Tree
///
/// The function code is written on filename.
/// If filename is 0, filename will be called nameoftree.C
///
/// The generated code includes the following:
///    - Identification of the original Tree and Input file name
///    - Connection of the Tree file
///    - Declaration of Tree variables
///    - Setting of branches addresses
///    - A skeleton for the entry loop
///
/// To use this function:
///    - connect your Tree file (eg: TFile f("myfile.root");)
///    - T->MakeCode("anal.C");
/// where T is the name of the Tree in file myfile.root
/// and anal.C the name of the file created by this function.
///
/// NOTE: Since the implementation of this function, a new and better
///       function TTree::MakeClass() has been developed.

Int_t TTreePlayer::MakeCode(const char *filename)
{
// Connect output file
   TString tfile;
   if (filename)
      tfile = filename;
   else
      tfile.Form("%s.C", fTree->GetName());
   FILE *fp = fopen(tfile, "w");
   if (!fp) {
      Error("MakeCode","cannot open output file %s", tfile.Data());
      return 3;
   }
   TString treefile;
   if (fTree->GetDirectory() && fTree->GetDirectory()->GetFile()) {
      treefile = fTree->GetDirectory()->GetFile()->GetName();
   } else {
      treefile = "Memory Directory";
   }
   // In the case of a chain, the GetDirectory information usually does
   // pertain to the Chain itself but to the currently loaded tree.
   // So we can not rely on it.
   Bool_t ischain = fTree->InheritsFrom(TChain::Class());

// Print header
   TObjArray *leaves = fTree->GetListOfLeaves();
   Int_t nleaves = leaves ? leaves->GetEntriesFast() : 0;
   TDatime td;
   fprintf(fp,"{\n");
   fprintf(fp,"//////////////////////////////////////////////////////////\n");
   fprintf(fp,"//   This file has been automatically generated \n");
   fprintf(fp,"//     (%s by ROOT version%s)\n",td.AsString(),gROOT->GetVersion());
   if (!ischain) {
      fprintf(fp,"//   from TTree %s/%s\n",fTree->GetName(),fTree->GetTitle());
      fprintf(fp,"//   found on file: %s\n",treefile.Data());
   } else {
      fprintf(fp,"//   from TChain %s/%s\n",fTree->GetName(),fTree->GetTitle());
   }
   fprintf(fp,"//////////////////////////////////////////////////////////\n");
   fprintf(fp,"\n");
   fprintf(fp,"\n");


// Reset and file connect
   fprintf(fp,"//Reset ROOT and connect tree file\n");
   fprintf(fp,"   gROOT->Reset();\n");
   if (ischain) {
      fprintf(fp,"\n#ifdef SINGLE_TREE\n");
      fprintf(fp,"   // The following code should be used if you want this code to access\n");
      fprintf(fp,"   // a single tree instead of a chain\n");
   }
   fprintf(fp,"   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(\"%s\");\n",treefile.Data());
   fprintf(fp,"   if (!f) {\n");
   fprintf(fp,"      f = new TFile(\"%s\");\n",treefile.Data());
   fprintf(fp,"   }\n");
   if (fTree->GetDirectory() != fTree->GetCurrentFile()) {
      fprintf(fp,"    TDirectory * dir = (TDirectory*)f->Get(\"%s\");\n",fTree->GetDirectory()->GetPath());
      fprintf(fp,"    dir->GetObject(\"%s\",tree);\n\n",fTree->GetName());
   } else {
      fprintf(fp,"    f->GetObject(\"%s\",tree);\n\n",fTree->GetName());
   }
   if (ischain) {
      fprintf(fp,"#else // SINGLE_TREE\n\n");
      fprintf(fp,"   // The following code should be used if you want this code to access a chain\n");
      fprintf(fp,"   // of trees.\n");
      fprintf(fp,"   TChain *%s = new TChain(\"%s\",\"%s\");\n",
                 fTree->GetName(),fTree->GetName(),fTree->GetTitle());
      {
         R__LOCKGUARD(gROOTMutex);
         TIter next(((TChain*)fTree)->GetListOfFiles());
         TChainElement *element;
         while ((element = (TChainElement*)next())) {
            fprintf(fp,"   %s->Add(\"%s/%s\");\n",fTree->GetName(),element->GetTitle(),element->GetName());
         }
      }
      fprintf(fp,"#endif // SINGLE_TREE\n\n");
   }

// First loop on all leaves to generate type declarations
   fprintf(fp,"//Declaration of leaves types\n");
   Int_t len, l;
   TLeaf *leafcount;
   TLeafObject *leafobj;
   char *bname;
   const char *headOK  = "   ";
   const char *headcom = " //";
   const char *head;
   char branchname[1024];
   for (l=0;l<nleaves;l++) {
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      len = leaf->GetLen();
      leafcount =leaf->GetLeafCount();
      TBranch *branch = leaf->GetBranch();
      if (branch->GetListOfBranches()->GetEntriesFast() > 0) continue;

      if ( branch->GetNleaves() > 1) {
         // More than one leaf for the branch we need to distinguish them
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
         strlcat(branchname,".",sizeof(branchname));
         strlcat(branchname,leaf->GetTitle(),sizeof(branchname));
         if (leafcount) {
            // remove any dimension in title
            char *dim =  (char*)strstr(branchname,"[");
            if (dim) dim[0] = 0;
         }
      } else {
         if (leafcount) strlcpy(branchname,branch->GetName(),sizeof(branchname));
         else           strlcpy(branchname,leaf->GetTitle(),sizeof(branchname));
      }
      char *twodim = (char*)strstr(leaf->GetTitle(),"][");
      bname = branchname;
      while (*bname) {
         if (*bname == '.') *bname='_';
         if (*bname == ',') *bname='_';
         if (*bname == ':') *bname='_';
         if (*bname == '<') *bname='_';
         if (*bname == '>') *bname='_';
         bname++;
      }
      if (branch->IsA() == TBranchObject::Class()) {
         leafobj = (TLeafObject*)leaf;
         if (leafobj->GetClass()) head = headOK;
         else                     head = headcom;
         fprintf(fp,"%s%-15s *%s = 0;\n",head,leafobj->GetTypeName(), leafobj->GetName());
         continue;
      }
      if (leafcount) {
         len = leafcount->GetMaximum();
         // Dimensions can be in the branchname for a split Object with a fix length C array.
         // Theses dimensions HAVE TO be placed after the dimension explicited by leafcount
         char *dimInName = (char*) strstr(branchname,"[");
         TString dimensions;
         if ( twodim || dimInName ) {
            if (dimInName) {
               dimensions = dimInName;
               dimInName[0] = 0; // terminate branchname before the array dimensions.
            }
            if (twodim) dimensions += (char*)(twodim+1);
         }
         if (dimensions.Length()) {
            fprintf(fp,"   %-15s %s[%d]%s;\n",leaf->GetTypeName(), branchname,len,dimensions.Data());
         } else {
            fprintf(fp,"   %-15s %s[%d];\n",leaf->GetTypeName(), branchname,len);
         }
      } else {
         if (strstr(branchname,"[")) len = 1;
         if (len < 2) fprintf(fp,"   %-15s %s;\n",leaf->GetTypeName(), branchname);
         else         fprintf(fp,"   %-15s %s[%d];\n",leaf->GetTypeName(), branchname,len);
      }
   }

// Second loop on all leaves to set the corresponding branch address
   fprintf(fp,"\n   // Set branch addresses.\n");
   for (l=0;l<nleaves;l++) {
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
      len = leaf->GetLen();
      leafcount =leaf->GetLeafCount();
      TBranch *branch = leaf->GetBranch();

      if ( branch->GetNleaves() > 1) {
         // More than one leaf for the branch we need to distinguish them
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
         strlcat(branchname,".",sizeof(branchname));
         strlcat(branchname,leaf->GetTitle(),sizeof(branchname));
         if (leafcount) {
            // remove any dimension in title
            char *dim =  (char*)strstr(branchname,"[");
            if (dim) dim[0] = 0;
         }
      } else {
         if (leafcount) strlcpy(branchname,branch->GetName(),sizeof(branchname));
         else           strlcpy(branchname,leaf->GetTitle(),sizeof(branchname));
      }
      bname = branchname;
      while (*bname) {
         if (*bname == '.') *bname='_';
         if (*bname == ',') *bname='_';
         if (*bname == ':') *bname='_';
         if (*bname == '<') *bname='_';
         if (*bname == '>') *bname='_';
         bname++;
      }
      char *brak = strstr(branchname,"[");
      if (brak) *brak = 0;
      head = headOK;
      if (branch->IsA() == TBranchObject::Class()) {
         strlcpy(branchname,branch->GetName(),sizeof(branchname));
         leafobj = (TLeafObject*)leaf;
         if (!leafobj->GetClass()) head = headcom;
      }
      if (leafcount) len = leafcount->GetMaximum()+1;
      if (len > 1 || brak) fprintf(fp,"%s%s->SetBranchAddress(\"%s\",%s);\n",head,fTree->GetName(),branch->GetName(),branchname);
      else                 fprintf(fp,"%s%s->SetBranchAddress(\"%s\",&%s);\n",head,fTree->GetName(),branch->GetName(),branchname);
   }

//Generate instructions to make the loop on entries
   fprintf(fp,"\n//     This is the loop skeleton\n");
   fprintf(fp,"//       To read only selected branches, Insert statements like:\n");
   fprintf(fp,"// %s->SetBranchStatus(\"*\",0);  // disable all branches\n",fTree->GetName());
   fprintf(fp,"// %s->SetBranchStatus(\"branchname\",1);  // activate branchname\n",GetName());
   fprintf(fp,"\n   Long64_t nentries = %s->GetEntries();\n",fTree->GetName());
   fprintf(fp,"\n   Long64_t nbytes = 0;\n");
   fprintf(fp,"//   for (Long64_t i=0; i<nentries;i++) {\n");
   fprintf(fp,"//      nbytes += %s->GetEntry(i);\n",fTree->GetName());
   fprintf(fp,"//   }\n");
   fprintf(fp,"}\n");

   printf("Macro: %s generated from Tree: %s\n",tfile.Data(), fTree->GetName());
   fclose(fp);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a skeleton analysis class for this Tree using TBranchProxy.
/// TBranchProxy is the base of a class hierarchy implementing an
/// indirect access to the content of the branches of a TTree.
///
/// "proxyClassname" is expected to be of the form:
/// ~~~{.cpp}
///    [path/]fileprefix
/// ~~~
/// The skeleton will then be generated in the file:
/// ~~~{.cpp}
///    fileprefix.h
/// ~~~
/// located in the current directory or in 'path/' if it is specified.
/// The class generated will be named 'fileprefix'.
/// If the fileprefix contains a period, the right side of the period
/// will be used as the extension (instead of 'h') and the left side
/// will be used as the classname.
///
/// "macrofilename" and optionally "cutfilename" are expected to point
/// to source file which will be included in by the generated skeletong.
/// Method of the same name as the file(minus the extension and path)
/// will be called by the generated skeleton's Process method as follow:
/// ~~~{.cpp}
///    [if (cutfilename())] htemp->Fill(macrofilename());
/// ~~~
/// "option" can be used select some of the optional features during
/// the code generation.  The possible options are:
/// -  nohist : indicates that the generated ProcessFill should not
///             fill the histogram.
///
/// 'maxUnrolling' controls how deep in the class hierarchy does the
/// system 'unroll' class that are not split.  'unrolling' a class
/// will allow direct access to its data members a class (this
/// emulates the behavior of TTreeFormula).
///
/// The main features of this skeleton are:
///
/// *  on-demand loading of branches
/// *  ability to use the 'branchname' as if it was a data member
/// *  protection against array out-of-bound
/// *  ability to use the branch data as object (when the user code is available)
///
/// For example with Event.root, if
/// ~~~{.cpp}
///    Double_t somepx = fTracks.fPx[2];
/// ~~~
/// is executed by one of the method of the skeleton,
/// somepx will be updated with the current value of fPx of the 3rd track.
///
/// Both macrofilename and the optional cutfilename are expected to be
/// the name of source files which contain at least a free standing
/// function with the signature:
/// ~~~{.cpp}
///     x_t macrofilename(); // i.e function with the same name as the file
/// ~~~
/// and
/// ~~~{.cpp}
///     y_t cutfilename();   // i.e function with the same name as the file
/// ~~~
/// x_t and y_t needs to be types that can convert respectively to a double
/// and a bool (because the skeleton uses:
/// ~~~{.cpp}
///     if (cutfilename()) htemp->Fill(macrofilename());
/// ~~~
/// This 2 functions are run in a context such that the branch names are
/// available as local variables of the correct (read-only) type.
///
/// Note that if you use the same 'variable' twice, it is more efficient
/// to 'cache' the value. For example
/// ~~~{.cpp}
///   Int_t n = fEventNumber; // Read fEventNumber
///   if (n<10 || n>10) { ... }
/// ~~~
/// is more efficient than
/// ~~~{.cpp}
///   if (fEventNumber<10 || fEventNumber>10)
/// ~~~
/// Access to TClonesArray.
///
/// If a branch (or member) is a TClonesArray (let's say fTracks), you
/// can access the TClonesArray itself by using ->:
/// ~~~{.cpp}
///    fTracks->GetLast();
/// ~~~
/// However this will load the full TClonesArray object and its content.
/// To quickly read the size of the TClonesArray use (note the dot):
/// ~~~{.cpp}
///    fTracks.GetEntries();
/// ~~~
/// This will read only the size from disk if the TClonesArray has been
/// split.
/// To access the content of the TClonesArray, use the [] operator:
/// ~~~
///    float px = fTracks[i].fPx; // fPx of the i-th track
/// ~~~
/// Warning:
///
/// The variable actually use for access are 'wrapper' around the
/// real data type (to add autoload for example) and hence getting to
/// the data involves the implicit call to a C++ conversion operator.
/// This conversion is automatic in most case.  However it is not invoked
/// in a few cases, in particular in variadic function (like printf).
/// So when using printf you should either explicitly cast the value or
/// use any intermediary variable:
/// ~~~{.cpp}
///      fprintf(stdout,"trs[%d].a = %d\n",i,(int)trs.a[i]);
/// ~~~
/// Also, optionally, the generated selector will also call methods named
/// macrofilename_methodname in each of 6 main selector methods if the method
/// macrofilename_methodname exist (Where macrofilename is stripped of its
/// extension).
///
/// Concretely, with the script named h1analysisProxy.C,
///
/// -  The method         calls the method (if it exist)
/// -  Begin           -> void h1analysisProxy_Begin(TTree*);
/// -  SlaveBegin      -> void h1analysisProxy_SlaveBegin(TTree*);
/// -  Notify          -> Bool_t h1analysisProxy_Notify();
/// -  Process         -> Bool_t h1analysisProxy_Process(Long64_t);
/// -  SlaveTerminate  -> void h1analysisProxy_SlaveTerminate();
/// -  Terminate       -> void h1analysisProxy_Terminate();
///
/// If a file name macrofilename.h (or .hh, .hpp, .hxx, .hPP, .hXX) exist
/// it is included before the declaration of the proxy class.  This can
/// be used in particular to insure that the include files needed by
/// the macro file are properly loaded.
///
/// The default histogram is accessible via the variable named 'htemp'.
///
/// If the library of the classes describing the data in the branch is
/// loaded, the skeleton will add the needed `include` statements and
/// give the ability to access the object stored in the branches.
///
/// To draw px using the file `hsimple.root (generated by the
/// hsimple.C tutorial), we need a file named hsimple.cxx:
///
/// ~~~{.cpp}
///     double hsimple() {
///        return px;
///     }
/// ~~~
/// MakeProxy can then be used indirectly via the TTree::Draw interface
/// as follow:
/// ~~~{.cpp}
///     new TFile("hsimple.root")
///     ntuple->Draw("hsimple.cxx");
/// ~~~
/// A more complete example is available in the tutorials directory:
/// h1analysisProxy.cxx , h1analysProxy.h and h1analysisProxyCut.C
/// which reimplement the selector found in h1analysis.C

Int_t TTreePlayer::MakeProxy(const char *proxyClassname,
                             const char *macrofilename, const char *cutfilename,
                             const char *option, Int_t maxUnrolling)
{
   if (macrofilename==0 || strlen(macrofilename)==0 ) {
      // We currently require a file name for the script
      Error("MakeProxy","A file name for the user script is required");
      return 0;
   }

   ROOT::Internal::TTreeProxyGenerator gp(fTree,macrofilename,cutfilename,proxyClassname,option,maxUnrolling);

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate skeleton selector class for this tree.
///
/// The following files are produced: classname.h and classname.C.
/// If classname is 0, the selector will be called "nameoftree".
/// The option can be used to specify the branches that will have a data member.
///    - If option is empty, readers will be generated for each leaf.
///    - If option is "@", readers will be generated for the topmost branches.
///    - Individual branches can also be picked by their name:
///       - "X" generates readers for leaves of X.
///       - "@X" generates a reader for X as a whole.
///       - "@X;Y" generates a reader for X as a whole and also readers for the
///         leaves of Y.
///    - For further examples see the figure below.
///
/// \image html ttree_makeselector_option_examples.png
///
/// The generated code in classname.h includes the following:
///    - Identification of the original Tree and Input file name
///    - Definition of selector class (data and functions)
///    - The following class functions:
///       - constructor and destructor
///       - void    Begin(TTree *tree)
///       - void    SlaveBegin(TTree *tree)
///       - void    Init(TTree *tree)
///       - Bool_t  Notify()
///       - Bool_t  Process(Long64_t entry)
///       - void    Terminate()
///       - void    SlaveTerminate()
///
/// The selector derives from TSelector.
/// The generated code in classname.C includes empty functions defined above.
///
/// To use this function:
///    - connect your Tree file (eg: `TFile f("myfile.root");`)
///    - `T->MakeSelector("myselect");`
///       where `T` is the name of the Tree in file `myfile.root`
///       and `myselect.h`, `myselect.C` the name of the files created by this
///       function.
///
/// In a ROOT session, you can do:
/// ~~~ {.cpp}
///    root > T->Process("myselect.C")
/// ~~~
Int_t TTreePlayer::MakeReader(const char *classname, Option_t *option)
{
   if (!classname) classname = fTree->GetName();

   ROOT::Internal::TTreeReaderGenerator gsr(fTree, classname, option);

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Interface to the Principal Components Analysis class.
///
/// Create an instance of TPrincipal
/// Fill it with the selected variables
///
/// -  if option "n" is specified, the TPrincipal object is filled with
///                 normalized variables.
/// -  If option "p" is specified, compute the principal components
/// -  If option "p" and "d" print results of analysis
/// -  If option "p" and "h" generate standard histograms
/// -  If option "p" and "c" generate code of conversion functions
///
/// return a pointer to the TPrincipal object. It is the user responsibility
/// to delete this object.
///
/// The option default value is "np"
///
/// See TTreePlayer::DrawSelect for explanation of the other parameters.

TPrincipal *TTreePlayer::Principal(const char *varexp, const char *selection, Option_t *option, Long64_t nentries, Long64_t firstentry)
{
   TTreeFormula **var;
   std::vector<TString> cnames;
   TString opt = option;
   opt.ToLower();
   TPrincipal *principal = 0;
   Long64_t entry,entryNumber;
   Int_t i,nch;
   Int_t ncols = 8;   // by default first 8 columns are printed only
   TObjArray *leaves = fTree->GetListOfLeaves();
   Int_t nleaves = leaves->GetEntriesFast();
   if (nleaves < ncols) ncols = nleaves;
   nch = varexp ? strlen(varexp) : 0;

   nentries = GetEntriesToProcess(firstentry, nentries);

//*-*- Compile selection expression if there is one
   TTreeFormula *select = 0;
   if (strlen(selection)) {
      select = new TTreeFormula("Selection",selection,fTree);
      if (!select) return principal;
      if (!select->GetNdim()) { delete select; return principal; }
      fFormulaList->Add(select);
   }
//*-*- if varexp is empty, take first 8 columns by default
   int allvar = 0;
   if (varexp && !strcmp(varexp, "*")) { ncols = nleaves; allvar = 1; }
   if (nch == 0 || allvar) {
      for (i=0;i<ncols;i++) {
         cnames.push_back( ((TLeaf*)leaves->At(i))->GetName() );
      }
//*-*- otherwise select only the specified columns
   } else {
      ncols = fSelector->SplitNames(varexp,cnames);
   }
   var = new TTreeFormula* [ncols];
   Double_t *xvars = new Double_t[ncols];

//*-*- Create the TreeFormula objects corresponding to each column
   for (i=0;i<ncols;i++) {
      var[i] = new TTreeFormula("Var1",cnames[i].Data(),fTree);
      fFormulaList->Add(var[i]);
   }

//*-*- Create a TreeFormulaManager to coordinate the formulas
   TTreeFormulaManager *manager=0;
   if (fFormulaList->LastIndex()>=0) {
      manager = new TTreeFormulaManager;
      for(i=0;i<=fFormulaList->LastIndex();i++) {
         manager->Add((TTreeFormula*)fFormulaList->At(i));
      }
      manager->Sync();
   }

//*-* Build the TPrincipal object
   if (opt.Contains("n")) principal = new TPrincipal(ncols, "n");
   else                   principal = new TPrincipal(ncols);

//*-*- loop on all selected entries
   fSelectedRows = 0;
   Int_t tnumber = -1;
   for (entry=firstentry;entry<firstentry+nentries;entry++) {
      entryNumber = fTree->GetEntryNumber(entry);
      if (entryNumber < 0) break;
      Long64_t localEntry = fTree->LoadTree(entryNumber);
      if (localEntry < 0) break;
      if (tnumber != fTree->GetTreeNumber()) {
         tnumber = fTree->GetTreeNumber();
         if (manager) manager->UpdateFormulaLeaves();
      }
      int ndata = 1;
      if (manager && manager->GetMultiplicity()) {
         ndata = manager->GetNdata();
      }

      for(int inst=0;inst<ndata;inst++) {
         Bool_t loaded = kFALSE;
         if (select) {
            if (select->EvalInstance(inst) == 0) {
               continue;
            }
         }

         if (inst==0) loaded = kTRUE;
         else if (!loaded) {
            // EvalInstance(0) always needs to be called so that
            // the proper branches are loaded.
            for (i=0;i<ncols;i++) {
               var[i]->EvalInstance(0);
            }
            loaded = kTRUE;
         }

         for (i=0;i<ncols;i++) {
            xvars[i] = var[i]->EvalInstance(inst);
         }
         principal->AddRow(xvars);
      }
   }

   //*-* some actions with principal ?
   if (opt.Contains("p")) {
      principal->MakePrincipals(); // Do the actual analysis
      if (opt.Contains("d")) principal->Print();
      if (opt.Contains("h")) principal->MakeHistograms();
      if (opt.Contains("c")) principal->MakeCode();
   }

//*-*- delete temporary objects
   fFormulaList->Clear();
   delete [] var;
   delete [] xvars;

   return principal;
}

////////////////////////////////////////////////////////////////////////////////
/// Process this tree executing the TSelector code in the specified filename.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.
///
/// The code in filename is loaded (interpreted or compiled, see below),
/// filename must contain a valid class implementation derived from TSelector,
/// where TSelector has the following member functions:
///
/// -  Begin():        called every time a loop on the tree starts,
///                    a convenient place to create your histograms.
/// -  SlaveBegin():   called after Begin(), when on PROOF called only on the
///                    slave servers.
/// -  Process():      called for each event, in this function you decide what
///                    to read and fill your histograms.
/// -  SlaveTerminate: called at the end of the loop on the tree, when on PROOF
///                    called only on the slave servers.
/// -  Terminate():    called at the end of the loop on the tree,
///                    a convenient place to draw/fit your histograms.
///
/// If filename is of the form file.C, the file will be interpreted.
/// If filename is of the form file.C++, the file file.C will be compiled
/// and dynamically loaded.
///
/// If filename is of the form file.C+, the file file.C will be compiled
/// and dynamically loaded. At next call, if file.C is older than file.o
/// and file.so, the file.C is not compiled, only file.so is loaded.
///
///  ### NOTE 1
///  It may be more interesting to invoke directly the other Process function
///  accepting a TSelector* as argument.eg
/// ~~~{.cpp}
///     MySelector *selector = (MySelector*)TSelector::GetSelector(filename);
///     selector->CallSomeFunction(..);
///     mytree.Process(selector,..);
/// ~~~
///  ### NOTE 2
///  One should not call this function twice with the same selector file
///  in the same script. If this is required, proceed as indicated in NOTE1,
///  by getting a pointer to the corresponding TSelector,eg
///#### workaround 1
/// ~~~{.cpp}
///void stubs1() {
///   TSelector *selector = TSelector::GetSelector("h1test.C");
///   TFile *f1 = new TFile("stubs_nood_le1.root");
///   TTree *h1 = (TTree*)f1->Get("h1");
///   h1->Process(selector);
///   TFile *f2 = new TFile("stubs_nood_le1_coarse.root");
///   TTree *h2 = (TTree*)f2->Get("h1");
///   h2->Process(selector);
///}
/// ~~~
///  or use ACLIC to compile the selector
///#### workaround 2
/// ~~~{.cpp}
///void stubs2() {
///   TFile *f1 = new TFile("stubs_nood_le1.root");
///   TTree *h1 = (TTree*)f1->Get("h1");
///   h1->Process("h1test.C+");
///   TFile *f2 = new TFile("stubs_nood_le1_coarse.root");
///   TTree *h2 = (TTree*)f2->Get("h1");
///   h2->Process("h1test.C+");
///}
/// ~~~

Long64_t TTreePlayer::Process(const char *filename,Option_t *option, Long64_t nentries, Long64_t firstentry)
{
   DeleteSelectorFromFile(); //delete previous selector if any

   // This might reloads the script and delete your option
   // string! so let copy it first:
   TString opt(option);
   TString file(filename);
   TSelector *selector = TSelector::GetSelector(file);
   if (!selector) return -1;

   fSelectorFromFile = selector;
   fSelectorClass    = selector->IsA();

   Long64_t nsel = Process(selector,opt,nentries,firstentry);
   return nsel;
}

////////////////////////////////////////////////////////////////////////////////
/// Process this tree executing the code in the specified selector.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.
///
///   The TSelector class has the following member functions:
///
/// -  Begin():        called every time a loop on the tree starts,
///                    a convenient place to create your histograms.
/// -  SlaveBegin():   called after Begin(), when on PROOF called only on the
///                    slave servers.
/// -  Process():      called for each event, in this function you decide what
///                    to read and fill your histograms.
/// -  SlaveTerminate: called at the end of the loop on the tree, when on PROOF
///                    called only on the slave servers.
/// -  Terminate():    called at the end of the loop on the tree,
///                    a convenient place to draw/fit your histograms.
///
///  If the Tree (Chain) has an associated EventList, the loop is on the nentries
///  of the EventList, starting at firstentry, otherwise the loop is on the
///  specified Tree entries.

Long64_t TTreePlayer::Process(TSelector *selector,Option_t *option, Long64_t nentries, Long64_t firstentry)
{
   nentries = GetEntriesToProcess(firstentry, nentries);

   TDirectory::TContext ctxt;

   fTree->SetNotify(selector);

   selector->SetOption(option);

   selector->Begin(fTree);       //<===call user initialization function
   selector->SlaveBegin(fTree);  //<===call user initialization function
   if (selector->Version() >= 2)
      selector->Init(fTree);
   selector->Notify();

   if (gMonitoringWriter)
      gMonitoringWriter->SendProcessingStatus("STARTED",kTRUE);

   Bool_t process = (selector->GetAbort() != TSelector::kAbortProcess &&
                    (selector->Version() != 0 || selector->GetStatus() != -1)) ? kTRUE : kFALSE;
   if (process) {

      Long64_t readbytesatstart = 0;
      readbytesatstart = TFile::GetFileBytesRead();

      //set the file cache
      TTreeCache *tpf = 0;
      TFile *curfile = fTree->GetCurrentFile();
      if (curfile && fTree->GetCacheSize() > 0) {
         tpf = (TTreeCache*)curfile->GetCacheRead(fTree);
         if (tpf)
            tpf->SetEntryRange(firstentry,firstentry+nentries);
         else {
            fTree->SetCacheSize(fTree->GetCacheSize());
            tpf = (TTreeCache*)curfile->GetCacheRead(fTree);
            if (tpf) tpf->SetEntryRange(firstentry,firstentry+nentries);
         }
      }

      //Create a timer to get control in the entry loop(s)
      TProcessEventTimer *timer = 0;
      Int_t interval = fTree->GetTimerInterval();
      if (!gROOT->IsBatch() && interval)
         timer = new TProcessEventTimer(interval);

      //loop on entries (elist or all entries)
      Long64_t entry, entryNumber, localEntry;

      Bool_t useCutFill = selector->Version() == 0;

      // force the first monitoring info
      if (gMonitoringWriter)
         gMonitoringWriter->SendProcessingProgress(0,0,kTRUE);

      //trying to set the first tree, because in the Draw function
      //the tree corresponding to firstentry has already been loaded,
      //so it is not set in the entry list
      fSelectorUpdate = selector;
      UpdateFormulaLeaves();

      for (entry=firstentry;entry<firstentry+nentries;entry++) {
         entryNumber = fTree->GetEntryNumber(entry);
         if (entryNumber < 0) break;
         if (timer && timer->ProcessEvents()) break;
         if (gROOT->IsInterrupted()) break;
         localEntry = fTree->LoadTree(entryNumber);
         if (localEntry < 0) break;
         if(useCutFill) {
            if (selector->ProcessCut(localEntry))
               selector->ProcessFill(localEntry); //<==call user analysis function
         } else {
            selector->Process(localEntry);        //<==call user analysis function
         }
         if (gMonitoringWriter)
            gMonitoringWriter->SendProcessingProgress((entry-firstentry),TFile::GetFileBytesRead()-readbytesatstart,kTRUE);
         if (selector->GetAbort() == TSelector::kAbortProcess) break;
         if (selector->GetAbort() == TSelector::kAbortFile) {
            // Skip to the next file.
            entry += fTree->GetTree()->GetEntries() - localEntry;
            // Reset the abort status.
            selector->ResetAbort();
         }
      }
      delete timer;
      //we must reset the cache
      {
         TFile *curfile2 = fTree->GetCurrentFile();
         if (curfile2 && fTree->GetCacheSize() > 0) {
            tpf = (TTreeCache*)curfile2->GetCacheRead(fTree);
            if (tpf) tpf->SetEntryRange(0,0);
         }
      }
   }

   process = (selector->GetAbort() != TSelector::kAbortProcess &&
             (selector->Version() != 0 || selector->GetStatus() != -1)) ? kTRUE : kFALSE;
   Long64_t res = (process) ? 0 : -1;
   if (process) {
      selector->SlaveTerminate();   //<==call user termination function
      selector->Terminate();        //<==call user termination function
      res = selector->GetStatus();
   }
   fTree->SetNotify(0); // Detach the selector from the tree.
   fSelectorUpdate = 0;
   if (gMonitoringWriter)
      gMonitoringWriter->SendProcessingStatus("DONE");

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// cleanup pointers in the player pointing to obj

void TTreePlayer::RecursiveRemove(TObject *obj)
{
   if (fHistogram == obj) fHistogram = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Loop on Tree and print entries passing selection. If varexp is 0 (or "")
/// then print only first 8 columns. If varexp = "*" print all columns.
/// Otherwise a columns selection can be made using "var1:var2:var3".
/// The function returns the number of entries passing the selection.
///
/// By default 50 rows are shown and you are asked for `<CR>`
/// to see the next 50 rows.
///
/// You can change the default number of rows to be shown before `<CR>`
/// via  mytree->SetScanField(maxrows) where maxrows is 50 by default.
/// if maxrows is set to 0 all rows of the Tree are shown.
///
/// This option is interesting when dumping the contents of a Tree to
/// an ascii file, eg from the command line.
/// ### with ROOT 5
/// ~~~{.cpp}
///   root [0] tree->SetScanField(0);
///   root [1] tree->Scan("*"); >tree.log
/// ~~~
/// ### with ROOT 6
/// ~~~{.cpp}
///   root [0] tree->SetScanField(0);
///   root [1] .> tree.log
///   tree->Scan("*");
///   .>
/// ~~~
///  will create a file tree.log
///
/// Arrays (within an entry) are printed in their linear forms.
/// If several arrays with multiple dimensions are printed together,
/// they will NOT be synchronized.  For example print
///   arr1[4][2] and arr2[2][3] will results in a printing similar to:
/// ~~~{.cpp}
/// ***********************************************
/// *    Row   * Instance *      arr1 *      arr2 *
/// ***********************************************
/// *        x *        0 * arr1[0][0]* arr2[0][0]*
/// *        x *        1 * arr1[0][1]* arr2[0][1]*
/// *        x *        2 * arr1[1][0]* arr2[0][2]*
/// *        x *        3 * arr1[1][1]* arr2[1][0]*
/// *        x *        4 * arr1[2][0]* arr2[1][1]*
/// *        x *        5 * arr1[2][1]* arr2[1][2]*
/// *        x *        6 * arr1[3][0]*           *
/// *        x *        7 * arr1[3][1]*           *
/// ~~~
/// However, if there is a selection criterion which is an array, then
/// all the formulas will be synchronized with the selection criterion
/// (see TTreePlayer::DrawSelect for more information).
///
/// The options string can contains the following parameters:
///
/// -  lenmax=dd
///       Where 'dd' is the maximum number of elements per array that should
///       be printed.  If 'dd' is 0, all elements are printed (this is the
///       default)
/// -  colsize=ss
///       Where 'ss' will be used as the default size for all the column
///       If this options is not specified, the default column size is 9
/// -  precision=pp
///       Where 'pp' will be used as the default 'precision' for the
///       printing format.
/// -  col=xxx
///       Where 'xxx' is colon (:) delimited list of printing format for
///       each column. The format string should follow the printf format
///       specification.  The value given will be prefixed by % and, if no
///       conversion specifier is given, will be suffixed by the letter g.
///       before being passed to fprintf.  If no format is specified for a
///       column, the default is used  (aka ${colsize}.${precision}g )
///
/// For example:
/// ~~~{.cpp}
///   tree->Scan("a:b:c","","colsize=30 precision=3 col=::20.10:#x:5ld");
/// ~~~
/// Will print 3 columns, the first 2 columns will be 30 characters long,
/// the third columns will be 20 characters long.  The printing format used
/// for the columns (assuming they are numbers) will be respectively:
/// ~~~ {.cpp}
/// %30.3g %30.3g %20.10g %#x %5ld
/// ~~~

Long64_t TTreePlayer::Scan(const char *varexp, const char *selection,
                           Option_t * option,
                           Long64_t nentries, Long64_t firstentry)
{

   TString opt = option;
   opt.ToLower();
   UInt_t ui;
   UInt_t lenmax = 0;
   UInt_t colDefaultSize = 9;
   UInt_t colPrecision = 9;
   std::vector<TString> colFormats;
   std::vector<Int_t> colSizes;

   if (opt.Contains("lenmax=")) {
      int start = opt.Index("lenmax=");
      int numpos = start + strlen("lenmax=");
      int numlen = 0;
      int len = opt.Length();
      while( (numpos+numlen<len) && isdigit(opt[numpos+numlen]) ) numlen++;
      TString num = opt(numpos,numlen);
      opt.Remove(start,strlen("lenmax")+numlen);

      lenmax = atoi(num.Data());
   }
   if (opt.Contains("colsize=")) {
      int start = opt.Index("colsize=");
      int numpos = start + strlen("colsize=");
      int numlen = 0;
      int len = opt.Length();
      while( (numpos+numlen<len) && isdigit(opt[numpos+numlen]) ) numlen++;
      TString num = opt(numpos,numlen);
      opt.Remove(start,strlen("size")+numlen);

      colDefaultSize = atoi(num.Data());
      colPrecision = colDefaultSize;
      if (colPrecision>18) colPrecision = 18;
   }
   if (opt.Contains("precision=")) {
      int start = opt.Index("precision=");
      int numpos = start + strlen("precision=");
      int numlen = 0;
      int len = opt.Length();
      while( (numpos+numlen<len) && isdigit(opt[numpos+numlen]) ) numlen++;
      TString num = opt(numpos,numlen);
      opt.Remove(start,strlen("precision")+numlen);

      colPrecision = atoi(num.Data());
   }
   TString defFormat = Form("%d.%d",colDefaultSize,colPrecision);
   if (opt.Contains("col=")) {
      int start = opt.Index("col=");
      int numpos = start + strlen("col=");
      int numlen = 0;
      int len = opt.Length();
      while( (numpos+numlen<len) &&
             (isdigit(opt[numpos+numlen])
              || opt[numpos+numlen] == 'c'
              || opt[numpos+numlen] == 'd'
              || opt[numpos+numlen] == 'i'
              || opt[numpos+numlen] == 'o'
              || opt[numpos+numlen] == 'x'
              || opt[numpos+numlen] == 'X'
              || opt[numpos+numlen] == 'u'
              || opt[numpos+numlen] == 'f'
              || opt[numpos+numlen] == 'e'
              || opt[numpos+numlen] == 'E'
              || opt[numpos+numlen] == 'g'
              || opt[numpos+numlen] == 'G'
              || opt[numpos+numlen] == 'l'
              || opt[numpos+numlen] == 'L'
              || opt[numpos+numlen] == 'h'
              || opt[numpos+numlen] == 's'
              || opt[numpos+numlen] == '#'
              || opt[numpos+numlen]=='.'
              || opt[numpos+numlen]==':')) numlen++;
      TString flist = opt(numpos,numlen);
      opt.Remove(start,strlen("col")+numlen);

      int i = 0;
      while(i<flist.Length() && flist[i]==':') {
         colFormats.push_back(defFormat);
         colSizes.push_back(colDefaultSize);
         ++i;
      }
      for(; i<flist.Length(); ++i) {
         int next = flist.Index(":",i);
         if (next==i) {
            colFormats.push_back(defFormat);
         } else if (next==kNPOS) {
            colFormats.push_back(flist(i,flist.Length()-i));
            i = flist.Length();
         } else {
            colFormats.push_back(flist(i,next-i));
            i = next;
         }
         UInt_t siz = atoi(colFormats[colFormats.size()-1].Data());
         colSizes.push_back( siz ? siz : colDefaultSize );
      }
   }

   TTreeFormula **var;
   std::vector<TString> cnames;
   TString onerow;
   Long64_t entry,entryNumber;
   Int_t i,nch;
   UInt_t ncols = 8;   // by default first 8 columns are printed only
   std::ofstream out;
   const char *fname = nullptr;
   TString fownname;
   if (fScanRedirect) {
      fTree->SetScanField(0);  // no page break if Scan is redirected
      fname = fScanFileName;
      if (!fname) fname = "";
      Int_t lenfile = strlen(fname);
      if (!lenfile) {
         fownname = fTree->GetName();
         fownname.Append("-scan.dat");
         fname = fownname.Data();
      }
      out.open(fname, std::ios::out);
      if (!out.good ()) {
         Error("Scan","Can not open file for redirection");
         return 0;
      }
   }
   TObjArray *leaves = fTree->GetListOfLeaves();
   if (leaves==0) return 0;
   UInt_t nleaves = leaves->GetEntriesFast();
   if (nleaves < ncols) ncols = nleaves;
   nch = varexp ? strlen(varexp) : 0;

   nentries = GetEntriesToProcess(firstentry, nentries);

//*-*- Compile selection expression if there is one
   TTreeFormula        *select  = 0;
   if (selection && strlen(selection)) {
      select = new TTreeFormula("Selection",selection,fTree);
      if (!select) return -1;
      if (!select->GetNdim()) { delete select; return -1; }
      fFormulaList->Add(select);
   }
//*-*- if varexp is empty, take first 8 columns by default
   int allvar = 0;
   if (varexp && !strcmp(varexp, "*")) { ncols = nleaves; allvar = 1; }
   if (nch == 0 || allvar) {
      UInt_t ncs = ncols;
      ncols = 0;
      for (ui=0;ui<ncs;++ui) {
         TLeaf *lf = (TLeaf*)leaves->At(ui);
         if (lf->GetBranch()->GetListOfBranches()->GetEntries() > 0) continue;
         cnames.push_back( lf->GetBranch()->GetMother()->GetName() );
         if (cnames[ncols] == lf->GetName() ) {
            // Already complete, let move on.
         } else if (cnames[ncols][cnames[ncols].Length()-1]=='.') {
            cnames[ncols] = lf->GetBranch()->GetName(); // name of branch already include mother's name
         } else {
            if (lf->GetBranch()->GetMother()->IsA()->InheritsFrom(TBranchElement::Class())) {
               TBranchElement *mother = (TBranchElement*)lf->GetBranch()->GetMother();
               if (mother->GetType() == 3 || mother->GetType() == 4) {
                  // The name of the mother branch is embedded in the sub-branch names.
                  cnames[ncols] = lf->GetBranch()->GetName();
                  ++ncols;
                  continue;
               }
            }
            if (!strchr(lf->GetBranch()->GetName() ,'[') ) {
               cnames[ncols].Append('.');
               cnames[ncols].Append( lf->GetBranch()->GetName() );
            }
         }
         if (lf->GetBranch()->IsA() == TBranch::Class() ||
             strcmp( lf->GetBranch()->GetName(), lf->GetName() ) != 0 ) {
            cnames[ncols].Append('.');
            cnames[ncols].Append( lf->GetName() );
         }
         ++ncols;
      }
//*-*- otherwise select only the specified columns
   } else {

      ncols = fSelector->SplitNames(varexp, cnames);

   }
   var = new TTreeFormula* [ncols];

   for(ui=colFormats.size();ui<ncols;++ui) {
      colFormats.push_back(defFormat);
      colSizes.push_back(colDefaultSize);
   }

//*-*- Create the TreeFormula objects corresponding to each column
   for (ui=0;ui<ncols;ui++) {
      var[ui] = new TTreeFormula("Var1",cnames[ui].Data(),fTree);
      fFormulaList->Add(var[ui]);
   }

//*-*- Create a TreeFormulaManager to coordinate the formulas
   TTreeFormulaManager *manager=0;
   Bool_t hasArray = kFALSE;
   Bool_t forceDim = kFALSE;
   if (fFormulaList->LastIndex()>=0) {
      if (select) {
         if (select->GetManager()->GetMultiplicity() > 0 ) {
            manager = new TTreeFormulaManager;
            for(i=0;i<=fFormulaList->LastIndex();i++) {
               manager->Add((TTreeFormula*)fFormulaList->At(i));
            }
            manager->Sync();
         }
      }
      for(i=0;i<=fFormulaList->LastIndex();i++) {
         TTreeFormula *form = ((TTreeFormula*)fFormulaList->At(i));
         switch( form->GetManager()->GetMultiplicity() ) {
            case  1:
            case  2:
               hasArray = kTRUE;
               forceDim = kTRUE;
               break;
            case -1:
               forceDim = kTRUE;
               break;
            case  0:
               break;
         }

      }
   }

//*-*- Print header
   onerow = "***********";
   if (hasArray) onerow += "***********";

   for (ui=0;ui<ncols;ui++) {
      TString starFormat = Form("*%%%d.%ds",colSizes[ui]+2,colSizes[ui]+2);
      onerow += Form(starFormat.Data(),var[ui]->PrintValue(-2));
   }
   if (fScanRedirect)
      out<<onerow.Data()<<"*"<<std::endl;
   else
      printf("%s*\n",onerow.Data());
   onerow = "*    Row   ";
   if (hasArray) onerow += "* Instance ";
   for (ui=0;ui<ncols;ui++) {
      TString numbFormat = Form("* %%%d.%ds ",colSizes[ui],colSizes[ui]);
      onerow += Form(numbFormat.Data(),var[ui]->PrintValue(-1));
   }
   if (fScanRedirect)
      out<<onerow.Data()<<"*"<<std::endl;
   else
      printf("%s*\n",onerow.Data());
   onerow = "***********";
   if (hasArray) onerow += "***********";
   for (ui=0;ui<ncols;ui++) {
      TString starFormat = Form("*%%%d.%ds",colSizes[ui]+2,colSizes[ui]+2);
      onerow += Form(starFormat.Data(),var[ui]->PrintValue(-2));
   }
   if (fScanRedirect)
      out<<onerow.Data()<<"*"<<std::endl;
   else
      printf("%s*\n",onerow.Data());
//*-*- loop on all selected entries
   fSelectedRows = 0;
   Int_t tnumber = -1;
   Bool_t exitloop = kFALSE;
   for (entry=firstentry;
        entry<(firstentry+nentries) && !exitloop;
        entry++) {
      entryNumber = fTree->GetEntryNumber(entry);
      if (entryNumber < 0) break;
      Long64_t localEntry = fTree->LoadTree(entryNumber);
      if (localEntry < 0) break;
      if (tnumber != fTree->GetTreeNumber()) {
         tnumber = fTree->GetTreeNumber();
         if (manager) manager->UpdateFormulaLeaves();
         else {
            for(i=0;i<=fFormulaList->LastIndex();i++) {
               ((TTreeFormula*)fFormulaList->At(i))->UpdateFormulaLeaves();
            }
         }
      }

      int ndata = 1;
      if (forceDim) {

         if (manager) {

            ndata = manager->GetNdata(kTRUE);

         } else {

            // let's print the max number of column
            for (ui=0;ui<ncols;ui++) {
               if (ndata < var[ui]->GetNdata() ) {
                  ndata = var[ui]->GetNdata();
               }
            }
            if (select && select->GetNdata()==0) ndata = 0;
         }

      }

      if (lenmax && ndata>(int)lenmax) ndata = lenmax;
      Bool_t loaded = kFALSE;
      for(int inst=0;inst<ndata;inst++) {
         if (select) {
            if (select->EvalInstance(inst) == 0) {
               continue;
            }
         }
         if (inst==0) loaded = kTRUE;
         else if (!loaded) {
            // EvalInstance(0) always needs to be called so that
            // the proper branches are loaded.
            for (ui=0;ui<ncols;ui++) {
               var[ui]->EvalInstance(0);
            }
            loaded = kTRUE;
         }
         onerow = Form("* %8lld ",entryNumber);
         if (hasArray) {
            onerow += Form("* %8d ",inst);
         }
         for (ui=0;ui<ncols;++ui) {
            TString numbFormat = Form("* %%%d.%ds ",colSizes[ui],colSizes[ui]);
            if (var[ui]->GetNdim()) onerow += Form(numbFormat.Data(),var[ui]->PrintValue(0,inst,colFormats[ui].Data()));
            else {
               TString emptyForm = Form("* %%%dc ",colSizes[ui]);
               onerow += Form(emptyForm.Data(),' ');
            }
         }
         fSelectedRows++;
         if (fScanRedirect)
            out<<onerow.Data()<<"*"<<std::endl;
         else
            printf("%s*\n",onerow.Data());
         if (fTree->GetScanField() > 0 && fSelectedRows > 0) {
            if (fSelectedRows%fTree->GetScanField() == 0) {
               fprintf(stderr,"Type <CR> to continue or q to quit ==> ");
               int answer, readch;
               readch = getchar();
               answer = readch;
               while (readch != '\n' && readch != EOF) readch = getchar();
               if (answer == 'q' || answer == 'Q') {
                  exitloop = kTRUE;
                  break;
               }
            }
         }
      }
   }
   onerow = "***********";
   if (hasArray) onerow += "***********";
   for (ui=0;ui<ncols;ui++) {
      TString starFormat = Form("*%%%d.%ds",colSizes[ui]+2,colSizes[ui]+2);
      onerow += Form(starFormat.Data(),var[ui]->PrintValue(-2));
   }
   if (fScanRedirect)
      out<<onerow.Data()<<"*"<<std::endl;
   else
      printf("%s*\n",onerow.Data());
   if (select) Printf("==> %lld selected %s", fSelectedRows,
                      fSelectedRows == 1 ? "entry" : "entries");
   if (fScanRedirect) printf("File <%s> created\n", fname);

//*-*- delete temporary objects
   fFormulaList->Clear();
   // The TTreeFormulaManager is deleted by the last TTreeFormula.
   delete [] var;
   return fSelectedRows;
}

////////////////////////////////////////////////////////////////////////////////
/// Loop on Tree and return TSQLResult object containing entries passing
/// selection. If varexp is 0 (or "") then print only first 8 columns.
/// If varexp = "*" print all columns. Otherwise a columns selection can
/// be made using "var1:var2:var3". In case of error 0 is returned otherwise
/// a TSQLResult object which must be deleted by the user.

TSQLResult *TTreePlayer::Query(const char *varexp, const char *selection,
                               Option_t *, Long64_t nentries, Long64_t firstentry)
{
   TTreeFormula **var;
   std::vector<TString> cnames;
   TString onerow;
   Long64_t entry,entryNumber;
   Int_t i,nch;
   Int_t ncols = 8;   // by default first 8 columns are printed only
   TObjArray *leaves = fTree->GetListOfLeaves();
   Int_t nleaves = leaves->GetEntriesFast();
   if (nleaves < ncols) ncols = nleaves;
   nch = varexp ? strlen(varexp) : 0;

   nentries = GetEntriesToProcess(firstentry, nentries);

   // compile selection expression if there is one
   TTreeFormula *select = 0;
   if (strlen(selection)) {
      select = new TTreeFormula("Selection",selection,fTree);
      if (!select) return 0;
      if (!select->GetNdim()) { delete select; return 0; }
      fFormulaList->Add(select);
   }

   // if varexp is empty, take first 8 columns by default
   int allvar = 0;
   if (varexp && !strcmp(varexp, "*")) { ncols = nleaves; allvar = 1; }
   if (nch == 0 || allvar) {
      for (i=0;i<ncols;i++) {
         cnames.push_back( ((TLeaf*)leaves->At(i))->GetName() );
      }
   } else {
      // otherwise select only the specified columns
      ncols = fSelector->SplitNames(varexp,cnames);
   }
   var = new TTreeFormula* [ncols];

   // create the TreeFormula objects corresponding to each column
   for (i=0;i<ncols;i++) {
      var[i] = new TTreeFormula("Var1",cnames[i].Data(),fTree);
      fFormulaList->Add(var[i]);
   }

   // fill header info into result object
   TTreeResult *res = new TTreeResult(ncols);
   for (i = 0; i < ncols; i++) {
      res->AddField(i, var[i]->PrintValue(-1));
   }

   //*-*- Create a TreeFormulaManager to coordinate the formulas
   TTreeFormulaManager *manager=0;
   if (fFormulaList->LastIndex()>=0) {
      manager = new TTreeFormulaManager;
      for(i=0;i<=fFormulaList->LastIndex();i++) {
         manager->Add((TTreeFormula*)fFormulaList->At(i));
      }
      manager->Sync();
   }

   // loop on all selected entries
   const char *aresult;
   Int_t len;
   char *arow = new char[ncols*50];
   fSelectedRows = 0;
   Int_t tnumber = -1;
   Int_t *fields = new Int_t[ncols];
   for (entry=firstentry;entry<firstentry+nentries;entry++) {
      entryNumber = fTree->GetEntryNumber(entry);
      if (entryNumber < 0) break;
      Long64_t localEntry = fTree->LoadTree(entryNumber);
      if (localEntry < 0) break;
      if (tnumber != fTree->GetTreeNumber()) {
         tnumber = fTree->GetTreeNumber();
         for (i=0;i<ncols;i++) var[i]->UpdateFormulaLeaves();
      }

      Int_t ndata = 1;
      if (manager && manager->GetMultiplicity()) {
         ndata = manager->GetNdata();
      }

      if (select) {
         select->GetNdata();
         if (select->EvalInstance(0) == 0) continue;
      }

      Bool_t loaded = kFALSE;
      for(int inst=0;inst<ndata;inst++) {
         if (select) {
            if (select->EvalInstance(inst) == 0) {
               continue;
            }
         }

         if (inst==0) loaded = kTRUE;
         else if (!loaded) {
            // EvalInstance(0) always needs to be called so that
            // the proper branches are loaded.
            for (i=0;i<ncols;i++) {
               var[i]->EvalInstance(0);
            }
            loaded = kTRUE;
         }
         for (i=0;i<ncols;i++) {
            aresult = var[i]->PrintValue(0,inst);
            len = strlen(aresult)+1;
            if (i == 0) {
               memcpy(arow,aresult,len);
               fields[i] = len;
            } else {
               memcpy(arow+fields[i-1],aresult,len);
               fields[i] = fields[i-1] + len;
            }
         }
         res->AddRow(new TTreeRow(ncols,fields,arow));
         fSelectedRows++;
      }
   }

   // delete temporary objects
   fFormulaList->Clear();
   // The TTreeFormulaManager is deleted by the last TTreeFormula.
   delete [] fields;
   delete [] arow;
   delete [] var;

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of entries to estimate variable limits.

void TTreePlayer::SetEstimate(Long64_t n)
{
   fSelector->SetEstimate(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Start the TTreeViewer on this TTree.
///
/// -  ww is the width of the canvas in pixels
/// -  wh is the height of the canvas in pixels

void TTreePlayer::StartViewer(Int_t ww, Int_t wh)
{
   if (!gApplication)
      TApplication::CreateApplication();
   // make sure that the Gpad and GUI libs are loaded
   TApplication::NeedGraphicsLibs();
   if (gApplication)
      gApplication->InitializeGraphics();
   if (gROOT->IsBatch()) {
      Warning("StartViewer", "The tree viewer cannot run in batch mode");
      return;
   }

   if (ww || wh) { }   // use unused variables
   TPluginHandler *h;
   if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualTreeViewer"))) {
      if (h->LoadPlugin() == -1)
         return;
      h->ExecPlugin(1,fTree);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Unbinned fit of one or more variable(s) from a Tree.
///
/// funcname is a TF1 function.
///
/// See TTree::Draw for explanations of the other parameters.
///
/// Fit the variable varexp using the function funcname using the
/// selection cuts given by selection.
///
/// The list of fit options is given in parameter option.
///
/// -  option = "Q" Quiet mode (minimum printing)
/// -  option = "V" Verbose mode (default is between Q and V)
/// -  option = "E" Perform better Errors estimation using Minos technique
/// -  option = "M" More. Improve fit results
/// -  option = "D" Draw the projected histogram with the fitted function
///             normalized to the number of selected rows
///             and multiplied by the bin width
///
/// You can specify boundary limits for some or all parameters via
/// ~~~{.cpp}
///        func->SetParLimits(p_number, parmin, parmax);
/// ~~~
/// if parmin>=parmax, the parameter is fixed
///
/// Note that you are not forced to fix the limits for all parameters.
/// For example, if you fit a function with 6 parameters, you can do:
/// ~~~{.cpp}
///     func->SetParameters(0,3.1,1.e-6,0.1,-8,100);
///     func->SetParLimits(4,-10,-4);
///     func->SetParLimits(5, 1,1);
/// ~~~
/// With this setup, parameters 0->3 can vary freely
/// -  Parameter 4 has boundaries [-10,-4] with initial value -8
/// -  Parameter 5 is fixed to 100.
///
///   For the fit to be meaningful, the function must be self-normalized.
///
/// i.e. It must have the same integral regardless of the parameter
/// settings.  Otherwise the fit will effectively just maximize the
/// area.
///
/// It is mandatory to have a normalization variable
/// which is fixed for the fit.  e.g.
/// ~~~{.cpp}
///     TF1* f1 = new TF1("f1", "gaus(0)/sqrt(2*3.14159)/[2]", 0, 5);
///     f1->SetParameters(1, 3.1, 0.01);
///     f1->SetParLimits(0, 1, 1); // fix the normalization parameter to 1
///     data->UnbinnedFit("f1", "jpsimass", "jpsipt>3.0");
/// ~~~
///
/// 1, 2 and 3 Dimensional fits are supported.
/// See also TTree::Fit
///
/// ### Return status
///
/// The function return the status of the fit in the following form
/// ~~~{.cpp}
///     fitResult = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult
/// ~~~
/// -  The fitResult is 0 is the fit is OK.
/// -  The fitResult is negative in case of an error not connected with the fit.
/// -  The number of entries used in the fit can be obtained via
/// ~~~{.cpp}
///     mytree.GetSelectedRows();
/// ~~~
/// -  If the number of selected entries is null the function returns -1
///
/// new implementation using new Fitter classes

Int_t TTreePlayer::UnbinnedFit(const char *funcname ,const char *varexp, const char *selection,Option_t *option ,Long64_t nentries, Long64_t firstentry)
{
   // function is given by name, find it in gROOT
   TF1* fitfunc = (TF1*)gROOT->GetFunction(funcname);
   if (!fitfunc) { Error("UnbinnedFit", "Unknown function: %s",funcname); return 0; }

   Int_t npar = fitfunc->GetNpar();
   if (npar <=0) { Error("UnbinnedFit", "Illegal number of parameters = %d",npar); return 0; }

   // Spin through the data to select out the events of interest
   // Make sure that the arrays V1,etc are created large enough to accommodate
   // all entries
   Long64_t oldEstimate = fTree->GetEstimate();
   Long64_t nent = fTree->GetEntriesFriend();
   fTree->SetEstimate(TMath::Min(nent,nentries));

   // build FitOptions
   TString opt = option;
   opt.ToUpper();
   Foption_t fitOption;
   if (opt.Contains("Q")) fitOption.Quiet   = 1;
   if (opt.Contains("V")){fitOption.Verbose = 1; fitOption.Quiet   = 0;}
   if (opt.Contains("E")) fitOption.Errors  = 1;
   if (opt.Contains("M")) fitOption.More    = 1;
   if (!opt.Contains("D")) fitOption.Nograph    = 1;  // what about 0
   // could add range and automatic normalization of functions and gradient

   TString drawOpt = "goff";
   if (!fitOption.Nograph) drawOpt = "";
   Long64_t nsel = DrawSelect(varexp, selection,drawOpt, nentries, firstentry);

   if (!fitOption.Nograph  && GetSelectedRows() <= 0 && GetDimension() > 4) {
      Info("UnbinnedFit","Ignore option D with more than 4 variables");
      nsel = DrawSelect(varexp, selection,"goff", nentries, firstentry);
   }

   //if no selected entries return
   Long64_t nrows = GetSelectedRows();

   if (nrows <= 0) {
      Error("UnbinnedFit", "Cannot fit: no entries selected");
      return -1;
   }

   // Check that function has same dimension as number of variables
   Int_t ndim = GetDimension();
   // do not check with TF1::GetNdim() since it returns 1 for TF1 classes created with
   // a C function with larger dimension


   // use pointer stored in the tree (not copy the data in)
   std::vector<double *> vlist(ndim);
   for (int i = 0; i < ndim; ++i)
      vlist[i] = fSelector->GetVal(i);

   // fill the fit data object
   // the object will be then managed by the fitted classes - however it will be invalid when the
   // data pointers (given by fSelector->GetVal() ) wil be invalidated
   ROOT::Fit::UnBinData * fitdata = new ROOT::Fit::UnBinData(nrows, ndim, vlist.begin());



   ROOT::Math::MinimizerOptions minOption;
   TFitResultPtr ret = ROOT::Fit::UnBinFit(fitdata,fitfunc, fitOption, minOption);

   //reset estimate
   fTree->SetEstimate(oldEstimate);

   //if option "D" is specified, draw the projected histogram
   //with the fitted function normalized to the number of selected rows
   //and multiplied by the bin width
   if (!fitOption.Nograph && fHistogram) {
      if (fHistogram->GetDimension() < 2) {
         TH1 *hf = (TH1*)fHistogram->Clone("unbinnedFit");
         hf->SetLineWidth(3);
         hf->Reset();
         Int_t nbins = fHistogram->GetXaxis()->GetNbins();
         Double_t norm = ((Double_t)nsel)*fHistogram->GetXaxis()->GetBinWidth(1);
         for (Int_t bin=1;bin<=nbins;bin++) {
            Double_t func = norm*fitfunc->Eval(hf->GetBinCenter(bin));
            hf->SetBinContent(bin,func);
         }
         fHistogram->GetListOfFunctions()->Add(hf,"lsame");
      }
      fHistogram->Draw();
   }


   return int(ret);

}

////////////////////////////////////////////////////////////////////////////////
/// this function is called by TChain::LoadTree when a new Tree is loaded.
/// Because Trees in a TChain may have a different list of leaves, one
/// must update the leaves numbers in the TTreeFormula used by the TreePlayer.

void TTreePlayer::UpdateFormulaLeaves()
{
   if (fSelector)  fSelector->Notify();
   if (fSelectorUpdate){
      //If the selector is writing into a TEntryList, the entry list's
      //sublists need to be changed according to the loaded tree
      if (fSelector==fSelectorUpdate) {
         //FIXME: should be more consistent with selector from file
         TObject *obj = fSelector->GetObject();
         if (obj){
            if (fSelector->GetObject()->InheritsFrom(TEntryList::Class())){
               ((TEntryList*)fSelector->GetObject())->SetTree(fTree->GetTree());
            }
         }
      }
      if (fSelectorFromFile==fSelectorUpdate) {
         TIter next(fSelectorFromFile->GetOutputList());
         TEntryList *elist=0;
         while ((elist=(TEntryList*)next())){
            if (elist->InheritsFrom(TEntryList::Class())){
               elist->SetTree(fTree->GetTree());
            }
         }
      }
   }

   if (fFormulaList->GetSize()) {
      TObjLink *lnk = fFormulaList->FirstLink();
      while (lnk) {
         lnk->GetObject()->Notify();
         lnk = lnk->Next();
      }
   }
}
