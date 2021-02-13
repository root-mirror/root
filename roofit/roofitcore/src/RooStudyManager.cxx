/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooStudyManager.cxx
\class RooStudyManager
\ingroup Roofitcore

RooStudyManager is a utility class to manage studies that consist of
repeated applications of generate-and-fit operations on a workspace

**/



#include "RooFit.h"
#include "Riostream.h"

#include "RooStudyManager.h"
#include "RooWorkspace.h"
#include "RooAbsStudy.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooStudyPackage.h"
#include "TFile.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TKey.h"
#include <string>
#include "TROOT.h"
#include "TSystem.h"

using namespace std ;

ClassImp(RooStudyManager);


////////////////////////////////////////////////////////////////////////////////

RooStudyManager::RooStudyManager(RooWorkspace& w)
{
  _pkg = new RooStudyPackage(w) ;
}



////////////////////////////////////////////////////////////////////////////////

RooStudyManager::RooStudyManager(RooWorkspace& w, RooAbsStudy& study)
{
  _pkg = new RooStudyPackage(w) ;
  _pkg->addStudy(study) ;
}


////////////////////////////////////////////////////////////////////////////////

RooStudyManager::RooStudyManager(const char* studyPackFileName)
{
  string pwd = gDirectory->GetName() ;
  TFile *f = new TFile(studyPackFileName) ;
  _pkg = dynamic_cast<RooStudyPackage*>(f->Get("studypack")) ;
  gDirectory->cd(Form("%s:",pwd.c_str())) ;
}



////////////////////////////////////////////////////////////////////////////////

void RooStudyManager::addStudy(RooAbsStudy& study)
{
  _pkg->addStudy(study) ;
}




////////////////////////////////////////////////////////////////////////////////

void RooStudyManager::run(Int_t nExperiments)
{
  _pkg->driver(nExperiments) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Open PROOF-Lite session

void RooStudyManager::runProof(Int_t nExperiments, const char* proofHost, Bool_t showGui)
{
  coutP(Generation) << "RooStudyManager::runProof(" << GetName() << ") opening PROOF session" << endl ;
  void* p = (void*) gROOT->ProcessLineFast(Form("TProof::Open(\"%s\")",proofHost)) ;

  // Check that PROOF initialization actually succeeeded
  if (p==0) {
    coutE(Generation) << "RooStudyManager::runProof(" << GetName() << ") ERROR initializing proof, aborting" << endl ;
    return ;
  }

  // Suppress GUI if so requested
  if (!showGui) {
    gROOT->ProcessLineFast(Form("((TProof*)0x%lx)->SetProgressDialog(0) ;",(ULong_t)p)) ;
  }

  // Propagate workspace to proof nodes
  coutP(Generation) << "RooStudyManager::runProof(" << GetName() << ") sending work package to PROOF servers" << endl ;
  gROOT->ProcessLineFast(Form("((TProof*)0x%lx)->AddInput((TObject*)0x%lx) ;",(ULong_t)p,(ULong_t)_pkg) ) ;

  // Run selector in parallel
  coutP(Generation) << "RooStudyManager::runProof(" << GetName() << ") starting PROOF processing of " << nExperiments << " experiments" << endl ;

  gROOT->ProcessLineFast(Form("((TProof*)0x%lx)->Process(\"RooProofDriverSelector\",%d) ;",(ULong_t)p,nExperiments)) ;

  // Aggregate results data
  coutP(Generation) << "RooStudyManager::runProof(" << GetName() << ") aggregating results data" << endl ;
  TList* olist = (TList*) gROOT->ProcessLineFast(Form("((TProof*)0x%lx)->GetOutputList()",(ULong_t)p)) ;
  aggregateData(olist) ;

  // cleaning up
  coutP(Generation) << "RooStudyManager::runProof(" << GetName() << ") cleaning up input list" << endl ;
  gROOT->ProcessLineFast(Form("((TProof*)0x%lx)->GetInputList()->Remove((TObject*)0x%lx) ;",(ULong_t)p,(ULong_t)_pkg) ) ;

}


////////////////////////////////////////////////////////////////////////////////
/// "Option_t *option" takes the parameters forwarded to gProof->Close(option).
///
/// This function is intended for scripts that run in loops
/// where it is essential to properly close all connections and delete
/// the TProof instance (frees ports).

void RooStudyManager::closeProof(Option_t *option)
{
  if (gROOT->GetListOfProofs()->LastIndex() != -1  &&  gROOT->ProcessLineFast("gProof;"))
  {
    gROOT->ProcessLineFast(Form("gProof->Close(\"%s\") ;",option)) ;
    gROOT->ProcessLineFast("gProof->CloseProgressDialog() ;") ;

    // CloseProgressDialog does not do anything when run without GUI. This detects
    // whether the proof instance is still there and deletes it if that is the case.
    if (gROOT->GetListOfProofs()->LastIndex() != -1  &&  gROOT->ProcessLineFast("gProof;")) {
      gROOT->ProcessLineFast("delete gProof ;") ;
    }
  } else {
    ooccoutI((TObject*)NULL,Generation) << "RooStudyManager: No global Proof objects. No connections closed." << endl ;
  }
}



////////////////////////////////////////////////////////////////////////////////

void RooStudyManager::prepareBatchInput(const char* studyName, Int_t nExpPerJob, Bool_t unifiedInput=kFALSE)
{
  TFile f(Form("study_data_%s.root",studyName),"RECREATE") ;
  _pkg->Write("studypack") ;
  f.Close() ;

  if (unifiedInput) {

    // Write header of driver script
    ofstream bdr(Form("study_driver_%s.sh",studyName)) ;
    bdr << "#!/bin/sh" << endl
        << Form("if [ ! -f study_data_%s.root ] ; then",studyName) << endl
        << "uudecode <<EOR" << endl ;
    bdr.close() ;

    // Write uuencoded ROOT file (base64) in driver script
    gSystem->Exec(Form("cat study_data_%s.root | uuencode -m study_data_%s.root >> study_driver_%s.sh",studyName,studyName,studyName)) ;

    // Write remainder of deriver script
    ofstream bdr2 (Form("study_driver_%s.sh",studyName),ios::app) ;
    bdr2 << "EOR" << endl
	 << "fi" << endl
	 << "root -l -b <<EOR" << endl
	 << Form("RooStudyPackage::processFile(\"%s\",%d) ;",studyName,nExpPerJob) << endl
	 << ".q" << endl
	 << "EOR" << endl ;
    // Remove binary input file
    gSystem->Unlink(Form("study_data_%s.root",studyName)) ;

    coutI(DataHandling) << "RooStudyManager::prepareBatchInput batch driver file is '" << Form("study_driver_%s.sh",studyName) << "," << endl
			<< "     input data files is embedded in driver script" << endl ;

  } else {

    ofstream bdr(Form("study_driver_%s.sh",studyName)) ;
    bdr << "#!/bin/sh" << endl
	<< "root -l -b <<EOR" << endl
	<< Form("RooStudyPackage::processFile(\"%s\",%d) ;",studyName,nExpPerJob) << endl
	<< ".q" << endl
	<< "EOR" << endl ;

    coutI(DataHandling) << "RooStudyManager::prepareBatchInput batch driver file is '" << Form("study_driver_%s.sh",studyName) << "," << endl
			<< "     input data file is " << Form("study_data_%s.root",studyName) << endl ;

  }
}




////////////////////////////////////////////////////////////////////////////////

void RooStudyManager::processBatchOutput(const char* filePat)
{
  list<string> flist ;
  expandWildCardSpec(filePat,flist) ;

  TList olist ;

  for (list<string>::iterator iter = flist.begin() ; iter!=flist.end() ; ++iter) {
    coutP(DataHandling) << "RooStudyManager::processBatchOutput() now reading file " << *iter << endl ;
    TFile f(iter->c_str()) ;

    TList* list = f.GetListOfKeys() ;
    TIterator* kiter = list->MakeIterator();

    TObject* obj ;
    TKey* key ;
    while((key=(TKey*)kiter->Next())) {
      obj = f.Get(key->GetName()) ;
      TObject* clone = obj->Clone(obj->GetName()) ;
      olist.Add(clone) ;
    }
    delete kiter ;
  }
  aggregateData(&olist) ;
  olist.Delete() ;
}


////////////////////////////////////////////////////////////////////////////////

void RooStudyManager::aggregateData(TList* olist)
{
  for (list<RooAbsStudy*>::iterator iter=_pkg->studies().begin() ; iter!=_pkg->studies().end() ; ++iter) {
    (*iter)->aggregateSummaryOutput(olist) ;
  }
}




////////////////////////////////////////////////////////////////////////////////
/// case with one single file

void RooStudyManager::expandWildCardSpec(const char* name, list<string>& result)
{
  if (!TString(name).MaybeWildcard()) {
    result.push_back(name) ;
    return ;
  }

   // wildcarding used in name
   TString basename(name);

   Int_t dotslashpos = -1;
   {
      Int_t next_dot = basename.Index(".root");
      while(next_dot>=0) {
         dotslashpos = next_dot;
         next_dot = basename.Index(".root",dotslashpos+1);
      }
      if (basename[dotslashpos+5]!='/') {
         // We found the 'last' .root in the name and it is not followed by
         // a '/', so the tree name is _not_ specified in the name.
         dotslashpos = -1;
      }
   }
   //Int_t dotslashpos = basename.Index(".root/");
   TString behind_dot_root;
   if (dotslashpos>=0) {
      // Copy the tree name specification
      behind_dot_root = basename(dotslashpos+6,basename.Length()-dotslashpos+6);
      // and remove it from basename
      basename.Remove(dotslashpos+5);
   }

   Int_t slashpos = basename.Last('/');
   TString directory;
   if (slashpos>=0) {
      directory = basename(0,slashpos); // Copy the directory name
      basename.Remove(0,slashpos+1);      // and remove it from basename
   } else {
      directory = gSystem->UnixPathName(gSystem->WorkingDirectory());
   }

   TString expand_directory = directory;
   gSystem->ExpandPathName(expand_directory);
   void *dir = gSystem->OpenDirectory(expand_directory.Data());

   if (dir) {
      //create a TList to store the file names (not yet sorted)
      TList l;
      TRegexp re(basename,kTRUE);
      const char *file;
      while ((file = gSystem->GetDirEntry(dir))) {
         if (!strcmp(file,".") || !strcmp(file,"..")) continue;
         TString s = file;
         if ( (basename!=file) && s.Index(re) == kNPOS) continue;
         l.Add(new TObjString(file));
      }
      gSystem->FreeDirectory(dir);
      //sort the files in alphanumeric order
      l.Sort();
      TIter next(&l);
      TObjString *obj;
      while ((obj = (TObjString*)next())) {
         file = obj->GetName();
         if (behind_dot_root.Length() != 0)
            result.push_back(Form("%s/%s/%s",directory.Data(),file,behind_dot_root.Data())) ;
         else
            result.push_back(Form("%s/%s",directory.Data(),file)) ;
      }
      l.Delete();
   }
}
