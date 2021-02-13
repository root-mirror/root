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
\file RooStudyPackage.cxx
\class RooStudyPackage
\ingroup Roofitcore

RooStudyPackage is a utility class to manage studies that consist of
repeated applications of generate-and-fit operations on a workspace

**/



#include "RooFit.h"
#include "Riostream.h"

#include "RooStudyPackage.h"
#include "RooWorkspace.h"
#include "RooAbsStudy.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "TFile.h"
#include "TRandom2.h"
#include "RooRandom.h"
#include "TMath.h"
#include "TEnv.h"

using namespace std ;

ClassImp(RooStudyPackage);
  ;



////////////////////////////////////////////////////////////////////////////////

RooStudyPackage::RooStudyPackage() : _ws(0)
{  
}



////////////////////////////////////////////////////////////////////////////////

RooStudyPackage::RooStudyPackage(RooWorkspace& w) : _ws(new RooWorkspace(w))
{  
}



////////////////////////////////////////////////////////////////////////////////

RooStudyPackage::RooStudyPackage(const RooStudyPackage& other) : TNamed(other), _ws(new RooWorkspace(*other._ws))
{      
  list<RooAbsStudy*>::const_iterator iter = other._studies.begin() ;
  for (;iter!=other._studies.end() ; ++iter) {
    _studies.push_back((*iter)->clone()) ;
  }
}



////////////////////////////////////////////////////////////////////////////////

void RooStudyPackage::addStudy(RooAbsStudy& study) 
{
  _studies.push_back(&study) ;
}



////////////////////////////////////////////////////////////////////////////////

void RooStudyPackage::driver(Int_t nExperiments)
{
  initialize() ;
  run(nExperiments) ;
  finalize() ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Make iterator over copy of studies attached to workspace

void RooStudyPackage::initialize() 
{
  for (list<RooAbsStudy*>::iterator iter=_studies.begin() ; iter!=_studies.end() ; ++iter) {
    (*iter)->attach(*_ws) ;
    (*iter)->initialize() ;
  }

}


////////////////////////////////////////////////////////////////////////////////

void RooStudyPackage::run(Int_t nExperiments) 
{
  // Run the requested number of experiments
  Int_t prescale = nExperiments>100 ? Int_t(nExperiments/100) : 1 ;
  for (Int_t i=0 ; i<nExperiments ; i++) {
    if (i%prescale==0) {
      coutP(Generation) << "RooStudyPackage::run(" << GetName() << ") processing experiment " << i << "/" << nExperiments << endl ;
    }    
    runOne() ;
  }
}



////////////////////////////////////////////////////////////////////////////////

void RooStudyPackage::runOne() 
{
  for (list<RooAbsStudy*>::iterator iter=_studies.begin() ; iter!=_studies.end() ; ++iter) {
    (*iter)->execute() ;
  }    
}




////////////////////////////////////////////////////////////////////////////////
/// Finalize all studies

void RooStudyPackage::finalize() 
{   
  for (list<RooAbsStudy*>::iterator iter=_studies.begin() ; iter!=_studies.end() ; ++iter) {
    (*iter)->finalize() ;
  }
}




////////////////////////////////////////////////////////////////////////////////

void RooStudyPackage::exportData(TList* olist, Int_t seqno)
{
  for (list<RooAbsStudy*>::iterator iter=_studies.begin() ; iter!=_studies.end() ; ++iter) {

    (*iter)->finalize() ;

    RooDataSet* summaryData = (*iter)->summaryData() ;
    if (summaryData) {
      summaryData->SetName(Form("%s_%d",summaryData->GetName(),seqno)) ;
      cout << "registering summary dataset: " ; summaryData->Print() ;
      olist->Add(summaryData) ;
    }

    RooLinkedList* detailedData = (*iter)->detailedData() ;
    if (detailedData && detailedData->GetSize()>0) {

      detailedData->SetName(Form("%s_%d",detailedData->GetName(),seqno)) ;
      cout << "registering detailed dataset " << detailedData->IsA()->GetName() << "::" 
	   << detailedData->GetName() << " with " << detailedData->GetSize() << " elements" << endl ;
      TIterator* diter = detailedData->MakeIterator() ;
      TNamed* dobj ;
      while((dobj=(TNamed*)diter->Next())) {
	dobj->SetName(Form("%s_%d",dobj->GetName(),seqno)) ;
      }
      delete diter ;
      olist->Add(detailedData) ;
      (*iter)->releaseDetailData() ;
    }
  }        
}



////////////////////////////////////////////////////////////////////////////////
/// Choose random seed for this process
/// in case pass a definite seed to have it deterministic
/// use also worker number

Int_t RooStudyPackage::initRandom()
{
  TRandom2 random(0); 
  //gRandom->SetSeed(0) ;
  Int_t seed = random.Integer(TMath::Limits<Int_t>::Max()) ;

  // get worker number 
  TString  worknumber = gEnv->GetValue("ProofServ.Ordinal","undef");
  int iworker = -1; 
  if (worknumber != "undef") 
     iworker = int( worknumber.Atof()*10 + 0.1);

  if (iworker >= 0)  {
     for (int i = 0; i <= iworker; ++i ) 
        seed = random.Integer( TMath::Limits<Int_t>::Max() ); 
  }

  RooRandom::randomGenerator()->SetSeed(seed) ;
  gRandom->SetSeed(seed) ;

  return seed ;
}



////////////////////////////////////////////////////////////////////////////////
/// Read in study package

void RooStudyPackage::processFile(const char* studyName, Int_t nexp) 
{
  string name_fin = Form("study_data_%s.root",studyName) ;
  TFile fin(name_fin.c_str()) ;
  RooStudyPackage* pkg = dynamic_cast<RooStudyPackage*>(fin.Get("studypack")) ;
  if (!pkg) {
    cout << "RooStudyPackage::processFile() ERROR input file " << name_fin << " does not contain a RooStudyPackage named 'studypack'" << endl ;
    return ;
  }
      
  // Initialize random seed
  Int_t seqno = pkg->initRandom() ;
  cout << "RooStudyPackage::processFile() Initial random seed for this run is " << seqno << endl ;

  // Run study
  pkg->driver(nexp) ;

  // Save result
  TList res ;
  pkg->exportData(&res,seqno) ;
  TFile fout(Form("study_result_%s_%d.root",studyName,seqno),"RECREATE") ;
  res.Write() ;
  fout.Close() ;    
}
