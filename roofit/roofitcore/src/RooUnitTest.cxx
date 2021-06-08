/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF, verkerke@nikhef.nl                         *
 *                                                                           *
 * Copyright (c) 2000-2011, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooUnitTest.cxx
\class RooUnitTest
\ingroup Roofitcore

RooUnit test is an abstract base class for unit regression tests for
RooFit and RooStats tests performed in stressRooFit and stressRooStats
Implementations of this class must implement abstract method testCode()
which defines the regression test to be performed. Inside testCode()
the regression test can define objects on which the regression is performed.
These are:
Object          | function
----------------|------------
   RooPlot      | regPlot()
   RooFitResult | regResult()
   Double_t     | regValue()
   RooTable     | regTable()
   TH1/2/3      | regTH()
   RooWorkspace | regWS()
**/

#include "RooFit.h"
#include "RooUnitTest.h"
#include "TDirectory.h"
#include "TClass.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooDouble.h"
#include "RooTrace.h"
#include "RooRandom.h"
#include <cmath>

ClassImp(RooUnitTest);

using namespace std;

TDirectory* RooUnitTest::gMemDir = nullptr;


////////////////////////////////////////////////////////////////////////////////

RooUnitTest::RooUnitTest(const char* name, TFile* refFile, Bool_t writeRef, Int_t verbose, int batchMode) : TNamed(name,name),
  			                         _refFile(refFile), _debug(kFALSE), _write(writeRef), _verb(verbose), _batchMode(batchMode)
{
}



////////////////////////////////////////////////////////////////////////////////

RooUnitTest::~RooUnitTest()
{
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::regPlot(RooPlot* frame, const char* refName)
{
  if (_refFile) {
    string refNameStr(refName) ;
    frame->SetName(refName) ;
    _regPlots.push_back(make_pair(frame,refNameStr)) ;
  } else {
    delete frame ;
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::regResult(RooFitResult* r, const char* refName)
{
  if (_refFile) {
    string refNameStr(refName) ;
    _regResults.push_back(make_pair(r,refNameStr)) ;
  } else {
    delete r ;
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::regValue(Double_t d, const char* refName)
{
  if (_refFile) {
    string refNameStr(refName) ;
    _regValues.push_back(make_pair(d,refNameStr)) ;
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::regTable(RooTable* t, const char* refName)
{
  if (_refFile) {
    string refNameStr(refName) ;
    _regTables.push_back(make_pair(t,refNameStr)) ;
  } else {
    delete t ;
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::regWS(RooWorkspace* ws, const char* refName)
{
  if (_refFile) {
    string refNameStr(refName) ;
    _regWS.push_back(make_pair(ws,refNameStr)) ;
  } else {
    delete ws ;
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::regTH(TH1* th, const char* refName)
{
  if (_refFile) {
    string refNameStr(refName) ;
    _regTH.push_back(make_pair(th,refNameStr)) ;
  } else {
    delete th ;
  }
}


////////////////////////////////////////////////////////////////////////////////

RooWorkspace* RooUnitTest::getWS(const char* refName)
{
  RooWorkspace* ws = dynamic_cast<RooWorkspace*>(_refFile->Get(refName)) ;
  if (!ws) {
    cout << "RooUnitTest ERROR: cannot retrieve RooWorkspace " << refName
	 << " from reference file, skipping " << endl ;
    return 0 ;
  }

  return ws ;
}


////////////////////////////////////////////////////////////////////////////////

Bool_t RooUnitTest::areTHidentical(TH1* htest, TH1* href)
{
  if (htest->GetDimension() != href->GetDimension()) {
    return kFALSE ;
  }

  // Use Kolmogorov distance as metric rather than probability
  // because we expect histograms to be identical rather
  // than drawn from the same parent distribution
  Double_t kmax = htest->KolmogorovTest(href,"M") ;

  if (kmax>htol()) {

    cout << "KS distances = " << kmax << endl ;

    Int_t ntest = htest->GetNbinsX() +2 ;
    Int_t nref  = href->GetNbinsX() +2 ;
    if (htest->GetDimension()>1) {
      ntest *= htest->GetNbinsY() + 2 ;
      nref *= href->GetNbinsY() + 2 ;
    }
    if (htest->GetDimension()>2) {
      ntest *= htest->GetNbinsZ() + 2 ;
      nref *= href->GetNbinsZ() + 2 ;
    }

    if (ntest != nref) {
      return kFALSE ;
    }

    for (Int_t i=0 ; i<ntest ; i++) {
      if (fabs(htest->GetBinContent(i)-href->GetBinContent(i))>htol()) {
	cout << "htest[" << i << "] = " << htest->GetBinContent(i) << " href[" << i << "] = " << href->GetBinContent(i) << endl;
      }
    }

    return kFALSE ;
  }

  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////

Bool_t RooUnitTest::runCompTests()
{
  Bool_t ret = kTRUE ;

  list<pair<RooPlot*, string> >::iterator iter = _regPlots.begin() ;
  while (iter!=_regPlots.end()) {

    if (!_write) {

      // Comparison mode

      // Retrieve benchmark
      RooPlot* bmark = dynamic_cast<RooPlot*>(_refFile->Get(iter->second.c_str())) ;
      if (!bmark) {
	cout << "RooUnitTest ERROR: cannot retrieve RooPlot " << iter->second << " from reference file, skipping " << endl ;
	ret = kFALSE ;
	++iter ;
	continue ;
      }

      if (_verb) {
	cout << "comparing RooPlot " << iter->first << " to benchmark " << iter->second << " = " << bmark << endl ;
	cout << "reference: " ; iter->first->Print() ;
	cout << "benchmark: " ; bmark->Print() ;
      }

      RooPlot* compPlot = _debug ? iter->first->emptyClone(Form("%s_comparison",iter->first->GetName())) : 0 ;
      Bool_t anyFail=kFALSE ;

      Stat_t nItems = iter->first->numItems() ;
      for (Stat_t i=0 ; i<nItems ; i++) {
	// coverity[NULL_RETURNS]
	TObject* obj = iter->first->getObject((Int_t)i) ;

	// Retrieve corresponding object from reference frame
	TObject* objRef = bmark->findObject(obj->GetName()) ;

	if (!objRef) {
	  cout << "RooUnitTest ERROR: cannot retrieve object " << obj->GetName() << " from reference  RooPlot " << iter->second << ", skipping" << endl ;
	  ret = kFALSE ;
	  break ;
	}

	// Histogram comparisons
	if (obj->IsA()==RooHist::Class()) {
	  RooHist* testHist = static_cast<RooHist*>(obj) ;
	  RooHist* refHist = static_cast<RooHist*>(objRef) ;
	  if (!testHist->isIdentical(*refHist,htol())) {
	    cout << "RooUnitTest ERROR: comparison of object " << obj->IsA()->GetName() << "::" << obj->GetName()
		 <<   " fails comparison with counterpart in reference RooPlot " << bmark->GetName() << endl ;

	    if (compPlot) {
	      compPlot->addPlotable((RooHist*)testHist->Clone(),"P") ;
	      compPlot->getAttLine()->SetLineColor(kRed) ;
	      compPlot->getAttMarker()->SetMarkerColor(kRed) ;
	      compPlot->getAttLine()->SetLineWidth(1) ;

	      compPlot->addPlotable((RooHist*)refHist->Clone(),"P") ;
	      compPlot->getAttLine()->SetLineColor(kBlue) ;
	      compPlot->getAttMarker()->SetMarkerColor(kBlue) ;
	      compPlot->getAttLine()->SetLineWidth(1) ;
	    }

	    anyFail=kTRUE ;
	    ret = kFALSE ;
	  }
	} else if (obj->IsA()==RooCurve::Class()) {
	  RooCurve* testCurve = static_cast<RooCurve*>(obj) ;
	  RooCurve* refCurve = static_cast<RooCurve*>(objRef) ;
	  if (!testCurve->isIdentical(*refCurve,ctol())) {
	    cout << "RooUnitTest ERROR: comparison of object " << obj->IsA()->GetName() << "::" << obj->GetName()
		 <<   " fails comparison with counterpart in reference RooPlot " << bmark->GetName() << endl ;

	    if (compPlot) {
	      compPlot->addPlotable((RooCurve*)testCurve->Clone()) ;
	      compPlot->getAttLine()->SetLineColor(kRed) ;
	      compPlot->getAttLine()->SetLineWidth(1) ;
	      compPlot->getAttLine()->SetLineStyle(kSolid) ;

	      compPlot->addPlotable((RooCurve*)refCurve->Clone()) ;
	      compPlot->getAttLine()->SetLineColor(kBlue) ;
	      compPlot->getAttLine()->SetLineWidth(1) ;
	      compPlot->getAttLine()->SetLineStyle(kDashed) ;
	    }

	    anyFail=kTRUE ;
	    ret = kFALSE ;
	  }

	}

      }

      if (anyFail && compPlot) {
	cout << "RooUnitTest INFO: writing comparison plot " << compPlot->GetName() << " of failed test to RooUnitTest_DEBUG.root" << endl ;
	TFile fdbg("RooUnitTest_DEBUG.root","UPDATE") ;
	compPlot->Write() ;
	fdbg.Close() ;
      } else {
	delete compPlot ;
      }

      // Delete RooPlot when comparison is finished to avoid noise in leak checking
      delete iter->first ;

    } else {

      // Writing mode

      cout <<"RooUnitTest: Writing reference RooPlot " << iter->first << " as benchmark " << iter->second << endl ;
      _refFile->cd() ;
      iter->first->Write(iter->second.c_str()) ;
      gMemDir->cd() ;
    }

    ++iter ;
  }


  list<pair<RooFitResult*, string> >::iterator iter2 = _regResults.begin() ;
  while (iter2!=_regResults.end()) {

    if (!_write) {

      // Comparison mode

     // Retrieve benchmark
      RooFitResult* bmark = dynamic_cast<RooFitResult*>(_refFile->Get(iter2->second.c_str())) ;
      if (!bmark) {
	cout << "RooUnitTest ERROR: cannot retrieve RooFitResult " << iter2->second << " from reference file, skipping " << endl ;
	++iter2 ;
	ret = kFALSE ;
	continue ;
      }

      if (_verb) {
	cout << "comparing RooFitResult " << iter2->first << " to benchmark " << iter2->second << " = " << bmark << endl ;
      }

      if (!iter2->first->isIdentical(*bmark,fptol(),fctol())) {
	cout << "RooUnitTest ERROR: comparison of object " << iter2->first->IsA()->GetName() << "::" << iter2->first->GetName()
	     << " from result " << iter2->second
	     <<   " fails comparison with counterpart in reference RooFitResult " << bmark->GetName() << endl ;
	ret = kFALSE ;
      }

      // Delete RooFitResult when comparison is finished to avoid noise in leak checking
      delete iter2->first ;


    } else {

      // Writing mode

      cout <<"RooUnitTest: Writing reference RooFitResult " << iter2->first << " as benchmark " << iter2->second << endl ;
      _refFile->cd() ;
      iter2->first->Write(iter2->second.c_str()) ;
      gMemDir->cd() ;
    }

    ++iter2 ;
  }

  list<pair<Double_t, string> >::iterator iter3 = _regValues.begin() ;
  while (iter3!=_regValues.end()) {

    if (!_write) {

      // Comparison mode

     // Retrieve benchmark
      RooDouble* ref = dynamic_cast<RooDouble*>(_refFile->Get(iter3->second.c_str())) ;
      if (!ref) {
	cout << "RooUnitTest ERROR: cannot retrieve RooDouble " << iter3->second << " from reference file, skipping " << endl ;
	++iter3 ;
	ret = kFALSE ;
	continue ;
      }

      if (_verb) {
	cout << "comparing value " << iter3->first << " to benchmark " << iter3->second << " = " << (Double_t)(*ref) << endl ;
      }

      if (fabs(iter3->first - (Double_t)(*ref))>vtol() ) {
	cout << "RooUnitTest ERROR: comparison of value " << iter3->first <<   " fails comparison with reference " << ref->GetName() << endl ;
	ret = kFALSE ;
      }


    } else {

      // Writing mode

      cout <<"RooUnitTest: Writing reference Double_t " << iter3->first << " as benchmark " << iter3->second << endl ;
      _refFile->cd() ;
      RooDouble* rd = new RooDouble(iter3->first) ;
      rd->Write(iter3->second.c_str()) ;
      gMemDir->cd() ;
    }

    ++iter3 ;
  }


  list<pair<RooTable*, string> >::iterator iter4 = _regTables.begin() ;
  while (iter4!=_regTables.end()) {

    if (!_write) {

      // Comparison mode

     // Retrieve benchmark
      RooTable* bmark = dynamic_cast<RooTable*>(_refFile->Get(iter4->second.c_str())) ;
      if (!bmark) {
	cout << "RooUnitTest ERROR: cannot retrieve RooTable " << iter4->second << " from reference file, skipping " << endl ;
	++iter4 ;
	ret = kFALSE ;
	continue ;
      }

      if (_verb) {
	cout << "comparing RooTable " << iter4->first << " to benchmark " << iter4->second << " = " << bmark << endl ;
      }

      if (!iter4->first->isIdentical(*bmark)) {
        cout << "RooUnitTest ERROR: comparison of object " << iter4->first->IsA()->GetName() << "::" << iter4->first->GetName()
	         <<   " fails comparison with counterpart in reference RooTable " << bmark->GetName() << endl ;
        if (_verb) {
          iter4->first->Print("V");
          bmark->Print("V");
        }
        ret = false;
      }

      // Delete RooTable when comparison is finished to avoid noise in leak checking
      delete iter4->first ;


    } else {

      // Writing mode

      cout <<"RooUnitTest: Writing reference RooTable " << iter4->first << " as benchmark " << iter4->second << endl ;
      _refFile->cd() ;
      iter4->first->Write(iter4->second.c_str()) ;
      gMemDir->cd() ;
    }

    ++iter4 ;
  }


  list<pair<RooWorkspace*, string> >::iterator iter5 = _regWS.begin() ;
  while (iter5!=_regWS.end()) {

    if (_write) {

      // Writing mode

      cout <<"RooUnitTest: Writing reference RooWorkspace " << iter5->first << " as benchmark " << iter5->second << endl ;
      _refFile->cd() ;
      iter5->first->Write(iter5->second.c_str()) ;
      gMemDir->cd() ;
    }

    ++iter5 ;
  }

  /////////////////
  list<pair<TH1*, string> >::iterator iter6 = _regTH.begin() ;
  while (iter6!=_regTH.end()) {

    if (!_write) {

      // Comparison mode

     // Retrieve benchmark
      TH1* bmark = dynamic_cast<TH1*>(_refFile->Get(iter6->second.c_str())) ;
      if (!bmark) {
	cout << "RooUnitTest ERROR: cannot retrieve TH1 " << iter6->second << " from reference file, skipping " << endl ;
	++iter6 ;
	ret = kFALSE ;
	continue ;
      }

      if (_verb) {
	cout << "comparing TH1 " << iter6->first << " to benchmark " << iter6->second << " = " << bmark << endl ;
      }

      if (!areTHidentical(iter6->first,bmark)) {
	// coverity[NULL_RETURNS]
	cout << "RooUnitTest ERROR: comparison of object " << iter6->first->IsA()->GetName() << "::" << iter6->first->GetName()
	     <<   " fails comparison with counterpart in reference TH1 " << bmark->GetName() << endl ;


      if (_debug) {
	cout << "RooUnitTest INFO: writing THx " << iter6->first->GetName() << " and " << bmark->GetName()
	     << " of failed test to RooUnitTest_DEBUG.root" << endl ;
	TFile fdbg("RooUnitTest_DEBUG.root","UPDATE") ;
	iter6->first->SetName(Form("%s_test",iter6->first->GetName())) ;
	iter6->first->Write() ;
	bmark->SetName(Form("%s_ref",bmark->GetName())) ;
	bmark->Write() ;
	fdbg.Close() ;
      }

	ret = kFALSE ;
      }

      // Delete TH1 when comparison is finished to avoid noise in leak checking
      delete iter6->first ;


    } else {

      // Writing mode

      cout <<"RooUnitTest: Writing reference TH1 " << iter6->first << " as benchmark " << iter6->second << endl ;
      _refFile->cd() ;
      iter6->first->Write(iter6->second.c_str()) ;
      gMemDir->cd() ;
    }

    ++iter6 ;
  }


  /////////////////

  return ret ;
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::setSilentMode()
{
  RooMsgService::instance().setSilentMode(kTRUE) ;
  for (Int_t i=0 ; i<RooMsgService::instance().numStreams() ; i++) {
    if (RooMsgService::instance().getStream(i).minLevel<RooFit::ERROR) {
      RooMsgService::instance().setStreamStatus(i,kFALSE) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooUnitTest::clearSilentMode()
{
  RooMsgService::instance().setSilentMode(kFALSE) ;
  for (Int_t i=0 ; i<RooMsgService::instance().numStreams() ; i++) {
    RooMsgService::instance().setStreamStatus(i,kTRUE) ;
  }
}



////////////////////////////////////////////////////////////////////////////////

Bool_t RooUnitTest::runTest()
{
  gMemDir->cd() ;

  if (_verb<2) {
    setSilentMode() ;
  } else {
    cout << "*** Begin of output of Unit Test at normal verbosity *************" << endl ;
  }

  RooMsgService::instance().clearErrorCount() ;

  // Reset random generator seed to make results independent of test ordering
  gRandom->SetSeed(12345) ;
  RooRandom::randomGenerator()->SetSeed(12345) ;

  RooTrace::callgrind_zero() ;
  if (!testCode()) return kFALSE ;
  RooTrace::callgrind_dump() ;

  if (_verb<2) {
    clearSilentMode() ;
  } else {
    cout << "*** End of output of Unit Test at normal verbosity ***************" << endl ;
  }

  if (RooMsgService::instance().errorCount()>0) {
    cout << "RooUnitTest: ERROR messages were logged, failing test" << endl ;
    return kFALSE ;
  }

  return runCompTests() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set gMemDir to memDir

void RooUnitTest::setMemDir(TDirectory* memDir) {
   gMemDir = memDir ;
}
