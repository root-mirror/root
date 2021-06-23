/****************************************************************************
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
\file RooTaskSpec.cxx
\class RooTaskSpec
\ingroup Roofitcore

RooTaskSpec is a precursor to an upgrade of the the multi-processor front-end for 
parallel calculation of RooAbsReal objects. The RooTaskSpec should return a table
containing the information that will be passed to the RooRealMPFE process. The
first calculation that is returned is the case being calculated since Binned and
Unbinned likelihood calcualtions need to be optimised for multiprocessing in very
different ways.
Several cases are possible:
Either the return value of createNLL() is a RooAbsTestStatisic or a RooAddition of
RooAbsTestStatistics and/or a RooConstraintSum. 
Case 1: return value is a RooAbsTestStatistic -> this class returns the tools to
analyse the  components of this likelihood.
Case 2: Multiple return values. The first of which is always a RooAbsTestStatistic
followed by either more RooAbsTestStatistics or a RooConstraintSum which also
signifies the last term. Here we only consider the first term since the others
represent odd and rare test cases. 
ToDo: Simple use demonstration.
**/


#include "RooTaskSpec.h"
// #include "Riostream.h"
// #include "RooFit.h"
// #include <sstream>
#include "RooAbsTestStatistic.h"
#include "RooAddition.h"
#include "RooAbsData.h"

using namespace std;
using namespace RooFit;

//ClassImp(RooTaskSpec);
RooTaskSpec::RooTaskSpec(){};

RooTaskSpec::RooTaskSpec(RooAbsTestStatistic* rats_nll){
//  cout << " NLL is a RooAbsTestStatistic (Case 1)" << endl ;
//  rats_nll->Print(); // WARNING: don't print MPFE values before they're fully initialized! Or make them dirty again afterwards.
  _fit_case = 1;
  _initialise(rats_nll);
}

RooTaskSpec::RooTaskSpec(RooAbsReal* nll){
//  cout <<"starting case2"<<endl;
  RooAddition* ra = dynamic_cast<RooAddition*>(nll) ;
  if (ra) {
//    cout <<"yes ra, printing RooAddition"<<endl;
    ra->Print();
//    cout <<"printed"<<endl;

    RooAbsTestStatistic* rats = dynamic_cast<RooAbsTestStatistic*>(ra->list().at(0)) ;
    if (!rats) {
//      cout << "ERROR: NLL is a RooAddition, but first element of addition is not a RooAbsTestStatistic!" << endl ;
      _fit_case = 0;
//      cout << "It is a "<<ra->list().at(0)<<endl;
    } else {
      _fit_case = 1;
//      cout << "NLL is a RooAddition (Case 2), first element is a RooAbsTestStatistic" << endl ;
      _initialise(rats);
    }
  }
  else {
//    cout<<"not ra, this would not cast to RooAddition"<<endl;
    nll->Print();
  }
}


void RooTaskSpec::_initialise (RooAbsTestStatistic* rats){
  // Check if nll is a AbsTestStatistic
  if (rats->numSimultaneous()==0){
    //    _set.add(_fill_task(0, rats));
    tasks.push_back(_fill_task(0, rats));
  }  else {
    for (Int_t i=0 ; i < rats->numSimultaneous() ; i++) {
//      cout << "SimComponent #" << i << " = " ;
      rats->simComponents()[i]->Print() ;
      RooAbsTestStatistic* comp = (RooAbsTestStatistic*) rats->simComponents()[i] ;
      tasks.push_back(_fill_task(i, comp));
    }
  }
}

RooTaskSpec::Task RooTaskSpec::_fill_task(Int_t n, RooAbsTestStatistic* rats){
  Bool_t b = rats->function().getAttribute("BinnedLikelihood") ;
  Task t;
  t.id = n;
  t.binned = b;
  if (b) {
    t.name = rats->function().GetName();
    t.entries = rats->data().numEntries();
    //   cout << "Binned Likelihood has probability model named " << rats->function().GetName()
    //	 << " and a binned dataset with " << rats->data().numEntries()
    //   << " bins with a total event count of " << rats->data().sumEntries() << endl ;
  } else {
    t.name = rats->function().GetName();
    t.entries = rats->data().numEntries();

    //    cout << "Unbinned likelihood has probability density model named " << rats->function().GetName()
    //    << " and an dataset with " << rats->data().numEntries() 
    //    << " events with a weight of " << rats->data().sumEntries() << endl ;
  } return t;
} 

