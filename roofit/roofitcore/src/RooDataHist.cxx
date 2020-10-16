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
\file RooDataHist.cxx
\class RooDataHist
\ingroup Roofitcore

The RooDataHist is a container class to hold N-dimensional binned data. Each bin's central
coordinates in N-dimensional space are represented by a RooArgSet containing RooRealVar, RooCategory
or RooStringVar objects, thus data can be binned in real and/or discrete dimensions.
**/

#include "RooDataHist.h"

#include "RooFit.h"
#include "Riostream.h"
#include "RooMsgService.h"
#include "RooDataHistSliceIter.h"
#include "RooAbsLValue.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "RooHistError.h"
#include "RooCategory.h"
#include "RooCmdConfig.h"
#include "RooLinkedListIter.h"
#include "RooTreeDataStore.h"
#include "RooVectorDataStore.h"
#include "RooTrace.h"
#include "RooFitLegacy/RooTreeData.h"
#include "RooHelpers.h"
#include "RooFormulaVar.h"
#include "RooFormula.h"
#include "RooUniformBinning.h"
#include "RooSpan.h"

#include "TH1.h"
#include "TTree.h"
#include "TBuffer.h"
#include "TMath.h"
#include "Math/Util.h"

using namespace std;

ClassImp(RooDataHist);


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooDataHist::RooDataHist() : _pbinvCacheMgr(0,10)
{
  _arrSize = 0 ;
  _wgt = 0 ;
  _errLo = 0 ;
  _errHi = 0 ;
  _sumw2 = 0 ;
  _binv = 0 ;
  _pbinv = 0 ;
  _curWeight = 0 ;
  _curIndex = -1 ;
  _binValid = 0 ;
  _curSumW2 = 0 ;
  _curVolume = 1 ;
  _curWgtErrHi = 0 ;
  _curWgtErrLo = 0 ;
  _cache_sum_valid = 0 ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of an empty data hist from a RooArgSet defining the dimensions
/// of the data space. The range and number of bins in each dimensions are taken
/// from getMin()getMax(),getBins() of each RooAbsArg representing that
/// dimension.
///
/// For real dimensions, the fit range and number of bins can be set independently
/// of the plot range and number of bins, but it is advisable to keep the
/// ratio of the plot bin width and the fit bin width an integer value.
/// For category dimensions, the fit ranges always comprises all defined states
/// and each state is always has its individual bin
///
/// To effective achive binning of real dimensions with variable bins sizes,
/// construct a RooThresholdCategory of the real dimension to be binned variably.
/// Set the thresholds at the desired bin boundaries, and construct the
/// data hist as function of the threshold category instead of the real variable.

RooDataHist::RooDataHist(const char *name, const char *title, const RooArgSet& vars, const char* binningName) : 
  RooAbsData(name,title,vars), _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) : 
                                         ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  
  initialize(binningName) ;

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;

  appendToDir(this,kTRUE) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a data hist from an existing data collection (binned or unbinned)
/// The RooArgSet 'vars' defines the dimensions of the histogram. 
/// The range and number of bins in each dimensions are taken
/// from getMin(), getMax(), getBins() of each argument passed.
///
/// For real dimensions, the fit range and number of bins can be set independently
/// of the plot range and number of bins, but it is advisable to keep the
/// ratio of the plot bin width and the fit bin width an integer value.
/// For category dimensions, the fit ranges always comprises all defined states
/// and each state is always has its individual bin
///
/// To effective achive binning of real dimensions with variable bins sizes,
/// construct a RooThresholdCategory of the real dimension to be binned variably.
/// Set the thresholds at the desired bin boundaries, and construct the
/// data hist as function of the threshold category instead of the real variable.
///
/// If the constructed data hist has less dimensions that in source data collection,
/// all missing dimensions will be projected.

RooDataHist::RooDataHist(const char *name, const char *title, const RooArgSet& vars, const RooAbsData& data, Double_t wgt) :
  RooAbsData(name,title,vars), _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) : 
                                         ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;

  initialize() ;
  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;

  add(data,(const RooFormulaVar*)0,wgt) ;
  appendToDir(this,kTRUE) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a data hist from a map of TH1,TH2 or TH3 that are collated into a x+1 dimensional
/// RooDataHist where the added dimension is a category that labels the input source as defined
/// in the histMap argument. The state names used in histMap must correspond to predefined states
/// 'indexCat'
///
/// The RooArgList 'vars' defines the dimensions of the histogram. 
/// The ranges and number of bins are taken from the input histogram and must be the same in all histograms

RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, RooCategory& indexCat, 
			 map<string,TH1*> histMap, Double_t wgt) :
  RooAbsData(name,title,RooArgSet(vars,&indexCat)), 
  _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) : 
                                         ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  
  importTH1Set(vars, indexCat, histMap, wgt, kFALSE) ;

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a data hist from a map of RooDataHists that are collated into a x+1 dimensional
/// RooDataHist where the added dimension is a category that labels the input source as defined
/// in the histMap argument. The state names used in histMap must correspond to predefined states
/// 'indexCat'
///
/// The RooArgList 'vars' defines the dimensions of the histogram. 
/// The ranges and number of bins are taken from the input histogram and must be the same in all histograms

RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, RooCategory& indexCat, 
			 map<string,RooDataHist*> dhistMap, Double_t wgt) :
  RooAbsData(name,title,RooArgSet(vars,&indexCat)), 
  _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) : 
                                         ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  
  importDHistSet(vars, indexCat, dhistMap, wgt) ;

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a data hist from an TH1,TH2 or TH3
/// The RooArgSet 'vars' defines the dimensions of the histogram. The ranges
/// and number of bins are taken from the input histogram, and the corresponding
/// values are set accordingly on the arguments in 'vars'

RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, const TH1* hist, Double_t wgt) :
  RooAbsData(name,title,vars), _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) : 
                                         ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;

  // Check consistency in number of dimensions
  if (vars.getSize() != hist->GetDimension()) {
    coutE(InputArguments) << "RooDataHist::ctor(" << GetName() << ") ERROR: dimension of input histogram must match "
			  << "number of dimension variables" << endl ;
    assert(0) ; 
  }

  importTH1(vars,*hist,wgt, kFALSE) ;

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a binned dataset from a RooArgSet defining the dimensions
/// of the data space. The range and number of bins in each dimensions are taken
/// from getMin() getMax(),getBins() of each RooAbsArg representing that
/// dimension.
///
/// <table>
/// <tr><th> Optional Argument <th> Effect
/// <tr><td> Import(TH1&, Bool_t impDens) <td> Import contents of the given TH1/2/3 into this binned dataset. The
///                                 ranges and binning of the binned dataset are automatically adjusted to
///                                 match those of the imported histogram. 
///
///                                 Please note: for TH1& with unequal binning _only_,
///                                 you should decide if you want to import the absolute bin content,
///                                 or the bin content expressed as density. The latter is default and will
///                                 result in the same histogram as the original TH1. For certain types of
///                                 bin contents (containing efficiencies, asymmetries, or ratio is general)
///                                 you should import the absolute value and set impDens to kFALSE
///                                 
///
/// <tr><td> Weight(Double_t)          <td> Apply given weight factor when importing histograms
///
/// <tr><td> Index(RooCategory&)       <td> Prepare import of multiple TH1/1/2/3 into a N+1 dimensional RooDataHist
///                              where the extra discrete dimension labels the source of the imported histogram
///                              If the index category defines states for which no histogram is be imported
///                              the corresponding bins will be left empty.
///                              
/// <tr><td> Import(const char*, TH1&) <td> Import a THx to be associated with the given state name of the index category
///                              specified in Index(). If the given state name is not yet defined in the index
///                              category it will be added on the fly. The import command can be specified
///                              multiple times. 
/// <tr><td> Import(map<string,TH1*>&) <td> As above, but allows specification of many imports in a single operation
///                              

RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3,
			 const RooCmdArg& arg4,const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) :
  RooAbsData(name,title,RooArgSet(vars,(RooAbsArg*)RooCmdConfig::decodeObjOnTheFly("RooDataHist::RooDataHist", "IndexCat",0,0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8))), 
  _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) : 
                                         ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;

  // Define configuration for this method
  RooCmdConfig pc(Form("RooDataHist::ctor(%s)",GetName())) ;
  pc.defineObject("impHist","ImportHisto",0) ;
  pc.defineInt("impDens","ImportHisto",0) ;
  pc.defineObject("indexCat","IndexCat",0) ;
  pc.defineObject("impSliceHist","ImportHistoSlice",0,0,kTRUE) ; // array
  pc.defineString("impSliceState","ImportHistoSlice",0,"",kTRUE) ; // array
  pc.defineObject("impSliceDHist","ImportDataHistSlice",0,0,kTRUE) ; // array
  pc.defineString("impSliceDState","ImportDataHistSlice",0,"",kTRUE) ; // array
  pc.defineDouble("weight","Weight",0,1) ; 
  pc.defineObject("dummy1","ImportDataHistSliceMany",0) ;
  pc.defineObject("dummy2","ImportHistoSliceMany",0) ;
  pc.defineMutex("ImportHisto","ImportHistoSlice","ImportDataHistSlice") ;
  pc.defineDependency("ImportHistoSlice","IndexCat") ;
  pc.defineDependency("ImportDataHistSlice","IndexCat") ;

  RooLinkedList l ;
  l.Add((TObject*)&arg1) ;  l.Add((TObject*)&arg2) ;  
  l.Add((TObject*)&arg3) ;  l.Add((TObject*)&arg4) ;
  l.Add((TObject*)&arg5) ;  l.Add((TObject*)&arg6) ;  
  l.Add((TObject*)&arg7) ;  l.Add((TObject*)&arg8) ;

  // Process & check varargs 
  pc.process(l) ;
  if (!pc.ok(kTRUE)) {
    assert(0) ;
    return ;
  }

  TH1* impHist = static_cast<TH1*>(pc.getObject("impHist")) ;
  Bool_t impDens = pc.getInt("impDens") ;
  Double_t initWgt = pc.getDouble("weight") ;
  const char* impSliceNames = pc.getString("impSliceState","",kTRUE) ;
  const RooLinkedList& impSliceHistos = pc.getObjectList("impSliceHist") ;
  RooCategory* indexCat = static_cast<RooCategory*>(pc.getObject("indexCat")) ;
  const char* impSliceDNames = pc.getString("impSliceDState","",kTRUE) ;
  const RooLinkedList& impSliceDHistos = pc.getObjectList("impSliceDHist") ;


  if (impHist) {
    
    // Initialize importing contents from TH1
    importTH1(vars,*impHist,initWgt, impDens) ;

  } else if (indexCat) {


    if (impSliceHistos.GetSize()>0) {

      // Initialize importing mapped set of TH1s
      map<string,TH1*> hmap ;
      TIterator* hiter = impSliceHistos.MakeIterator() ;
      for (const auto& token : RooHelpers::tokenise(impSliceNames, ",")) {
        auto histo = static_cast<TH1*>(hiter->Next());
        assert(histo);
        hmap[token] = histo;
      }
      importTH1Set(vars,*indexCat,hmap,initWgt,kFALSE) ;
    } else {

      // Initialize importing mapped set of RooDataHists
      map<string,RooDataHist*> dmap ;
      TIterator* hiter = impSliceDHistos.MakeIterator() ;
      for (const auto& token : RooHelpers::tokenise(impSliceDNames, ",")) {
        dmap[token] = (RooDataHist*) hiter->Next() ;
      }
      importDHistSet(vars,*indexCat,dmap,initWgt) ;
    }


  } else {

    // Initialize empty
    initialize() ;
    appendToDir(this,kTRUE) ;    

  }

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  TRACE_CREATE

}




////////////////////////////////////////////////////////////////////////////////
/// Import data from given TH1/2/3 into this RooDataHist

void RooDataHist::importTH1(const RooArgList& vars, const TH1& histo, Double_t wgt, Bool_t doDensityCorrection) 
{
  // Adjust binning of internal observables to match that of input THx
  Int_t offset[3]{0, 0, 0};
  adjustBinning(vars, histo, offset) ;
  
  // Initialize internal data structure
  initialize() ;
  appendToDir(this,kTRUE) ;

  // Define x,y,z as 1st, 2nd and 3rd observable
  RooRealVar* xvar = (RooRealVar*) _vars.find(vars.at(0)->GetName()) ;
  RooRealVar* yvar = (RooRealVar*) (vars.at(1) ? _vars.find(vars.at(1)->GetName()) : 0 ) ;
  RooRealVar* zvar = (RooRealVar*) (vars.at(2) ? _vars.find(vars.at(2)->GetName()) : 0 ) ;

  // Transfer contents
  Int_t xmin(0),ymin(0),zmin(0) ;
  RooArgSet vset(*xvar) ;
  Double_t volume = xvar->getMax()-xvar->getMin() ;
  xmin = offset[0] ;
  if (yvar) {
    vset.add(*yvar) ;
    ymin = offset[1] ;
    volume *= (yvar->getMax()-yvar->getMin()) ;
  }
  if (zvar) {
    vset.add(*zvar) ;
    zmin = offset[2] ;
    volume *= (zvar->getMax()-zvar->getMin()) ;
  }
  //Double_t avgBV = volume / numEntries() ;
//   cout << "average bin volume = " << avgBV << endl ;

  Int_t ix(0),iy(0),iz(0) ;
  for (ix=0 ; ix < xvar->getBins() ; ix++) {
    xvar->setBin(ix) ;
    if (yvar) {
      for (iy=0 ; iy < yvar->getBins() ; iy++) {
        yvar->setBin(iy) ;
        if (zvar) {
          for (iz=0 ; iz < zvar->getBins() ; iz++) {
            zvar->setBin(iz) ;
            Double_t bv = doDensityCorrection ? binVolume(vset) : 1;
            add(vset,bv*histo.GetBinContent(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,bv*TMath::Power(histo.GetBinError(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,2)) ;
          }
        } else {
          Double_t bv = doDensityCorrection ? binVolume(vset) : 1;
          add(vset,bv*histo.GetBinContent(ix+1+xmin,iy+1+ymin)*wgt,bv*TMath::Power(histo.GetBinError(ix+1+xmin,iy+1+ymin)*wgt,2)) ;
        }
      }
    } else {
      Double_t bv = doDensityCorrection ? binVolume(vset) : 1 ;
      add(vset,bv*histo.GetBinContent(ix+1+xmin)*wgt,bv*TMath::Power(histo.GetBinError(ix+1+xmin)*wgt,2)) ;	    
    }
  }  

}





////////////////////////////////////////////////////////////////////////////////
/// Import data from given set of TH1/2/3 into this RooDataHist. The category indexCat labels the sources
/// in the constructed RooDataHist. The stl map provides the mapping between the indexCat state labels
/// and the import source

void RooDataHist::importTH1Set(const RooArgList& vars, RooCategory& indexCat, map<string,TH1*> hmap, Double_t wgt, Bool_t doDensityCorrection) 
{
  RooCategory* icat = (RooCategory*) _vars.find(indexCat.GetName()) ;

  TH1* histo(0) ;  
  Bool_t init(kFALSE) ;
  for (map<string,TH1*>::iterator hiter = hmap.begin() ; hiter!=hmap.end() ; ++hiter) {
    // Store pointer to first histogram from which binning specification will be taken
    if (!histo) {
      histo = hiter->second ;
    }
    // Define state labels in index category (both in provided indexCat and in internal copy in dataset)
    if (!indexCat.lookupType(hiter->first.c_str())) {
      indexCat.defineType(hiter->first) ;
      coutI(InputArguments) << "RooDataHist::importTH1Set(" << GetName() << ") defining state \"" << hiter->first << "\" in index category " << indexCat.GetName() << endl ;
    }
    if (!icat->lookupType(hiter->first.c_str())) {	
      icat->defineType(hiter->first) ;
    }
  }

  // Check consistency in number of dimensions
  if (histo && (vars.getSize() != histo->GetDimension())) {
    coutE(InputArguments) << "RooDataHist::ctor(" << GetName() << ") ERROR: dimension of input histogram must match "
			  << "number of continuous variables" << endl ;
    assert(0) ; 
  }
  
  // Copy bins and ranges from THx to dimension observables
  Int_t offset[3] ;
  adjustBinning(vars,*histo,offset) ;
  
  // Initialize internal data structure
  if (!init) {
    initialize() ;
    appendToDir(this,kTRUE) ;
    init = kTRUE ;
  }

  // Define x,y,z as 1st, 2nd and 3rd observable
  RooRealVar* xvar = (RooRealVar*) _vars.find(vars.at(0)->GetName()) ;
  RooRealVar* yvar = (RooRealVar*) (vars.at(1) ? _vars.find(vars.at(1)->GetName()) : 0 ) ;
  RooRealVar* zvar = (RooRealVar*) (vars.at(2) ? _vars.find(vars.at(2)->GetName()) : 0 ) ;

  // Transfer contents
  Int_t xmin(0),ymin(0),zmin(0) ;
  RooArgSet vset(*xvar) ;
  Double_t volume = xvar->getMax()-xvar->getMin() ;
  xmin = offset[0] ;
  if (yvar) {
    vset.add(*yvar) ;
    ymin = offset[1] ;
    volume *= (yvar->getMax()-yvar->getMin()) ;
  }
  if (zvar) {
    vset.add(*zvar) ;
    zmin = offset[2] ;
    volume *= (zvar->getMax()-zvar->getMin()) ;
  }
  Double_t avgBV = volume / numEntries() ;
  
  Int_t ic(0),ix(0),iy(0),iz(0) ;
  for (ic=0 ; ic < icat->numBins(0) ; ic++) {
    icat->setBin(ic) ;
    histo = hmap[icat->getCurrentLabel()] ;
    for (ix=0 ; ix < xvar->getBins() ; ix++) {
      xvar->setBin(ix) ;
      if (yvar) {
	for (iy=0 ; iy < yvar->getBins() ; iy++) {
	  yvar->setBin(iy) ;
	  if (zvar) {
	    for (iz=0 ; iz < zvar->getBins() ; iz++) {
	      zvar->setBin(iz) ;
	      Double_t bv = doDensityCorrection ? binVolume(vset)/avgBV : 1;
	      add(vset,bv*histo->GetBinContent(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,bv*TMath::Power(histo->GetBinError(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,2)) ;
	    }
	  } else {
	    Double_t bv = doDensityCorrection ? binVolume(vset)/avgBV : 1;
	    add(vset,bv*histo->GetBinContent(ix+1+xmin,iy+1+ymin)*wgt,bv*TMath::Power(histo->GetBinError(ix+1+xmin,iy+1+ymin)*wgt,2)) ;
	  }
	}
      } else {
	Double_t bv = doDensityCorrection ? binVolume(vset)/avgBV : 1;
	add(vset,bv*histo->GetBinContent(ix+1+xmin)*wgt,bv*TMath::Power(histo->GetBinError(ix+1+xmin)*wgt,2)) ;	    
      }
    }  
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Import data from given set of TH1/2/3 into this RooDataHist. The category indexCat labels the sources
/// in the constructed RooDataHist. The stl map provides the mapping between the indexCat state labels
/// and the import source

void RooDataHist::importDHistSet(const RooArgList& /*vars*/, RooCategory& indexCat, std::map<std::string,RooDataHist*> dmap, Double_t initWgt) 
{
  RooCategory* icat = (RooCategory*) _vars.find(indexCat.GetName()) ;

  for (map<string,RooDataHist*>::iterator diter = dmap.begin() ; diter!=dmap.end() ; ++diter) {

    // Define state labels in index category (both in provided indexCat and in internal copy in dataset)
    if (!indexCat.lookupType(diter->first.c_str())) {
      indexCat.defineType(diter->first) ;
      coutI(InputArguments) << "RooDataHist::importDHistSet(" << GetName() << ") defining state \"" << diter->first << "\" in index category " << indexCat.GetName() << endl ;
    }
    if (!icat->lookupType(diter->first.c_str())) {	
      icat->defineType(diter->first) ;
    }
  }

  initialize() ;
  appendToDir(this,kTRUE) ;  


  for (map<string,RooDataHist*>::iterator diter = dmap.begin() ; diter!=dmap.end() ; ++diter) {

    RooDataHist* dhist = diter->second ;

    icat->setLabel(diter->first.c_str()) ;

    // Transfer contents
    for (Int_t i=0 ; i<dhist->numEntries() ; i++) {
      _vars = *dhist->get(i) ;
      add(_vars,dhist->weight()*initWgt, pow(dhist->weightError(SumW2),2) ) ;
    }

  }
}

////////////////////////////////////////////////////////////////////////////////
/// Helper doing the actual work of adjustBinning().

void RooDataHist::_adjustBinning(RooRealVar &theirVar, const TAxis &axis,
    RooRealVar *ourVar, Int_t *offset)
{
  if (!dynamic_cast<RooRealVar*>(ourVar)) {
    coutE(InputArguments) << "RooDataHist::adjustBinning(" << GetName() << ") ERROR: dimension " << ourVar->GetName() << " must be real" << endl ;
    assert(0) ;
  }

  const double xlo = theirVar.getMin();
  const double xhi = theirVar.getMax();

  if (axis.GetXbins()->GetArray()) {
    RooBinning xbins(axis.GetNbins(), axis.GetXbins()->GetArray());

    const double tolerance = 1e-6 * xbins.averageBinWidth();
    
    // Adjust xlo/xhi to nearest boundary
    const double xloAdj = xbins.binLow(xbins.binNumber(xlo + tolerance));
    const double xhiAdj = xbins.binHigh(xbins.binNumber(xhi - tolerance));
    xbins.setRange(xloAdj, xhiAdj);

    theirVar.setBinning(xbins);

    if (true || fabs(xloAdj - xlo) > tolerance || fabs(xhiAdj - xhi) > tolerance) {
      coutI(DataHandling)<< "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << ourVar->GetName() << " expanded to nearest bin boundaries: [" << xlo << "," << xhi << "] --> [" << xloAdj << "," << xhiAdj << "]" << endl;
    }

    ourVar->setBinning(xbins);

    if (offset) {
      *offset = xbins.rawBinNumber(xloAdj + tolerance);
    }
  } else {
    RooBinning xbins(axis.GetXmin(), axis.GetXmax());
    xbins.addUniform(axis.GetNbins(), axis.GetXmin(), axis.GetXmax());

    const double tolerance = 1e-6 * xbins.averageBinWidth();

    // Adjust xlo/xhi to nearest boundary
    const double xloAdj = xbins.binLow(xbins.binNumber(xlo + tolerance));
    const double xhiAdj = xbins.binHigh(xbins.binNumber(xhi - tolerance));
    xbins.setRange(xloAdj, xhiAdj);
    theirVar.setRange(xloAdj, xhiAdj);

    if (fabs(xloAdj - xlo) > tolerance || fabs(xhiAdj - xhi) > tolerance) {
      coutI(DataHandling)<< "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << ourVar->GetName() << " expanded to nearest bin boundaries: [" << xlo << "," << xhi << "] --> [" << xloAdj << "," << xhiAdj << "]" << endl;
    }

    RooUniformBinning xbins2(xloAdj, xhiAdj, xbins.numBins());
    ourVar->setBinning(xbins2);

    if (offset) {
      *offset = xbins.rawBinNumber(xloAdj + tolerance);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Adjust binning specification on first and optionally second and third
/// observable to binning in given reference TH1. Used by constructors
/// that import data from an external TH1.
/// Both the variables in vars and in this RooDataHist are adjusted.
/// @param vars List with variables that are supposed to have their binning adjusted.
/// @param href Reference histogram that dictates the binning
/// @param offset If not nullptr, a possible bin count offset for the axes x,y,z is saved here as Int_t[3]

void RooDataHist::adjustBinning(const RooArgList& vars, const TH1& href, Int_t* offset)
{
  auto xvar = static_cast<RooRealVar*>( _vars.find(*vars.at(0)) );
  _adjustBinning(*static_cast<RooRealVar*>(vars.at(0)), *href.GetXaxis(), xvar, offset ? &offset[0] : nullptr);

  if (vars.at(1)) {
    auto yvar = static_cast<RooRealVar*>(_vars.find(*vars.at(1)));
    if (yvar)
      _adjustBinning(*static_cast<RooRealVar*>(vars.at(1)), *href.GetYaxis(), yvar, offset ? &offset[1] : nullptr);
  }

  if (vars.at(2)) {
    auto zvar = static_cast<RooRealVar*>(_vars.find(*vars.at(2)));
    if (zvar)
      _adjustBinning(*static_cast<RooRealVar*>(vars.at(2)), *href.GetZaxis(), zvar, offset ? &offset[2] : nullptr);
  }

}





////////////////////////////////////////////////////////////////////////////////
/// Initialization procedure: allocate weights array, calculate
/// multipliers needed for N-space to 1-dim array jump table,
/// and fill the internal tree with all bin center coordinates

void RooDataHist::initialize(const char* binningName, Bool_t fillTree)
{

  // Save real dimensions of dataset separately
  for (const auto real : _vars) {
    if (dynamic_cast<RooAbsReal*>(real)) _realVars.add(*real);
  }

  // Fill array of LValue pointers to variables
  for (const auto rvarg : _vars) {
    if (binningName) {
      RooRealVar* rrv = dynamic_cast<RooRealVar*>(rvarg); 
      if (rrv) {
	rrv->setBinning(rrv->getBinning(binningName));
      }
    }
    // coverity[FORWARD_NULL]
    _lvvars.push_back(dynamic_cast<RooAbsLValue*>(rvarg));    
    // coverity[FORWARD_NULL]
    const RooAbsBinning* binning = dynamic_cast<RooAbsLValue*>(rvarg)->getBinningPtr(0);
    _lvbins.push_back(binning ? binning->clone() : 0);
  }

  
  // Allocate coefficients array
  _idxMult.resize(_vars.getSize()) ;

  _arrSize = 1 ;
  Int_t n(0), i ;
  for (const auto var : _vars) {
    auto arg = dynamic_cast<const RooAbsLValue*>(var);
    
    // Calculate sub-index multipliers for master index
    for (i=0 ; i<n ; i++) {
      _idxMult[i] *= arg->numBins() ;
    }
    _idxMult[n++] = 1 ;

    // Calculate dimension of weight array
    _arrSize *= arg->numBins() ;
  }  

  // Allocate and initialize weight array if necessary
  if (!_wgt) {
    _wgt = new Double_t[_arrSize] ;
    _errLo = new Double_t[_arrSize] ;
    _errHi = new Double_t[_arrSize] ;
    _sumw2 = new Double_t[_arrSize] ;
    _binv = new Double_t[_arrSize] ;

    // Refill array pointers in data store when reading
    // from Streamer
    if (fillTree==kFALSE) {
      _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;      
    }
    
    for (i=0 ; i<_arrSize ; i++) {
      _wgt[i] = 0 ;
      _errLo[i] = -1 ;
      _errHi[i] = -1 ;
      _sumw2[i] = 0 ;
    }
  }

  if (!fillTree) return ;

  // Fill TTree with bin center coordinates
  // Calculate plot bins of components from master index

  Int_t ibin ;
  for (ibin=0 ; ibin<_arrSize ; ibin++) {
    Int_t j(0), idx(0), tmp(ibin) ;
    Double_t theBinVolume(1) ;
    for (auto arg2 : _vars) {
      idx  = tmp / _idxMult[j] ;
      tmp -= idx*_idxMult[j++] ;
      auto arglv = dynamic_cast<RooAbsLValue*>(arg2);
      arglv->setBin(idx) ;
      theBinVolume *= arglv->getBinWidth(idx) ;
//       cout << "init: bin width at idx=" << idx << " = " << arglv->getBinWidth(idx) << " binv[" << idx << "] = " << theBinVolume << endl ;
    }
    _binv[ibin] = theBinVolume ;
//     cout << "binv[" << ibin << "] = " << theBinVolume << endl ;
    fill() ;
  }


}


////////////////////////////////////////////////////////////////////////////////

void RooDataHist::checkBinBounds() const
{
  if (!_binbounds.empty()) return;
  for (std::vector<const RooAbsBinning*>::const_iterator it = _lvbins.begin();
      _lvbins.end() != it; ++it) {
    _binbounds.push_back(std::vector<Double_t>());
    if (*it) {
      std::vector<Double_t>& bounds = _binbounds.back();
      bounds.reserve(2 * (*it)->numBins());
      for (Int_t i = 0; i < (*it)->numBins(); ++i) {
	bounds.push_back((*it)->binLow(i));
	bounds.push_back((*it)->binHigh(i));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooDataHist::RooDataHist(const RooDataHist& other, const char* newname) :
  RooAbsData(other,newname), RooDirItem(), _idxMult(other._idxMult), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(other._pbinvCacheMgr,0), _cache_sum_valid(0)
{
  Int_t i ;

  // Allocate and initialize weight array 
  _arrSize = other._arrSize ;
  _wgt = new Double_t[_arrSize] ;
  _errLo = new Double_t[_arrSize] ;
  _errHi = new Double_t[_arrSize] ;
  _binv = new Double_t[_arrSize] ;
  _sumw2 = new Double_t[_arrSize] ;
  for (i=0 ; i<_arrSize ; i++) {
    _wgt[i] = other._wgt[i] ;
    _errLo[i] = other._errLo[i] ;
    _errHi[i] = other._errHi[i] ;
    _sumw2[i] = other._sumw2[i] ;
    _binv[i] = other._binv[i] ;
  }  

  // Save real dimensions of dataset separately
  for (const auto arg : _vars) {
    if (dynamic_cast<RooAbsReal*>(arg) != nullptr) _realVars.add(*arg) ;
  }

  // Fill array of LValue pointers to variables
  for (const auto rvarg : _vars) {
    // coverity[FORWARD_NULL]
    _lvvars.push_back(dynamic_cast<RooAbsLValue*>(rvarg)) ;
    // coverity[FORWARD_NULL]
    const RooAbsBinning* binning = dynamic_cast<RooAbsLValue*>(rvarg)->getBinningPtr(0) ;
    _lvbins.push_back(binning ? binning->clone() : 0) ;    
  }

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;

 appendToDir(this,kTRUE) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a data hist from (part of) an existing data hist. The dimensions
/// of the data set are defined by the 'vars' RooArgSet, which can be identical
/// to 'dset' dimensions, or a subset thereof. Reduced dimensions will be projected
/// in the output data hist. The optional 'cutVar' formula variable can used to 
/// select the subset of bins to be copied.
///
/// For most uses the RooAbsData::reduce() wrapper function, which uses this constructor, 
/// is the most convenient way to create a subset of an existing data  

RooDataHist::RooDataHist(const char* name, const char* title, RooDataHist* h, const RooArgSet& varSubset, 
			 const RooFormulaVar* cutVar, const char* cutRange, Int_t nStart, Int_t nStop, Bool_t copyCache) :
  RooAbsData(name,title,varSubset),
  _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
{
  // Initialize datastore
  _dstore = new RooTreeDataStore(name,title,*h->_dstore,_vars,cutVar,cutRange,nStart,nStop,copyCache) ;
  
  initialize(0,kFALSE) ;

  _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;

  // Copy weight array etc
  Int_t i ;
  for (i=0 ; i<_arrSize ; i++) {
    _wgt[i] = h->_wgt[i] ;
    _errLo[i] = h->_errLo[i] ;
    _errHi[i] = h->_errHi[i] ;
    _sumw2[i] = h->_sumw2[i] ;
    _binv[i] = h->_binv[i] ;
  }  


  appendToDir(this,kTRUE) ;
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////
/// Construct a clone of this dataset that contains only the cached variables

RooAbsData* RooDataHist::cacheClone(const RooAbsArg* newCacheOwner, const RooArgSet* newCacheVars, const char* newName) 
{
  checkInit() ;

  RooDataHist* dhist = new RooDataHist(newName?newName:GetName(),GetTitle(),this,*get(),0,0,0,2000000000,kTRUE) ; 

  RooArgSet* selCacheVars = (RooArgSet*) newCacheVars->selectCommon(dhist->_cachedVars) ;
  dhist->attachCache(newCacheOwner, *selCacheVars) ;
  delete selCacheVars ;

  return dhist ;
}



////////////////////////////////////////////////////////////////////////////////
/// Implementation of RooAbsData virtual method that drives the RooAbsData::reduce() methods

RooAbsData* RooDataHist::reduceEng(const RooArgSet& varSubset, const RooFormulaVar* cutVar, const char* cutRange, 
    std::size_t nStart, std::size_t nStop, Bool_t /*copyCache*/)
{
  checkInit() ;
  RooArgSet* myVarSubset = (RooArgSet*) _vars.selectCommon(varSubset) ;
  RooDataHist *rdh = new RooDataHist(GetName(), GetTitle(), *myVarSubset) ;
  delete myVarSubset ;

  RooFormulaVar* cloneVar = 0;
  RooArgSet* tmp(0) ;
  if (cutVar) {
    // Deep clone cutVar and attach clone to this dataset
    tmp = (RooArgSet*) RooArgSet(*cutVar).snapshot() ;
    if (!tmp) {
      coutE(DataHandling) << "RooDataHist::reduceEng(" << GetName() << ") Couldn't deep-clone cut variable, abort," << endl ;
      return 0 ;
    }
    cloneVar = (RooFormulaVar*) tmp->find(*cutVar) ;
    cloneVar->attachDataSet(*this) ;
  }

  Double_t lo,hi ;
  const std::size_t nevt = nStop < static_cast<std::size_t>(numEntries()) ? nStop : static_cast<std::size_t>(numEntries());
  for (auto i=nStart; i<nevt ; i++) {
    const RooArgSet* row = get(i) ;

    Bool_t doSelect(kTRUE) ;
    if (cutRange) {
      for (auto arg : *get()) {
        if (!arg->inRange(cutRange)) {
          doSelect = kFALSE ;
          break ;
        }
      }
    }
    if (!doSelect) continue ;

    if (!cloneVar || cloneVar->getVal()) {
      weightError(lo,hi,SumW2) ;
      rdh->add(*row,weight(),lo*lo) ;
    }
  }

  if (cloneVar) {
    delete tmp ;
  } 

  return rdh ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooDataHist::~RooDataHist() 
{
  if (_wgt) delete[] _wgt ;
  if (_errLo) delete[] _errLo ;
  if (_errHi) delete[] _errHi ;
  if (_sumw2) delete[] _sumw2 ;
  if (_binv) delete[] _binv ;
  if (_binValid) delete[] _binValid ;
  vector<const RooAbsBinning*>::iterator iter = _lvbins.begin() ;
  while(iter!=_lvbins.end()) {
    delete *iter ;
    ++iter ;
  }

   removeFromDir(this) ;
  TRACE_DESTROY
}




////////////////////////////////////////////////////////////////////////////////
/// Calculate bin number of the given coordinates. If only a subset of the internal
/// coordinates are passed, the missing coordinates are taken at their current value.
/// \param[in] coord Variables that are representing the coordinates.
/// \param[in] fast If the variables in `coord` and the ones of the data hist have the
/// same size and layout, see RooAbsCollection::hasSameLayout(), fast can be set to skip
/// checking that all variables are present in `coord`.
Int_t RooDataHist::getIndex(const RooArgSet& coord, Bool_t fast)
{
  checkInit() ;
  if (fast) {
    _vars.assignFast(coord,kFALSE) ;
  } else {
    _vars.assignValueOnly(coord) ;
  }
  return calcTreeIndex() ;
}




////////////////////////////////////////////////////////////////////////////////
/// Calculate the index for the weights array corresponding to 
/// to the bin enclosing the current coordinates of the internal argset

Int_t RooDataHist::calcTreeIndex() const 
{
  int masterIdx(0);
  for (unsigned int i=0; i < _lvvars.size(); ++i) {
    const RooAbsLValue*  lvvar = _lvvars[i];
    const RooAbsBinning* binning = _lvbins[i];
    masterIdx += _idxMult[i] * lvvar->getBin(binning);
  }

  return masterIdx ;
}



////////////////////////////////////////////////////////////////////////////////
/// Debug stuff, should go...

void RooDataHist::dump2() 
{  
  cout << "_arrSize = " << _arrSize << endl ;
  for (Int_t i=0 ; i<_arrSize ; i++) {
    cout << "wgt[" << i << "] = " << _wgt[i] << "sumw2[" << i << "] = " << _sumw2[i] << " vol[" << i << "] = " << _binv[i] << endl ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Back end function to plotting functionality. Plot RooDataHist on given
/// frame in mode specified by plot options 'o'. The main purpose of
/// this function is to match the specified binning on 'o' to the
/// internal binning of the plot observable in this RooDataHist.
/// \see RooAbsData::plotOn() for plotting options.
RooPlot *RooDataHist::plotOn(RooPlot *frame, PlotOpt o) const 
{
  checkInit() ;
  if (o.bins) return RooAbsData::plotOn(frame,o) ;

  if(0 == frame) {
    coutE(InputArguments) << ClassName() << "::" << GetName() << ":plotOn: frame is null" << endl;
    return 0;
  }
  RooAbsRealLValue *var= (RooAbsRealLValue*) frame->getPlotVar();
  if(0 == var) {
    coutE(InputArguments) << ClassName() << "::" << GetName()
	 << ":plotOn: frame does not specify a plot variable" << endl;
    return 0;
  }

  RooRealVar* dataVar = (RooRealVar*) _vars.find(*var) ;
  if (!dataVar) {
    coutE(InputArguments) << ClassName() << "::" << GetName()
	 << ":plotOn: dataset doesn't contain plot frame variable" << endl;
    return 0;
  }

  o.bins = &dataVar->getBinning() ;
  o.correctForBinWidth = kFALSE ;
  return RooAbsData::plotOn(frame,o) ;
}




////////////////////////////////////////////////////////////////////////////////

Double_t RooDataHist::weightSquared() const {
  return _curSumW2 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the weight at given coordinates with optional
/// interpolation. If intOrder is zero, the weight
/// for the bin enclosing the coordinates
/// contained in 'bin' is returned. For higher values,
/// the result is interpolated in the real dimensions 
/// of the dataset

Double_t RooDataHist::weight(const RooArgSet& bin, Int_t intOrder, Bool_t correctForBinSize, Bool_t cdfBoundaries, Bool_t oneSafe) 
{
  //cout << "RooDataHist::weight(" << bin << "," << intOrder << "," << correctForBinSize << "," << cdfBoundaries << "," << oneSafe << ")" << endl ;

  checkInit() ;

  // Handle illegal intOrder values
  if (intOrder<0) {
    coutE(InputArguments) << "RooDataHist::weight(" << GetName() << ") ERROR: interpolation order must be positive" << endl ;
    return 0 ;
  }

  // Handle no-interpolation case
  if (intOrder==0) {    
    _vars.assignValueOnly(bin,oneSafe) ;
    Int_t idx = calcTreeIndex() ;
    //cout << "intOrder 0, idx = " << idx << endl ;
    if (correctForBinSize) {
      //calculatePartialBinVolume(*get()) ;
      //cout << "binw[" << idx << "] = " << _wgt[idx] <<  " / " << _binv[idx] << endl ;
      return get_wgt(idx) / _binv[idx];
    } else {
      //cout << "binw[" << idx << "] = " << _wgt[idx] << endl ;
      return get_wgt(idx);
    }
  }

  // Handle all interpolation cases
  _vars.assignValueOnly(bin) ;

  Double_t wInt(0) ;
  if (_realVars.getSize()==1) {

    // 1-dimensional interpolation
    const auto real = static_cast<RooRealVar*>(_realVars[static_cast<std::size_t>(0)]);
    const RooAbsBinning* binning = real->getBinningPtr(0) ;
    wInt = interpolateDim(*real,binning,((RooAbsReal*)bin.find(*real))->getVal(), intOrder, correctForBinSize, cdfBoundaries) ;
    
  } else if (_realVars.getSize()==2) {

    // 2-dimensional interpolation
    const auto realX = static_cast<RooRealVar*>(_realVars[static_cast<std::size_t>(0)]);
    const auto realY = static_cast<RooRealVar*>(_realVars[static_cast<std::size_t>(1)]);
    Double_t xval = ((RooAbsReal*)bin.find(*realX))->getVal() ;
    Double_t yval = ((RooAbsReal*)bin.find(*realY))->getVal() ;
    
    Int_t ybinC = realY->getBin() ;
    Int_t ybinLo = ybinC-intOrder/2 - ((yval<realY->getBinning().binCenter(ybinC))?1:0) ;
    Int_t ybinM = realY->numBins() ;
    
    Int_t i ;
    Double_t yarr[10] ;
    Double_t xarr[10] ;
    const RooAbsBinning* binning = realX->getBinningPtr(0) ;
    for (i=ybinLo ; i<=intOrder+ybinLo ; i++) {
      Int_t ibin ;
      if (i>=0 && i<ybinM) {
	// In range
	ibin = i ;
	realY->setBin(ibin) ;
	xarr[i-ybinLo] = realY->getVal() ;
      } else if (i>=ybinM) {
	// Overflow: mirror
	ibin = 2*ybinM-i-1 ;
	realY->setBin(ibin) ;
	xarr[i-ybinLo] = 2*realY->getMax()-realY->getVal() ;
      } else {
	// Underflow: mirror
	ibin = -i -1;
	realY->setBin(ibin) ;
	xarr[i-ybinLo] = 2*realY->getMin()-realY->getVal() ;
      }
      yarr[i-ybinLo] = interpolateDim(*realX,binning,xval,intOrder,correctForBinSize,kFALSE) ;	
    }

    if (gDebug>7) {
      cout << "RooDataHist interpolating data is" << endl ;
      cout << "xarr = " ;
      for (int q=0; q<=intOrder ; q++) cout << xarr[q] << " " ;
      cout << " yarr = " ;
      for (int q=0; q<=intOrder ; q++) cout << yarr[q] << " " ;
      cout << endl ;
    }
    wInt = RooMath::interpolate(xarr,yarr,intOrder+1,yval) ;
    
  } else {

    // Higher dimensional scenarios not yet implemented
    coutE(InputArguments) << "RooDataHist::weight(" << GetName() << ") interpolation in " 
	 << _realVars.getSize() << " dimensions not yet implemented" << endl ;
    return weight(bin,0) ;

  }

  // Cut off negative values
//   if (wInt<=0) {
//     wInt=0 ; 
//   }

  //cout << "RooDataHist wInt = " << wInt << endl ;
  return wInt ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return the error on current weight

void RooDataHist::weightError(Double_t& lo, Double_t& hi, ErrorType etype) const 
{ 
  checkInit() ;

  switch (etype) {

  case Auto:
    throw string(Form("RooDataHist::weightError(%s) error type Auto not allowed here",GetName())) ;
    break ;

  case Expected:
    throw string(Form("RooDataHist::weightError(%s) error type Expected not allowed here",GetName())) ;
    break ;

  case Poisson:
    if (_curWgtErrLo>=0) {
      // Weight is preset or precalculated    
      lo = _curWgtErrLo ;
      hi = _curWgtErrHi ;
      return ;
    }
    
    // Calculate poisson errors
    Double_t ym,yp ;  
    RooHistError::instance().getPoissonInterval(Int_t(weight()+0.5),ym,yp,1) ;
    _curWgtErrLo = weight()-ym ;
    _curWgtErrHi = yp-weight() ;
    _errLo[_curIndex] = _curWgtErrLo ;
    _errHi[_curIndex] = _curWgtErrHi ;
    lo = _curWgtErrLo ;
    hi = _curWgtErrHi ;
    return ;

  case SumW2:
    lo = sqrt(_curSumW2) ;
    hi = sqrt(_curSumW2) ;
    return ;

  case None:
    lo = 0 ;
    hi = 0 ;
    return ;
  }
}


// wve adjust for variable bin sizes

////////////////////////////////////////////////////////////////////////////////
/// Perform boundary safe 'intOrder'-th interpolation of weights in dimension 'dim'
/// at current value 'xval'

Double_t RooDataHist::interpolateDim(RooRealVar& dim, const RooAbsBinning* binning, Double_t xval, Int_t intOrder, Bool_t correctForBinSize, Bool_t cdfBoundaries) 
{
  // Fill workspace arrays spanning interpolation area
  Int_t fbinC = dim.getBin(*binning) ;
  Int_t fbinLo = fbinC-intOrder/2 - ((xval<binning->binCenter(fbinC))?1:0) ;
  Int_t fbinM = dim.numBins(*binning) ;


  Int_t i ;
  Double_t yarr[10] ;
  Double_t xarr[10] ;
  for (i=fbinLo ; i<=intOrder+fbinLo ; i++) {
    Int_t ibin ;
    if (i>=0 && i<fbinM) {
      // In range
      ibin = i ;
      dim.setBinFast(ibin,*binning) ;
      //cout << "INRANGE: dim.getVal(ibin=" << ibin << ") = " << dim.getVal() << endl ;
      xarr[i-fbinLo] = dim.getVal() ;
      Int_t idx = calcTreeIndex();
      yarr[i - fbinLo] = get_wgt(idx);
      if (correctForBinSize) yarr[i-fbinLo] /=  _binv[idx] ;
    } else if (i>=fbinM) {
      // Overflow: mirror
      ibin = 2*fbinM-i-1 ;
      dim.setBinFast(ibin,*binning) ;
      //cout << "OVERFLOW: dim.getVal(ibin=" << ibin << ") = " << dim.getVal() << endl ;
      if (cdfBoundaries) {	
	xarr[i-fbinLo] = dim.getMax()+1e-10*(i-fbinM+1) ;
	yarr[i-fbinLo] = 1.0 ;
      } else {
	Int_t idx = calcTreeIndex() ;      
	xarr[i-fbinLo] = 2*dim.getMax()-dim.getVal() ;
   yarr[i - fbinLo] = get_wgt(idx);
   if (correctForBinSize)
      yarr[i - fbinLo] /= _binv[idx];
      }
    } else {
      // Underflow: mirror
      ibin = -i - 1 ;
      dim.setBinFast(ibin,*binning) ;
      //cout << "UNDERFLOW: dim.getVal(ibin=" << ibin << ") = " << dim.getVal() << endl ;
      if (cdfBoundaries) {
	xarr[i-fbinLo] = dim.getMin()-ibin*(1e-10) ; ;
	yarr[i-fbinLo] = 0.0 ;
      } else {
	Int_t idx = calcTreeIndex() ;      
	xarr[i-fbinLo] = 2*dim.getMin()-dim.getVal() ;
   yarr[i - fbinLo] = get_wgt(idx);
   if (correctForBinSize)
      yarr[i - fbinLo] /= _binv[idx];
      }
    }
    //cout << "ibin = " << ibin << endl ;
  }
//   for (int k=0 ; k<=intOrder ; k++) {
//     cout << "k=" << k << " x = " << xarr[k] << " y = " << yarr[k] << endl ;
//   }
  dim.setBinFast(fbinC,*binning) ;
  Double_t ret = RooMath::interpolate(xarr,yarr,intOrder+1,xval) ;
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Increment the weight of the bin enclosing the coordinates given
/// by 'row' by the specified amount. Add the sum of weights squared
/// for the bin by 'sumw2' rather than wgt^2

void RooDataHist::add(const RooArgSet& row, Double_t wgt, Double_t sumw2) 
{
  checkInit() ;

//   cout << "RooDataHist::add() accepted coordinate is " << endl ;
//   _vars.Print("v") ;

  _vars = row ;
  Int_t idx = calcTreeIndex() ;
  _wgt[idx] += wgt ;  
  _sumw2[idx] += (sumw2>0?sumw2:wgt*wgt) ;
  _errLo[idx] = -1 ;
  _errHi[idx] = -1 ;

  _cache_sum_valid = kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set the weight and errors of the bin enclosing the coordinates
/// given by `row`. Associate errors [wgtErrLo,wgtErrHi] with the event weight on this bin.
void RooDataHist::set(const RooArgSet& row, Double_t wgt, Double_t wgtErrLo, Double_t wgtErrHi) 
{
  checkInit() ;

  _vars = row ;
  Int_t idx = calcTreeIndex() ;
  _wgt[idx] = wgt ;  
  _errLo[idx] = wgtErrLo ;  
  _errHi[idx] = wgtErrHi ;  

  _cache_sum_valid = kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set the weight and weight error of the bin enclosing the current (i.e. last-used)
/// coordinates.
void RooDataHist::set(Double_t wgt, Double_t wgtErr) 
{
  checkInit() ;
  
  if (_curIndex<0) {
    _curIndex = calcTreeIndex() ;
  }

  _wgt[_curIndex] = wgt ;  
  _errLo[_curIndex] = wgtErr ;  
  _errHi[_curIndex] = wgtErr ;  
  _sumw2[_curIndex] = wgtErr*wgtErr ;

  _cache_sum_valid = kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set the weight and weight error of the bin enclosing the coordinates
/// given by `row`.
void RooDataHist::set(const RooArgSet& row, Double_t wgt, Double_t wgtErr) 
{
  checkInit() ;

  _vars = row ;
  Int_t idx = calcTreeIndex() ;
  _wgt[idx] = wgt ;  
  _errLo[idx] = wgtErr ;  
  _errHi[idx] = wgtErr ;  
  _sumw2[idx] = wgtErr*wgtErr ;

  _cache_sum_valid = kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add all data points contained in 'dset' to this data set with given weight.
/// Optional cut string expression selects the data points to be added and can
/// reference any variable contained in this data set

void RooDataHist::add(const RooAbsData& dset, const char* cut, Double_t wgt) 
{  
  RooFormulaVar cutVar("select",cut,*dset.get()) ;
  add(dset,&cutVar,wgt) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add all data points contained in 'dset' to this data set with given weight.
/// Optional RooFormulaVar pointer selects the data points to be added.

void RooDataHist::add(const RooAbsData& dset, const RooFormulaVar* cutVar, Double_t wgt) 
{
  checkInit() ;

  RooFormulaVar* cloneVar = 0;
  RooArgSet* tmp(0) ;
  if (cutVar) {
    // Deep clone cutVar and attach clone to this dataset
    tmp = (RooArgSet*) RooArgSet(*cutVar).snapshot() ;
    if (!tmp) {
      coutE(DataHandling) << "RooDataHist::add(" << GetName() << ") Couldn't deep-clone cut variable, abort," << endl ;
      return ;
    }

    cloneVar = (RooFormulaVar*) tmp->find(*cutVar) ;
    cloneVar->attachDataSet(dset) ;
  }


  Int_t i ;
  for (i=0 ; i<dset.numEntries() ; i++) {
    const RooArgSet* row = dset.get(i) ;
    if (!cloneVar || cloneVar->getVal()) {
       add(*row,wgt*dset.weight(), wgt*wgt*dset.weightSquared()) ;
    }
  }

  if (cloneVar) {
    delete tmp ;
  } 

  _cache_sum_valid = kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the sum of the weights of all hist bins.
///
/// If correctForBinSize is specified, the sum of weights
/// is multiplied by the N-dimensional bin volume,
/// making the return value the integral over the function
/// represented by this histogram

Double_t RooDataHist::sum(Bool_t correctForBinSize, Bool_t inverseBinCor) const 
{
  checkInit() ;

  // Check if result was cached
  Int_t cache_code = 1 + (correctForBinSize?1:0) + ((correctForBinSize&&inverseBinCor)?1:0) ;
  if (_cache_sum_valid==cache_code) {
    return _cache_sum ;
  }

  Int_t i ;
  Double_t total(0), carry(0);
  for (i=0 ; i<_arrSize ; i++) {
    
    Double_t theBinVolume = correctForBinSize ? (inverseBinCor ? 1/_binv[i] : _binv[i]) : 1.0 ;
    // cout << "total += " << _wgt[i] << "*" << theBinVolume << endl ;
    // Double_t y = _wgt[i]*theBinVolume - carry;
    Double_t y = get_wgt(i) * theBinVolume - carry;
    Double_t t = total + y;
    carry = (t - total) - y;
    total = t;
  }

  // Store result in cache
  _cache_sum_valid=cache_code ;
  _cache_sum = total ;

  return total ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the sum of the weights of a multi-dimensional slice of the histogram
/// by summing only over the dimensions specified in sumSet.
///   
/// The coordinates of all other dimensions are fixed to those given in sliceSet
///
/// If correctForBinSize is specified, the sum of weights
/// is multiplied by the M-dimensional bin volume, (M = N(sumSet)),
/// making the return value the integral over the function
/// represented by this histogram

Double_t RooDataHist::sum(const RooArgSet& sumSet, const RooArgSet& sliceSet, Bool_t correctForBinSize, Bool_t inverseBinCor)
{
  checkInit() ;

  RooArgSet varSave ;
  varSave.addClone(_vars) ;

  RooArgSet* sliceOnlySet = new RooArgSet(sliceSet) ;
  sliceOnlySet->remove(sumSet,kTRUE,kTRUE) ;

  _vars = *sliceOnlySet ;
  calculatePartialBinVolume(*sliceOnlySet) ;
  delete sliceOnlySet ;
  
  // Calculate mask and refence plot bins for non-iterating variables
  Bool_t* mask = new Bool_t[_vars.getSize()] ;
  Int_t*  refBin = new Int_t[_vars.getSize()] ;

  for (unsigned int i = 0; i < _vars.size(); ++i) {
    const auto arg = _vars[i];

    if (sumSet.find(*arg)) {
      mask[i] = kFALSE ;
    } else {
      mask[i] = kTRUE ;
      refBin[i] = dynamic_cast<RooAbsLValue*>(arg)->getBin();
    }
  }
    
  // Loop over entire data set, skipping masked entries
  Double_t total(0), carry(0);
  Int_t ibin ;
  for (ibin=0 ; ibin<_arrSize ; ibin++) {

    Int_t idx(0), tmp(ibin), ivar(0) ;
    Bool_t skip(kFALSE) ;

    // Check if this bin belongs in selected slice
    for (unsigned int i = 0; !skip && i < _vars.size(); ++i) {
      idx  = tmp / _idxMult[ivar] ;
      tmp -= idx*_idxMult[ivar] ;
      if (mask[ivar] && idx!=refBin[ivar]) skip=kTRUE ;
      ivar++ ;
    }
    
    if (!skip) {
      Double_t theBinVolume = correctForBinSize ? (inverseBinCor ? 1/(*_pbinv)[_vars.size()] : (*_pbinv)[_vars.size()] ) : 1.0 ;
      //       cout << "adding bin[" << ibin << "] to sum wgt = " << _wgt[ibin] << " binv = " << theBinVolume << endl ;
      // Double_t y = _wgt[ibin]*theBinVolume - carry;
      Double_t y = get_wgt(ibin) * theBinVolume - carry;
      Double_t t = total + y;
      carry = (t - total) - y;
      total = t;
    }
  }

  delete[] mask ;
  delete[] refBin ;

  _vars = varSave ;

  return total ;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the sum of the weights of a multi-dimensional slice of the histogram
/// by summing only over the dimensions specified in sumSet.
///   
/// The coordinates of all other dimensions are fixed to those given in sliceSet
///
/// If correctForBinSize is specified, the sum of weights
/// is multiplied by the M-dimensional bin volume, (M = N(sumSet)),
/// or the fraction of it that falls inside the range rangeName,
/// making the return value the integral over the function
/// represented by this histogram.
///
/// If correctForBinSize is not specified, the weights are multiplied by the
/// fraction of the bin volume that falls inside the range, i.e. a factor of
/// binVolumeInRange/totalBinVolume.

Double_t RooDataHist::sum(const RooArgSet& sumSet, const RooArgSet& sliceSet,
	Bool_t correctForBinSize, Bool_t inverseBinCor,
	const std::map<const RooAbsArg*, std::pair<Double_t, Double_t> >& ranges)
{
  checkInit();
  checkBinBounds();
  RooArgSet varSave;
  varSave.addClone(_vars);
  {
    RooArgSet sliceOnlySet(sliceSet);
    sliceOnlySet.remove(sumSet,kTRUE,kTRUE);
    _vars = sliceOnlySet;
  }

  // Calculate mask and reference plot bins for non-iterating variables,
  // and get ranges for iterating variables
  std::vector<bool> mask(_vars.getSize());
  std::vector<Int_t> refBin(_vars.getSize());
  std::vector<Double_t> rangeLo(_vars.getSize(), -std::numeric_limits<Double_t>::infinity());
  std::vector<Double_t> rangeHi(_vars.getSize(), +std::numeric_limits<Double_t>::infinity());

  for (std::size_t i = 0; i < _vars.size(); ++i) {
    const auto arg = _vars[i];
    RooAbsArg* sumsetv = sumSet.find(*arg);
    RooAbsArg* slicesetv = sliceSet.find(*arg);
    mask[i] = !sumsetv;
    if (mask[i]) {
      auto argLV = dynamic_cast<const RooAbsLValue*>(arg);
      assert(argLV);
      refBin[i] = argLV->getBin();
    }

	auto it = ranges.find(sumsetv ? sumsetv : slicesetv);
    if (ranges.end() != it) {
      rangeLo[i] = it->second.first;
      rangeHi[i] = it->second.second;
    }
  }

  // Loop over entire data set, skipping masked entries
  Double_t total(0), carry(0);
  for (Int_t ibin = 0; ibin < _arrSize; ++ibin) {
    // Check if this bin belongs in selected slice
    Bool_t skip(kFALSE);
    for (int ivar = 0, tmp = ibin; !skip && ivar < int(_vars.size()); ++ivar) {
      const Int_t idx = tmp / _idxMult[ivar];
      tmp -= idx*_idxMult[ivar];
      if (mask[ivar] && idx!=refBin[ivar]) skip=kTRUE;
    }

    if (skip) continue;

    // work out bin volume
    Double_t theBinVolume = 1.;
    for (Int_t ivar = 0, tmp = ibin; ivar < (int)_vars.size(); ++ivar) {
      const Int_t idx = tmp / _idxMult[ivar];
      tmp -= idx*_idxMult[ivar];
      if (_binbounds[ivar].empty()) continue;
      const Double_t binLo = _binbounds[ivar][2 * idx];
      const Double_t binHi = _binbounds[ivar][2 * idx + 1];
      if (binHi < rangeLo[ivar] || binLo > rangeHi[ivar]) {
        // bin is outside of allowed range - effective bin volume is zero
        theBinVolume = 0.;
        break;
      }
      theBinVolume *= 
          (std::min(rangeHi[ivar], binHi) - std::max(rangeLo[ivar], binLo));
    }
    const Double_t corrPartial = theBinVolume / _binv[ibin];
    if (0. == corrPartial) continue;
    const Double_t corr = correctForBinSize ? (inverseBinCor ? 1. / _binv[ibin] : _binv[ibin] ) : 1.0;
    //cout << "adding bin[" << ibin << "] to sum wgt = " << _wgt[ibin] << " binv = " << theBinVolume << " _binv[" << ibin << "] " << _binv[ibin] << endl;
    // const Double_t y = _wgt[ibin] * corr * corrPartial - carry;
    const Double_t y = get_wgt(ibin) * corr * corrPartial - carry;
    const Double_t t = total + y;
    carry = (t - total) - y;
    total = t;
  }

  _vars = varSave;

  return total;
}



////////////////////////////////////////////////////////////////////////////////
/// Fill the transient cache with partial bin volumes with up-to-date
/// values for the partial volume specified by observables 'dimSet'

void RooDataHist::calculatePartialBinVolume(const RooArgSet& dimSet) const 
{
  // Allocate cache if not yet existing
  vector<Double_t> *pbinv = _pbinvCacheMgr.getObj(&dimSet) ;
  if (pbinv) {
    _pbinv = pbinv ;
    return ;
  }

  pbinv = new vector<Double_t>(_arrSize) ;

  // Calculate plot bins of components from master index
  Bool_t* selDim = new Bool_t[_vars.getSize()] ;
  Int_t i(0) ;
  for (const auto v : _vars) {
    selDim[i++] = dimSet.find(*v) ? kTRUE : kFALSE ;
  }

  // Recalculate partial bin volume cache
  Int_t ibin ;
  for (ibin=0 ; ibin<_arrSize ; ibin++) {
    Int_t j(0), idx(0), tmp(ibin) ;
    Double_t theBinVolume(1) ;
    for (const auto absArg : _vars) {
      auto arg = dynamic_cast<const RooAbsLValue*>(absArg);
      if (!arg)
        break;

      idx  = tmp / _idxMult[j] ;
      tmp -= idx*_idxMult[j++] ;
      if (selDim[j-1]) {
        theBinVolume *= arg->getBinWidth(idx) ;
      }
    }
    (*pbinv)[ibin] = theBinVolume ;
  }

  delete[] selDim ;

  // Put in cache (which takes ownership) 
  _pbinvCacheMgr.setObj(&dimSet,pbinv) ;

  // Publicize the array
  _pbinv = pbinv ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the number of bins

Int_t RooDataHist::numEntries() const 
{
  return RooAbsData::numEntries() ;
}



////////////////////////////////////////////////////////////////////////////////

Double_t RooDataHist::sumEntries() const 
{
  Int_t i ;
  Double_t n(0), carry(0);
  for (i=0 ; i<_arrSize ; i++) {
    if (!_binValid || _binValid[i]) {
       // Double_t y = _wgt[i] - carry;
       Double_t y = get_wgt(i) - carry;
       Double_t t = n + y;
       carry = (t - n) - y;
       n = t;
    }
  }
  return n ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the sum of weights in all entries matching cutSpec (if specified)
/// and in named range cutRange (if specified)
/// Return the

Double_t RooDataHist::sumEntries(const char* cutSpec, const char* cutRange) const
{
  checkInit() ;

  if (cutSpec==0 && cutRange==0) {
    return sumEntries();
  } else {
    
    // Setup RooFormulaVar for cutSpec if it is present
    RooFormula* select = 0 ;
    if (cutSpec) {
      select = new RooFormula("select",cutSpec,*get()) ;
    }
    
    // Otherwise sum the weights in the event
    Double_t sumw(0), carry(0);
    Int_t i ;
    for (i=0 ; i<numEntries() ; i++) {
      get(i) ;
      if (select && select->eval()==0.) continue ;
      if (cutRange && !_vars.allInRange(cutRange)) continue ;

      if (!_binValid || _binValid[i]) {
	Double_t y = weight() - carry;
	Double_t t = sumw + y;
	carry = (t - sumw) - y;
	sumw = t;
      }
    }
    
    if (select) delete select ;
    
    return sumw ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Reset all bin weights to zero

void RooDataHist::reset() 
{
  // WVE DO NOT CALL RooTreeData::reset() for binned
  // datasets as this will delete the bin definitions

  Int_t i ;
  for (i=0 ; i<_arrSize ; i++) {
    _wgt[i] = 0. ;
    _errLo[i] = -1 ;
    _errHi[i] = -1 ;
  }
  _curWeight = 0 ;
  _curWgtErrLo = -1 ;
  _curWgtErrHi = -1 ;
  _curVolume = 1 ;

  _cache_sum_valid = kFALSE ;

}



////////////////////////////////////////////////////////////////////////////////
/// Return an argset with the bin center coordinates for 
/// bin sequential number 'masterIdx'. For iterative use.

const RooArgSet* RooDataHist::get(Int_t masterIdx) const  
{
  checkInit() ;
  _curWeight = _wgt[masterIdx] ;
  _curWgtErrLo = _errLo[masterIdx] ;
  _curWgtErrHi = _errHi[masterIdx] ;
  _curSumW2 = _sumw2[masterIdx] ;
  _curVolume = _binv[masterIdx] ; 
  _curIndex  = masterIdx ;
  return RooAbsData::get(masterIdx) ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Return a RooArgSet with center coordinates of the bin
/// enclosing the point 'coord'

const RooArgSet* RooDataHist::get(const RooArgSet& coord) const
{
  ((RooDataHist*)this)->_vars = coord ;
  return get(calcTreeIndex()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the volume of the bin enclosing coordinates 'coord'

Double_t RooDataHist::binVolume(const RooArgSet& coord) 
{
  checkInit() ;
  ((RooDataHist*)this)->_vars = coord ;
  return _binv[calcTreeIndex()] ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set all the event weight of all bins to the specified value

void RooDataHist::setAllWeights(Double_t value) 
{
  for (Int_t i=0 ; i<_arrSize ; i++) {
    _wgt[i] = value ;
  }

  _cache_sum_valid = kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Create an iterator over all bins in a slice defined by the subset of observables
/// listed in sliceArg. The position of the slice is given by otherArgs

TIterator* RooDataHist::sliceIterator(RooAbsArg& sliceArg, const RooArgSet& otherArgs) 
{
  // Update to current position
  _vars = otherArgs ;
  _curIndex = calcTreeIndex() ;
  
  RooAbsArg* intArg = _vars.find(sliceArg) ;
  if (!intArg) {
    coutE(InputArguments) << "RooDataHist::sliceIterator() variable " << sliceArg.GetName() << " is not part of this RooDataHist" << endl ;
    return 0 ;
  }
  return new RooDataHistSliceIter(*this,*intArg) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Change the name of the RooDataHist

void RooDataHist::SetName(const char *name) 
{
  if (_dir) _dir->GetList()->Remove(this);
  TNamed::SetName(name) ;
  if (_dir) _dir->GetList()->Add(this);
}


////////////////////////////////////////////////////////////////////////////////
/// Change the title of this RooDataHist

void RooDataHist::SetNameTitle(const char *name, const char* title) 
{
  if (_dir) _dir->GetList()->Remove(this);
  TNamed::SetNameTitle(name,title) ;
  if (_dir) _dir->GetList()->Add(this);
}


////////////////////////////////////////////////////////////////////////////////
/// Print value of the dataset, i.e. the sum of weights contained in the dataset

void RooDataHist::printValue(ostream& os) const 
{
  os << numEntries() << " bins (" << sumEntries() << " weights)" ;
}




////////////////////////////////////////////////////////////////////////////////
/// Print argument of dataset, i.e. the observable names

void RooDataHist::printArgs(ostream& os) const 
{
  os << "[" ;    
  Bool_t first(kTRUE) ;
  for (const auto arg : _vars) {
    if (first) {
      first=kFALSE ;
    } else {
      os << "," ;
    }
    os << arg->GetName() ;
  }
  os << "]" ;
}



////////////////////////////////////////////////////////////////////////////////
/// Cache the datahist entries with bin centers that are inside/outside the
/// current observable definitio

void RooDataHist::cacheValidEntries() 
{
  checkInit() ;

  if (!_binValid) {
    _binValid = new Bool_t[_arrSize] ;
  }
  TIterator* iter = _vars.createIterator() ;
  RooAbsArg* arg ;
  for (Int_t i=0 ; i<_arrSize ; i++) {
    get(i) ;
    _binValid[i] = kTRUE ;
    iter->Reset() ;
    while((arg=(RooAbsArg*)iter->Next())) {
      // coverity[CHECKED_RETURN]
      _binValid[i] &= arg->inRange(0) ;      
    }
  }
  delete iter ;

}


////////////////////////////////////////////////////////////////////////////////
/// Return true if currently loaded coordinate is considered valid within
/// the current range definitions of all observables

Bool_t RooDataHist::valid() const 
{
  // If caching is enabled, use the precached result
  if (_binValid) {
    return _binValid[_curIndex] ;
  }

  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Returns true if datasets contains entries with a non-integer weight

Bool_t RooDataHist::isNonPoissonWeighted() const
{
  for (int i=0 ; i<numEntries() ; i++) {
    if (fabs(_wgt[i]-Int_t(_wgt[i]))>1e-10) return kTRUE ;
  }
  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Print the details on the dataset contents

void RooDataHist::printMultiline(ostream& os, Int_t content, Bool_t verbose, TString indent) const 
{
  RooAbsData::printMultiline(os,content,verbose,indent) ;  

  os << indent << "Binned Dataset " << GetName() << " (" << GetTitle() << ")" << endl ;
  os << indent << "  Contains " << numEntries() << " bins with a total weight of " << sumEntries() << endl;
  
  if (!verbose) {
    os << indent << "  Observables " << _vars << endl ;
  } else {
    os << indent << "  Observables: " ;
    _vars.printStream(os,kName|kValue|kExtras|kTitle,kVerbose,indent+"  ") ;
  }

  if(verbose) {
    if (_cachedVars.getSize()>0) {
      os << indent << "  Caches " << _cachedVars << endl ;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooDataHist.

void RooDataHist::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {

     UInt_t R__s, R__c;
     Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
     
      if (R__v>2) {

	R__b.ReadClassBuffer(RooDataHist::Class(),this,R__v,R__s,R__c);
	initialize(0,kFALSE) ;

      } else {	

	// Legacy dataset conversion happens here. Legacy RooDataHist inherits from RooTreeData
	// which in turn inherits from RooAbsData. Manually stream RooTreeData contents on 
	// file here and convert it into a RooTreeDataStore which is installed in the 
	// new-style RooAbsData base class
	
	// --- This is the contents of the streamer code of RooTreeData version 2 ---
	UInt_t R__s1, R__c1;
	Version_t R__v1 = R__b.ReadVersion(&R__s1, &R__c1); if (R__v1) { }
	
	RooAbsData::Streamer(R__b);
	TTree* X_tree(0) ; R__b >> X_tree;
	RooArgSet X_truth ; X_truth.Streamer(R__b);
	TString X_blindString ; X_blindString.Streamer(R__b);
	R__b.CheckByteCount(R__s1, R__c1, RooTreeData::Class());
	// --- End of RooTreeData-v1 streamer
	
	// Construct RooTreeDataStore from X_tree and complete initialization of new-style RooAbsData
	_dstore = new RooTreeDataStore(X_tree,_vars) ;
	_dstore->SetName(GetName()) ;
	_dstore->SetTitle(GetTitle()) ;
	_dstore->checkInit() ;       
	
	RooDirItem::Streamer(R__b);
	R__b >> _arrSize;
	delete [] _wgt;
	_wgt = new Double_t[_arrSize];
	R__b.ReadFastArray(_wgt,_arrSize);
	delete [] _errLo;
	_errLo = new Double_t[_arrSize];
	R__b.ReadFastArray(_errLo,_arrSize);
	delete [] _errHi;
	_errHi = new Double_t[_arrSize];
	R__b.ReadFastArray(_errHi,_arrSize);
	delete [] _sumw2;
	_sumw2 = new Double_t[_arrSize];
	R__b.ReadFastArray(_sumw2,_arrSize);
	delete [] _binv;
	_binv = new Double_t[_arrSize];
	R__b.ReadFastArray(_binv,_arrSize);
	_realVars.Streamer(R__b);
	R__b >> _curWeight;
	R__b >> _curWgtErrLo;
	R__b >> _curWgtErrHi;
	R__b >> _curSumW2;
	R__b >> _curVolume;
	R__b >> _curIndex;
	R__b.CheckByteCount(R__s, R__c, RooDataHist::IsA());

      }
      
   } else {

      R__b.WriteClassBuffer(RooDataHist::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return event weights of all events in range [first, first+len).
/// If no contiguous structure of weights is stored, an empty batch can be returned.
/// This indicates that the weight is constant. Use weight() to retrieve it.
RooSpan<const double> RooDataHist::getWeightBatch(std::size_t first, std::size_t len) const {
  if (first > static_cast<std::size_t>(_arrSize))
    return {};

  return {&_wgt[first], std::min(_arrSize - first + len, len)};
}


////////////////////////////////////////////////////////////////////////////////
/// Write information to retrieve data columns into `evalData.spans`.
/// All spans belonging to variables of this dataset are overwritten. Spans to other
/// variables remain intact.
/// \param[out] evalData Store references to all data batches in this struct's `spans`.
/// The key to retrieve an item is the pointer of the variable that owns the data.
/// \param first Index of first event that ends up in the batch.
/// \param len   Number of events in each batch.
void RooDataHist::getBatches(BatchHelpers::RunContext& evalData, std::size_t begin, std::size_t len) const {
  for (auto&& batch : store()->getBatches(begin, len).spans) {
    evalData.spans[batch.first] = std::move(batch.second);
  }
}

