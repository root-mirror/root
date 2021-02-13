// @(#)root/roostats:$Id$
// Author: George Lewis, Kyle Cranmer
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOBARLOWBEESTONLL
#define ROOBARLOWBEESTONLL

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include <map>
#include <set>
#include <string>
#include <vector>

class RooMinuit ;

namespace RooStats{
  namespace HistFactory{

class RooBarlowBeestonLL : public RooAbsReal {
public:

  RooBarlowBeestonLL() ;
  RooBarlowBeestonLL(const char *name, const char *title, RooAbsReal& nll /*, const RooArgSet& observables*/);
  RooBarlowBeestonLL(const RooBarlowBeestonLL& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBarlowBeestonLL(*this,newname); }
  virtual ~RooBarlowBeestonLL() ;

  // A simple class to store the 
  // necessary objects for a 
  // single gamma in a single channel
  class BarlowCache {
  public:
    BarlowCache() : hasStatUncert(false), gamma(NULL), 
		    observables(NULL), bin_center(NULL), 
		    tau(NULL), nom_pois_mean(NULL),
		    sumPdf(NULL),  nData(-1) {}
    bool hasStatUncert;
    RooRealVar* gamma;
    RooArgSet* observables;
    RooArgSet* bin_center; // Snapshot
    RooRealVar* tau;
    RooAbsReal* nom_pois_mean;
    RooAbsReal* sumPdf;
    double nData;
    double binVolume;
    void SetBinCenter() const;
    /*
    // Restore original values and constant status of observables
    TIterator* iter = obsSetOrig->createIterator() ;
    RooRealVar* var ;
    while((var=(RooRealVar*)iter->Next())) {
    RooRealVar* target = (RooRealVar*) _obs.find(var->GetName()) ;
    target->setVal(var->getVal()) ;
    target->setConstant(var->isConstant()) ;
    }
     */

  };
  

  void initializeBarlowCache();

  RooArgSet* getParameters(const RooArgSet* depList, Bool_t stripDisconnected=kTRUE) const;

  // void setAlwaysStartFromMin(Bool_t flag) { _startFromMin = flag ; }
  // Bool_t alwaysStartFromMin() const { return _startFromMin ; }

  //RooMinuit* minuit() { return _minuit ; }
  RooAbsReal& nll() { return const_cast<RooAbsReal&>(_nll.arg()) ; }
  // const RooArgSet& bestFitParams() const ;
  // const RooArgSet& bestFitObs() const ;

  //  virtual RooAbsReal* createProfile(const RooArgSet& paramsOfInterest) ;
  
  virtual Bool_t redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) ;

  // void clearAbsMin() { _absMinValid = kFALSE ; }

  // Int_t numEval() const { return _neval ; }

  void setPdf(RooAbsPdf* pdf) { _pdf = pdf; }
  void setDataset(RooAbsData* data) { _data = data; }
  
  //void FactorizePdf(const RooArgSet &observables, RooAbsPdf &pdf, 
  //	    RooArgList &obsTerms, RooArgList &constraints) const;


protected:

  // void validateAbsMin() const ;


  RooRealProxy _nll ;    // Input -log(L) function
  /*
  RooSetProxy _obs ;     // Parameters of profile likelihood
  RooSetProxy _par ;     // Marginialized parameters of likelihood
  */
  RooAbsPdf* _pdf;
  RooAbsData* _data;
  mutable std::map< std::string, std::vector< BarlowCache > > _barlowCache;
  mutable std::set< std::string > _statUncertParams;
  // Bool_t _startFromMin ; // Always start minimization for global minimum?

  /*
  TIterator* _piter ; //! Iterator over profile likelihood parameters to be minimized 
  TIterator* _oiter ; //! Iterator of profile likelihood output parameter(s)
  */

  // mutable RooMinuit* _minuit ; //! Internal minuit instance

  // mutable Bool_t _absMinValid ; // flag if absmin is up-to-date
  // mutable Double_t _absMin ; // absolute minimum of -log(L)
  // mutable RooArgSet _paramAbsMin ; // Parameter values at absolute minimum
  // mutable RooArgSet _obsAbsMin ; // Observable values at absolute minimum
  mutable std::map<std::string,bool> _paramFixed ; // Parameter constant status at last time of use
  // mutable Int_t _neval ; // Number evaluations used in last minimization
  Double_t evaluate() const ;
  //Double_t evaluate_bad() const ;


private:

  ClassDef(RooStats::HistFactory::RooBarlowBeestonLL,0) // Real-valued function representing a Barlow-Beeston minimized profile likelihood of external (likelihood) function
};

  }
}
 
#endif
