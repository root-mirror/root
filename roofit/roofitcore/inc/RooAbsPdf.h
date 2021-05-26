/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsPdf.h,v 1.90 2007/07/21 21:32:52 wouter Exp $
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
#ifndef ROO_ABS_PDF
#define ROO_ABS_PDF

#include "RooAbsReal.h"
//#include "RooRealIntegral.h"
#include "RooNameSet.h"
#include "RooObjCacheManager.h"
#include "RooCmdArg.h"

class RooDataSet;
class RooDataHist ;
class RooArgSet ;
class RooAbsGenContext ;
class RooFitResult ;
class RooExtendPdf ;
class RooCategory ;
class TPaveText;
class TH1F;
class TH2F;
class TList ;
class RooLinkedList ;
class RooNumGenConfig ;
class RooRealIntegral ;
namespace RooBatchCompute {
struct RunContext;
}

class RooAbsPdf : public RooAbsReal {
public:

  // Constructors, assignment etc
  RooAbsPdf() ;
  RooAbsPdf(const char *name, const char *title=0) ;
  RooAbsPdf(const char *name, const char *title, Double_t minVal, Double_t maxVal) ;
  // RooAbsPdf(const RooAbsPdf& other, const char* name=0);
  virtual ~RooAbsPdf();

  // Toy MC generation

  ////////////////////////////////////////////////////////////////////////////////
  /// See RooAbsPdf::generate(const RooArgSet&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&)
  /// \param[in] nEvents How many events to generate
  RooDataSet *generate(const RooArgSet &whatVars, Int_t nEvents, const RooCmdArg& arg1,
                       const RooCmdArg& arg2=RooCmdArg::none(), const RooCmdArg& arg3=RooCmdArg::none(),
                       const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none()) {
    return generate(whatVars,RooFit::NumEvents(nEvents),arg1,arg2,arg3,arg4,arg5) ;
  }
  RooDataSet *generate(const RooArgSet &whatVars,
                       const RooCmdArg& arg1=RooCmdArg::none(),const RooCmdArg& arg2=RooCmdArg::none(),
                       const RooCmdArg& arg3=RooCmdArg::none(),const RooCmdArg& arg4=RooCmdArg::none(),
                       const RooCmdArg& arg5=RooCmdArg::none(),const RooCmdArg& arg6=RooCmdArg::none()) ;
  RooDataSet *generate(const RooArgSet &whatVars, Double_t nEvents = 0, Bool_t verbose=kFALSE, Bool_t autoBinned=kTRUE,
             const char* binnedTag="", Bool_t expectedData=kFALSE, Bool_t extended = kFALSE) const;
  RooDataSet *generate(const RooArgSet &whatVars, const RooDataSet &prototype, Int_t nEvents= 0,
             Bool_t verbose=kFALSE, Bool_t randProtoOrder=kFALSE, Bool_t resampleProto=kFALSE) const;


  class GenSpec {
  public:
    virtual ~GenSpec() ;
    GenSpec() { _genContext = 0 ; _protoData = 0 ; _init = kFALSE ; _extended=kFALSE, _nGen=0 ; _randProto = kFALSE ; _resampleProto=kFALSE ; }
  private:
    GenSpec(RooAbsGenContext* context, const RooArgSet& whatVars, RooDataSet* protoData, Int_t nGen, Bool_t extended,
       Bool_t randProto, Bool_t resampleProto, TString dsetName, Bool_t init=kFALSE) ;
    GenSpec(const GenSpec& other) ;

    friend class RooAbsPdf ;
    RooAbsGenContext* _genContext ;
    RooArgSet _whatVars ;
    RooDataSet* _protoData ;
    Int_t _nGen ;
    Bool_t _extended ;
    Bool_t _randProto ;
    Bool_t _resampleProto ;
    TString _dsetName ;
    Bool_t _init ;
    ClassDef(GenSpec,0) // Generation specification
  } ;

  ///Prepare GenSpec configuration object for efficient generation of multiple datasets from identical specification.
  GenSpec* prepareMultiGen(const RooArgSet &whatVars,
            const RooCmdArg& arg1=RooCmdArg::none(),const RooCmdArg& arg2=RooCmdArg::none(),
            const RooCmdArg& arg3=RooCmdArg::none(),const RooCmdArg& arg4=RooCmdArg::none(),
            const RooCmdArg& arg5=RooCmdArg::none(),const RooCmdArg& arg6=RooCmdArg::none()) ;
  ///Generate according to GenSpec obtained from prepareMultiGen().
  RooDataSet* generate(GenSpec&) const ;


  ////////////////////////////////////////////////////////////////////////////////
  /// As RooAbsPdf::generateBinned(const RooArgSet&, const RooCmdArg&,const RooCmdArg&, const RooCmdArg&,const RooCmdArg&, const RooCmdArg&,const RooCmdArg&)
  /// \param[in] nEvents How many events to generate
  virtual RooDataHist *generateBinned(const RooArgSet &whatVars, Double_t nEvents, const RooCmdArg& arg1,
               const RooCmdArg& arg2=RooCmdArg::none(), const RooCmdArg& arg3=RooCmdArg::none(),
               const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none()) const {
    return generateBinned(whatVars,RooFit::NumEvents(nEvents),arg1,arg2,arg3,arg4,arg5);
  }
  virtual RooDataHist *generateBinned(const RooArgSet &whatVars,
               const RooCmdArg& arg1=RooCmdArg::none(),const RooCmdArg& arg2=RooCmdArg::none(),
               const RooCmdArg& arg3=RooCmdArg::none(),const RooCmdArg& arg4=RooCmdArg::none(),
               const RooCmdArg& arg5=RooCmdArg::none(),const RooCmdArg& arg6=RooCmdArg::none()) const;
  virtual RooDataHist *generateBinned(const RooArgSet &whatVars, Double_t nEvents, Bool_t expectedData=kFALSE, Bool_t extended=kFALSE) const;

  virtual RooDataSet* generateSimGlobal(const RooArgSet& whatVars, Int_t nEvents) ;

  ///Helper calling plotOn(RooPlot*, RooLinkedList&) const
  virtual RooPlot* plotOn(RooPlot* frame,
           const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(),
           const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(),
           const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(),
           const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none(),
           const RooCmdArg& arg9=RooCmdArg::none(), const RooCmdArg& arg10=RooCmdArg::none()
              ) const {
    return RooAbsReal::plotOn(frame,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10) ;
  }
  virtual RooPlot* plotOn(RooPlot* frame, RooLinkedList& cmdList) const ;

  /// Add a box with parameter values (and errors) to the specified frame
  virtual RooPlot* paramOn(RooPlot* frame,
                           const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(),
                           const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(),
                           const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(),
                           const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;

  virtual RooPlot* paramOn(RooPlot* frame, const RooAbsData* data, const char *label= "", Int_t sigDigits = 2,
            Option_t *options = "NELU", Double_t xmin=0.65,
            Double_t xmax = 0.9, Double_t ymax = 0.9) ;

  // Built-in generator support
  virtual Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  virtual void initGenerator(Int_t code) ;
  virtual void generateEvent(Int_t code);
  virtual Bool_t isDirectGenSafe(const RooAbsArg& arg) const ;

  // Configuration of MC generators used for this pdf
  const RooNumGenConfig* getGeneratorConfig() const ;
  static RooNumGenConfig* defaultGeneratorConfig()  ;
  RooNumGenConfig* specialGeneratorConfig() const ;
  RooNumGenConfig* specialGeneratorConfig(Bool_t createOnTheFly) ;
  void setGeneratorConfig() ;
  void setGeneratorConfig(const RooNumGenConfig& config) ;

  // -log(L) fits to binned and unbinned data
  virtual RooFitResult* fitTo(RooAbsData& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),
                              const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),
                              const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;
  virtual RooFitResult* fitTo(RooAbsData& data, const RooLinkedList& cmdList) ;

  virtual RooAbsReal* createNLL(RooAbsData& data, const RooLinkedList& cmdList) ;
  virtual RooAbsReal* createNLL(RooAbsData& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),
            const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),
            const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;

  // Chi^2 fits to histograms
  using RooAbsReal::chi2FitTo ;
  using RooAbsReal::createChi2 ;
  virtual RooFitResult* chi2FitTo(RooDataHist& data, const RooLinkedList& cmdList) ;
  virtual RooAbsReal* createChi2(RooDataHist& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),
             const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),
             const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;

  // Chi^2 fits to X-Y datasets
  virtual RooAbsReal* createChi2(RooDataSet& data, const RooLinkedList& cmdList) ;





  // Constraint management
  virtual RooArgSet* getConstraints(const RooArgSet& /*observables*/, RooArgSet& /*constrainedParams*/, Bool_t /*stripDisconnected*/) const {
    // Interface to retrieve constraint terms on this pdf. Default implementation returns null
    return 0 ;
  }
  virtual RooArgSet* getAllConstraints(const RooArgSet& observables, RooArgSet& constrainedParams, Bool_t stripDisconnected=kTRUE) const ;

  // Project p.d.f into lower dimensional p.d.f
  virtual RooAbsPdf* createProjection(const RooArgSet& iset) ;

  // Create cumulative density function from p.d.f
  RooAbsReal* createCdf(const RooArgSet& iset, const RooArgSet& nset=RooArgSet()) ;
  RooAbsReal* createCdf(const RooArgSet& iset, const RooCmdArg& arg1, const RooCmdArg& arg2=RooCmdArg::none(),
         const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(),
         const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(),
         const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;
  RooAbsReal* createScanCdf(const RooArgSet& iset, const RooArgSet& nset, Int_t numScanBins, Int_t intOrder) ;

  // Function evaluation support
  virtual Double_t getValV(const RooArgSet* set=0) const ;
  virtual Double_t getLogVal(const RooArgSet* set=0) const ;

  RooSpan<const double> getValues(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const;
  RooSpan<const double> getLogValBatch(std::size_t begin, std::size_t batchSize,
      const RooArgSet* normSet = nullptr) const;
  RooSpan<const double> getLogProbabilities(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet = nullptr) const;

  /// \copydoc getNorm(const RooArgSet*) const
  Double_t getNorm(const RooArgSet& nset) const {
    return getNorm(&nset) ;
  }
  virtual Double_t getNorm(const RooArgSet* set=0) const ;

  virtual void resetErrorCounters(Int_t resetValue=10) ;
  void setTraceCounter(Int_t value, Bool_t allNodes=kFALSE) ;
private:
  Bool_t traceEvalPdf(Double_t value) const;

public:

  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;

  /// Shows if a PDF is self-normalized, which means that no attempt is made to add a normalization term.
  /// Always returns false, unless a PDF overrides this function.
  virtual Bool_t selfNormalized() const {
    return kFALSE ;
  }

  // Support for extended maximum likelihood, switched off by default
  enum ExtendMode { CanNotBeExtended, CanBeExtended, MustBeExtended } ;
  /// Returns ability of PDF to provide extended likelihood terms. Possible
  /// answers are in the enumerator RooAbsPdf::ExtendMode.
  /// This default implementation always returns CanNotBeExtended.
  virtual ExtendMode extendMode() const { return CanNotBeExtended; }
  /// If true, PDF can provide extended likelihood term.
  inline Bool_t canBeExtended() const {
    return (extendMode() != CanNotBeExtended) ;
  }
  /// If true PDF must provide extended likelihood term.
  inline Bool_t mustBeExtended() const {
    return (extendMode() == MustBeExtended) ;
  }
  virtual Double_t expectedEvents(const RooArgSet* nset) const ;
  /// Return expected number of events to be used in calculation of extended likelihood.
  virtual Double_t expectedEvents(const RooArgSet& nset) const {
    return expectedEvents(&nset) ;
  }

  // Printing interface (human readable)
  virtual void printValue(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;

  static void verboseEval(Int_t stat) ;
  static int verboseEval() ;

  virtual Double_t extendedTerm(Double_t observedEvents, const RooArgSet* nset=0) const ;

  void setNormRange(const char* rangeName) ;
  const char* normRange() const {
    return _normRange.Length()>0 ? _normRange.Data() : 0 ;
  }
  void setNormRangeOverride(const char* rangeName) ;

  const RooAbsReal* getNormIntegral(const RooArgSet& nset) const { return getNormObj(0,&nset,0) ; }

protected:

public:
  virtual const RooAbsReal* getNormObj(const RooArgSet* set, const RooArgSet* iset, const TNamed* rangeName=0) const ;

  virtual RooAbsGenContext* binnedGenContext(const RooArgSet &vars, Bool_t verbose= kFALSE) const ;

  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0,
                                  const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;

  virtual RooAbsGenContext* autoGenContext(const RooArgSet &vars, const RooDataSet* prototype=0, const RooArgSet* auxProto=0,
                  Bool_t verbose=kFALSE, Bool_t autoBinned=kTRUE, const char* binnedTag="") const ;

private:

  RooDataSet *generate(RooAbsGenContext& context, const RooArgSet& whatVars, const RooDataSet* prototype,
             Double_t nEvents, Bool_t verbose, Bool_t randProtoOrder, Bool_t resampleProto, Bool_t skipInit=kFALSE,
             Bool_t extended=kFALSE) const ;

  // Implementation version
  virtual RooPlot* paramOn(RooPlot* frame, const RooArgSet& params, Bool_t showConstants=kFALSE,
                           const char *label= "", Int_t sigDigits = 2, Option_t *options = "NELU", Double_t xmin=0.65,
            Double_t xmax= 0.99,Double_t ymax=0.95, const RooCmdArg* formatCmd=0) ;

  void logBatchComputationErrors(RooSpan<const double>& outputs, std::size_t begin) const;

protected:
  virtual RooPlot *plotOn(RooPlot *frame, PlotOpt o) const;

  friend class RooEffGenContext ;
  friend class RooAddGenContext ;
  friend class RooProdGenContext ;
  friend class RooSimGenContext ;
  friend class RooSimSplitGenContext ;
  friend class RooConvGenContext ;
  friend class RooSimultaneous ;
  friend class RooAddGenContextOrig ;
  friend class RooProdPdf ;
  friend class RooMCStudy ;

  Int_t* randomizeProtoOrder(Int_t nProto,Int_t nGen,Bool_t resample=kFALSE) const ;

  friend class RooExtendPdf ;
  // This also forces the definition of a copy ctor in derived classes
  RooAbsPdf(const RooAbsPdf& other, const char* name = 0);

  friend class RooRealIntegral ;
  static Int_t _verboseEval ;

  virtual Bool_t syncNormalization(const RooArgSet* dset, Bool_t adjustProxies=kTRUE) const ;

  friend class RooAbsAnaConvPdf ;
  mutable Double_t _rawValue ;
  mutable RooAbsReal* _norm   ;      //! Normalization integral (owned by _normMgr)
  mutable RooArgSet* _normSet ;      //! Normalization set with for above integral

  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem(RooAbsReal& norm) : _norm(&norm) {} ;
    virtual ~CacheElem() ;
    virtual RooArgList containedArgs(Action) { return RooArgList(*_norm) ; }
    RooAbsReal* _norm ;
  } ;
  mutable RooObjCacheManager _normMgr ; // The cache manager

  friend class CacheElem ; // Cache needs to be able to clear _norm pointer

  virtual Bool_t redirectServersHook(const RooAbsCollection&, Bool_t, Bool_t, Bool_t) {
    // Hook function intercepting redirectServer calls. Discard current normalization
    // object if any server is redirected

    // Object is own by _normCacheManager that will delete object as soon as cache
    // is sterilized by server redirect
    _norm = 0 ;
    return kFALSE ;
  } ;


  mutable Int_t _errorCount ;        // Number of errors remaining to print
  mutable Int_t _traceCount ;        // Number of traces remaining to print
  mutable Int_t _negCount ;          // Number of negative probabilities remaining to print

  Bool_t _selectComp ;               // Component selection flag for RooAbsPdf::plotCompOn

  RooNumGenConfig* _specGeneratorConfig ; //! MC generator configuration specific for this object

  TString _normRange ; // Normalization range
  static TString _normRangeOverride ;

private:
  template<class Minimizer>
  int calculateAsymptoticCorrectedCovMatrix(Minimizer& minimizer, RooAbsData const& data);

  template<class Minimizer>
  int calculateSumW2CorrectedCovMatrix(Minimizer& minimizer, RooAbsReal const& nll) const;

  ClassDef(RooAbsPdf,4) // Abstract PDF with normalization support
};


#endif
