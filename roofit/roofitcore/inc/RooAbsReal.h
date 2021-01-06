/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsReal.h,v 1.75 2007/07/13 21:50:24 wouter Exp $
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
#ifndef ROO_ABS_REAL
#define ROO_ABS_REAL

#include "RooAbsArg.h"
#include "RooCmdArg.h"
#include "RooCurve.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooGlobalFunc.h"
#include "RooSpan.h"

#include <map>

class RooArgList ;
class RooDataSet ;
class RooPlot;
class RooRealVar;
class RooAbsFunc;
class RooAbsCategoryLValue ;
class RooCategory ;
class RooLinkedList ;
class RooNumIntConfig ;
class RooDataHist ;
class RooFunctor ;
class RooGenFunction ;
class RooMultiGenFunction ;
class RooFitResult ;
class RooAbsMoment ;
class RooDerivative ;
class RooVectorDataStore ;
namespace RooBatchCompute{
struct RunContext;
typedef std::unordered_map<const RooAbsReal*,RooSpan<const double>> DataMap;
}
namespace rbc = RooBatchCompute;
struct TreeReadBuffer; /// A space to attach TBranches

class TH1;
class TH1F;
class TH2F;
class TH3F;

#include <list>
#include <string>
#include <iostream>
#include <sstream>

class RooAbsReal : public RooAbsArg {
public:
  using value_type = double;

  // Constructors, assignment etc
  RooAbsReal() ;
  RooAbsReal(const char *name, const char *title, const char *unit= "") ;
  RooAbsReal(const char *name, const char *title, Double_t minVal, Double_t maxVal, 
	     const char *unit= "") ;
  RooAbsReal(const RooAbsReal& other, const char* name=0);
  RooAbsReal& operator=(const RooAbsReal& other);
  virtual ~RooAbsReal();




  //////////////////////////////////////////////////////////////////////////////////
  /// Evaluate object. Returns either cached value or triggers a recalculation.
  /// The recalculation happens by calling getValV(), which in the end calls the
  /// virtual evaluate() functions of the respective PDFs.
  /// \param[in] normalisationSet getValV() reacts differently depending on the value of the normalisation set.
  /// If the set is `nullptr`, an unnormalised value is returned.
  /// \note The normalisation is arbitrary, because it is up to the implementation
  /// of the PDF to e.g. leave out normalisation constants for speed reasons. The range
  /// of the variables is also ignored.
  ///
  /// To normalise the result properly, a RooArgSet has to be passed, which contains
  /// the variables to normalise over.
  /// These are integrated over their current ranges to compute the normalisation constant,
  /// and the unnormalised result is divided by this value.
  inline Double_t getVal(const RooArgSet* normalisationSet = nullptr) const {
    // Sometimes, the calling code uses an empty RooArgSet to request evaluation
    // without normalization set instead of following the `nullptr` convention.
    // To remove this ambiguity which might not always be correctly handled in
    // downstream code, we set `normalisationSet` to nullptr if it is pointing
    // to an empty set.
    if(normalisationSet && normalisationSet->empty()) {
      normalisationSet = nullptr;
    }
#ifdef ROOFIT_CHECK_CACHED_VALUES
    return _DEBUG_getVal(normalisationSet);
#else

#ifndef _WIN32
    return (_fast && !_inhibitDirty) ? _value : getValV(normalisationSet) ;
#else
    return (_fast && !inhibitDirty()) ? _value : getValV(normalisationSet) ;
#endif

#endif
  }

  /// Like getVal(const RooArgSet*), but always requires an argument for normalisation.
  inline  Double_t getVal(const RooArgSet& normalisationSet) const {
    // Sometimes, the calling code uses an empty RooArgSet to request evaluation
    // without normalization set instead of following the `nullptr` convention.
    // To remove this ambiguity which might not always be correctly handled in
    // downstream code, we set `normalisationSet` to nullptr if it is an empty set.
    return _fast ? _value : getValV(normalisationSet.empty() ? nullptr : &normalisationSet) ;
  }

  virtual Double_t getValV(const RooArgSet* normalisationSet = nullptr) const ;

  /// \deprecated getValBatch() has been removed in favour of the faster getValues(). If your code is affected
  /// by this change, please consult the release notes for ROOT 6.24 for guidance on how to make this transition.
  /// https://root.cern/doc/v624/release-notes.html
#ifndef R__MACOSX
  virtual RooSpan<const double> getValBatch(std::size_t /*begin*/, std::size_t /*maxSize*/, const RooArgSet* /*normSet*/ = nullptr) = delete;
#else
  //AppleClang in MacOS10.14 has a linker bug and fails to link programs that create objects of classes containing virtual deleted methods.
  //This can be safely deleted when MacOS10.14 is no longer supported by ROOT. See https://reviews.llvm.org/D37830
  virtual RooSpan<const double> getValBatch(std::size_t /*begin*/, std::size_t /*maxSize*/, const RooArgSet* /*normSet*/ = nullptr) final {
    throw std::logic_error("Deprecated getValBatch() has been removed in favour of the faster getValues(). If your code is affected by this change, please consult the release notes for ROOT 6.24 for guidance on how to make this transition. https://root.cern/doc/v624/release-notes.html");
  }
#endif
  /// by this change, please consult the release notes for ROOT 6.24 for guidance on how to make this transition.
  virtual RooSpan<const double> getValues(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet = nullptr) const;

  Double_t getPropagatedError(const RooFitResult &fr, const RooArgSet &nset = RooArgSet()) const;

  Bool_t operator==(Double_t value) const ;
  virtual Bool_t operator==(const RooAbsArg& other) const;
  virtual Bool_t isIdentical(const RooAbsArg& other, Bool_t assumeSameType=kFALSE) const;


  inline const Text_t *getUnit() const { 
    // Return string with unit description
    return _unit.Data(); 
  }
  inline void setUnit(const char *unit) { 
    // Set unit description to given string
    _unit= unit; 
  }
  TString getTitle(Bool_t appendUnit= kFALSE) const;

  // Lightweight interface adaptors (caller takes ownership)
  RooAbsFunc *bindVars(const RooArgSet &vars, const RooArgSet* nset=0, Bool_t clipInvalid=kFALSE) const;

  // Create a fundamental-type object that can hold our value.
  RooAbsArg *createFundamental(const char* newname=0) const;

  // Analytical integration support
  virtual Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const { 
    // Interface to force RooRealIntegral to offer given observable for internal integration
    // even if this is deemed unsafe. This default implementation returns always flase
    return kFALSE ; 
  }
  virtual void forceNumInt(Bool_t flag=kTRUE) { 
    // If flag is true, all advertised analytical integrals will be ignored
    // and all integrals are calculated numerically
    _forceNumInt = flag ; 
  }
  Bool_t getForceNumInt() const { return _forceNumInt ; }

  // Chi^2 fits to histograms
  virtual RooFitResult* chi2FitTo(RooDataHist& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),  
                              const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),  
                              const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;
  virtual RooFitResult* chi2FitTo(RooDataHist& data, const RooLinkedList& cmdList) ;

  virtual RooAbsReal* createChi2(RooDataHist& data, const RooLinkedList& cmdList) ;
  virtual RooAbsReal* createChi2(RooDataHist& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),  
				 const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),  
				 const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;

  // Chi^2 fits to X-Y datasets
  virtual RooFitResult* chi2FitTo(RooDataSet& xydata, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),  
                              const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),  
                              const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;
  virtual RooFitResult* chi2FitTo(RooDataSet& xydata, const RooLinkedList& cmdList) ;

  virtual RooAbsReal* createChi2(RooDataSet& data, const RooLinkedList& cmdList) ;
  virtual RooAbsReal* createChi2(RooDataSet& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),  
				   const RooCmdArg& arg3=RooCmdArg::none(),  const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg& arg5=RooCmdArg::none(),  
				   const RooCmdArg& arg6=RooCmdArg::none(),  const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;


  virtual RooAbsReal* createProfile(const RooArgSet& paramsOfInterest) ;


  RooAbsReal* createIntegral(const RooArgSet& iset, const RooCmdArg& arg1, const RooCmdArg& arg2=RooCmdArg::none(),
                             const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(), 
			     const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(), 
			     const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) const ;

  /// Create integral over observables in iset in range named rangeName.
  RooAbsReal* createIntegral(const RooArgSet& iset, const char* rangeName) const { 
    return createIntegral(iset,0,0,rangeName) ; 
  }
  /// Create integral over observables in iset in range named rangeName with integrand normalized over observables in nset
  RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet& nset, const char* rangeName=0) const { 
    return createIntegral(iset,&nset,0,rangeName) ; 
  }
  /// Create integral over observables in iset in range named rangeName with integrand normalized over observables in nset while
  /// using specified configuration for any numeric integration.
  RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet& nset, const RooNumIntConfig& cfg, const char* rangeName=0) const {
    return createIntegral(iset,&nset,&cfg,rangeName) ; 
  }
  /// Create integral over observables in iset in range named rangeName using specified configuration for any numeric integration.
  RooAbsReal* createIntegral(const RooArgSet& iset, const RooNumIntConfig& cfg, const char* rangeName=0) const { 
    return createIntegral(iset,0,&cfg,rangeName) ; 
  }
  virtual RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet* nset=0, const RooNumIntConfig* cfg=0, const char* rangeName=0) const ;  


  void setParameterizeIntegral(const RooArgSet& paramVars) ;

  // Create running integrals
  RooAbsReal* createRunningIntegral(const RooArgSet& iset, const RooArgSet& nset=RooArgSet()) ;
  RooAbsReal* createRunningIntegral(const RooArgSet& iset, const RooCmdArg& arg1, const RooCmdArg& arg2=RooCmdArg::none(),
			const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(), 
			const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(), 
			const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;
  RooAbsReal* createIntRI(const RooArgSet& iset, const RooArgSet& nset=RooArgSet()) ;
  RooAbsReal* createScanRI(const RooArgSet& iset, const RooArgSet& nset, Int_t numScanBins, Int_t intOrder) ;

  
  // Optimized accept/reject generator support
  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;
  virtual Int_t minTrialSamples(const RooArgSet& /*arGenObs*/) const { return 0 ; }


  // Plotting options
  void setPlotLabel(const char *label);
  const char *getPlotLabel() const;

  virtual Double_t defaultErrorLevel() const { 
    // Return default level for MINUIT error analysis
    return 1.0 ; 
  }

  const RooNumIntConfig* getIntegratorConfig() const ;
  RooNumIntConfig* getIntegratorConfig() ;
  static RooNumIntConfig* defaultIntegratorConfig()  ;
  RooNumIntConfig* specialIntegratorConfig() const ;
  RooNumIntConfig* specialIntegratorConfig(Bool_t createOnTheFly) ;
  void setIntegratorConfig() ;
  void setIntegratorConfig(const RooNumIntConfig& config) ;

  virtual void fixAddCoefNormalization(const RooArgSet& addNormSet=RooArgSet(),Bool_t force=kTRUE) ;
  virtual void fixAddCoefRange(const char* rangeName=0,Bool_t force=kTRUE) ;

  virtual void preferredObservableScanOrder(const RooArgSet& obs, RooArgSet& orderedObs) const ;

  // User entry point for plotting
  virtual RooPlot* plotOn(RooPlot* frame, 
			  const RooCmdArg& arg1=RooCmdArg(), const RooCmdArg& arg2=RooCmdArg(),
			  const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
			  const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),
			  const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg(),
			  const RooCmdArg& arg9=RooCmdArg(), const RooCmdArg& arg10=RooCmdArg()
              ) const ;


  enum ScaleType { Raw, Relative, NumEvent, RelativeExpected } ;

  // Forwarder function for backward compatibility
  virtual RooPlot *plotSliceOn(RooPlot *frame, const RooArgSet& sliceSet, Option_t* drawOptions="L", 
			       Double_t scaleFactor=1.0, ScaleType stype=Relative, const RooAbsData* projData=0) const;

  // Fill an existing histogram
  TH1 *fillHistogram(TH1 *hist, const RooArgList &plotVars,
		     Double_t scaleFactor= 1, const RooArgSet *projectedVars= 0, Bool_t scaling=kTRUE,
		     const RooArgSet* condObs=0, Bool_t setError=kTRUE) const;

  // Create 1,2, and 3D histograms from and fill it
  TH1 *createHistogram(const char* varNameList, Int_t xbins=0, Int_t ybins=0, Int_t zbins=0) const ;
  TH1* createHistogram(const char *name, const RooAbsRealLValue& xvar, RooLinkedList& argList) const ;
  TH1 *createHistogram(const char *name, const RooAbsRealLValue& xvar,
                       const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(), 
                       const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(), 
                       const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(), 
                       const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) const ;

  // Fill a RooDataHist
  RooDataHist* fillDataHist(RooDataHist *hist, const RooArgSet* nset, Double_t scaleFactor,
			    Bool_t correctForBinVolume=kFALSE, Bool_t showProgress=kFALSE) const ;

  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

  // Printing interface (human readable)
  virtual void printValue(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;

  inline void setCachedValue(double value, bool notifyClients = true) final;

  // Evaluation error logging 
  class EvalError {
  public:
    EvalError() { }
    EvalError(const EvalError& other) : _msg(other._msg), _srvval(other._srvval) { }
    void setMessage(const char* tmp) { std::string s(tmp); s.swap(_msg); }
    void setServerValues(const char* tmp) { std::string s(tmp); s.swap(_srvval); }
    std::string _msg;
    std::string _srvval;
  } ;

  enum ErrorLoggingMode { PrintErrors, CollectErrors, CountErrors, Ignore } ;
  static ErrorLoggingMode evalErrorLoggingMode() ;
  static void setEvalErrorLoggingMode(ErrorLoggingMode m) ;
  void logEvalError(const char* message, const char* serverValueString=0) const ;
  static void logEvalError(const RooAbsReal* originator, const char* origName, const char* message, const char* serverValueString=0) ;
  static void printEvalErrors(std::ostream&os=std::cout, Int_t maxPerNode=10000000) ;
  static Int_t numEvalErrors() ;
  static Int_t numEvalErrorItems() ;

   
  typedef std::map<const RooAbsArg*,std::pair<std::string,std::list<EvalError> > >::const_iterator EvalErrorIter ; 
  static EvalErrorIter evalErrorIter() ;

  static void clearEvalErrorLog() ;
  
  /// Tests if the distribution is binned. Unless overridden by derived classes, this always returns false.
  virtual Bool_t isBinnedDistribution(const RooArgSet& /*obs*/) const { return kFALSE ; }
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const;

  RooGenFunction* iGenFunction(RooRealVar& x, const RooArgSet& nset=RooArgSet()) ;
  RooMultiGenFunction* iGenFunction(const RooArgSet& observables, const RooArgSet& nset=RooArgSet()) ;

  RooFunctor* functor(const RooArgList& obs, const RooArgList& pars=RooArgList(), const RooArgSet& nset=RooArgSet()) const ;
  TF1* asTF(const RooArgList& obs, const RooArgList& pars=RooArgList(), const RooArgSet& nset=RooArgSet()) const ;

  RooDerivative* derivative(RooRealVar& obs, Int_t order=1, Double_t eps=0.001) ;
  RooDerivative* derivative(RooRealVar& obs, const RooArgSet& normSet, Int_t order, Double_t eps=0.001) ; 

  RooAbsMoment* moment(RooRealVar& obs, Int_t order, Bool_t central, Bool_t takeRoot) ;
  RooAbsMoment* moment(RooRealVar& obs, const RooArgSet& normObs, Int_t order, Bool_t central, Bool_t takeRoot, Bool_t intNormObs) ;

  RooAbsMoment* mean(RooRealVar& obs) { return moment(obs,1,kFALSE,kFALSE) ; }
  RooAbsMoment* mean(RooRealVar& obs, const RooArgSet& nset) { return moment(obs,nset,1,kFALSE,kFALSE,kTRUE) ; }
  RooAbsMoment* sigma(RooRealVar& obs) { return moment(obs,2,kTRUE,kTRUE) ; }
  RooAbsMoment* sigma(RooRealVar& obs, const RooArgSet& nset) { return moment(obs,nset,2,kTRUE,kTRUE,kTRUE) ; }

  Double_t findRoot(RooRealVar& x, Double_t xmin, Double_t xmax, Double_t yval) ;


  virtual Bool_t setData(RooAbsData& /*data*/, Bool_t /*cloneData*/=kTRUE) { return kTRUE ; }

  virtual void enableOffsetting(Bool_t) {} ;
  virtual Bool_t isOffsetting() const { return kFALSE ; }
  virtual Double_t offset() const { return 0 ; }
  
  static void setHideOffset(Bool_t flag);
  static Bool_t hideOffset() ;

protected:
  // Hook for objects with normalization-dependent parameters interperetation
  virtual void selectNormalization(const RooArgSet* depSet=0, Bool_t force=kFALSE) ;
  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;

  // Helper functions for plotting
  Bool_t plotSanityChecks(RooPlot* frame) const ;
  void makeProjectionSet(const RooAbsArg* plotVar, const RooArgSet* allVars, 
			 RooArgSet& projectedVars, Bool_t silent) const ;

  TString integralNameSuffix(const RooArgSet& iset, const RooArgSet* nset=0, const char* rangeName=0, Bool_t omitEmpty=kFALSE) const ;


  Bool_t isSelectedComp() const ;

  
 public:
  const RooAbsReal* createPlotProjection(const RooArgSet& depVars, const RooArgSet& projVars, RooArgSet*& cloneSet) const ;
  const RooAbsReal *createPlotProjection(const RooArgSet &dependentVars, const RooArgSet *projectedVars,
				         RooArgSet *&cloneSet, const char* rangeName=0, const RooArgSet* condObs=0) const;
  virtual void computeBatch(double*, size_t, rbc::DataMap&) const {throw std::runtime_error("computeBatch not implemented in this derived class"); }

 protected:

  RooFitResult* chi2FitDriver(RooAbsReal& fcn, RooLinkedList& cmdList) ;

  void plotOnCompSelect(RooArgSet* selNodes) const ;
  RooPlot* plotOnWithErrorBand(RooPlot* frame,const RooFitResult& fr, Double_t Z, const RooArgSet* params, const RooLinkedList& argList, Bool_t method1) const ;

  // Support interface for subclasses to advertise their analytic integration
  // and generator capabilities in their analyticalIntegral() and generateEvent()
  // implementations.
  Bool_t matchArgs(const RooArgSet& allDeps, RooArgSet& numDeps, 
		   const RooArgProxy& a) const ;
  Bool_t matchArgs(const RooArgSet& allDeps, RooArgSet& numDeps, 
		   const RooArgProxy& a, const RooArgProxy& b) const ;
  Bool_t matchArgs(const RooArgSet& allDeps, RooArgSet& numDeps, 
		   const RooArgProxy& a, const RooArgProxy& b, const RooArgProxy& c) const ;
  Bool_t matchArgs(const RooArgSet& allDeps, RooArgSet& numDeps, 
		   const RooArgProxy& a, const RooArgProxy& b, 		   
		   const RooArgProxy& c, const RooArgProxy& d) const ;

  Bool_t matchArgs(const RooArgSet& allDeps, RooArgSet& numDeps, 
		   const RooArgSet& set) const ;


  RooAbsReal* createIntObj(const RooArgSet& iset, const RooArgSet* nset, const RooNumIntConfig* cfg, const char* rangeName) const ;
  void findInnerMostIntegration(const RooArgSet& allObs, RooArgSet& innerObs, const char* rangeName) const ;


  // Internal consistency checking (needed by RooDataSet)
  /// Check if current value is valid.
  virtual bool isValid() const { return isValidReal(_value); }
  /// Interface function to check if given value is a valid value for this object. Returns true unless overridden.
  virtual bool isValidReal(double /*value*/, bool printError = false) const { (void)printError; return true; }


  // Function evaluation and error tracing
  Double_t traceEval(const RooArgSet* set) const ;

  /// Evaluate this PDF / function / constant. Needs to be overridden by all derived classes.
  virtual Double_t evaluate() const = 0;

  /// \deprecated evaluateBatch() has been removed in favour of the faster evaluateSpan(). If your code is affected
  /// by this change, please consult the release notes for ROOT 6.24 for guidance on how to make this transition.
  /// https://root.cern/doc/v624/release-notes.html
#ifndef R__MACOSX
  virtual RooSpan<double> evaluateBatch(std::size_t /*begin*/, std::size_t /*maxSize*/) = delete;
#else
  //AppleClang in MacOS10.14 has a linker bug and fails to link programs that create objects of classes containing virtual deleted methods.
  //This can be safely deleted when MacOS10.14 is no longer supported by ROOT. See https://reviews.llvm.org/D37830
  virtual RooSpan<double> evaluateBatch(std::size_t /*begin*/, std::size_t /*maxSize*/) final {
    throw std::logic_error("Deprecated evaluatedBatch() has been removed in favour of the faster evaluateSpan(). If your code is affected by this change, please consult the release notes for ROOT 6.24 for guidance on how to make this transition. https://root.cern/doc/v624/release-notes.html");
  }
#endif

  virtual RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const;

  //---------- Interface to access batch data ---------------------------
  //
  friend class BatchInterfaceAccessor;
  
 private:
  void checkBatchComputation(const RooBatchCompute::RunContext& evalData, std::size_t evtNo, const RooArgSet* normSet = nullptr, double relAccuracy = 1.E-13) const;

  /// Debug version of getVal(), which is slow and does error checking.
  Double_t _DEBUG_getVal(const RooArgSet* normalisationSet) const;

  //--------------------------------------------------------------------

 protected:
  // Hooks for RooDataSet interface
  friend class RooRealIntegral ;
  friend class RooVectorDataStore ;
  virtual void syncCache(const RooArgSet* set=0) { getVal(set) ; }
  virtual void copyCache(const RooAbsArg* source, Bool_t valueOnly=kFALSE, Bool_t setValDirty=kTRUE) ;
  virtual void attachToTree(TTree& t, Int_t bufSize=32000) ;
  virtual void attachToVStore(RooVectorDataStore& vstore) ;
  virtual void setTreeBranchStatus(TTree& t, Bool_t active) ;
  virtual void fillTreeBranch(TTree& t) ;

  friend class RooRealBinding ;
  Double_t _plotMin ;       // Minimum of plot range
  Double_t _plotMax ;       // Maximum of plot range
  Int_t    _plotBins ;      // Number of plot bins
  mutable Double_t _value ; // Cache for current value of object
  TString  _unit ;          // Unit for objects value
  TString  _label ;         // Plot label for objects value
  Bool_t   _forceNumInt ;   // Force numerical integration if flag set

  friend class RooAbsPdf ;
  friend class RooAbsAnaConvPdf ;

  RooNumIntConfig* _specIntegratorConfig ; // Numeric integrator configuration specific for this object

  friend class RooDataProjBinding ;
  friend class RooAbsOptGoodnessOfFit ;
  
  struct PlotOpt {
   PlotOpt() : drawOptions("L"), scaleFactor(1.0), stype(Relative), projData(0), binProjData(kFALSE), projSet(0), precision(1e-3), 
               shiftToZero(kFALSE),projDataSet(0),normRangeName(0),rangeLo(0),rangeHi(0),postRangeFracScale(kFALSE),wmode(RooCurve::Extended),
               projectionRangeName(0),curveInvisible(kFALSE), curveName(0),addToCurveName(0),addToWgtSelf(1.),addToWgtOther(1.),
               numCPU(1),interleave(RooFit::Interleave),curveNameSuffix(""), numee(10), eeval(0), doeeval(kFALSE), progress(kFALSE), errorFR(0) {} ;
   Option_t* drawOptions ;
   Double_t scaleFactor ;	 
   ScaleType stype ;
   const RooAbsData* projData ;
   Bool_t binProjData ;
   const RooArgSet* projSet ;
   Double_t precision ;
   Bool_t shiftToZero ;
   const RooArgSet* projDataSet ;
   const char* normRangeName ;
   Double_t rangeLo ;
   Double_t rangeHi ;
   Bool_t postRangeFracScale ;
   RooCurve::WingMode wmode ;
   const char* projectionRangeName ;
   Bool_t curveInvisible ;
   const char* curveName ;
   const char* addToCurveName ;
   Double_t addToWgtSelf ;
   Double_t addToWgtOther ;
   Int_t    numCPU ;
   RooFit::MPSplit interleave ;
   const char* curveNameSuffix ; 
   Int_t    numee ;
   Double_t eeval ;
   Bool_t   doeeval ;
   Bool_t progress ;
   const RooFitResult* errorFR ;
  } ;

  // Plot implementation functions
  virtual RooPlot *plotOn(RooPlot* frame, PlotOpt o) const;

public:
  // PlotOn with command list
  virtual RooPlot* plotOn(RooPlot* frame, RooLinkedList& cmdList) const ;

 protected:
  virtual RooPlot *plotAsymOn(RooPlot *frame, const RooAbsCategoryLValue& asymCat, PlotOpt o) const;


private:

  static ErrorLoggingMode _evalErrorMode ;
  static std::map<const RooAbsArg*,std::pair<std::string,std::list<EvalError> > > _evalErrorList ;
  static Int_t _evalErrorCount ;

  Bool_t matchArgsByName(const RooArgSet &allArgs, RooArgSet &matchedArgs, const TList &nameList) const;

  std::unique_ptr<TreeReadBuffer> _treeReadBuffer; //! A buffer for reading values from trees

protected:


  friend class RooRealSumPdf ;
  friend class RooRealSumFunc;
  friend class RooAddPdf ;
  friend class RooAddModel ;
  void selectComp(Bool_t flag) { 
    // If flag is true, only selected component will be included in evaluates of RooAddPdf components
    _selectComp = flag ; 
  }
  static void globalSelectComp(Bool_t flag) ;
  Bool_t _selectComp ;               //! Component selection flag for RooAbsPdf::plotCompOn
  static Bool_t _globalSelectComp ;  // Global activation switch for component selection
  // This struct can be used to flip the global switch to select components.
  // Doing this with RAII prevents forgetting to reset the state.
  struct GlobalSelectComponentRAII {
      GlobalSelectComponentRAII(bool state) :
      _oldState{_globalSelectComp} {
        if (state != RooAbsReal::_globalSelectComp)
          RooAbsReal::_globalSelectComp = state;
      }

      ~GlobalSelectComponentRAII() {
        if (RooAbsReal::_globalSelectComp != _oldState)
          RooAbsReal::_globalSelectComp = _oldState;
      }

      bool _oldState;
  };


  mutable RooArgSet* _lastNSet ; //!
  static Bool_t _hideOffset ; // Offset hiding flag

  ClassDef(RooAbsReal,2) // Abstract real-valued variable
};


/// Helper class to access a batch-related part of RooAbsReal's interface, which should not leak to the outside world.
class BatchInterfaceAccessor {
  public:
    static void checkBatchComputation(const RooAbsReal& theReal, const RooBatchCompute::RunContext& evalData, std::size_t evtNo,
        const RooArgSet* normSet = nullptr, double relAccuracy = 1.E-13) {
      theReal.checkBatchComputation(evalData, evtNo, normSet, relAccuracy);
    }
};


////////////////////////////////////////////////////////////////////////////////
/// Overwrite the value stored in this object's cache.
/// This can be used to fake a computation that resulted in `value`.
/// \param[in] value Value to write.
/// \param[in] setValDirty If true, notify users of this object that its value changed.
/// This is the default.
void RooAbsReal::setCachedValue(double value, bool notifyClients) {
  _value = value;

  if (notifyClients) {
    setValueDirty();
    _valueDirty = false;
  }
}


#endif
