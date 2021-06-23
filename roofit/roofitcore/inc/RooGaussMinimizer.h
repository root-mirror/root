/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_GAUSS_MINIMIZER
#define ROO_GAUSS_MINIMIZER

#include "TObject.h"
#include "TStopwatch.h"
#include <fstream>
#include "TMatrixDSymfwd.h"


#include "Fit/Fitter.h"
#include "RooMinimizerFcn.h"
#include "RooGaussMinimizerFcn.h"

class RooAbsReal ;
class RooFitResult ;
class RooArgList ;
class RooRealVar ;
class RooArgSet ;
class TH2F ;
class RooPlot ;

class RooGaussMinimizer : public TObject {
public:

  RooGaussMinimizer(RooAbsReal& function) ;
  virtual ~RooGaussMinimizer() ;

  enum Strategy { Speed=0, Balance=1, Robustness=2 } ;
  enum PrintLevel { None=-1, Reduced=0, Normal=1, ExtraForProblem=2, Maximum=3 } ;
  void optimizeConst(Int_t flag) ;

  RooFitResult* fit(const char* options) ;

  Int_t migrad() ;
  Int_t hesse() ;
  Int_t minos() ;
  Int_t minos(const RooArgSet& minosParamList) ;

  Int_t minimize(const char* type, const char* alg=0) ;

  static void cleanup() ;

  void saveStatus(const char* label, Int_t status) { _statusHistory.push_back(std::pair<std::string,int>(label,status)) ; }

  Int_t evalCounter() const { return fitterFcn()->evalCounter() ; }
  void zeroEvalCount() { fitterFcn()->zeroEvalCount() ; }

  ROOT::Fit::Fitter* fitter() ;
  const ROOT::Fit::Fitter* fitter() const ;
  
protected:

  friend class RooAbsPdf ;

  inline Int_t getNPar() const { return fitterFcn()->NDim() ; }
  inline std::ofstream* logfile() { return fitterFcn()->GetLogFile(); }
  inline Double_t& maxFCN() { return fitterFcn()->GetMaxFCN() ; }
  
  const RooGaussMinimizerFcn* fitterFcn() const {  return ( fitter()->GetFCN() ? (dynamic_cast<RooGaussMinimizerFcn*>(fitter()->GetFCN())) : _fcn ) ; }
  RooGaussMinimizerFcn* fitterFcn() { return ( fitter()->GetFCN() ? (dynamic_cast<RooGaussMinimizerFcn*>(fitter()->GetFCN())) : _fcn ) ; }

private:

  Int_t       _printLevel ;
  Int_t       _status ;
  Bool_t      _optConst ;
  Bool_t      _profile ;
  RooAbsReal* _func ;

  Bool_t      _verbose ;
  TStopwatch  _timer ;
  TStopwatch  _cumulTimer ;
  Bool_t      _profileStart ;

  TMatrixDSym* _extV ;

  RooGaussMinimizerFcn *_fcn;
  std::string _minimizerType;

  static ROOT::Fit::Fitter *_theFitter ;

  std::vector<std::pair<std::string,int> > _statusHistory ;

  RooGaussMinimizer(const RooGaussMinimizer&) ;
	
} ;

#endif
