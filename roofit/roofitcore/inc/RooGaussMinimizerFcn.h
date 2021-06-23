/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       * 
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOGAUSSMINIMIZER

#ifndef ROO_GAUSS_MINIMIZER_FCN
#define ROO_GAUSS_MINIMIZER_FCN

#include "Math/IFunction.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"

#include "TMatrixDSym.h"

#include "RooAbsReal.h"
#include "RooArgList.h"

#include <iostream>
#include <fstream>

class RooGaussMinimizer;


class RooGaussMinimizerFcn : public ROOT::Math::IMultiGradFunction {

 public:

  RooGaussMinimizerFcn(RooAbsReal *funct, RooGaussMinimizer *context, 
	       bool verbose = false);
  RooGaussMinimizerFcn(const RooGaussMinimizerFcn& other);
  virtual ~RooGaussMinimizerFcn();

  virtual ROOT::Math::IMultiGradFunction* Clone() const;
  virtual unsigned int NDim() const { return _nDim; }

  RooArgList* GetFloatParamList() { return _floatParamList; }
  RooArgList* GetConstParamList() { return _constParamList; }
  RooArgList* GetInitFloatParamList() { return _initFloatParamList; }
  RooArgList* GetInitConstParamList() { return _initConstParamList; }

  void SetEvalErrorWall(Bool_t flag) { _doEvalErrorWall = flag ; }
  void SetPrintEvalErrors(Int_t numEvalErrors) { _printEvalErrors = numEvalErrors ; }
  Bool_t SetLogFile(const char* inLogfile);
  std::ofstream* GetLogFile() { return _logfile; }
  void SetVerbose(Bool_t flag=kTRUE) { _verbose = flag ; }

  Double_t& GetMaxFCN() { return _maxFCN; }
  Int_t GetNumInvalidNLL() { return _numBadNLL; }

  Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters, 
		     Bool_t optConst, Bool_t verbose);
  void BackProp(const ROOT::Fit::FitResult &results);  
  void ApplyCovarianceMatrix(TMatrixDSym& V); 

  Int_t evalCounter() const { return _evalCounter ; }
  void zeroEvalCount() { _evalCounter = 0 ; }


 private:
  
  Double_t GetPdfParamVal(Int_t index);
  Double_t GetPdfParamErr(Int_t index);
  void SetPdfParamErr(Int_t index, Double_t value);
  void ClearPdfParamAsymErr(Int_t index);
  void SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal);

  inline Bool_t SetPdfParamVal(const Int_t &index, const Double_t &value) const;


  virtual double DoEval(const double * x) const;  
  void updateFloatVec() ;

  virtual double DoDerivative(const double *x, unsigned int icoord) const;

private:

  mutable Int_t _evalCounter ;
  
  RooAbsReal *_funct;
  RooGaussMinimizer *_context;

  mutable double _maxFCN;
  mutable int _numBadNLL;
  mutable int _printEvalErrors;
  Bool_t _doEvalErrorWall;

  int _nDim;
  std::ofstream *_logfile;
  bool _verbose;

  RooArgList* _floatParamList;
  std::vector<RooAbsArg*> _floatParamVec ;
  RooArgList* _constParamList;
  RooArgList* _initFloatParamList;
  RooArgList* _initConstParamList;

};

#endif
#endif
