/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

#ifndef ROO_MINIMIZER_FCN
#define ROO_MINIMIZER_FCN

#include "Math/IFunction.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"

#include "RooAbsReal.h"
#include "RooArgList.h"

#include <fstream>
#include <vector>

class RooMinimizer;
template<typename T> class TMatrixTSym;
using TMatrixDSym = TMatrixTSym<double>;

class RooMinimizerFcn : public ROOT::Math::IBaseFunctionMultiDim {

 public:

  RooMinimizerFcn(RooAbsReal *funct, RooMinimizer *context,
          bool verbose = false);
  RooMinimizerFcn(const RooMinimizerFcn& other);
  virtual ~RooMinimizerFcn();

  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const;
  virtual unsigned int NDim() const { return _nDim; }

  RooArgList* GetFloatParamList() { return _floatParamList; }
  RooArgList* GetConstParamList() { return _constParamList; }
  RooArgList* GetInitFloatParamList() { return _initFloatParamList; }
  RooArgList* GetInitConstParamList() { return _initConstParamList; }

  void SetEvalErrorWall(Bool_t flag) { _doEvalErrorWall = flag ; }
  /// Try to recover from invalid function values. When invalid function values are encountered,
  /// a penalty term is returned to the minimiser to make it back off. This sets the strength of this penalty.
  /// \note A strength of zero is equivalent to a constant penalty (= the gradient vanishes, ROOT < 6.24).
  /// Positive values lead to a gradient pointing away from the undefined regions. Use ~10 to force the minimiser
  /// away from invalid function values.
  void SetRecoverFromNaNStrength(double strength) { _recoverFromNaNStrength = strength; }
  void SetPrintEvalErrors(Int_t numEvalErrors) { _printEvalErrors = numEvalErrors ; }
  Bool_t SetLogFile(const char* inLogfile);
  std::ofstream* GetLogFile() { return _logfile; }
  void SetVerbose(Bool_t flag=kTRUE) { _verbose = flag ; }

  Double_t& GetMaxFCN() { return _maxFCN; }
  Int_t GetNumInvalidNLL() const { return _numBadNLL; }

  Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters,
           Bool_t optConst, Bool_t verbose);
  void BackProp(const ROOT::Fit::FitResult &results);
  void ApplyCovarianceMatrix(TMatrixDSym& V);

  Int_t evalCounter() const { return _evalCounter ; }
  void zeroEvalCount() { _evalCounter = 0 ; }
  /// Return a possible offset that's applied to the function to separate invalid function values from valid ones.
  double getOffset() const { return _funcOffset; }

 private:
  void SetPdfParamErr(Int_t index, Double_t value);
  void ClearPdfParamAsymErr(Int_t index);
  void SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal);

  Bool_t SetPdfParamVal(int index, double value) const;
  void printEvalErrors() const;

  virtual double DoEval(const double * x) const;


  RooAbsReal *_funct;
  const RooMinimizer *_context;

  mutable double _maxFCN;
  mutable double _funcOffset{0.};
  double _recoverFromNaNStrength{10.};
  mutable int _numBadNLL;
  mutable int _printEvalErrors;
  mutable int _evalCounter{0};
  int _nDim;

  RooArgList* _floatParamList;
  RooArgList* _constParamList;
  RooArgList* _initFloatParamList;
  RooArgList* _initConstParamList;

  std::ofstream *_logfile;
  bool _doEvalErrorWall{true};
  bool _verbose;

};

#endif
#endif
