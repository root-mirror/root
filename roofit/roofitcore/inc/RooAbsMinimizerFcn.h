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

#ifndef ROO_ABS_MINIMIZER_FCN
#define ROO_ABS_MINIMIZER_FCN

#include "Math/IFunction.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"

#include "TMatrixDSym.h"

#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooRealVar.h"

#include <iostream>
#include <fstream>
#include <string>

// forward declaration
class RooMinimizer;

class RooAbsMinimizerFcn {

public:
   RooAbsMinimizerFcn(RooArgList paramList, RooMinimizer *context, bool verbose = false);
   RooAbsMinimizerFcn(const RooAbsMinimizerFcn &other);
   virtual ~RooAbsMinimizerFcn();

   // inform Minuit through its parameter_settings vector of RooFit parameter properties
   Bool_t synchronize_parameter_settings(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose);
   // same, but can be overridden to e.g. also include gradient strategy synchronization in subclasses:
   virtual Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose);

   // used for export to RooFitResult from Minimizer:
   RooArgList *GetFloatParamList();
   RooArgList *GetConstParamList();
   RooArgList *GetInitFloatParamList();
   RooArgList *GetInitConstParamList();
   Int_t GetNumInvalidNLL() const;

   // need access from Minimizer:
   void SetEvalErrorWall(Bool_t flag);
   /// Try to recover from invalid function values. When invalid function values are encountered,
   /// a penalty term is returned to the minimiser to make it back off. This sets the strength of this penalty.
   /// \note A strength of zero is equivalent to a constant penalty (= the gradient vanishes, ROOT < 6.24).
   /// Positive values lead to a gradient pointing away from the undefined regions. Use ~10 to force the minimiser
   /// away from invalid function values.
   void SetRecoverFromNaNStrength(double strength) { _recoverFromNaNStrength = strength; }
   void SetPrintEvalErrors(Int_t numEvalErrors);
   Double_t &GetMaxFCN();
   Int_t evalCounter() const;
   void zeroEvalCount();
   /// Return a possible offset that's applied to the function to separate invalid function values from valid ones.
   double getOffset() const { return _funcOffset; }
   void SetVerbose(Bool_t flag = kTRUE);

   // put Minuit results back into RooFit objects:
   void BackProp(const ROOT::Fit::FitResult &results);

   // used in several Minimizer functions:
   virtual std::string getFunctionName() const = 0;
   virtual std::string getFunctionTitle() const = 0;

   // set different external covariance matrix
   void ApplyCovarianceMatrix(TMatrixDSym &V);

   Bool_t SetLogFile(const char *inLogfile);
   std::ofstream *GetLogFile() { return _logfile; }

   unsigned int get_nDim() const { return _nDim; }

   virtual void setOptimizeConst(Int_t flag) = 0;

   bool getOptConst();
   std::vector<double> get_parameter_values() const;

   Bool_t SetPdfParamVal(int index, double value) const;

protected:
   virtual void optimizeConstantTerms(bool constStatChange, bool constValChange) = 0;

   // used in BackProp (Minuit results -> RooFit) and ApplyCovarianceMatrix
   void SetPdfParamErr(Int_t index, Double_t value);
   void ClearPdfParamAsymErr(Int_t index);
   void SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal);

   void printEvalErrors() const;

   // members
   RooMinimizer *_context;

   // the following four are mutable because DoEval is const (in child classes)
   // Reset the *largest* negative log-likelihood value we have seen so far:
   mutable double _maxFCN = -std::numeric_limits<double>::infinity();
   mutable double _funcOffset{0.};
   double _recoverFromNaNStrength{10.};
   mutable int _numBadNLL = 0;
   mutable int _printEvalErrors = 10;
   mutable int _evalCounter{0};

   unsigned int _nDim = 0;

   Bool_t _optConst = kFALSE;

   RooArgList *_floatParamList;
   RooArgList *_constParamList;
   RooArgList *_initFloatParamList;
   RooArgList *_initConstParamList;

   std::ofstream *_logfile = nullptr;
   bool _doEvalErrorWall{true};
   bool _verbose;
};

#endif
#endif
