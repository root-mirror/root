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

// Global helper functions

#include "RooFit.h"

#include "RooGlobalFunc.h"
#include "RooCategory.h"
#include "RooRealConstant.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooNumIntConfig.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "TH1.h"

using namespace std;

namespace RooFit {

  // RooAbsReal::plotOn arguments
  RooCmdArg DrawOption(const char* opt)            {
     return RooCmdArg("DrawOption", 0, 0, 0, 0, opt, nullptr, nullptr, nullptr);
  }
  RooCmdArg Slice(const RooArgSet& sliceSet)       {
     return RooCmdArg("SliceVars", 0, 0, 0, 0, nullptr, nullptr, &sliceSet, nullptr);
  }
  RooCmdArg Slice(RooCategory& cat, const char* label) {
     return RooCmdArg("SliceCat", 0, 0, 0, 0, label, nullptr, &cat, nullptr);
  }

  RooCmdArg Project(const RooArgSet& projSet)      {
     return RooCmdArg("Project", 0, 0, 0, 0, nullptr, nullptr, &projSet, nullptr);
  }
  RooCmdArg ProjWData(const RooArgSet& projSet, 
                      const RooAbsData& projData,
                      Bool_t binData)              {
     return RooCmdArg("ProjData", binData, 0, 0, 0, nullptr, nullptr, &projSet, &projData);
  }
  RooCmdArg ProjWData(const RooAbsData& projData,
                      Bool_t binData)              {
     return RooCmdArg("ProjData", binData, 0, 0, 0, nullptr, nullptr, nullptr, &projData);
  }
  RooCmdArg Asymmetry(const RooCategory& cat)      {
     return RooCmdArg("Asymmetry", 0, 0, 0, 0, nullptr, nullptr, &cat, nullptr);
  }
  RooCmdArg Precision(Double_t prec)               {
     return RooCmdArg("Precision", 0, 0, prec, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ShiftToZero()                          {
     return RooCmdArg("ShiftToZero", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Normalization(Double_t scaleFactor)    {
     return RooCmdArg("Normalization", RooAbsReal::Relative, 0, scaleFactor, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Range(const char* rangeName, Bool_t adjustNorm)   {
     return RooCmdArg("RangeWithName", adjustNorm, 0, 0, 0, rangeName, nullptr, nullptr, nullptr);
  }
  RooCmdArg Range(Double_t lo, Double_t hi, Bool_t adjustNorm){
     return RooCmdArg("Range", adjustNorm, 0, lo, hi, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg NormRange(const char* rangeNameList)   {
     return RooCmdArg("NormRange", 0, 0, 0, 0, rangeNameList, nullptr, nullptr, nullptr);
  }
  RooCmdArg VLines()                               {
     return RooCmdArg("VLines", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg LineColor(Color_t color)               {
     return RooCmdArg("LineColor", color, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg LineStyle(Style_t style)               {
     return RooCmdArg("LineStyle", style, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg LineWidth(Width_t width)               {
     return RooCmdArg("LineWidth", width, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg FillColor(Color_t color)               {
     return RooCmdArg("FillColor", color, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg FillStyle(Style_t style)               {
     return RooCmdArg("FillStyle", style, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ProjectionRange(const char* rangeName) {
     return RooCmdArg("ProjectionRange", 0, 0, 0, 0, rangeName, nullptr, nullptr, nullptr);
  }
  RooCmdArg Name(const char* name)                 {
     return RooCmdArg("Name", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg Invisible()                            {
     return RooCmdArg("Invisible", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg AddTo(const char* name, double wgtSel, double wgtOther) {
     return RooCmdArg("AddTo", 0, 0, wgtSel, wgtOther, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg EvalErrorValue(Double_t val)           {
     return RooCmdArg("EvalErrorValue", 1, 0, val, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg MoveToBack()                           {
     return RooCmdArg("MoveToBack", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg VisualizeError(const RooFitResult& fitres, Double_t Z, Bool_t EVmethod)  {
     return RooCmdArg("VisualizeError", EVmethod, 0, Z, 0, nullptr, nullptr, &fitres, nullptr);
  }
  RooCmdArg VisualizeError(const RooFitResult& fitres, const RooArgSet& param, Double_t Z, Bool_t EVmethod) 
                                                                  {
                                                                     return RooCmdArg("VisualizeError", EVmethod, 0, Z,
                                                                                      0, nullptr, nullptr, &fitres,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      &param);
                                                                  }
                                                                  RooCmdArg VisualizeError(const RooDataSet &paramData,
                                                                                           Double_t Z)
                                                                  {
                                                                     return RooCmdArg("VisualizeErrorData", 0, 0, Z, 0,
                                                                                      nullptr, nullptr, &paramData,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg ShowProgress()
                                                                  {
                                                                     return RooCmdArg("ShowProgress", 1, 0, 0, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }

                                                                  // RooAbsPdf::plotOn arguments
                                                                  RooCmdArg Components(const RooArgSet &compSet)
                                                                  {
                                                                     return RooCmdArg("SelectCompSet", 0, 0, 0, 0,
                                                                                      nullptr, nullptr, &compSet,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg Components(const char *compSpec)
                                                                  {
                                                                     return RooCmdArg("SelectCompSpec", 0, 0, 0, 0,
                                                                                      compSpec, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg Normalization(Double_t scaleFactor,
                                                                                          Int_t scaleType)
                                                                  {
                                                                     return RooCmdArg("Normalization", scaleType, 0,
                                                                                      scaleFactor, 0, nullptr, nullptr,
                                                                                      nullptr, nullptr);
                                                                  }

                                                                  // RooAbsData::plotOn arguments
                                                                  RooCmdArg Cut(const char *cutSpec)
                                                                  {
                                                                     return RooCmdArg("CutSpec", 0, 0, 0, 0, cutSpec,
                                                                                      nullptr, nullptr, nullptr);
                                                                  }
                                                                  RooCmdArg Cut(const RooFormulaVar &cutVar)
                                                                  {
                                                                     return RooCmdArg("CutVar", 0, 0, 0, 0, nullptr,
                                                                                      nullptr, &cutVar, nullptr);
                                                                  }
                                                                  RooCmdArg Binning(const RooAbsBinning &binning)
                                                                  {
                                                                     return RooCmdArg("Binning", 0, 0, 0, 0, nullptr,
                                                                                      nullptr, &binning, nullptr);
                                                                  }
                                                                  RooCmdArg Binning(const char *binningName)
                                                                  {
                                                                     return RooCmdArg("BinningName", 0, 0, 0, 0,
                                                                                      binningName, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg Binning(Int_t nBins, Double_t xlo,
                                                                                    Double_t xhi)
                                                                  {
                                                                     return RooCmdArg("BinningSpec", nBins, 0, xlo, xhi,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg MarkerStyle(Style_t style)
                                                                  {
                                                                     return RooCmdArg("MarkerStyle", style, 0, 0, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg MarkerSize(Size_t size)
                                                                  {
                                                                     return RooCmdArg("MarkerSize", 0, 0, size, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg MarkerColor(Color_t color)
                                                                  {
                                                                     return RooCmdArg("MarkerColor", color, 0, 0, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg CutRange(const char *rangeName)
                                                                  {
                                                                     return RooCmdArg("CutRange", 0, 0, 0, 0, rangeName,
                                                                                      nullptr, nullptr, nullptr);
                                                                  }
                                                                  RooCmdArg AddTo(const char *name)
                                                                  {
                                                                     return RooCmdArg("AddTo", 0, 0, 0, 0, name,
                                                                                      nullptr, nullptr, nullptr);
                                                                  }
                                                                  RooCmdArg XErrorSize(Double_t width)
                                                                  {
                                                                     return RooCmdArg("XErrorSize", 0, 0, width, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg RefreshNorm()
                                                                  {
                                                                     return RooCmdArg("RefreshNorm", 1, 0, 0, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }
                                                                  RooCmdArg Efficiency(const RooCategory &cat)
                                                                  {
                                                                     return RooCmdArg("Efficiency", 0, 0, 0, 0, nullptr,
                                                                                      nullptr, &cat, nullptr);
                                                                  }
                                                                  RooCmdArg Rescale(Double_t factor)
                                                                  {
                                                                     return RooCmdArg("Rescale", 0, 0, factor, 0,
                                                                                      nullptr, nullptr, nullptr,
                                                                                      nullptr);
                                                                  }

                                                                  // RooDataHist::ctor arguments
                                                                  RooCmdArg Weight(Double_t wgt)
                                                                  {
                                                                     return RooCmdArg("Weight", 0, 0, wgt, 0, nullptr,
                                                                                      nullptr, nullptr, nullptr);
                                                                  }
                                                                  RooCmdArg Index(RooCategory &icat)
                                                                  {
                                                                     return RooCmdArg("IndexCat", 0, 0, 0, 0, nullptr,
                                                                                      nullptr, &icat, nullptr);
                                                                  }
                                                                  RooCmdArg Import(const char *state, TH1 &histo)
                                                                  {
                                                                     return RooCmdArg("ImportHistoSlice", 0, 0, 0, 0,
                                                                                      state, nullptr, &histo, nullptr);
                                                                  }
                                                                  RooCmdArg Import(const char *state,
                                                                                   RooDataHist &dhist)
                                                                  {
                                                                     return RooCmdArg("ImportDataHistSlice", 0, 0, 0, 0,
                                                                                      state, nullptr, &dhist, nullptr);
                                                                  }
                                                                  RooCmdArg Import(TH1 &histo, Bool_t importDensity)
                                                                  {
                                                                     return RooCmdArg("ImportHisto", importDensity, 0,
                                                                                      0, 0, nullptr, nullptr, &histo,
                                                                                      nullptr);
                                                                  }

                                                                  RooCmdArg Import(
                                                                     const std::map<std::string, RooDataHist *> &arg) {
                                                                     RooCmdArg container("ImportDataHistSliceMany", 0,
                                                                                         0, 0, 0, nullptr, nullptr,
                                                                                         nullptr, nullptr);
                                                                     std::map<std::string,
                                                                              RooDataHist *>::const_iterator iter;
                                                                     for (iter = arg.begin(); iter != arg.end();
                                                                          ++iter) {
                                                                        container.addArg(Import(iter->first.c_str(),
                                                                                                *(iter->second)));
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }
  RooCmdArg Import(const std::map<std::string,TH1*>& arg) {
     RooCmdArg container("ImportHistoSliceMany", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
     std::map<std::string, TH1 *>::const_iterator iter;
     for (iter = arg.begin(); iter != arg.end(); ++iter) {
        container.addArg(Import(iter->first.c_str(), *(iter->second)));
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }

  
  // RooDataSet::ctor arguments
  RooCmdArg WeightVar(const char* name, Bool_t reinterpretAsWeight) {
     return RooCmdArg("WeightVarName", reinterpretAsWeight, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg WeightVar(const RooRealVar& arg, Bool_t reinterpretAsWeight)  {
     return RooCmdArg("WeightVar", reinterpretAsWeight, 0, 0, 0, nullptr, nullptr, &arg, nullptr);
  }
  RooCmdArg Link(const char* state, RooAbsData& data)   {
     return RooCmdArg("LinkDataSlice", 0, 0, 0, 0, state, nullptr, &data, nullptr);
  }
  RooCmdArg Import(const char* state, RooDataSet& data) {
     return RooCmdArg("ImportDataSlice", 0, 0, 0, 0, state, nullptr, &data, nullptr);
  }
  RooCmdArg Import(RooDataSet& data)                    {
     return RooCmdArg("ImportData", 0, 0, 0, 0, nullptr, nullptr, &data, nullptr);
  }
  RooCmdArg Import(TTree& tree)                         {
     return RooCmdArg("ImportTree", 0, 0, 0, 0, nullptr, nullptr, reinterpret_cast<TObject *>(&tree), nullptr);
  }
  RooCmdArg ImportFromFile(const char* fname, const char* tname){
     return RooCmdArg("ImportFromFile", 0, 0, 0, 0, fname, tname, nullptr, nullptr);
  }
  RooCmdArg StoreError(const RooArgSet& aset)           {
     return RooCmdArg("StoreError", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &aset);
  }
  RooCmdArg StoreAsymError(const RooArgSet& aset)       {
     return RooCmdArg("StoreAsymError", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &aset);
  }
  RooCmdArg OwnLinked()                                 {
     return RooCmdArg("OwnLinked", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  RooCmdArg Import(const std::map<std::string,RooDataSet*>& arg) {
     RooCmdArg container("ImportDataSliceMany", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
     std::map<std::string, RooDataSet *>::const_iterator iter;
     for (iter = arg.begin(); iter != arg.end(); ++iter) {
        container.addArg(Import(iter->first.c_str(), *(iter->second)));
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }
  RooCmdArg Link(const std::map<std::string,RooAbsData*>& arg) {
     RooCmdArg container("LinkDataSliceMany", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
     std::map<std::string, RooAbsData *>::const_iterator iter;
     for (iter = arg.begin(); iter != arg.end(); ++iter) {
        container.addArg(Link(iter->first.c_str(), *(iter->second)));
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }
 

  // RooChi2Var::ctor arguments
  RooCmdArg Extended(Bool_t flag) {
     return RooCmdArg("Extended", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg DataError(Int_t etype) {
     return RooCmdArg("DataError", (Int_t)etype, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg NumCPU(Int_t nCPU, Int_t interleave)   {
     return RooCmdArg("NumCPU", nCPU, interleave, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooAbsCollection::printLatex arguments
  RooCmdArg Columns(Int_t ncol)                           {
     return RooCmdArg("Columns", ncol, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg OutputFile(const char* fileName)              {
     return RooCmdArg("OutputFile", 0, 0, 0, 0, fileName, nullptr, nullptr, nullptr);
  }
  RooCmdArg Sibling(const RooAbsCollection& sibling)      {
     return RooCmdArg("Sibling", 0, 0, 0, 0, nullptr, nullptr, &sibling, nullptr);
  }
  RooCmdArg Format(const char* format, Int_t sigDigit)    {
     return RooCmdArg("Format", sigDigit, 0, 0, 0, format, nullptr, nullptr, nullptr);
  }
  RooCmdArg Format(const char* what, const RooCmdArg& arg1,const RooCmdArg& arg2,const RooCmdArg& arg3,const RooCmdArg& arg4,
                   const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) {
     RooCmdArg ret("FormatArgs", 0, 0, 0, 0, what, nullptr, nullptr, nullptr);
     ret.addArg(arg1);
     ret.addArg(arg2);
     ret.addArg(arg3);
     ret.addArg(arg4);
     ret.addArg(arg5);
     ret.addArg(arg6);
     ret.addArg(arg7);
     ret.addArg(arg8);
     ret.setProcessRecArgs(kFALSE);
     return ret;
  }
  
  // RooAbsRealLValue::frame arguments
  RooCmdArg Title(const char* name) {
     return RooCmdArg("Title", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg Bins(Int_t nbin)        {
     return RooCmdArg("Bins", nbin, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg AutoSymRange(const RooAbsData& data, Double_t marginFactor) {
     return RooCmdArg("AutoRange", 1, 0, marginFactor, 0, nullptr, nullptr, &data, nullptr);
  }
  RooCmdArg AutoRange(const RooAbsData& data, Double_t marginFactor) {
     return RooCmdArg("AutoRange", 0, 0, marginFactor, 0, nullptr, nullptr, &data, nullptr);
  }

  // RooAbsData::reduce arguments
  RooCmdArg SelectVars(const RooArgSet& vars)     {
     return RooCmdArg("SelectVars", 0, 0, 0, 0, nullptr, nullptr, &vars, nullptr);
  }
  RooCmdArg EventRange(Int_t nStart, Int_t nStop) {
     return RooCmdArg("EventRange", nStart, nStop, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooAbsPdf::fitTo arguments
  RooCmdArg FitOptions(const char* opts) {
     return RooCmdArg("FitOptions", 0, 0, 0, 0, opts, nullptr, nullptr, nullptr);
  }
  RooCmdArg Optimize(Int_t flag)         {
     return RooCmdArg("Optimize", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Verbose(Bool_t flag)         {
     return RooCmdArg("Verbose", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Save(Bool_t flag)            {
     return RooCmdArg("Save", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Timer(Bool_t flag)           {
     return RooCmdArg("Timer", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg PrintLevel(Int_t level)      {
     return RooCmdArg("PrintLevel", level, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Warnings(Bool_t flag)        {
     return RooCmdArg("Warnings", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Strategy(Int_t code)         {
     return RooCmdArg("Strategy", code, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg InitialHesse(Bool_t flag)    {
     return RooCmdArg("InitialHesse", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Hesse(Bool_t flag)           {
     return RooCmdArg("Hesse", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Minos(Bool_t flag)           {
     return RooCmdArg("Minos", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Minos(const RooArgSet& minosArgs)            {
     return RooCmdArg("Minos", kTRUE, 0, 0, 0, nullptr, nullptr, &minosArgs, nullptr);
  }
  RooCmdArg ConditionalObservables(const RooArgSet& set) {
     return RooCmdArg("ProjectedObservables", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &set);
  }
  RooCmdArg ProjectedObservables(const RooArgSet& set)   {
     return RooCmdArg("ProjectedObservables", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &set);
  }
  RooCmdArg SplitRange(Bool_t flag)                      {
     return RooCmdArg("SplitRange", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg SumCoefRange(const char* rangeName)          {
     return RooCmdArg("SumCoefRange", 0, 0, 0, 0, rangeName, nullptr, nullptr, nullptr);
  }
  RooCmdArg Constrain(const RooArgSet& params)           {
     return RooCmdArg("Constrain", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &params);
  }
  RooCmdArg GlobalObservables(const RooArgSet& globs)    {
     return RooCmdArg("GlobalObservables", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &globs);
  }
  RooCmdArg GlobalObservablesTag(const char* tagName)    {
     return RooCmdArg("GlobalObservablesTag", 0, 0, 0, 0, tagName, nullptr, nullptr, nullptr);
  }
  RooCmdArg Constrained()                                {
     return RooCmdArg("Constrained", kTRUE, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ExternalConstraints(const RooArgSet& cpdfs)  {
     return RooCmdArg("ExternalConstraints", 0, 0, 0, 0, nullptr, nullptr, &cpdfs, nullptr, nullptr, nullptr, &cpdfs);
  }
  RooCmdArg PrintEvalErrors(Int_t numErrors)             {
     return RooCmdArg("PrintEvalErrors", numErrors, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg EvalErrorWall(Bool_t flag)                   {
     return RooCmdArg("EvalErrorWall", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg SumW2Error(Bool_t flag)                      {
     return RooCmdArg("SumW2Error", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg CloneData(Bool_t flag)                       {
     return RooCmdArg("CloneData", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Integrate(Bool_t flag)                       {
     return RooCmdArg("Integrate", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Minimizer(const char* type, const char* alg) {
     return RooCmdArg("Minimizer", 0, 0, 0, 0, type, alg, nullptr, nullptr);
  }
  RooCmdArg Offset(Bool_t flag)                          {
     return RooCmdArg("OffsetLikelihood", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooAbsPdf::paramOn arguments
  RooCmdArg Label(const char* str) {
     return RooCmdArg("Label", 0, 0, 0, 0, str, nullptr, nullptr, nullptr);
  }
  RooCmdArg Layout(Double_t xmin, Double_t xmax, Double_t ymin) {
     return RooCmdArg("Layout", Int_t(ymin * 10000), 0, xmin, xmax, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Parameters(const RooArgSet& params) {
     return RooCmdArg("Parameters", 0, 0, 0, 0, nullptr, nullptr, &params, nullptr);
  }
  RooCmdArg ShowConstants(Bool_t flag) {
     return RooCmdArg("ShowConstants", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooTreeData::statOn arguments
  RooCmdArg What(const char* str) {
     return RooCmdArg("What", 0, 0, 0, 0, str, nullptr, nullptr, nullptr);
  }

  // RooProdPdf::ctor arguments
  RooCmdArg Conditional(const RooArgSet& pdfSet, const RooArgSet& depSet, Bool_t depsAreCond) {
     return RooCmdArg("Conditional", depsAreCond, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                      &pdfSet, &depSet);
  };

  // RooAbsPdf::generate arguments
  RooCmdArg ProtoData(const RooDataSet& protoData, Bool_t randomizeOrder, Bool_t resample) 
                                         {
                                            return RooCmdArg("PrototypeData", randomizeOrder, resample, 0, 0, nullptr,
                                                             nullptr, &protoData, nullptr);
                                         }
                                         RooCmdArg NumEvents(Int_t numEvents)   {
                                            return RooCmdArg("NumEvents", numEvents, 0, 0, 0, nullptr, nullptr, nullptr,
                                                             nullptr);
                                         }
                                         RooCmdArg NumEvents(Double_t numEvents)
                                         {
                                            return RooCmdArg("NumEventsD", 0, 0, numEvents, 0, nullptr, nullptr,
                                                             nullptr, nullptr);
                                         }
                                         RooCmdArg ExpectedData(Bool_t flag)    {
                                            return RooCmdArg("ExpectedData", flag, 0, 0, 0, nullptr, nullptr, nullptr,
                                                             nullptr);
                                         }
                                         RooCmdArg Asimov(Bool_t flag)          { return ExpectedData(flag) ; }
  RooCmdArg AutoBinned(Bool_t flag)      {
     return RooCmdArg("AutoBinned", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg GenBinned(const char* tag)   {
     return RooCmdArg("GenBinned", 0, 0, 0, 0, tag, nullptr, nullptr, nullptr);
  }
  RooCmdArg AllBinned()                  {
     return RooCmdArg("GenBinned", 0, 0, 0, 0, "*", nullptr, nullptr, nullptr);
  }

  // RooAbsRealLValue::createHistogram arguments
  RooCmdArg YVar(const RooAbsRealLValue& var, const RooCmdArg& arg)       {
     return RooCmdArg("YVar", 0, 0, 0, 0, nullptr, nullptr, &var, nullptr, &arg);
  }
  RooCmdArg ZVar(const RooAbsRealLValue& var, const RooCmdArg& arg)       {
     return RooCmdArg("ZVar", 0, 0, 0, 0, nullptr, nullptr, &var, nullptr, &arg);
  }
  RooCmdArg AxisLabel(const char* name)                                   {
     return RooCmdArg("AxisLabel", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg Scaling(Bool_t flag)                                          {
     return RooCmdArg("Scaling", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooAbsReal::createHistogram arguments
  RooCmdArg IntrinsicBinning(Bool_t flag)                                 {
     return RooCmdArg("IntrinsicBinning", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooAbsData::createHistogram arguments
  RooCmdArg AutoSymBinning(Int_t nbins, Double_t marginFactor) {
     return RooCmdArg("AutoRangeData", 1, nbins, marginFactor, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg AutoBinning(Int_t nbins, Double_t marginFactor) {
     return RooCmdArg("AutoRangeData", 0, nbins, marginFactor, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooAbsReal::fillHistogram arguments
  RooCmdArg IntegratedObservables(const RooArgSet& intObs) {
     return RooCmdArg("IntObs", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &intObs, nullptr);
  };

  // RooAbsReal::createIntegral arguments
  RooCmdArg NormSet(const RooArgSet& nset)           {
     return RooCmdArg("NormSet", 0, 0, 0, 0, nullptr, nullptr, &nset, nullptr);
  }
  RooCmdArg NumIntConfig(const RooNumIntConfig& cfg) {
     return RooCmdArg("NumIntConfig", 0, 0, 0, 0, nullptr, nullptr, &cfg, nullptr);
  }

  // RooMCStudy::ctor arguments
  RooCmdArg Silence(Bool_t flag) {
     return RooCmdArg("Silence", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg FitModel(RooAbsPdf& pdf) {
     return RooCmdArg("FitModel", 0, 0, 0, 0, nullptr, nullptr, &pdf, nullptr);
  }
  RooCmdArg FitOptions(const RooCmdArg& arg1 ,const RooCmdArg& arg2, const RooCmdArg& arg3,
                       const RooCmdArg& arg4, const RooCmdArg& arg5, const RooCmdArg& arg6) {
     RooCmdArg ret("FitOptArgs", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
     ret.addArg(arg1);
     ret.addArg(arg2);
     ret.addArg(arg3);
     ret.addArg(arg4);
     ret.addArg(arg5);
     ret.addArg(arg6);
     ret.setProcessRecArgs(kFALSE);
     return ret; 
  }
  RooCmdArg Binned(Bool_t flag)               {
     return RooCmdArg("Binned", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg BootStrapData(const RooDataSet& dset) {
     return RooCmdArg("BootStrapData", 0, 0, 0, 0, nullptr, nullptr, &dset, nullptr);
  }

  // RooMCStudy::plot* arguments
  RooCmdArg Frame(const RooCmdArg& arg1,const RooCmdArg& arg2,
                  const RooCmdArg& arg3,const RooCmdArg& arg4,
                  const RooCmdArg& arg5,const RooCmdArg& arg6) {
     RooCmdArg ret("FrameArgs", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
     ret.addArg(arg1);
     ret.addArg(arg2);
     ret.addArg(arg3);
     ret.addArg(arg4);
     ret.addArg(arg5);
     ret.addArg(arg6);
     ret.setProcessRecArgs(kFALSE);
     return ret;
  }
  RooCmdArg FrameBins(Int_t nbins)                 {
     return RooCmdArg("Bins", nbins, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg FrameRange(Double_t xlo, Double_t xhi) {
     return RooCmdArg("Range", 0, 0, xlo, xhi, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg FitGauss(Bool_t flag)                  {
     return RooCmdArg("FitGauss", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooRealVar::format arguments
  RooCmdArg ShowName(Bool_t flag)             {
     return RooCmdArg("ShowName", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ShowValue(Bool_t flag)            {
     return RooCmdArg("ShowValue", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ShowError(Bool_t flag)            {
     return RooCmdArg("ShowError", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ShowAsymError(Bool_t flag)        {
     return RooCmdArg("ShowAsymError", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ShowUnit(Bool_t flag)             {
     return RooCmdArg("ShowUnit", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg AutoPrecision(Int_t ndigit)   {
     return RooCmdArg("AutoPrecision", ndigit, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg FixedPrecision(Int_t ndigit)  {
     return RooCmdArg("FixedPrecision", ndigit, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg TLatexStyle(Bool_t flag)      {
     return RooCmdArg("TLatexStyle", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg LatexStyle(Bool_t flag)       {
     return RooCmdArg("LatexStyle", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg LatexTableStyle(Bool_t flag)  {
     return RooCmdArg("LatexTableStyle", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg VerbatimName(Bool_t flag)     {
     return RooCmdArg("VerbatimName", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooMsgService::addReportingStream arguments
  RooCmdArg Topic(Int_t topic)              {
     return RooCmdArg("Topic", topic, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ObjectName(const char* name)    {
     return RooCmdArg("ObjectName", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg ClassName(const char* name)     {
     return RooCmdArg("ClassName", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg BaseClassName(const char* name) {
     return RooCmdArg("BaseClassName", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg TagName(const char* name)     {
     return RooCmdArg("LabelName", 0, 0, 0, 0, name, nullptr, nullptr, nullptr);
  }
  RooCmdArg OutputStream(ostream &os)
  {
     return RooCmdArg("OutputStream", 0, 0, 0, 0, nullptr, nullptr, reinterpret_cast<TObject *>(&os), nullptr);
  }
  RooCmdArg Prefix(Bool_t flag)          {
     return RooCmdArg("Prefix", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg Color(Color_t color)         {
     return RooCmdArg("Color", color, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooWorkspace::import() arguments
  RooCmdArg RenameConflictNodes(const char* suffix, Bool_t ro) {
     return RooCmdArg("RenameConflictNodes", ro, 0, 0, 0, suffix, nullptr, nullptr, nullptr);
  }
  RooCmdArg RecycleConflictNodes(Bool_t flag)               {
     return RooCmdArg("RecycleConflictNodes", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg RenameAllNodes(const char* suffix)              {
     return RooCmdArg("RenameAllNodes", 0, 0, 0, 0, suffix, nullptr, nullptr, nullptr);
  }
  RooCmdArg RenameAllVariables(const char* suffix)          {
     return RooCmdArg("RenameAllVariables", 0, 0, 0, 0, suffix, nullptr, nullptr, nullptr);
  }
  RooCmdArg RenameAllVariablesExcept(const char* suffix, const char* except)    {
     return RooCmdArg("RenameAllVariables", 0, 0, 0, 0, suffix, except, nullptr, nullptr);
  }
  RooCmdArg RenameVariable(const char* in, const char* out) {
     return RooCmdArg("RenameVar", 0, 0, 0, 0, in, out, nullptr, nullptr);
  }
  RooCmdArg Rename(const char* suffix)                      {
     return RooCmdArg("Rename", 0, 0, 0, 0, suffix, nullptr, nullptr, nullptr);
  }
  RooCmdArg Embedded(Bool_t flag)                           {
     return RooCmdArg("Embedded", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg NoRecursion(Bool_t flag)                        {
     return RooCmdArg("NoRecursion", flag, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  // RooSimCloneTool::build() arguments
  RooCmdArg SplitParam(const char* varname, const char* catname)         {
     return RooCmdArg("SplitParam", 0, 0, 0, 0, varname, catname, nullptr, nullptr);
  }
  RooCmdArg SplitParam(const RooRealVar& var, const RooAbsCategory& cat) {
     return RooCmdArg("SplitParam", 0, 0, 0, 0, var.GetName(), cat.GetName(), nullptr, nullptr);
  }
  RooCmdArg SplitParamConstrained(const char* varname, const char* catname, const char* rsname)        {
     return RooCmdArg("SplitParamConstrained", 0, 0, 0, 0, varname, catname, nullptr, nullptr, nullptr, rsname);
  }
  RooCmdArg SplitParamConstrained(const RooRealVar& var, const RooAbsCategory& cat, const char* rsname) {
     return RooCmdArg("SplitParamConstrained", 0, 0, 0, 0, var.GetName(), cat.GetName(), nullptr, nullptr, nullptr,
                      rsname);
  }
  RooCmdArg Restrict(const char* catName, const char* stateNameList) {
     return RooCmdArg("Restrict", 0, 0, 0, 0, catName, stateNameList, nullptr, nullptr);
  }

  // RooAbsPdf::createCdf() arguments
  RooCmdArg SupNormSet(const RooArgSet& nset) {
     return RooCmdArg("SupNormSet", 0, 0, 0, 0, nullptr, nullptr, &nset, nullptr);
  }
  RooCmdArg ScanParameters(Int_t nbins,Int_t intOrder) {
     return RooCmdArg("ScanParameters", nbins, intOrder, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ScanNumCdf() {
     return RooCmdArg("ScanNumCdf", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ScanAllCdf() {
     return RooCmdArg("ScanAllCdf", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }
  RooCmdArg ScanNoCdf() {
     return RooCmdArg("ScanNoCdf", 1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
  }

  RooCmdArg MultiArg(const RooCmdArg& arg1,const RooCmdArg& arg2,const RooCmdArg& arg3,const RooCmdArg& arg4,
                     const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) {
     RooCmdArg ret("MultiArg", 0, 0, 0, 0, nullptr, nullptr, nullptr, nullptr);
     ret.addArg(arg1);
     ret.addArg(arg2);
     ret.addArg(arg3);
     ret.addArg(arg4);
     ret.addArg(arg5);
     ret.addArg(arg6);
     ret.addArg(arg7);
     ret.addArg(arg8);
     ret.setProcessRecArgs(kTRUE, kFALSE);
     return ret;
  }

  RooConstVar& RooConst(Double_t val) { return RooRealConstant::value(val) ; }

 
} // End namespace RooFit

namespace RooFitShortHand {

RooArgSet S(const RooAbsArg& v1) { return RooArgSet(v1) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2) { return RooArgSet(v1,v2) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3) { return RooArgSet(v1,v2,v3) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4) { return RooArgSet(v1,v2,v3,v4) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5) 
          { return RooArgSet(v1,v2,v3,v4,v5) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6) { return RooArgSet(v1,v2,v3,v4,v5,v6) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6, const RooAbsArg& v7) { return RooArgSet(v1,v2,v3,v4,v5,v6,v7) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8) { return RooArgSet(v1,v2,v3,v4,v5,v6,v7,v8) ; }
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8, const RooAbsArg& v9) 
          { return RooArgSet(v1,v2,v3,v4,v5,v6,v7,v8,v9) ; }

RooArgList L(const RooAbsArg& v1) { return RooArgList(v1) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2) { return RooArgList(v1,v2) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3) { return RooArgList(v1,v2,v3) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4) { return RooArgList(v1,v2,v3,v4) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5) 
           { return RooArgList(v1,v2,v3,v4,v5) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6) { return RooArgList(v1,v2,v3,v4,v5,v6) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6, const RooAbsArg& v7) { return RooArgList(v1,v2,v3,v4,v5,v6,v7) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8) { return RooArgList(v1,v2,v3,v4,v5,v6,v7,v8) ; }
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8, const RooAbsArg& v9) 
           { return RooArgList(v1,v2,v3,v4,v5,v6,v7,v8,v9) ; }

RooConstVar& C(Double_t value) { return RooFit::RooConst(value) ; }

} // End namespace Shorthand

