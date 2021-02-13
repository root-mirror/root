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

#include "RooGlobalFunc.h"

#include "RooFit.h"
#include "RooCategory.h"
#include "RooRealConstant.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooNumIntConfig.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooFormulaVar.h"
#include "RooHelpers.h"
#include "TH1.h"

using namespace std;

namespace RooFit {

  // RooAbsReal::plotOn arguments
  RooCmdArg DrawOption(const char* opt)            { return RooCmdArg("DrawOption",0,0,0,0,opt,0,0,0) ; }
  RooCmdArg Slice(const RooArgSet& sliceSet)       { return RooCmdArg("SliceVars",0,0,0,0,0,0,&sliceSet,0) ; }
  RooCmdArg Slice(RooCategory& cat, const char* label) { return RooCmdArg("SliceCat",0,0,0,0,label,0,&cat,0) ; }

  RooCmdArg Project(const RooArgSet& projSet)      { return RooCmdArg("Project",0,0,0,0,0,0,&projSet,0) ; }
  RooCmdArg ProjWData(const RooArgSet& projSet, 
                      const RooAbsData& projData,
                      Bool_t binData)              { return RooCmdArg("ProjData",binData,0,0,0,0,0,&projSet,&projData) ; }
  RooCmdArg ProjWData(const RooAbsData& projData,
                      Bool_t binData)              { return RooCmdArg("ProjData",binData,0,0,0,0,0,0,&projData) ; }
  RooCmdArg Asymmetry(const RooCategory& cat)      { return RooCmdArg("Asymmetry",0,0,0,0,0,0,&cat,0) ; }
  RooCmdArg Precision(Double_t prec)               { return RooCmdArg("Precision",0,0,prec,0,0,0,0,0) ; }
  RooCmdArg ShiftToZero()                          { return RooCmdArg("ShiftToZero",1,0,0,0,0,0,0,0) ; }
  RooCmdArg Normalization(Double_t scaleFactor)    { return RooCmdArg("Normalization",RooAbsReal::Relative,0,scaleFactor,0,0,0,0,0) ; }
  RooCmdArg Range(const char* rangeName, Bool_t adjustNorm)   { return RooCmdArg("RangeWithName",adjustNorm,0,0,0,rangeName,0,0,0) ; }
  RooCmdArg Range(Double_t lo, Double_t hi, Bool_t adjustNorm){ return RooCmdArg("Range",adjustNorm,0,lo,hi,0,0,0,0) ; }
  RooCmdArg NormRange(const char* rangeNameList)   { return RooCmdArg("NormRange",0,0,0,0,rangeNameList,0,0,0) ; }
  RooCmdArg VLines()                               { return RooCmdArg("VLines",1,0,0,0,0,0,0,0) ; } 
  RooCmdArg LineColor(Color_t color)               { return RooCmdArg("LineColor",color,0,0,0,0,0,0,0) ; }
  RooCmdArg LineStyle(Style_t style)               { return RooCmdArg("LineStyle",style,0,0,0,0,0,0,0) ; }
  RooCmdArg LineWidth(Width_t width)               { return RooCmdArg("LineWidth",width,0,0,0,0,0,0,0) ; }
  RooCmdArg FillColor(Color_t color)               { return RooCmdArg("FillColor",color,0,0,0,0,0,0,0) ; }
  RooCmdArg FillStyle(Style_t style)               { return RooCmdArg("FillStyle",style,0,0,0,0,0,0,0) ; }
  RooCmdArg ProjectionRange(const char* rangeName) { return RooCmdArg("ProjectionRange",0,0,0,0,rangeName,0,0,0) ; }
  RooCmdArg Name(const char* name)                 { return RooCmdArg("Name",0,0,0,0,name,0,0,0) ; }
  RooCmdArg Invisible(bool inv)                    { return RooCmdArg("Invisible",inv,0,0,0,0,0,0,0) ; }
  RooCmdArg AddTo(const char* name, double wgtSel, double wgtOther) { return RooCmdArg("AddTo",0,0,wgtSel,wgtOther,name,0,0,0) ; }
  RooCmdArg EvalErrorValue(Double_t val)           { return RooCmdArg("EvalErrorValue",1,0,val,0,0,0,0,0) ; }
  RooCmdArg MoveToBack()                           { return RooCmdArg("MoveToBack",1,0,0,0,0,0,0,0) ; }
  RooCmdArg VisualizeError(const RooFitResult& fitres, Double_t Z, Bool_t EVmethod)  { return RooCmdArg("VisualizeError",EVmethod,0,Z,0,0,0,&fitres,0) ; }
  RooCmdArg VisualizeError(const RooFitResult& fitres, const RooArgSet& param, Double_t Z, Bool_t EVmethod) 
                                                                  { return RooCmdArg("VisualizeError",EVmethod,0,Z,0,0,0,&fitres,0,0,0,&param) ; }
  RooCmdArg VisualizeError(const RooDataSet& paramData, Double_t Z) { return RooCmdArg("VisualizeErrorData",0,0,Z,0,0,0,&paramData,0) ; }
  RooCmdArg ShowProgress()                         { return RooCmdArg("ShowProgress",1,0,0,0,0,0,0,0) ; }
  
  // RooAbsPdf::plotOn arguments
  RooCmdArg Components(const RooArgSet& compSet) { return RooCmdArg("SelectCompSet",0,0,0,0,0,0,&compSet,0) ; }
  RooCmdArg Components(const char* compSpec) { return RooCmdArg("SelectCompSpec",0,0,0,0,compSpec,0,0,0) ; }
  RooCmdArg Normalization(Double_t scaleFactor, Int_t scaleType) 
                                                   { return RooCmdArg("Normalization",scaleType,0,scaleFactor,0,0,0,0,0) ; }
  
  // RooAbsData::plotOn arguments
  RooCmdArg Cut(const char* cutSpec)              { return RooCmdArg("CutSpec",0,0,0,0,cutSpec,0,0,0) ; }
  RooCmdArg Cut(const RooFormulaVar& cutVar)      { return RooCmdArg("CutVar",0,0,0,0,0,0,&cutVar,0) ; }
  RooCmdArg Binning(const RooAbsBinning& binning) { return RooCmdArg("Binning",0,0,0,0,0,0,&binning,0) ;}
  RooCmdArg Binning(const char* binningName) { return RooCmdArg("BinningName",0,0,0,0,binningName,0,0,0) ;}
  RooCmdArg Binning(Int_t nBins, Double_t xlo, Double_t xhi) { return RooCmdArg("BinningSpec",nBins,0,xlo,xhi,0,0,0,0) ;}
  RooCmdArg MarkerStyle(Style_t style)            { return RooCmdArg("MarkerStyle",style,0,0,0,0,0,0,0) ; }
  RooCmdArg MarkerSize(Size_t size)               { return RooCmdArg("MarkerSize",0,0,size,0,0,0,0,0) ; }
  RooCmdArg MarkerColor(Color_t color)            { return RooCmdArg("MarkerColor",color,0,0,0,0,0,0,0) ; }
  RooCmdArg CutRange(const char* rangeName)       { return RooCmdArg("CutRange",0,0,0,0,rangeName,0,0,0) ; }
  RooCmdArg AddTo(const char* name)               { return RooCmdArg("AddTo",0,0,0,0,name,0,0,0) ; }
  RooCmdArg XErrorSize(Double_t width)            { return RooCmdArg("XErrorSize",0,0,width,0,0,0,0,0) ; }
  RooCmdArg RefreshNorm()                         { return RooCmdArg("RefreshNorm",1,0,0,0,0,0,0,0) ; }
  RooCmdArg Efficiency(const RooCategory& cat)    { return RooCmdArg("Efficiency",0,0,0,0,0,0,&cat,0) ; }
  RooCmdArg Rescale(Double_t factor)              { return RooCmdArg("Rescale",0,0,factor,0,0,0,0,0) ; }

  // RooDataHist::ctor arguments
  RooCmdArg Weight(Double_t wgt)                          { return RooCmdArg("Weight",0,0,wgt,0,0,0,0,0) ; }
  RooCmdArg Index(RooCategory& icat)                      { return RooCmdArg("IndexCat",0,0,0,0,0,0,&icat,0) ; }
  RooCmdArg Import(const char* state, TH1& histo)         { return RooCmdArg("ImportHistoSlice",0,0,0,0,state,0,&histo,0) ; }
  RooCmdArg Import(const char* state, RooDataHist& dhist) { return RooCmdArg("ImportDataHistSlice",0,0,0,0,state,0,&dhist,0) ; }
  RooCmdArg Import(TH1& histo, Bool_t importDensity)      { return RooCmdArg("ImportHisto",importDensity,0,0,0,0,0,&histo,0) ; }

  RooCmdArg Import(const std::map<std::string,RooDataHist*>& arg) {
    RooCmdArg container("ImportDataHistSliceMany",0,0,0,0,0,0,0,0) ; 
    std::map<std::string,RooDataHist*>::const_iterator iter ;
    for (iter = arg.begin() ; iter!=arg.end() ; ++iter) {
      container.addArg(Import(iter->first.c_str(),*(iter->second))) ;
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }
  RooCmdArg Import(const std::map<std::string,TH1*>& arg) {
    RooCmdArg container("ImportHistoSliceMany",0,0,0,0,0,0,0,0) ; 
    std::map<std::string,TH1*>::const_iterator iter ;
    for (iter = arg.begin() ; iter!=arg.end() ; ++iter) {
      container.addArg(Import(iter->first.c_str(),*(iter->second))) ;
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }

  
  // RooDataSet::ctor arguments
  RooCmdArg WeightVar(const char* name, Bool_t reinterpretAsWeight) { return RooCmdArg("WeightVarName",reinterpretAsWeight,0,0,0,name,0,0,0) ; }
  RooCmdArg WeightVar(const RooRealVar& arg, Bool_t reinterpretAsWeight)  { return RooCmdArg("WeightVar",reinterpretAsWeight,0,0,0,0,0,&arg,0) ; }
  RooCmdArg Link(const char* state, RooAbsData& data)   { return RooCmdArg("LinkDataSlice",0,0,0,0,state,0,&data,0) ;} 
  RooCmdArg Import(const char* state, RooDataSet& data) { return RooCmdArg("ImportDataSlice",0,0,0,0,state,0,&data,0) ; }
  RooCmdArg Import(RooDataSet& data)                    { return RooCmdArg("ImportData",0,0,0,0,0,0,&data,0) ; }
  RooCmdArg Import(TTree& tree)                         { return RooCmdArg("ImportTree",0,0,0,0,0,0,reinterpret_cast<TObject*>(&tree),0) ; }
  RooCmdArg ImportFromFile(const char* fname, const char* tname){ return RooCmdArg("ImportFromFile",0,0,0,0,fname,tname,0,0) ; }
  RooCmdArg StoreError(const RooArgSet& aset)           { return RooCmdArg("StoreError",0,0,0,0,0,0,0,0,0,0,&aset) ; }
  RooCmdArg StoreAsymError(const RooArgSet& aset)       { return RooCmdArg("StoreAsymError",0,0,0,0,0,0,0,0,0,0,&aset) ; }
  RooCmdArg OwnLinked()                                 { return RooCmdArg("OwnLinked",1,0,0,0,0,0,0,0,0,0,0) ; }

  RooCmdArg Import(const std::map<std::string,RooDataSet*>& arg) {
    RooCmdArg container("ImportDataSliceMany",0,0,0,0,0,0,0,0) ; 
    std::map<std::string,RooDataSet*>::const_iterator iter ;
    for (iter = arg.begin() ; iter!=arg.end() ; ++iter) {
      container.addArg(Import(iter->first.c_str(),*(iter->second))) ;
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }
  RooCmdArg Link(const std::map<std::string,RooAbsData*>& arg) {
    RooCmdArg container("LinkDataSliceMany",0,0,0,0,0,0,0,0) ; 
    std::map<std::string,RooAbsData*>::const_iterator iter ;
    for (iter = arg.begin() ; iter!=arg.end() ; ++iter) {
      container.addArg(Link(iter->first.c_str(),*(iter->second))) ;
    }
    container.setProcessRecArgs(kTRUE,kFALSE) ;
    return container ;
  }
 

  // RooChi2Var::ctor / RooNLLVar arguments
  RooCmdArg Extended(Bool_t flag) { return RooCmdArg("Extended",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg DataError(Int_t etype) { return RooCmdArg("DataError",(Int_t)etype,0,0,0,0,0,0,0) ; }
  RooCmdArg NumCPU(Int_t nCPU, Int_t interleave)   { return RooCmdArg("NumCPU",nCPU,interleave,0,0,0,0,0,0) ; }
  RooCmdArg BatchMode(bool flag) { return RooCmdArg("BatchMode", flag); }
  /// Integrate the PDF over bins. Improves accuracy for binned fits. Switch off using `0.` as argument. \see RooAbsPdf::fitTo().
  RooCmdArg IntegrateBins(double precision) { return RooCmdArg("IntegrateBins", 0, 0, precision); }

  // RooAbsCollection::printLatex arguments
  RooCmdArg Columns(Int_t ncol)                           { return RooCmdArg("Columns",ncol,0,0,0,0,0,0,0) ; }
  RooCmdArg OutputFile(const char* fileName)              { return RooCmdArg("OutputFile",0,0,0,0,fileName,0,0,0) ; }
  RooCmdArg Sibling(const RooAbsCollection& sibling)      { return RooCmdArg("Sibling",0,0,0,0,0,0,&sibling,0) ; }
  RooCmdArg Format(const char* format, Int_t sigDigit)    { return RooCmdArg("Format",sigDigit,0,0,0,format,0,0,0) ; }
  RooCmdArg Format(const char* what, const RooCmdArg& arg1,const RooCmdArg& arg2,const RooCmdArg& arg3,const RooCmdArg& arg4,
                   const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) {
            RooCmdArg ret("FormatArgs",0,0,0,0,what,0,0,0) ; ret.addArg(arg1) ; ret.addArg(arg2) ; 
            ret.addArg(arg3) ; ret.addArg(arg4) ; ret.addArg(arg5) ; ret.addArg(arg6) ; ret.addArg(arg7) ; ret.addArg(arg8) ;
            ret.setProcessRecArgs(kFALSE) ; return ret ;
  }
  
  // RooAbsRealLValue::frame arguments
  RooCmdArg Title(const char* name) { return RooCmdArg("Title",0,0,0,0,name,0,0,0) ; }
  RooCmdArg Bins(Int_t nbin)        { return RooCmdArg("Bins",nbin,0,0,0,0,0,0,0) ; }
  RooCmdArg AutoSymRange(const RooAbsData& data, Double_t marginFactor) { return RooCmdArg("AutoRange",1,0,marginFactor,0,0,0,&data,0) ; }
  RooCmdArg AutoRange(const RooAbsData& data, Double_t marginFactor) { return RooCmdArg("AutoRange",0,0,marginFactor,0,0,0,&data,0) ; }
  
  // RooAbsData::reduce arguments
  RooCmdArg SelectVars(const RooArgSet& vars)     { return RooCmdArg("SelectVars",0,0,0,0,0,0,&vars,0) ; }
  RooCmdArg EventRange(Int_t nStart, Int_t nStop) { return RooCmdArg("EventRange",nStart,nStop,0,0,0,0,0,0) ; }
  
  // RooAbsPdf::fitTo arguments
  RooCmdArg PrefitDataFraction(Double_t data_ratio)  { return RooCmdArg("Prefit",0,0,data_ratio,0,nullptr,nullptr,nullptr,nullptr) ; }
  RooCmdArg FitOptions(const char* opts) { return RooCmdArg("FitOptions",0,0,0,0,opts,0,0,0) ; }
  RooCmdArg Optimize(Int_t flag)         { return RooCmdArg("Optimize",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Verbose(Bool_t flag)         { return RooCmdArg("Verbose",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Save(Bool_t flag)            { return RooCmdArg("Save",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Timer(Bool_t flag)           { return RooCmdArg("Timer",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg PrintLevel(Int_t level)      { return RooCmdArg("PrintLevel",level,0,0,0,0,0,0,0) ; }
  RooCmdArg Warnings(Bool_t flag)        { return RooCmdArg("Warnings",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Strategy(Int_t code)         { return RooCmdArg("Strategy",code,0,0,0,0,0,0,0) ; }
  RooCmdArg InitialHesse(Bool_t flag)    { return RooCmdArg("InitialHesse",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Hesse(Bool_t flag)           { return RooCmdArg("Hesse",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Minos(Bool_t flag)           { return RooCmdArg("Minos",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Minos(const RooArgSet& minosArgs)            { return RooCmdArg("Minos",kTRUE,0,0,0,0,0,&minosArgs,0) ; }
  RooCmdArg ConditionalObservables(const RooArgSet& set) { return RooCmdArg("ProjectedObservables",0,0,0,0,0,0,0,0,0,0,&set) ; }
  RooCmdArg ProjectedObservables(const RooArgSet& set)   { return RooCmdArg("ProjectedObservables",0,0,0,0,0,0,0,0,0,0,&set) ; }
  RooCmdArg SplitRange(Bool_t flag)                      { return RooCmdArg("SplitRange",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg SumCoefRange(const char* rangeName)          { return RooCmdArg("SumCoefRange",0,0,0,0,rangeName,0,0,0) ; }
  RooCmdArg Constrain(const RooArgSet& params)           { return RooCmdArg("Constrain",0,0,0,0,0,0,0,0,0,0,&params) ; }
  RooCmdArg GlobalObservables(const RooArgSet& globs)    { return RooCmdArg("GlobalObservables",0,0,0,0,0,0,0,0,0,0,&globs) ; }
  RooCmdArg GlobalObservablesTag(const char* tagName)    { return RooCmdArg("GlobalObservablesTag",0,0,0,0,tagName,0,0,0) ; }
//  RooCmdArg Constrained()                                { return RooCmdArg("Constrained",kTRUE,0,0,0,0,0,0,0) ; }
  RooCmdArg ExternalConstraints(const RooArgSet& cpdfs)  { return RooCmdArg("ExternalConstraints",0,0,0,0,0,0,&cpdfs,0,0,0,&cpdfs) ; }
  RooCmdArg PrintEvalErrors(Int_t numErrors)             { return RooCmdArg("PrintEvalErrors",numErrors,0,0,0,0,0,0,0) ; }
  RooCmdArg EvalErrorWall(Bool_t flag)                   { return RooCmdArg("EvalErrorWall",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg SumW2Error(Bool_t flag)                      { return RooCmdArg("SumW2Error",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg AsymptoticError(Bool_t flag)                      { return RooCmdArg("AsymptoticError",flag,0,0,0,0,0,0,0) ; }  
  RooCmdArg CloneData(Bool_t flag)                       { return RooCmdArg("CloneData",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Integrate(Bool_t flag)                       { return RooCmdArg("Integrate",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Minimizer(const char* type, const char* alg) { return RooCmdArg("Minimizer",0,0,0,0,type,alg,0,0) ; }
  RooCmdArg Offset(Bool_t flag)                          { return RooCmdArg("OffsetLikelihood",flag,0,0,0,0,0,0,0) ; }
  /// When parameters are chosen such that a PDF is undefined, try to indicate to the minimiser how to leave this region.
  /// \param strength Strength of hints for minimiser. Set to zero to switch off.
  RooCmdArg RecoverFromUndefinedRegions(double strength) { return RooCmdArg("RecoverFromUndefinedRegions",0,0,strength,0,0,0,0,0) ; }

  
  // RooAbsPdf::paramOn arguments
  RooCmdArg Label(const char* str) { return RooCmdArg("Label",0,0,0,0,str,0,0,0) ; }
  RooCmdArg Layout(Double_t xmin, Double_t xmax, Double_t ymin) { return RooCmdArg("Layout",Int_t(ymin*10000),0,xmin,xmax,0,0,0,0) ; }
  RooCmdArg Parameters(const RooArgSet& params) { return RooCmdArg("Parameters",0,0,0,0,0,0,&params,0) ; }
  RooCmdArg ShowConstants(Bool_t flag) { return RooCmdArg("ShowConstants",flag,0,0,0,0,0,0,0) ; }

  // RooTreeData::statOn arguments
  RooCmdArg What(const char* str) { return RooCmdArg("What",0,0,0,0,str,0,0,0) ; }

  // RooProdPdf::ctor arguments
  RooCmdArg Conditional(const RooArgSet& pdfSet, const RooArgSet& depSet, Bool_t depsAreCond) { return RooCmdArg("Conditional",depsAreCond,0,0,0,0,0,0,0,0,0,&pdfSet,&depSet) ; } ;
  
  // RooAbsPdf::generate arguments
  RooCmdArg ProtoData(const RooDataSet& protoData, Bool_t randomizeOrder, Bool_t resample) 
                                         { return RooCmdArg("PrototypeData",randomizeOrder,resample,0,0,0,0,&protoData,0) ; }
  RooCmdArg NumEvents(Int_t numEvents)   { return RooCmdArg("NumEvents",numEvents,0,0,0,0,0,0,0) ; }
  RooCmdArg NumEvents(Double_t numEvents)   { return RooCmdArg("NumEventsD",0,0,numEvents,0,0,0,0,0) ; }
  RooCmdArg ExpectedData(Bool_t flag)    { return RooCmdArg("ExpectedData",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Asimov(Bool_t flag)          { return ExpectedData(flag) ; }
  RooCmdArg AutoBinned(Bool_t flag)      { return RooCmdArg("AutoBinned",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg GenBinned(const char* tag)   { return RooCmdArg("GenBinned",0,0,0,0,tag,0,0,0) ; }
  RooCmdArg AllBinned()                  { return RooCmdArg("GenBinned",0,0,0,0,"*",0,0,0) ; }

  
  // RooAbsRealLValue::createHistogram arguments
  RooCmdArg YVar(const RooAbsRealLValue& var, const RooCmdArg& arg)       { return RooCmdArg("YVar",0,0,0,0,0,0,&var,0,&arg) ; }
  RooCmdArg ZVar(const RooAbsRealLValue& var, const RooCmdArg& arg)       { return RooCmdArg("ZVar",0,0,0,0,0,0,&var,0,&arg) ; }
  RooCmdArg AxisLabel(const char* name)                                   { return RooCmdArg("AxisLabel",0,0,0,0,name,0,0,0) ; }
  RooCmdArg Scaling(Bool_t flag)                                          { return RooCmdArg("Scaling",flag,0,0,0,0,0,0,0) ; }

  // RooAbsReal::createHistogram arguments
  RooCmdArg IntrinsicBinning(Bool_t flag)                                 { return RooCmdArg("IntrinsicBinning",flag,0,0,0,0,0,0,0) ; }

  // RooAbsData::createHistogram arguments
  RooCmdArg AutoSymBinning(Int_t nbins, Double_t marginFactor) { return RooCmdArg("AutoRangeData",1,nbins,marginFactor,0,0,0,0,0) ; }
  RooCmdArg AutoBinning(Int_t nbins, Double_t marginFactor) { return RooCmdArg("AutoRangeData",0,nbins,marginFactor,0,0,0,0,0) ; }

  // RooAbsReal::fillHistogram arguments
  RooCmdArg IntegratedObservables(const RooArgSet& intObs) {  return RooCmdArg("IntObs",0,0,0,0,0,0,0,0,0,0,&intObs,0) ; } ;  
 
  // RooAbsReal::createIntegral arguments
  RooCmdArg NormSet(const RooArgSet& nset)           { return RooCmdArg("NormSet",0,0,0,0,0,0,&nset,0) ; }
  RooCmdArg NumIntConfig(const RooNumIntConfig& cfg) { return RooCmdArg("NumIntConfig",0,0,0,0,0,0,&cfg,0) ; }

  // RooMCStudy::ctor arguments
  RooCmdArg Silence(Bool_t flag) { return RooCmdArg("Silence",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg FitModel(RooAbsPdf& pdf) { return RooCmdArg("FitModel",0,0,0,0,0,0,&pdf,0) ; }
  RooCmdArg FitOptions(const RooCmdArg& arg1 ,const RooCmdArg& arg2, const RooCmdArg& arg3,
                       const RooCmdArg& arg4, const RooCmdArg& arg5, const RooCmdArg& arg6) {
             RooCmdArg ret("FitOptArgs",0,0,0,0,0,0,0,0) ; ret.addArg(arg1) ; ret.addArg(arg2) ; 
             ret.addArg(arg3) ; ret.addArg(arg4) ; ret.addArg(arg5) ; ret.addArg(arg6) ; 
             ret.setProcessRecArgs(kFALSE) ; return ret ; 
  }
  RooCmdArg Binned(Bool_t flag)               { return RooCmdArg("Binned",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg BootStrapData(const RooDataSet& dset) { return RooCmdArg("BootStrapData",0,0,0,0,0,0,&dset,0) ; }


  // RooMCStudy::plot* arguments
  RooCmdArg Frame(const RooCmdArg& arg1,const RooCmdArg& arg2,
                  const RooCmdArg& arg3,const RooCmdArg& arg4,
                  const RooCmdArg& arg5,const RooCmdArg& arg6) {
            RooCmdArg ret("FrameArgs",0,0,0,0,0,0,0,0) ; ret.addArg(arg1) ; ret.addArg(arg2) ; 
            ret.addArg(arg3) ; ret.addArg(arg4) ; ret.addArg(arg5) ; ret.addArg(arg6) ; 
            ret.setProcessRecArgs(kFALSE) ; return ret ;
  }
  RooCmdArg FrameBins(Int_t nbins)                 { return RooCmdArg("Bins",nbins,0,0,0,0,0,0,0) ; }
  RooCmdArg FrameRange(Double_t xlo, Double_t xhi) { return RooCmdArg("Range",0,0,xlo,xhi,0,0,0,0) ; }
  RooCmdArg FitGauss(Bool_t flag)                  { return RooCmdArg("FitGauss",flag,0,0,0,0,0,0,0) ; }

  // RooRealVar::format arguments
  RooCmdArg ShowName(Bool_t flag)             { return RooCmdArg("ShowName",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg ShowValue(Bool_t flag)            { return RooCmdArg("ShowValue",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg ShowError(Bool_t flag)            { return RooCmdArg("ShowError",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg ShowAsymError(Bool_t flag)        { return RooCmdArg("ShowAsymError",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg ShowUnit(Bool_t flag)             { return RooCmdArg("ShowUnit",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg AutoPrecision(Int_t ndigit)   { return RooCmdArg("AutoPrecision",ndigit,0,0,0,0,0,0,0) ; }
  RooCmdArg FixedPrecision(Int_t ndigit)  { return RooCmdArg("FixedPrecision",ndigit,0,0,0,0,0,0,0) ; }
  RooCmdArg TLatexStyle(Bool_t flag)      { return RooCmdArg("TLatexStyle",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg LatexStyle(Bool_t flag)       { return RooCmdArg("LatexStyle",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg LatexTableStyle(Bool_t flag)  { return RooCmdArg("LatexTableStyle",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg VerbatimName(Bool_t flag)     { return RooCmdArg("VerbatimName",flag,0,0,0,0,0,0,0) ; }

  // RooMsgService::addReportingStream arguments
  RooCmdArg Topic(Int_t topic)              { return RooCmdArg("Topic",topic,0,0,0,0,0,0,0) ; }
  RooCmdArg ObjectName(const char* name)    { return RooCmdArg("ObjectName",0,0,0,0,name,0,0,0) ; }
  RooCmdArg ClassName(const char* name)     { return RooCmdArg("ClassName",0,0,0,0,name,0,0,0) ; }
  RooCmdArg BaseClassName(const char* name) { return RooCmdArg("BaseClassName",0,0,0,0,name,0,0,0) ; }
  RooCmdArg TagName(const char* name)     { return RooCmdArg("LabelName",0,0,0,0,name,0,0,0) ; }
  RooCmdArg OutputStream(ostream& os)    { return RooCmdArg("OutputStream",0,0,0,0,0,0,new RooHelpers::WrapIntoTObject<ostream>(os),0) ; }
  RooCmdArg Prefix(Bool_t flag)          { return RooCmdArg("Prefix",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg Color(Color_t color)         { return RooCmdArg("Color",color,0,0,0,0,0,0,0) ; }


  // RooWorkspace::import() arguments
  RooCmdArg RenameConflictNodes(const char* suffix, Bool_t ro) { return RooCmdArg("RenameConflictNodes",ro,0,0,0,suffix,0,0,0) ; }
  RooCmdArg RecycleConflictNodes(Bool_t flag)               { return RooCmdArg("RecycleConflictNodes",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg RenameAllNodes(const char* suffix)              { return RooCmdArg("RenameAllNodes",0,0,0,0,suffix,0,0,0) ; }
  RooCmdArg RenameAllVariables(const char* suffix)          { return RooCmdArg("RenameAllVariables",0,0,0,0,suffix,0,0,0) ; }
  RooCmdArg RenameAllVariablesExcept(const char* suffix, const char* except)    { return RooCmdArg("RenameAllVariables",0,0,0,0,suffix,except,0,0) ; }
  RooCmdArg RenameVariable(const char* in, const char* out) { return RooCmdArg("RenameVar",0,0,0,0,in,out,0,0) ; }
  RooCmdArg Rename(const char* suffix)                      { return RooCmdArg("Rename",0,0,0,0,suffix,0,0,0) ; }
  RooCmdArg Embedded(Bool_t flag)                           { return RooCmdArg("Embedded",flag,0,0,0,0,0,0,0) ; }
  RooCmdArg NoRecursion(Bool_t flag)                        { return RooCmdArg("NoRecursion",flag,0,0,0,0,0,0,0) ; }

  // RooSimCloneTool::build() arguments
  RooCmdArg SplitParam(const char* varname, const char* catname)         { return RooCmdArg("SplitParam",0,0,0,0,varname,catname,0,0) ; }
  RooCmdArg SplitParam(const RooRealVar& var, const RooAbsCategory& cat) { return RooCmdArg("SplitParam",0,0,0,0,var.GetName(),cat.GetName(),0,0) ; }
  RooCmdArg SplitParamConstrained(const char* varname, const char* catname, const char* rsname)        { return RooCmdArg("SplitParamConstrained",0,0,0,0,varname,catname,0,0,0,rsname) ; }
  RooCmdArg SplitParamConstrained(const RooRealVar& var, const RooAbsCategory& cat, const char* rsname) { return RooCmdArg("SplitParamConstrained",0,0,0,0,var.GetName(),cat.GetName(),0,0,0,rsname) ; }
  RooCmdArg Restrict(const char* catName, const char* stateNameList) { return RooCmdArg("Restrict",0,0,0,0,catName,stateNameList,0,0) ; }

  // RooAbsPdf::createCdf() arguments
  RooCmdArg SupNormSet(const RooArgSet& nset) { return RooCmdArg("SupNormSet",0,0,0,0,0,0,&nset,0) ; } 
  RooCmdArg ScanParameters(Int_t nbins,Int_t intOrder) { return RooCmdArg("ScanParameters",nbins,intOrder,0,0,0,0,0,0) ; }
  RooCmdArg ScanNumCdf() { return RooCmdArg("ScanNumCdf",1,0,0,0,0,0,0,0) ; }
  RooCmdArg ScanAllCdf() { return RooCmdArg("ScanAllCdf",1,0,0,0,0,0,0,0) ; }
  RooCmdArg ScanNoCdf() { return RooCmdArg("ScanNoCdf",1,0,0,0,0,0,0,0) ; }


  RooCmdArg MultiArg(const RooCmdArg& arg1,const RooCmdArg& arg2,const RooCmdArg& arg3,const RooCmdArg& arg4,
                     const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) {
    RooCmdArg ret("MultiArg",0,0,0,0,0,0,0,0) ; ret.addArg(arg1) ; ret.addArg(arg2) ; 
    ret.addArg(arg3) ; ret.addArg(arg4) ; ret.addArg(arg5) ; ret.addArg(arg6) ; ret.addArg(arg7) ; ret.addArg(arg8) ;
    ret.setProcessRecArgs(kTRUE,kFALSE) ; return ret ;
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
