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
\file RooAbsData.cxx
\class RooAbsData
\ingroup Roofitcore

RooAbsData is the common abstract base class for binned and unbinned
datasets. The abstract interface defines plotting and tabulating entry
points for its contents and provides an iterator over its elements
(bins for binned data sets, data points for unbinned datasets).
**/

#include "RooAbsData.h"

#include "TBuffer.h"
#include "TClass.h"
#include "TMath.h"
#include "TTree.h"

#include "RooFormulaVar.h"
#include "RooCmdConfig.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"
#include "RooMultiCategory.h"
#include "Roo1DTable.h"
#include "RooAbsDataStore.h"
#include "RooVectorDataStore.h"
#include "RooTreeDataStore.h"
#include "RooDataHist.h"
#include "RooCompositeDataStore.h"
#include "RooCategory.h"
#include "RooTrace.h"
#include "RooUniformBinning.h"

#include "RooRealVar.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooHelpers.h"

#include "TMatrixDSym.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "Math/Util.h"

#include <iostream>
#include <memory>


using namespace std;

ClassImp(RooAbsData);
;

static std::map<RooAbsData*,int> _dcc ;

RooAbsData::StorageType RooAbsData::defaultStorageType=RooAbsData::Vector ;

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::setDefaultStorageType(RooAbsData::StorageType s)
{
   if (RooAbsData::Composite == s) {
      cout << "Composite storage is not a valid *default* storage type." << endl;
   } else {
      defaultStorageType = s;
   }
}

////////////////////////////////////////////////////////////////////////////////

RooAbsData::StorageType RooAbsData::getDefaultStorageType( )
{
  return defaultStorageType;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::claimVars(RooAbsData* data)
{
  _dcc[data]++ ;
  //cout << "RooAbsData(" << data << ") claim incremented to " << _dcc[data] << endl ;
}

////////////////////////////////////////////////////////////////////////////////
/// If return value is true variables can be deleted

Bool_t RooAbsData::releaseVars(RooAbsData* data)
{
  if (_dcc[data]>0) {
    _dcc[data]-- ;
  }

  //cout << "RooAbsData(" << data << ") claim decremented to " << _dcc[data] << endl ;
  return (_dcc[data]==0) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooAbsData::RooAbsData()
{
  claimVars(this) ;
  _dstore = 0 ;
  storageType = defaultStorageType;

  RooTrace::create(this) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a set of variables. Only fundamental elements of vars
/// (RooRealVar,RooCategory etc) are stored as part of the dataset

RooAbsData::RooAbsData(std::string_view name, std::string_view title, const RooArgSet& vars, RooAbsDataStore* dstore) :
  TNamed(TString{name},TString{title}),
  _vars("Dataset Variables"),
  _cachedVars("Cached Variables"),
  _dstore(dstore)
{
   if (dynamic_cast<RooTreeDataStore *>(dstore)) {
      storageType = RooAbsData::Tree;
   } else if (dynamic_cast<RooVectorDataStore *>(dstore)) {
      storageType = RooAbsData::Vector;
   } else {
      storageType = RooAbsData::Composite;
   }
   // cout << "created dataset " << this << endl ;
   claimVars(this);

   // clone the fundamentals of the given data set into internal buffer
   for (const auto var : vars) {
      if (!var->isFundamental()) {
         coutE(InputArguments) << "RooAbsDataStore::initialize(" << GetName()
                               << "): Data set cannot contain non-fundamental types, ignoring " << var->GetName()
                               << endl;
         throw std::invalid_argument(std::string("Only fundamental variables can be placed into datasets. This is violated for ") + var->GetName());
      } else {
         _vars.addClone(*var);
      }
   }

   // reconnect any parameterized ranges to internal dataset observables
   for (auto var : _vars) {
      var->attachArgs(_vars);
   }

   RooTrace::create(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAbsData::RooAbsData(const RooAbsData& other, const char* newname) :
  TNamed(newname?newname:other.GetName(),other.GetTitle()),
  RooPrintable(other), _vars(),
  _cachedVars("Cached Variables")
{
  //cout << "created dataset " << this << endl ;
  claimVars(this) ;
  _vars.addClone(other._vars) ;

  // reconnect any parameterized ranges to internal dataset observables
  for (auto var : _vars) {
    var->attachArgs(_vars);
  }


  if (other._ownedComponents.size()>0) {

    // copy owned components here

    map<string,RooAbsDataStore*> smap ;
    for (auto& itero : other._ownedComponents) {
      RooAbsData* dclone = (RooAbsData*) itero.second->Clone();
      _ownedComponents[itero.first] = dclone;
      smap[itero.first] = dclone->store();
    }

    RooCategory* idx = (RooCategory*) _vars.find(*((RooCompositeDataStore*)other.store())->index()) ;
    _dstore = new RooCompositeDataStore(newname?newname:other.GetName(),other.GetTitle(),_vars,*idx,smap) ;
    storageType = RooAbsData::Composite;

  } else {

    // Convert to vector store if default is vector
    _dstore = other._dstore->clone(_vars,newname?newname:other.GetName()) ;
    storageType = other.storageType;
  }

  RooTrace::create(this) ;
}

RooAbsData& RooAbsData::operator=(const RooAbsData& other) {
  TNamed::operator=(other);
  RooPrintable::operator=(other);

  claimVars(this);
  _vars.Clear();
  _vars.addClone(other._vars);

  // reconnect any parameterized ranges to internal dataset observables
  for (const auto var : _vars) {
    var->attachDataSet(*this) ;
  }


  if (other._ownedComponents.size()>0) {

    // copy owned components here

    map<string,RooAbsDataStore*> smap ;
    for (auto& itero : other._ownedComponents) {
      RooAbsData* dclone = (RooAbsData*) itero.second->Clone();
      _ownedComponents[itero.first] = dclone;
      smap[itero.first] = dclone->store();
    }

    RooCategory* idx = (RooCategory*) _vars.find(*((RooCompositeDataStore*)other.store())->index()) ;
    _dstore = new RooCompositeDataStore(GetName(), GetTitle(), _vars, *idx, smap);
    storageType = RooAbsData::Composite;

  } else {

    // Convert to vector store if default is vector
    _dstore = other._dstore->clone(_vars);
    storageType = other.storageType;
  }

  return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAbsData::~RooAbsData()
{
  if (releaseVars(this)) {
    // will cause content to be deleted subsequently in dtor
  } else {
    _vars.releaseOwnership() ;
  }

  // delete owned contents.
  delete _dstore ;

  // Delete owned dataset components
  for(map<std::string,RooAbsData*>::iterator iter = _ownedComponents.begin() ; iter!= _ownedComponents.end() ; ++iter) {
    delete iter->second ;
  }

  RooTrace::destroy(this) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert tree-based storage to vector-based storage

void RooAbsData::convertToVectorStore()
{
   if (storageType == RooAbsData::Tree) {
      RooVectorDataStore *newStore = new RooVectorDataStore(*(RooTreeDataStore *)_dstore, _vars, GetName());
      delete _dstore;
      _dstore = newStore;
      storageType = RooAbsData::Vector;
   }
}

////////////////////////////////////////////////////////////////////////////////

Bool_t RooAbsData::changeObservableName(const char* from, const char* to)
{
  Bool_t ret =  _dstore->changeObservableName(from,to) ;

  RooAbsArg* tmp = _vars.find(from) ;
  if (tmp) {
    tmp->SetName(to) ;
  }
  return ret ;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::fill()
{
  _dstore->fill() ;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooAbsData::numEntries() const
{
  return nullptr != _dstore ? _dstore->numEntries() : 0;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::reset()
{
  _dstore->reset() ;
}

////////////////////////////////////////////////////////////////////////////////

const RooArgSet* RooAbsData::get(Int_t index) const
{
  checkInit() ;
  return _dstore->get(index) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method -- Cache given set of functions with data

void RooAbsData::cacheArgs(const RooAbsArg* cacheOwner, RooArgSet& varSet, const RooArgSet* nset, Bool_t skipZeroWeights)
{
  _dstore->cacheArgs(cacheOwner,varSet,nset,skipZeroWeights) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method -- Remove cached function values

void RooAbsData::resetCache()
{
  _dstore->resetCache() ;
  _cachedVars.removeAll() ;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method -- Attach dataset copied with cache contents to copied instances of functions

void RooAbsData::attachCache(const RooAbsArg* newOwner, const RooArgSet& cachedVars)
{
  _dstore->attachCache(newOwner, cachedVars) ;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::setArgStatus(const RooArgSet& set, Bool_t active)
{
  _dstore->setArgStatus(set,active) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Control propagation of dirty flags from observables in dataset

void RooAbsData::setDirtyProp(Bool_t flag)
{
  _dstore->setDirtyProp(flag) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a reduced copy of this dataset. The caller takes ownership of the returned dataset
///
/// The following optional named arguments are accepted
/// <table>
/// <tr><td> `SelectVars(const RooArgSet& vars)`   <td> Only retain the listed observables in the output dataset
/// <tr><td> `Cut(const char* expression)`   <td> Only retain event surviving the given cut expression
/// <tr><td> `Cut(const RooFormulaVar& expr)`   <td> Only retain event surviving the given cut formula
/// <tr><td> `CutRange(const char* name)`   <td> Only retain events inside range with given name. Multiple CutRange
///     arguments may be given to select multiple ranges
/// <tr><td> `EventRange(int lo, int hi)`   <td> Only retain events with given sequential event numbers
/// <tr><td> `Name(const char* name)`   <td> Give specified name to output dataset
/// <tr><td> `Title(const char* name)`   <td> Give specified title to output dataset
/// </table>

RooAbsData* RooAbsData::reduce(const RooCmdArg& arg1,const RooCmdArg& arg2,const RooCmdArg& arg3,const RooCmdArg& arg4,
                const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8)
{
  // Define configuration for this method
  RooCmdConfig pc(Form("RooAbsData::reduce(%s)",GetName())) ;
  pc.defineString("name","Name",0,"") ;
  pc.defineString("title","Title",0,"") ;
  pc.defineString("cutRange","CutRange",0,"") ;
  pc.defineString("cutSpec","CutSpec",0,"") ;
  pc.defineObject("cutVar","CutVar",0,0) ;
  pc.defineInt("evtStart","EventRange",0,0) ;
  pc.defineInt("evtStop","EventRange",1,std::numeric_limits<int>::max()) ;
  pc.defineObject("varSel","SelectVars",0,0) ;
  pc.defineMutex("CutVar","CutSpec") ;

  // Process & check varargs
  pc.process(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) ;
  if (!pc.ok(kTRUE)) {
    return 0 ;
  }

  // Extract values from named arguments
  const char* cutRange = pc.getString("cutRange",0,kTRUE) ;
  const char* cutSpec = pc.getString("cutSpec",0,kTRUE) ;
  RooFormulaVar* cutVar = static_cast<RooFormulaVar*>(pc.getObject("cutVar",0)) ;
  Int_t nStart = pc.getInt("evtStart",0) ;
  Int_t nStop = pc.getInt("evtStop",std::numeric_limits<int>::max()) ;
  RooArgSet* varSet = static_cast<RooArgSet*>(pc.getObject("varSel")) ;
  const char* name = pc.getString("name",0,kTRUE) ;
  const char* title = pc.getString("title",0,kTRUE) ;

  // Make sure varSubset doesn't contain any variable not in this dataset
  RooArgSet varSubset ;
  if (varSet) {
    varSubset.add(*varSet) ;
    for (const auto arg : varSubset) {
      if (!_vars.find(arg->GetName())) {
        coutW(InputArguments) << "RooAbsData::reduce(" << GetName() << ") WARNING: variable "
            << arg->GetName() << " not in dataset, ignored" << endl ;
        varSubset.remove(*arg) ;
      }
    }
  } else {
    varSubset.add(*get()) ;
  }

  RooAbsData* ret = 0 ;
  if (cutSpec) {

    RooFormulaVar cutVarTmp(cutSpec,cutSpec,*get()) ;
    ret =  reduceEng(varSubset,&cutVarTmp,cutRange,nStart,nStop,kFALSE) ;

  } else if (cutVar) {

    ret = reduceEng(varSubset,cutVar,cutRange,nStart,nStop,kFALSE) ;

  } else {

    ret = reduceEng(varSubset,0,cutRange,nStart,nStop,kFALSE) ;

  }

  if (!ret) return 0 ;

  if (name) {
    ret->SetName(name) ;
  }
  if (title) {
    ret->SetTitle(title) ;
  }

  return ret ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the data set by applying the given cut on the data points.
/// The cut expression can refer to any variable in the data set. For cuts involving
/// other variables, such as intermediate formula objects, use the equivalent
/// reduce method specifying the as a RooFormulVar reference.

RooAbsData* RooAbsData::reduce(const char* cut)
{
  RooFormulaVar cutVar(cut,cut,*get()) ;
  return reduceEng(*get(),&cutVar,0,0,std::numeric_limits<std::size_t>::max(),kFALSE) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the data set by applying the given cut on the data points.
/// The 'cutVar' formula variable is used to select the subset of data points to be
/// retained in the reduced data collection.

RooAbsData* RooAbsData::reduce(const RooFormulaVar& cutVar)
{
  return reduceEng(*get(),&cutVar,0,0,std::numeric_limits<std::size_t>::max(),kFALSE) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the data set by applying the given cut on the data points
/// and reducing the dimensions to the specified set.
///
/// The cut expression can refer to any variable in the data set. For cuts involving
/// other variables, such as intermediate formula objects, use the equivalent
/// reduce method specifying the as a RooFormulVar reference.

RooAbsData* RooAbsData::reduce(const RooArgSet& varSubset, const char* cut)
{
  // Make sure varSubset doesn't contain any variable not in this dataset
  RooArgSet varSubset2(varSubset) ;
  for (const auto arg : varSubset) {
    if (!_vars.find(arg->GetName())) {
      coutW(InputArguments) << "RooAbsData::reduce(" << GetName() << ") WARNING: variable "
             << arg->GetName() << " not in dataset, ignored" << endl ;
      varSubset2.remove(*arg) ;
    }
  }

  if (cut && strlen(cut)>0) {
    RooFormulaVar cutVar(cut, cut, *get(), false);
    return reduceEng(varSubset2,&cutVar,0,0,std::numeric_limits<std::size_t>::max(),kFALSE) ;
  }
  return reduceEng(varSubset2,0,0,0,std::numeric_limits<std::size_t>::max(),kFALSE) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the data set by applying the given cut on the data points
/// and reducing the dimensions to the specified set.
///
/// The 'cutVar' formula variable is used to select the subset of data points to be
/// retained in the reduced data collection.

RooAbsData* RooAbsData::reduce(const RooArgSet& varSubset, const RooFormulaVar& cutVar)
{
  // Make sure varSubset doesn't contain any variable not in this dataset
  RooArgSet varSubset2(varSubset) ;
  TIterator* iter = varSubset.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!_vars.find(arg->GetName())) {
      coutW(InputArguments) << "RooAbsData::reduce(" << GetName() << ") WARNING: variable "
             << arg->GetName() << " not in dataset, ignored" << endl ;
      varSubset2.remove(*arg) ;
    }
  }
  delete iter ;

  return reduceEng(varSubset2,&cutVar,0,0,std::numeric_limits<std::size_t>::max(),kFALSE) ;
}


RooPlot* RooAbsData::plotOn(RooPlot* frame, const RooCmdArg& arg1, const RooCmdArg& arg2,
             const RooCmdArg& arg3, const RooCmdArg& arg4, const RooCmdArg& arg5,
             const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8) const
{
  RooLinkedList l ;
  l.Add((TObject*)&arg1) ;  l.Add((TObject*)&arg2) ;
  l.Add((TObject*)&arg3) ;  l.Add((TObject*)&arg4) ;
  l.Add((TObject*)&arg5) ;  l.Add((TObject*)&arg6) ;
  l.Add((TObject*)&arg7) ;  l.Add((TObject*)&arg8) ;
  return plotOn(frame,l) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and fill a ROOT histogram TH1,TH2 or TH3 with the values of this dataset for the variables with given names
/// The range of each observable that is histogrammed is always automatically calculated from the distribution in
/// the dataset. The number of bins can be controlled using the [xyz]bins parameters. For a greater degree of control
/// use the createHistogram() method below with named arguments
///
/// The caller takes ownership of the returned histogram

TH1 *RooAbsData::createHistogram(const char* varNameList, Int_t xbins, Int_t ybins, Int_t zbins) const
{
  // Parse list of variable names
  const auto varNames = RooHelpers::tokenise(varNameList, ",:");
  RooLinkedList argList;
  RooRealVar* vars[3] = {nullptr, nullptr, nullptr};

  for (unsigned int i = 0; i < varNames.size(); ++i) {
    if (i >= 3) {
      coutW(InputArguments) << "RooAbsData::createHistogram(" << GetName() << "): Can only create 3-dimensional histograms. Variable "
          << i << " " << varNames[i] << " unused." << std::endl;
      continue;
    }

    vars[i] = static_cast<RooRealVar*>( get()->find(varNames[i].data()) );
    if (!vars[i]) {
      coutE(InputArguments) << "RooAbsData::createHistogram(" << GetName() << ") ERROR: dataset does not contain an observable named " << varNames[i] << std::endl;
      return nullptr;
    }
  }

  if (!vars[0]) {
    coutE(InputArguments) << "RooAbsData::createHistogram(" << GetName() << "): No variable to be histogrammed in list '" << varNameList << "'" << std::endl;
    return nullptr;
  }

  if (xbins<=0  || !vars[0]->hasMax() || !vars[0]->hasMin() ) {
    argList.Add(RooFit::AutoBinning(xbins==0?vars[0]->numBins():abs(xbins)).Clone()) ;
  } else {
    argList.Add(RooFit::Binning(xbins).Clone()) ;
  }

  if (vars[1]) {
    if (ybins<=0 || !vars[1]->hasMax() || !vars[1]->hasMin() ) {
      argList.Add(RooFit::YVar(*vars[1],RooFit::AutoBinning(ybins==0?vars[1]->numBins():abs(ybins))).Clone()) ;
    } else {
      argList.Add(RooFit::YVar(*vars[1],RooFit::Binning(ybins)).Clone()) ;
    }
  }
  if (vars[2]) {
    if (zbins<=0 || !vars[2]->hasMax() || !vars[2]->hasMin() ) {
      argList.Add(RooFit::ZVar(*vars[2],RooFit::AutoBinning(zbins==0?vars[2]->numBins():abs(zbins))).Clone()) ;
    } else {
      argList.Add(RooFit::ZVar(*vars[2],RooFit::Binning(zbins)).Clone()) ;
    }
  }


  // Call implementation function
  TH1* result = createHistogram(GetName(), *vars[0], argList);

  // Delete temporary list of RooCmdArgs
  argList.Delete() ;

  return result ;
}


TH1 *RooAbsData::createHistogram(const char *name, const RooAbsRealLValue& xvar,
             const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
             const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8) const
{
  RooLinkedList l ;
  l.Add((TObject*)&arg1) ;  l.Add((TObject*)&arg2) ;
  l.Add((TObject*)&arg3) ;  l.Add((TObject*)&arg4) ;
  l.Add((TObject*)&arg5) ;  l.Add((TObject*)&arg6) ;
  l.Add((TObject*)&arg7) ;  l.Add((TObject*)&arg8) ;

  return createHistogram(name,xvar,l) ;
}

////////////////////////////////////////////////////////////////////////////////
///
/// This function accepts the following arguments
///
/// \param[in] name Name of the ROOT histogram
/// \param[in] xvar Observable to be mapped on x axis of ROOT histogram
/// \return Histogram now owned by user.
///
/// <table>
/// <tr><td> `AutoBinning(Int_t nbins, Double_y margin)`   <td> Automatically calculate range with given added fractional margin, set binning to nbins
/// <tr><td> `AutoSymBinning(Int_t nbins, Double_y margin)`   <td> Automatically calculate range with given added fractional margin,
///     with additional constraint that mean of data is in center of range, set binning to nbins
/// <tr><td> `Binning(const char* name)`   <td> Apply binning with given name to x axis of histogram
/// <tr><td> `Binning(RooAbsBinning& binning)`   <td> Apply specified binning to x axis of histogram
/// <tr><td> `Binning(int nbins, double lo, double hi)`   <td> Apply specified binning to x axis of histogram
///
/// <tr><td> `YVar(const RooAbsRealLValue& var,...)`   <td> Observable to be mapped on y axis of ROOT histogram
/// <tr><td> `ZVar(const RooAbsRealLValue& var,...)`   <td> Observable to be mapped on z axis of ROOT histogram
/// </table>
///
/// The YVar() and ZVar() arguments can be supplied with optional Binning() Auto(Sym)Range() arguments to control the binning of the Y and Z axes, e.g.
/// ```
/// createHistogram("histo",x,Binning(-1,1,20), YVar(y,Binning(-1,1,30)), ZVar(z,Binning("zbinning")))
/// ```
///
/// The caller takes ownership of the returned histogram

TH1 *RooAbsData::createHistogram(const char *name, const RooAbsRealLValue& xvar, const RooLinkedList& argListIn) const
{
  RooLinkedList argList(argListIn) ;

  // Define configuration for this method
  RooCmdConfig pc(Form("RooAbsData::createHistogram(%s)",GetName())) ;
  pc.defineString("cutRange","CutRange",0,"",kTRUE) ;
  pc.defineString("cutString","CutSpec",0,"") ;
  pc.defineObject("yvar","YVar",0,0) ;
  pc.defineObject("zvar","ZVar",0,0) ;
  pc.allowUndefined() ;

  // Process & check varargs
  pc.process(argList) ;
  if (!pc.ok(kTRUE)) {
    return 0 ;
  }

  const char* cutSpec = pc.getString("cutString",0,kTRUE) ;
  const char* cutRange = pc.getString("cutRange",0,kTRUE) ;

  RooArgList vars(xvar) ;
  RooAbsArg* yvar = static_cast<RooAbsArg*>(pc.getObject("yvar")) ;
  if (yvar) {
    vars.add(*yvar) ;
  }
  RooAbsArg* zvar = static_cast<RooAbsArg*>(pc.getObject("zvar")) ;
  if (zvar) {
    vars.add(*zvar) ;
  }

  pc.stripCmdList(argList,"CutRange,CutSpec") ;

  // Swap Auto(Sym)RangeData with a Binning command
  RooLinkedList ownedCmds ;
  RooCmdArg* autoRD = (RooCmdArg*) argList.find("AutoRangeData") ;
  if (autoRD) {
    Double_t xmin,xmax ;
    if (!getRange((RooRealVar&)xvar,xmin,xmax,autoRD->getDouble(0),autoRD->getInt(0))) {
       RooCmdArg* bincmd = (RooCmdArg*) RooFit::Binning(autoRD->getInt(1),xmin,xmax).Clone() ;
       ownedCmds.Add(bincmd) ;
       argList.Replace(autoRD,bincmd) ;
    }
  }

  if (yvar) {
    RooCmdArg* autoRDY = (RooCmdArg*) ((RooCmdArg*)argList.find("YVar"))->subArgs().find("AutoRangeData") ;
    if (autoRDY) {
      Double_t ymin,ymax ;
      if (!getRange((RooRealVar&)(*yvar),ymin,ymax,autoRDY->getDouble(0),autoRDY->getInt(0))) {
         RooCmdArg* bincmd = (RooCmdArg*) RooFit::Binning(autoRDY->getInt(1),ymin,ymax).Clone() ;
         //ownedCmds.Add(bincmd) ;
         ((RooCmdArg*)argList.find("YVar"))->subArgs().Replace(autoRDY,bincmd) ;
      }
      delete autoRDY ;
    }
  }

  if (zvar) {
    RooCmdArg* autoRDZ = (RooCmdArg*) ((RooCmdArg*)argList.find("ZVar"))->subArgs().find("AutoRangeData") ;
    if (autoRDZ) {
      Double_t zmin,zmax ;
      if (!getRange((RooRealVar&)(*zvar),zmin,zmax,autoRDZ->getDouble(0),autoRDZ->getInt(0))) {
         RooCmdArg* bincmd = (RooCmdArg*) RooFit::Binning(autoRDZ->getInt(1),zmin,zmax).Clone() ;
         //ownedCmds.Add(bincmd) ;
         ((RooCmdArg*)argList.find("ZVar"))->subArgs().Replace(autoRDZ,bincmd) ;
      }
      delete autoRDZ ;
    }
  }


  TH1* histo = xvar.createHistogram(name,argList) ;
  fillHistogram(histo,vars,cutSpec,cutRange) ;

  ownedCmds.Delete() ;

  return histo ;
}

////////////////////////////////////////////////////////////////////////////////
/// Construct table for product of categories in catSet

Roo1DTable* RooAbsData::table(const RooArgSet& catSet, const char* cuts, const char* opts) const
{
  RooArgSet catSet2 ;

  string prodName("(") ;
  TIterator* iter = catSet.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (dynamic_cast<RooAbsCategory*>(arg)) {
      RooAbsCategory* varsArg = dynamic_cast<RooAbsCategory*>(_vars.find(arg->GetName())) ;
      if (varsArg != 0) catSet2.add(*varsArg) ;
      else catSet2.add(*arg) ;
      if (prodName.length()>1) {
   prodName += " x " ;
      }
      prodName += arg->GetName() ;
    } else {
      coutW(InputArguments) << "RooAbsData::table(" << GetName() << ") non-RooAbsCategory input argument " << arg->GetName() << " ignored" << endl ;
    }
  }
  prodName += ")" ;
  delete iter ;

  RooMultiCategory tmp(prodName.c_str(),prodName.c_str(),catSet2) ;
  return table(tmp,cuts,opts) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Print name of dataset

void RooAbsData::printName(ostream& os) const
{
  os << GetName() ;
}

////////////////////////////////////////////////////////////////////////////////
/// Print title of dataset

void RooAbsData::printTitle(ostream& os) const
{
  os << GetTitle() ;
}

////////////////////////////////////////////////////////////////////////////////
/// Print class name of dataset

void RooAbsData::printClassName(ostream& os) const
{
  os << IsA()->GetName() ;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::printMultiline(ostream& os, Int_t contents, Bool_t verbose, TString indent) const
{
  _dstore->printMultiline(os,contents,verbose,indent) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Define default print options, for a given print style

Int_t RooAbsData::defaultPrintContents(Option_t* /*opt*/) const
{
  return kName|kClassName|kArgs|kValue ;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate standardized moment.
///
/// \param[in] var Variable to be used for calculating the moment.
/// \param[in] order Order of the moment.
/// \param[in] cutSpec  If specified, the moment is calculated on the subset of the data which pass the C++ cut specification expression 'cutSpec'
/// \param[in] cutRange If specified, calculate inside the range named 'cutRange' (also applies cut spec)
/// \return \f$ \frac{\left< \left( X - \left< X \right> \right)^n \right>}{\sigma^n} \f$,  where n = order.

Double_t RooAbsData::standMoment(const RooRealVar &var, Double_t order, const char* cutSpec, const char* cutRange) const
{
  // Hardwire invariant answer for first and second moment
  if (order==1) return 0 ;
  if (order==2) return 1 ;

  return moment(var,order,cutSpec,cutRange) / TMath::Power(sigma(var,cutSpec,cutRange),order) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate moment of requested order.
///
/// \param[in] var Variable to be used for calculating the moment.
/// \param[in] order Order of the moment.
/// \param[in] cutSpec  If specified, the moment is calculated on the subset of the data which pass the C++ cut specification expression 'cutSpec'
/// \param[in] cutRange If specified, calculate inside the range named 'cutRange' (also applies cut spec)
/// \return \f$ \left< \left( X - \left< X \right> \right)^n \right> \f$ of order \f$n\f$.
///

Double_t RooAbsData::moment(const RooRealVar& var, Double_t order, const char* cutSpec, const char* cutRange) const
{
  Double_t offset = order>1 ? moment(var,1,cutSpec,cutRange) : 0 ;
  return moment(var,order,offset,cutSpec,cutRange) ;

}

////////////////////////////////////////////////////////////////////////////////
/// Return the 'order'-ed moment of observable 'var' in this dataset. If offset is non-zero it is subtracted
/// from the values of 'var' prior to the moment calculation. If cutSpec and/or cutRange are specified
/// the moment is calculated on the subset of the data which pass the C++ cut specification expression 'cutSpec'
/// and/or are inside the range named 'cutRange'

Double_t RooAbsData::moment(const RooRealVar& var, Double_t order, Double_t offset, const char* cutSpec, const char* cutRange) const
{
  // Lookup variable in dataset
  auto arg = _vars.find(var.GetName());
  if (!arg) {
    coutE(InputArguments) << "RooDataSet::moment(" << GetName() << ") ERROR: unknown variable: " << var.GetName() << endl ;
    return 0;
  }

  auto varPtr = dynamic_cast<const RooRealVar*>(arg);
  // Check if found variable is of type RooRealVar
  if (!varPtr) {
    coutE(InputArguments) << "RooDataSet::moment(" << GetName() << ") ERROR: variable " << var.GetName() << " is not of type RooRealVar" << endl ;
    return 0;
  }

  // Check if dataset is not empty
  if(sumEntries(cutSpec, cutRange) == 0.) {
    coutE(InputArguments) << "RooDataSet::moment(" << GetName() << ") WARNING: empty dataset" << endl ;
    return 0;
  }

  // Setup RooFormulaVar for cutSpec if it is present
  std::unique_ptr<RooFormula> select;
  if (cutSpec) {
    select.reset(new RooFormula("select",cutSpec,*get()));
  }


  // Calculate requested moment
  ROOT::Math::KahanSum<double> sum;
  for(Int_t index= 0; index < numEntries(); index++) {
    const RooArgSet* vars = get(index) ;
    if (select && select->eval()==0) continue ;
    if (cutRange && vars->allInRange(cutRange)) continue ;

    sum += weight() * TMath::Power(varPtr->getVal() - offset,order);
  }

  return sum/sumEntries(cutSpec, cutRange);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method to check if given RooRealVar maps to a RooRealVar in this dataset

RooRealVar* RooAbsData::dataRealVar(const char* methodname, const RooRealVar& extVar) const
{
  // Lookup variable in dataset
  RooRealVar *xdata = (RooRealVar*) _vars.find(extVar.GetName());
  if(!xdata) {
    coutE(InputArguments) << "RooDataSet::" << methodname << "(" << GetName() << ") ERROR: variable : " << extVar.GetName() << " is not in data" << endl ;
    return 0;
  }
  // Check if found variable is of type RooRealVar
  if (!dynamic_cast<RooRealVar*>(xdata)) {
    coutE(InputArguments) << "RooDataSet::" << methodname << "(" << GetName() << ") ERROR: variable : " << extVar.GetName() << " is not of type RooRealVar in data" << endl ;
    return 0;
  }
  return xdata;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method to calculate single correlation and covariance elements

Double_t RooAbsData::corrcov(const RooRealVar &x, const RooRealVar &y, const char* cutSpec, const char* cutRange, Bool_t corr) const
{
  // Lookup variable in dataset
  RooRealVar *xdata = dataRealVar(corr?"correlation":"covariance",x) ;
  RooRealVar *ydata = dataRealVar(corr?"correlation":"covariance",y) ;
  if (!xdata||!ydata) return 0 ;

  // Check if dataset is not empty
  if(sumEntries(cutSpec, cutRange) == 0.) {
    coutW(InputArguments) << "RooDataSet::" << (corr?"correlation":"covariance") << "(" << GetName() << ") WARNING: empty dataset, returning zero" << endl ;
    return 0;
  }

  // Setup RooFormulaVar for cutSpec if it is present
  RooFormula* select = cutSpec ? new RooFormula("select",cutSpec,*get()) : 0 ;

  // Calculate requested moment
  Double_t xysum(0),xsum(0),ysum(0),x2sum(0),y2sum(0);
  const RooArgSet* vars ;
  for(Int_t index= 0; index < numEntries(); index++) {
    vars = get(index) ;
    if (select && select->eval()==0) continue ;
    if (cutRange && vars->allInRange(cutRange)) continue ;

    xysum += weight()*xdata->getVal()*ydata->getVal() ;
    xsum += weight()*xdata->getVal() ;
    ysum += weight()*ydata->getVal() ;
    if (corr) {
      x2sum += weight()*xdata->getVal()*xdata->getVal() ;
      y2sum += weight()*ydata->getVal()*ydata->getVal() ;
    }
  }

  // Normalize entries
  xysum/=sumEntries(cutSpec, cutRange) ;
  xsum/=sumEntries(cutSpec, cutRange) ;
  ysum/=sumEntries(cutSpec, cutRange) ;
  if (corr) {
    x2sum/=sumEntries(cutSpec, cutRange) ;
    y2sum/=sumEntries(cutSpec, cutRange) ;
  }

  // Cleanup
  if (select) delete select ;

  // Return covariance or correlation as requested
  if (corr) {
    return (xysum-xsum*ysum)/(sqrt(x2sum-(xsum*xsum))*sqrt(y2sum-(ysum*ysum))) ;
  } else {
    return (xysum-xsum*ysum);
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Return covariance matrix from data for given list of observables

TMatrixDSym* RooAbsData::corrcovMatrix(const RooArgList& vars, const char* cutSpec, const char* cutRange, Bool_t corr) const
{
  RooArgList varList ;
  TIterator* iter = vars.createIterator() ;
  RooRealVar* var ;
  while((var=(RooRealVar*)iter->Next())) {
    RooRealVar* datavar = dataRealVar("covarianceMatrix",*var) ;
    if (!datavar) {
      delete iter ;
      return 0 ;
    }
    varList.add(*datavar) ;
  }
  delete iter ;


  // Check if dataset is not empty
  if(sumEntries(cutSpec, cutRange) == 0.) {
    coutW(InputArguments) << "RooDataSet::covariance(" << GetName() << ") WARNING: empty dataset, returning zero" << endl ;
    return 0;
  }

  // Setup RooFormulaVar for cutSpec if it is present
  RooFormula* select = cutSpec ? new RooFormula("select",cutSpec,*get()) : 0 ;

  iter = varList.createIterator() ;
  TIterator* iter2 = varList.createIterator() ;

  TMatrixDSym xysum(varList.getSize()) ;
  vector<double> xsum(varList.getSize()) ;
  vector<double> x2sum(varList.getSize()) ;

  // Calculate <x_i> and <x_i y_j>
  for(Int_t index= 0; index < numEntries(); index++) {
    const RooArgSet* dvars = get(index) ;
    if (select && select->eval()==0) continue ;
    if (cutRange && dvars->allInRange(cutRange)) continue ;

    RooRealVar* varx, *vary ;
    iter->Reset() ;
    Int_t ix=0,iy=0 ;
    while((varx=(RooRealVar*)iter->Next())) {
      xsum[ix] += weight()*varx->getVal() ;
      if (corr) {
   x2sum[ix] += weight()*varx->getVal()*varx->getVal() ;
      }

      *iter2=*iter ; iy=ix ;
      vary=varx ;
      while(vary) {
   xysum(ix,iy) += weight()*varx->getVal()*vary->getVal() ;
   xysum(iy,ix) = xysum(ix,iy) ;
   iy++ ;
   vary=(RooRealVar*)iter2->Next() ;
      }
      ix++ ;
    }

  }

  // Normalize sums
  for (Int_t ix=0 ; ix<varList.getSize() ; ix++) {
    xsum[ix] /= sumEntries(cutSpec, cutRange) ;
    if (corr) {
      x2sum[ix] /= sumEntries(cutSpec, cutRange) ;
    }
    for (Int_t iy=0 ; iy<varList.getSize() ; iy++) {
      xysum(ix,iy) /= sumEntries(cutSpec, cutRange) ;
    }
  }

  // Calculate covariance matrix
  TMatrixDSym* C = new TMatrixDSym(varList.getSize()) ;
  for (Int_t ix=0 ; ix<varList.getSize() ; ix++) {
    for (Int_t iy=0 ; iy<varList.getSize() ; iy++) {
      (*C)(ix,iy) = xysum(ix,iy)-xsum[ix]*xsum[iy] ;
      if (corr) {
   (*C)(ix,iy) /= sqrt((x2sum[ix]-(xsum[ix]*xsum[ix]))*(x2sum[iy]-(xsum[iy]*xsum[iy]))) ;
      }
    }
  }

  if (select) delete select ;
  delete iter ;
  delete iter2 ;

  return C ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a RooRealVar containing the mean of observable 'var' in
/// this dataset.  If cutSpec and/or cutRange are specified the
/// moment is calculated on the subset of the data which pass the C++
/// cut specification expression 'cutSpec' and/or are inside the
/// range named 'cutRange'

RooRealVar* RooAbsData::meanVar(const RooRealVar &var, const char* cutSpec, const char* cutRange) const
{
  // Create a new variable with appropriate strings. The error is calculated as
  // RMS/Sqrt(N) which is generally valid.

  // Create holder variable for mean
  TString name(var.GetName()),title("Mean of ") ;
  name.Append("Mean");
  title.Append(var.GetTitle());
  RooRealVar *meanv= new RooRealVar(name,title,0) ;
  meanv->setConstant(kFALSE) ;

  // Adjust plot label
  TString label("<") ;
  label.Append(var.getPlotLabel());
  label.Append(">");
  meanv->setPlotLabel(label.Data());

  // fill in this variable's value and error
  Double_t meanVal=moment(var,1,0,cutSpec,cutRange) ;
  Double_t N(sumEntries(cutSpec,cutRange)) ;

  Double_t rmsVal= sqrt(moment(var,2,meanVal,cutSpec,cutRange)*N/(N-1));
  meanv->setVal(meanVal) ;
  meanv->setError(N > 0 ? rmsVal/sqrt(N) : 0);

  return meanv;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a RooRealVar containing the RMS of observable 'var' in
/// this dataset.  If cutSpec and/or cutRange are specified the
/// moment is calculated on the subset of the data which pass the C++
/// cut specification expression 'cutSpec' and/or are inside the
/// range named 'cutRange'

RooRealVar* RooAbsData::rmsVar(const RooRealVar &var, const char* cutSpec, const char* cutRange) const
{
  // Create a new variable with appropriate strings. The error is calculated as
  // RMS/(2*Sqrt(N)) which is only valid if the variable has a Gaussian distribution.

  // Create RMS value holder
  TString name(var.GetName()),title("RMS of ") ;
  name.Append("RMS");
  title.Append(var.GetTitle());
  RooRealVar *rms= new RooRealVar(name,title,0) ;
  rms->setConstant(kFALSE) ;

  // Adjust plot label
  TString label(var.getPlotLabel());
  label.Append("_{RMS}");
  rms->setPlotLabel(label);

  // Fill in this variable's value and error
  Double_t meanVal(moment(var,1,0,cutSpec,cutRange)) ;
  Double_t N(sumEntries(cutSpec, cutRange));
  Double_t rmsVal= sqrt(moment(var,2,meanVal,cutSpec,cutRange)*N/(N-1));
  rms->setVal(rmsVal) ;
  rms->setError(rmsVal/sqrt(2*N));

  return rms;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a box with statistics information to the specified frame. By default a box with the
/// event count, mean and rms of the plotted variable is added.
///
/// The following optional named arguments are accepted
/// <table>
/// <tr><td> `What(const char* whatstr)`   <td> Controls what is printed: "N" = count, "M" is mean, "R" is RMS.
/// <tr><td> `Format(const char* optStr)`   <td> \deprecated Classing parameter formatting options, provided for backward compatibility
///
/// <tr><td> `Format(const char* what,...)`   <td> Parameter formatting options.
///   <table>
///   <tr><td> const char* what <td> Controls what is shown:
///     - "N" adds name
///     - "E" adds error
///     - "A" shows asymmetric error
///     - "U" shows unit
///     - "H" hides the value
///   <tr><td> `FixedPrecision(int n)`   <td> Controls precision, set fixed number of digits
///   <tr><td> `AutoPrecision(int n)`   <td> Controls precision. Number of shown digits is calculated from error + n specified additional digits (1 is sensible default)
///   <tr><td> `VerbatimName(Bool_t flag)`   <td> Put variable name in a \\verb+   + clause.
///   </table>
/// <tr><td> `Label(const chat* label)`   <td> Add header label to parameter box
/// <tr><td> `Layout(Double_t xmin, Double_t xmax, Double_t ymax)`   <td> Specify relative position of left,right side of box and top of box. Position of
///     bottom of box is calculated automatically from number lines in box
/// <tr><td> `Cut(const char* expression)`   <td> Apply given cut expression to data when calculating statistics
/// <tr><td> `CutRange(const char* rangeName)`   <td> Only consider events within given range when calculating statistics. Multiple
///     CutRange() argument may be specified to combine ranges.
///
/// </table>

RooPlot* RooAbsData::statOn(RooPlot* frame, const RooCmdArg& arg1, const RooCmdArg& arg2,
             const RooCmdArg& arg3, const RooCmdArg& arg4, const RooCmdArg& arg5,
             const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8)
{
  // Stuff all arguments in a list
  RooLinkedList cmdList;
  cmdList.Add(const_cast<RooCmdArg*>(&arg1)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg2)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg3)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg4)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg5)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg6)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg7)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg8)) ;

  // Select the pdf-specific commands
  RooCmdConfig pc(Form("RooTreeData::statOn(%s)",GetName())) ;
  pc.defineString("what","What",0,"MNR") ;
  pc.defineString("label","Label",0,"") ;
  pc.defineDouble("xmin","Layout",0,0.65) ;
  pc.defineDouble("xmax","Layout",1,0.99) ;
  pc.defineInt("ymaxi","Layout",0,Int_t(0.95*10000)) ;
  pc.defineString("formatStr","Format",0,"NELU") ;
  pc.defineInt("sigDigit","Format",0,2) ;
  pc.defineInt("dummy","FormatArgs",0,0) ;
  pc.defineString("cutRange","CutRange",0,"",kTRUE) ;
  pc.defineString("cutString","CutSpec",0,"") ;
  pc.defineMutex("Format","FormatArgs") ;

  // Process and check varargs
  pc.process(cmdList) ;
  if (!pc.ok(kTRUE)) {
    return frame ;
  }

  const char* label = pc.getString("label") ;
  Double_t xmin = pc.getDouble("xmin") ;
  Double_t xmax = pc.getDouble("xmax") ;
  Double_t ymax = pc.getInt("ymaxi") / 10000. ;
  const char* formatStr = pc.getString("formatStr") ;
  Int_t sigDigit = pc.getInt("sigDigit") ;
  const char* what = pc.getString("what") ;

  const char* cutSpec = pc.getString("cutString",0,kTRUE) ;
  const char* cutRange = pc.getString("cutRange",0,kTRUE) ;

  if (pc.hasProcessed("FormatArgs")) {
    RooCmdArg* formatCmd = static_cast<RooCmdArg*>(cmdList.FindObject("FormatArgs")) ;
    return statOn(frame,what,label,0,0,xmin,xmax,ymax,cutSpec,cutRange,formatCmd) ;
  } else {
    return statOn(frame,what,label,sigDigit,formatStr,xmin,xmax,ymax,cutSpec,cutRange) ;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation back-end of statOn() method with named arguments

RooPlot* RooAbsData::statOn(RooPlot* frame, const char* what, const char *label, Int_t sigDigits,
              Option_t *options, Double_t xmin, Double_t xmax, Double_t ymax,
              const char* cutSpec, const char* cutRange, const RooCmdArg* formatCmd)
{
  Bool_t showLabel= (label != 0 && strlen(label) > 0);

  TString whatStr(what) ;
  whatStr.ToUpper() ;
  Bool_t showN = whatStr.Contains("N") ;
  Bool_t showR = whatStr.Contains("R") ;
  Bool_t showM = whatStr.Contains("M") ;
  Int_t nPar= 0;
  if (showN) nPar++ ;
  if (showR) nPar++ ;
  if (showM) nPar++ ;

  // calculate the box's size
  Double_t dy(0.06), ymin(ymax-nPar*dy);
  if(showLabel) ymin-= dy;

  // create the box and set its options
  TPaveText *box= new TPaveText(xmin,ymax,xmax,ymin,"BRNDC");
  if(!box) return 0;
  box->SetName(Form("%s_statBox",GetName())) ;
  box->SetFillColor(0);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextSize(0.04F);
  box->SetFillStyle(1001);

  // add formatted text for each statistic
  RooRealVar N("N","Number of Events",sumEntries(cutSpec,cutRange));
  N.setPlotLabel("Entries") ;
  RooRealVar *meanv= meanVar(*(RooRealVar*)frame->getPlotVar(),cutSpec,cutRange);
  meanv->setPlotLabel("Mean") ;
  RooRealVar *rms= rmsVar(*(RooRealVar*)frame->getPlotVar(),cutSpec,cutRange);
  rms->setPlotLabel("RMS") ;
  TString *rmsText, *meanText, *NText ;
  if (options) {
    rmsText= rms->format(sigDigits,options);
    meanText= meanv->format(sigDigits,options);
    NText= N.format(sigDigits,options);
  } else {
    rmsText= rms->format(*formatCmd);
    meanText= meanv->format(*formatCmd);
    NText= N.format(*formatCmd);
  }
  if (showR) box->AddText(rmsText->Data());
  if (showM) box->AddText(meanText->Data());
  if (showN) box->AddText(NText->Data());

  // cleanup heap memory
  delete NText;
  delete meanText;
  delete rmsText;
  delete meanv;
  delete rms;

  // add the optional label if specified
  if(showLabel) box->AddText(label);

  frame->addObject(box) ;
  return frame ;
}

////////////////////////////////////////////////////////////////////////////////
/// Loop over columns of our tree data and fill the input histogram. Returns a pointer to the
/// input histogram, or zero in case of an error. The input histogram can be any TH1 subclass, and
/// therefore of arbitrary dimension. Variables are matched with the (x,y,...) dimensions of the input
/// histogram according to the order in which they appear in the input plotVars list.

TH1 *RooAbsData::fillHistogram(TH1 *hist, const RooArgList &plotVars, const char *cuts, const char* cutRange) const
{
  // Do we have a valid histogram to use?
  if(0 == hist) {
    coutE(InputArguments) << ClassName() << "::" << GetName() << ":fillHistogram: no valid histogram to fill" << endl;
    return 0;
  }

  // Check that the number of plotVars matches the input histogram's dimension
  Int_t hdim= hist->GetDimension();
  if(hdim != plotVars.getSize()) {
    coutE(InputArguments) << ClassName() << "::" << GetName() << ":fillHistogram: plotVars has the wrong dimension" << endl;
    return 0;
  }

  // Check that the plot variables are all actually RooAbsReal's and print a warning if we do not
  // explicitly depend on one of them. Clone any variables that we do not contain directly and
  // redirect them to use our event data.
  RooArgSet plotClones,localVars;
  for(Int_t index= 0; index < plotVars.getSize(); index++) {
    const RooAbsArg *var= plotVars.at(index);
    const RooAbsReal *realVar= dynamic_cast<const RooAbsReal*>(var);
    if(0 == realVar) {
      coutE(InputArguments) << ClassName() << "::" << GetName() << ":fillHistogram: cannot plot variable \"" << var->GetName()
      << "\" of type " << var->ClassName() << endl;
      return 0;
    }
    RooAbsArg *found= _vars.find(realVar->GetName());
    if(!found) {
      RooAbsArg *clone= plotClones.addClone(*realVar,kTRUE); // do not complain about duplicates
      assert(0 != clone);
      if(!clone->dependsOn(_vars)) {
        coutE(InputArguments) << ClassName() << "::" << GetName()
            << ":fillHistogram: Data does not contain the variable '" << realVar->GetName() << "'." << endl;
        return nullptr;
      }
      else {
   clone->recursiveRedirectServers(_vars);
      }
      localVars.add(*clone);
    }
    else {
      localVars.add(*found);
    }
  }

  // Create selection formula if selection cuts are specified
  std::unique_ptr<RooFormula> select;
  if (cuts != nullptr && strlen(cuts) > 0) {
    select.reset(new RooFormula(cuts, cuts, _vars, false));
    if (!select || !select->ok()) {
      coutE(InputArguments) << ClassName() << "::" << GetName() << ":fillHistogram: invalid cuts \"" << cuts << "\"" << endl;
      return 0 ;
    }
  }

  // Lookup each of the variables we are binning in our tree variables
  const RooAbsReal *xvar = 0;
  const RooAbsReal *yvar = 0;
  const RooAbsReal *zvar = 0;
  switch(hdim) {
  case 3:
    zvar= dynamic_cast<RooAbsReal*>(localVars.find(plotVars.at(2)->GetName()));
    assert(0 != zvar);
    // fall through to next case...
  case 2:
    yvar= dynamic_cast<RooAbsReal*>(localVars.find(plotVars.at(1)->GetName()));
    assert(0 != yvar);
    // fall through to next case...
  case 1:
    xvar= dynamic_cast<RooAbsReal*>(localVars.find(plotVars.at(0)->GetName()));
    assert(0 != xvar);
    break;
  default:
    coutE(InputArguments) << ClassName() << "::" << GetName() << ":fillHistogram: cannot fill histogram with "
    << hdim << " dimensions" << endl;
    break;
  }

  // Parse cutRange specification
  const auto cutVec = RooHelpers::tokenise(cutRange ? cutRange : "", ",", /*returnEmptyToken =*/ false);

  // Loop over events and fill the histogram
  if (hist->GetSumw2()->fN==0) {
    hist->Sumw2() ;
  }
  Int_t nevent= numEntries() ; //(Int_t)_tree->GetEntries();
  for(Int_t i=0; i < nevent; ++i) {

    //Int_t entryNumber= _tree->GetEntryNumber(i);
    //if (entryNumber<0) break;
    get(i);

    // Apply expression based selection criteria
    if (select && select->eval()==0) {
      continue ;
    }


    // Apply range based selection criteria
    Bool_t selectByRange = kTRUE ;
    if (cutRange) {
      for (const auto arg : _vars) {
        Bool_t selectThisArg = kFALSE ;
        UInt_t icut ;
        for (icut=0 ; icut<cutVec.size() ; icut++) {
          if (arg->inRange(cutVec[icut].c_str())) {
            selectThisArg = kTRUE ;
            break ;
          }
        }
        if (!selectThisArg) {
          selectByRange = kFALSE ;
          break ;
        }
      }
    }

    if (!selectByRange) {
      // Go to next event in loop over events
      continue ;
    }

    Int_t bin(0);
    switch(hdim) {
    case 1:
      bin= hist->FindBin(xvar->getVal());
      hist->Fill(xvar->getVal(),weight()) ;
      break;
    case 2:
      bin= hist->FindBin(xvar->getVal(),yvar->getVal());
      static_cast<TH2*>(hist)->Fill(xvar->getVal(),yvar->getVal(),weight()) ;
      break;
    case 3:
      bin= hist->FindBin(xvar->getVal(),yvar->getVal(),zvar->getVal());
      static_cast<TH3*>(hist)->Fill(xvar->getVal(),yvar->getVal(),zvar->getVal(),weight()) ;
      break;
    default:
      assert(hdim < 3);
      break;
    }


    Double_t error2 = TMath::Power(hist->GetBinError(bin),2)-TMath::Power(weight(),2)  ;
    Double_t we = weightError(RooAbsData::SumW2) ;
    if (we==0) we = weight() ;
    error2 += TMath::Power(we,2) ;


//     Double_t we = weightError(RooAbsData::SumW2) ;
//     Double_t error2(0) ;
//     if (we==0) {
//       we = weight() ; //sqrt(weight()) ;
//       error2 = TMath::Power(hist->GetBinError(bin),2)-TMath::Power(weight(),2) + TMath::Power(we,2) ;
//     } else {
//       error2 = TMath::Power(hist->GetBinError(bin),2)-TMath::Power(weight(),2) + TMath::Power(we,2) ;
//     }
    //hist->AddBinContent(bin,weight());
    hist->SetBinError(bin,sqrt(error2)) ;

    //cout << "RooTreeData::fillHistogram() bin = " << bin << " weight() = " << weight() << " we = " << we << endl ;

  }

  return hist;
}

////////////////////////////////////////////////////////////////////////////////
/// Split dataset into subsets based on states of given splitCat in this dataset.
/// A TList of RooDataSets is returned in which each RooDataSet is named
/// after the state name of splitCat of which it contains the dataset subset.
/// The observables splitCat itself is no longer present in the sub datasets.
/// If createEmptyDataSets is kFALSE (default) this method only creates datasets for states
/// which have at least one entry The caller takes ownership of the returned list and its contents

TList* RooAbsData::split(const RooAbsCategory& splitCat, Bool_t createEmptyDataSets) const
{
  // Sanity check
  if (!splitCat.dependsOn(*get())) {
    coutE(InputArguments) << "RooTreeData::split(" << GetName() << ") ERROR category " << splitCat.GetName()
    << " doesn't depend on any variable in this dataset" << endl ;
    return 0 ;
  }

  // Clone splitting category and attach to self
  RooAbsCategory* cloneCat =0;
  RooArgSet* cloneSet = 0;
  if (splitCat.isDerived()) {
    cloneSet = (RooArgSet*) RooArgSet(splitCat).snapshot(kTRUE) ;
    if (!cloneSet) {
      coutE(InputArguments) << "RooTreeData::split(" << GetName() << ") Couldn't deep-clone splitting category, abort." << endl ;
      return 0 ;
    }
    cloneCat = (RooAbsCategory*) cloneSet->find(splitCat.GetName()) ;
    cloneCat->attachDataSet(*this) ;
  } else {
    cloneCat = dynamic_cast<RooAbsCategory*>(get()->find(splitCat.GetName())) ;
    if (!cloneCat) {
      coutE(InputArguments) << "RooTreeData::split(" << GetName() << ") ERROR category " << splitCat.GetName()
      << " is fundamental and does not appear in this dataset" << endl ;
      return 0 ;
    }
  }

  // Split a dataset in a series of subsets, each corresponding
  // to a state of splitCat
  TList* dsetList = new TList ;

  // Construct set of variables to be included in split sets = full set - split category
  RooArgSet subsetVars(*get()) ;
  if (splitCat.isDerived()) {
    RooArgSet* vars = splitCat.getVariables() ;
    subsetVars.remove(*vars,kTRUE,kTRUE) ;
    delete vars ;
  } else {
    subsetVars.remove(splitCat,kTRUE,kTRUE) ;
  }

  // Add weight variable explicitly if dataset has weights, but no top-level weight
  // variable exists (can happen with composite datastores)
  Bool_t addWV(kFALSE) ;
  RooRealVar newweight("weight","weight",-1e9,1e9) ;
  if (isWeighted() && !IsA()->InheritsFrom(RooDataHist::Class())) {
    subsetVars.add(newweight) ;
    addWV = kTRUE ;
  }

  // If createEmptyDataSets is true, prepopulate with empty sets corresponding to all states
  if (createEmptyDataSets) {
    for (const auto& nameIdx : *cloneCat) {
      RooAbsData* subset = emptyClone(nameIdx.first.c_str(), nameIdx.first.c_str(), &subsetVars,(addWV?"weight":0)) ;
      dsetList->Add((RooAbsArg*)subset) ;
    }
  }


  // Loop over dataset and copy event to matching subset
  const bool propWeightSquared = isWeighted();
  for (Int_t i = 0; i < numEntries(); ++i) {
    const RooArgSet* row =  get(i);
    RooAbsData* subset = (RooAbsData*) dsetList->FindObject(cloneCat->getCurrentLabel());
    if (!subset) {
      subset = emptyClone(cloneCat->getCurrentLabel(),cloneCat->getCurrentLabel(),&subsetVars,(addWV?"weight":0));
      dsetList->Add((RooAbsArg*)subset);
    }
    if (!propWeightSquared) {
   subset->add(*row, weight());
    } else {
   subset->add(*row, weight(), weightSquared());
    }
  }

  delete cloneSet;
  return dsetList;
}

////////////////////////////////////////////////////////////////////////////////
/// Plot dataset on specified frame.
///
/// By default:
/// - An unbinned dataset will use the default binning of the target frame.
/// - A binned dataset will retain its intrinsic binning.
///
/// The following optional named arguments can be used to modify the behaviour:
///
/// <table>
/// <tr><th> <th> Data representation options
/// <tr><td> `Asymmetry(const RooCategory& c)`   <td> Show the asymmetry of the data in given two-state category [F(+)-F(-)] / [F(+)+F(-)].
///     Category must have two states with indices -1 and +1 or three states with indices -1,0 and +1.
/// <tr><td> `Efficiency(const RooCategory& c)`   <td> Show the efficiency F(acc)/[F(acc)+F(rej)]. Category must have two states with indices 0 and 1
/// <tr><td> `DataError(RooAbsData::EType)`   <td> Select the type of error drawn:
///    - `Auto(default)` results in Poisson for unweighted data and SumW2 for weighted data
///    - `Poisson` draws asymmetric Poisson confidence intervals.
///    - `SumW2` draws symmetric sum-of-weights error ( \f$ \left( \sum w \right)^2 / \sum\left(w^2\right) \f$ )
///    - `None` draws no error bars
/// <tr><td> `Binning(int nbins, double xlo, double xhi)`   <td> Use specified binning to draw dataset
/// <tr><td> `Binning(const RooAbsBinning&)`   <td>  Use specified binning to draw dataset
/// <tr><td> `Binning(const char* name)`   <td>  Use binning with specified name to draw dataset
/// <tr><td> `RefreshNorm(Bool_t flag)`   <td> Force refreshing for PDF normalization information in frame.
///     If set, any subsequent PDF will normalize to this dataset, even if it is
///     not the first one added to the frame. By default only the 1st dataset
///     added to a frame will update the normalization information
/// <tr><td> `Rescale(Double_t f)`   <td> Rescale drawn histogram by given factor.
/// <tr><td> `Cut(const char*)`      <td> Only plot entries that pass the given cut.
///                                       Apart from cutting in continuous variables `Cut("x>5")`, this can also be used to plot a specific
///                                       category state. Use something like `Cut("myCategory == myCategory::stateA")`, where
///                                       `myCategory` resolves to the state number for a given entry and
///                                       `myCategory::stateA` resolves to the state number of the state named "stateA".
///
/// <tr><td> `CutRange(const char*)` <td> Only plot data from given range. Separate multiple ranges with ",".
/// \note This often requires passing the normalisation when plotting the PDF because RooFit does not save
/// how many events were being plotted (it will only work for cutting slices out of uniformly distributed variables).
/// ```
///  data->plotOn(frame01, CutRange("SB1"));
///  const double nData = data->sumEntries("", "SB1");
///  // Make clear that the target normalisation is nData. The enumerator NumEvent
///  // is needed to switch between relative and absolute scaling.
///  model.plotOn(frame01, Normalization(nData, RooAbsReal::NumEvent),
///    ProjectionRange("SB1"));
/// ```
///
/// <tr><th> <th> Histogram drawing options
/// <tr><td> `DrawOption(const char* opt)`   <td> Select ROOT draw option for resulting TGraph object
/// <tr><td> `LineStyle(Int_t style)`   <td> Select line style by ROOT line style code, default is solid
/// <tr><td> `LineColor(Int_t color)`   <td> Select line color by ROOT color code, default is black
/// <tr><td> `LineWidth(Int_t width)`   <td> Select line with in pixels, default is 3
/// <tr><td> `MarkerStyle(Int_t style)`   <td> Select the ROOT marker style, default is 21
/// <tr><td> `MarkerColor(Int_t color)`   <td> Select the ROOT marker color, default is black
/// <tr><td> `MarkerSize(Double_t size)`   <td> Select the ROOT marker size
/// <tr><td> `FillStyle(Int_t style)`   <td> Select fill style, default is filled.
/// <tr><td> `FillColor(Int_t color)`   <td> Select fill color by ROOT color code
/// <tr><td> `XErrorSize(Double_t frac)`   <td> Select size of X error bar as fraction of the bin width, default is 1
///
///
/// <tr><th> <th> Misc. other options
/// <tr><td> `Name(const chat* name)`   <td> Give curve specified name in frame. Useful if curve is to be referenced later
/// <tr><td> `Invisible()`   <td> Add curve to frame, but do not display. Useful in combination AddTo()
/// <tr><td> `AddTo(const char* name, double_t wgtSelf, double_t wgtOther)`   <td> Add constructed histogram to already existing histogram with given name and relative weight factors
/// </table>

RooPlot* RooAbsData::plotOn(RooPlot* frame, const RooLinkedList& argList) const
{
  // New experimental plotOn() with varargs...

  // Define configuration for this method
  RooCmdConfig pc(Form("RooAbsData::plotOn(%s)",GetName())) ;
  pc.defineString("drawOption","DrawOption",0,"P") ;
  pc.defineString("cutRange","CutRange",0,"",kTRUE) ;
  pc.defineString("cutString","CutSpec",0,"") ;
  pc.defineString("histName","Name",0,"") ;
  pc.defineObject("cutVar","CutVar",0) ;
  pc.defineObject("binning","Binning",0) ;
  pc.defineString("binningName","BinningName",0,"") ;
  pc.defineInt("nbins","BinningSpec",0,100) ;
  pc.defineDouble("xlo","BinningSpec",0,0) ;
  pc.defineDouble("xhi","BinningSpec",1,1) ;
  pc.defineObject("asymCat","Asymmetry",0) ;
  pc.defineObject("effCat","Efficiency",0) ;
  pc.defineInt("lineColor","LineColor",0,-999) ;
  pc.defineInt("lineStyle","LineStyle",0,-999) ;
  pc.defineInt("lineWidth","LineWidth",0,-999) ;
  pc.defineInt("markerColor","MarkerColor",0,-999) ;
  pc.defineInt("markerStyle","MarkerStyle",0,-999) ;
  pc.defineDouble("markerSize","MarkerSize",0,-999) ;
  pc.defineInt("fillColor","FillColor",0,-999) ;
  pc.defineInt("fillStyle","FillStyle",0,-999) ;
  pc.defineInt("errorType","DataError",0,(Int_t)RooAbsData::Auto) ;
  pc.defineInt("histInvisible","Invisible",0,0) ;
  pc.defineInt("refreshFrameNorm","RefreshNorm",0,1) ;
  pc.defineString("addToHistName","AddTo",0,"") ;
  pc.defineDouble("addToWgtSelf","AddTo",0,1.) ;
  pc.defineDouble("addToWgtOther","AddTo",1,1.) ;
  pc.defineDouble("xErrorSize","XErrorSize",0,1.) ;
  pc.defineDouble("scaleFactor","Rescale",0,1.) ;
  pc.defineMutex("DataError","Asymmetry","Efficiency") ;
  pc.defineMutex("Binning","BinningName","BinningSpec") ;

  // Process & check varargs
  pc.process(argList) ;
  if (!pc.ok(kTRUE)) {
    return frame ;
  }

  PlotOpt o ;

  // Extract values from named arguments
  o.drawOptions = pc.getString("drawOption") ;
  o.cuts = pc.getString("cutString") ;
  if (pc.hasProcessed("Binning")) {
    o.bins = (RooAbsBinning*) pc.getObject("binning") ;
  } else if (pc.hasProcessed("BinningName")) {
    o.bins = &frame->getPlotVar()->getBinning(pc.getString("binningName")) ;
  } else if (pc.hasProcessed("BinningSpec")) {
    Double_t xlo = pc.getDouble("xlo") ;
    Double_t xhi = pc.getDouble("xhi") ;
    o.bins = new RooUniformBinning((xlo==xhi)?frame->getPlotVar()->getMin():xlo,
               (xlo==xhi)?frame->getPlotVar()->getMax():xhi,pc.getInt("nbins")) ;
  }
  const RooAbsCategoryLValue* asymCat = (const RooAbsCategoryLValue*) pc.getObject("asymCat") ;
  const RooAbsCategoryLValue* effCat = (const RooAbsCategoryLValue*) pc.getObject("effCat") ;
  o.etype = (RooAbsData::ErrorType) pc.getInt("errorType") ;
  o.histInvisible = pc.getInt("histInvisible") ;
  o.xErrorSize = pc.getDouble("xErrorSize") ;
  o.cutRange = pc.getString("cutRange",0,kTRUE) ;
  o.histName = pc.getString("histName",0,kTRUE) ;
  o.addToHistName = pc.getString("addToHistName",0,kTRUE) ;
  o.addToWgtSelf = pc.getDouble("addToWgtSelf") ;
  o.addToWgtOther = pc.getDouble("addToWgtOther") ;
  o.refreshFrameNorm = pc.getInt("refreshFrameNorm") ;
  o.scaleFactor = pc.getDouble("scaleFactor") ;

  // Map auto error type to actual type
  if (o.etype == Auto) {
    o.etype = isNonPoissonWeighted() ? SumW2 : Poisson ;
    if (o.etype == SumW2) {
      coutI(InputArguments) << "RooAbsData::plotOn(" << GetName()
             << ") INFO: dataset has non-integer weights, auto-selecting SumW2 errors instead of Poisson errors" << endl ;
    }
  }

  if (o.addToHistName && !frame->findObject(o.addToHistName,RooHist::Class())) {
    coutE(InputArguments) << "RooAbsData::plotOn(" << GetName() << ") cannot find existing histogram " << o.addToHistName
           << " to add to in RooPlot" << endl ;
    return frame ;
  }

  RooPlot* ret ;
  if (!asymCat && !effCat) {
    ret = plotOn(frame,o) ;
  } else if (asymCat) {
    ret = plotAsymOn(frame,*asymCat,o) ;
  } else {
    ret = plotEffOn(frame,*effCat,o) ;
  }

  Int_t lineColor   = pc.getInt("lineColor") ;
  Int_t lineStyle   = pc.getInt("lineStyle") ;
  Int_t lineWidth   = pc.getInt("lineWidth") ;
  Int_t markerColor = pc.getInt("markerColor") ;
  Int_t markerStyle = pc.getInt("markerStyle") ;
  Size_t markerSize  = pc.getDouble("markerSize") ;
  Int_t fillColor = pc.getInt("fillColor") ;
  Int_t fillStyle = pc.getInt("fillStyle") ;
  if (lineColor!=-999) ret->getAttLine()->SetLineColor(lineColor) ;
  if (lineStyle!=-999) ret->getAttLine()->SetLineStyle(lineStyle) ;
  if (lineWidth!=-999) ret->getAttLine()->SetLineWidth(lineWidth) ;
  if (markerColor!=-999) ret->getAttMarker()->SetMarkerColor(markerColor) ;
  if (markerStyle!=-999) ret->getAttMarker()->SetMarkerStyle(markerStyle) ;
  if (markerSize!=-999) ret->getAttMarker()->SetMarkerSize(markerSize) ;
  if (fillColor!=-999) ret->getAttFill()->SetFillColor(fillColor) ;
  if (fillStyle!=-999) ret->getAttFill()->SetFillStyle(fillStyle) ;

  if (pc.hasProcessed("BinningSpec")) {
    delete o.bins ;
  }

  return ret ;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and fill a histogram of the frame's variable and append it to the frame.
/// The frame variable must be one of the data sets dimensions.
///
/// The plot range and the number of plot bins is determined by the parameters
/// of the plot variable of the frame (RooAbsReal::setPlotRange(), RooAbsReal::setPlotBins()).
///
/// The optional cut string expression can be used to select the events to be plotted.
/// The cut specification may refer to any variable contained in the data set.
///
/// The drawOptions are passed to the TH1::Draw() method.
/// \see RooAbsData::plotOn(RooPlot*,const RooLinkedList&) const
RooPlot *RooAbsData::plotOn(RooPlot *frame, PlotOpt o) const
{
  if(0 == frame) {
    coutE(Plotting) << ClassName() << "::" << GetName() << ":plotOn: frame is null" << endl;
    return 0;
  }
  RooAbsRealLValue *var= (RooAbsRealLValue*) frame->getPlotVar();
  if(0 == var) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotOn: frame does not specify a plot variable" << endl;
    return 0;
  }

  // create and fill a temporary histogram of this variable
  TString histName(GetName());
  histName.Append("_plot");
  std::unique_ptr<TH1> hist;
  if (o.bins) {
    hist.reset( var->createHistogram(histName.Data(), RooFit::AxisLabel("Events"), RooFit::Binning(*o.bins)) );
  } else if (!frame->getPlotVar()->getBinning().isUniform()) {
    hist.reset( var->createHistogram(histName.Data(), RooFit::AxisLabel("Events"),
        RooFit::Binning(frame->getPlotVar()->getBinning())) );
  } else {
    hist.reset( var->createHistogram(histName.Data(), "Events",
        frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax(), frame->GetNbinsX()) );
  }

  // Keep track of sum-of-weights error
  hist->Sumw2() ;

  if(0 == fillHistogram(hist.get(), RooArgList(*var),o.cuts,o.cutRange)) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotOn: fillHistogram() failed" << endl;
    return 0;
  }

  // If frame has no predefined bin width (event density) it will be adjusted to
  // our histograms bin width so we should force that bin width here
  Double_t nomBinWidth ;
  if (frame->getFitRangeNEvt()==0 && o.bins) {
    nomBinWidth = o.bins->averageBinWidth() ;
  } else {
    nomBinWidth = o.bins ? frame->getFitRangeBinW() : 0 ;
  }

  // convert this histogram to a RooHist object on the heap
  RooHist *graph= new RooHist(*hist,nomBinWidth,1,o.etype,o.xErrorSize,o.correctForBinWidth,o.scaleFactor);
  if(0 == graph) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotOn: unable to create a RooHist object" << endl;
    return 0;
  }

  // If the dataset variable has a wide range than the plot variable,
  // calculate the number of entries in the dataset in the plot variable fit range
  RooAbsRealLValue* dataVar = (RooAbsRealLValue*) _vars.find(var->GetName()) ;
  Double_t nEnt(sumEntries()) ;
  if (dataVar->getMin()<var->getMin() || dataVar->getMax()>var->getMax()) {
    RooAbsData* tmp = ((RooAbsData*)this)->reduce(*var) ;
    nEnt = tmp->sumEntries() ;
    delete tmp ;
  }

  // Store the number of entries before the cut, if any was made
  if ((o.cuts && strlen(o.cuts)) || o.cutRange) {
    coutI(Plotting) << "RooTreeData::plotOn: plotting " << hist->GetSumOfWeights() << " events out of " << nEnt << " total events" << endl ;
    graph->setRawEntries(nEnt) ;
  }

  // Add self to other hist if requested
  if (o.addToHistName) {
    RooHist* otherGraph = static_cast<RooHist*>(frame->findObject(o.addToHistName,RooHist::Class())) ;

    if (!graph->hasIdenticalBinning(*otherGraph)) {
      coutE(Plotting) << "RooTreeData::plotOn: ERROR Histogram to be added to, '" << o.addToHistName << "',has different binning" << endl ;
      delete graph ;
      return frame ;
    }

    RooHist* sumGraph = new RooHist(*graph,*otherGraph,o.addToWgtSelf,o.addToWgtOther,o.etype) ;
    delete graph ;
    graph = sumGraph ;
  }

  // Rename graph if requested
  if (o.histName) {
    graph->SetName(o.histName) ;
  } else {
    TString hname(Form("h_%s",GetName())) ;
    if (o.cutRange && strlen(o.cutRange)>0) {
      hname.Append(Form("_CutRange[%s]",o.cutRange)) ;
    }
    if (o.cuts && strlen(o.cuts)>0) {
      hname.Append(Form("_Cut[%s]",o.cuts)) ;
    }
    graph->SetName(hname.Data()) ;
  }

  // initialize the frame's normalization setup, if necessary
  frame->updateNormVars(_vars);


  // add the RooHist to the specified plot
  frame->addPlotable(graph,o.drawOptions,o.histInvisible,o.refreshFrameNorm);

  return frame;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and fill a histogram with the asymmetry N[+] - N[-] / ( N[+] + N[-] ),
/// where N(+/-) is the number of data points with asymCat=+1 and asymCat=-1
/// as function of the frames variable. The asymmetry category 'asymCat' must
/// have exactly 2 (or 3) states defined with index values +1,-1 (and 0)
///
/// The plot range and the number of plot bins is determined by the parameters
/// of the plot variable of the frame (RooAbsReal::setPlotRange(), RooAbsReal::setPlotBins())
///
/// The optional cut string expression can be used to select the events to be plotted.
/// The cut specification may refer to any variable contained in the data set
///
/// The drawOptions are passed to the TH1::Draw() method

RooPlot* RooAbsData::plotAsymOn(RooPlot* frame, const RooAbsCategoryLValue& asymCat, PlotOpt o) const
{
  if(0 == frame) {
    coutE(Plotting) << ClassName() << "::" << GetName() << ":plotAsymOn: frame is null" << endl;
    return 0;
  }
  RooAbsRealLValue *var= (RooAbsRealLValue*) frame->getPlotVar();
  if(0 == var) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotAsymOn: frame does not specify a plot variable" << endl;
    return 0;
  }

  // create and fill temporary histograms of this variable for each state
  TString hist1Name(GetName()),hist2Name(GetName());
  hist1Name.Append("_plot1");
  std::unique_ptr<TH1> hist1, hist2;
  hist2Name.Append("_plot2");

  if (o.bins) {
    hist1.reset( var->createHistogram(hist1Name.Data(), "Events", *o.bins) );
    hist2.reset( var->createHistogram(hist2Name.Data(), "Events", *o.bins) );
  } else {
    hist1.reset( var->createHistogram(hist1Name.Data(), "Events",
            frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax(),
            frame->GetNbinsX()) );
    hist2.reset( var->createHistogram(hist2Name.Data(), "Events",
            frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax(),
            frame->GetNbinsX()) );
  }

  assert(hist1 && hist2);

  TString cuts1,cuts2 ;
  if (o.cuts && strlen(o.cuts)) {
    cuts1 = Form("(%s)&&(%s>0)",o.cuts,asymCat.GetName());
    cuts2 = Form("(%s)&&(%s<0)",o.cuts,asymCat.GetName());
  } else {
    cuts1 = Form("(%s>0)",asymCat.GetName());
    cuts2 = Form("(%s<0)",asymCat.GetName());
  }

  if(! fillHistogram(hist1.get(), RooArgList(*var),cuts1.Data(),o.cutRange) ||
     ! fillHistogram(hist2.get(), RooArgList(*var),cuts2.Data(),o.cutRange)) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotAsymOn: createHistogram() failed" << endl;
    return 0;
  }

  // convert this histogram to a RooHist object on the heap
  RooHist *graph= new RooHist(*hist1,*hist2,0,1,o.etype,o.xErrorSize,kFALSE,o.scaleFactor);
  graph->setYAxisLabel(Form("Asymmetry in %s",asymCat.GetName())) ;

  // initialize the frame's normalization setup, if necessary
  frame->updateNormVars(_vars);

  // Rename graph if requested
  if (o.histName) {
    graph->SetName(o.histName) ;
  } else {
    TString hname(Form("h_%s_Asym[%s]",GetName(),asymCat.GetName())) ;
    if (o.cutRange && strlen(o.cutRange)>0) {
      hname.Append(Form("_CutRange[%s]",o.cutRange)) ;
    }
    if (o.cuts && strlen(o.cuts)>0) {
      hname.Append(Form("_Cut[%s]",o.cuts)) ;
    }
    graph->SetName(hname.Data()) ;
  }

  // add the RooHist to the specified plot
  frame->addPlotable(graph,o.drawOptions,o.histInvisible,o.refreshFrameNorm);

  return frame;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and fill a histogram with the efficiency N[1] / ( N[1] + N[0] ),
/// where N(1/0) is the number of data points with effCat=1 and effCat=0
/// as function of the frames variable. The efficiency category 'effCat' must
/// have exactly 2 +1 and 0.
///
/// The plot range and the number of plot bins is determined by the parameters
/// of the plot variable of the frame (RooAbsReal::setPlotRange(), RooAbsReal::setPlotBins())
///
/// The optional cut string expression can be used to select the events to be plotted.
/// The cut specification may refer to any variable contained in the data set
///
/// The drawOptions are passed to the TH1::Draw() method

RooPlot* RooAbsData::plotEffOn(RooPlot* frame, const RooAbsCategoryLValue& effCat, PlotOpt o) const
{
  if(0 == frame) {
    coutE(Plotting) << ClassName() << "::" << GetName() << ":plotEffOn: frame is null" << endl;
    return 0;
  }
  RooAbsRealLValue *var= (RooAbsRealLValue*) frame->getPlotVar();
  if(0 == var) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotEffOn: frame does not specify a plot variable" << endl;
    return 0;
  }

  // create and fill temporary histograms of this variable for each state
  TString hist1Name(GetName()),hist2Name(GetName());
  hist1Name.Append("_plot1");
  std::unique_ptr<TH1> hist1, hist2;
  hist2Name.Append("_plot2");

  if (o.bins) {
    hist1.reset( var->createHistogram(hist1Name.Data(), "Events", *o.bins) );
    hist2.reset( var->createHistogram(hist2Name.Data(), "Events", *o.bins) );
  } else {
    hist1.reset( var->createHistogram(hist1Name.Data(), "Events",
            frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax(),
            frame->GetNbinsX()) );
    hist2.reset( var->createHistogram(hist2Name.Data(), "Events",
            frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax(),
            frame->GetNbinsX()) );
  }

  assert(hist1 && hist2);

  TString cuts1,cuts2 ;
  if (o.cuts && strlen(o.cuts)) {
    cuts1 = Form("(%s)&&(%s==1)",o.cuts,effCat.GetName());
    cuts2 = Form("(%s)&&(%s==0)",o.cuts,effCat.GetName());
  } else {
    cuts1 = Form("(%s==1)",effCat.GetName());
    cuts2 = Form("(%s==0)",effCat.GetName());
  }

  if(! fillHistogram(hist1.get(), RooArgList(*var),cuts1.Data(),o.cutRange) ||
     ! fillHistogram(hist2.get(), RooArgList(*var),cuts2.Data(),o.cutRange)) {
    coutE(Plotting) << ClassName() << "::" << GetName()
    << ":plotEffOn: createHistogram() failed" << endl;
    return 0;
  }

  // convert this histogram to a RooHist object on the heap
  RooHist *graph= new RooHist(*hist1,*hist2,0,1,o.etype,o.xErrorSize,kTRUE);
  graph->setYAxisLabel(Form("Efficiency of %s=%s", effCat.GetName(), effCat.lookupName(1).c_str()));

  // initialize the frame's normalization setup, if necessary
  frame->updateNormVars(_vars);

  // Rename graph if requested
  if (o.histName) {
    graph->SetName(o.histName) ;
  } else {
    TString hname(Form("h_%s_Eff[%s]",GetName(),effCat.GetName())) ;
    if (o.cutRange && strlen(o.cutRange)>0) {
      hname.Append(Form("_CutRange[%s]",o.cutRange)) ;
    }
    if (o.cuts && strlen(o.cuts)>0) {
      hname.Append(Form("_Cut[%s]",o.cuts)) ;
    }
    graph->SetName(hname.Data()) ;
  }

  // add the RooHist to the specified plot
  frame->addPlotable(graph,o.drawOptions,o.histInvisible,o.refreshFrameNorm);

  return frame;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and fill a 1-dimensional table for given category column
/// This functions is the equivalent of plotOn() for category dimensions.
///
/// The optional cut string expression can be used to select the events to be tabulated
/// The cut specification may refer to any variable contained in the data set
///
/// The option string is currently not used

Roo1DTable* RooAbsData::table(const RooAbsCategory& cat, const char* cuts, const char* /*opts*/) const
{
  // First see if var is in data set
  RooAbsCategory* tableVar = (RooAbsCategory*) _vars.find(cat.GetName()) ;
  RooArgSet *tableSet = 0;
  Bool_t ownPlotVar(kFALSE) ;
  if (!tableVar) {
    if (!cat.dependsOn(_vars)) {
      coutE(Plotting) << "RooTreeData::Table(" << GetName() << "): Argument " << cat.GetName()
      << " is not in dataset and is also not dependent on data set" << endl ;
      return 0 ;
    }

    // Clone derived variable
    tableSet = (RooArgSet*) RooArgSet(cat).snapshot(kTRUE) ;
    if (!tableSet) {
      coutE(Plotting) << "RooTreeData::table(" << GetName() << ") Couldn't deep-clone table category, abort." << endl ;
      return 0 ;
    }
    tableVar = (RooAbsCategory*) tableSet->find(cat.GetName()) ;
    ownPlotVar = kTRUE ;

    //Redirect servers of derived clone to internal ArgSet representing the data in this set
    tableVar->recursiveRedirectServers(_vars) ;
  }

  TString tableName(GetName()) ;
  if (cuts && strlen(cuts)) {
    tableName.Append("(") ;
    tableName.Append(cuts) ;
    tableName.Append(")") ;
  }
  Roo1DTable* table2 = tableVar->createTable(tableName) ;

  // Make cut selector if cut is specified
  RooFormulaVar* cutVar = 0;
  if (cuts && strlen(cuts)) {
    cutVar = new RooFormulaVar("cutVar",cuts,_vars) ;
  }

  // Dump contents
  Int_t nevent= numEntries() ;
  for(Int_t i=0; i < nevent; ++i) {
    get(i);

    if (cutVar && cutVar->getVal()==0) continue ;

    table2->fill(*tableVar,weight()) ;
  }

  if (ownPlotVar) delete tableSet ;
  if (cutVar) delete cutVar ;

  return table2 ;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill Doubles 'lowest' and 'highest' with the lowest and highest value of
/// observable 'var' in this dataset. If the return value is kTRUE and error
/// occurred

Bool_t RooAbsData::getRange(const RooAbsRealLValue& var, Double_t& lowest, Double_t& highest, Double_t marginFrac, Bool_t symMode) const
{
  // Lookup variable in dataset
  const auto arg = _vars.find(var.GetName());
  if (!arg) {
    coutE(InputArguments) << "RooDataSet::getRange(" << GetName() << ") ERROR: unknown variable: " << var.GetName() << endl ;
    return kTRUE;
  }

  auto varPtr = dynamic_cast<const RooRealVar*>(arg);
  // Check if found variable is of type RooRealVar
  if (!varPtr) {
    coutE(InputArguments) << "RooDataSet::getRange(" << GetName() << ") ERROR: variable " << var.GetName() << " is not of type RooRealVar" << endl ;
    return kTRUE;
  }

  // Check if dataset is not empty
  if(sumEntries() == 0.) {
    coutE(InputArguments) << "RooDataSet::getRange(" << GetName() << ") WARNING: empty dataset" << endl ;
    return kTRUE;
  }

  // Look for highest and lowest value
  lowest = RooNumber::infinity() ;
  highest = -RooNumber::infinity() ;
  for (Int_t i=0 ; i<numEntries() ; i++) {
    get(i) ;
    if (varPtr->getVal()<lowest) {
      lowest = varPtr->getVal() ;
    }
    if (varPtr->getVal()>highest) {
      highest = varPtr->getVal() ;
    }
  }

  if (marginFrac>0) {
    if (symMode==kFALSE) {

      Double_t margin = marginFrac*(highest-lowest) ;
      lowest -= margin ;
      highest += margin ;
      if (lowest<var.getMin()) lowest = var.getMin() ;
      if (highest>var.getMax()) highest = var.getMax() ;

    } else {

      Double_t mom1 = moment(*varPtr,1) ;
      Double_t delta = ((highest-mom1)>(mom1-lowest)?(highest-mom1):(mom1-lowest))*(1+marginFrac) ;
      lowest = mom1-delta ;
      highest = mom1+delta ;
      if (lowest<var.getMin()) lowest = var.getMin() ;
      if (highest>var.getMax()) highest = var.getMax() ;

    }
  }

  return kFALSE ;
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare dataset for use with cached constant terms listed in
/// 'cacheList' of expression 'arg'. Deactivate tree branches
/// for any dataset observable that is either not used at all,
/// or is used exclusively by cached branch nodes.

void RooAbsData::optimizeReadingWithCaching(RooAbsArg& arg, const RooArgSet& cacheList, const RooArgSet& keepObsList)
{
  RooArgSet pruneSet ;

  // Add unused observables in this dataset to pruneSet
  pruneSet.add(*get()) ;
  RooArgSet* usedObs = arg.getObservables(*this) ;
  pruneSet.remove(*usedObs,kTRUE,kTRUE) ;

  // Add observables exclusively used to calculate cached observables to pruneSet
  TIterator* vIter = get()->createIterator() ;
  RooAbsArg *var ;
  while ((var=(RooAbsArg*) vIter->Next())) {
    if (allClientsCached(var,cacheList)) {
      pruneSet.add(*var) ;
    }
  }
  delete vIter ;


  if (pruneSet.getSize()!=0) {

    // Go over all used observables and check if any of them have parameterized
    // ranges in terms of pruned observables. If so, remove those observable
    // from the pruning list
    TIterator* uIter = usedObs->createIterator() ;
    RooAbsArg* obs ;
    while((obs=(RooAbsArg*)uIter->Next())) {
      RooRealVar* rrv = dynamic_cast<RooRealVar*>(obs) ;
      if (rrv && !rrv->getBinning().isShareable()) {
   RooArgSet depObs ;
   RooAbsReal* loFunc = rrv->getBinning().lowBoundFunc() ;
   RooAbsReal* hiFunc = rrv->getBinning().highBoundFunc() ;
   if (loFunc) {
     loFunc->leafNodeServerList(&depObs,0,kTRUE) ;
   }
   if (hiFunc) {
     hiFunc->leafNodeServerList(&depObs,0,kTRUE) ;
   }
   if (depObs.getSize()>0) {
     pruneSet.remove(depObs,kTRUE,kTRUE) ;
   }
      }
    }
    delete uIter ;
  }


  // Remove all observables in keep list from prune list
  pruneSet.remove(keepObsList,kTRUE,kTRUE) ;

  if (pruneSet.getSize()!=0) {

    // Deactivate tree branches here
    cxcoutI(Optimization) << "RooTreeData::optimizeReadingForTestStatistic(" << GetName() << "): Observables " << pruneSet
             << " in dataset are either not used at all, orserving exclusively p.d.f nodes that are now cached, disabling reading of these observables for TTree" << endl ;
    setArgStatus(pruneSet,kFALSE) ;
  }

  delete usedObs ;

}

////////////////////////////////////////////////////////////////////////////////
/// Utility function that determines if all clients of object 'var'
/// appear in given list of cached nodes.

Bool_t RooAbsData::allClientsCached(RooAbsArg* var, const RooArgSet& cacheList)
{
  Bool_t ret(kTRUE), anyClient(kFALSE) ;

  for (const auto client : var->valueClients()) {
    anyClient = kTRUE ;
    if (!cacheList.find(client->GetName())) {
      // If client is not cached recurse
      ret &= allClientsCached(client,cacheList) ;
    }
  }

  return anyClient?ret:kFALSE ;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::attachBuffers(const RooArgSet& extObs)
{
  _dstore->attachBuffers(extObs) ;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::resetBuffers()
{
  _dstore->resetBuffers() ;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t RooAbsData::canSplitFast() const
{
  if (_ownedComponents.size()>0) {
    return kTRUE ;
  }
  return kFALSE ;
}

////////////////////////////////////////////////////////////////////////////////

RooAbsData* RooAbsData::getSimData(const char* name)
{
  map<string,RooAbsData*>::iterator i = _ownedComponents.find(name) ;
  if (i==_ownedComponents.end()) return 0 ;
  return i->second ;
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::addOwnedComponent(const char* idxlabel, RooAbsData& data)
{
  _ownedComponents[idxlabel]= &data ;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooAbsData.

void RooAbsData::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RooAbsData::Class(),this);

      // Convert on the fly to vector storage if that the current working default
      if (defaultStorageType==RooAbsData::Vector) {
   convertToVectorStore() ;
      }

   } else {
      R__b.WriteClassBuffer(RooAbsData::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////

void RooAbsData::checkInit() const
{
  _dstore->checkInit() ;
}

////////////////////////////////////////////////////////////////////////////////
/// Forward draw command to data store

void RooAbsData::Draw(Option_t* option)
{
  if (_dstore) _dstore->Draw(option) ;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t RooAbsData::hasFilledCache() const
{
  return _dstore->hasFilledCache() ;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TTree which stores the data. Returns a nullpointer
/// if vector-based storage is used. The RooAbsData remains owner of the tree.
/// GetClonedTree() can be used to get a tree even if the internal storage does not use one.

const TTree *RooAbsData::tree() const
{
   if (storageType == RooAbsData::Tree) {
      return _dstore->tree();
   } else {
      coutW(InputArguments) << "RooAbsData::tree(" << GetName() << ") WARNING: is not of StorageType::Tree. "
                            << "Use GetClonedTree() instead or convert to tree storage." << endl;
      return (TTree *)nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a clone of the TTree which stores the data or create such a tree
/// if vector storage is used. The user is responsible for deleting the tree

TTree *RooAbsData::GetClonedTree() const
{
   if (storageType == RooAbsData::Tree) {
      auto tmp = const_cast<TTree *>(_dstore->tree());
      return tmp->CloneTree();
   } else {
      RooTreeDataStore buffer(GetName(), GetTitle(), *get(), *_dstore);
      return buffer.tree().CloneTree();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert vector-based storage to tree-based storage

void RooAbsData::convertToTreeStore()
{
   if (storageType != RooAbsData::Tree) {
      RooTreeDataStore *newStore = new RooTreeDataStore(GetName(), GetTitle(), _vars, *_dstore);
      delete _dstore;
      _dstore = newStore;
      storageType = RooAbsData::Tree;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If one of the TObject we have a referenced to is deleted, remove the
/// reference.

void RooAbsData::RecursiveRemove(TObject *obj)
{
  for(auto &iter : _ownedComponents) {
    if (iter.second == obj) {
      iter.second = nullptr;
    }
  }
}
