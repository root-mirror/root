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

//////////////////////////////////////////////////////////////////////////////
/// \class RooFitResult
/// RooFitResult is a container class to hold the input and output
/// of a PDF fit to a dataset. It contains:
///
///   * Values of all constant parameters
///   * Initial and final values of floating parameters with error
///   * Correlation matrix and global correlation coefficients
///   * NLL and EDM at mininum
///
/// No references to the fitted PDF and dataset are stored
///

#include "RooFit.h"

#include <iostream>
#include <iomanip>

#include "TBuffer.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TMarker.h"
#include "TLine.h"
#include "TBox.h"
#include "TGaxis.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TClass.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooEllipse.h"
#include "RooRandom.h"
#include "RooMsgService.h"
#include "TH2D.h"
#include "TText.h"
#include "TMatrixDSym.h"
#include "RooMultiVarGaussian.h"


using namespace std;

ClassImp(RooFitResult);



////////////////////////////////////////////////////////////////////////////////
/// Constructor with name and title

RooFitResult::RooFitResult(const char* name, const char* title) : 
  TNamed(name,title), _constPars(0), _initPars(0), _finalPars(0), _globalCorr(0), _randomPars(0), _Lt(0),
  _CM(0), _VM(0), _GC(0)
{  
  if (name) appendToDir(this,kTRUE) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooFitResult::RooFitResult(const RooFitResult& other) : 
  TNamed(other),
  RooPrintable(other),
  RooDirItem(other),
  _status(other._status),
  _covQual(other._covQual),
  _numBadNLL(other._numBadNLL),
  _minNLL(other._minNLL),
  _edm(other._edm),
  _globalCorr(0),
  _randomPars(0),
  _Lt(0),
  _CM(0),
  _VM(0),
  _GC(0),
  _statusHistory(other._statusHistory)
{
  _constPars = (RooArgList*) other._constPars->snapshot() ;
  _initPars = (RooArgList*) other._initPars->snapshot() ;
  _finalPars = (RooArgList*) other._finalPars->snapshot() ;
  if (other._randomPars) _randomPars = (RooArgList*) other._randomPars->snapshot() ;
  if (other._Lt) _Lt = new TMatrix(*other._Lt);
  if (other._VM) _VM = new TMatrixDSym(*other._VM) ;
  if (other._CM) _CM = new TMatrixDSym(*other._CM) ;
  if (other._GC) _GC = new TVectorD(*other._GC) ;

  if (GetName())
    appendToDir(this, kTRUE);
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooFitResult::~RooFitResult() 
{
  if (_constPars) delete _constPars ;
  if (_initPars)  delete _initPars ;
  if (_finalPars) delete _finalPars ;
  if (_globalCorr) delete _globalCorr;
  if (_randomPars) delete _randomPars;
  if (_Lt) delete _Lt;
  if (_CM) delete _CM ;
  if (_VM) delete _VM ;
  if (_GC) delete _GC ;
  
  _corrMatrix.RemoveAll();
  _corrMatrix.Delete();

  removeFromDir(this) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the list of constant parameters

void RooFitResult::setConstParList(const RooArgList& list) 
{
  if (_constPars) delete _constPars ;
  _constPars = (RooArgList*) list.snapshot() ;
  TIterator* iter = _constPars->createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    RooRealVar* rrv = dynamic_cast<RooRealVar*>(arg) ;
    if (rrv) {
      rrv->deleteSharedProperties() ;
    }
  }
  delete iter ;
}



////////////////////////////////////////////////////////////////////////////////
/// Fill the list of initial values of the floating parameters 

void RooFitResult::setInitParList(const RooArgList& list)
{
  if (_initPars) delete _initPars ;
  _initPars = (RooArgList*) list.snapshot() ;
  TIterator* iter = _initPars->createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    RooRealVar* rrv = dynamic_cast<RooRealVar*>(arg) ;
    if (rrv) {
      rrv->deleteSharedProperties() ;
    }
  }
  delete iter ;
}



////////////////////////////////////////////////////////////////////////////////
/// Fill the list of final values of the floating parameters 

void RooFitResult::setFinalParList(const RooArgList& list)
{
  if (_finalPars) delete _finalPars ;
  _finalPars = (RooArgList*) list.snapshot() ;

  TIterator* iter = _finalPars->createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    RooRealVar* rrv = dynamic_cast<RooRealVar*>(arg) ;
    if (rrv) {
      rrv->deleteSharedProperties() ;
    }
  }
  delete iter ;
}



////////////////////////////////////////////////////////////////////////////////

Int_t RooFitResult::statusCodeHistory(UInt_t icycle) const
{ 
  if (icycle>=_statusHistory.size()) {
    coutE(InputArguments) << "RooFitResult::statusCodeHistory(" << GetName() 
			  << " ERROR request for status history slot " 
			  << icycle << " exceeds history count of " << _statusHistory.size() << endl ;
  }
  return _statusHistory[icycle].second ; 
}



////////////////////////////////////////////////////////////////////////////////

const char* RooFitResult::statusLabelHistory(UInt_t icycle) const 
{ 
  if (icycle>=_statusHistory.size()) {
    coutE(InputArguments) << "RooFitResult::statusLabelHistory(" << GetName() 
			  << " ERROR request for status history slot " 
			  << icycle << " exceeds history count of " << _statusHistory.size() << endl ;
  }
  return _statusHistory[icycle].first.c_str() ; 
}



////////////////////////////////////////////////////////////////////////////////
/// Add objects to a 2D plot that represent the fit results for the
/// two named parameters.  The input frame with the objects added is
/// returned, or zero in case of an error.  Which objects are added
/// are determined by the options string which should be a concatenation
/// of the following (not case sensitive):
///
/// *  M - a marker at the best fit result
/// *  E - an error ellipse calculated at 1-sigma using the error matrix at the minimum
/// *  1 - the 1-sigma error bar for parameter 1
/// *  2 - the 1-sigma error bar for parameter 2
/// *  B - the bounding box for the error ellipse
/// *  H - a line and horizontal axis for reading off the correlation coefficient
/// *  V - a line and vertical axis for reading off the correlation coefficient
/// *  A - draw axes for reading off the correlation coefficients with the H or V options
///
/// You can change the attributes of objects in the returned RooPlot using the
/// various `RooPlot::getAttXxx(name)` member functions, e.g.
/// ```
///   plot->getAttLine("contour")->SetLineStyle(kDashed);
/// ```
/// Use `plot->Print()` for a list of all objects and their names (unfortunately most
/// of the ROOT builtin graphics objects like TLine are unnamed). Drag the left mouse
/// button along the labels of either axis button to interactively zoom in a plot.

RooPlot *RooFitResult::plotOn(RooPlot *frame, const char *parName1, const char *parName2,
			      const char *options) const 
{
  // lookup the input parameters by name: we require that they were floated in our fit
  const RooRealVar *par1= dynamic_cast<const RooRealVar*>(floatParsFinal().find(parName1));
  if(0 == par1) {
    coutE(InputArguments) << "RooFitResult::correlationPlot: parameter not floated in fit: " << parName1 << endl;
    return 0;
  }
  const RooRealVar *par2= dynamic_cast<const RooRealVar*>(floatParsFinal().find(parName2));
  if(0 == par2) {
    coutE(InputArguments) << "RooFitResult::correlationPlot: parameter not floated in fit: " << parName2 << endl;
    return 0;
  }

  // options are not case sensitive
  TString opt(options);
  opt.ToUpper();

  // lookup the 2x2 covariance matrix elements for these variables
  Double_t x1= par1->getVal();
  Double_t x2= par2->getVal();
  Double_t s1= par1->getError();
  Double_t s2= par2->getError();
  Double_t rho= correlation(parName1, parName2);

  // add a 1-sigma error ellipse, if requested
  if(opt.Contains("E")) {
    RooEllipse *contour= new RooEllipse("contour",x1,x2,s1,s2,rho);
    contour->SetLineWidth(2) ;
    frame->addPlotable(contour);
  }

  // add the error bar for parameter 1, if requested
  if(opt.Contains("1")) {
    TLine *hline= new TLine(x1-s1,x2,x1+s1,x2);
    hline->SetLineColor(kRed);
    frame->addObject(hline);
  }

  if(opt.Contains("2")) {
    TLine *vline= new TLine(x1,x2-s2,x1,x2+s2);
    vline->SetLineColor(kRed);
    frame->addObject(vline);
  }

  if(opt.Contains("B")) {
    TBox *box= new TBox(x1-s1,x2-s2,x1+s1,x2+s2);
    box->SetLineStyle(kDashed);
    box->SetLineColor(kRed);
    box->SetFillStyle(0);
    frame->addObject(box);
  }

  if(opt.Contains("H")) {
    TLine *line= new TLine(x1-rho*s1,x2-s2,x1+rho*s1,x2+s2);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2) ;
    frame->addObject(line);
    if(opt.Contains("A")) {
      TGaxis *axis= new TGaxis(x1-s1,x2-s2,x1+s1,x2-s2,-1.,+1.,502,"-=");
      axis->SetLineColor(kBlue);
      frame->addObject(axis);
    }
  }

  if(opt.Contains("V")) {
    TLine *line= new TLine(x1-s1,x2-rho*s2,x1+s1,x2+rho*s2);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2) ;
    frame->addObject(line);
    if(opt.Contains("A")) {
      TGaxis *axis= new TGaxis(x1-s1,x2-s2,x1-s1,x2+s2,-1.,+1.,502,"-=");
      axis->SetLineColor(kBlue);
      frame->addObject(axis);
    }
  }

  // add a marker at the fitted value, if requested
  if(opt.Contains("M")) {
    TMarker *marker= new TMarker(x1,x2,20);
    marker->SetMarkerColor(kBlack);
    frame->addObject(marker);
  }

  return frame;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a list of floating parameter values that are perturbed from the final
/// fit values by random amounts sampled from the covariance matrix. The returned
/// object is overwritten with each call and belongs to the RooFitResult. Uses
/// the "square root method" to decompose the covariance matrix, which makes inverting
/// it unnecessary.

const RooArgList& RooFitResult::randomizePars() const 
{
  Int_t nPar= _finalPars->getSize();
  if(0 == _randomPars) { // first-time initialization
    assert(0 != _finalPars);
    // create the list of random values to fill
    _randomPars= (RooArgList*)_finalPars->snapshot();
    // calculate the elements of the upper-triangular matrix L that gives Lt*L = C
    // where Lt is the transpose of L (the "square-root method")
    TMatrix L(nPar,nPar);
    for(Int_t iPar= 0; iPar < nPar; iPar++) {
      // calculate the diagonal term first
      L(iPar,iPar)= covariance(iPar,iPar);
      for(Int_t k= 0; k < iPar; k++) {
	Double_t tmp= L(k,iPar);
	L(iPar,iPar)-= tmp*tmp;
      }
      L(iPar,iPar)= sqrt(L(iPar,iPar));
      // then the off-diagonal terms
      for(Int_t jPar= iPar+1; jPar < nPar; jPar++) {
	L(iPar,jPar)= covariance(iPar,jPar);
	for(Int_t k= 0; k < iPar; k++) {
	  L(iPar,jPar)-= L(k,iPar)*L(k,jPar);
	}
	L(iPar,jPar)/= L(iPar,iPar);
      }
    }
    // remember Lt
    _Lt= new TMatrix(TMatrix::kTransposed,L);
  }
  else {
    // reset to the final fit values
    *_randomPars= *_finalPars;
  }

  // create a vector of unit Gaussian variables
  TVector g(nPar);
  for(Int_t k= 0; k < nPar; k++) g(k)= RooRandom::gaussian();
  // multiply this vector by Lt to introduce the appropriate correlations
  g*= (*_Lt);
  // add the mean value offsets and store the results
  TIterator *iter= _randomPars->createIterator();
  RooRealVar *par(0);
  Int_t index(0);
  while((0 != (par= (RooRealVar*)iter->Next()))) {
    par->setVal(par->getVal() + g(index++));
  }
  delete iter;

  return *_randomPars;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the correlation between parameters 'par1' and 'par2'

Double_t RooFitResult::correlation(const char* parname1, const char* parname2) const 
{
  Int_t idx1 = _finalPars->index(parname1) ;
  Int_t idx2 = _finalPars->index(parname2) ;
  if (idx1<0) {
    coutE(InputArguments) << "RooFitResult::correlation(" << GetName() << ") parameter " << parname1 << " is not a floating fit parameter" << endl ;
    return 0 ;
  }
  if (idx2<0) {
    coutE(InputArguments) << "RooFitResult::correlation(" << GetName() << ") parameter " << parname2 << " is not a floating fit parameter" << endl ;
    return 0 ;
  }
  return correlation(idx1,idx2) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the set of correlation coefficients of parameter 'par' with
/// all other floating parameters

const RooArgList* RooFitResult::correlation(const char* parname) const 
{
  if (_globalCorr==0) {
    fillLegacyCorrMatrix() ;
  }

  RooAbsArg* arg = _initPars->find(parname) ;
  if (!arg) {
    coutE(InputArguments) << "RooFitResult::correlation: variable " << parname << " not a floating parameter in fit" << endl ;
    return 0 ;
  }    
  return (RooArgList*)_corrMatrix.At(_initPars->index(arg)) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the global correlation of the named parameter

Double_t RooFitResult::globalCorr(const char* parname) 
{
  if (_globalCorr==0) {
    fillLegacyCorrMatrix() ;
  }

  RooAbsArg* arg = _initPars->find(parname) ;
  if (!arg) {
    coutE(InputArguments) << "RooFitResult::globalCorr: variable " << parname << " not a floating parameter in fit" << endl ;
    return 0 ;
  }    

  if (_globalCorr) {
    return ((RooAbsReal*)_globalCorr->at(_initPars->index(arg)))->getVal() ;
  } else {
    return 1.0 ; 
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Return the list of all global correlations

const RooArgList* RooFitResult::globalCorr() 
{
  if (_globalCorr==0) {
    fillLegacyCorrMatrix() ;
  }

  return _globalCorr ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return a correlation matrix element addressed with numeric indices.

Double_t RooFitResult::correlation(Int_t row, Int_t col) const 
{
  return (*_CM)(row,col) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the covariance matrix element addressed with numeric indices.

Double_t RooFitResult::covariance(Int_t row, Int_t col) const 
{
  return (*_VM)(row,col) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print fit result to stream 'os'. In Verbose mode, the constant parameters and
/// the initial and final values of the floating parameters are printed. 
/// Standard mode only the final values of the floating parameters are printed

void RooFitResult::printMultiline(ostream& os, Int_t /*contents*/, Bool_t verbose, TString indent) const
{

  os << endl 
     << indent << "  RooFitResult: minimized FCN value: " << _minNLL << ", estimated distance to minimum: " << _edm << endl
     << indent << "                covariance matrix quality: " ;
  switch(_covQual) {
  case -1 : os << "Unknown, matrix was externally provided" ; break ;
  case 0  : os << "Not calculated at all" ; break ;
  case 1  : os << "Approximation only, not accurate" ; break ;
  case 2  : os << "Full matrix, but forced positive-definite" ; break ;
  case 3  : os << "Full, accurate covariance matrix" ; break ;
  }
  os << endl ; 
  os << indent << "                Status : " ;
  for (vector<pair<string,int> >::const_iterator iter = _statusHistory.begin() ; iter != _statusHistory.end() ; ++iter) {
    os << iter->first << "=" << iter->second << " " ;
  }
  os << endl << endl ;;

  Int_t i ;
  if (verbose) {
    if (_constPars->getSize()>0) {
      os << indent << "    Constant Parameter    Value     " << endl
	 << indent << "  --------------------  ------------" << endl ;

      for (i=0 ; i<_constPars->getSize() ; i++) {
	os << indent << "  " << setw(20) << ((RooAbsArg*)_constPars->at(i))->GetName()
	   << "  " << setw(12) << Form("%12.4e",((RooRealVar*)_constPars->at(i))->getVal())
	   << endl ;
      }

      os << endl ;
    }

    // Has any parameter asymmetric errors?
    Bool_t doAsymErr(kFALSE) ;
    for (i=0 ; i<_finalPars->getSize() ; i++) {
      if (((RooRealVar*)_finalPars->at(i))->hasAsymError()) {
	doAsymErr=kTRUE ;
	break ;
      }
    }

    if (doAsymErr) {
      os << indent << "    Floating Parameter  InitialValue    FinalValue (+HiError,-LoError)    GblCorr." << endl
	 << indent << "  --------------------  ------------  ----------------------------------  --------" << endl ;
    } else {
      os << indent << "    Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr." << endl
	 << indent << "  --------------------  ------------  --------------------------  --------" << endl ;
    }

    for (i=0 ; i<_finalPars->getSize() ; i++) {
      os << indent << "  "    << setw(20) << ((RooAbsArg*)_finalPars->at(i))->GetName() ;
      os << indent << "  "    << setw(12) << Form("%12.4e",((RooRealVar*)_initPars->at(i))->getVal())
	 << indent << "  "    << setw(12) << Form("%12.4e",((RooRealVar*)_finalPars->at(i))->getVal()) ;
      
      if (((RooRealVar*)_finalPars->at(i))->hasAsymError()) {
	os << setw(21) << Form(" (+%8.2e,-%8.2e)",((RooRealVar*)_finalPars->at(i))->getAsymErrorHi(),
	                       -1*((RooRealVar*)_finalPars->at(i))->getAsymErrorLo()) ;
      } else {
	Double_t err = ((RooRealVar*)_finalPars->at(i))->getError() ;
	os << (doAsymErr?"        ":"") << " +/- " << setw(9)  << Form("%9.2e",err) ;
      }

      if (_globalCorr) {
	os << "  "    << setw(8)  << Form("%8.6f" ,((RooRealVar*)_globalCorr->at(i))->getVal()) ;
      } else {
	os << "  <none>" ;
      } 

      os << endl ;
    }

  } else {
    os << indent << "    Floating Parameter    FinalValue +/-  Error   " << endl
       << indent << "  --------------------  --------------------------" << endl ;

    for (i=0 ; i<_finalPars->getSize() ; i++) {
      Double_t err = ((RooRealVar*)_finalPars->at(i))->getError() ;
      os << indent << "  "    << setw(20) << ((RooAbsArg*)_finalPars->at(i))->GetName()
	 << "  "    << setw(12) << Form("%12.4e",((RooRealVar*)_finalPars->at(i))->getVal())
	 << " +/- " << setw(9)  << Form("%9.2e",err)
	 << endl ;
    }
  }
  

  os << endl ;
}


////////////////////////////////////////////////////////////////////////////////
/// Function called by RooMinimizer

void RooFitResult::fillCorrMatrix(const std::vector<double>& globalCC, const TMatrixDSym& corrs, const TMatrixDSym& covs)
{
  // Sanity check
  if (globalCC.empty() || corrs.GetNoElements() < 1 || covs.GetNoElements() < 1) {
    coutI(Minimization) << "RooFitResult::fillCorrMatrix: number of floating parameters is zero, correlation matrix not filled" << endl ;
    return ;
  }

  if (!_initPars) {
    coutE(Minimization) << "RooFitResult::fillCorrMatrix: ERROR: list of initial parameters must be filled first" << endl ;
    return ;
  }

  // Delete eventual prevous correlation data holders
  if (_CM) delete _CM ;
  if (_VM) delete _VM ;
  if (_GC) delete _GC ;

  // Build holding arrays for correlation coefficients
  _CM = new TMatrixDSym(corrs) ;
  _VM = new TMatrixDSym(covs) ;
  _GC = new TVectorD(_CM->GetNcols()) ;
  for(int i=0 ; i<_CM->GetNcols() ; i++) {
    (*_GC)[i] = globalCC[i] ;
  }
  //fillLegacyCorrMatrix() ;
}





////////////////////////////////////////////////////////////////////////////////
/// Sanity check

void RooFitResult::fillLegacyCorrMatrix() const 
{
  if (!_CM) return ;

  // Delete eventual prevous correlation data holders
  if (_globalCorr) delete _globalCorr ;
  _corrMatrix.Delete();

  // Build holding arrays for correlation coefficients
  _globalCorr = new RooArgList("globalCorrelations") ;

  TIterator* vIter = _initPars->createIterator() ;
  RooAbsArg* arg ;
  Int_t idx(0) ;
  while((arg=(RooAbsArg*)vIter->Next())) {
    // Create global correlation value holder
    TString gcName("GC[") ;
    gcName.Append(arg->GetName()) ;
    gcName.Append("]") ;
    TString gcTitle(arg->GetTitle()) ;
    gcTitle.Append(" Global Correlation") ;
    _globalCorr->addOwned(*(new RooRealVar(gcName.Data(),gcTitle.Data(),0.))) ;

    // Create array with correlation holders for this parameter
    TString name("C[") ;
    name.Append(arg->GetName()) ;
    name.Append(",*]") ;
    RooArgList* corrMatrixRow = new RooArgList(name.Data()) ;
    _corrMatrix.Add(corrMatrixRow) ;
    TIterator* vIter2 = _initPars->createIterator() ;
    RooAbsArg* arg2 ;
    while((arg2=(RooAbsArg*)vIter2->Next())) {

      TString cName("C[") ;
      cName.Append(arg->GetName()) ;
      cName.Append(",") ;
      cName.Append(arg2->GetName()) ;
      cName.Append("]") ;
      TString cTitle("Correlation between ") ;
      cTitle.Append(arg->GetName()) ;
      cTitle.Append(" and ") ;
      cTitle.Append(arg2->GetName()) ;
      corrMatrixRow->addOwned(*(new RooRealVar(cName.Data(),cTitle.Data(),0.))) ;      
    }
    delete vIter2 ;
    idx++ ;
  }
  delete vIter ;

  TIterator *gcIter = _globalCorr->createIterator() ;
  TIterator *parIter = _finalPars->createIterator() ;
  RooRealVar* gcVal = 0;
  for (unsigned int i = 0; i < (unsigned int)_CM->GetNcols() ; ++i) {

    // Find the next global correlation slot to fill, skipping fixed parameters
    gcVal = (RooRealVar*) gcIter->Next() ;
    gcVal->setVal((*_GC)(i)) ; // WVE FIX THIS 

    // Fill a row of the correlation matrix
    TIterator* cIter = ((RooArgList*)_corrMatrix.At(i))->createIterator() ;
    for (unsigned int it = 0; it < (unsigned int)_CM->GetNcols() ; ++it) {
      RooRealVar* cVal = (RooRealVar*) cIter->Next() ;
      double value = (*_CM)(i,it) ;
      cVal->setVal(value);      
      (*_CM)(i,it) = value;
    }
    delete cIter ;
  }

  delete gcIter ;
  delete parIter ;

}





////////////////////////////////////////////////////////////////////////////////
/// Internal utility method to extract the correlation matrix and the
/// global correlation coefficients from the MINUIT memory buffer and
/// fill the internal arrays.

void RooFitResult::fillCorrMatrix()
{
  // Sanity check
  if (gMinuit->fNpar < 1) {
    coutI(Minimization) << "RooFitResult::fillCorrMatrix: number of floating parameters is zero, correlation matrix not filled" << endl ;
    return ;
  }

  if (!_initPars) {
    coutE(Minimization) << "RooFitResult::fillCorrMatrix: ERROR: list of initial parameters must be filled first" << endl ;
    return ;
  }

  // Delete eventual prevous correlation data holders
  if (_CM) delete _CM ;
  if (_VM) delete _VM ;
  if (_GC) delete _GC ;

  // Build holding arrays for correlation coefficients
  _CM = new TMatrixDSym(_initPars->getSize()) ;
  _VM = new TMatrixDSym(_initPars->getSize()) ;
  _GC = new TVectorD(_initPars->getSize()) ;

  // Extract correlation information for MINUIT (code taken from TMinuit::mnmatu() )

  // WVE: This code directly manipulates minuit internal workspace, 
  //      if TMinuit code changes this may need updating
  Int_t ndex, i, j, m, n, it /* nparm,id,ix */ ;
  Int_t ndi, ndj /*, iso, isw2, isw5*/;
  for (i = 1; i <= gMinuit->fNpar; ++i) {
    ndi = i*(i + 1) / 2;
    for (j = 1; j <= gMinuit->fNpar; ++j) {
      m    = TMath::Max(i,j);
      n    = TMath::Min(i,j);
      ndex = m*(m-1) / 2 + n;
      ndj  = j*(j + 1) / 2;
      gMinuit->fMATUvline[j-1] = gMinuit->fVhmat[ndex-1] / TMath::Sqrt(TMath::Abs(gMinuit->fVhmat[ndi-1]*gMinuit->fVhmat[ndj-1]));
    }

    (*_GC)(i-1) = gMinuit->fGlobcc[i-1] ;

    // Fill a row of the correlation matrix
    for (it = 1; it <= gMinuit->fNpar ; ++it) {
      (*_CM)(i-1,it-1) = gMinuit->fMATUvline[it-1] ;
    }
  }

  for (int ii=0 ; ii<_finalPars->getSize() ; ii++) {
    for (int jj=0 ; jj<_finalPars->getSize() ; jj++) {
      (*_VM)(ii,jj) = (*_CM)(ii,jj) * ((RooRealVar*)_finalPars->at(ii))->getError() * ((RooRealVar*)_finalPars->at(jj))->getError() ;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooFitResult::fillPrefitCorrMatrix()
{

   // Delete eventual prevous correlation data holders
   if (_CM)
      delete _CM;
   if (_VM)
      delete _VM;
   if (_GC)
      delete _GC;

   // Build holding arrays for correlation coefficients
   _CM = new TMatrixDSym(_initPars->getSize());
   _VM = new TMatrixDSym(_initPars->getSize());
   _GC = new TVectorD(_initPars->getSize());

   for (int ii = 0; ii < _finalPars->getSize(); ii++) {
      (*_CM)(ii, ii) = 1;
      (*_VM)(ii, ii) = ((RooRealVar *)_finalPars->at(ii))->getError() * ((RooRealVar *)_finalPars->at(ii))->getError();
      (*_GC)(ii) = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if this fit result is identical to other within tolerances.
/// \param[in] other Fit result to test against.
/// \param[in] tol **Relative** tolerance for parameters and NLL.
/// \param[in] tolCorr **absolute** tolerance for correlation coefficients.
/// \param[in] verbose Ignored.

Bool_t RooFitResult::isIdentical(const RooFitResult& other, Double_t tol, Double_t tolCorr, Bool_t /*verbose*/) const 
{
  Bool_t ret = kTRUE ;
  auto deviation = [tol](const double left, const double right){
    if (right != 0.)
      return fabs((left - right)/right) >= tol;
    else
      return fabs(left) >= tol;
  };

  auto errMsg = [](std::string msgHead, const RooAbsReal* tv, const RooAbsReal* ov) {
    cout << "RooFitResult::isIdentical: " << msgHead << " " << tv->GetName() << " differs in value:\t"
        << tv->getVal() << " vs.\t" << ov->getVal()
        << "\t(" << (tv->getVal()-ov->getVal())/ov->getVal() << ")" << endl;
  };

  if (deviation(_minNLL, other._minNLL)) {
    cout << "RooFitResult::isIdentical: minimized value of -log(L) is different " << _minNLL << " vs. " << other._minNLL << endl ;
    ret = kFALSE ;
  }

  for (Int_t i=0 ; i<_constPars->getSize() ; i++) {
    auto tv = static_cast<const RooAbsReal*>(_constPars->at(i));
    auto ov = static_cast<const RooAbsReal*>(other._constPars->find(tv->GetName())) ;
    if (!ov) {
      cout << "RooFitResult::isIdentical: cannot find constant parameter " << _constPars->at(i)->GetName() << " in reference" << endl ;
      ret = kFALSE ;
    }
    if (ov && deviation(tv->getVal(), ov->getVal())) {
      errMsg("constant parameter", tv, ov);
      ret = kFALSE ;
    }
  }

  for (Int_t i=0 ; i<_initPars->getSize() ; i++) {
    auto ov = static_cast<const RooAbsReal*>(other._initPars->find(_initPars->at(i)->GetName())) ;
    auto tv = static_cast<const RooAbsReal*>(_initPars->at(i));
    if (!ov) {
      cout << "RooFitResult::isIdentical: cannot find initial parameter " << _initPars->at(i)->GetName() << " in reference" << endl ;
      ret = kFALSE ;
    }
    if (ov && deviation(tv->getVal(), ov->getVal())) {
      errMsg("initial parameter", tv, ov);
      ret = kFALSE ;
    }
  }

  for (Int_t i=0 ; i<_finalPars->getSize() ; i++) {
    auto tv = static_cast<const RooAbsReal*>(_finalPars->at(i));
    auto ov = static_cast<const RooAbsReal*>(other._finalPars->find(tv->GetName())) ;
    if (!ov) {
      cout << "RooFitResult::isIdentical: cannot find final parameter " << tv->GetName() << " in reference" << endl ;
      ret = kFALSE ;
    }
    if (ov && deviation(tv->getVal(), ov->getVal())) {
      errMsg("final parameter", tv, ov);
      ret = kFALSE ;
    }
  }

  auto deviationCorr = [tolCorr](const double left, const double right){
    return fabs(left - right) >= tolCorr;
  };

  // Only examine correlations for cases with >1 floating parameter
  if (_finalPars->getSize()>1) {

    fillLegacyCorrMatrix() ;
    other.fillLegacyCorrMatrix() ;

    for (Int_t i=0 ; i<_globalCorr->getSize() ; i++) {
      auto tv = static_cast<const RooAbsReal*>(_globalCorr->at(i));
      auto ov = static_cast<const RooAbsReal*>(other._globalCorr->find(_globalCorr->at(i)->GetName())) ;
      if (!ov) {
        cout << "RooFitResult::isIdentical: cannot find global correlation coefficient " << tv->GetName() << " in reference" << endl ;
        ret = kFALSE ;
      }
      if (ov && deviationCorr(tv->getVal(), ov->getVal())) {
        errMsg("global correlation coefficient", tv, ov);
        ret = kFALSE ;
      }
    }

    for (Int_t j=0 ; j<_corrMatrix.GetSize() ; j++) {
      RooArgList* row = (RooArgList*) _corrMatrix.At(j) ;
      RooArgList* orow = (RooArgList*) other._corrMatrix.At(j) ;
      for (Int_t i=0 ; i<row->getSize() ; i++) {
        auto tv = static_cast<const RooAbsReal*>(row->at(i));
        auto ov = static_cast<const RooAbsReal*>(orow->find(tv->GetName())) ;
        if (!ov) {
          cout << "RooFitResult::isIdentical: cannot find correlation coefficient " << tv->GetName() << " in reference" << endl ;
          ret = kFALSE ;
        }
        if (ov && deviationCorr(tv->getVal(), ov->getVal())) {
          errMsg("correlation coefficient", tv, ov);
          ret = kFALSE ;
        }
      }
    }
  }    

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Import the results of the last fit performed by gMinuit, interpreting
/// the fit parameters as the given varList of parameters.

RooFitResult* RooFitResult::lastMinuitFit(const RooArgList& varList) 
{
  // Verify length of supplied varList
  if (varList.getSize()>0 && varList.getSize()!=gMinuit->fNu) {
    oocoutE((TObject*)0,InputArguments) << "RooFitResult::lastMinuitFit: ERROR: supplied variable list must be either empty " << endl 
					<< "                             or match the number of variables of the last fit (" << gMinuit->fNu << ")" << endl ;
    return 0 ;
  }

  // Verify that all members of varList are of type RooRealVar
  TIterator* iter = varList.createIterator() ;
  RooAbsArg* arg  ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!dynamic_cast<RooRealVar*>(arg)) {
      oocoutE((TObject*)0,InputArguments) << "RooFitResult::lastMinuitFit: ERROR: variable '" << arg->GetName() << "' is not of type RooRealVar" << endl ;
      return 0 ;
    }
  }
  delete iter ;

  RooFitResult* r = new RooFitResult("lastMinuitFit","Last MINUIT fit") ;

  // Extract names of fit parameters from MINUIT 
  // and construct corresponding RooRealVars
  RooArgList constPars("constPars") ;
  RooArgList floatPars("floatPars") ;

  Int_t i ;
  for (i = 1; i <= gMinuit->fNu; ++i) {
    if (gMinuit->fNvarl[i-1] < 0) continue;
    Int_t l = gMinuit->fNiofex[i-1];
    TString varName(gMinuit->fCpnam[i-1]) ;
    Bool_t isConst(l==0) ;
    
    Double_t xlo = gMinuit->fAlim[i-1];
    Double_t xhi = gMinuit->fBlim[i-1];
    Double_t xerr = gMinuit->fWerr[l-1];
    Double_t xval = gMinuit->fU[i-1] ;

    RooRealVar* var ;
    if (varList.getSize()==0) {

      if ((xlo<xhi) && !isConst) {
	var = new RooRealVar(varName,varName,xval,xlo,xhi) ;
      } else {
	var = new RooRealVar(varName,varName,xval) ;
      }
      var->setConstant(isConst) ;
    } else {

      var = (RooRealVar*) varList.at(i-1)->Clone() ;
      var->setConstant(isConst) ;
      var->setVal(xval) ;
      if (xlo<xhi) {
	var->setRange(xlo,xhi) ;
      }
      if (varName.CompareTo(var->GetName())) {
	oocoutI((TObject*)0,Eval) << "RooFitResult::lastMinuitFit: fit parameter '" << varName 
				  << "' stored in variable '" << var->GetName() << "'" << endl ;
      }

    }

    if (isConst) {
      constPars.addOwned(*var) ;
    } else {
      var->setError(xerr) ;
      floatPars.addOwned(*var) ;
    }
  }

  Int_t icode,npari,nparx ;
  Double_t fmin,edm,errdef ;
  gMinuit->mnstat(fmin,edm,errdef,npari,nparx,icode) ;
  
  r->setConstParList(constPars) ;
  r->setInitParList(floatPars) ;
  r->setFinalParList(floatPars) ;
  r->setMinNLL(fmin) ;
  r->setEDM(edm) ; 
  r->setCovQual(icode) ;
  r->setStatus(gMinuit->fStatus) ;
  r->fillCorrMatrix() ;

  return r ;
}



////////////////////////////////////////////////////////////////////////////////
/// Import the results of the last fit performed by gMinuit, interpreting
/// the fit parameters as the given varList of parameters.

RooFitResult *RooFitResult::prefitResult(const RooArgList &paramList)
{
   // Verify that all members of varList are of type RooRealVar
   TIterator *iter = paramList.createIterator();
   RooAbsArg *arg;
   while ((arg = (RooAbsArg *)iter->Next())) {
      if (!dynamic_cast<RooRealVar *>(arg)) {
         oocoutE((TObject *)0, InputArguments) << "RooFitResult::lastMinuitFit: ERROR: variable '" << arg->GetName()
                                               << "' is not of type RooRealVar" << endl;
         return 0;
      }
   }

   RooFitResult *r = new RooFitResult("lastMinuitFit", "Last MINUIT fit");

   // Extract names of fit parameters from MINUIT
   // and construct corresponding RooRealVars
   RooArgList constPars("constPars");
   RooArgList floatPars("floatPars");

   iter->Reset();
   while ((arg = (RooAbsArg *)iter->Next())) {
      if (arg->isConstant()) {
         constPars.addClone(*arg);
      } else {
         floatPars.addClone(*arg);
      }
   }
   delete iter;

   r->setConstParList(constPars);
   r->setInitParList(floatPars);
   r->setFinalParList(floatPars);
   r->setMinNLL(0);
   r->setEDM(0);
   r->setCovQual(0);
   r->setStatus(0);
   r->fillPrefitCorrMatrix();

   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Store externally provided correlation matrix in this RooFitResult ;

void RooFitResult::setCovarianceMatrix(TMatrixDSym& V) 
{
  // Delete any previous matrices
  if (_VM) {
    delete _VM ;
  }
  if (_CM) {
    delete _CM ;
  }
  
  // Clone input covariance matrix ;
  _VM = (TMatrixDSym*) V.Clone() ;

  // Now construct correlation matrix from it
  _CM = (TMatrixDSym*) _VM->Clone() ;
  for (Int_t i=0 ; i<_CM->GetNrows() ; i++) {
    for (Int_t j=0 ; j<_CM->GetNcols() ; j++) {
      if (i!=j) {
	(*_CM)(i,j) = (*_CM)(i,j) / sqrt((*_CM)(i,i)*(*_CM)(j,j)) ;
      }
    }
  }
  for (Int_t i=0 ; i<_CM->GetNrows() ; i++) {
    (*_CM)(i,i) = 1.0 ;
  }

  _covQual = -1 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return TH2D of correlation matrix 

TH2* RooFitResult::correlationHist(const char* name) const 
{
  Int_t n = _CM->GetNcols() ;

  TH2D* hh = new TH2D(name,name,n,0,n,n,0,n) ;
  
  for (Int_t i = 0 ; i<n ; i++) {
    for (Int_t j = 0 ; j<n; j++) {
      hh->Fill(i+0.5,n-j-0.5,(*_CM)(i,j)) ;
    }
    hh->GetXaxis()->SetBinLabel(i+1,_finalPars->at(i)->GetName()) ;
    hh->GetYaxis()->SetBinLabel(n-i,_finalPars->at(i)->GetName()) ;    
  }
  hh->SetMinimum(-1) ;
  hh->SetMaximum(+1) ;


  return hh ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return covariance matrix 

const TMatrixDSym& RooFitResult::covarianceMatrix() const 
{
  return *_VM ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return a reduced covariance matrix (Note that Vred _is_ a simple sub-matrix of V,
/// row/columns are ordered to matched the convention given in input argument 'params'

TMatrixDSym RooFitResult::reducedCovarianceMatrix(const RooArgList& params) const 
{
  const TMatrixDSym& V = covarianceMatrix() ;


  // Make sure that all given params were floating parameters in the represented fit
  RooArgList params2 ;
  TIterator* iter = params.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (_finalPars->find(arg->GetName())) {
      params2.add(*arg) ;
    } else {
      coutW(InputArguments) << "RooFitResult::reducedCovarianceMatrix(" << GetName() << ") WARNING input variable " 
			    << arg->GetName() << " was not a floating parameters in fit result and is ignored" << endl ;
    }
  }
  delete iter ;
   
   // fix for bug ROOT-8044
   // use same order given bby vector params
   vector<int> indexMap(params2.getSize());
   for (int i=0 ; i<params2.getSize() ; i++) {
      indexMap[i] = _finalPars->index(params2[i].GetName());
      assert(indexMap[i] < V.GetNrows());
   }
   
   TMatrixDSym Vred(indexMap.size());
   for (int i = 0; i < Vred.GetNrows(); ++i) {
      for (int j = 0; j < Vred.GetNcols(); ++j) {
         Vred(i,j) = V( indexMap[i], indexMap[j]);
      }
   }
   return Vred;
}



////////////////////////////////////////////////////////////////////////////////
/// Return a reduced covariance matrix, which is calculated as
/// \f[
///   V_\mathrm{red} = \bar{V_{22}} = V_{11} - V_{12} \cdot V_{22}^{-1} \cdot V_{21},
/// \f]
/// where \f$ V_{11},V_{12},V_{21},V_{22} \f$ represent a block decomposition of the covariance matrix into observables that
/// are propagated (labeled by index '1') and that are not propagated (labeled by index '2'), and \f$ \bar{V_{22}} \f$
/// is the Shur complement of \f$ V_{22} \f$, calculated as shown above.
///
/// (Note that \f$ V_\mathrm{red} \f$ is *not* a simple sub-matrix of \f$ V \f$)

TMatrixDSym RooFitResult::conditionalCovarianceMatrix(const RooArgList& params) const 
{
  const TMatrixDSym& V = covarianceMatrix() ;

  // Handle case where V==Vred here
  if (V.GetNcols()==params.getSize()) {
    return V ;
  }

  Double_t det = V.Determinant() ;

  if (det<=0) {
    coutE(Eval) << "RooFitResult::conditionalCovarianceMatrix(" << GetName() << ") ERROR: covariance matrix is not positive definite (|V|=" 
		<< det << ") cannot reduce it" << endl ;
    throw string("RooFitResult::conditionalCovarianceMatrix() ERROR, input covariance matrix is not positive definite") ;
  }

  // Make sure that all given params were floating parameters in the represented fit
  RooArgList params2 ;
  TIterator* iter = params.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (_finalPars->find(arg->GetName())) {
      params2.add(*arg) ;
    } else {
      coutW(InputArguments) << "RooFitResult::conditionalCovarianceMatrix(" << GetName() << ") WARNING input variable " 
			    << arg->GetName() << " was not a floating parameters in fit result and is ignored" << endl ;
    }
  }
  delete iter ;

  // Need to order params in vector in same order as in covariance matrix
  RooArgList params3 ;
  iter = _finalPars->createIterator() ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (params2.find(arg->GetName())) {
      params3.add(*arg) ;
    }
  }
  delete iter ;

  // Find (subset) of parameters that are stored in the covariance matrix
  vector<int> map1, map2 ;
  for (int i=0 ; i<_finalPars->getSize() ; i++) {
    if (params3.find(_finalPars->at(i)->GetName())) {
      map1.push_back(i) ;
    } else {
      map2.push_back(i) ;
    }
  }

  // Rearrange matrix in block form with 'params' first and 'others' last
  // (preserving relative order) 
  TMatrixDSym S11, S22 ;
  TMatrixD S12, S21 ;
  RooMultiVarGaussian::blockDecompose(V,map1,map2,S11,S12,S21,S22) ;

  // Constructed conditional matrix form         -1
  // F(X1|X2) --> CovI --> S22bar = S11 - S12 S22  S21
  
  // Do eigenvalue decomposition
  TMatrixD S22Inv(TMatrixD::kInverted,S22) ;
  TMatrixD S22bar =  S11 - S12 * (S22Inv * S21) ;

  // Convert explicitly to symmetric form
  TMatrixDSym Vred(S22bar.GetNcols()) ;
  for (int i=0 ; i<Vred.GetNcols() ; i++) {
    for (int j=i ; j<Vred.GetNcols() ; j++) {
      Vred(i,j) = (S22bar(i,j) + S22bar(j,i))/2 ;
      Vred(j,i) = Vred(i,j) ;
    }
  }

  return Vred ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return correlation matrix ;

const TMatrixDSym& RooFitResult::correlationMatrix() const 
{
  return *_CM ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return a p.d.f that represents the fit result as a multi-variate probability densisty
/// function on the floating fit parameters, including correlations

RooAbsPdf* RooFitResult::createHessePdf(const RooArgSet& params) const
{
  const TMatrixDSym& V = covarianceMatrix() ;
  Double_t det = V.Determinant() ;

  if (det<=0) {
    coutE(Eval) << "RooFitResult::createHessePdf(" << GetName() << ") ERROR: covariance matrix is not positive definite (|V|=" 
		<< det << ") cannot construct p.d.f" << endl ;
    return 0 ;
  }

  // Make sure that all given params were floating parameters in the represented fit
  RooArgList params2 ;
  TIterator* iter = params.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (_finalPars->find(arg->GetName())) {
      params2.add(*arg) ;
    } else {
      coutW(InputArguments) << "RooFitResult::createHessePdf(" << GetName() << ") WARNING input variable " 
			    << arg->GetName() << " was not a floating parameters in fit result and is ignored" << endl ;
    }
  }
  delete iter ;

  // Need to order params in vector in same order as in covariance matrix
  RooArgList params3 ;
  iter = _finalPars->createIterator() ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (params2.find(arg->GetName())) {
      params3.add(*arg) ;
    }
  }
  delete iter ;


  // Handle special case of representing full covariance matrix here
  if (params3.getSize()==_finalPars->getSize()) {

    RooArgList mu ;
    for (Int_t i=0 ; i<_finalPars->getSize() ; i++) {
      RooRealVar* parclone = (RooRealVar*) _finalPars->at(i)->Clone(Form("%s_centralvalue",_finalPars->at(i)->GetName())) ;
      parclone->setConstant(kTRUE) ;
      mu.add(*parclone) ;      
    }

    string name  = Form("pdf_%s",GetName()) ;
    string title = Form("P.d.f of %s",GetTitle()) ;
    
    // Create p.d.f.
    RooAbsPdf* mvg = new RooMultiVarGaussian(name.c_str(),title.c_str(),params3,mu,V) ; 
    mvg->addOwnedComponents(mu) ;
    return  mvg ;
  }

  //                                       -> ->
  // Handle case of conditional p.d.f. MVG(p1|p2) here

  // Find (subset) of parameters that are stored in the covariance matrix
  vector<int> map1, map2 ;
  for (int i=0 ; i<_finalPars->getSize() ; i++) {
    if (params3.find(_finalPars->at(i)->GetName())) {
      map1.push_back(i) ;
    } else {
      map2.push_back(i) ;
    }
  }

  // Rearrange matrix in block form with 'params' first and 'others' last
  // (preserving relative order) 
  TMatrixDSym S11, S22 ;
  TMatrixD S12, S21 ;
  RooMultiVarGaussian::blockDecompose(V,map1,map2,S11,S12,S21,S22) ;

  // Calculate offset vectors mu1 and mu2
  RooArgList mu1 ;
  for (UInt_t i=0 ; i<map1.size() ; i++) {
    RooRealVar* parclone = (RooRealVar*) _finalPars->at(map1[i])->Clone(Form("%s_centralvalue",_finalPars->at(map1[i])->GetName())) ;
    parclone->setConstant(kTRUE) ;
    mu1.add(*parclone) ;      
  }

  // Constructed conditional matrix form         -1
  // F(X1|X2) --> CovI --> S22bar = S11 - S12 S22  S21
  
  // Do eigenvalue decomposition
  TMatrixD S22Inv(TMatrixD::kInverted,S22) ;
  TMatrixD S22bar =  S11 - S12 * (S22Inv * S21) ;

  // Convert explicitly to symmetric form
  TMatrixDSym Vred(S22bar.GetNcols()) ;
  for (int i=0 ; i<Vred.GetNcols() ; i++) {
    for (int j=i ; j<Vred.GetNcols() ; j++) {
      Vred(i,j) = (S22bar(i,j) + S22bar(j,i))/2 ;
      Vred(j,i) = Vred(i,j) ;
    }
  }
  string name  = Form("pdf_%s",GetName()) ;
  string title = Form("P.d.f of %s",GetTitle()) ;

  // Create p.d.f.
  RooAbsPdf* ret =  new RooMultiVarGaussian(name.c_str(),title.c_str(),params3,mu1,Vred) ;
  ret->addOwnedComponents(mu1) ;  
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Change name of RooFitResult object

void RooFitResult::SetName(const char *name) 
{
  if (_dir) _dir->GetList()->Remove(this);
  TNamed::SetName(name) ;
  if (_dir) _dir->GetList()->Add(this);
}


////////////////////////////////////////////////////////////////////////////////
/// Change name and title of RooFitResult object

void RooFitResult::SetNameTitle(const char *name, const char* title) 
{
  if (_dir) _dir->GetList()->Remove(this);
  TNamed::SetNameTitle(name,title) ;
  if (_dir) _dir->GetList()->Add(this);
}


////////////////////////////////////////////////////////////////////////////////
/// Print name of fit result

void RooFitResult::printName(ostream& os) const 
{
  os << GetName() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print title of fit result

void RooFitResult::printTitle(ostream& os) const 
{
  os << GetTitle() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print class name of fit result

void RooFitResult::printClassName(ostream& os) const 
{
  os << IsA()->GetName() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print arguments of fit result, i.e. the parameters of the fit

void RooFitResult::printArgs(ostream& os) const 
{
  os << "[constPars=" << *_constPars << ",floatPars=" << *_finalPars << "]" ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print the value of the fit result, i.e.g the status, minimized FCN, edm and covariance quality code

void RooFitResult::printValue(ostream& os) const 
{
  os << "(status=" << _status << ",FCNmin=" << _minNLL << ",EDM=" << _edm << ",covQual=" << _covQual << ")" ;
}


////////////////////////////////////////////////////////////////////////////////
/// Configure default contents to be printed

Int_t RooFitResult::defaultPrintContents(Option_t* /*opt*/) const 
{
  return kName|kClassName|kArgs|kValue ;
}


////////////////////////////////////////////////////////////////////////////////
/// Configure mapping of Print() arguments to RooPrintable print styles

RooPrintable::StyleOption RooFitResult::defaultPrintStyle(Option_t* opt) const 
{
  if (!opt || strlen(opt)==0) {
    return kStandard ;
  }
  return RooPrintable::defaultPrintStyle(opt) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooFitResult.

void RooFitResult::Streamer(TBuffer &R__b)
{
  if (R__b.IsReading()) {
    UInt_t R__s, R__c;
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c);     
    if (R__v>3) {    
      R__b.ReadClassBuffer(RooFitResult::Class(),this,R__v,R__s,R__c);
      RooAbsArg::ioStreamerPass2Finalize();
      _corrMatrix.SetOwner();
    } else {
      // backward compatibitily streaming 
      TNamed::Streamer(R__b);
      RooPrintable::Streamer(R__b);
      RooDirItem::Streamer(R__b);
      R__b >> _status;
      R__b >> _covQual;
      R__b >> _numBadNLL;
      R__b >> _minNLL;
      R__b >> _edm;
      R__b >> _constPars;
      R__b >> _initPars;
      R__b >> _finalPars;
      R__b >> _globalCorr;
      _corrMatrix.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooFitResult::IsA());

      // Now fill new-style covariance and correlation matrix information
      // from legacy form
      _CM = new TMatrixDSym(_finalPars->getSize()) ;
      _VM = new TMatrixDSym(_CM->GetNcols()) ;
      _GC = new TVectorD(_CM->GetNcols()) ;
      
      TIterator *gcIter = _globalCorr->createIterator() ;
      TIterator *parIter = _finalPars->createIterator() ;
      RooRealVar* gcVal = 0;
      for (unsigned int i = 0; i < (unsigned int)_CM->GetNcols() ; ++i) {
	
	// Find the next global correlation slot to fill, skipping fixed parameters
	gcVal = (RooRealVar*) gcIter->Next() ;
	(*_GC)(i) = gcVal->getVal() ;
	
	// Fill a row of the correlation matrix
	TIterator* cIter = ((RooArgList*)_corrMatrix.At(i))->createIterator() ;
	for (unsigned int it = 0; it < (unsigned int)_CM->GetNcols() ; ++it) {
	  RooRealVar* cVal = (RooRealVar*) cIter->Next() ;	  
	  double value = cVal->getVal() ;
	  (*_CM)(it,i) = value ;
	  (*_CM)(i,it) = value;
	  (*_VM)(it,i) = value*((RooRealVar*)_finalPars->at(i))->getError()*((RooRealVar*)_finalPars->at(it))->getError() ;
	  (*_VM)(i,it) = (*_VM)(it,i) ;
	}
	delete cIter ;
      }
      
      delete gcIter ;
      delete parIter ;                 
    }

   } else {
      R__b.WriteClassBuffer(RooFitResult::Class(),this);
   }
}

