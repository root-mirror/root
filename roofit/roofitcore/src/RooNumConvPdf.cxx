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
\file RooNumConvPdf.cxx
\class RooNumConvPdf
\ingroup Roofitcore

Numeric 1-dimensional convolution operator PDF. This class can convolve any PDF
with any other PDF using a straightforward numeric calculation of the
convolution integral
This class should be used as last resort as numeric convolution calculated
this way is computationally intensive and prone to stability fitting problems. 
<b>The preferred way to compute numeric convolutions is RooFFTConvPdf</b>,
which calculates convolutions using Fourier Transforms (requires external free
FFTW3 package)
RooNumConvPdf implements reasonable defaults that should convolve most
functions reasonably well, but results strongly depend on the shape of your
input PDFS so always check your result.
The default integration engine for the numeric convolution is the
adaptive Gauss-Kronrod method, which empirically seems the most robust
for this task. You can override the convolution integration settings via
the RooNumIntConfig object reference returned by the convIntConfig() member
function
By default the numeric convolution is integrated from -infinity to
+infinity through a <pre>x -> 1/x</pre> coordinate transformation of the
tails. For convolution with a very small bandwidth it may be
advantageous (for both CPU consumption and stability) if the
integration domain is limited to a finite range. The function
setConvolutionWindow(mean,width,scale) allows to set a sliding
window around the x value to be calculated taking a RooAbsReal
expression for an offset and a width to be taken around the x
value. These input expression can be RooFormulaVars or other
function objects although the 3d 'scale' argument 'scale'
multiplies the width RooAbsReal expression given in the 2nd
argument, allowing for an appropriate window definition for most
cases without need for a RooFormulaVar object: e.g. a Gaussian
resolution PDF do setConvolutionWindow(gaussMean,gaussSigma,5)
Note that for a 'wide' Gaussian the -inf to +inf integration
may converge more quickly than that over a finite range!
The default numeric precision is 1e-7, i.e. the global default for
numeric integration but you should experiment with this value to
see if it is sufficient for example by studying the number of function
calls that MINUIT needs to fit your function as function of the
convolution precision. 
**/

#include "RooFit.h"

#include "Riostream.h"
#include "TH2F.h"
#include "RooNumConvPdf.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCustomizer.h"
#include "RooConvIntegrandBinding.h"
#include "RooNumIntFactory.h"
#include "RooGenContext.h"
#include "RooConvGenContext.h"



using namespace std;

ClassImp(RooNumConvPdf);




////////////////////////////////////////////////////////////////////////////////

RooNumConvPdf::RooNumConvPdf() :
  _init(kFALSE),
  _conv(0)
{  
}




//_____________________________________________________________________________R
RooNumConvPdf::RooNumConvPdf(const char *name, const char *title, RooRealVar& convVar, RooAbsPdf& inPdf, RooAbsPdf& resmodel) : 
  RooAbsPdf(name,title), 
  _init(kFALSE),
  _conv(0),
  _origVar("!origVar","Original Convolution variable",this,convVar),
  _origPdf("!origPdf","Original Input PDF",this,inPdf),
  _origModel("!origModel","Original Resolution model",this,resmodel)
{
  // Constructor of convolution operator PDF
  // 
  // convVar  :  convolution variable (on which both pdf and resmodel should depend)
  // pdf      :  input 'physics' pdf
  // resmodel :  input 'resultion' pdf
  //
  // output is pdf(x) (X) resmodel(x) = Int [ pdf(x') resmodel (x-x') ] dx'
  //
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooNumConvPdf::RooNumConvPdf(const RooNumConvPdf& other, const char* name) :
  RooAbsPdf(other,name), 
  _init(kFALSE),
  _origVar("!origVar",this,other._origVar),
  _origPdf("!origPdf",this,other._origPdf),
  _origModel("!origModel",this,other._origModel)
{
  // Make temporary clone of original convolution to preserve configuration information
  // This information will be propagated to a newly create convolution in a subsequent
  // call to initialize() 
  if (other._conv) {
    _conv = new RooNumConvolution(*other._conv,Form("%s_CONV",name?name:GetName())) ;
  } else {
    _conv = 0 ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooNumConvPdf::~RooNumConvPdf() 
{
  if (_init) {
    delete _conv ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate and return value of p.d.f

Double_t RooNumConvPdf::evaluate() const
{
  if (!_init) initialize() ;

  return _conv->evaluate() ;
}



////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of object

void RooNumConvPdf::initialize() const
{
  // Save pointer to any prototype convolution object (only present if this object is made through
  // a copy constructor) 
  RooNumConvolution* protoConv = _conv ;

  // Optionally pass along configuration data from prototype object
  _conv = new RooNumConvolution(Form("%s_CONV",GetName()),GetTitle(),var(),pdf(),model(),protoConv) ;

  // Delete prototype object now
  if (protoConv) {
    delete protoConv ;  
  }

  _init = kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return appropriate generator context for this convolved p.d.f. If both pdf and resolution
/// model support internal generation return and optimization convolution generation context
/// that uses a smearing algorithm. Otherwise return a standard accept/reject sampling
/// context on the convoluted shape.

RooAbsGenContext* RooNumConvPdf::genContext(const RooArgSet &vars, const RooDataSet *prototype, 
					    const RooArgSet* auxProto, Bool_t verbose) const 
{
  if (!_init) initialize() ;

  // Check if physics PDF and resolution model can both directly generate the convolution variable
  RooArgSet* modelDep = _conv->model().getObservables(&vars) ;
  modelDep->remove(_conv->var(),kTRUE,kTRUE) ;
  Int_t numAddDep = modelDep->getSize() ;
  delete modelDep ;

  RooArgSet dummy ;
  Bool_t pdfCanDir = (((RooAbsPdf&)_conv->pdf()).getGenerator(_conv->var(),dummy) != 0 && \
		      ((RooAbsPdf&)_conv->pdf()).isDirectGenSafe(_conv->var())) ;
  Bool_t resCanDir = (((RooAbsPdf&)_conv->model()).getGenerator(_conv->var(),dummy) !=0  && 
		      ((RooAbsPdf&)_conv->model()).isDirectGenSafe(_conv->var())) ;

  if (numAddDep>0 || !pdfCanDir || !resCanDir) {
    // Any resolution model with more dependents than the convolution variable
    // or pdf or resmodel do not support direct generation
    return new RooGenContext(*this,vars,prototype,auxProto,verbose) ;
  } 
  
  // Any other resolution model: use specialized generator context
  return new RooConvGenContext(*this,vars,prototype,auxProto,verbose) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooNumConvPdf to more intuitively reflect the contents of the
/// product operator construction

void RooNumConvPdf::printMetaArgs(ostream& os) const 
{
  os << _origPdf.arg().GetName() << "(" << _origVar.arg().GetName() << ") (*) " << _origModel.arg().GetName() << "(" << _origVar.arg().GetName() << ") " ;
}
