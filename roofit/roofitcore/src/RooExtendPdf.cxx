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

/** \class RooExtendPdf
RooExtendPdf is a wrapper around an existing PDF that adds a 
parameteric extended likelihood term to the PDF, optionally divided by a 
fractional term from a partial normalization of the PDF:
\f[
      n_\mathrm{Expected} = N \quad \text{or} \quad n_\mathrm{Expected} = N / \mathrm{frac} 
\f]
where \f$ N \f$ is supplied as a RooAbsReal to RooExtendPdf.
The fractional term is defined as
\f[
    \mathrm{frac} = \frac{\int_\mathrm{cutRegion[x]} \mathrm{pdf}(x,y) \; \mathrm{d}x \mathrm{d}y}{
      \int_\mathrm{normRegion[x]} \mathrm{pdf}(x,y) \; \mathrm{d}x \mathrm{d}y}
\f]

where \f$ x \f$ is the set of dependents involved in the selection region and \f$ y \f$
is the set of remaining dependents.

\f$ \mathrm{cutRegion}[x] \f$ is a limited integration range that is contained in
the nominal integration range \f$ \mathrm{normRegion}[x] \f$.
*/

#include "RooFit.h"
#include "Riostream.h"

#include "RooExtendPdf.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooNameReg.h"
#include "RooMsgService.h"



using namespace std;

ClassImp(RooExtendPdf);
;


RooExtendPdf::RooExtendPdf() : _rangeName(0)
{
  // Default constructor
}

/// Constructor. The ExtendPdf behaves identical to the supplied input pdf,
/// but adds an extended likelihood term. expectedEvents() will return 
/// `norm` if `rangeName` remains empty. If `rangeName` is not empty,
/// `norm` will refer to this range, and expectedEvents will return the
/// total number of events over the full range of the observables.
/// \param[in] name   Name of the pdf
/// \param[in] title  Title of the pdf (for plotting)
/// \param[in] pdf    The pdf to be extended
/// \param[in] norm   Expected number of events
/// \param[in] rangeName  If given, the number of events denoted by `norm` is interpreted as
/// the number of events in this range only
RooExtendPdf::RooExtendPdf(const char *name, const char *title, RooAbsPdf& pdf,
			   RooAbsReal& norm, const char* rangeName) :
  RooAbsPdf(name,title),
  _pdf("pdf", "PDF", this, pdf),
  _n("n","Normalization",this,norm),
  _rangeName(RooNameReg::ptr(rangeName))
{

  // Copy various setting from pdf
  setUnit(_pdf->getUnit()) ;
  setPlotLabel(_pdf->getPlotLabel()) ;
}



RooExtendPdf::RooExtendPdf(const RooExtendPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _pdf("pdf",this,other._pdf),
  _n("n",this,other._n),
  _rangeName(other._rangeName)
{
  // Copy constructor
}


RooExtendPdf::~RooExtendPdf() 
{
  // Destructor

}


/// Return the number of expected events over the full range of all variables.
/// `norm`, the variable set as normalisation constant in the constructor,
/// will yield the number of events in the range set in the constructor. That is, the function returns
/// \f[
///     N = \mathrm{norm} \; \cdot \; \frac{\int_{(x_F,y_F)} \mathrm{pdf}(x,y) }{\int_{(x_C,y_F)} \mathrm{pdf}(x,y)}
/// \f]
/// Where \f$ x \f$ is the set of dependents with a restricted range (defined by `rangeName` in the constructor),
/// and \f$ y \f$ are the other dependents. \f$ x_C \f$ is the integration
/// of \f$ x \f$ over the restricted range, and \f$ x_F \f$ is the integration of
/// \f$ x \f$ over the full range. `norm` is the number of events given as parameter to the constructor.
///
/// If the nested PDF can be extended, \f$ N \f$ is further scaled by its expected number of events.
Double_t RooExtendPdf::expectedEvents(const RooArgSet* nset) const 
{
  const RooAbsPdf& pdf = *_pdf;

  if (_rangeName && (!nset || nset->getSize()==0)) {
    coutW(InputArguments) << "RooExtendPdf::expectedEvents(" << GetName() << ") WARNING: RooExtendPdf needs non-null normalization set to calculate fraction in range " 
			  << _rangeName << ".  Results may be nonsensical" << endl ;  
  }

  Double_t nExp = _n ;

  // Optionally multiply with fractional normalization
  if (_rangeName) {

    double fracInt;
    {
      GlobalSelectComponentRAII globalSelComp(true);
      fracInt = pdf.getNormObj(nset,nset,_rangeName)->getVal();
    }


    if ( fracInt == 0. || _n == 0.) {
      coutW(Eval) << "RooExtendPdf(" << GetName() << ") WARNING: nExpected = " << _n << " / " 
		  << fracInt << " for nset = " << (nset?*nset:RooArgSet()) << endl ;
    }

    nExp /= fracInt ;    


    // cout << "RooExtendPdf::expectedEvents(" << GetName() << ") fracInt = " << fracInt << " _n = " << _n << " nExpect = " << nExp << endl ;

  }

  // Multiply with original Nexpected, if defined
  if (pdf.canBeExtended()) nExp *= pdf.expectedEvents(nset) ;

  return nExp ;
}



