/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   AH, Andreas Hoecker,  Orsay,            hoecker@slac.stanford.edu       *
 *   SL, Sandrine Laplace, Orsay,            laplace@slac.stanford.edu       *
 *   JS, Jan Stark,        Paris,            stark@slac.stanford.edu         *
 *   WV, Wouter Verkerke,  UC Santa Barbara, verkerke@slac.stanford.edu      *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California,         *
 *                          IN2P3. All rights reserved.                      *
 *                                                                           *
 * History                                                                   *
 *   Nov-2001   WV Created initial version                                   *
 *   Dec-2001   SL mischarge correction, direct CPV                          *
 *   Jan-2002   AH built dedicated generator + code cleaning                 *
 *   Mar-2002   JS committed debugged version to CVS                         *
 *   Apr-2002   AH allow prompt (ie, non-Pdf) mischarge treatment            *
 *   May-2002   JS Changed the set of CP parameters (mathematically equiv.)  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooNonCPEigenDecay
    \ingroup Roofit

Time-dependent RooAbsAnaConvPdf for CP violating decays
to Non-CP eigenstates (eg, \f$ B_0 \rightarrow \rho^\pm \pi^\mp\f$).
For a description of the physics model see the
BaBar Physics Book, section 6.5.2.3 .
The set of CP parameters used in this class is equivalent to
the one used in the Physics Book, but it is not exactly the
same. Starting from the set in the BaBar Book, in order to
get the parameters used here you have to change the sign of both
\f$a_c^+\f$ and \f$a_c^-\f$, and then substitute:
\f[
   a_s^Q = S + Q \cdot \delta S \\
   a_c^Q = C + Q \cdot \delta C
\f]
where Q denotes the charge of the \f$\rho\f$ meson.
**/

#include "RooFit.h"

#include "Riostream.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooNonCPEigenDecay.h"
#include "TMath.h"
#include "RooRealIntegral.h"

using namespace std;

ClassImp(RooNonCPEigenDecay);

#define Debug_RooNonCPEigenDecay 1


////////////////////////////////////////////////////////////////////////////////

RooNonCPEigenDecay::RooNonCPEigenDecay( const char *name, const char *title,
               RooRealVar&     t,
               RooAbsCategory& tag,
               RooAbsReal&     tau,
               RooAbsReal&     dm,
               RooAbsReal&     avgW,
               RooAbsReal&     delW,
               RooAbsCategory& rhoQ,
               RooAbsReal&     correctQ,
               RooAbsReal&     wQ,
               RooAbsReal&     acp,
               RooAbsReal&     C,
               RooAbsReal&     delC,
               RooAbsReal&     S,
               RooAbsReal&     delS,
               const RooResolutionModel& model,
               DecayType       type )
  : RooAbsAnaConvPdf( name, title, model, t ),
  _acp      ( "acp",      "acp",                this, acp      ),
  _avgC        ( "C",        "C",                  this, C        ),
  _delC     ( "delC",     "delC",               this, delC     ),
  _avgS        ( "S",        "S",                  this, S        ),
  _delS     ( "delS",     "delS",               this, delS     ),
  _avgW     ( "avgW",     "Average mistag rate",this, avgW     ),
  _delW     ( "delW",     "Shift mistag rate",  this, delW     ),
  _t        ( "t",        "time",               this, t        ),
  _tau      ( "tau",      "decay time",         this, tau      ),
  _dm       ( "dm",       "mixing frequency",   this, dm       ),
  _tag      ( "tag",      "CP state",           this, tag      ),
  _rhoQ     ( "rhoQ",     "Charge of the rho",  this, rhoQ     ),
  _correctQ ( "correctQ", "correction of rhoQ", this, correctQ ),
  _wQ       ( "wQ",       "mischarge",          this, wQ       ),
  _genB0Frac     ( 0 ),
  _genRhoPlusFrac( 0 ),
  _type     ( type )
{
  // Constructor
  switch(type) {
  case SingleSided:
    _basisExp = declareBasis( "exp(-@0/@1)",            RooArgList( tau     ) );
    _basisSin = declareBasis( "exp(-@0/@1)*sin(@0*@2)", RooArgList( tau, dm ) );
    _basisCos = declareBasis( "exp(-@0/@1)*cos(@0*@2)", RooArgList( tau, dm ) );
    break;
  case Flipped:
    _basisExp = declareBasis( "exp(@0)/@1)",            RooArgList( tau     ) );
    _basisSin = declareBasis( "exp(@0/@1)*sin(@0*@2)",  RooArgList( tau, dm ) );
    _basisCos = declareBasis( "exp(@0/@1)*cos(@0*@2)",  RooArgList( tau, dm ) );
    break;
  case DoubleSided:
    _basisExp = declareBasis( "exp(-abs(@0)/@1)",            RooArgList( tau     ) );
    _basisSin = declareBasis( "exp(-abs(@0)/@1)*sin(@0*@2)", RooArgList( tau, dm ) );
    _basisCos = declareBasis( "exp(-abs(@0)/@1)*cos(@0*@2)", RooArgList( tau, dm ) );
    break;
  }
}

////////////////////////////////////////////////////////////////////////////////

RooNonCPEigenDecay::RooNonCPEigenDecay( const char *name, const char *title,
               RooRealVar&     t,
               RooAbsCategory& tag,
               RooAbsReal&     tau,
               RooAbsReal&     dm,
               RooAbsReal&     avgW,
               RooAbsReal&     delW,
               RooAbsCategory& rhoQ,
               RooAbsReal&     correctQ,
               RooAbsReal&     acp,
               RooAbsReal&     C,
               RooAbsReal&     delC,
               RooAbsReal&     S,
               RooAbsReal&     delS,
               const RooResolutionModel& model,
               DecayType       type )
  : RooAbsAnaConvPdf( name, title, model, t ),
  _acp      ( "acp",      "acp",                this, acp      ),
  _avgC        ( "C",        "C",                  this, C        ),
  _delC     ( "delC",     "delC",               this, delC     ),
  _avgS        ( "S",        "S",                  this, S        ),
  _delS     ( "delS",     "delS",               this, delS     ),
  _avgW     ( "avgW",     "Average mistag rate",this, avgW     ),
  _delW     ( "delW",     "Shift mistag rate",  this, delW     ),
  _t        ( "t",        "time",               this, t        ),
  _tau      ( "tau",      "decay time",         this, tau      ),
  _dm       ( "dm",       "mixing frequency",   this, dm       ),
  _tag      ( "tag",      "CP state",           this, tag      ),
  _rhoQ     ( "rhoQ",     "Charge of the rho",  this, rhoQ     ),
  _correctQ ( "correctQ", "correction of rhoQ", this, correctQ ),
  _genB0Frac     ( 0 ),
  _genRhoPlusFrac( 0 ),
  _type     ( type )
{
  // dummy mischarge (must be set to zero!)
  _wQ = RooRealProxy( "wQ", "mischarge", this, *(new RooRealVar( "wQ", "wQ", 0 )) );

  switch(type) {
  case SingleSided:
    _basisExp = declareBasis( "exp(-@0/@1)",            RooArgList( tau     ) );
    _basisSin = declareBasis( "exp(-@0/@1)*sin(@0*@2)", RooArgList( tau, dm ) );
    _basisCos = declareBasis( "exp(-@0/@1)*cos(@0*@2)", RooArgList( tau, dm ) );
    break;
  case Flipped:
    _basisExp = declareBasis( "exp(@0)/@1)",            RooArgList( tau     ) );
    _basisSin = declareBasis( "exp(@0/@1)*sin(@0*@2)",  RooArgList( tau, dm ) );
    _basisCos = declareBasis( "exp(@0/@1)*cos(@0*@2)",  RooArgList( tau, dm ) );
    break;
  case DoubleSided:
    _basisExp = declareBasis( "exp(-abs(@0)/@1)",            RooArgList( tau     ) );
    _basisSin = declareBasis( "exp(-abs(@0)/@1)*sin(@0*@2)", RooArgList( tau, dm ) );
    _basisCos = declareBasis( "exp(-abs(@0)/@1)*cos(@0*@2)", RooArgList( tau, dm ) );
    break;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooNonCPEigenDecay::RooNonCPEigenDecay( const RooNonCPEigenDecay& other, const char* name )
  : RooAbsAnaConvPdf( other, name ),
  _acp      ( "acp",      this, other._acp      ),
  _avgC        ( "C",        this, other._avgC        ),
  _delC     ( "delC",     this, other._delC     ),
  _avgS        ( "S",        this, other._avgS        ),
  _delS     ( "delS",     this, other._delS     ),
  _avgW     ( "avgW",     this, other._avgW     ),
  _delW     ( "delW",     this, other._delW     ),
  _t        ( "t",        this, other._t        ),
  _tau      ( "tau",      this, other._tau      ),
  _dm       ( "dm",       this, other._dm       ),
  _tag      ( "tag",      this, other._tag      ),
  _rhoQ     ( "rhoQ",     this, other._rhoQ     ),
  _correctQ ( "correctQ", this, other._correctQ ),
  _wQ       ( "wQ",       this, other._wQ       ),
  _genB0Frac     ( other._genB0Frac      ),
  _genRhoPlusFrac( other._genRhoPlusFrac ),
  _type          ( other._type           ),
  _basisExp      ( other._basisExp       ),
  _basisSin      ( other._basisSin       ),
  _basisCos      ( other._basisCos       )
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooNonCPEigenDecay::~RooNonCPEigenDecay( void )
{
}

////////////////////////////////////////////////////////////////////////////////
///  - B0    : _tag  == -1
///  - B0bar : _tag  == +1
///  - rho+  : _rhoQ == +1
///  - rho-  : _rhoQ == -1
///  - the charge correction factor "_correctQ" serves to implement mis-charges

Double_t RooNonCPEigenDecay::coefficient( Int_t basisIndex ) const
{
  Int_t rhoQc = _rhoQ * int(_correctQ);
  assert( rhoQc == 1 || rhoQc == -1 );

  Double_t a_sin_p = _avgS + _delS;
  Double_t a_sin_m = _avgS - _delS;
  Double_t a_cos_p = _avgC + _delC;
  Double_t a_cos_m = _avgC - _delC;

  if (basisIndex == _basisExp) {
    if (rhoQc == -1 || rhoQc == +1)
      return (1 + rhoQc*_acp*(1 - 2*_wQ))*(1 + 0.5*_tag*(2*_delW));
    else
      return 1;
  }

  if (basisIndex == _basisSin) {

    if (rhoQc == -1)

      return - ((1 - _acp)*a_sin_m*(1 - _wQ) + (1 + _acp)*a_sin_p*_wQ)*(1 - 2*_avgW)*_tag;

    else if (rhoQc == +1)

      return - ((1 + _acp)*a_sin_p*(1 - _wQ) + (1 - _acp)*a_sin_m*_wQ)*(1 - 2*_avgW)*_tag;

    else
       return - _tag*((a_sin_p + a_sin_m)/2)*(1 - 2*_avgW);
  }

  if (basisIndex == _basisCos) {

    if ( rhoQc == -1)

      return + ((1 - _acp)*a_cos_m*(1 - _wQ) + (1 + _acp)*a_cos_p*_wQ)*(1 - 2*_avgW)*_tag;

    else if (rhoQc == +1)

      return + ((1 + _acp)*a_cos_p*(1 - _wQ) + (1 - _acp)*a_cos_m*_wQ)*(1 - 2*_avgW)*_tag;

    else
      return _tag*((a_cos_p + a_cos_m)/2)*(1 - 2*_avgW);
  }

  return 0;
}

// advertise analytical integration

////////////////////////////////////////////////////////////////////////////////

Int_t RooNonCPEigenDecay::getCoefAnalyticalIntegral( Int_t /*code*/, RooArgSet& allVars,
                       RooArgSet& analVars, const char* rangeName ) const
{
  if (rangeName) return 0 ;

  if (matchArgs( allVars, analVars, _tag, _rhoQ )) return 3;
  if (matchArgs( allVars, analVars, _rhoQ       )) return 2;
  if (matchArgs( allVars, analVars, _tag        )) return 1;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// correct for the right/wrong charge...

Double_t RooNonCPEigenDecay::coefAnalyticalIntegral( Int_t basisIndex,
                       Int_t code, const char* /*rangeName*/ ) const
{
  Int_t rhoQc = _rhoQ*int(_correctQ);

  Double_t a_sin_p = _avgS + _delS;
  Double_t a_sin_m = _avgS - _delS;
  Double_t a_cos_p = _avgC + _delC;
  Double_t a_cos_m = _avgC - _delC;

  switch(code) {

    // No integration
  case 0: return coefficient(basisIndex);

    // Integration over 'tag'
  case 1:
    if (basisIndex == _basisExp) return 2*(1 + rhoQc*_acp*(1 - 2*_wQ));
    if (basisIndex == _basisSin || basisIndex==_basisCos) return 0;
    assert( kFALSE );

    // Integration over 'rhoQ'
  case 2:
    if (basisIndex == _basisExp) return 2*(1 + 0.5*_tag*(2.*_delW));

    if (basisIndex == _basisSin)

      return - ( (1 - _acp)*a_sin_m + (1 + _acp)*a_sin_p )*(1 - 2*_avgW)*_tag;

    if (basisIndex == _basisCos)

      return + ( (1 - _acp)*a_cos_m + (1 + _acp)*a_cos_p )*(1 - 2*_avgW)*_tag;

    assert( kFALSE );

    // Integration over 'tag' and 'rhoQ'
  case 3:
    if (basisIndex == _basisExp) return 2*2; // for both: tag and charge
    if (basisIndex == _basisSin || basisIndex==_basisCos) return 0;
    assert( kFALSE );

  default:
    assert( kFALSE );
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooNonCPEigenDecay::getGenerator( const RooArgSet& directVars,
               RooArgSet&       generateVars, Bool_t staticInitOK ) const
{
  if (staticInitOK) {
    if (matchArgs( directVars, generateVars, _t, _tag, _rhoQ )) return 4;
    if (matchArgs( directVars, generateVars, _t, _rhoQ       )) return 3;
    if (matchArgs( directVars, generateVars, _t, _tag        )) return 2;
  }
  if (matchArgs( directVars, generateVars, _t              )) return 1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void RooNonCPEigenDecay::initGenerator( Int_t code )
{
  if (code == 2 || code == 4) {
    // Calculate the fraction of mixed events to generate
    Double_t sumInt1 = RooRealIntegral( "sumInt1", "sum integral1", *this,
               RooArgSet( _t.arg(), _tag.arg(), _rhoQ.arg() )
                  ).getVal();
    _tag = -1;
    Double_t b0Int1 = RooRealIntegral( "mixInt1", "mix integral1", *this,
                   RooArgSet( _t.arg(), _rhoQ.arg() )
                 ).getVal();
    _genB0Frac = b0Int1/sumInt1;

    if (Debug_RooNonCPEigenDecay == 1)
      cout << "     o RooNonCPEigenDecay::initgenerator: genB0Frac     : "
      << _genB0Frac
      << ", tag dilution: " << (1 - 2*_avgW)
      << endl;
  }

  if (code == 3 || code == 4) {
    // Calculate the fraction of positive rho's to generate
    Double_t sumInt2 = RooRealIntegral( "sumInt2", "sum integral2", *this,
               RooArgSet( _t.arg(), _tag.arg(), _rhoQ.arg() )
                  ).getVal();
    _rhoQ = 1;
    Double_t b0Int2 = RooRealIntegral( "mixInt2", "mix integral2", *this,
                   RooArgSet( _t.arg(), _tag.arg() )
                 ).getVal();
    _genRhoPlusFrac = b0Int2/sumInt2;

    if (Debug_RooNonCPEigenDecay == 1)
      cout << "     o RooNonCPEigenDecay::initgenerator: genRhoPlusFrac: "
      << _genRhoPlusFrac << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooNonCPEigenDecay::generateEvent( Int_t code )
{
  // Generate delta-t dependent
  while (kTRUE) {

    // B flavor and rho charge (we do not use the integrated weights)
    if (code != 1) {
      if (code != 3) _tag  = (RooRandom::uniform()<=0.5) ? -1 : +1;
      if (code != 2) _rhoQ = (RooRandom::uniform()<=0.5) ?  1 : -1;
    }

    // opposite charge?
    // Int_t rhoQc = _rhoQ*int(_correctQ);

    Double_t a_sin_p = _avgS + _delS;
    Double_t a_sin_m = _avgS - _delS;
    Double_t a_cos_p = _avgC + _delC;
    Double_t a_cos_m = _avgC - _delC;

    // maximum probability density
    double a1 = 1 + sqrt(TMath::Power(a_cos_m, 2) + TMath::Power(a_sin_m, 2));
    double a2 = 1 + sqrt(TMath::Power(a_cos_p, 2) + TMath::Power(a_sin_p, 2));

    Double_t maxAcceptProb = (1.10 + TMath::Abs(_acp)) * (a1 > a2 ? a1 : a2);
    // The 1.10 in the above line is a security feature to prevent crashes close to the limit at 1.00

    Double_t rand = RooRandom::uniform();
    Double_t tval(0);

    switch(_type) {

    case SingleSided:
      tval = -_tau*log(rand);
      break;

    case Flipped:
      tval = +_tau*log(rand);
      break;

    case DoubleSided:
      tval = (rand<=0.5) ? -_tau*log(2*rand) : +_tau*log(2*(rand-0.5));
      break;
    }

    // get coefficients
    Double_t expC = coefficient( _basisExp );
    Double_t sinC = coefficient( _basisSin );
    Double_t cosC = coefficient( _basisCos );

    // probability density
    Double_t acceptProb  = expC + sinC*sin(_dm*tval) + cosC*cos(_dm*tval);

    // sanity check...
    assert( acceptProb <= maxAcceptProb );

    // hit or miss...
    Bool_t accept = maxAcceptProb*RooRandom::uniform() < acceptProb ? kTRUE : kFALSE;

    if (accept && tval<_t.max() && tval>_t.min()) {
      _t = tval;
      break;
    }
  }
}
