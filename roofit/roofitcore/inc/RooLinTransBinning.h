/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooLinTransBinning.h,v 1.9 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_LIN_TRANS_BINNING
#define ROO_LIN_TRANS_BINNING

#include "Rtypes.h"
#include "RooAbsBinning.h"

class RooLinTransBinning : public RooAbsBinning {
public:

  RooLinTransBinning(const char* name=0) : RooAbsBinning(name) { }
  RooLinTransBinning(const RooAbsBinning& input, Double_t slope=1.0, Double_t offset=0.0, const char* name=0);
  RooLinTransBinning(const RooLinTransBinning&, const char* name=0);
  virtual RooAbsBinning* clone(const char* name=0) const { return new RooLinTransBinning(*this,name) ; }
  virtual ~RooLinTransBinning() ;

  virtual Int_t numBoundaries() const { return _input->numBoundaries() ; }
  virtual Int_t binNumber(Double_t x) const { return _input->binNumber(invTrans(x)) ; }
  virtual Double_t binCenter(Int_t bin) const { return trans(_input->binCenter(binTrans(bin))) ; }
  virtual Double_t binWidth(Int_t bin) const { return _slope*_input->binWidth(binTrans(bin)) ; }
  virtual Double_t binLow(Int_t bin) const { if (_slope>0) return trans(_input->binLow(binTrans(bin))) ; else return trans(_input->binHigh(binTrans(bin))) ; }
  virtual Double_t binHigh(Int_t bin) const { if (_slope>0) return trans(_input->binHigh(binTrans(bin))) ; else return trans(_input->binLow(binTrans(bin))) ; }

  virtual void setRange(Double_t xlo, Double_t xhi) ;
  virtual void setMin(Double_t xlo) { setRange(xlo,highBound()) ; }
  virtual void setMax(Double_t xhi) { setRange(lowBound(),xhi) ; }

  virtual Double_t lowBound() const { if (_slope>0) return trans(_input->lowBound()) ; else return trans(_input->highBound()) ; }
  virtual Double_t highBound() const { if (_slope>0) return trans(_input->highBound()) ; else return trans(_input->lowBound()) ; }
  virtual Double_t averageBinWidth() const { return _slope*_input->averageBinWidth() ; }

  virtual Double_t* array() const ;

  void updateInput(const RooAbsBinning& input, Double_t slope=1.0, Double_t offset=0.0) ;

protected:

  inline Int_t binTrans(Int_t bin) const { if (_slope>0) return bin ; else return numBins()-bin-1 ; }
  inline Double_t trans(Double_t x) const { return x*_slope + _offset ; }
  inline Double_t invTrans(Double_t x) const { if (_slope==0.) return 0 ; return (x-_offset)/_slope ; }

  Double_t _slope{0.};   // Slope of transformation
  Double_t _offset{0.};  // Offset of transformation
  RooAbsBinning* _input{nullptr};    // Input binning
  mutable Double_t *_array{nullptr}; //! Array of transformed bin boundaries

  ClassDef(RooLinTransBinning,1) // Linear transformation of binning specification
};

#endif
