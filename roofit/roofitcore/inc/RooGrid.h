/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooGrid.h,v 1.10 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_GRID
#define ROO_GRID

#include "TObject.h"
#include "RooPrintable.h"

class RooAbsFunc;

class RooGrid : public TObject, public RooPrintable {
public:
  RooGrid() ;
  RooGrid(const RooAbsFunc &function);
  virtual ~RooGrid();

  // Printing interface
  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const;

  inline virtual void Print(Option_t *options= 0) const {
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  inline Bool_t isValid() const { return _valid; }
  inline UInt_t getDimension() const { return _dim; }
  inline Double_t getVolume() const { return _vol; }
  inline UInt_t getNBins() const { return _bins; }
  inline UInt_t getNBoxes() const { return _boxes; }
  inline void setNBoxes(UInt_t boxes) { _boxes= boxes; }

  inline Double_t *createPoint() const { return _valid ? new Double_t[_dim] : 0; }
  inline UInt_t *createIndexVector() const { return _valid ? new UInt_t[_dim] : 0; }

  Bool_t initialize(const RooAbsFunc &function);
  void resize(UInt_t bins);
  void resetValues();
  void generatePoint(const UInt_t box[], Double_t x[], UInt_t bin[],
           Double_t &vol, Bool_t useQuasiRandom= kTRUE) const;
  void accumulate(const UInt_t bin[], Double_t amount);
  void refine(Double_t alpha= 1.5);

  void firstBox(UInt_t box[]) const;
  Bool_t nextBox(UInt_t box[]) const;

  enum { maxBins = 50 }; // must be even

  // Accessor for the j-th normalized grid point along the i-th dimension
public:
  inline Double_t coord(Int_t i, Int_t j) const { return _xi[i*_dim + j]; }
  inline Double_t value(Int_t i,Int_t j) const { return _d[i*_dim + j]; }
protected:
  inline Double_t& coord(Int_t i, Int_t j) { return _xi[i*_dim + j]; }
  inline Double_t& value(Int_t i,Int_t j) { return _d[i*_dim + j]; }
  inline Double_t& newCoord(Int_t i) { return _xin[i]; }

protected:

  Bool_t _valid;              // Is configuration valid
  UInt_t _dim,_bins,_boxes;   // Number of dimensions, bins and boxes
  Double_t _vol;              // Volume

  Double_t *_xl;     //! Internal workspace
  Double_t *_xu;     //! Internal workspace
  Double_t *_delx;   //! Internal workspace
  Double_t *_d;      //! Internal workspace
  Double_t *_xi;     //! Internal workspace
  Double_t *_xin;    //! Internal workspace
  Double_t *_weight; //! Internal workspace

  ClassDef(RooGrid,1) // Utility class for RooMCIntegrator holding a multi-dimensional grid
};

#endif


