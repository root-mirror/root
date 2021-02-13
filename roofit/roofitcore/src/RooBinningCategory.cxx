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
\file RooBinningCategory.cxx
\class RooBinningCategory
\ingroup Roofitcore

Class RooBinningCategory provides a real-to-category mapping defined
by a series of thresholds. It evaluates the value of `inputVar` passed in the
constructor, and converts this into a bin number using a binning defined for
the inputVar. The name of this binning is passed in the constructor.
**/


#include "RooBinningCategory.h"

#include "RooFit.h"
#include "Riostream.h"
#include "RooStreamParser.h"

using namespace std;

ClassImp(RooBinningCategory);



////////////////////////////////////////////////////////////////////////////////
/// Constructor with input function to be mapped and name and index of default
/// output state of unmapped values

RooBinningCategory::RooBinningCategory(const char *name, const char *title, RooAbsRealLValue& inputVar, 
					   const char* binningName, const char* catTypeName) :
  RooAbsCategory(name, title), _inputVar("inputVar","Input category",this,inputVar), _bname(binningName)
{
  initialize(catTypeName) ;

}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooBinningCategory::RooBinningCategory(const RooBinningCategory& other, const char *name) :
  RooAbsCategory(other,name), _inputVar("inputVar",this,other._inputVar), _bname(other._bname)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooBinningCategory::~RooBinningCategory() 
{
}




////////////////////////////////////////////////////////////////////////////////
/// Iterator over all bins in input variable and define corresponding state labels

void RooBinningCategory::initialize(const char* catTypeName)
{
  const int nbins = _inputVar->getBinning(_bname.Length() > 0 ? _bname.Data() : nullptr).numBins();
  for (Int_t i=0 ; i<nbins ; i++) {
    string name = catTypeName!=0 ? Form("%s%d",catTypeName,i)
            : (_bname.Length()>0 ? Form("%s_%s_bin%d",_inputVar.arg().GetName(),_bname.Data(),i) 
            : Form("%s_bin%d",_inputVar.arg().GetName(),i)) ;
    defineState(name,i);
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Calculate and return the value of the mapping function

RooAbsCategory::value_type RooBinningCategory::evaluate() const
{
  Int_t ibin = _inputVar->getBin(_bname.Length() > 0 ? _bname.Data() : nullptr);

  if (!hasIndex(ibin)) {
    string name = (_bname.Length()>0) ? Form("%s_%s_bin%d",_inputVar.arg().GetName(),_bname.Data(),ibin) 
	                              : Form("%s_bin%d",_inputVar.arg().GetName(),ibin) ;
    const_cast<RooBinningCategory*>(this)->defineState(name,ibin);
  }

  return ibin;
}




////////////////////////////////////////////////////////////////////////////////
/// Print info about this threshold category to the specified stream. In addition to the info
/// from RooAbsCategory::printStream() we add:
///
///  Standard : input category
///     Shape : default value
///   Verbose : list of thresholds

void RooBinningCategory::printMultiline(ostream& os, Int_t content, Bool_t verbose, TString indent) const
{
   RooAbsCategory::printMultiline(os,content,verbose,indent);

   if (verbose) {
     os << indent << "--- RooBinningCategory ---" << endl
	<< indent << "  Maps from " ;
     _inputVar.arg().printStream(os,kName|kValue,kSingleLine);
   }
}


