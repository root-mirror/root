// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_LININTERPVAR
#define ROOSTATS_LININTERPVAR

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include <vector>

class RooRealVar;
class RooArgList ;

namespace RooStats{
namespace HistFactory{

  class LinInterpVar : public RooAbsReal {
  public:

    LinInterpVar() ;
    LinInterpVar(const char *name, const char *title,
                 const RooArgList& _paramList, double nominal, std::vector<double> low, std::vector<double> high);

    LinInterpVar(const char *name, const char *title);
    LinInterpVar(const LinInterpVar&, const char*);

    virtual TObject* clone(const char* newname) const { return new LinInterpVar(*this, newname); }
    virtual ~LinInterpVar() ;


  protected:

    RooListProxy _paramList ;
    double _nominal;
    std::vector<double> _low;
    std::vector<double> _high;
    
    TIterator* _paramIter ;  //! do not persist

    Double_t evaluate() const;

    ClassDef(RooStats::HistFactory::LinInterpVar,1) // Piecewise linear interpolation
  };
}
}

#endif
