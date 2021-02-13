/*****************************************************************************
 * Project: RooStats
 * Package: RooStats
 *    File: $Id$
 * author: Kyle Cranmer
 *****************************************************************************/
#ifndef JEFFREYSPRIOR
#define JEFFREYSPRIOR

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooJeffreysPrior : public RooAbsPdf {
public:

  RooJeffreysPrior() { };
  RooJeffreysPrior(const char *name, const char *title, RooAbsPdf& nominal, const RooArgList& paramSet, const RooArgList& obsSet) ;
  virtual ~RooJeffreysPrior() ;

  RooJeffreysPrior(const RooJeffreysPrior& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooJeffreysPrior(*this, newname); }

  const RooArgList& lowList() const { return _obsSet ; }
  const RooArgList& paramList() const { return _paramSet ; }

protected:

  RooTemplateProxy<RooAbsPdf> _nominal;    // Proxy to the PDF for this prior.
  RooListProxy _obsSet ;   // Observables of the PDF.
  RooListProxy _paramSet ; // Parameters of the PDF.

  Double_t evaluate() const;

private:
  struct CacheElem : public RooAbsCacheElement {
  public:
      virtual ~CacheElem() = default;
      // Payload
      std::unique_ptr<RooAbsPdf> _pdf;
      std::unique_ptr<RooArgSet> _pdfVariables;

      virtual RooArgList containedArgs(Action) override {
        RooArgList list(*_pdf);
        list.add(*_pdfVariables, true);
        return list;
      }
  };
  mutable RooObjCacheManager _cacheMgr; //!

  ClassDef(RooJeffreysPrior,2) // Sum of RooAbsReal objects
};

#endif
