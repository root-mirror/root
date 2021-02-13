/** \class RooJeffreysPrior
\ingroup Roofit

Implementation of Jeffrey's prior. This class estimates the fisher information matrix by generating
a binned Asimov dataset from the supplied PDFs, fitting it, retrieving the covariance matrix and inverting
it. It returns the square root of the determinant of this matrix.
Numerical integration is used to normalise. Since each integration step requires fits to be run,
evaluating complicated PDFs may take long.

Check the tutorial rs302_JeffreysPriorDemo.C for a demonstration with a simple PDF.
**/

#include "RooJeffreysPrior.h"

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "TMatrixDSym.h"
#include "RooDataHist.h"
#include "RooNumIntConfig.h"
#include "RooRealVar.h"
#include "RooHelpers.h"

using namespace std;

ClassImp(RooJeffreysPrior);

using namespace RooFit;

////////////////////////////////////////////////////////////////////////////////
/// Construct a new JeffreysPrior.
/// \param[in] name     Name of this object.
/// \param[in] title    Title (for plotting)
/// \param[in] nominal  The PDF to base Jeffrey's prior on.
/// \param[in] paramSet Parameters of the PDF.
/// \param[in] obsSet   Observables of the PDF.

RooJeffreysPrior::RooJeffreysPrior(const char* name, const char* title,
              RooAbsPdf& nominal,
              const RooArgList& paramSet,
              const RooArgList& obsSet) :
  RooAbsPdf(name, title),
  _nominal("nominal","nominal",this, nominal, false, false),
  _obsSet("!obsSet","Observables",this, false, false),
  _paramSet("!paramSet","Parameters",this),
  _cacheMgr(this, 1, true, false)
{
  for (const auto comp : obsSet) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "RooJeffreysPrior::ctor(" << GetName() << ") ERROR: component " << comp->GetName()
             << " in observable list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _obsSet.add(*comp) ;
  }

  for (const auto comp : paramSet) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "RooJeffreysPrior::ctor(" << GetName() << ") ERROR: component " << comp->GetName()
             << " in parameter list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _paramSet.add(*comp) ;
  }

  // use a different integrator by default.
  if(paramSet.getSize()==1)
    this->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")  ;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooJeffreysPrior::RooJeffreysPrior(const RooJeffreysPrior& other, const char* name) :
  RooAbsPdf(other, name),
  _nominal("!nominal",this,other._nominal),
  _obsSet("!obsSet",this,other._obsSet),
  _paramSet("!paramSet",this,other._paramSet),
  _cacheMgr(this, 1, true, false)
{

}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooJeffreysPrior::~RooJeffreysPrior()
{

}

////////////////////////////////////////////////////////////////////////////////
/// Calculate and return current value of self

Double_t RooJeffreysPrior::evaluate() const
{
  RooHelpers::LocalChangeMsgLevel msgLvlRAII(RooFit::WARNING);


  CacheElem* cacheElm = (CacheElem*) _cacheMgr.getObj(nullptr);
  if (!cacheElm) {
    //Internally, we have to enlarge the range of fit parameters to make
    //fits converge even if we are close to the limit of a parameter. Therefore, we clone the pdf and its
    //observables here. If something happens to the external PDF, the cache is wiped,
    //and we start to clone again.
    auto& pdf = _nominal.arg();
    RooAbsPdf* clonePdf = static_cast<RooAbsPdf*>(pdf.cloneTree());
    auto vars = clonePdf->getParameters(_obsSet);
    for (auto varTmp : *vars) {
      auto& var = static_cast<RooRealVar&>(*varTmp);
      auto range = var.getRange();
      double span = range.second - range.first;
      var.setRange(range.first - 0.1*span, range.second + 0.1 * span);
    }

    cacheElm = new CacheElem;
    cacheElm->_pdf.reset(clonePdf);
    cacheElm->_pdfVariables.reset(vars);

    _cacheMgr.setObj(nullptr, cacheElm);
  }

  auto& cachedPdf = *cacheElm->_pdf;
  auto& pdfVars = *cacheElm->_pdfVariables;
  pdfVars = _paramSet;

  std::unique_ptr<RooDataHist> data( cachedPdf.generateBinned(_obsSet,ExpectedData()) );
  std::unique_ptr<RooFitResult> res( cachedPdf.fitTo(*data, Save(),PrintLevel(-1),Minos(false),SumW2Error(false)) );
  TMatrixDSym cov = res->covarianceMatrix();
  cov.Invert();

  return sqrt(cov.Determinant());
}

