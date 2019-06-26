/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

/** \class RooHistConstraint
 * \ingroup Roofit
 * The RooHistConstraint implements constraint terms for a binned PDF with statistical uncertainties.
 * Following the Barlow-Beeston method, it adds Poisson constraints for each bin that
 * constrain the statistical uncertainty of the template histogram.
 *
 * It can therefore be used to estimate the Monte Carlo uncertainty of a fit.
 *
 * Check also the tutorial rf709_BarlowBeeston.C
 *
 */

#include "Riostream.h"

#include "RooHistConstraint.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooParamHistFunc.h"
#include "RooRealVar.h"
#include <math.h>
#include "TMath.h"

using namespace std;

ClassImp(RooHistConstraint);

////////////////////////////////////////////////////////////////////////////////
/// Create a new RooHistConstraint.
/// \param[in] name Name of the PDF. This is used to identify it in a likelihood model.
/// \param[in] title Title for plotting etc.
/// \param[in] phfSet Set of parametrised histogram functions (RooParamHistFunc).
/// \param[in] threshold Threshold (bin content) up to which statistcal uncertainties are taken into account.
RooHistConstraint::RooHistConstraint(const char *name, const char *title,
    const RooArgSet& phfSet, Int_t threshold) :
  RooAbsPdf(name,title),
  _gamma("gamma","gamma",this),
  _nominal("nominal","nominal",this),
  _relParam(kTRUE)
{
  // Implementing constraint on sum of RooParamHists
  //
  // Step 1 - Create new gamma parameters for sum
  // Step 2 - Replace entries in gamma listproxy of components with new sum components
  // Step 3 - Implement constraints in terms of gamma sum parameters


  if (phfSet.getSize()==1) {

    RooParamHistFunc* phf = dynamic_cast<RooParamHistFunc*>(phfSet.first()) ;

    if (!phf) {
      coutE(InputArguments) << "RooHistConstraint::ctor(" << GetName()
                 << ") ERROR: input object must be a RooParamHistFunc" << endl ;
      throw std::string("RooHistConstraint::ctor ERROR incongruent input arguments") ;
    }

    // Now populate nominal with parameters
    RooArgSet allVars ;
    for (Int_t i=0 ; i<phf->_dh.numEntries() ; i++) {
      phf->_dh.get(i) ;
      if (phf->_dh.weight()<threshold && phf->_dh.weight() != 0.) {
        const char* vname = Form("%s_nominal_bin_%i",GetName(),i) ;
        RooRealVar* var = new RooRealVar(vname,vname,0,1.E30) ;
        var->setVal(phf->_dh.weight()) ;
        var->setConstant(true);
        allVars.add(*var) ;
        _nominal.add(*var) ;

        RooRealVar* gam = (RooRealVar*) phf->_p.at(i) ;
        if (var->getVal()>0) {
          gam->setConstant(false);
        }
        _gamma.add(*gam) ;
      }
    }

    addOwnedComponents(allVars) ;

    return ;
  }



  Int_t nbins(-1) ;
  vector<RooParamHistFunc*> phvec ;
  RooArgSet gammaSet ;
  string bin0_name ;
  for (const auto arg : phfSet) {

    RooParamHistFunc* phfComp = dynamic_cast<RooParamHistFunc*>(arg) ;
    if (phfComp) {
      phvec.push_back(phfComp) ;
      if (nbins==-1) {
        nbins = phfComp->_p.getSize() ;
        bin0_name = phfComp->_p.at(0)->GetName() ;
        gammaSet.add(phfComp->_p) ;
      } else {
        if (phfComp->_p.getSize()!=nbins) {
          coutE(InputArguments) << "RooHistConstraint::ctor(" << GetName()
                << ") ERROR: incongruent input arguments: all input RooParamHistFuncs should have same #bins" << endl ;
          throw std::string("RooHistConstraint::ctor ERROR incongruent input arguments") ;
        }
        if (bin0_name != phfComp->_p.at(0)->GetName()) {
          coutE(InputArguments) << "RooHistConstraint::ctor(" << GetName()
                << ") ERROR: incongruent input arguments: all input RooParamHistFuncs should have the same bin parameters.\n"
                << "Previously found " << bin0_name << ", now found " << phfComp->_p.at(0)->GetName() << ".\n"
                << "Check that the right RooParamHistFuncs have been passed to this RooHistConstraint." << std::endl;
          throw std::string("RooHistConstraint::ctor ERROR incongruent input arguments") ;
        }

      }
    } else {
      coutW(InputArguments) << "RooHistConstraint::ctor(" << GetName()
                 << ") WARNING: ignoring input argument " << arg->GetName() << " which is not of type RooParamHistFunc" << endl;
    }
  }

  _gamma.add(gammaSet) ;

  // Now populate nominal and nominalErr with parameters
  RooArgSet allVars ;
  for (Int_t i=0 ; i<nbins ; i++) {

    Double_t sumVal(0) ;
    for (const auto phfunc : phvec) {
      sumVal += phfunc->getNominal(i);
    }

    if (sumVal<threshold && sumVal != 0.) {

      const char* vname = Form("%s_nominal_bin_%i",GetName(),i) ;
      RooRealVar* var = new RooRealVar(vname,vname,0,1000) ;

      Double_t sumVal2(0) ;
      for (vector<RooParamHistFunc*>::iterator iter = phvec.begin() ; iter != phvec.end() ; ++iter) {
        sumVal2 += (*iter)->getNominal(i) ;
      }
      var->setVal(sumVal2) ;
      var->setConstant(kTRUE) ;

      vname = Form("%s_nominal_error_bin_%i",GetName(),i) ;
      RooRealVar* vare = new RooRealVar(vname,vname,0,1000) ;

      Double_t sumErr2(0) ;
      for (vector<RooParamHistFunc*>::iterator iter = phvec.begin() ; iter != phvec.end() ; ++iter) {
        sumErr2 += pow((*iter)->getNominalError(i),2) ;
      }
      vare->setVal(sqrt(sumErr2)) ;
      vare->setConstant(kTRUE) ;

      allVars.add(RooArgSet(*var,*vare)) ;
      _nominal.add(*var) ;
      //      _nominalErr.add(*vare) ;

      ((RooRealVar*)_gamma.at(i))->setConstant(kFALSE) ;

    }
  }
  addOwnedComponents(allVars) ;
}

////////////////////////////////////////////////////////////////////////////////

 RooHistConstraint::RooHistConstraint(const RooHistConstraint& other, const char* name) :
   RooAbsPdf(other,name),
   _gamma("gamma",this,other._gamma),
   _nominal("nominal",this,other._nominal),
   _relParam(other._relParam)
 {
 }

////////////////////////////////////////////////////////////////////////////////

 Double_t RooHistConstraint::evaluate() const
 {
   double prod(1);

   for (unsigned int i=0; i < _nominal.size(); ++i) {
     const auto& gamma = static_cast<const RooAbsReal&>(_gamma[i]);
     const auto& nominal = static_cast<const RooAbsReal&>(_nominal[i]);

     double gamVal = gamma.getVal();
     if (_relParam)
       gamVal *= nominal.getVal();

     const double pois = TMath::Poisson(nominal.getVal(),gamVal);
     prod *= pois;
   }

   return prod;
 }

////////////////////////////////////////////////////////////////////////////////

Double_t RooHistConstraint::getLogVal(const RooArgSet* /*set*/) const
{
   double sum = 0.;
   for (unsigned int i=0; i < _nominal.size(); ++i) {
     const auto& gamma = static_cast<const RooAbsReal&>(_gamma[i]);
     const auto& nominal = static_cast<const RooAbsReal&>(_nominal[i]);
     double gam = gamma.getVal();
     Int_t  nom = (Int_t) nominal.getVal();

     if (_relParam)
       gam *= nom;

     if (gam>0) {
       const double logPoisson = nom * log(gam) - gam - std::lgamma(nom+1);
       sum += logPoisson ;
     } else if (nom>0) {
       cerr << "ERROR in RooHistConstraint: gam=0 and nom>0" << endl ;
     }
   }

   return sum ;
}


