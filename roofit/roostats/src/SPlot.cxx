// @(#)root/roostats:$Id$
// Author: Kyle Cranmer   28/07/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*****************************************************************************
 * Project: RooStats
 * Package: RooFit/RooStats
 *
 * Authors:
 *   Original code from M. Pivk as part of MLFit package from BaBar.
 *   Modifications:
 *     Giacinto Piacquadio, Maurizio Pierini: modifications for new RooFit version
 *     George H. Lewis, Kyle Cranmer: generalized for weighted events
 *
 * Porting to RooStats (with permission) by Kyle Cranmer, July 2008
 *
 *****************************************************************************/


/** \class RooStats::SPlot
    \ingroup Roostats

  A class to calculate "sWeights" used to create an "sPlot".
  An sPlot can reweight a dataset to show different components (e.g. signal / background),
  but it doesn't use cuts, and therefore doesn't have to sort events into signal/background (or other) categories.
  Instead of *assigning* a category to each event in the dataset, all events are *weighted*.
  To compute the weights, a PDF with different components is analysed, and the weights are added
  to the dataset. When plotting the dataset with the weights of the signal or
  background components, the data looks like "signal", but all events in the dataset are used.

  The result is similar to a likelihood projection plot, but without cuts.

  \note SPlot needs to fit the pdf to the data once, so make sure that all relevant fit arguments such as
  the fit range are passed in the constructor.

  The code is based on
  ``SPlot: A statistical tool to unfold data distributions,''
  Nucl. Instrum. Meth. A 555, 356 (2005) [arXiv:physics/0402083].

  ### Creating an SPlot
  To use this class, you first must have a pdf that includes
  yield parameters for (possibly several) different species, for example a signal and background
  yield. Those yields must be of type RooRealVar / RooLinearVar (or anything that derives from
  RooAbsRealLValue). This is necessary because
  RooStats needs to be able to set the yields to 0 and 1 to probe the PDF. After
  constructing the s weights, the yields will be restored to their original values.

  To create an instance of the SPlot, supply a data set, the pdf to analyse,
  and a list which parameters of the pdf are yields. The SPlot will calculate SWeights, and
  include these as columns in the RooDataSet. The dataset will have two additional columns
  for every yield with name "`<varname>`":
  - `L_<varname>` is the the likelihood for each event, *i.e.*, the pdf evaluated for the given value of the variable "varname".
  - `<varname>_sw` is the value of the sWeight for the variable "varname" for each event.

  In SPlot::SPlot(), one can choose whether columns should be added to an existing dataset or whether a copy of the dataset
  should be created.

  ### Plotting s-weighted data
  After computing the s weights, create a new dataset that uses the s weights of the variable of interest for weighting.
  If the yield parameter for signal was e.g. "signalYield", the dataset can be constructed as follows:
  ~~~{.cpp}
  RooDataSet data_signal("<name>", "<title>", <dataWithSWeights>, <variables>, 0, "signalYield_sw");
  ~~~

  A complete tutorial with an extensive model is rs301_splot.C

  #### Using ratios as yield parameters
  As mentioned, RooStats needs to be able to modify the yield parameters. That means that they have to be a RooRealVar
  of a RooLinearVar. This allows using ratio parameters as in the following example:
  ~~~{.cpp}
  RooRealVar x("x", "observable", 0, 0, 20);
  RooRealVar m("m", "mean", 5., -10, 10);
  RooRealVar s("s", "sigma", 2., 0, 10);
  RooGaussian gaus("gaus", "gaus", x, m, s);

  RooRealVar a("a", "exp", -0.2, -10., 0.);
  RooExponential ex("ex", "ex", x, a);

  RooRealVar common("common", "common scale", 3., 0, 10);
  RooRealVar r1("r1", "ratio of signal events", 0.3, 0, 10);
  RooRealVar r2("r2", "ratio of background events", 0.5, 0, 10);
  RooLinearVar c1("c1", "c1", r1, common, RooFit::RooConst(0.));
  RooLinearVar c2("c2", "c2", r2, common, RooFit::RooConst(0.));

  RooAddPdf sum("sum", "sum", RooArgSet(gaus, ex), RooArgSet(c1, c2));
  auto data = sum.generate(x, 1000);

  RooStats::SPlot splot("splot", "splot", *data, &sum, RooArgSet(c1, c2));
  ~~~
*/

#include <vector>
#include <map>

#include "RooStats/SPlot.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"
#include "RooStats/RooStatsUtils.h"


#include "TMatrixD.h"


ClassImp(RooStats::SPlot); ;

using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

SPlot::~SPlot()
{
   if(TestBit(kOwnData) && fSData)
      delete fSData;

}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

SPlot::SPlot():
  TNamed()
{
  RooArgList Args;

  fSWeightVars = Args;

  fSData = NULL;

}

////////////////////////////////////////////////////////////////////////////////

SPlot::SPlot(const char* name, const char* title):
  TNamed(name, title)
{
  RooArgList Args;

  fSWeightVars = Args;

  fSData = NULL;

}

////////////////////////////////////////////////////////////////////////////////
///Constructor from a RooDataSet
///No sWeighted variables are present

SPlot::SPlot(const char* name, const char* title, const RooDataSet &data):
  TNamed(name, title)
{
  RooArgList Args;

  fSWeightVars = Args;

  fSData = (RooDataSet*) &data;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor from another SPlot

SPlot::SPlot(const SPlot &other):
  TNamed(other)
{
  RooArgList Args = (RooArgList) other.GetSWeightVars();

  fSWeightVars.addClone(Args);

  fSData = (RooDataSet*) other.GetSDataSet();

}

////////////////////////////////////////////////////////////////////////////////
///Construct a new SPlot instance, calculate sWeights, and include them
///in the RooDataSet held by this instance.
///
/// The constructor automatically calls AddSWeight() to add s weights to the dataset.
/// These can be retrieved later using GetSWeight() or GetSDataSet().
///\param[in] name Name of the instance.
///\param[in] title Title of the instance.
///\param[in] data Dataset to fit to.
///\param[in] pdf PDF to compute s weights for.
///\param[in] yieldsList List of parameters in `pdf` that are yields. These must be RooRealVar or RooLinearVar, since RooStats will need to modify their values.
///\param[in] projDeps Don't normalise over these parameters when calculating the sWeights. Will be passed on to AddSWeight().
///\param[in] useWeights Include weights of the input data in calculation of s weights.
///\param[in] cloneData Make a clone of the incoming data before adding weights.
///\param[in] newName New name for the data.
///\param[in] argX Additional arguments for the fitting step in AddSWeight().
SPlot::SPlot(const char* name, const char* title, RooDataSet& data, RooAbsPdf* pdf,
        const RooArgList &yieldsList, const RooArgSet &projDeps,
        bool useWeights, bool cloneData, const char* newName,
        const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8):
  TNamed(name, title)
{
  if(cloneData == 1) {
    fSData = (RooDataSet*) data.Clone(newName);
    SetBit(kOwnData);
  }
  else
    fSData = (RooDataSet*) &data;

  // Add check that yieldsList contains all RooRealVar / RooAbsRealLValue
  for (const auto arg : yieldsList) {
    if (!dynamic_cast<const RooAbsRealLValue*>(arg)) {
      coutE(InputArguments) << "SPlot::SPlot(" << GetName() << ") input argument "
             << arg->GetName() << " is not of type RooRealVar (or RooLinearVar)."
             << "\nRooStats must be able to set it to 0 and to 1 to probe the PDF." << endl ;
      throw std::invalid_argument(Form("SPlot::SPlot(%s) input argument %s is not of type RooRealVar/RooLinearVar",GetName(),arg->GetName())) ;
    }
  }

  //Construct a new SPlot class,
  //calculate sWeights, and include them
  //in the RooDataSet of this class.

  this->AddSWeight(pdf, yieldsList, projDeps, useWeights, arg5, arg6, arg7, arg8);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dataset (if not passed in constructor).
RooDataSet* SPlot::SetSData(RooDataSet* data)
{
  if(data)    {
    fSData = (RooDataSet*) data;
    return fSData;
  }  else
    return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve s-weighted data.
/// It does **not** automatically call AddSWeight(). This needs to be done manually.
RooDataSet* SPlot::GetSDataSet() const
{
  return fSData;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve an s weight.
/// \param[in] numEvent Event number to retrieve s weight for.
/// \param[in] sVariable The yield parameter to retrieve the s weight for.
Double_t SPlot::GetSWeight(Int_t numEvent, const char* sVariable) const
{
  if(numEvent > fSData->numEntries() )
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  if(numEvent < 0)
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  Double_t totalYield = 0;

  std::string varname(sVariable);
  varname += "_sw";


  if(fSWeightVars.find(sVariable) )
    {
      RooArgSet Row(*fSData->get(numEvent));
      totalYield += Row.getRealValue(sVariable);

      return totalYield;
    }

  if( fSWeightVars.find(varname.c_str())  )
    {

      RooArgSet Row(*fSData->get(numEvent));
      totalYield += Row.getRealValue(varname.c_str() );

      return totalYield;
    }

  else
    coutE(InputArguments) << "InputVariable not in list of sWeighted variables" << endl;

  return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Sum the SWeights for a particular event.
/// This sum should equal the total weight of that event.
/// This method is intended to be used as a check.

Double_t SPlot::GetSumOfEventSWeight(Int_t numEvent) const
{
  if(numEvent > fSData->numEntries() )
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  if(numEvent < 0)
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  Int_t numSWeightVars = this->GetNumSWeightVars();

  Double_t eventSWeight = 0;

  RooArgSet Row(*fSData->get(numEvent));

  for (Int_t i = 0; i < numSWeightVars; i++)
    eventSWeight += Row.getRealValue(fSWeightVars.at(i)->GetName() );

  return  eventSWeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Sum the SWeights for a particular species over all events.
/// This should equal the total (weighted) yield of that species.
/// This method is intended as a check.

Double_t SPlot::GetYieldFromSWeight(const char* sVariable) const
{

  Double_t totalYield = 0;

  std::string varname(sVariable);
  varname += "_sw";


  if(fSWeightVars.find(sVariable) )
    {
      for(Int_t i=0; i < fSData->numEntries(); i++)
   {
     RooArgSet Row(*fSData->get(i));
     totalYield += Row.getRealValue(sVariable);
   }

      return totalYield;
    }

  if( fSWeightVars.find(varname.c_str())  )
    {
      for(Int_t i=0; i < fSData->numEntries(); i++)
   {
     RooArgSet Row(*fSData->get(i));
     totalYield += Row.getRealValue(varname.c_str() );
   }

      return totalYield;
    }

  else
    coutE(InputArguments) << "InputVariable not in list of sWeighted variables" << endl;

  return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a RooArgList containing all paramters that have s weights.

RooArgList SPlot::GetSWeightVars() const
{

  RooArgList Args = fSWeightVars;

  return  Args;

}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of SWeights
/// In other words, return the number of
/// species that we are trying to extract.

Int_t SPlot::GetNumSWeightVars() const
{
  RooArgList Args = fSWeightVars;

  return Args.getSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Method which adds the sWeights to the dataset.
///
/// The SPlot will contain two new variables for each yield parameter:
/// - `L_<varname>` is the the likelihood for each event, *i.e.*, the pdf evaluated for the a given value of the variable "varname".
/// - `<varname>_sw` is the value of the sWeight for the variable "varname" for each event.
///
/// Find Parameters in the PDF to be considered fixed when calculating the SWeights
/// and be sure to NOT include the yields in that list.
///
/// After fixing non-yield parameters, this function will start a fit by calling
/// ```
/// pdf->fitTo(*fSData, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1)).
/// ```
/// One can pass additional arguments to `fitTo`, such as `RooFit::Range("fitrange")`, as `arg5`, `arg6`, `arg7`, `arg8`.
///
/// \note A `RooFit::Range` may be necessary to get expected results if you initially fit in a range
/// and/or called `pdf->fixCoefRange("fitrange")` on `pdf`.
/// Pass `arg5`, `arg6`, `arg7`, `arg8` AT YOUR OWN RISK.
///
/// \param[in] pdf PDF to fit to data to compute s weights.
/// \param[in] yieldsTmp Yields to use to compute s weights.
/// \param[in] projDeps These will not be normalized over when calculating the sWeights,
/// and will be considered parameters, not observables.
/// \param[in] includeWeights Include weights of the input data in calculation of s weights.
/// \param[in] argX Optional additional arguments for the fitting step.
void SPlot::AddSWeight( RooAbsPdf* pdf, const RooArgList &yieldsTmp,
         const RooArgSet &projDeps, bool includeWeights,
         const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8)
{

  // Find Parameters in the PDF to be considered fixed when calculating the SWeights
  // and be sure to NOT include the yields in that list
  RooArgList* constParameters = (RooArgList*)pdf->getParameters(fSData) ;
  for (unsigned int i=0; i < constParameters->size(); ++i) {
    // Need a counting loop since collection is being modified
    auto& par = (*constParameters)[i];
    if (std::any_of(yieldsTmp.begin(), yieldsTmp.end(), [&](const RooAbsArg* yield){ return yield->dependsOn(par); })) {
      constParameters->remove(par, kTRUE, kTRUE);
      --i;
    }
  }


  // Set these parameters constant and store them so they can later
  // be set to not constant
  std::vector<RooAbsRealLValue*> constVarHolder;

  for(Int_t i = 0; i < constParameters->getSize(); i++)
  {
    RooAbsRealLValue* varTemp = static_cast<RooAbsRealLValue*>( constParameters->at(i) );
    if(varTemp &&  varTemp->isConstant() == 0 )
    {
      varTemp->setConstant();
      constVarHolder.push_back(varTemp);
    }
  }

  // Fit yields to the data with all other variables held constant
  // This is necessary because SPlot assumes the yields minimise -Log(likelihood)
  pdf->fitTo(*fSData, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), arg5, arg6, arg7, arg8);

  // Hold the value of the fitted yields
  std::vector<double> yieldsHolder;

  for(Int_t i = 0; i < yieldsTmp.getSize(); i++)
    yieldsHolder.push_back(static_cast<RooAbsReal*>(yieldsTmp.at(i))->getVal());

  const Int_t nspec = yieldsTmp.getSize();
  RooArgList yields = *(RooArgList*)yieldsTmp.snapshot(kFALSE);

  if (RooMsgService::instance().isActive(this, RooFit::InputArguments, RooFit::DEBUG)) {
    coutI(InputArguments) << "Printing Yields" << endl;
    yields.Print();
  }

  // The list of variables to normalize over when calculating PDF values.

  RooArgSet vars(*fSData->get() );
  vars.remove(projDeps, kTRUE, kTRUE);

  // Attach data set

  // const_cast<RooAbsPdf*>(pdf)->attachDataSet(*fSData);

  pdf->attachDataSet(*fSData);

  // first calculate the pdf values for all species and all events
  std::vector<RooAbsRealLValue*> yieldvars ;
  RooArgSet pdfServers;
  pdf->treeNodeServerList(&pdfServers);

  std::vector<Double_t> yieldvalues ;
  for (Int_t k = 0; k < nspec; ++k) {
    auto thisyield = static_cast<const RooAbsReal*>(yields.at(k)) ;
    auto yieldinpdf = static_cast<RooAbsRealLValue*>(pdfServers.find(thisyield->GetName()));
    assert(pdf->dependsOn(*yieldinpdf));

    if (yieldinpdf) {
      coutI(InputArguments)<< "yield in pdf: " << yieldinpdf->GetName() << " " << thisyield->getVal() << endl;

      yieldvars.push_back(yieldinpdf) ;
      yieldvalues.push_back(thisyield->getVal()) ;
    }
  }

  Int_t numevents = fSData->numEntries() ;




  // set all yield to zero
  for(Int_t m=0; m<nspec; ++m) {
    auto theVar = static_cast<RooAbsRealLValue*>(yieldvars[m]);
    theVar->setVal(0) ;

    //Check that range of yields is at least (0,1), and fix otherwise
    if (theVar->getMin() > 0) {
      coutE(InputArguments)  << "Yield variables need to have a range that includes at least [0, 1]. Minimum for "
          << theVar->GetName() << " is " << theVar->getMin() << std::endl;
      if (RooRealVar* realVar = dynamic_cast<RooRealVar*>(theVar)) {
        coutE(InputArguments)  << "Setting min range to 0" << std::endl;
        realVar->setMin(0);
      } else {
        throw std::invalid_argument(std::string("Yield variable ") + theVar->GetName() + " must have a range that includes 0.");
      }
    }

    if (theVar->getMax() < 1) {
      coutW(InputArguments)  << "Yield variables need to have a range that includes at least [0, 1]. Maximum for "
          << theVar->GetName() << " is " << theVar->getMax() << std::endl;
      if (RooRealVar* realVar = dynamic_cast<RooRealVar*>(theVar)) {
        coutE(InputArguments)  << "Setting max range to 1" << std::endl;
        realVar->setMax(1);
      } else {
        throw std::invalid_argument(std::string("Yield variable ") + theVar->GetName() + " must have a range that includes 1.");
      }
    }
  }


  // For every event and for every species,
  // calculate the value of the component pdf for that specie
  // by setting the yield of that specie to 1
  // and all others to 0.  Evaluate the pdf for each event
  // and store the values.

  RooArgSet * pdfvars = pdf->getVariables();
  std::vector<std::vector<Double_t> > pdfvalues(numevents,std::vector<Double_t>(nspec,0)) ;

  for (Int_t ievt = 0; ievt <numevents; ievt++)
  {
    //WVE: FIX THIS PART, EVALUATION PROGRESS!!

    RooStats::SetParameters(fSData->get(ievt), pdfvars);

    for(Int_t k = 0; k < nspec; ++k) {
      auto theVar = static_cast<RooAbsRealLValue*>(yieldvars[k]);

      // set this yield to 1
      theVar->setVal( 1 ) ;
      // evaluate the pdf
      Double_t f_k = pdf->getVal(&vars) ;
      pdfvalues[ievt][k] = f_k ;
      if( !(f_k>1 || f_k<1) )
        coutW(InputArguments) << "Strange pdf value: " << ievt << " " << k << " " << f_k << std::endl ;
      theVar->setVal( 0 ) ;
    }
  }
  delete pdfvars;

  // check that the likelihood normalization is fine
  std::vector<Double_t> norm(nspec,0) ;
  for (Int_t ievt = 0; ievt <numevents ; ievt++)
    {
      Double_t dnorm(0) ;
      for(Int_t k=0; k<nspec; ++k) dnorm += yieldvalues[k] * pdfvalues[ievt][k] ;
      for(Int_t j=0; j<nspec; ++j) norm[j] += pdfvalues[ievt][j]/dnorm ;
    }

  coutI(Contents) << "likelihood norms: "  ;

  for(Int_t k=0; k<nspec; ++k)  coutI(Contents) << norm[k] << " " ;
  coutI(Contents) << std::endl ;

  // Make a TMatrixD to hold the covariance matrix.
  TMatrixD covInv(nspec, nspec);
  for (Int_t i = 0; i < nspec; i++) for (Int_t j = 0; j < nspec; j++) covInv(i,j) = 0;

  coutI(Contents) << "Calculating covariance matrix";


  // Calculate the inverse covariance matrix, using weights
  for (Int_t ievt = 0; ievt < numevents; ++ievt)
    {

      fSData->get(ievt) ;

      // Calculate contribution to the inverse of the covariance
      // matrix. See BAD 509 V2 eqn. 15

      // Sum for the denominator
      Double_t dsum(0);
      for(Int_t k = 0; k < nspec; ++k)
   dsum += pdfvalues[ievt][k] * yieldvalues[k] ;

      for(Int_t n=0; n<nspec; ++n)
   for(Int_t j=0; j<nspec; ++j)
     {
       if(includeWeights)
         covInv(n,j) +=  fSData->weight()*pdfvalues[ievt][n]*pdfvalues[ievt][j]/(dsum*dsum) ;
       else
         covInv(n,j) +=                   pdfvalues[ievt][n]*pdfvalues[ievt][j]/(dsum*dsum) ;
     }

      //ADDED WEIGHT ABOVE

    }

  // Covariance inverse should now be computed!

  // Invert to get the covariance matrix
  if (covInv.Determinant() <=0)
    {
      coutE(Eval) << "SPlot Error: covariance matrix is singular; I can't invert it!" << std::endl;
      covInv.Print();
      return;
    }

  TMatrixD covMatrix(TMatrixD::kInverted,covInv);

  //check cov normalization
  if (RooMsgService::instance().isActive(this, RooFit::Eval, RooFit::DEBUG)) {
    coutI(Eval) << "Checking Likelihood normalization:  " << std::endl;
    coutI(Eval) << "Yield of specie  Sum of Row in Matrix   Norm" << std::endl;
    for(Int_t k=0; k<nspec; ++k)
    {
      Double_t covnorm(0) ;
      for(Int_t m=0; m<nspec; ++m) covnorm += covInv[k][m]*yieldvalues[m] ;
      Double_t sumrow(0) ;
      for(Int_t m = 0; m < nspec; ++m) sumrow += covMatrix[k][m] ;
      coutI(Eval)  << yieldvalues[k] << " " << sumrow << " " << covnorm << endl ;
    }
  }

  // calculate for each event the sWeight (BAD 509 V2 eq. 21)
  coutI(Eval) << "Calculating sWeight" << std::endl;
  std::vector<RooRealVar*> sweightvec ;
  std::vector<RooRealVar*> pdfvec ;
  RooArgSet sweightset ;

  // Create and label the variables
  // used to store the SWeights

  fSWeightVars.Clear();

  for(Int_t k=0; k<nspec; ++k)
    {
       std::string wname = std::string(yieldvars[k]->GetName()) + "_sw";
       RooRealVar* var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
       sweightvec.push_back( var) ;
       sweightset.add(*var) ;
       fSWeightVars.add(*var);

       wname = "L_" + std::string(yieldvars[k]->GetName());
       var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
       pdfvec.push_back( var) ;
       sweightset.add(*var) ;
    }

  // Create and fill a RooDataSet
  // with the SWeights

  RooDataSet* sWeightData = new RooDataSet("dataset", "dataset with sWeights", sweightset);

  for(Int_t ievt = 0; ievt < numevents; ++ievt)
    {

      fSData->get(ievt) ;

      // sum for denominator
      Double_t dsum(0);
      for(Int_t k = 0; k < nspec; ++k)   dsum +=  pdfvalues[ievt][k] * yieldvalues[k] ;
      // covariance weighted pdf for each specief
      for(Int_t n=0; n<nspec; ++n)
   {
     Double_t nsum(0) ;
     for(Int_t j=0; j<nspec; ++j) nsum += covMatrix(n,j) * pdfvalues[ievt][j] ;


     //Add the sWeights here!!
     //Include weights,
     //ie events weights are absorbed into sWeight


     if(includeWeights) sweightvec[n]->setVal(fSData->weight() * nsum/dsum) ;
     else  sweightvec[n]->setVal( nsum/dsum) ;

     pdfvec[n]->setVal( pdfvalues[ievt][n] ) ;

     if( !(fabs(nsum/dsum)>=0 ) )
       {
         coutE(Contents) << "error: " << nsum/dsum << endl ;
         return;
       }
   }

      sWeightData->add(sweightset) ;
    }


  // Add the SWeights to the original data set

  fSData->merge(sWeightData);

  delete sWeightData;

  //Restore yield values

  for(Int_t i = 0; i < yieldsTmp.getSize(); i++)
    static_cast<RooAbsRealLValue*>(yieldsTmp.at(i))->setVal(yieldsHolder.at(i));

  //Make any variables that were forced to constant no longer constant

  for(Int_t i=0; i < (Int_t) constVarHolder.size(); i++)
    constVarHolder.at(i)->setConstant(kFALSE);

  return;

}
