// @(#)root/roostats:$Id$
// Author: Sven Kreiss    June 2010
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::NuisanceParametersSampler
    \ingroup Roostats

Helper class for ToyMCSampler. Handles all of the nuisance parameter related
functions. Once instantiated, it gives a new nuisance parameter point
at each call to nextPoint(...).
*/

/** \class RooStats::ToyMCSampler
    \ingroup Roostats

ToyMCSampler is an implementation of the TestStatSampler interface.
It generates Toy Monte Carlo for a given parameter point and evaluates a
TestStatistic.

For parallel runs, ToyMCSampler can be given an instance of ProofConfig
and then run in parallel using proof or proof-lite. Internally, it uses
ToyMCStudy with the RooStudyManager.
*/

#include "RooStats/ToyMCSampler.h"

#include "RooMsgService.h"

#include "RooDataHist.h"

#include "RooRealVar.h"

#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRandom.h"

#include "RooStudyManager.h"
#include "RooStats/ToyMCStudy.h"
#include "RooStats/DetailedOutputAggregator.h"
#include "RooStats/RooStatsUtils.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"

#include "TMath.h"


using namespace RooFit;
using namespace std;


ClassImp(RooStats::ToyMCSampler);

namespace RooStats {

////////////////////////////////////////////////////////////////////////////////
/// Assigns new nuisance parameter point to members of nuisPoint.
/// nuisPoint can be more objects than just the nuisance
/// parameters.

void NuisanceParametersSampler::NextPoint(RooArgSet& nuisPoint, Double_t& weight) {

   // check whether to get new set of nuisanceParPoints
   if (fIndex >= fNToys) {
      Refresh();
      fIndex = 0;
   }

   // get value
   nuisPoint =  *fPoints->get(fIndex++);
   weight = fPoints->weight();

   // check whether result will have any influence
   if(fPoints->weight() == 0.0) {
      oocoutI((TObject*)NULL,Generation) << "Weight 0 encountered. Skipping." << endl;
      NextPoint(nuisPoint, weight);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the initial set of nuisance parameter points. It also refills the
/// set with new parameter points if called repeatedly. This helps with
/// adaptive sampling as the required number of nuisance parameter points
/// might increase during the run.

void NuisanceParametersSampler::Refresh() {

   if (!fPrior || !fParams) return;

   if (fPoints) delete fPoints;

   if (fExpected) {
      // UNDER CONSTRUCTION
      oocoutI((TObject*)NULL,InputArguments) << "Using expected nuisance parameters." << endl;

      int nBins = fNToys;

      // From FeldmanCousins.cxx:
      // set nbins for the POI
      TIter it2 = fParams->createIterator();
      RooRealVar *myarg2;
      while ((myarg2 = dynamic_cast<RooRealVar*>(it2.Next()))) {
        myarg2->setBins(nBins);
      }


      fPoints = fPrior->generate(
         *fParams,
         AllBinned(),
         ExpectedData(),
         NumEvents(1) // for Asimov set, this is only a scale factor
      );
      if(fPoints->numEntries() != fNToys) {
         fNToys = fPoints->numEntries();
         oocoutI((TObject*)NULL,InputArguments) <<
            "Adjusted number of toys to number of bins of nuisance parameters: " << fNToys << endl;
      }

/*
      // check
      TCanvas *c1 = new TCanvas;
      RooPlot *p = dynamic_cast<RooRealVar*>(fParams->first())->frame();
      fPoints->plotOn(p);
      p->Draw();
      for(int x=0; x < fPoints->numEntries(); x++) {
         fPoints->get(x)->Print("v");
         cout << fPoints->weight() << endl;
      }
*/

   }else{
      oocoutI((TObject*)NULL,InputArguments) << "Using randomized nuisance parameters." << endl;

      fPoints = fPrior->generate(*fParams, fNToys);
   }
}

Bool_t ToyMCSampler::fgAlwaysUseMultiGen = kFALSE ;

////////////////////////////////////////////////////////////////////////////////

void ToyMCSampler::SetAlwaysUseMultiGen(Bool_t flag) { fgAlwaysUseMultiGen = flag ; }

////////////////////////////////////////////////////////////////////////////////
/// Proof constructor. Do not use.

ToyMCSampler::ToyMCSampler() : fSamplingDistName("SD"), fNToys(1)
{

   fPdf = NULL;
   fParametersForTestStat = NULL;
   fPriorNuisance = NULL;
   fNuisancePars = NULL;
   fObservables = NULL;
   fGlobalObservables = NULL;

   fSize = 0.05;
   fNEvents = 0;
   fGenerateBinned = kFALSE;
   fGenerateBinnedTag = "";
   fGenerateAutoBinned = kTRUE;
   fExpectedNuisancePar = kFALSE;

   fToysInTails = 0.0;
   fMaxToys = RooNumber::infinity();
   fAdaptiveLowLimit = -RooNumber::infinity();
   fAdaptiveHighLimit = RooNumber::infinity();

   fProtoData = NULL;

   fProofConfig = NULL;
   fNuisanceParametersSampler = NULL;

   _allVars = NULL ;
   _gs1 = NULL ;
   _gs2 = NULL ;
   _gs3 = NULL ;
   _gs4 = NULL ;

   //suppress messages for num integration of Roofit
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

   fUseMultiGen = kFALSE ;
}

////////////////////////////////////////////////////////////////////////////////

ToyMCSampler::ToyMCSampler(TestStatistic &ts, Int_t ntoys) : fSamplingDistName(ts.GetVarName().Data()), fNToys(ntoys)
{
   fPdf = NULL;
   fParametersForTestStat = NULL;
   fPriorNuisance = NULL;
   fNuisancePars = NULL;
   fObservables = NULL;
   fGlobalObservables = NULL;

   fSize = 0.05;
   fNEvents = 0;
   fGenerateBinned = kFALSE;
   fGenerateBinnedTag = "";
   fGenerateAutoBinned = kTRUE;
   fExpectedNuisancePar = kFALSE;

   fToysInTails = 0.0;
   fMaxToys = RooNumber::infinity();
   fAdaptiveLowLimit = -RooNumber::infinity();
   fAdaptiveHighLimit = RooNumber::infinity();

   fProtoData = NULL;

   fProofConfig = NULL;
   fNuisanceParametersSampler = NULL;

   _allVars = NULL ;
   _gs1 = NULL ;
   _gs2 = NULL ;
   _gs3 = NULL ;
   _gs4 = NULL ;

   //suppress messages for num integration of Roofit
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

   fUseMultiGen = kFALSE ;

   AddTestStatistic(&ts);
}

////////////////////////////////////////////////////////////////////////////////

ToyMCSampler::~ToyMCSampler() {
   if(fNuisanceParametersSampler) delete fNuisanceParametersSampler;

   ClearCache();
}

////////////////////////////////////////////////////////////////////////////////
/// only checks, no guessing/determination (do this in calculators,
/// e.g. using ModelConfig::GuessObsAndNuisance(...))

Bool_t ToyMCSampler::CheckConfig(void) {
   bool goodConfig = true;

   if(fTestStatistics.size() == 0 || fTestStatistics[0] == NULL) { ooccoutE((TObject*)NULL,InputArguments) << "Test statistic not set." << endl; goodConfig = false; }
   if(!fObservables) { ooccoutE((TObject*)NULL,InputArguments) << "Observables not set." << endl; goodConfig = false; }
   if(!fParametersForTestStat) { ooccoutE((TObject*)NULL,InputArguments) << "Parameter values used to evaluate the test statistic are not set." << endl; goodConfig = false; }
   if(!fPdf) { ooccoutE((TObject*)NULL,InputArguments) << "Pdf not set." << endl; goodConfig = false; }


   //ooccoutI((TObject*)NULL,InputArguments) << "ToyMCSampler configuration:" << endl;
   //ooccoutI((TObject*)NULL,InputArguments) << "Pdf from SetPdf: " << fPdf << endl;
   // for( unsigned int i=0; i < fTestStatistics.size(); i++ ) {
   //   ooccoutI((TObject*)NULL,InputArguments) << "test statistics["<<i<<"]: " << fTestStatistics[i] << endl;
   // }
   //ooccoutI((TObject*)NULL,InputArguments) << endl;

   return goodConfig;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate all test statistics, returning result and any detailed output.
/// PDF parameter values are saved in case they are modified by
/// TestStatistic::Evaluate (eg. SimpleLikelihoodRatioTestStat).

RooArgList* ToyMCSampler::EvaluateAllTestStatistics(RooAbsData& data, const RooArgSet& poi) {
   DetailedOutputAggregator detOutAgg;
   const RooArgList* allTS = EvaluateAllTestStatistics(data, poi, detOutAgg);
   if (!allTS) return 0;
   // no need to delete allTS, it is deleted in destructor of detOutAgg
   return  dynamic_cast<RooArgList*>(allTS->snapshot());
}

////////////////////////////////////////////////////////////////////////////////

const RooArgList* ToyMCSampler::EvaluateAllTestStatistics(RooAbsData& data, const RooArgSet& poi, DetailedOutputAggregator& detOutAgg) {
   RooArgSet *allVars = fPdf ? fPdf->getVariables() : nullptr;
   RooArgSet *saveAll = allVars ? allVars->snapshot() : nullptr;
   for( unsigned int i = 0; i < fTestStatistics.size(); i++ ) {
      if( fTestStatistics[i] == NULL ) continue;
      TString name( TString::Format("%s_TS%u", fSamplingDistName.c_str(), i) );
      std::unique_ptr<RooArgSet> parForTS(poi.snapshot());
      RooRealVar ts( name, fTestStatistics[i]->GetVarName(), fTestStatistics[i]->Evaluate( data, *parForTS ) );
      RooArgList tset(ts);
      detOutAgg.AppendArgSet(&tset);
      if (const RooArgSet* detOut = fTestStatistics[i]->GetDetailedOutput()) {
        name.Append("_");
        detOutAgg.AppendArgSet(detOut, name);
      }
      if (saveAll) *allVars = *saveAll;  // restore values, perhaps modified by fTestStatistics[i]->Evaluate()
   }
   delete saveAll;
   delete allVars;
   return detOutAgg.GetAsArgList();
}

////////////////////////////////////////////////////////////////////////////////

SamplingDistribution* ToyMCSampler::GetSamplingDistribution(RooArgSet& paramPointIn) {
   if(fTestStatistics.size() > 1) {
      oocoutW((TObject*)NULL, InputArguments) << "Multiple test statistics defined, but only one distribution will be returned." << endl;
      for( unsigned int i=0; i < fTestStatistics.size(); i++ ) {
         oocoutW((TObject*)NULL, InputArguments) << " \t test statistic: " << fTestStatistics[i] << endl;
      }
   }

   RooDataSet* r = GetSamplingDistributions(paramPointIn);
   if(r == NULL || r->numEntries() == 0) {
      oocoutW((TObject*)NULL, Generation) << "no sampling distribution generated" << endl;
      return NULL;
   }

   SamplingDistribution* samp = new SamplingDistribution( r->GetName(), r->GetTitle(), *r );
   delete r;
   return samp;
}

////////////////////////////////////////////////////////////////////////////////
/// Use for serial and parallel runs.

RooDataSet* ToyMCSampler::GetSamplingDistributions(RooArgSet& paramPointIn)
{

   // ======= S I N G L E   R U N ? =======
   if(!fProofConfig)
      return GetSamplingDistributionsSingleWorker(paramPointIn);

   // ======= P A R A L L E L   R U N =======
   if (!CheckConfig()){
      oocoutE((TObject*)NULL, InputArguments)
         << "Bad COnfiguration in ToyMCSampler "
         << endl;
      return nullptr;
   }

   // turn adaptive sampling off if given
   if(fToysInTails) {
      fToysInTails = 0;
      oocoutW((TObject*)NULL, InputArguments)
         << "Adaptive sampling in ToyMCSampler is not supported for parallel runs."
         << endl;
   }

   // adjust number of toys on the slaves to keep the total number of toys constant
   Int_t totToys = fNToys;
   fNToys = (int)ceil((double)fNToys / (double)fProofConfig->GetNExperiments()); // round up

   // create the study instance for parallel processing
   ToyMCStudy* toymcstudy = new ToyMCStudy ;
   toymcstudy->SetToyMCSampler(*this);
   toymcstudy->SetParamPoint(paramPointIn);
   toymcstudy->SetRandomSeed(RooRandom::randomGenerator()->Integer(TMath::Limits<unsigned int>::Max() ) );

   // temporary workspace for proof to avoid messing with TRef
   RooWorkspace w(fProofConfig->GetWorkspace());
   RooStudyManager studymanager(w, *toymcstudy);
   studymanager.runProof(fProofConfig->GetNExperiments(), fProofConfig->GetHost(), fProofConfig->GetShowGui());

   RooDataSet* output = toymcstudy->merge();

   // reset the number of toys
   fNToys = totToys;

   delete toymcstudy;
   return output;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the main function for serial runs. It is called automatically
/// from inside GetSamplingDistribution when no ProofConfig is given.
/// You should not call this function yourself. This function should
/// be used by ToyMCStudy on the workers (ie. when you explicitly want
/// a serial run although ProofConfig is present).
///

RooDataSet* ToyMCSampler::GetSamplingDistributionsSingleWorker(RooArgSet& paramPointIn)
{
  // Make sure the cache is clear. It is important to clear it here, because
  // the cache might be invalid even when just the firstPOI was changed, for which
  // no accessor has to be called. (Fixes a bug when ToyMCSampler is
  // used with the Neyman Construction)
   ClearCache();

   if (!CheckConfig()){
      oocoutE((TObject*)NULL, InputArguments)
         << "Bad COnfiguration in ToyMCSampler "
         << endl;
      return nullptr;
   }

   // important to cache the paramPoint b/c test statistic might
   // modify it from event to event
   RooArgSet *paramPoint = (RooArgSet*) paramPointIn.snapshot();
   RooArgSet *allVars = fPdf->getVariables();
   RooArgSet *saveAll = (RooArgSet*) allVars->snapshot();


   DetailedOutputAggregator detOutAgg;

   // counts the number of toys in the limits set for adaptive sampling
   // (taking weights into account; always on first test statistic)
   Double_t toysInTails = 0.0;

   for (Int_t i = 0; i < fMaxToys; ++i) {
      // need to check at the beginning for case that zero toys are requested
      if (toysInTails >= fToysInTails  &&  i+1 > fNToys) break;

      // status update
      if ( i% 500 == 0 && i>0 ) {
         oocoutP((TObject*)0,Generation) << "generated toys: " << i << " / " << fNToys;
         if (fToysInTails) ooccoutP((TObject*)0,Generation) << " (tails: " << toysInTails << " / " << fToysInTails << ")" << std::endl;
         else ooccoutP((TObject*)0,Generation) << endl;
      }

      // TODO: change this treatment to keep track of all values so that the threshold
      // for adaptive sampling is counted for all distributions and not just the
      // first one.
      Double_t valueFirst = -999.0, weight = 1.0;

      // set variables to requested parameter point
      *allVars = *saveAll; // important for example for SimpleLikelihoodRatioTestStat

      RooAbsData* toydata = GenerateToyData(*paramPoint, weight);
      if (i == 0 && !fPdf->canBeExtended() && dynamic_cast<RooSimultaneous*>(fPdf)) {
        const RooArgSet* toySet = toydata->get();
        if (std::none_of(toySet->begin(), toySet->end(), [](const RooAbsArg* arg){
          return dynamic_cast<const RooAbsCategory*>(arg) != nullptr;
        }))
          oocoutE((TObject*)nullptr, Generation) << "ToyMCSampler: Generated toy data didn't contain a category variable, although"
            " a simultaneous PDF is in use. To generate events for a simultaneous PDF, all components need to be"
            " extended. Otherwise, the number of events to generate per component cannot be determined." << std::endl;
      }

      *allVars = *fParametersForTestStat;

      const RooArgList* allTS = EvaluateAllTestStatistics(*toydata, *fParametersForTestStat, detOutAgg);
      if (allTS->getSize() > Int_t(fTestStatistics.size()))
        detOutAgg.AppendArgSet( fGlobalObservables, "globObs_" );
      if (RooRealVar* firstTS = dynamic_cast<RooRealVar*>(allTS->first()))
         valueFirst = firstTS->getVal();

      delete toydata;

      // check for nan
      if(valueFirst != valueFirst) {
         oocoutW((TObject*)NULL, Generation) << "skip: " << valueFirst << ", " << weight << endl;
         continue;
      }

      detOutAgg.CommitSet(weight);

      // adaptive sampling checks
      if (valueFirst <= fAdaptiveLowLimit  ||  valueFirst >= fAdaptiveHighLimit) {
         if(weight >= 0.) toysInTails += weight;
         else toysInTails += 1.;
      }
   }

   // clean up
   *allVars = *saveAll;
   delete saveAll;
   delete allVars;
   delete paramPoint;

   return detOutAgg.GetAsDataSet(fSamplingDistName, fSamplingDistName);
}

////////////////////////////////////////////////////////////////////////////////

void ToyMCSampler::GenerateGlobalObservables(RooAbsPdf& pdf) const {


   if(!fGlobalObservables  ||  fGlobalObservables->getSize()==0) {
      ooccoutE((TObject*)NULL,InputArguments) << "Global Observables not set." << endl;
      return;
   }


   if (fUseMultiGen || fgAlwaysUseMultiGen) {

      // generate one set of global observables and assign it
      // has problem for sim pdfs
      RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>( &pdf );
      if (!simPdf) {
         RooDataSet *one = pdf.generate(*fGlobalObservables, 1);

         const RooArgSet *values = one->get(0);
         if (!_allVars) {
            _allVars = pdf.getVariables();
         }
         *_allVars = *values;
         delete one;

      } else {

         if (_pdfList.size() == 0) {
            RooCategory& channelCat = (RooCategory&)simPdf->indexCat();
            int nCat = channelCat.numTypes();
            for (int i=0; i < nCat; ++i){
               channelCat.setIndex(i);
               RooAbsPdf* pdftmp = simPdf->getPdf(channelCat.getCurrentLabel());
               assert(pdftmp);
               RooArgSet* globtmp = pdftmp->getObservables(*fGlobalObservables);
               RooAbsPdf::GenSpec* gs = pdftmp->prepareMultiGen(*globtmp, NumEvents(1));
               _pdfList.push_back(pdftmp);
               _obsList.push_back(globtmp);
               _gsList.push_back(gs);
            }
         }

         list<RooArgSet*>::iterator oiter = _obsList.begin();
         list<RooAbsPdf::GenSpec*>::iterator giter = _gsList.begin();
         for (list<RooAbsPdf*>::iterator iter = _pdfList.begin(); iter != _pdfList.end(); ++iter, ++giter, ++oiter) {
            //RooDataSet* tmp = (*iter)->generate(**oiter,1) ;
            RooDataSet* tmp = (*iter)->generate(**giter);
            **oiter = *tmp->get(0);
            delete tmp;
         }
      }


   } else {

      // not using multigen for global observables
      RooDataSet* one = pdf.generateSimGlobal( *fGlobalObservables, 1 );
      const RooArgSet *values = one->get(0);
      RooArgSet* allVars = pdf.getVariables();
      *allVars = *values;
      delete allVars;
      delete one;

   }
}

////////////////////////////////////////////////////////////////////////////////
/// This method generates a toy data set for the given parameter point taking
/// global observables into account.
/// The values of the generated global observables remain in the pdf's variables.
/// They have to have those values for the subsequent evaluation of the
/// test statistics.

RooAbsData* ToyMCSampler::GenerateToyData(RooArgSet& paramPoint, double& weight, RooAbsPdf& pdf) const {

   if(!fObservables) {
      ooccoutE((TObject*)NULL,InputArguments) << "Observables not set." << endl;
      return NULL;
   }

   // assign input paramPoint
   RooArgSet* allVars = fPdf->getVariables();
   *allVars = paramPoint;


   // create nuisance parameter points
   if(!fNuisanceParametersSampler && fPriorNuisance && fNuisancePars) {
      fNuisanceParametersSampler = new NuisanceParametersSampler(fPriorNuisance, fNuisancePars, fNToys, fExpectedNuisancePar);
      if ((fUseMultiGen || fgAlwaysUseMultiGen) &&  fNuisanceParametersSampler )
         oocoutI((TObject*)NULL,InputArguments) << "Cannot use multigen when nuisance parameters vary for every toy" << endl;
   }

   // generate global observables
   RooArgSet observables(*fObservables);
   if(fGlobalObservables  &&  fGlobalObservables->getSize()) {
      observables.remove(*fGlobalObservables);
      GenerateGlobalObservables(pdf);
   }

   // save values to restore later.
   // but this must remain after(!) generating global observables
   const RooArgSet* saveVars = (const RooArgSet*)allVars->snapshot();

   if(fNuisanceParametersSampler) { // use nuisance parameters?
      // Construct a set of nuisance parameters that has the parameters
      // in the input paramPoint removed. Therefore, no parameter in
      // paramPoint is randomized.
      // Therefore when a parameter is given (should be held fixed),
      // but is also in the list of nuisance parameters, the parameter
      // will be held fixed. This is useful for debugging to hold single
      // parameters fixed although under "normal" circumstances it is
      // randomized.
      RooArgSet allVarsMinusParamPoint(*allVars);
      allVarsMinusParamPoint.remove(paramPoint, kFALSE, kTRUE); // match by name

      // get nuisance parameter point and weight
      fNuisanceParametersSampler->NextPoint(allVarsMinusParamPoint, weight);


   }else{
      weight = 1.0;
   }

   RooAbsData *data = Generate(pdf, observables);

   // We generated the data with the randomized nuisance parameter (if hybrid)
   // but now we are setting the nuisance parameters back to where they were.
   *allVars = *saveVars;
   delete allVars;
   delete saveVars;

   return data;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the generate function to use in the context of the ToyMCSampler
/// instead of the standard RooAbsPdf::generate(...).
/// It takes into account whether the number of events is given explicitly
/// or whether it should use the expected number of events. It also takes
/// into account the option to generate a binned data set (*i.e.* RooDataHist).

RooAbsData* ToyMCSampler::Generate(RooAbsPdf &pdf, RooArgSet &observables, const RooDataSet* protoData, int forceEvents) const {

  if(fProtoData) {
    protoData = fProtoData;
    forceEvents = protoData->numEntries();
  }

  RooAbsData *data = NULL;
  int events = forceEvents;
  if(events == 0) events = fNEvents;

  // cannot use multigen when the nuisance parameters change for every toy
  bool useMultiGen = (fUseMultiGen || fgAlwaysUseMultiGen) && !fNuisanceParametersSampler;

  if (events == 0) {
    if (pdf.canBeExtended() && pdf.expectedEvents(observables) > 0) {
      if(fGenerateBinned) {
        if(protoData) data = pdf.generate(observables, AllBinned(), Extended(), ProtoData(*protoData, true, true));
        else          data = pdf.generate(observables, AllBinned(), Extended());
      } else {
        if (protoData) {
          if (useMultiGen) {
            if (!_gs2) { _gs2 = pdf.prepareMultiGen(observables, Extended(), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag), ProtoData(*protoData, true, true)) ; }
            data = pdf.generate(*_gs2) ;
          } else {
            data = pdf.generate                    (observables, Extended(), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag), ProtoData(*protoData, true, true));
          }
        } else {
          if (useMultiGen) {
            if (!_gs1) { _gs1 = pdf.prepareMultiGen(observables, Extended(), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag) ) ; }
            data = pdf.generate(*_gs1) ;
          } else {
            data = pdf.generate                    (observables, Extended(), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag) );
          }

        }
      }
    } else {
      oocoutE((TObject*)0,InputArguments)
                << "ToyMCSampler: Error : pdf is not extended and number of events per toy is zero"
                << endl;
    }
  } else {
    if (fGenerateBinned) {
      if(protoData) data = pdf.generate(observables, events, AllBinned(), ProtoData(*protoData, true, true));
      else          data = pdf.generate(observables, events, AllBinned());
    } else {
      if (protoData) {
        if (useMultiGen) {
          if (!_gs3) { _gs3 = pdf.prepareMultiGen(observables, NumEvents(events), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag), ProtoData(*protoData, true, true)); }
          data = pdf.generate(*_gs3) ;
        } else {
          data = pdf.generate                    (observables, NumEvents(events), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag), ProtoData(*protoData, true, true));
        }
      } else {
        if (useMultiGen) {
          if (!_gs4) { _gs4 = pdf.prepareMultiGen(observables, NumEvents(events), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag)); }
          data = pdf.generate(*_gs4) ;
        } else {
          data = pdf.generate                    (observables, NumEvents(events), AutoBinned(fGenerateAutoBinned), GenBinned(fGenerateBinnedTag));
        }
      }
    }
  }

  // in case of number counting print observables
  // if (data->numEntries() == 1) {
  //    std::cout << "generate observables : ";
  //    RooStats::PrintListContent(*data->get(0), std::cout);
  // }

  return data;
}

////////////////////////////////////////////////////////////////////////////////
/// Extended interface to append to sampling distribution more samples

SamplingDistribution* ToyMCSampler::AppendSamplingDistribution(
   RooArgSet& allParameters,
   SamplingDistribution* last,
   Int_t additionalMC)
{
   Int_t tmp = fNToys;
   fNToys = additionalMC;
   SamplingDistribution* newSamples = GetSamplingDistribution(allParameters);
   fNToys = tmp;

   if(last){
     last->Add(newSamples);
     delete newSamples;
     return last;
   }

   return newSamples;
}

////////////////////////////////////////////////////////////////////////////////
/// clear the cache obtained from the pdf used for speeding the toy and global observables generation
/// needs to be called every time the model pdf (fPdf) changes

void ToyMCSampler::ClearCache() {

   if (_gs1) delete _gs1;
   _gs1 = 0;
   if (_gs2) delete _gs2;
   _gs2 = 0;
   if (_gs3) delete _gs3;
   _gs3 = 0;
   if (_gs4) delete _gs4;
   _gs4 = 0;

   // no need to delete the _pdfList since it is managed by the RooSimultaneous object
   if (_pdfList.size() > 0) {
      std::list<RooArgSet*>::iterator oiter = _obsList.begin();
      for (std::list<RooAbsPdf::GenSpec*>::iterator giter  = _gsList.begin(); giter != _gsList.end(); ++giter, ++oiter) {
         delete *oiter;
         delete *giter;
      }
      _pdfList.clear();
      _obsList.clear();
      _gsList.clear();
   }

   //LM: is this set really needed ??
   if (_allVars) delete _allVars;
   _allVars = 0;

}

} // end namespace RooStats
