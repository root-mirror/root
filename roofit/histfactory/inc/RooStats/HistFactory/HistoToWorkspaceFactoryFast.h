// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HISTOTOWORKSPACEFACTORYFAST
#define ROOSTATS_HISTOTOWORKSPACEFACTORYFAST

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TObject.h>
#include <TH1.h>
#include <TDirectory.h>

#include "RooStats/HistFactory/Systematics.h"
class ParamHistFunc;


namespace RooStats{
  namespace HistFactory{

    // Forward Declarations FTW
    class Measurement;
    class Channel;
    class Sample;

    class HistoToWorkspaceFactoryFast: public TObject {
    
    public:


      HistoToWorkspaceFactoryFast();
      HistoToWorkspaceFactoryFast(  RooStats::HistFactory::Measurement& Meas );
      virtual ~HistoToWorkspaceFactoryFast();

      static void ConfigureWorkspaceForMeasurement( const std::string& ModelName, 
						    RooWorkspace* ws_single, 
						    Measurement& measurement );
    
      RooWorkspace* MakeSingleChannelModel( Measurement& measurement, Channel& channel );
      RooWorkspace*  MakeCombinedModel(std::vector<std::string>, std::vector<RooWorkspace*>);
    
      static RooWorkspace* MakeCombinedModel( Measurement& measurement );
      static void PrintCovarianceMatrix(RooFitResult* result, RooArgSet* params, 
					std::string filename);

      void SetFunctionsToPreprocess(std::vector<std::string> lines) { fPreprocessFunctions=lines; }

    protected:

       void AddConstraintTerms(RooWorkspace* proto, Measurement& measurement, std::string prefix, std::string interpName,
			      std::vector<OverallSys>& systList, 			 
			      std::vector<std::string>& likelihoodTermNames, 
			      std::vector<std::string>& totSystTermNames);

      std::string AddNormFactor(RooWorkspace* proto, std::string& channel, 
				std::string& sigmaEpsilon, Sample& sample, bool doRatio);

      void AddMultiVarGaussConstraint(RooWorkspace* proto, std::string prefix, 
				      int lowBin, int highBin, 
				      std::vector<std::string>& likelihoodTermNames);
    
      void AddPoissonTerms(RooWorkspace* proto, std::string prefix, std::string obsPrefix, 
			   std::string expPrefix, int lowBin, int highBin,
			   std::vector<std::string>& likelihoodTermNames);
    
      static void EditSyst(RooWorkspace* proto, const char* pdfNameChar, 
			   std::map<std::string,double> gammaSyst, 
			   std::map<std::string,double> uniformSyst, 
			   std::map<std::string,double> logNormSyst, 
			   std::map<std::string,double> noSyst);

      void LinInterpWithConstraint(RooWorkspace* proto, const TH1* nominal, std::vector<HistoSys>,
				   std::string prefix, std::string productPrefix, 
				   std::string systTerm, 
				   std::vector<std::string>& likelihoodTermNames);

      RooWorkspace* MakeSingleChannelWorkspace(Measurement& measurement, Channel& channel);

      void MakeTotalExpected(RooWorkspace* proto, std::string totName, 
			     std::vector<std::string>& syst_x_expectedPrefixNames,
			     std::vector<std::string>& normByNames);
    
      RooDataSet* MergeDataSets(RooWorkspace* combined,
				std::vector<RooWorkspace*> wspace_vec, 
				std::vector<std::string> channel_names, 
				std::string dataSetName,
				RooArgList obsList,
				RooCategory* channelCat);

      void ProcessExpectedHisto(const TH1* hist, RooWorkspace* proto, std::string prefix,
				std::string productPrefix, std::string systTerm );

      void SetObsToExpected(RooWorkspace* proto, std::string obsPrefix, std::string expPrefix, 
			    int lowBin, int highBin);

      std::unique_ptr<TH1> MakeScaledUncertaintyHist(const std::string& Name,
				     std::vector< std::pair<const TH1*, const TH1*> > HistVec ) const;

      TH1* MakeAbsolUncertaintyHist( const std::string& Name, const TH1* Hist );

      RooArgList createStatConstraintTerms( RooWorkspace* proto, 
					    std::vector<std::string>& constraintTerms, 
					    ParamHistFunc& paramHist, const TH1* uncertHist,
					    Constraint::Type type, Double_t minSigma );

      void ConfigureHistFactoryDataset(RooDataSet* obsData, TH1* nominal, RooWorkspace* proto,
				       std::vector<std::string> obsNameVec);
    
      std::vector<std::string> fSystToFix;
      std::map<std::string, double> fParamValues;
      double fNomLumi;
      double fLumiError;
      int fLowBin; 
      int fHighBin;    

    private:
    
      void GuessObsNameVec(const TH1* hist);
    
      std::vector<std::string> fObsNameVec;
      std::string fObsName;
      std::vector<std::string> fPreprocessFunctions;
    
      ClassDef(RooStats::HistFactory::HistoToWorkspaceFactoryFast,3)
    };
  
  }
}

#endif
