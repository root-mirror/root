// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss, Eckhard von Toerne
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodTMlpANN                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodTMlpANN
\ingroup TMVA

This is the TMVA TMultiLayerPerceptron interface class. It provides the
training and testing the ROOT internal MLP class in the TMVA framework.

Available learning methods:<br>

  - Stochastic
  - Batch
  - SteepestDescent
  - RibierePolak
  - FletcherReeves
  - BFGS

See the TMultiLayerPerceptron class description
for details on this ANN.
*/

#include "TMVA/MethodTMlpANN.h"

#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Tools.h"

#include "TLeaf.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMultiLayerPerceptron.h"
#include "ThreadLocalStorage.h"

#include <cstdlib>
#include <iostream>
#include <fstream>


using std::atoi;

// some additional TMlpANN options
const Bool_t EnforceNormalization__=kTRUE;

REGISTER_METHOD(TMlpANN)

ClassImp(TMVA::MethodTMlpANN);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

   TMVA::MethodTMlpANN::MethodTMlpANN( const TString& jobName,
                                       const TString& methodTitle,
                                       DataSetInfo& theData,
                                       const TString& theOption) :
   TMVA::MethodBase( jobName, Types::kTMlpANN, methodTitle, theData, theOption),
   fMLP(0),
   fLocalTrainingTree(0),
   fNcycles(100),
   fValidationFraction(0.5),
   fLearningMethod( "" )
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodTMlpANN::MethodTMlpANN( DataSetInfo& theData,
                                    const TString& theWeightFile) :
   TMVA::MethodBase( Types::kTMlpANN, theData, theWeightFile),
   fMLP(0),
   fLocalTrainingTree(0),
   fNcycles(100),
   fValidationFraction(0.5),
   fLearningMethod( "" )
{
}

////////////////////////////////////////////////////////////////////////////////
/// TMlpANN can handle classification with 2 classes

Bool_t TMVA::MethodTMlpANN::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses,
                                             UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// default initialisations

void TMVA::MethodTMlpANN::Init( void )
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodTMlpANN::~MethodTMlpANN( void )
{
   if (fMLP) delete fMLP;
}

////////////////////////////////////////////////////////////////////////////////
/// translates options from option string into TMlpANN language

void TMVA::MethodTMlpANN::CreateMLPOptions( TString layerSpec )
{
   fHiddenLayer = ":";

   while (layerSpec.Length()>0) {
      TString sToAdd="";
      if (layerSpec.First(',')<0) {
         sToAdd = layerSpec;
         layerSpec = "";
      }
      else {
         sToAdd = layerSpec(0,layerSpec.First(','));
         layerSpec = layerSpec(layerSpec.First(',')+1,layerSpec.Length());
      }
      int nNodes = 0;
      if (sToAdd.BeginsWith("N")) { sToAdd.Remove(0,1); nNodes = GetNvar(); }
      nNodes += atoi(sToAdd);
      fHiddenLayer = Form( "%s%i:", (const char*)fHiddenLayer, nNodes );
   }

   // set input vars
   std::vector<TString>::iterator itrVar    = (*fInputVars).begin();
   std::vector<TString>::iterator itrVarEnd = (*fInputVars).end();
   fMLPBuildOptions = "";
   for (; itrVar != itrVarEnd; ++itrVar) {
      if (EnforceNormalization__) fMLPBuildOptions += "@";
      TString myVar = *itrVar; ;
      fMLPBuildOptions += myVar;
      fMLPBuildOptions += ",";
   }
   fMLPBuildOptions.Chop(); // remove last ","

   // prepare final options for MLP kernel
   fMLPBuildOptions += fHiddenLayer;
   fMLPBuildOptions += "type";

   Log() << kINFO << "Use " << fNcycles << " training cycles" << Endl;
   Log() << kINFO << "Use configuration (nodes per hidden layer): " << fHiddenLayer << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string
///
/// know options:
///
///  - NCycles       <integer>    Number of training cycles (too many cycles could overtrain the network)
///  - HiddenLayers  <string>     Layout of the hidden layers (nodes per layer)
///     * specifications for each hidden layer are separated by comma
///     * for each layer the number of nodes can be either absolut (simply a number)
///          or relative to the number of input nodes to the neural net (N)
///     * there is always a single node in the output layer
///
///   example: a net with 6 input nodes and "Hiddenlayers=N-1,N-2" has 6,5,4,1 nodes in the
///   layers 1,2,3,4, respectively

void TMVA::MethodTMlpANN::DeclareOptions()
{
   DeclareOptionRef( fNcycles    = 200,       "NCycles",      "Number of training cycles" );
   DeclareOptionRef( fLayerSpec  = "N,N-1",   "HiddenLayers", "Specification of hidden layer architecture (N stands for number of variables; any integers may also be used)" );

   DeclareOptionRef( fValidationFraction = 0.5, "ValidationFraction",
                     "Fraction of events in training tree used for cross validation" );

   DeclareOptionRef( fLearningMethod = "Stochastic", "LearningMethod", "Learning method" );
   AddPreDefVal( TString("Stochastic") );
   AddPreDefVal( TString("Batch") );
   AddPreDefVal( TString("SteepestDescent") );
   AddPreDefVal( TString("RibierePolak") );
   AddPreDefVal( TString("FletcherReeves") );
   AddPreDefVal( TString("BFGS") );
}

////////////////////////////////////////////////////////////////////////////////
/// builds the neural network as specified by the user

void TMVA::MethodTMlpANN::ProcessOptions()
{
   CreateMLPOptions(fLayerSpec);

   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kFATAL << "Mechanism to ignore events with negative weights in training not available for method"
            << GetMethodTypeName()
            << " --> please remove \"IgnoreNegWeightsInTraining\" option from booking string."
            << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the value of the neural net for the current event

Double_t TMVA::MethodTMlpANN::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   const Event* ev = GetEvent();
   TTHREAD_TLS_DECL_ARG(Double_t*, d, new Double_t[Data()->GetNVariables()]);

   for (UInt_t ivar = 0; ivar<Data()->GetNVariables(); ivar++) {
      d[ivar] = (Double_t)ev->GetValue(ivar);
   }
   Double_t mvaVal = fMLP->Evaluate(0,d);

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return mvaVal;
}

////////////////////////////////////////////////////////////////////////////////
/// performs TMlpANN training
/// available learning methods:
///
///  - TMultiLayerPerceptron::kStochastic
///  - TMultiLayerPerceptron::kBatch
///  - TMultiLayerPerceptron::kSteepestDescent
///  - TMultiLayerPerceptron::kRibierePolak
///  - TMultiLayerPerceptron::kFletcherReeves
///  - TMultiLayerPerceptron::kBFGS
///
/// TMultiLayerPerceptron wants test and training tree at once
/// so merge the training and testing trees from the MVA factory first:

void TMVA::MethodTMlpANN::Train( void )
{
   Int_t type;
   Float_t weight;
   const Long_t basketsize = 128000;
   Float_t* vArr = new Float_t[GetNvar()];

   TTree *localTrainingTree = new TTree( "TMLPtrain", "Local training tree for TMlpANN" );
   localTrainingTree->Branch( "type",       &type,        "type/I",        basketsize );
   localTrainingTree->Branch( "weight",     &weight,      "weight/F",      basketsize );

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      const char* myVar = GetInternalVarName(ivar).Data();
      localTrainingTree->Branch( myVar, &vArr[ivar], Form("Var%02i/F", ivar), basketsize );
   }

   for (UInt_t ievt=0; ievt<Data()->GetNEvents(); ievt++) {
      const Event *ev = GetEvent(ievt);
      for (UInt_t i=0; i<GetNvar(); i++) {
         vArr[i] = ev->GetValue( i );
      }
      type   = DataInfo().IsSignal( ev ) ? 1 : 0;
      weight = ev->GetWeight();
      localTrainingTree->Fill();
   }

   // These are the event lists for the mlp train method
   // first events in the tree are for training
   // the rest for internal testing (cross validation)...
   // NOTE: the training events are ordered: first part is signal, second part background
   TString trainList = "Entry$<";
   trainList += 1.0-fValidationFraction;
   trainList += "*";
   trainList += (Int_t)Data()->GetNEvtSigTrain();
   trainList += " || (Entry$>";
   trainList += (Int_t)Data()->GetNEvtSigTrain();
   trainList += " && Entry$<";
   trainList += (Int_t)(Data()->GetNEvtSigTrain() + (1.0 - fValidationFraction)*Data()->GetNEvtBkgdTrain());
   trainList += ")";
   TString testList  = TString("!(") + trainList + ")";

   // print the requirements
   Log() << kHEADER << "Requirement for training   events: \"" << trainList << "\"" << Endl;
   Log() << kINFO << "Requirement for validation events: \"" << testList << "\"" << Endl;

   // localTrainingTree->Print();

   // create NN
   if (fMLP != 0) { delete fMLP; fMLP = 0; }
   fMLP = new TMultiLayerPerceptron( fMLPBuildOptions.Data(),
                                     localTrainingTree,
                                     trainList,
                                     testList );
   fMLP->SetEventWeight( "weight" );

   // set learning method
   TMultiLayerPerceptron::ELearningMethod learningMethod = TMultiLayerPerceptron::kStochastic;

   fLearningMethod.ToLower();
   if      (fLearningMethod == "stochastic"      ) learningMethod = TMultiLayerPerceptron::kStochastic;
   else if (fLearningMethod == "batch"           ) learningMethod = TMultiLayerPerceptron::kBatch;
   else if (fLearningMethod == "steepestdescent" ) learningMethod = TMultiLayerPerceptron::kSteepestDescent;
   else if (fLearningMethod == "ribierepolak"    ) learningMethod = TMultiLayerPerceptron::kRibierePolak;
   else if (fLearningMethod == "fletcherreeves"  ) learningMethod = TMultiLayerPerceptron::kFletcherReeves;
   else if (fLearningMethod == "bfgs"            ) learningMethod = TMultiLayerPerceptron::kBFGS;
   else {
      Log() << kFATAL << "Unknown Learning Method: \"" << fLearningMethod << "\"" << Endl;
   }
   fMLP->SetLearningMethod( learningMethod );

   // train NN
   fMLP->Train(fNcycles, "" ); //"text,update=50" );

   // write weights to File;
   // this is not nice, but fMLP gets deleted at the end of Train()
   delete localTrainingTree;
   delete [] vArr;
}

////////////////////////////////////////////////////////////////////////////////
/// write weights to xml file

void TMVA::MethodTMlpANN::AddWeightsXMLTo( void* parent ) const
{
   // first the architecture
   void *wght = gTools().AddChild(parent, "Weights");
   void* arch = gTools().AddChild( wght, "Architecture" );
   gTools().AddAttr( arch, "BuildOptions", fMLPBuildOptions.Data() );

   // dump weights first in temporary txt file, read from there into xml
   const TString tmpfile=GetWeightFileDir()+"/TMlp.nn.weights.temp";
   fMLP->DumpWeights( tmpfile.Data() );
   std::ifstream inf( tmpfile.Data() );
   char temp[256];
   TString data("");
   void *ch=NULL;
   while (inf.getline(temp,256)) {
      TString dummy(temp);
      //std::cout << dummy << std::endl; // remove annoying debug printout with std::cout
      if (dummy.BeginsWith('#')) {
         if (ch!=0) gTools().AddRawLine( ch, data.Data() );
         dummy = dummy.Strip(TString::kLeading, '#');
         dummy = dummy(0,dummy.First(' '));
         ch = gTools().AddChild(wght, dummy);
         data.Resize(0);
         continue;
      }
      data += (dummy + " ");
   }
   if (ch != 0) gTools().AddRawLine( ch, data.Data() );

   inf.close();
}

////////////////////////////////////////////////////////////////////////////////
/// rebuild temporary textfile from xml weightfile and load this
/// file into MLP

void  TMVA::MethodTMlpANN::ReadWeightsFromXML( void* wghtnode )
{
   void* ch = gTools().GetChild(wghtnode);
   gTools().ReadAttr( ch, "BuildOptions", fMLPBuildOptions );

   ch = gTools().GetNextChild(ch);
   const TString fname = GetWeightFileDir()+"/TMlp.nn.weights.temp";
   std::ofstream fout( fname.Data() );
   double temp1=0,temp2=0;
   while (ch) {
      const char* nodecontent = gTools().GetContent(ch);
      std::stringstream content(nodecontent);
      if (strcmp(gTools().GetName(ch),"input")==0) {
         fout << "#input normalization" << std::endl;
         while ((content >> temp1) &&(content >> temp2)) {
            fout << temp1 << " " << temp2 << std::endl;
         }
      }
      if (strcmp(gTools().GetName(ch),"output")==0) {
         fout << "#output normalization" << std::endl;
         while ((content >> temp1) &&(content >> temp2)) {
            fout << temp1 << " " << temp2 << std::endl;
         }
      }
      if (strcmp(gTools().GetName(ch),"neurons")==0) {
         fout << "#neurons weights" << std::endl;
         while (content >> temp1) {
            fout << temp1 << std::endl;
         }
      }
      if (strcmp(gTools().GetName(ch),"synapses")==0) {
         fout << "#synapses weights" ;
         while (content >> temp1) {
            fout << std::endl << temp1 ;
         }
      }
      ch = gTools().GetNextChild(ch);
   }
   fout.close();;

   // Here we create a dummy tree necessary to create a minimal NN
   // to be used for testing, evaluation and application
   TTHREAD_TLS_DECL_ARG(Double_t*, d, new Double_t[Data()->GetNVariables()]);
   TTHREAD_TLS(Int_t) type;

   gROOT->cd();
   TTree * dummyTree = new TTree("dummy","Empty dummy tree", 1);
   for (UInt_t ivar = 0; ivar<Data()->GetNVariables(); ivar++) {
      TString vn = DataInfo().GetVariableInfo(ivar).GetInternalName();
      dummyTree->Branch(Form("%s",vn.Data()), d+ivar, Form("%s/D",vn.Data()));
   }
   dummyTree->Branch("type", &type, "type/I");

   if (fMLP != 0) { delete fMLP; fMLP = 0; }
   fMLP = new TMultiLayerPerceptron( fMLPBuildOptions.Data(), dummyTree );
   fMLP->LoadWeights( fname );
}

////////////////////////////////////////////////////////////////////////////////
/// read weights from stream
/// since the MLP can not read from the stream, we
/// 1st: write the weights to temporary file

void  TMVA::MethodTMlpANN::ReadWeightsFromStream( std::istream& istr )
{
   std::ofstream fout( "./TMlp.nn.weights.temp" );
   fout << istr.rdbuf();
   fout.close();
   // 2nd: load the weights from the temporary file into the MLP
   // the MLP is already build
   Log() << kINFO << "Load TMLP weights into " << fMLP << Endl;

   Double_t* d = new Double_t[Data()->GetNVariables()] ;
   Int_t type;
   gROOT->cd();
   TTree * dummyTree = new TTree("dummy","Empty dummy tree", 1);
   for (UInt_t ivar = 0; ivar<Data()->GetNVariables(); ivar++) {
      TString vn = DataInfo().GetVariableInfo(ivar).GetLabel();
      dummyTree->Branch(Form("%s",vn.Data()), d+ivar, Form("%s/D",vn.Data()));
   }
   dummyTree->Branch("type", &type, "type/I");

   if (fMLP != 0) { delete fMLP; fMLP = 0; }
   fMLP = new TMultiLayerPerceptron( fMLPBuildOptions.Data(), dummyTree );

   fMLP->LoadWeights( "./TMlp.nn.weights.temp" );
   // here we can delete the temporary file
   // how?
   delete [] d;
}

////////////////////////////////////////////////////////////////////////////////
/// create reader class for classifier -> overwrites base class function
/// create specific class for TMultiLayerPerceptron

void TMVA::MethodTMlpANN::MakeClass( const TString& theClassFileName ) const
{
   // the default consists of
   TString classFileName = "";
   if (theClassFileName == "")
      classFileName = GetWeightFileDir() + "/" + GetJobName() + "_" + GetMethodName() + ".class";
   else
      classFileName = theClassFileName;

   classFileName.ReplaceAll(".class","");
   Log() << kINFO << "Creating specific (TMultiLayerPerceptron) standalone response class: " << classFileName << Endl;
   fMLP->Export( classFileName.Data() );
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response
/// nothing to do here - all taken care of by TMultiLayerPerceptron

void TMVA::MethodTMlpANN::MakeClassSpecific( std::ostream& /*fout*/, const TString& /*className*/ ) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodTMlpANN::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "This feed-forward multilayer perceptron neural network is the " << Endl;
   Log() << "standard implementation distributed with ROOT (class TMultiLayerPerceptron)." << Endl;
   Log() << Endl;
   Log() << "Detailed information is available here:" << Endl;
   if (gConfig().WriteOptionsReference()) {
      Log() << "<a href=\"http://root.cern.ch/root/html/TMultiLayerPerceptron.html\">";
      Log() << "http://root.cern.ch/root/html/TMultiLayerPerceptron.html</a>" << Endl;
   }
   else Log() << "http://root.cern.ch/root/html/TMultiLayerPerceptron.html" << Endl;
   Log() << Endl;
}
