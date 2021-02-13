// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Matt Jachowski, Jan Therhaag, Jiahang Zhong

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodANNBase                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Artificial neural network base class for the discrimination of signal     *
 *      from background.                                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Krzysztof Danielowski <danielow@cern.ch>       - IFJ & AGH, Poland        *
 *      Andreas Hoecker       <Andreas.Hocker@cern.ch> - CERN, Switzerland        *
 *      Matt Jachowski        <jachowski@stanford.edu> - Stanford University, USA *
 *      Kamil Kraszewski      <kalq@cern.ch>           - IFJ & UJ, Poland         *
 *      Maciej Kruk           <mkruk@cern.ch>          - IFJ & AGH, Poland        *
 *      Peter Speckmayer      <peter.speckmayer@cern.ch> - CERN, Switzerland      *
 *      Joerg Stelzer         <stelzer@cern.ch>        - DESY, Germany            *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Jiahang Zhong         <Jiahang.Zhong@cern.ch>  - Academia Sinica, Taipei  *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodANNBase
\ingroup TMVA

Base class for all TMVA methods using artificial neural networks.

*/

#include "TMVA/MethodBase.h"

#include "TMVA/Configurable.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/MethodANNBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/TNeuron.h"
#include "TMVA/TSynapse.h"
#include "TMVA/TActivationChooser.h"
#include "TMVA/TActivationTanh.h"
#include "TMVA/Types.h"
#include "TMVA/Tools.h"
#include "TMVA/TNeuronInputChooser.h"
#include "TMVA/Ranking.h"
#include "TMVA/Version.h"

#include "TString.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1.h"
#include "TMath.h"
#include "TMatrixT.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#if __cplusplus > 199711L
#include <atomic>
#endif


using std::vector;

ClassImp(TMVA::MethodANNBase);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor
/// Note: Right now it is an option to choose the neuron input function,
/// but only the input function "sum" leads to weight convergence --
/// otherwise the weights go to nan and lead to an ABORT.

TMVA::MethodANNBase::MethodANNBase( const TString& jobName,
                                    Types::EMVA methodType,
                                    const TString& methodTitle,
                                    DataSetInfo& theData,
                                    const TString& theOption )
: TMVA::MethodBase( jobName, methodType, methodTitle, theData, theOption)
   , fEstimator(kMSE)
   , fUseRegulator(kFALSE)
   , fRandomSeed(0)
{
   InitANNBase();

   DeclareOptions();
}

////////////////////////////////////////////////////////////////////////////////
/// construct the Method from the weight file

TMVA::MethodANNBase::MethodANNBase( Types::EMVA methodType,
                                    DataSetInfo& theData,
                                    const TString& theWeightFile)
   : TMVA::MethodBase( methodType, theData, theWeightFile)
   , fEstimator(kMSE)
   , fUseRegulator(kFALSE)
   , fRandomSeed(0)
{
   InitANNBase();

   DeclareOptions();
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string
/// here the options valid for ALL MVA methods are declared.
///
/// know options:
///
///  - NCycles=xx              :the number of training cycles
///  - Normalize=kTRUE,kFALSe  :if normalised in put variables should be used
///  - HiddenLayser="N-1,N-2"  :the specification of the hidden layers
///  - NeuronType=sigmoid,tanh,radial,linar  : the type of activation function
///    used at the neuron

void TMVA::MethodANNBase::DeclareOptions()
{
   DeclareOptionRef( fNcycles    = 500,       "NCycles",         "Number of training cycles" );
   DeclareOptionRef( fLayerSpec  = "N,N-1",   "HiddenLayers",    "Specification of hidden layer architecture" );
   DeclareOptionRef( fNeuronType = "sigmoid", "NeuronType",      "Neuron activation function type" );
   DeclareOptionRef( fRandomSeed = 1, "RandomSeed", "Random seed for initial synapse weights (0 means unique seed for each run; default value '1')");

   DeclareOptionRef(fEstimatorS="MSE", "EstimatorType",
                    "MSE (Mean Square Estimator) for Gaussian Likelihood or CE(Cross-Entropy) for Bernoulli Likelihood" ); //zjh
   AddPreDefVal(TString("MSE"));  //zjh
   AddPreDefVal(TString("CE"));   //zjh


   TActivationChooser aChooser;
   std::vector<TString>* names = aChooser.GetAllActivationNames();
   Int_t nTypes = names->size();
   for (Int_t i = 0; i < nTypes; i++)
      AddPreDefVal(names->at(i));
   delete names;

   DeclareOptionRef(fNeuronInputType="sum", "NeuronInputType","Neuron input function type");
   TNeuronInputChooser iChooser;
   names = iChooser.GetAllNeuronInputNames();
   nTypes = names->size();
   for (Int_t i = 0; i < nTypes; i++) AddPreDefVal(names->at(i));
   delete names;
}


////////////////////////////////////////////////////////////////////////////////
/// do nothing specific at this moment

void TMVA::MethodANNBase::ProcessOptions()
{
   if      ( DoRegression() || DoMulticlass())  fEstimatorS = "MSE";    //zjh
   else    fEstimatorS = "CE" ;                                         //hhv
   if      (fEstimatorS == "MSE" )  fEstimator = kMSE;
   else if (fEstimatorS == "CE")    fEstimator = kCE;      //zjh
   std::vector<Int_t>* layout = ParseLayoutString(fLayerSpec);
   BuildNetwork(layout);
   delete layout;
}

////////////////////////////////////////////////////////////////////////////////
/// parse layout specification string and return a vector, each entry
/// containing the number of neurons to go in each successive layer

std::vector<Int_t>* TMVA::MethodANNBase::ParseLayoutString(TString layerSpec)
{
   std::vector<Int_t>* layout = new std::vector<Int_t>();
   layout->push_back((Int_t)GetNvar());
   while(layerSpec.Length()>0) {
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
      if (sToAdd.BeginsWith("n") || sToAdd.BeginsWith("N")) { sToAdd.Remove(0,1); nNodes = GetNvar(); }
      nNodes += atoi(sToAdd);
      layout->push_back(nNodes);
   }
   if( DoRegression() )
      layout->push_back( DataInfo().GetNTargets() );  // one output node for each target
   else if( DoMulticlass() )
      layout->push_back( DataInfo().GetNClasses() );  // one output node for each class
   else
      layout->push_back(1);  // one output node (for signal/background classification)

   int n = 0;
   for( std::vector<Int_t>::iterator it = layout->begin(); it != layout->end(); ++it ){
      n++;
   }

   return layout;
}

////////////////////////////////////////////////////////////////////////////////
/// initialize ANNBase object

void TMVA::MethodANNBase::InitANNBase()
{
   fNetwork         = NULL;
   frgen            = NULL;
   fActivation      = NULL;
   fOutput          = NULL; //zjh
   fIdentity        = NULL;
   fInputCalculator = NULL;
   fSynapses        = NULL;
   fEstimatorHistTrain = NULL;
   fEstimatorHistTest  = NULL;

   // reset monitoring histogram vectors
   fEpochMonHistS.clear();
   fEpochMonHistB.clear();
   fEpochMonHistW.clear();

   // these will be set in BuildNetwork()
   fInputLayer = NULL;
   fOutputNeurons.clear();

   frgen = new TRandom3(fRandomSeed);

   fSynapses = new TObjArray();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodANNBase::~MethodANNBase()
{
   DeleteNetwork();
}

////////////////////////////////////////////////////////////////////////////////
/// delete/clear network

void TMVA::MethodANNBase::DeleteNetwork()
{
   if (fNetwork != NULL) {
      TObjArray *layer;
      Int_t numLayers = fNetwork->GetEntriesFast();
      for (Int_t i = 0; i < numLayers; i++) {
         layer = (TObjArray*)fNetwork->At(i);
         DeleteNetworkLayer(layer);
      }
      delete fNetwork;
   }

   if (frgen != NULL)            delete frgen;
   if (fActivation != NULL)      delete fActivation;
   if (fOutput != NULL)          delete fOutput;  //zjh
   if (fIdentity != NULL)        delete fIdentity;
   if (fInputCalculator != NULL) delete fInputCalculator;
   if (fSynapses != NULL)        delete fSynapses;

   fNetwork         = NULL;
   frgen            = NULL;
   fActivation      = NULL;
   fOutput          = NULL; //zjh
   fIdentity        = NULL;
   fInputCalculator = NULL;
   fSynapses        = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// delete a network layer

void TMVA::MethodANNBase::DeleteNetworkLayer( TObjArray*& layer )
{
   TNeuron* neuron;
   Int_t numNeurons = layer->GetEntriesFast();
   for (Int_t i = 0; i < numNeurons; i++) {
      neuron = (TNeuron*)layer->At(i);
      neuron->DeletePreLinks();
      delete neuron;
   }
   delete layer;
}

////////////////////////////////////////////////////////////////////////////////
/// build network given a layout (number of neurons in each layer)
/// and optional weights array

void TMVA::MethodANNBase::BuildNetwork( std::vector<Int_t>* layout, std::vector<Double_t>* weights, Bool_t fromFile )
{
   if (fEstimatorS == "MSE")  fEstimator = kMSE;    //zjh
   else if (fEstimatorS == "CE")    fEstimator = kCE;      //zjh
   else Log()<<kWARNING<<"fEstimator="<<fEstimator<<"\tfEstimatorS="<<fEstimatorS<<Endl;
   if (fEstimator!=kMSE && fEstimator!=kCE) Log()<<kWARNING<<"Estimator type unspecified \t"<<Endl; //zjh


   Log() << kHEADER << "Building Network. " << Endl;

   DeleteNetwork();
   InitANNBase();

   // set activation and input functions
   TActivationChooser aChooser;
   fActivation = aChooser.CreateActivation(fNeuronType);
   fIdentity   = aChooser.CreateActivation("linear");
   if (fEstimator==kMSE)  fOutput = aChooser.CreateActivation("linear");  //zjh
   else if (fEstimator==kCE)   fOutput = aChooser.CreateActivation("sigmoid"); //zjh
   TNeuronInputChooser iChooser;
   fInputCalculator = iChooser.CreateNeuronInput(fNeuronInputType);

   fNetwork = new TObjArray();
   fRegulatorIdx.clear();
   fRegulators.clear();
   BuildLayers( layout, fromFile );

   // cache input layer and output neuron for fast access
   fInputLayer   = (TObjArray*)fNetwork->At(0);
   TObjArray* outputLayer = (TObjArray*)fNetwork->At(fNetwork->GetEntriesFast()-1);
   fOutputNeurons.clear();
   for (Int_t i = 0; i < outputLayer->GetEntries(); i++) {
      fOutputNeurons.push_back( (TNeuron*)outputLayer->At(i) );
   }

   if (weights == NULL) InitWeights();
   else                 ForceWeights(weights);
}

////////////////////////////////////////////////////////////////////////////////
/// build the network layers

void TMVA::MethodANNBase::BuildLayers( std::vector<Int_t>* layout, Bool_t fromFile )
{
   TObjArray* curLayer;
   TObjArray* prevLayer = NULL;

   Int_t numLayers = layout->size();

   for (Int_t i = 0; i < numLayers; i++) {
      curLayer = new TObjArray();
      BuildLayer(layout->at(i), curLayer, prevLayer, i, numLayers, fromFile);
      prevLayer = curLayer;
      fNetwork->Add(curLayer);
   }

   // cache pointers to synapses for fast access, the order matters
   for (Int_t i = 0; i < numLayers; i++) {
      TObjArray* layer = (TObjArray*)fNetwork->At(i);
      Int_t numNeurons = layer->GetEntriesFast();
      if (i!=0 && i!=numLayers-1) fRegulators.push_back(0.);  //zjh
      for (Int_t j = 0; j < numNeurons; j++) {
         if (i==0) fRegulators.push_back(0.);//zjh
         TNeuron* neuron = (TNeuron*)layer->At(j);
         Int_t numSynapses = neuron->NumPostLinks();
         for (Int_t k = 0; k < numSynapses; k++) {
            TSynapse* synapse = neuron->PostLinkAt(k);
            fSynapses->Add(synapse);
            fRegulatorIdx.push_back(fRegulators.size()-1);//zjh
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// build a single layer with neurons and synapses connecting this
/// layer to the previous layer

void TMVA::MethodANNBase::BuildLayer( Int_t numNeurons, TObjArray* curLayer,
                                      TObjArray* prevLayer, Int_t layerIndex,
                                      Int_t numLayers, Bool_t fromFile )
{
   TNeuron* neuron;
   for (Int_t j = 0; j < numNeurons; j++) {
      if (fromFile && (layerIndex != numLayers-1) && (j==numNeurons-1)){
         neuron = new TNeuron();
         neuron->SetActivationEqn(fIdentity);
         neuron->SetBiasNeuron();
         neuron->ForceValue(1.0);
         curLayer->Add(neuron);
      }
      else {
         neuron = new TNeuron();
         neuron->SetInputCalculator(fInputCalculator);

         // input layer
         if (layerIndex == 0) {
            neuron->SetActivationEqn(fIdentity);
            neuron->SetInputNeuron();
         }
         else {
            // output layer
            if (layerIndex == numLayers-1) {
               neuron->SetOutputNeuron();
               neuron->SetActivationEqn(fOutput);     //zjh
            }
            // hidden layers
            else neuron->SetActivationEqn(fActivation);
            AddPreLinks(neuron, prevLayer);
         }

         curLayer->Add(neuron);
      }
   }

   // add bias neutron (except to output layer)
   if(!fromFile){
      if (layerIndex != numLayers-1) {
         neuron = new TNeuron();
         neuron->SetActivationEqn(fIdentity);
         neuron->SetBiasNeuron();
         neuron->ForceValue(1.0);
         curLayer->Add(neuron);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// add synapses connecting a neuron to its preceding layer

void TMVA::MethodANNBase::AddPreLinks(TNeuron* neuron, TObjArray* prevLayer)
{
   TSynapse* synapse;
   int numNeurons = prevLayer->GetEntriesFast();
   TNeuron* preNeuron;

   for (Int_t i = 0; i < numNeurons; i++) {
      preNeuron = (TNeuron*)prevLayer->At(i);
      synapse = new TSynapse();
      synapse->SetPreNeuron(preNeuron);
      synapse->SetPostNeuron(neuron);
      preNeuron->AddPostLink(synapse);
      neuron->AddPreLink(synapse);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// initialize the synapse weights randomly

void TMVA::MethodANNBase::InitWeights()
{
  PrintMessage("Initializing weights");

   // init synapse weights
   Int_t numSynapses = fSynapses->GetEntriesFast();
   TSynapse* synapse;
   for (Int_t i = 0; i < numSynapses; i++) {
      synapse = (TSynapse*)fSynapses->At(i);
      synapse->SetWeight(4.0*frgen->Rndm() - 2.0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// force the synapse weights

void TMVA::MethodANNBase::ForceWeights(std::vector<Double_t>* weights)
{
   PrintMessage("Forcing weights");

   Int_t numSynapses = fSynapses->GetEntriesFast();
   TSynapse* synapse;
   for (Int_t i = 0; i < numSynapses; i++) {
      synapse = (TSynapse*)fSynapses->At(i);
      synapse->SetWeight(weights->at(i));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// force the input values of the input neurons
/// force the value for each input neuron

void TMVA::MethodANNBase::ForceNetworkInputs( const Event* ev, Int_t ignoreIndex)
{
   Double_t x;
   TNeuron* neuron;

   //   const Event* ev = GetEvent();
   for (UInt_t j = 0; j < GetNvar(); j++) {

      x = (j != (UInt_t)ignoreIndex)?ev->GetValue(j):0;

      neuron = GetInputNeuron(j);
      neuron->ForceValue(x);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// calculate input values to each neuron

void TMVA::MethodANNBase::ForceNetworkCalculations()
{
   TObjArray* curLayer;
   TNeuron* neuron;
   Int_t numLayers = fNetwork->GetEntriesFast();
   Int_t numNeurons;

   for (Int_t i = 0; i < numLayers; i++) {
      curLayer = (TObjArray*)fNetwork->At(i);
      numNeurons = curLayer->GetEntriesFast();

      for (Int_t j = 0; j < numNeurons; j++) {
         neuron = (TNeuron*) curLayer->At(j);
         neuron->CalculateValue();
         neuron->CalculateActivationValue();

      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print messages, turn off printing by setting verbose and debug flag appropriately

void TMVA::MethodANNBase::PrintMessage(TString message, Bool_t force) const
{
   if (Verbose() || Debug() || force) Log() << kINFO << message << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// wait for keyboard input, for debugging

void TMVA::MethodANNBase::WaitForKeyboard()
{
   std::string dummy;
   Log() << kINFO << "***Type anything to continue (q to quit): ";
   std::getline(std::cin, dummy);
   if (dummy == "q" || dummy == "Q") {
      PrintMessage( "quit" );
      delete this;
      exit(0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print network representation, for debugging

void TMVA::MethodANNBase::PrintNetwork() const
{
   if (!Debug()) return;

   Log() << kINFO << Endl;
   PrintMessage( "Printing network " );
   Log() << kINFO << "-------------------------------------------------------------------" << Endl;

   TObjArray* curLayer;
   Int_t numLayers = fNetwork->GetEntriesFast();

   for (Int_t i = 0; i < numLayers; i++) {

      curLayer = (TObjArray*)fNetwork->At(i);
      Int_t numNeurons = curLayer->GetEntriesFast();

      Log() << kINFO << "Layer #" << i << " (" << numNeurons << " neurons):" << Endl;
      PrintLayer( curLayer );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print a single layer, for debugging

void TMVA::MethodANNBase::PrintLayer(TObjArray* layer) const
{
   Int_t numNeurons = layer->GetEntriesFast();
   TNeuron* neuron;

   for (Int_t j = 0; j < numNeurons; j++) {
      neuron = (TNeuron*) layer->At(j);
      Log() << kINFO << "\tNeuron #" << j << " (LinksIn: " << neuron->NumPreLinks()
            << " , LinksOut: " << neuron->NumPostLinks() << ")" << Endl;
      PrintNeuron( neuron );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print a neuron, for debugging

void TMVA::MethodANNBase::PrintNeuron(TNeuron* neuron) const
{
   Log() << kINFO
         << "\t\tValue:\t"     << neuron->GetValue()
         << "\t\tActivation: " << neuron->GetActivationValue()
         << "\t\tDelta: "      << neuron->GetDelta() << Endl;
   Log() << kINFO << "\t\tActivationEquation:\t";
   neuron->PrintActivationEqn();
   Log() << kINFO << "\t\tLinksIn:" << Endl;
   neuron->PrintPreLinks();
   Log() << kINFO << "\t\tLinksOut:" << Endl;
   neuron->PrintPostLinks();
}

////////////////////////////////////////////////////////////////////////////////
/// get the mva value generated by the NN

Double_t TMVA::MethodANNBase::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   TNeuron* neuron;

   TObjArray* inputLayer = (TObjArray*)fNetwork->At(0);

   const Event * ev = GetEvent();

   for (UInt_t i = 0; i < GetNvar(); i++) {
      neuron = (TNeuron*)inputLayer->At(i);
      neuron->ForceValue( ev->GetValue(i) );
   }
   ForceNetworkCalculations();

   // check the output of the network
   TObjArray* outputLayer = (TObjArray*)fNetwork->At( fNetwork->GetEntriesFast()-1 );
   neuron = (TNeuron*)outputLayer->At(0);

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return neuron->GetActivationValue();
}

////////////////////////////////////////////////////////////////////////////////
/// get the regression value generated by the NN

const std::vector<Float_t> &TMVA::MethodANNBase::GetRegressionValues()
{
   TNeuron* neuron;

   TObjArray* inputLayer = (TObjArray*)fNetwork->At(0);

   const Event * ev = GetEvent();

   for (UInt_t i = 0; i < GetNvar(); i++) {
      neuron = (TNeuron*)inputLayer->At(i);
      neuron->ForceValue( ev->GetValue(i) );
   }
   ForceNetworkCalculations();

   // check the output of the network
   TObjArray* outputLayer = (TObjArray*)fNetwork->At( fNetwork->GetEntriesFast()-1 );

   if (fRegressionReturnVal == NULL) fRegressionReturnVal = new std::vector<Float_t>();
   fRegressionReturnVal->clear();

   Event * evT = new Event(*ev);
   UInt_t ntgts = outputLayer->GetEntriesFast();
   for (UInt_t itgt = 0; itgt < ntgts; itgt++) {
      evT->SetTarget(itgt,((TNeuron*)outputLayer->At(itgt))->GetActivationValue());
   }

   const Event* evT2 = GetTransformationHandler().InverseTransform( evT );
   for (UInt_t itgt = 0; itgt < ntgts; itgt++) {
      fRegressionReturnVal->push_back( evT2->GetTarget(itgt) );
   }

   delete evT;

   return *fRegressionReturnVal;
}

////////////////////////////////////////////////////////////////////////////////
/// get the multiclass classification values generated by the NN

const std::vector<Float_t> &TMVA::MethodANNBase::GetMulticlassValues()
{
   TNeuron* neuron;

   TObjArray* inputLayer = (TObjArray*)fNetwork->At(0);

   const Event * ev = GetEvent();

   for (UInt_t i = 0; i < GetNvar(); i++) {
      neuron = (TNeuron*)inputLayer->At(i);
      neuron->ForceValue( ev->GetValue(i) );
   }
   ForceNetworkCalculations();

   // check the output of the network

   if (fMulticlassReturnVal == NULL) fMulticlassReturnVal = new std::vector<Float_t>();
   fMulticlassReturnVal->clear();
   std::vector<Float_t> temp;

   UInt_t nClasses = DataInfo().GetNClasses();
   for (UInt_t icls = 0; icls < nClasses; icls++) {
      temp.push_back(GetOutputNeuron( icls )->GetActivationValue() );
   }

   for(UInt_t iClass=0; iClass<nClasses; iClass++){
      Double_t norm = 0.0;
      for(UInt_t j=0;j<nClasses;j++){
         if(iClass!=j)
            norm+=exp(temp[j]-temp[iClass]);
      }
      (*fMulticlassReturnVal).push_back(1.0/(1.0+norm));
   }



   return *fMulticlassReturnVal;
}


////////////////////////////////////////////////////////////////////////////////
/// create XML description of ANN classifier

void TMVA::MethodANNBase::AddWeightsXMLTo( void* parent ) const
{
   Int_t numLayers = fNetwork->GetEntriesFast();
   void* wght = gTools().xmlengine().NewChild(parent, 0, "Weights");
   void* xmlLayout = gTools().xmlengine().NewChild(wght, 0, "Layout");
   gTools().xmlengine().NewAttr(xmlLayout, 0, "NLayers", gTools().StringFromInt(fNetwork->GetEntriesFast()) );
   TString weights = "";
   for (Int_t i = 0; i < numLayers; i++) {
      TObjArray* layer = (TObjArray*)fNetwork->At(i);
      Int_t numNeurons = layer->GetEntriesFast();
      void* layerxml = gTools().xmlengine().NewChild(xmlLayout, 0, "Layer");
      gTools().xmlengine().NewAttr(layerxml, 0, "Index",    gTools().StringFromInt(i) );
      gTools().xmlengine().NewAttr(layerxml, 0, "NNeurons", gTools().StringFromInt(numNeurons) );
      for (Int_t j = 0; j < numNeurons; j++) {
         TNeuron* neuron = (TNeuron*)layer->At(j);
         Int_t numSynapses = neuron->NumPostLinks();
         void* neuronxml = gTools().AddChild(layerxml, "Neuron");
         gTools().AddAttr(neuronxml, "NSynapses", gTools().StringFromInt(numSynapses) );
         if(numSynapses==0) continue;
         std::stringstream s("");
         s.precision( 16 );
         for (Int_t k = 0; k < numSynapses; k++) {
            TSynapse* synapse = neuron->PostLinkAt(k);
            s << std::scientific << synapse->GetWeight() << " ";
         }
         gTools().AddRawLine( neuronxml, s.str().c_str() );
      }
   }

   // if inverse hessian exists, write inverse hessian to weight file
   if( fInvHessian.GetNcols()>0 ){
      void* xmlInvHessian = gTools().xmlengine().NewChild(wght, 0, "InverseHessian");

      // get the matrix dimensions
      Int_t nElements = fInvHessian.GetNoElements();
      Int_t nRows     = fInvHessian.GetNrows();
      Int_t nCols     = fInvHessian.GetNcols();
      gTools().xmlengine().NewAttr(xmlInvHessian, 0, "NElements", gTools().StringFromInt(nElements) );
      gTools().xmlengine().NewAttr(xmlInvHessian, 0, "NRows", gTools().StringFromInt(nRows) );
      gTools().xmlengine().NewAttr(xmlInvHessian, 0, "NCols", gTools().StringFromInt(nCols) );

      // read in the matrix elements
      Double_t* elements = new Double_t[nElements+10];
      fInvHessian.GetMatrix2Array( elements );

      // store the matrix elements row-wise
      Int_t index = 0;
      for( Int_t row = 0; row < nRows; ++row ){
         void* xmlRow = gTools().xmlengine().NewChild(xmlInvHessian, 0, "Row");
         gTools().xmlengine().NewAttr(xmlRow, 0, "Index", gTools().StringFromInt(row) );

         // create the rows
         std::stringstream s("");
         s.precision( 16 );
         for( Int_t col = 0; col < nCols; ++col ){
            s << std::scientific << (*(elements+index)) << " ";
            ++index;
         }
         gTools().xmlengine().AddRawLine( xmlRow, s.str().c_str() );
      }
      delete[] elements;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// read MLP from xml weight file

void TMVA::MethodANNBase::ReadWeightsFromXML( void* wghtnode )
{
   // build the layout first
   Bool_t fromFile = kTRUE;
   std::vector<Int_t>* layout = new std::vector<Int_t>();

   void* xmlLayout = NULL;
   xmlLayout = gTools().GetChild(wghtnode, "Layout");
   if( !xmlLayout )
      xmlLayout = wghtnode;

   UInt_t nLayers;
   gTools().ReadAttr( xmlLayout, "NLayers", nLayers );
   layout->resize( nLayers );

   void* ch = gTools().xmlengine().GetChild(xmlLayout);
   UInt_t index;
   UInt_t nNeurons;
   while (ch) {
      gTools().ReadAttr( ch, "Index",   index   );
      gTools().ReadAttr( ch, "NNeurons", nNeurons );
      layout->at(index) = nNeurons;
      ch = gTools().GetNextChild(ch);
   }

   BuildNetwork( layout, NULL, fromFile );
   // use 'slow' (exact) TanH if processing old weigh file to ensure 100% compatible results
   // otherwise use the new default, the 'tast tanh' approximation
   if (GetTrainingTMVAVersionCode() < TMVA_VERSION(4,2,1) && fActivation->GetExpression().Contains("tanh")){
      TActivationTanh* act = dynamic_cast<TActivationTanh*>( fActivation );
      if (act) act->SetSlow();
   }

   // fill the weights of the synapses
   UInt_t nSyn;
   Float_t weight;
   ch = gTools().xmlengine().GetChild(xmlLayout);
   UInt_t iLayer = 0;
   while (ch) {  // layers
      TObjArray* layer = (TObjArray*)fNetwork->At(iLayer);
      gTools().ReadAttr( ch, "Index",   index   );
      gTools().ReadAttr( ch, "NNeurons", nNeurons );

      void* nodeN = gTools().GetChild(ch);
      UInt_t iNeuron = 0;
      while( nodeN ){ // neurons
         TNeuron *neuron = (TNeuron*)layer->At(iNeuron);
         gTools().ReadAttr( nodeN, "NSynapses", nSyn );
         if( nSyn > 0 ){
            const char* content = gTools().GetContent(nodeN);
            std::stringstream s(content);
            for (UInt_t iSyn = 0; iSyn<nSyn; iSyn++) { // synapses

               TSynapse* synapse = neuron->PostLinkAt(iSyn);
               s >> weight;
               //Log() << kWARNING << neuron << " " << weight <<  Endl;
               synapse->SetWeight(weight);
            }
         }
         nodeN = gTools().GetNextChild(nodeN);
         iNeuron++;
      }
      ch = gTools().GetNextChild(ch);
      iLayer++;
   }

   delete layout;

   void* xmlInvHessian = NULL;
   xmlInvHessian = gTools().GetChild(wghtnode, "InverseHessian");
   if( !xmlInvHessian )
      // no inverse hessian available
      return;

   fUseRegulator = kTRUE;

   Int_t nElements = 0;
   Int_t nRows     = 0;
   Int_t nCols     = 0;
   gTools().ReadAttr( xmlInvHessian, "NElements", nElements );
   gTools().ReadAttr( xmlInvHessian, "NRows", nRows );
   gTools().ReadAttr( xmlInvHessian, "NCols", nCols );

   // adjust the matrix dimensions
   fInvHessian.ResizeTo( nRows, nCols );

   // prepare an array to read in the values
   Double_t* elements;
   if (nElements > std::numeric_limits<int>::max()-100){
      Log() << kFATAL << "you tried to read a hessian matrix with " << nElements << " elements, --> too large, guess s.th. went wrong reading from the weight file" << Endl;
      return;
   } else {
      elements = new Double_t[nElements+10];
   }



   void* xmlRow = gTools().xmlengine().GetChild(xmlInvHessian);
   Int_t row = 0;
   index = 0;
   while (xmlRow) {  // rows
      gTools().ReadAttr( xmlRow, "Index",   row   );

      const char* content = gTools().xmlengine().GetNodeContent(xmlRow);

      std::stringstream s(content);
      for (Int_t iCol = 0; iCol<nCols; iCol++) { // columns
         s >> (*(elements+index));
         ++index;
      }
      xmlRow = gTools().xmlengine().GetNext(xmlRow);
      ++row;
   }

   fInvHessian.SetMatrixArray( elements );

   delete[] elements;
}

////////////////////////////////////////////////////////////////////////////////
/// destroy/clear the network then read it back in from the weights file

void TMVA::MethodANNBase::ReadWeightsFromStream( std::istream & istr)
{
   // delete network so we can reconstruct network from scratch

   TString dummy;

   // synapse weights
   Double_t weight;
   std::vector<Double_t>* weights = new std::vector<Double_t>();
   istr>> dummy;
   while (istr>> dummy >> weight) weights->push_back(weight); // use w/ slower write-out

   ForceWeights(weights);


   delete weights;
}

////////////////////////////////////////////////////////////////////////////////
/// compute ranking of input variables by summing function of weights

const TMVA::Ranking* TMVA::MethodANNBase::CreateRanking()
{
   // create the ranking object
   fRanking = new Ranking( GetName(), "Importance" );

   TNeuron*  neuron;
   TSynapse* synapse;
   Double_t  importance, avgVal;
   TString varName;

   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) {

      neuron = GetInputNeuron(ivar);
      Int_t numSynapses = neuron->NumPostLinks();
      importance = 0;
      varName = GetInputVar(ivar); // fix this line

      // figure out average value of variable i
      Double_t meanS, meanB, rmsS, rmsB, xmin, xmax;
      Statistics( TMVA::Types::kTraining, varName,
                  meanS, meanB, rmsS, rmsB, xmin, xmax );

      avgVal = (TMath::Abs(meanS) + TMath::Abs(meanB))/2.0;
      double meanrms = (TMath::Abs(rmsS) + TMath::Abs(rmsB))/2.;
      if (avgVal<meanrms) avgVal = meanrms;
      if (IsNormalised()) avgVal = 0.5*(1 + gTools().NormVariable( avgVal, GetXmin( ivar ), GetXmax( ivar )));

      for (Int_t j = 0; j < numSynapses; j++) {
         synapse = neuron->PostLinkAt(j);
         importance += synapse->GetWeight() * synapse->GetWeight();
      }

      importance *= avgVal * avgVal;

      fRanking->AddRank( Rank( varName, importance ) );
   }

   return fRanking;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodANNBase::CreateWeightMonitoringHists( const TString& bulkname,
                                                       std::vector<TH1*>* hv ) const
{
   TH2F* hist;
   Int_t numLayers = fNetwork->GetEntriesFast();

   for (Int_t i = 0; i < numLayers-1; i++) {

      TObjArray* layer1 = (TObjArray*)fNetwork->At(i);
      TObjArray* layer2 = (TObjArray*)fNetwork->At(i+1);
      Int_t numNeurons1 = layer1->GetEntriesFast();
      Int_t numNeurons2 = layer2->GetEntriesFast();

      TString name = Form("%s%i%i", bulkname.Data(), i, i+1);
      hist = new TH2F(name + "", name + "",
                      numNeurons1, 0, numNeurons1, numNeurons2, 0, numNeurons2);

      for (Int_t j = 0; j < numNeurons1; j++) {

         TNeuron* neuron = (TNeuron*)layer1->At(j);
         Int_t numSynapses = neuron->NumPostLinks();

         for (Int_t k = 0; k < numSynapses; k++) {

            TSynapse* synapse = neuron->PostLinkAt(k);
            hist->SetBinContent(j+1, k+1, synapse->GetWeight());

         }
      }

      if (hv) hv->push_back( hist );
      else {
         hist->Write();
         delete hist;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// write histograms to file

void TMVA::MethodANNBase::WriteMonitoringHistosToFile() const
{
   PrintMessage(Form("Write special histos to file: %s", BaseDir()->GetPath()), kTRUE);

   if (fEstimatorHistTrain) fEstimatorHistTrain->Write();
   if (fEstimatorHistTest ) fEstimatorHistTest ->Write();

   // histograms containing weights for architecture plotting (used in macro "network.cxx")
   CreateWeightMonitoringHists( "weights_hist" );

   // now save all the epoch-wise monitoring information
#if __cplusplus > 199711L
   static std::atomic<int> epochMonitoringDirectoryNumber{0};
#else
   static int epochMonitoringDirectoryNumber = 0;
#endif
   int epochVal = epochMonitoringDirectoryNumber++;
   TDirectory* epochdir = NULL;
   if( epochVal == 0 )
      epochdir = BaseDir()->mkdir( "EpochMonitoring" );
   else
      epochdir = BaseDir()->mkdir( Form("EpochMonitoring_%4d",epochVal) );

   epochdir->cd();
   for (std::vector<TH1*>::const_iterator it = fEpochMonHistS.begin(); it != fEpochMonHistS.end(); ++it) {
      (*it)->Write();
      delete (*it);
   }
   for (std::vector<TH1*>::const_iterator it = fEpochMonHistB.begin(); it != fEpochMonHistB.end(); ++it) {
      (*it)->Write();
      delete (*it);
   }
   for (std::vector<TH1*>::const_iterator it = fEpochMonHistW.begin(); it != fEpochMonHistW.end(); ++it) {
      (*it)->Write();
      delete (*it);
   }
   BaseDir()->cd();
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response

void TMVA::MethodANNBase::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   Int_t numLayers = fNetwork->GetEntries();

   fout << std::endl;
   fout << "   double ActivationFnc(double x) const;" << std::endl;
   fout << "   double OutputActivationFnc(double x) const;" << std::endl;     //zjh
   fout << std::endl;
   int numNodesFrom = -1;
   for (Int_t lIdx = 0; lIdx < numLayers; lIdx++) {
      int numNodesTo = ((TObjArray*)fNetwork->At(lIdx))->GetEntries();
      if (numNodesFrom<0) { numNodesFrom=numNodesTo; continue; }
      fout << "   double fWeightMatrix" << lIdx-1  << "to" << lIdx << "[" << numNodesTo << "][" << numNodesFrom << "];";
      fout << "   // weight matrix from layer " << lIdx-1  << " to " << lIdx << std::endl;
      numNodesFrom = numNodesTo;
   }
   fout << std::endl;
   fout << "};" << std::endl;

   fout << std::endl;

   fout << "inline void " << className << "::Initialize()" << std::endl;
   fout << "{" << std::endl;
   fout << "   // build network structure" << std::endl;

   for (Int_t i = 0; i < numLayers-1; i++) {
      fout << "   // weight matrix from layer " << i  << " to " << i+1 << std::endl;
      TObjArray* layer = (TObjArray*)fNetwork->At(i);
      Int_t numNeurons = layer->GetEntriesFast();
      for (Int_t j = 0; j < numNeurons; j++) {
         TNeuron* neuron = (TNeuron*)layer->At(j);
         Int_t numSynapses = neuron->NumPostLinks();
         for (Int_t k = 0; k < numSynapses; k++) {
            TSynapse* synapse = neuron->PostLinkAt(k);
            fout << "   fWeightMatrix" << i  << "to" << i+1 << "[" << k << "][" << j << "] = " << synapse->GetWeight() << ";" << std::endl;
         }
      }
   }

   fout << "}" << std::endl;
   fout << std::endl;

   // writing of the GetMvaValue__ method
   fout << "inline double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   if (inputValues.size() != (unsigned int)" << ((TObjArray *)fNetwork->At(0))->GetEntries() - 1 << ") {"
        << std::endl;
   fout << "      std::cout << \"Input vector needs to be of size \" << "
        << ((TObjArray *)fNetwork->At(0))->GetEntries() - 1 << " << std::endl;" << std::endl;
   fout << "      return 0;" << std::endl;
   fout << "   }" << std::endl;
   fout << std::endl;
   for (Int_t lIdx = 1; lIdx < numLayers; lIdx++) {
      TObjArray *layer = (TObjArray *)fNetwork->At(lIdx);
      int numNodes = layer->GetEntries();
      fout << "   std::array<double, " << numNodes << "> fWeights" << lIdx << " {{}};" << std::endl;
   }
   for (Int_t lIdx = 1; lIdx < numLayers - 1; lIdx++) {
      fout << "   fWeights" << lIdx << ".back() = 1.;" << std::endl;
   }
   fout << std::endl;
   for (Int_t i = 0; i < numLayers - 1; i++) {
      fout << "   // layer " << i << " to " << i + 1 << std::endl;
      if (i + 1 == numLayers - 1) {
         fout << "   for (int o=0; o<" << ((TObjArray *)fNetwork->At(i + 1))->GetEntries() << "; o++) {" << std::endl;
      } else {
         fout << "   for (int o=0; o<" << ((TObjArray *)fNetwork->At(i + 1))->GetEntries() - 1 << "; o++) {"
              << std::endl;
      }
      if (0 == i) {
         fout << "      std::array<double, " << ((TObjArray *)fNetwork->At(i))->GetEntries()
              << "> buffer; // no need to initialise" << std::endl;
         fout << "      for (int i = 0; i<" << ((TObjArray *)fNetwork->At(i))->GetEntries() << " - 1; i++) {"
              << std::endl;
         fout << "         buffer[i] = fWeightMatrix" << i << "to" << i + 1 << "[o][i] * inputValues[i];" << std::endl;
         fout << "      } // loop over i" << std::endl;
         fout << "      buffer.back() = fWeightMatrix" << i << "to" << i + 1 << "[o]["
              << ((TObjArray *)fNetwork->At(i))->GetEntries() - 1 << "];" << std::endl;
      } else {
         fout << "      std::array<double, " << ((TObjArray *)fNetwork->At(i))->GetEntries()
              << "> buffer; // no need to initialise" << std::endl;
         fout << "      for (int i=0; i<" << ((TObjArray *)fNetwork->At(i))->GetEntries() << "; i++) {" << std::endl;
         fout << "         buffer[i] = fWeightMatrix" << i << "to" << i + 1 << "[o][i] * fWeights" << i << "[i];"
              << std::endl;
         fout << "      } // loop over i" << std::endl;
      }
      fout << "      for (int i=0; i<" << ((TObjArray *)fNetwork->At(i))->GetEntries() << "; i++) {" << std::endl;
      if (fNeuronInputType == "sum") {
         fout << "         fWeights" << i + 1 << "[o] += buffer[i];" << std::endl;
      } else if (fNeuronInputType == "sqsum") {
         fout << "         fWeights" << i + 1 << "[o] += buffer[i]*buffer[i];" << std::endl;
      } else { // fNeuronInputType == TNeuronInputChooser::kAbsSum
         fout << "         fWeights" << i + 1 << "[o] += fabs(buffer[i]);" << std::endl;
      }
      fout << "      } // loop over i" << std::endl;
      fout << "    } // loop over o" << std::endl;
      if (i + 1 == numLayers - 1) {
         fout << "   for (int o=0; o<" << ((TObjArray *)fNetwork->At(i + 1))->GetEntries() << "; o++) {" << std::endl;
      } else {
         fout << "   for (int o=0; o<" << ((TObjArray *)fNetwork->At(i + 1))->GetEntries() - 1 << "; o++) {"
              << std::endl;
      }
      if (i+1 != numLayers-1) // in the last layer no activation function is applied
         fout << "      fWeights" << i + 1 << "[o] = ActivationFnc(fWeights" << i + 1 << "[o]);" << std::endl;
      else
         fout << "      fWeights" << i + 1 << "[o] = OutputActivationFnc(fWeights" << i + 1 << "[o]);"
              << std::endl; // zjh
      fout << "   } // loop over o" << std::endl;
   }
   fout << std::endl;
   fout << "   return fWeights" << numLayers - 1 << "[0];" << std::endl;
   fout << "}" << std::endl;

   fout << std::endl;
   TString fncName = className+"::ActivationFnc";
   fActivation->MakeFunction(fout, fncName);
   fncName = className+"::OutputActivationFnc";   //zjh
   fOutput->MakeFunction(fout, fncName);//zjh

   fout << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear()" << std::endl;
   fout << "{" << std::endl;
   fout << "}" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// who the hell makes such strange Debug flags that even use "global pointers"..

Bool_t TMVA::MethodANNBase::Debug() const
{
   return fgDEBUG;
}
