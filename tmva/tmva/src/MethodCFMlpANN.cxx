// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate Data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::MethodCFMlpANN                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Xavier Prudent  <prudent@lapp.in2p3.fr>  - LAPP, France                   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodCFMlpANN
\ingroup TMVA

Interface to Clermond-Ferrand artificial neural network


The CFMlpANN belong to the class of Multilayer Perceptrons (MLP), which are
feed-forward networks according to the following propagation schema:

\image html tmva_mlp.png Schema for artificial neural network.

The input layer contains as many neurons as input variables used in the MVA.
The output layer contains two neurons for the signal and background
event classes. In between the input and output layers are a variable number
of <i>k</i> hidden layers with arbitrary numbers of neurons. (While the
structure of the input and output layers is determined by the problem, the
hidden layers can be configured by the user through the option string
of the method booking.)

As indicated in the sketch, all neuron inputs to a layer are linear
combinations of the neuron output of the previous layer. The transfer
from input to output within a neuron is performed by means of an "activation
function". In general, the activation function of a neuron can be
zero (deactivated), one (linear), or non-linear. The above example uses
a sigmoid activation function. The transfer function of the output layer
is usually linear. As a consequence: an ANN without hidden layer should
give identical discrimination power as a linear discriminant analysis (Fisher).
In case of one hidden layer, the ANN computes a linear combination of
sigmoid.

The learning method used by the CFMlpANN is only stochastic.
*/


#include "TMVA/MethodCFMlpANN.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCFMlpANN_def.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TMatrix.h"
#include "TMath.h"

#include <cstdlib>
#include <iostream>
#include <string>



REGISTER_METHOD(CFMlpANN)

using std::stringstream;
using std::make_pair;
using std::atoi;

ClassImp(TMVA::MethodCFMlpANN);



////////////////////////////////////////////////////////////////////////////////
/// standard constructor
///
/// option string: "n_training_cycles:n_hidden_layers"
///
/// default is:  n_training_cycles = 5000, n_layers = 4
///
/// * note that the number of hidden layers in the NN is:
///   n_hidden_layers = n_layers - 2
///
/// * since there is one input and one output layer. The number of
///   nodes (neurons) is predefined to be:
///
///     n_nodes[i] = nvars + 1 - i (where i=1..n_layers)
///
///   with nvars being the number of variables used in the NN.
///
/// Hence, the default case is:
///
///     n_neurons(layer 1 (input)) : nvars
///     n_neurons(layer 2 (hidden)): nvars-1
///     n_neurons(layer 3 (hidden)): nvars-1
///     n_neurons(layer 4 (out))   : 2
///
/// This artificial neural network usually needs a relatively large
/// number of cycles to converge (8000 and more). Overtraining can
/// be efficiently tested by comparing the signal and background
/// output of the NN for the events that were used for training and
/// an independent data sample (with equal properties). If the separation
/// performance is significantly better for the training sample, the
/// NN interprets statistical effects, and is hence overtrained. In
/// this case, the number of cycles should be reduced, or the size
/// of the training sample increased.

TMVA::MethodCFMlpANN::MethodCFMlpANN( const TString& jobName,
                                      const TString& methodTitle,
                                      DataSetInfo& theData,
                                      const TString& theOption  ) :
   TMVA::MethodBase( jobName, Types::kCFMlpANN, methodTitle, theData, theOption),
   fData(0),
   fClass(0),
   fNlayers(0),
   fNcycles(0),
   fNodes(0),
   fYNN(0),
   MethodCFMlpANN_nsel(0)
{
   MethodCFMlpANN_Utils::SetLogger(&Log());
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodCFMlpANN::MethodCFMlpANN( DataSetInfo& theData,
                                      const TString& theWeightFile):
   TMVA::MethodBase( Types::kCFMlpANN, theData, theWeightFile),
   fData(0),
   fClass(0),
   fNlayers(0),
   fNcycles(0),
   fNodes(0),
   fYNN(0),
   MethodCFMlpANN_nsel(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// CFMlpANN can handle classification with 2 classes

Bool_t TMVA::MethodCFMlpANN::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string
/// know options: NCycles=xx              :the number of training cycles
///               HiddenLayser="N-1,N-2"  :the specification of the hidden layers

void TMVA::MethodCFMlpANN::DeclareOptions()
{
   DeclareOptionRef( fNcycles  =3000,      "NCycles",      "Number of training cycles" );
   DeclareOptionRef( fLayerSpec="N,N-1",   "HiddenLayers", "Specification of hidden layer architecture" );
}

////////////////////////////////////////////////////////////////////////////////
/// decode the options in the option string

void TMVA::MethodCFMlpANN::ProcessOptions()
{
   fNodes = new Int_t[20]; // number of nodes per layer (maximum 20 layers)
   fNlayers = 2;
   Int_t currentHiddenLayer = 1;
   TString layerSpec(fLayerSpec);
   while(layerSpec.Length()>0) {
      TString sToAdd = "";
      if (layerSpec.First(',')<0) {
         sToAdd = layerSpec;
         layerSpec = "";
      }
      else {
         sToAdd = layerSpec(0,layerSpec.First(','));
         layerSpec = layerSpec(layerSpec.First(',')+1,layerSpec.Length());
      }
      Int_t nNodes = 0;
      if (sToAdd.BeginsWith("N") || sToAdd.BeginsWith("n")) { sToAdd.Remove(0,1); nNodes = GetNvar(); }
      nNodes += atoi(sToAdd);
      fNodes[currentHiddenLayer++] = nNodes;
      fNlayers++;
   }
   fNodes[0]          = GetNvar(); // number of input nodes
   fNodes[fNlayers-1] = 2;         // number of output nodes

   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kFATAL << "Mechanism to ignore events with negative weights in training not yet available for method: "
            << GetMethodTypeName()
            << " --> please remove \"IgnoreNegWeightsInTraining\" option from booking string."
            << Endl;
   }

   Log() << kINFO << "Use configuration (nodes per layer): in=";
   for (Int_t i=0; i<fNlayers-1; i++) Log() << kINFO << fNodes[i] << ":";
   Log() << kINFO << fNodes[fNlayers-1] << "=out" << Endl;

   // some info
   Log() << "Use " << fNcycles << " training cycles" << Endl;

   Int_t nEvtTrain = Data()->GetNTrainingEvents();

   // note that one variable is type
   if (nEvtTrain>0) {

      // Data LUT
      fData  = new TMatrix( nEvtTrain, GetNvar() );
      fClass = new std::vector<Int_t>( nEvtTrain );

      // ---- fill LUTs

      UInt_t ivar;
      for (Int_t ievt=0; ievt<nEvtTrain; ievt++) {
         const Event * ev = GetEvent(ievt);

         // identify signal and background events
         (*fClass)[ievt] = DataInfo().IsSignal(ev) ? 1 : 2;

         // use normalized input Data
         for (ivar=0; ivar<GetNvar(); ivar++) {
            (*fData)( ievt, ivar ) = ev->GetValue(ivar);
         }
      }

      //Log() << kVERBOSE << Data()->GetNEvtSigTrain() << " Signal and "
      //        << Data()->GetNEvtBkgdTrain() << " background" << " events in trainingTree" << Endl;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// default initialisation called by all constructors

void TMVA::MethodCFMlpANN::Init( void )
{
   // CFMlpANN prefers normalised input variables
   SetNormalised( kTRUE );

   // initialize dimensions
   MethodCFMlpANN_nsel = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodCFMlpANN::~MethodCFMlpANN( void )
{
   delete fData;
   delete fClass;
   delete[] fNodes;

   if (fYNN!=0) {
      for (Int_t i=0; i<fNlayers; i++) delete[] fYNN[i];
      delete[] fYNN;
      fYNN=0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// training of the Clement-Ferrand NN classifier

void TMVA::MethodCFMlpANN::Train( void )
{
   Double_t dumDat(0);
   Int_t ntrain(Data()->GetNTrainingEvents());
   Int_t ntest(0);
   Int_t nvar(GetNvar());
   Int_t nlayers(fNlayers);
   Int_t *nodes = new Int_t[nlayers];
   Int_t ncycles(fNcycles);

   for (Int_t i=0; i<nlayers; i++) nodes[i] = fNodes[i]; // full copy of class member

   if (fYNN != 0) {
      for (Int_t i=0; i<fNlayers; i++) delete[] fYNN[i];
      delete[] fYNN;
      fYNN = 0;
   }
   fYNN = new Double_t*[nlayers];
   for (Int_t layer=0; layer<nlayers; layer++)
      fYNN[layer] = new Double_t[fNodes[layer]];

   // please check
#ifndef R__WIN32
   Train_nn( &dumDat, &dumDat, &ntrain, &ntest, &nvar, &nlayers, nodes, &ncycles );
#else
   Log() << kWARNING << "<Train> sorry CFMlpANN does not run on Windows" << Endl;
#endif

   delete [] nodes;

   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// returns CFMlpANN output (normalised within [0,1])

Double_t TMVA::MethodCFMlpANN::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   Bool_t isOK = kTRUE;

   const Event* ev = GetEvent();

   // copy of input variables
   std::vector<Double_t> inputVec( GetNvar() );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) inputVec[ivar] = ev->GetValue(ivar);

   Double_t myMVA = EvalANN( inputVec, isOK );
   if (!isOK) Log() << kFATAL << "EvalANN returns (!isOK) for event " << Endl;

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return myMVA;
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates NN value as function of input variables

Double_t TMVA::MethodCFMlpANN::EvalANN( std::vector<Double_t>& inVar, Bool_t& isOK )
{
   // hardcopy of input variables (necessary because they are update later)
   Double_t* xeev = new Double_t[GetNvar()];
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) xeev[ivar] = inVar[ivar];

   // ---- now apply the weights: get NN output
   isOK = kTRUE;
   for (UInt_t jvar=0; jvar<GetNvar(); jvar++) {

      if (fVarn_1.xmax[jvar] < xeev[jvar]) xeev[jvar] = fVarn_1.xmax[jvar];
      if (fVarn_1.xmin[jvar] > xeev[jvar]) xeev[jvar] = fVarn_1.xmin[jvar];
      if (fVarn_1.xmax[jvar] == fVarn_1.xmin[jvar]) {
         isOK = kFALSE;
         xeev[jvar] = 0;
      }
      else {
         xeev[jvar] = xeev[jvar] - ((fVarn_1.xmax[jvar] + fVarn_1.xmin[jvar])/2);
         xeev[jvar] = xeev[jvar] / ((fVarn_1.xmax[jvar] - fVarn_1.xmin[jvar])/2);
      }
   }

   NN_ava( xeev );

   Double_t retval = 0.5*(1.0 + fYNN[fParam_1.layerm-1][0]);

   delete [] xeev;

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// auxiliary functions

void  TMVA::MethodCFMlpANN::NN_ava( Double_t* xeev )
{
   for (Int_t ivar=0; ivar<fNeur_1.neuron[0]; ivar++) fYNN[0][ivar] = xeev[ivar];

   for (Int_t layer=1; layer<fParam_1.layerm; layer++) {
      for (Int_t j=1; j<=fNeur_1.neuron[layer]; j++) {

         Double_t x = Ww_ref(fNeur_1.ww, layer+1,j); // init with the bias layer

         for (Int_t k=1; k<=fNeur_1.neuron[layer-1]; k++) { // neurons of originating layer
            x += fYNN[layer-1][k-1]*W_ref(fNeur_1.w, layer+1, j, k);
         }
         fYNN[layer][j-1] = NN_fonc( layer, x );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// activation function

Double_t TMVA::MethodCFMlpANN::NN_fonc( Int_t i, Double_t u ) const
{
   Double_t f(0);

   if      (u/fDel_1.temp[i] >  170) f = +1;
   else if (u/fDel_1.temp[i] < -170) f = -1;
   else {
      Double_t yy = TMath::Exp(-u/fDel_1.temp[i]);
      f  = (1 - yy)/(1 + yy);
   }

   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// read back the weight from the training from file (stream)

void TMVA::MethodCFMlpANN::ReadWeightsFromStream( std::istream & istr )
{
   TString var;

   // read number of variables and classes
   UInt_t nva(0), lclass(0);
   istr >> nva >> lclass;

   if (GetNvar() != nva) // wrong file
      Log() << kFATAL << "<ReadWeightsFromFile> mismatch in number of variables" << Endl;

   // number of output classes must be 2
   if (lclass != 2) // wrong file
      Log() << kFATAL << "<ReadWeightsFromFile> mismatch in number of classes" << Endl;

   // check that we are not at the end of the file
   if (istr.eof( ))
      Log() << kFATAL << "<ReadWeightsFromStream> reached EOF prematurely " << Endl;

   // read extrema of input variables
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++)
      istr >> fVarn_1.xmax[ivar] >> fVarn_1.xmin[ivar];

   // read number of layers (sum of: input + output + hidden)
   istr >> fParam_1.layerm;

   if (fYNN != 0) {
      for (Int_t i=0; i<fNlayers; i++) delete[] fYNN[i];
      delete[] fYNN;
      fYNN = 0;
   }
   fYNN = new Double_t*[fParam_1.layerm];
   for (Int_t layer=0; layer<fParam_1.layerm; layer++) {
      // read number of neurons for each layer
      // coverity[tainted_data_argument]
      istr >> fNeur_1.neuron[layer];
      fYNN[layer] = new Double_t[fNeur_1.neuron[layer]];
   }

   // to read dummy lines
   const Int_t nchar( 100 );
   char* dumchar = new char[nchar];

   // read weights
   for (Int_t layer=1; layer<=fParam_1.layerm-1; layer++) {

      Int_t nq = fNeur_1.neuron[layer]/10;
      Int_t nr = fNeur_1.neuron[layer] - nq*10;

      Int_t kk(0);
      if (nr==0) kk = nq;
      else       kk = nq+1;

      for (Int_t k=1; k<=kk; k++) {
         Int_t jmin = 10*k - 9;
         Int_t jmax = 10*k;
         if (fNeur_1.neuron[layer]<jmax) jmax = fNeur_1.neuron[layer];
         for (Int_t j=jmin; j<=jmax; j++) {
            istr >> Ww_ref(fNeur_1.ww, layer+1, j);
         }
         for (Int_t i=1; i<=fNeur_1.neuron[layer-1]; i++) {
            for (Int_t j=jmin; j<=jmax; j++) {
               istr >> W_ref(fNeur_1.w, layer+1, j, i);
            }
         }
         // skip two empty lines
         istr.getline( dumchar, nchar );
      }
   }

   for (Int_t layer=0; layer<fParam_1.layerm; layer++) {

      // skip 2 empty lines
      istr.getline( dumchar, nchar );
      istr.getline( dumchar, nchar );

      istr >> fDel_1.temp[layer];
   }

   // sanity check
   if ((Int_t)GetNvar() != fNeur_1.neuron[0]) {
      Log() << kFATAL << "<ReadWeightsFromFile> mismatch in zeroth layer:"
            << GetNvar() << " " << fNeur_1.neuron[0] << Endl;
   }

   fNlayers = fParam_1.layerm;
   delete[] dumchar;
}

////////////////////////////////////////////////////////////////////////////////
/// data interface function

Int_t TMVA::MethodCFMlpANN::DataInterface( Double_t* /*tout2*/, Double_t*  /*tin2*/,
                                           Int_t* /* icode*/, Int_t*  /*flag*/,
                                           Int_t*  /*nalire*/, Int_t* nvar,
                                           Double_t* xpg, Int_t* iclass, Int_t* ikend )
{
   // icode and ikend are dummies needed to match f2c mlpl3 functions
   *ikend = 0;


   // sanity checks
   if (0 == xpg) {
      Log() << kFATAL << "ERROR in MethodCFMlpANN_DataInterface zero pointer xpg" << Endl;
   }
   if (*nvar != (Int_t)this->GetNvar()) {
      Log() << kFATAL << "ERROR in MethodCFMlpANN_DataInterface mismatch in num of variables: "
            << *nvar << " " << this->GetNvar() << Endl;
   }

   // fill variables
   *iclass = (int)this->GetClass( MethodCFMlpANN_nsel );
   for (UInt_t ivar=0; ivar<this->GetNvar(); ivar++)
      xpg[ivar] = (double)this->GetData( MethodCFMlpANN_nsel, ivar );

   ++MethodCFMlpANN_nsel;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// write weights to xml file

void TMVA::MethodCFMlpANN::AddWeightsXMLTo( void* parent ) const
{
   void *wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr(wght,"NVars",fParam_1.nvar);
   gTools().AddAttr(wght,"NClasses",fParam_1.lclass);
   gTools().AddAttr(wght,"NLayers",fParam_1.layerm);
   void* minmaxnode = gTools().AddChild(wght, "VarMinMax");
   stringstream s;
   s.precision( 16 );
   for (Int_t ivar=0; ivar<fParam_1.nvar; ivar++)
      s << std::scientific << fVarn_1.xmin[ivar] <<  " " << fVarn_1.xmax[ivar] <<  " ";
   gTools().AddRawLine( minmaxnode, s.str().c_str() );
   void* neurons = gTools().AddChild(wght, "NNeurons");
   stringstream n;
   n.precision( 16 );
   for (Int_t layer=0; layer<fParam_1.layerm; layer++)
      n << std::scientific << fNeur_1.neuron[layer] << " ";
   gTools().AddRawLine( neurons, n.str().c_str() );
   for (Int_t layer=1; layer<fParam_1.layerm; layer++) {
      void* layernode = gTools().AddChild(wght, "Layer"+gTools().StringFromInt(layer));
      gTools().AddAttr(layernode,"NNeurons",fNeur_1.neuron[layer]);
      void* neuronnode=NULL;
      for (Int_t neuron=0; neuron<fNeur_1.neuron[layer]; neuron++) {
         neuronnode = gTools().AddChild(layernode,"Neuron"+gTools().StringFromInt(neuron));
         stringstream weights;
         weights.precision( 16 );
         weights << std::scientific << Ww_ref(fNeur_1.ww, layer+1, neuron+1);
         for (Int_t i=0; i<fNeur_1.neuron[layer-1]; i++) {
            weights << " " << std::scientific << W_ref(fNeur_1.w, layer+1, neuron+1, i+1);
         }
         gTools().AddRawLine( neuronnode, weights.str().c_str() );
      }
   }
   void* tempnode = gTools().AddChild(wght, "LayerTemp");
   stringstream temp;
   temp.precision( 16 );
   for (Int_t layer=0; layer<fParam_1.layerm; layer++) {
      temp << std::scientific << fDel_1.temp[layer] << " ";
   }
   gTools().AddRawLine(tempnode, temp.str().c_str() );
}
////////////////////////////////////////////////////////////////////////////////
/// read weights from xml file

void TMVA::MethodCFMlpANN::ReadWeightsFromXML( void* wghtnode )
{
   gTools().ReadAttr( wghtnode, "NLayers",fParam_1.layerm );
   void* minmaxnode = gTools().GetChild(wghtnode);
   const char* minmaxcontent = gTools().GetContent(minmaxnode);
   stringstream content(minmaxcontent);
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++)
      content >> fVarn_1.xmin[ivar] >> fVarn_1.xmax[ivar];
   if (fYNN != 0) {
      for (Int_t i=0; i<fNlayers; i++) delete[] fYNN[i];
      delete[] fYNN;
      fYNN = 0;
   }
   fYNN = new Double_t*[fParam_1.layerm];
   void *layernode=gTools().GetNextChild(minmaxnode);
   const char* neuronscontent = gTools().GetContent(layernode);
   stringstream ncontent(neuronscontent);
   for (Int_t layer=0; layer<fParam_1.layerm; layer++) {
      // read number of neurons for each layer;
      // coverity[tainted_data_argument]
      ncontent >> fNeur_1.neuron[layer];
      fYNN[layer] = new Double_t[fNeur_1.neuron[layer]];
   }
   for (Int_t layer=1; layer<fParam_1.layerm; layer++) {
      layernode=gTools().GetNextChild(layernode);
      void* neuronnode=NULL;
      neuronnode = gTools().GetChild(layernode);
      for (Int_t neuron=0; neuron<fNeur_1.neuron[layer]; neuron++) {
         const char* neuronweights = gTools().GetContent(neuronnode);
         stringstream weights(neuronweights);
         weights >> Ww_ref(fNeur_1.ww, layer+1, neuron+1);
         for (Int_t i=0; i<fNeur_1.neuron[layer-1]; i++) {
            weights >> W_ref(fNeur_1.w, layer+1, neuron+1, i+1);
         }
         neuronnode=gTools().GetNextChild(neuronnode);
      }
   }
   void* tempnode=gTools().GetNextChild(layernode);
   const char* temp = gTools().GetContent(tempnode);
   stringstream t(temp);
   for (Int_t layer=0; layer<fParam_1.layerm; layer++) {
      t >> fDel_1.temp[layer];
   }
   fNlayers = fParam_1.layerm;
}

////////////////////////////////////////////////////////////////////////////////
/// write the weights of the neural net

void TMVA::MethodCFMlpANN::PrintWeights( std::ostream & o ) const
{
   // write number of variables and classes
   o << "Number of vars " << fParam_1.nvar << std::endl;
   o << "Output nodes   " << fParam_1.lclass << std::endl;

   // write extrema of input variables
   for (Int_t ivar=0; ivar<fParam_1.nvar; ivar++)
      o << "Var " << ivar << " [" << fVarn_1.xmin[ivar] << " - " << fVarn_1.xmax[ivar] << "]" << std::endl;

   // write number of layers (sum of: input + output + hidden)
   o << "Number of layers " << fParam_1.layerm << std::endl;

   o << "Nodes per layer ";
   for (Int_t layer=0; layer<fParam_1.layerm; layer++)
      // write number of neurons for each layer
      o << fNeur_1.neuron[layer] << "     ";
   o << std::endl;

   // write weights
   for (Int_t layer=1; layer<=fParam_1.layerm-1; layer++) {

      Int_t nq = fNeur_1.neuron[layer]/10;
      Int_t nr = fNeur_1.neuron[layer] - nq*10;

      Int_t kk(0);
      if (nr==0) kk = nq;
      else       kk = nq+1;

      for (Int_t k=1; k<=kk; k++) {
         Int_t jmin = 10*k - 9;
         Int_t jmax = 10*k;
         Int_t i, j;
         if (fNeur_1.neuron[layer]<jmax) jmax = fNeur_1.neuron[layer];
         for (j=jmin; j<=jmax; j++) {

            //o << fNeur_1.ww[j*max_nLayers_ + layer - 6] << "   ";
            o << Ww_ref(fNeur_1.ww, layer+1, j) << "   ";

         }
         o << std::endl;
         //for (i=1; i<=fNeur_1.neuron[layer-1]; i++) {
         for (i=1; i<=fNeur_1.neuron[layer-1]; i++) {
            for (j=jmin; j<=jmax; j++) {
               //               o << fNeur_1.w[(i*max_nNodes_ + j)*max_nLayers_ + layer - 186] << "   ";
               o << W_ref(fNeur_1.w, layer+1, j, i) << "   ";
            }
            o << std::endl;
         }

         // skip two empty lines
         o << std::endl;
      }
   }
   for (Int_t layer=0; layer<fParam_1.layerm; layer++) {
      o << "Del.temp in layer " << layer << " :  " << fDel_1.temp[layer] << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodCFMlpANN::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   // write specific classifier response
   fout << "   // not implemented for class: \"" << className << "\"" << std::endl;
   fout << "};" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response for header

void TMVA::MethodCFMlpANN::MakeClassSpecificHeader( std::ostream& , const TString&  ) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodCFMlpANN::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
}
