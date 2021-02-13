// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodDNN                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      NeuralNetwork                                                             *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer      <peter.speckmayer@gmx.at>  - CERN, Switzerland      *
 *      Simon Pfreundschuh    <s.pfreundschuh@gmail.com> - CERN, Switzerland      *
 *                                                                                *
 * Copyright (c) 2005-2015:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

//#pragma once

#ifndef ROOT_TMVA_MethodDNN
#define ROOT_TMVA_MethodDNN

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodDNN                                                             //
//                                                                      //
// Neural Network implementation                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TMVA/MethodBase.h"
#include "TMVA/NeuralNet.h"

#include "TMVA/Tools.h"

#include "TMVA/DNN/Net.h"
#include "TMVA/DNN/Minimizers.h"
#include "TMVA/DNN/Architectures/Reference.h"

#ifdef R__HAS_TMVACPU
#define DNNCPU
#endif
#ifdef R__HAS_TMVAGPU
//#define DNNCUDA
#endif

#ifdef DNNCPU
#include "TMVA/DNN/Architectures/Cpu.h"
#endif

#ifdef DNNCUDA
#include "TMVA/DNN/Architectures/Cuda.h"
#endif

namespace TMVA {

class MethodDNN : public MethodBase
{
   friend struct TestMethodDNNValidationSize;

   using Architecture_t = DNN::TReference<Float_t>;
   using Net_t          = DNN::TNet<Architecture_t>;
   using Matrix_t       = typename Architecture_t::Matrix_t;
   using Scalar_t       = typename Architecture_t::Scalar_t;

private:
   using LayoutVector_t   = std::vector<std::pair<int, DNN::EActivationFunction>>;
   using KeyValueVector_t = std::vector<std::map<TString, TString>>;

   struct TTrainingSettings
   {
       size_t                batchSize;
       size_t                testInterval;
       size_t                convergenceSteps;
       DNN::ERegularization  regularization;
       Double_t              learningRate;
       Double_t              momentum;
       Double_t              weightDecay;
       std::vector<Double_t> dropoutProbabilities;
       bool                  multithreading;
   };

   // the option handling methods
   void DeclareOptions();
   void ProcessOptions();

   UInt_t GetNumValidationSamples();

   // general helper functions
   void     Init();

   Net_t                fNet;
   DNN::EInitialization fWeightInitialization;
   DNN::EOutputFunction fOutputFunction;

   TString                        fLayoutString;
   TString                        fErrorStrategy;
   TString                        fTrainingStrategyString;
   TString                        fWeightInitializationString;
   TString                        fArchitectureString;
   TString                        fValidationSize;
   LayoutVector_t                 fLayout;
   std::vector<TTrainingSettings> fTrainingSettings;
   bool                           fResume;

   KeyValueVector_t fSettings;

   ClassDef(MethodDNN,0); // neural network

   static inline void WriteMatrixXML(void *parent, const char *name,
                                     const TMatrixT<Double_t> &X);
   static inline void ReadMatrixXML(void *xml, const char *name,
                                    TMatrixT<Double_t> &X);
protected:

   void MakeClassSpecific( std::ostream&, const TString& ) const;
   void GetHelpMessage() const;

public:

   // Standard Constructors
   MethodDNN(const TString& jobName,
             const TString&  methodTitle,
             DataSetInfo& theData,
             const TString& theOption);
   MethodDNN(DataSetInfo& theData,
             const TString& theWeightFile);
   virtual ~MethodDNN();

   virtual Bool_t HasAnalysisType(Types::EAnalysisType type,
                                  UInt_t numberClasses,
                                  UInt_t numberTargets );
   LayoutVector_t   ParseLayoutString(TString layerSpec);
   KeyValueVector_t ParseKeyValueString(TString parseString,
                                      TString blockDelim,
                                      TString tokenDelim);
   void Train();
   void TrainGpu();
   void TrainCpu();

   virtual Double_t GetMvaValue( Double_t* err=0, Double_t* errUpper=0 );
   virtual const std::vector<Float_t>& GetRegressionValues();
   virtual const std::vector<Float_t>& GetMulticlassValues();

   using MethodBase::ReadWeightsFromStream;

   // write weights to stream
   void AddWeightsXMLTo     ( void* parent ) const;

   // read weights from stream
   void ReadWeightsFromStream( std::istream & i );
   void ReadWeightsFromXML   ( void* wghtnode );

   // ranking of input variables
   const Ranking* CreateRanking();

};

inline void MethodDNN::WriteMatrixXML(void *parent,
                                      const char *name,
                                      const TMatrixT<Double_t> &X)
{
   std::stringstream matrixStringStream("");
   matrixStringStream.precision( 16 );

   for (size_t i = 0; i < (size_t) X.GetNrows(); i++)
   {
      for (size_t j = 0; j < (size_t) X.GetNcols(); j++)
      {
         matrixStringStream << std::scientific << X(i,j) << " ";
      }
   }
   std::string s = matrixStringStream.str();
   void* matxml = gTools().xmlengine().NewChild(parent, 0, name);
   gTools().xmlengine().NewAttr(matxml, 0, "rows",
                                gTools().StringFromInt((int)X.GetNrows()));
   gTools().xmlengine().NewAttr(matxml, 0, "cols",
                                gTools().StringFromInt((int)X.GetNcols()));
   gTools().xmlengine().AddRawLine (matxml, s.c_str());
}

inline void MethodDNN::ReadMatrixXML(void *xml,
                                     const char *name,
                                     TMatrixT<Double_t> &X)
{
   void *matrixXML = gTools().GetChild(xml, name);
   size_t rows, cols;
   gTools().ReadAttr(matrixXML, "rows", rows);
   gTools().ReadAttr(matrixXML, "cols", cols);

   const char * matrixString = gTools().xmlengine().GetNodeContent(matrixXML);
   std::stringstream matrixStringStream(matrixString);

   for (size_t i = 0; i < rows; i++)
   {
      for (size_t j = 0; j < cols; j++)
      {
         matrixStringStream >> X(i,j);
      }
   }
}
} // namespace TMVA

#endif
