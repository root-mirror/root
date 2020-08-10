// @(#)root/tmva/tmva/dnn:$Id$
// Author: Vladimir Ilievski, Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing MethodDL with DNN for various optimizers                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ravi Kiran S           <sravikiran0606@gmail.com>  - CERN, Switzerland    *
 *      Vladimir Ilievski      <ilievski.vladimir@live.com>  - CERN, Switzerland  *
 *                                                                                *
 * Copyright (c) 2005-2018:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef TMVA_TEST_DNN_TEST_METHOD_DL_OPTIMIZATION_H
#define TMVA_TEST_DNN_TEST_METHOD_DL_OPTIMIZATION_H

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"

#include "TMVA/MethodDL.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Config.h"

#include <iostream>

/** Testing the entire pipeline of the Method DL, when only a Multilayer Percepton
 *  is constructed with various optimizers. */
//______________________________________________________________________________
void testMethodDL_DNN(TString architectureStr, TString optimizerStr)
{

   ROOT::EnableImplicitMT(1);
   TMVA::Config::Instance();

   TFile *input(0);
   // TString fname = "tmva_class_example.root";
   TString fname = "http://root.cern/files/tmva_class_example.root";
   TFile::SetCacheFileDir(".");
   TString fopt = "CACHEREAD";
   input = TFile::Open(fname, fopt);

   // Register the training and test trees
   TTree *signalTree = (TTree *)input->Get("TreeS");
   TTree *background = (TTree *)input->Get("TreeB");

   TString outfileName = TString::Format("TMVA_DNN_%s_%s.root",optimizerStr.Data(),architectureStr.Data());
   TFile *outputFile = TFile::Open(outfileName, "RECREATE");

   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   dataloader->AddVariable("myvar1 := var1+var2", 'F');
   dataloader->AddVariable("myvar2 := var1-var2", "Expression 2", "", 'F');
   dataloader->AddVariable("var3", "Variable 3", "units", 'F');
   dataloader->AddVariable("var4", "Variable 4", "units", 'F');
   dataloader->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
   dataloader->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');

   // global event weights per tree
   Double_t signalWeight = 1.0;
   Double_t backgroundWeight = 1.0;

   // Add signal and background trees
   dataloader->AddSignalTree(signalTree, signalWeight);
   dataloader->AddBackgroundTree(background, backgroundWeight);

   dataloader->SetBackgroundWeightExpression("weight");

   TCut mycuts = "";
   TCut mycutb = "";

   dataloader->PrepareTrainingAndTestTree(
      mycuts, mycutb, "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V");

   // Input Layout
   TString inputLayoutString("InputLayout=1|1|4");

   // Batch Layout
   TString batchLayoutString("BatchLayout=256|1|4");

   // General layout.
   TString layoutString("Layout=RESHAPE|1|1|4|FLAT,DENSE|128|TANH,DENSE|128|TANH,DENSE|128|TANH,DENSE|1|LINEAR");

   // Training strategies.
   TString training0("LearningRate=1e-2,Optimizer=" + optimizerStr +
                     ",Momentum=0.9,Repetitions=1,"
                     "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                     "WeightDecay=1e-4,Regularization=L2,"
                     "DropConfig=0.0+0.5+0.5+0.5,MaxEpochs=100");
   // Having only one phase of training to reduce the test time.
   // TString training1("LearningRate=1e-2,Optimizer=" + optimizerStr +
   //                   ",Momentum=0.9,Repetitions=1,"
   //                   "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
   //                   "WeightDecay=1e-4,Regularization=L2,"
   //                   "DropConfig=0.0+0.0+0.0+0.0,Multithreading=True");
   // TString training2("LearningRate=1e-3,Optimizer=" + optimizerStr +
   //                   ",Momentum=0.9,Repetitions=1,"
   //                   "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
   //                   "WeightDecay=1e-4,Regularization=L2,"
   //                   "DropConfig=0.0+0.0+0.0+0.0,Multithreading=True");
   TString trainingStrategyString("TrainingStrategy=");
   trainingStrategyString += training0;

   // General Options.
   TString dnnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:"
                      "WeightInitialization=XAVIERUNIFORM");

   // Concatenate all option strings
   dnnOptions.Append(":");
   dnnOptions.Append(inputLayoutString);

   dnnOptions.Append(":");
   dnnOptions.Append(batchLayoutString);

   dnnOptions.Append(":");
   dnnOptions.Append(layoutString);

   dnnOptions.Append(":");
   dnnOptions.Append(trainingStrategyString);

   dnnOptions.Append(":Architecture=");
   dnnOptions.Append(architectureStr);

   // create factory
   TMVA::Factory *factory =
      new TMVA::Factory("TMVAClassification", outputFile,
                        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

   TString methodTitle = TString::Format("DL_%s_%s",optimizerStr.Data(),architectureStr.Data());
   factory->BookMethod(dataloader, TMVA::Types::kDL, methodTitle, dnnOptions);

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
}

#endif
