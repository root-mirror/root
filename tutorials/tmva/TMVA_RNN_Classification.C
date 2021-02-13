/// \file
/// \ingroup tutorial_tmva
/// \notebook
///  TMVA Classification Example Using a Recurrent Neural Network
///
/// This is an example of using a RNN in TMVA. We do classification using a toy time dependent data set
/// that is generated when running this example macro
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta
/***

    # TMVA Classification Example Using a Recurrent Neural Network

    This is an example of using a RNN in TMVA.
    We do the classification using a toy data set containing a time series of data sample ntimes
    and with dimension ndim that is generated when running the provided function `MakeTimeData (nevents, ntime, ndim)`


**/

#include<TROOT.h>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Config.h"
#include "TMVA/MethodDL.h"


#include "TFile.h"
#include "TTree.h"

///  Helper function to generate the time data set
///  make some time data but not of fixed length.
///  use a poisson with mu = 5 and troncated at 10
///
void MakeTimeData(int n, int ntime, int ndim )
{

   // const int ntime = 10;
   // const int ndim = 30; // number of dim/time
   TString fname = TString::Format("time_data_t%d_d%d.root", ntime, ndim);
   std::vector<TH1 *> v1(ntime);
   std::vector<TH1 *> v2(ntime);
   int i = 0;
   for (int i = 0; i < ntime; ++i) {
      v1[i] = new TH1D(TString::Format("h1_%d", i), "h1", ndim, 0, 10);
      v2[i] = new TH1D(TString::Format("h2_%d", i), "h2", ndim, 0, 10);
   }

   auto f1 = new TF1("f1", "gaus");
   auto f2 = new TF1("f2", "gaus");

   TTree sgn("sgn", "sgn");
   TTree bkg("bkg", "bkg");
   TFile f(fname, "RECREATE");

   std::vector<std::vector<float>> x1(ntime);
   std::vector<std::vector<float>> x2(ntime);

   for (int i = 0; i < ntime; ++i) {
      x1[i] = std::vector<float>(ndim);
      x2[i] = std::vector<float>(ndim);
   }

   for (auto i = 0; i < ntime; i++) {
      bkg.Branch(Form("vars_time%d", i), "std::vector<float>", &x1[i]);
      sgn.Branch(Form("vars_time%d", i), "std::vector<float>", &x2[i]);
   }

   sgn.SetDirectory(&f);
   bkg.SetDirectory(&f);
   gRandom->SetSeed(0);

   std::vector<double> mean1(ntime);
   std::vector<double> mean2(ntime);
   std::vector<double> sigma1(ntime);
   std::vector<double> sigma2(ntime);
   for (int j = 0; j < ntime; ++j) {
      mean1[j] = 5. + 0.2 * sin(TMath::Pi() * j / double(ntime));
      mean2[j] = 5. + 0.2 * cos(TMath::Pi() * j / double(ntime));
      sigma1[j] = 4 + 0.3 * sin(TMath::Pi() * j / double(ntime));
      sigma2[j] = 4 + 0.3 * cos(TMath::Pi() * j / double(ntime));
   }
   for (int i = 0; i < n; ++i) {

      if (i % 1000 == 0)
         std::cout << "Generating  event ... " << i << std::endl;

      for (int j = 0; j < ntime; ++j) {
         auto h1 = v1[j];
         auto h2 = v2[j];
         h1->Reset();
         h2->Reset();

         f1->SetParameters(1, mean1[j], sigma1[j]);
         f2->SetParameters(1, mean2[j], sigma2[j]);

         h1->FillRandom("f1", 1000);
         h2->FillRandom("f2", 1000);

         for (int k = 0; k < ndim; ++k) {
            // std::cout << j*10+k << "   ";
            x1[j][k] = h1->GetBinContent(k + 1) + gRandom->Gaus(0, 10);
            x2[j][k] = h2->GetBinContent(k + 1) + gRandom->Gaus(0, 10);
         }
      }
      // std::cout << std::endl;
      sgn.Fill();
      bkg.Fill();

      if (n == 1) {
         auto c1 = new TCanvas();
         c1->Divide(ntime, 2);
         for (int j = 0; j < ntime; ++j) {
            c1->cd(j + 1);
            v1[j]->Draw();
         }
         for (int j = 0; j < ntime; ++j) {
            c1->cd(ntime + j + 1);
            v2[j]->Draw();
         }
         gPad->Update();
      }
   }
   if (n > 1) {
      sgn.Write();
      bkg.Write();
      sgn.Print();
      bkg.Print();
      f.Close();
   }
}
/// macro for performing a classification using a Recurrent Neural Network
/// @param use_type
///    use_type = 0    use Simple RNN network
///    use_type = 1    use LSTM network
///    use_type = 2    use GRU
///    use_type = 3    build 3 different networks with RNN, LSTM and GRU

void TMVA_RNN_Classification(int use_type = 1)
{

   const int ninput = 30;
   const int ntime = 10;
   const int batchSize = 100;
   const int maxepochs = 20;

   int nTotEvts = 10000; // total events to be generated for signal or background

   bool useKeras = true;


   bool useTMVA_RNN = true;
   bool useTMVA_DNN = true;
   bool useTMVA_BDT = false;

   std::vector<std::string> rnn_types = {"RNN", "LSTM", "GRU"};
   std::vector<bool> use_rnn_type = {1, 1, 1};
   if (use_type >=0 && use_type < 3) {
      use_rnn_type = {0,0,0};
      use_rnn_type[use_type] = 1;
   }
   bool useGPU = true;   // use GPU for TMVA if available

#ifndef R__HAS_TMVAGPU
   useGPU = false;
#ifndef R__HAS_TMVACPU
   Warning("TMVA_RNN_Classification", "TMVA is not build with GPU or CPU multi-thread support. Cannot use TMVA Deep Learning for RNN");
   useTMVA_RNN = false;
#endif
#endif


   TString archString = (useGPU) ? "GPU" : "CPU";

   bool writeOutputFile = true;



   const char *rnn_type = "RNN";

#ifdef R__HAS_PYMVA
   TMVA::PyMethodBase::PyInitialize();
#else
   useKeras = false;
#endif

   int num_threads = 0;   // use by default all threads 
   // do enable MT running
   if (num_threads >= 0) {
      ROOT::EnableImplicitMT(num_threads);
      if (num_threads > 0) gSystem->Setenv("OMP_NUM_THREADS", TString::Format("%d",num_threads));
   }
   else
      gSystem->Setenv("OMP_NUM_THREADS", "1");

   TMVA::Config::Instance();

   std::cout << "Running with nthreads  = " << ROOT::GetThreadPoolSize() << std::endl;

   TString inputFileName = "time_data_t10_d30.root";

   bool fileExist = !gSystem->AccessPathName(inputFileName);

   // if file does not exists create it
   if (!fileExist) {
      MakeTimeData(nTotEvts,ntime, ninput);
   }


   auto inputFile = TFile::Open(inputFileName);
   if (!inputFile) {
      Error("TMVA_RNN_Classification", "Error opening input file %s - exit", inputFileName.Data());
      return;
   }


   std::cout << "--- RNNClassification  : Using input file: " << inputFile->GetName() << std::endl;

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName(TString::Format("data_RNN_%s.root", archString.Data()));
   TFile *outputFile = nullptr;
   if (writeOutputFile) outputFile = TFile::Open(outfileName, "RECREATE");

   /**
    ## Declare Factory

    Create the Factory class. Later you can choose the methods
    whose performance you'd like to investigate.

    The factory is the major TMVA object you have to interact with. Here is the list of parameters you need to
pass

    - The first argument is the base of the name of all the output
    weightfiles in the directory weight/ that will be created with the
    method parameters

    - The second argument is the output file for the training results

    - The third argument is a string option defining some general configuration for the TMVA session.
      For example all TMVA output can be suppressed by removing the "!" (not) in front of the "Silent" argument in
the option string

    **/

   // Creating the factory object
   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=None:!Correlations:"
                                              "AnalysisType=Classification:ModelPersistence");
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   TTree *signalTree = (TTree *)inputFile->Get("sgn");
   TTree *background = (TTree *)inputFile->Get("bkg");

   const int nvar = ninput * ntime;

   /// add variables - use new AddVariablesArray function
   for (auto i = 0; i < ntime; i++) {
      dataloader->AddVariablesArray(Form("vars_time%d", i), ninput);
   }

   dataloader->AddSignalTree(signalTree, 1.0);
   dataloader->AddBackgroundTree(background, 1.0);

   // check given input
   auto &datainfo = dataloader->GetDataSetInfo();
   auto vars = datainfo.GetListOfVariables();
   std::cout << "number of variables is " << vars.size() << std::endl;
   for (auto &v : vars)
      std::cout << v << ",";
   std::cout << std::endl;

   int nTrainSig = 0.8 * nTotEvts;
   int nTrainBkg = 0.8 *  nTotEvts;

   // build the string options for DataLoader::PrepareTrainingAndTestTree
   TString prepareOptions = TString::Format("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:SplitSeed=100:NormMode=NumEvents:!V:!CalcCorrelations", nTrainSig, nTrainBkg);

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = "";

   dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, prepareOptions);

   std::cout << "prepared DATA LOADER " << std::endl;

   /**
       ## Book TMVA  recurrent models

      Book the different types of recurrent models in TMVA  (SimpleRNN, LSTM or GRU)

 **/

   if (useTMVA_RNN) {

      for (int i = 0; i < 3; ++i) {

         if (!use_rnn_type[i])
            continue;

         const char *rnn_type = rnn_types[i].c_str();

         /// define the inputlayout string for RNN
         /// the input data should be organize as   following:
         //// input layout for RNN:    time x ndim

         TString inputLayoutString = TString::Format("InputLayout=%d|%d", ntime, ninput);

         /// Define RNN layer layout
         ///  it should be   LayerType (RNN or LSTM or GRU) |  number of units | number of inputs | time steps | remember output (typically no=0 | return full sequence
         TString rnnLayout = TString::Format("%s|10|%d|%d|0|1", rnn_type, ninput, ntime);

         /// add after RNN a reshape layer (needed top flatten the output) and a dense layer with 64 units and a last one
         /// Note the last layer is linear because  when using Crossentropy a Sigmoid is applied already
         TString layoutString = TString("Layout=") + rnnLayout + TString(",RESHAPE|FLAT,DENSE|64|TANH,LINEAR");

         /// Defining Training strategies. Different training strings can be concatenate. Use however only one
         TString trainingString1 = TString::Format("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                                             "ConvergenceSteps=5,BatchSize=%d,TestRepetitions=1,"
                                             "WeightDecay=1e-2,Regularization=None,MaxEpochs=%d,"
                                             "Optimizer=ADAM,DropConfig=0.0+0.+0.+0.",
                                             batchSize,maxepochs);

         TString trainingStrategyString("TrainingStrategy=");
         trainingStrategyString += trainingString1; // + "|" + trainingString2

         /// Define the full RNN Noption string adding the final options for all network
         TString rnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                            "WeightInitialization=XAVIERUNIFORM:ValidationSize=0.2:RandomSeed=1234");

         rnnOptions.Append(":");
         rnnOptions.Append(inputLayoutString);
         rnnOptions.Append(":");
         rnnOptions.Append(layoutString);
         rnnOptions.Append(":");
         rnnOptions.Append(trainingStrategyString);
         rnnOptions.Append(":");
         rnnOptions.Append(TString::Format("Architecture=%s", archString.Data()));

         TString rnnName = "TMVA_" + TString(rnn_type);
         factory->BookMethod(dataloader, TMVA::Types::kDL, rnnName, rnnOptions);

      }
   }

   /**
      ## Book TMVA  fully connected dense layer  models

   **/

   if (useTMVA_DNN) {
      // Method DL with Dense Layer
      TString inputLayoutString = TString::Format("InputLayout=1|1|%d", ntime * ninput);

      TString layoutString("Layout=DENSE|64|TANH,DENSE|TANH|64,DENSE|TANH|64,LINEAR");
      // Training strategies.
      TString trainingString1("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                        "ConvergenceSteps=10,BatchSize=256,TestRepetitions=1,"
                        "WeightDecay=1e-4,Regularization=None,MaxEpochs=20"
                        "DropConfig=0.0+0.+0.+0.,Optimizer=ADAM");
      TString trainingStrategyString("TrainingStrategy=");
      trainingStrategyString += trainingString1; // + "|" + trainingString2

      // General Options.
      TString dnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                         "WeightInitialization=XAVIER:RandomSeed=0");

      dnnOptions.Append(":");
      dnnOptions.Append(inputLayoutString);
      dnnOptions.Append(":");
      dnnOptions.Append(layoutString);
      dnnOptions.Append(":");
      dnnOptions.Append(trainingStrategyString);
      dnnOptions.Append(":");
      dnnOptions.Append(archString);

      TString dnnName = "TMVA_DNN";
      factory->BookMethod(dataloader, TMVA::Types::kDL, dnnName, dnnOptions);
   }

   /**
    ## Book Keras recurrent models

     Book the different types of recurrent models in Keras  (SimpleRNN, LSTM or GRU)

   **/

   if (useKeras) {

      for (int i = 0; i < 3; i++) {

         if (use_rnn_type[i]) {

            TString modelName = TString::Format("model_%s.h5", rnn_types[i].c_str());
            TString trainedModelName = TString::Format("trained_model_%s.h5", rnn_types[i].c_str());

            Info("TMVA_RNN_Classification", "Building recurrent keras model using a %s layer", rnn_types[i].c_str());
            // create python script which can be executed
            // create 2 conv2d layer + maxpool + dense
            TMacro m;
            m.AddLine("import tensorflow");
            m.AddLine("from tensorflow.keras.models import Sequential");
            m.AddLine("from tensorflow.keras.optimizers import Adam");
            m.AddLine("from tensorflow.keras.layers import Input, Dense, Dropout, Flatten, SimpleRNN, GRU, LSTM, Reshape, "
                      "BatchNormalization");
            m.AddLine("");
            m.AddLine("model = Sequential() ");
            m.AddLine("model.add(Reshape((10, 30), input_shape = (10*30, )))");
            // add recurrent neural network depending on type / Use option to return the full output
            if (rnn_types[i] == "LSTM")
               m.AddLine("model.add(LSTM(units=10, return_sequences=True) )");
            else if (rnn_types[i] == "GRU")
               m.AddLine("model.add(GRU(units=10, return_sequences=True) )");
            else
               m.AddLine("model.add(SimpleRNN(units=10, return_sequences=True) )");

            // m.AddLine("model.add(BatchNormalization())");
            m.AddLine("model.add(Flatten())"); // needed if returning the full time output sequence
            m.AddLine("model.add(Dense(64, activation = 'tanh')) ");
            m.AddLine("model.add(Dense(2, activation = 'sigmoid')) ");
            m.AddLine(
               "model.compile(loss = 'binary_crossentropy', optimizer = Adam(lr = 0.001), metrics = ['accuracy'])");
            m.AddLine(TString::Format("modelName = '%s'", modelName.Data()));
            m.AddLine("model.save(modelName)");
            m.AddLine("model.summary()");

            m.SaveSource("make_rnn_model.py");
            // execute
            gSystem->Exec("python make_rnn_model.py");

            if (gSystem->AccessPathName(modelName)) {
               Warning("TMVA_RNN_Classification", "Error creating Keras recurrennt model file - Skip using Keras");
               useKeras = false;
            } else {
               // book PyKeras method only if Keras model could be created
               Info("TMVA_RNN_Classification", "Booking Keras %s model", rnn_types[i].c_str());
               factory->BookMethod(dataloader, TMVA::Types::kPyKeras,
                                   TString::Format("PyKeras_%s", rnn_types[i].c_str()),
                                   TString::Format("!H:!V:VarTransform=None:FilenameModel=%s:tf.keras:"
                                                   "FilenameTrainedModel=%s:GpuOptions=allow_growth=True:"
                                                   "NumEpochs=%d:BatchSize=%d",
                                                   modelName.Data(), trainedModelName.Data(), maxepochs, batchSize));
            }
         }
      }
   }

   // use BDT in case not using Keras or TMVA DL
   if (!useKeras || !useTMVA_BDT)
      useTMVA_BDT = true;

   /**
         ## Book TMVA BDT
   **/

   if (useTMVA_BDT) {

      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                          "!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:"
                          "BaggedSampleFraction=0.5:nCuts=20:"
                          "MaxDepth=2");

   }

   /// Train all methods
   factory->TrainAllMethods();

   std::cout << "nthreads  = " << ROOT::GetThreadPoolSize() << std::endl;

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // check method

   // plot ROC curve
   auto c1 = factory->GetROCCurve(dataloader);
   c1->Draw();

   if (outputFile) outputFile->Close();
}
