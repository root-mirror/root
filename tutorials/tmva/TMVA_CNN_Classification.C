/// \file
/// \ingroup tutorial_tmva
/// \notebook
///  TMVA Classification Example Using a Convolutional Neural Network
///
/// This is an example of using a CNN in TMVA. We do classification using a toy image data set
/// that is generated when running the example macro
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta

/***

    # TMVA Classification Example Using a Convolutional Neural Network


**/

/// Helper function to create input images data
/// we create a signal and background 2D histograms from 2d gaussians
/// with a location (means in X and Y)  different for each event
/// The difference between signal and background is in the gaussian width.
/// The width for the bakground gaussian is slightly larger than the signal width by few % values
///
///
void MakeImagesTree(int n, int nh, int nw)
{

   // image size (nh x nw)
   const int ntot = nh * nw;
   const TString fileOutName = TString::Format("images_data_%dx%d.root", nh, nw);

   const int nRndmEvts = 10000; // number of events we use to fill each image
   double delta_sigma = 0.1;    // 5% difference in the sigma
   double pixelNoise = 5;

   double sX1 = 3;
   double sY1 = 3;
   double sX2 = sX1 + delta_sigma;
   double sY2 = sY1 - delta_sigma;

   auto h1 = new TH2D("h1", "h1", nh, 0, 10, nw, 0, 10);
   auto h2 = new TH2D("h2", "h2", nh, 0, 10, nw, 0, 10);

   auto f1 = new TF2("f1", "xygaus");
   auto f2 = new TF2("f2", "xygaus");
   TTree sgn("sig_tree", "signal_tree");
   TTree bkg("bkg_tree", "bakground_tree");

   TFile f(fileOutName, "RECREATE");

   std::vector<float> x1(ntot);
   std::vector<float> x2(ntot);

   // create signal and background trees with a single branch
   // an std::vector<float> of size nh x nw containing the image data

   std::vector<float> *px1 = &x1;
   std::vector<float> *px2 = &x2;

   bkg.Branch("vars", "std::vector<float>", &px1);
   sgn.Branch("vars", "std::vector<float>", &px2);

   // std::cout << "create tree " << std::endl;

   sgn.SetDirectory(&f);
   bkg.SetDirectory(&f);

   f1->SetParameters(1, 5, sX1, 5, sY1);
   f2->SetParameters(1, 5, sX2, 5, sY2);
   gRandom->SetSeed(0);
   std::cout << "Filling ROOT tree " << std::endl;
   for (int i = 0; i < n; ++i) {
      if (i % 1000 == 0)
         std::cout << "Generating image event ... " << i << std::endl;
      h1->Reset();
      h2->Reset();
      // generate random means in range [3,7] to be not too much on the border
      f1->SetParameter(1, gRandom->Uniform(3, 7));
      f1->SetParameter(3, gRandom->Uniform(3, 7));
      f2->SetParameter(1, gRandom->Uniform(3, 7));
      f2->SetParameter(3, gRandom->Uniform(3, 7));

      h1->FillRandom("f1", nRndmEvts);
      h2->FillRandom("f2", nRndmEvts);

      for (int k = 0; k < nh; ++k) {
         for (int l = 0; l < nw; ++l) {
            int m = k * nw + l;
            // add some noise in each bin
            x1[m] = h1->GetBinContent(k + 1, l + 1) + gRandom->Gaus(0, pixelNoise);
            x2[m] = h2->GetBinContent(k + 1, l + 1) + gRandom->Gaus(0, pixelNoise);
         }
      }
      sgn.Fill();
      bkg.Fill();
   }
   sgn.Write();
   bkg.Write();

   Info("MakeImagesTree", "Signal and background tree with images data written to the file %s", f.GetName());
   sgn.Print();
   bkg.Print();
   f.Close();
}

void TMVA_CNN_Classification(std::vector<bool> opt = {1, 1, 1, 1, 1})
{

   bool useTMVACNN = (opt.size() > 0) ? opt[0] : false;
   bool useKerasCNN = (opt.size() > 1) ? opt[1] : false;
   bool useTMVADNN = (opt.size() > 2) ? opt[2] : false;
   bool useTMVABDT = (opt.size() > 3) ? opt[3] : false;
   bool usePyTorchCNN = (opt.size() > 4) ? opt[4] : false;
#ifndef R__HAS_TMVACPU
#ifndef R__HAS_TMVAGPU
   Warning("TMVA_CNN_Classification",
           "TMVA is not build with GPU or CPU multi-thread support. Cannot use TMVA Deep Learning for CNN");
   useTMVACNN = false;
#endif
#endif

   bool writeOutputFile = true;

   int num_threads = 0;  // use default threads

   TMVA::Tools::Instance();

   // do enable MT running
   if (num_threads >= 0) {
      ROOT::EnableImplicitMT(num_threads);
      if (num_threads > 0) gSystem->Setenv("OMP_NUM_THREADS", TString::Format("%d",num_threads));
   }
   else
      gSystem->Setenv("OMP_NUM_THREADS", "1");

   std::cout << "Running with nthreads  = " << ROOT::GetThreadPoolSize() << std::endl;

#ifdef R__HAS_PYMVA
   gSystem->Setenv("KERAS_BACKEND", "tensorflow");
   // for using Keras
   TMVA::PyMethodBase::PyInitialize();
#else
   useKerasCNN = false;
#endif

   TFile *outputFile = nullptr;
   if (writeOutputFile)
      outputFile = TFile::Open("TMVA_CNN_ClassificationOutput.root", "RECREATE");

   /***
       ## Create TMVA Factory

    Create the Factory class. Later you can choose the methods
    whose performance you'd like to investigate.

    The factory is the major TMVA object you have to interact with. Here is the list of parameters you need to pass

    - The first argument is the base of the name of all the output
    weightfiles in the directory weight/ that will be created with the
    method parameters

    - The second argument is the output file for the training results

    - The third argument is a string option defining some general configuration for the TMVA session.
      For example all TMVA output can be suppressed by removing the "!" (not) in front of the "Silent" argument in the
   option string

    - note that we disable any pre-transformation of the input variables and we avoid computing correlations between
   input variables
   ***/

   TMVA::Factory factory(
      "TMVA_CNN_Classification", outputFile,
      "!V:ROC:!Silent:Color:AnalysisType=Classification:Transformations=None:!Correlations");

   /***

       ## Declare DataLoader(s)

       The next step is to declare the DataLoader class that deals with input variables

       Define the input variables that shall be used for the MVA training
       note that you may also use variable expressions, which can be parsed by TTree::Draw( "expression" )]

       In this case the input data consists of an image of 16x16 pixels. Each single pixel is a branch in a ROOT TTree

   **/

   TMVA::DataLoader *loader = new TMVA::DataLoader("dataset");

   /***

       ## Setup Dataset(s)

       Define input data file and signal and background trees

   **/

   int imgSize = 16 * 16;
   TString inputFileName = "images_data_16x16.root";

   bool fileExist = !gSystem->AccessPathName(inputFileName);

   // if file does not exists create it
   if (!fileExist) {
      MakeImagesTree(5000, 16, 16);
   }

   // TString inputFileName = "tmva_class_example.root";

   auto inputFile = TFile::Open(inputFileName);
   if (!inputFile) {
      Error("TMVA_CNN_Classification", "Error opening input file %s - exit", inputFileName.Data());
      return;
   }

   // --- Register the training and test trees

   TTree *signalTree = (TTree *)inputFile->Get("sig_tree");
   TTree *backgroundTree = (TTree *)inputFile->Get("bkg_tree");

   int nEventsSig = signalTree->GetEntries();
   int nEventsBkg = backgroundTree->GetEntries();

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   loader->AddSignalTree(signalTree, signalWeight);
   loader->AddBackgroundTree(backgroundTree, backgroundWeight);

   /// add event variables (image)
   /// use new method (from ROOT 6.20 to add a variable array for all image data)
   loader->AddVariablesArray("vars", imgSize);

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   // loader->SetBackgroundWeightExpression( "weight" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //    loader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // It is possible also to specify the number of training and testing events,
   // note we disable the computation of the correlation matrix of the input variables

   int nTrainSig = 0.8 * nEventsSig;
   int nTrainBkg = 0.8 * nEventsBkg;

   // build the string options for DataLoader::PrepareTrainingAndTestTree
   TString prepareOptions = TString::Format(
      "nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:SplitSeed=100:NormMode=NumEvents:!V:!CalcCorrelations",
      nTrainSig, nTrainBkg);

   loader->PrepareTrainingAndTestTree(mycuts, mycutb, prepareOptions);

   /***

       DataSetInfo              : [dataset] : Added class "Signal"
       : Add Tree sig_tree of type Signal with 10000 events
       DataSetInfo              : [dataset] : Added class "Background"
       : Add Tree bkg_tree of type Background with 10000 events



   **/

   // signalTree->Print();

   /****
        # Booking Methods

        Here we book the TMVA methods. We book a Boosted Decision Tree method (BDT)

   **/

   // Boosted Decision Trees
   if (useTMVABDT) {
      factory.BookMethod(loader, TMVA::Types::kBDT, "BDT",
                         "!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:"
                         "UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
   }
   /**

      #### Booking Deep Neural Network

      Here we book the DNN of TMVA. See the example TMVA_Higgs_Classification.C for a detailed description of the
   options

   **/

   if (useTMVADNN) {

      TString layoutString(
         "Layout=DENSE|100|RELU,BNORM,DENSE|100|RELU,BNORM,DENSE|100|RELU,BNORM,DENSE|100|RELU,DENSE|1|LINEAR");

      // Training strategies
      // one can catenate several training strings with different parameters (e.g. learning rates or regularizations
      // parameters) The training string must be concatenates with the `|` delimiter
      TString trainingString1("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
                              "ConvergenceSteps=5,BatchSize=100,TestRepetitions=1,"
                              "MaxEpochs=20,WeightDecay=1e-4,Regularization=None,"
                              "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.");

      TString trainingStrategyString("TrainingStrategy=");
      trainingStrategyString += trainingString1; // + "|" + trainingString2 + ....

      // Build now the full DNN Option string

      TString dnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                         "WeightInitialization=XAVIER");
      dnnOptions.Append(":");
      dnnOptions.Append(layoutString);
      dnnOptions.Append(":");
      dnnOptions.Append(trainingStrategyString);

      TString dnnMethodName = "TMVA_DNN_CPU";
// use GPU if available
#ifdef R__HAS_TMVAGPU
      dnnOptions += ":Architecture=GPU";
      dnnMethodName = "TMVA_DNN_GPU";
#elif defined(R__HAS_TMVACPU)
      dnnOptions += ":Architecture=CPU";
#endif

      factory.BookMethod(loader, TMVA::Types::kDL, dnnMethodName, dnnOptions);
   }

   /***
    ### Book Convolutional Neural Network in TMVA

    For building a CNN one needs to define

    -  Input Layout :  number of channels (in this case = 1)  | image height | image width
    -  Batch Layout :  batch size | number of channels | image size = (height*width)

    Then one add Convolutional layers and MaxPool layers.

    -  For Convolutional layer the option string has to be:
       - CONV | number of units | filter height | filter width | stride height | stride width | padding height | paddig
   width | activation function

       - note in this case we are using a filer 3x3 and padding=1 and stride=1 so we get the output dimension of the
   conv layer equal to the input

      - note we use after the first convolutional layer a batch normalization layer. This seems to help significatly the
   convergence

     - For the MaxPool layer:
        - MAXPOOL  | pool height | pool width | stride height | stride width

    The RESHAPE layer is needed to flatten the output before the Dense layer


    Note that to run the CNN is required to have CPU  or GPU support

   ***/

   if (useTMVACNN) {

      TString inputLayoutString("InputLayout=1|16|16");

      // Batch Layout
      TString layoutString("Layout=CONV|10|3|3|1|1|1|1|RELU,BNORM,CONV|10|3|3|1|1|1|1|RELU,MAXPOOL|2|2|1|1,"
                           "RESHAPE|FLAT,DENSE|100|RELU,DENSE|1|LINEAR");

      // Training strategies.
      TString trainingString1("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
                              "ConvergenceSteps=5,BatchSize=100,TestRepetitions=1,"
                              "MaxEpochs=20,WeightDecay=1e-4,Regularization=None,"
                              "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.0");

      TString trainingStrategyString("TrainingStrategy=");
      trainingStrategyString +=
         trainingString1; // + "|" + trainingString2 + "|" + trainingString3; for concatenating more training strings

      // Build full CNN Options.
      TString cnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                         "WeightInitialization=XAVIER");

      cnnOptions.Append(":");
      cnnOptions.Append(inputLayoutString);
      cnnOptions.Append(":");
      cnnOptions.Append(layoutString);
      cnnOptions.Append(":");
      cnnOptions.Append(trainingStrategyString);

      //// New DL (CNN)
      TString cnnMethodName = "TMVA_CNN_CPU";
// use GPU if available
#ifdef R__HAS_TMVAGPU
      cnnOptions += ":Architecture=GPU";
      cnnMethodName = "TMVA_CNN_GPU";
#else
      cnnOptions += ":Architecture=CPU";
      cnnMethodName = "TMVA_CNN_CPU";
#endif

      factory.BookMethod(loader, TMVA::Types::kDL, cnnMethodName, cnnOptions);
   }

   /**
      ### Book Convolutional Neural Network in Keras using a generated model

   **/

   if (useKerasCNN) {

      Info("TMVA_CNN_Classification", "Building convolutional keras model");
      // create python script which can be executed
      // crceate 2 conv2d layer + maxpool + dense
      TMacro m;
      m.AddLine("import tensorflow");
      m.AddLine("from tensorflow.keras.models import Sequential");
      m.AddLine("from tensorflow.keras.optimizers import Adam");
      m.AddLine(
         "from tensorflow.keras.layers import Input, Dense, Dropout, Flatten, Conv2D, MaxPooling2D, Reshape, BatchNormalization");
      m.AddLine("");
      m.AddLine("model = Sequential() ");
      m.AddLine("model.add(Reshape((16, 16, 1), input_shape = (256, )))");
      m.AddLine("model.add(Conv2D(10, kernel_size = (3, 3), kernel_initializer = 'glorot_normal',activation = "
                "'relu', padding = 'same'))");
      m.AddLine("model.add(BatchNormalization())");
      m.AddLine("model.add(Conv2D(10, kernel_size = (3, 3), kernel_initializer = 'glorot_normal',activation = "
                "'relu', padding = 'same'))");
      // m.AddLine("model.add(BatchNormalization())");
      m.AddLine("model.add(MaxPooling2D(pool_size = (2, 2), strides = (1,1))) ");
      m.AddLine("model.add(Flatten())");
      m.AddLine("model.add(Dense(256, activation = 'relu')) ");
      m.AddLine("model.add(Dense(2, activation = 'sigmoid')) ");
      m.AddLine("model.compile(loss = 'binary_crossentropy', optimizer = Adam(lr = 0.001), metrics = ['accuracy'])");
      m.AddLine("model.save('model_cnn.h5')");
      m.AddLine("model.summary()");

      m.SaveSource("make_cnn_model.py");
      // execute
      gSystem->Exec("python make_cnn_model.py");

      if (gSystem->AccessPathName("model_cnn.h5")) {
         Warning("TMVA_CNN_Classification", "Error creating Keras model file - skip using Keras");
      } else {
         // book PyKeras method only if Keras model could be created
         Info("TMVA_CNN_Classification", "Booking tf.Keras CNN model");
         factory.BookMethod(
            loader, TMVA::Types::kPyKeras, "PyKeras",
            "H:!V:VarTransform=None:FilenameModel=model_cnn.h5:tf.keras:"
            "FilenameTrainedModel=trained_model_cnn.h5:NumEpochs=20:BatchSize=100:"
            "GpuOptions=allow_growth=True"); // needed for RTX NVidia card and to avoid TF allocates all GPU memory
      }
   }

   if (usePyTorchCNN) {

      Info("TMVA_CNN_Classification", "Using Convolutional PyTorch Model");
      TString pyTorchFileName = gROOT->GetTutorialDir() + TString("/tmva/PyTorch_Generate_CNN_Model.py");
      // check that pytorch can be imported and file defining the model and used later when booking the method is existing
      if (gSystem->Exec("python -c 'import torch'")  || gSystem->AccessPathName(pyTorchFileName) ) {
         Warning("TMVA_CNN_Classification", "PyTorch is not installed or model building file is not existing - skip using PyTorch");
      }
      else {
         // book PyTorch method only if PyTorch model could be created
         Info("TMVA_CNN_Classification", "Booking PyTorch CNN model");
         TString methodOpt = "H:!V:VarTransform=None:FilenameModel=PyTorchModelCNN.pt:"
                             "FilenameTrainedModel=PyTorchTrainedModelCNN.pt:NumEpochs=20:BatchSize=100";
         methodOpt += TString(":UserCode=") + pyTorchFileName;
         factory.BookMethod(loader, TMVA::Types::kPyTorch, "PyTorch", methodOpt);
      }
   }


   ////  ## Train Methods

   factory.TrainAllMethods();

   /// ## Test and Evaluate Methods

   factory.TestAllMethods();

   factory.EvaluateAllMethods();

   /// ## Plot ROC Curve

   auto c1 = factory.GetROCCurve(loader);
   c1->Draw();

   // close outputfile to save output file
   outputFile->Close();
}
