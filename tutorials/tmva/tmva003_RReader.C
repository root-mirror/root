/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This tutorial shows how to apply with the modern interfaces models saved in
/// TMVA XML files.
///
/// \macro_code
/// \macro_output
///
/// \date July 2019
/// \author Stefan Wunsch

using namespace TMVA::Experimental;

void train(const std::string &filename)
{
   // Create factory
   auto output = TFile::Open("TMVA.root", "RECREATE");
   auto factory = new TMVA::Factory("tmva003",
           output, "!V:!DrawProgressBar:AnalysisType=Classification");

   // Open trees with signal and background events
   auto data = TFile::Open(filename.c_str());
   auto signal = (TTree *)data->Get("TreeS");
   auto background = (TTree *)data->Get("TreeB");

   // Add variables and register the trees with the dataloader
   auto dataloader = new TMVA::DataLoader("tmva003_BDT");
   const std::vector<std::string> variables = {"var1", "var2", "var3", "var4"};
   for (const auto &var : variables) {
      dataloader->AddVariable(var);
   }
   dataloader->AddSignalTree(signal, 1.0);
   dataloader->AddBackgroundTree(background, 1.0);
   dataloader->PrepareTrainingAndTestTree("", "");

   // Train a TMVA method
   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!V:!H:NTrees=300:MaxDepth=2");
   factory->TrainAllMethods();
}

void tmva003_RReader()
{
   // First, let's train a model with TMVA.
   const std::string filename = "http://root.cern/files/tmva_class_example.root";
   train(filename);

   // Next, we load the model from the TMVA XML file.
   RReader model("tmva003_BDT/weights/tmva003_BDT.weights.xml");

   // In case you need a reminder of the names and order of the variables during
   // training, you can ask the model for it.
   auto variables = model.GetVariableNames();

   // The model can now be applied in different scenarios:
   // 1) Event-by-event inference
   // 2) Batch inference on data of multiple events
   // 3) Inference as part of an RDataFrame graph

   // 1) Event-by-event inference
   // The event-by-event inference takes the values of the variables as a std::vector<float>.
   // Note that the return value is as well a std::vector<float> since the reader
   // is also capable to process models with multiple outputs.
   auto prediction = model.Compute({0.5, 1.0, -0.2, 1.5});
   std::cout << "Single-event inference: " << prediction[0] << "\n\n";

   // 2) Batch inference on data of multiple events
   // For batch inference, the data needs to be structured as a matrix. For this
   // purpose, TMVA makes use of the RTensor class. For convenience, we use RDataFrame
   // and the AsTensor utility to make the read-out from the ROOT file.
   ROOT::RDataFrame df("TreeS", filename);
   auto df2 = df.Range(3); // Read only a small subset of the dataset
   auto x = AsTensor<float>(df2, variables);
   auto y = model.Compute(x);

   std::cout << "RTensor input for inference on data of multiple events:\n" << x << "\n\n";
   std::cout << "Prediction performed on multiple events: " << y << "\n\n";

   // 3) Perform inference as part of an RDataFrame graph
   // We write a small lambda function that performs for us the inference on
   // a dataframe to omit code duplication.
   auto make_histo = [&](const std::string &treename) {
      ROOT::RDataFrame df(treename, filename);
      auto df2 = df.Define("y", Compute<4, float>(model), variables);
      return df2.Histo1D({treename.c_str(), ";BDT score;N_{Events}", 30, -0.5, 0.5}, "y");
   };

   auto sig = make_histo("TreeS");
   auto bkg = make_histo("TreeB");

   // Make plot
   gStyle->SetOptStat(0);
   auto c = new TCanvas("", "", 800, 800);

   sig->SetLineColor(kRed);
   bkg->SetLineColor(kBlue);
   sig->SetLineWidth(2);
   bkg->SetLineWidth(2);
   bkg->Draw("HIST");
   sig->Draw("HIST SAME");

   TLegend legend(0.7, 0.7, 0.89, 0.89);
   legend.SetBorderSize(0);
   legend.AddEntry("TreeS", "Signal", "l");
   legend.AddEntry("TreeB", "Background", "l");
   legend.Draw();

   c->DrawClone();
}
