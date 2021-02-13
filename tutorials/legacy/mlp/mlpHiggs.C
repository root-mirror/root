/// \file
/// \ingroup tutorial_mlp
/// \notebook
/// Example of a Multi Layer Perceptron
/// For a LEP search for invisible Higgs boson, a neural network
/// was used to separate the signal from the background passing
/// some selection cuts. Here is a simplified version of this network,
/// taking into account only WW events.
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Christophe Delaere

void mlpHiggs(Int_t ntrain=100) {
   const char *fname = "mlpHiggs.root";
   TFile *input = 0;
   if (!gSystem->AccessPathName(fname)) {
      input = TFile::Open(fname);
   } else if (!gSystem->AccessPathName(Form("%s/legacy/mlp/%s", TROOT::GetTutorialDir().Data(), fname))) {
      input = TFile::Open(Form("%s/legacy/mlp/%s", TROOT::GetTutorialDir().Data(), fname));
   } else {
      printf("accessing %s file from http://root.cern.ch/files\n",fname);
      input = TFile::Open(Form("http://root.cern.ch/files/%s",fname));
   }
   if (!input) return;

   TTree *sig_filtered = (TTree *) input->Get("sig_filtered");
   TTree *bg_filtered = (TTree *) input->Get("bg_filtered");
   TTree *simu = new TTree("MonteCarlo", "Filtered Monte Carlo Events");
   Float_t ptsumf, qelep, nch, msumf, minvis, acopl, acolin;
   Int_t type;
   sig_filtered->SetBranchAddress("ptsumf", &ptsumf);
   sig_filtered->SetBranchAddress("qelep",  &qelep);
   sig_filtered->SetBranchAddress("nch",    &nch);
   sig_filtered->SetBranchAddress("msumf",  &msumf);
   sig_filtered->SetBranchAddress("minvis", &minvis);
   sig_filtered->SetBranchAddress("acopl",  &acopl);
   sig_filtered->SetBranchAddress("acolin", &acolin);
   bg_filtered->SetBranchAddress("ptsumf", &ptsumf);
   bg_filtered->SetBranchAddress("qelep",  &qelep);
   bg_filtered->SetBranchAddress("nch",    &nch);
   bg_filtered->SetBranchAddress("msumf",  &msumf);
   bg_filtered->SetBranchAddress("minvis", &minvis);
   bg_filtered->SetBranchAddress("acopl",  &acopl);
   bg_filtered->SetBranchAddress("acolin", &acolin);
   simu->Branch("ptsumf", &ptsumf, "ptsumf/F");
   simu->Branch("qelep",  &qelep,  "qelep/F");
   simu->Branch("nch",    &nch,    "nch/F");
   simu->Branch("msumf",  &msumf,  "msumf/F");
   simu->Branch("minvis", &minvis, "minvis/F");
   simu->Branch("acopl",  &acopl,  "acopl/F");
   simu->Branch("acolin", &acolin, "acolin/F");
   simu->Branch("type",   &type,   "type/I");
   type = 1;
   Int_t i;
   for (i = 0; i < sig_filtered->GetEntries(); i++) {
      sig_filtered->GetEntry(i);
      simu->Fill();
   }
   type = 0;
   for (i = 0; i < bg_filtered->GetEntries(); i++) {
      bg_filtered->GetEntry(i);
      simu->Fill();
   }
   // Build and train the NN ptsumf is used as a weight since we are primarily
   // interested  by high pt events.
   // The datasets used here are the same as the default ones.
   TMultiLayerPerceptron *mlp =
      new TMultiLayerPerceptron("@msumf,@ptsumf,@acolin:5:3:type",
                                "ptsumf",simu,"Entry$%2","(Entry$+1)%2");
   mlp->Train(ntrain, "text,graph,update=10");
   mlp->Export("test","python");
   // Use TMLPAnalyzer to see what it looks for
   TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
   mlpa_canvas->Divide(2,2);
   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
   mlpa_canvas->cd(1);
   // shows how each variable influences the network
   ana.DrawDInputs();
   mlpa_canvas->cd(2);
   // shows the network structure
   mlp->Draw();
   mlpa_canvas->cd(3);
   // draws the resulting network
   ana.DrawNetwork(0,"type==1","type==0");
   mlpa_canvas->cd(4);
   // Use the NN to plot the results for each sample
   // This will give approx. the same result as DrawNetwork.
   // All entries are used, while DrawNetwork focuses on
   // the test sample. Also the xaxis range is manually set.
   TH1F *bg = new TH1F("bgh", "NN output", 50, -.5, 1.5);
   TH1F *sig = new TH1F("sigh", "NN output", 50, -.5, 1.5);
   bg->SetDirectory(0);
   sig->SetDirectory(0);
   Double_t params[3];
   for (i = 0; i < bg_filtered->GetEntries(); i++) {
      bg_filtered->GetEntry(i);
      params[0] = msumf;
      params[1] = ptsumf;
      params[2] = acolin;
      bg->Fill(mlp->Evaluate(0, params));
   }
   for (i = 0; i < sig_filtered->GetEntries(); i++) {
      sig_filtered->GetEntry(i);
      params[0] = msumf;
      params[1] = ptsumf;
      params[2] = acolin;
      sig->Fill(mlp->Evaluate(0,params));
   }
   bg->SetLineColor(kBlue);
   bg->SetFillStyle(3008);   bg->SetFillColor(kBlue);
   sig->SetLineColor(kRed);
   sig->SetFillStyle(3003); sig->SetFillColor(kRed);
   bg->SetStats(0);
   sig->SetStats(0);
   bg->Draw();
   sig->Draw("same");
   TLegend *legend = new TLegend(.75, .80, .95, .95);
   legend->AddEntry(bg, "Background (WW)");
   legend->AddEntry(sig, "Signal (Higgs)");
   legend->Draw();
   mlpa_canvas->cd(0);
   delete input;
}
