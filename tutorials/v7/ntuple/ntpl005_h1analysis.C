/// \file
/// \ingroup tutorial_ntuple
/// \notebook
/// Multiple RNTuple files containing data from DESY are chained.
/// This tutorial illustrates how multiple RNTuple files with identical fields can be chained and be treated
/// like a single file.
/// For reading, the tutorial shows the use of an ntuple View, which selectively accesses specific fields.
/// If a view is used for reading, there is no need to load the dictionary for custom-made objects.
/// The advantage of a view is that it directly accesses RNTuple's data buffers without making an additional
/// memory copy. It also illustrates the usage of CollectionViews.
///
/// \macro_image
/// \macro_code
///
/// \date October 2019
/// \author The ROOT Team

// NOTE: The RNTuple classes are experimental at this point.
// Functionality, interface, and data format is still subject to changes.
// Do not use for real data!

#include "ntpl005_headers/convert.hxx"
#include "ntpl005_headers/h1event.hxx"

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>

#include <TCanvas.h>
#include <TF1.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>

#include <cstdint>

// Import classes from experimental namespace for the time being
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleReader = ROOT::Experimental::RNTupleReader;

const std::string_view ntuplename = "eventdata";
// files to convert from TTree to RNTuple and read in.
const std::vector<std::string> ntupleFiles { "dstarmb.root", "dstarp1a.root", "dstarp1b.root", "dstarp2.root"};
const Double_t dxbin = (0.17-0.13)/40;   // Bin-width
const Double_t sigma = 0.0012;


Double_t fdm5(Double_t *xx, Double_t *par)
{
   Double_t x = xx[0];
   if (x <= 0.13957) return 0;
   Double_t xp3 = (x-par[3])*(x-par[3]);
   Double_t res = dxbin*(par[0]*TMath::Power(x-0.13957, par[1])
       + par[2] / 2.5066/par[4]*TMath::Exp(-xp3/2/par[4]/par[4]));
   return res;
}


Double_t fdm2(Double_t *xx, Double_t *par)
{
   Double_t x = xx[0];
   if (x <= 0.13957) return 0;
   Double_t xp3 = (x-0.1454)*(x-0.1454);
   Double_t res = dxbin*(par[0]*TMath::Power(x-0.13957, 0.25)
       + par[1] / 2.5066/sigma*TMath::Exp(-xp3/2/sigma/sigma));
   return res;
}

void ntpl005_h1analysis()
{
   // converts file .root file containing TTree to .root file containing ntuple
   for (const auto& f: ntupleFiles) {
      convert(f);
   }
   
   auto emptyModel = RNTupleModel::Create();
   // Workaround: RNTuple should only lazily load the default model
   // Note here, that a std::vector of filenames (ntupleFiles) is passed here to chain multiple files.
   auto ntuple = RNTupleReader::Open(std::move(emptyModel), "h42", ntupleFiles);

   auto hdmd = new TH1F("hdmd","dm_d",40,0.13,0.17);
   auto h2   = new TH2F("h2","ptD0 vs dm_d",30,0.135,0.165,30,-3,6);

   // Create a view for every field from which data is read.
   auto dm_dView = ntuple->GetView<float>("event.dm_d");
   auto rpd0_tView = ntuple->GetView<float>("event.rpd0_t");
   auto ptd0_dView = ntuple->GetView<float>("event.ptd0_d");

   auto ptds_dView = ntuple->GetView<float>("event.ptds_d");
   auto etads_dView = ntuple->GetView<float>("event.etads_d");
   auto ikView = ntuple->GetView<std::int32_t>("event.ik");
   auto ipiView = ntuple->GetView<std::int32_t>("event.ipi");
   auto ipisView = ntuple->GetView<std::int32_t>("event.ipis");
   auto md0_dView = ntuple->GetView<float>("event.md0_d");
   
   // CollectionViews allow to iterate through entries of its subfields
   auto trackView = ntuple->GetViewCollection("event.tracks");
   auto nhitrpView = ntuple->GetView<std::int32_t>("event.tracks.H1Event::Track.nhitrp");
   auto rstartView = ntuple->GetView<float>("event.tracks.H1Event::Track.rstart");
   auto rendView = ntuple->GetView<float>("event.tracks.H1Event::Track.rend");
   auto nlhkView = ntuple->GetView<float>("event.tracks.H1Event::Track.nlhk");
   auto nlhpiView = ntuple->GetView<float>("event.tracks.H1Event::Track.nlhpi");
   auto njetsView = ntuple->GetView<std::vector<H1Event::Jet>>("event.jets");

   for (auto i : ntuple->GetViewRange()) {
      auto ik = ikView(i) - 1;
      auto ipi = ipiView(i) - 1;
      auto ipis = ipisView(i) - 1;

      if (TMath::Abs(md0_dView(i)-1.8646) >= 0.04) continue;
      if (ptds_dView(i) <= 2.5) continue;
      if (TMath::Abs(etads_dView(i)) >= 1.5) continue;

      // nhitrpView(*trackView.GetViewRange(i).begin()) returns the value of the first element of the i-th entry.
      // nhitrpView(i) would not do the same, because here a single entry contains a vector of elements. (In this
      // case the field: "event.tracks" contains a std::vector of H1Event::Track, which is the parent field of
      // "event.tracks.H1Event::Track.nhitrp")
      if (nhitrpView(*trackView.GetViewRange(i).begin()+ik) * nhitrpView(*trackView.GetViewRange(i).begin()+ipi) <= 1) continue;
      if (rendView(*trackView.GetViewRange(i).begin()+ik) - rstartView(*trackView.GetViewRange(i).begin()+ik) <= 22) continue;
      if (rendView(*trackView.GetViewRange(i).begin()+ipi) - rstartView(*trackView.GetViewRange(i).begin()+ipi) <= 22) continue;
      if (nlhkView(*trackView.GetViewRange(i).begin()+ik) <= 0.1) continue;
      if (nlhpiView(*trackView.GetViewRange(i).begin()+ipi) <= 0.1) continue;
      if (nlhpiView(*trackView.GetViewRange(i).begin()+ipis) <= 0.1) continue;
      if (njetsView(i).size() < 1) continue;

      hdmd->Fill(dm_dView(i));
      h2->Fill(dm_dView(i),rpd0_tView(i)/0.029979*1.8646/ptd0_dView(i));
   }

   gStyle->SetOptFit();
   TCanvas *c1 = new TCanvas("c1","h1analysis analysis",10,10,800,600);
   c1->SetBottomMargin(0.15);
   hdmd->GetXaxis()->SetTitle("m_{K#pi#pi} - m_{K#pi}[GeV/c^{2}]");
   hdmd->GetXaxis()->SetTitleOffset(1.4);

   TF1 *f5 = new TF1("f5",fdm5,0.139,0.17,5);
   f5->SetParameters(1000000, .25, 2000, .1454, .001);
   hdmd->Fit("f5","lr");
   hdmd->DrawCopy();

   gStyle->SetOptFit(0);
   gStyle->SetOptStat(1100);
   TCanvas *c2 = new TCanvas("c2","tauD0",100,100,800,600);
   c2->SetGrid();
   c2->SetBottomMargin(0.15);

   TF1 *f2 = new TF1("f2",fdm2,0.139,0.17,2);
   f2->SetParameters(10000, 10);
   h2->FitSlicesX(f2,0,-1,1,"qln");
   TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");
   h2_1->GetXaxis()->SetTitle("#tau[ps]");
   h2_1->SetMarkerStyle(21);
   h2_1->DrawCopy();
   c2->Update();
   TLine *line = new TLine(0,0,0,c2->GetUymax());
   line->Draw();
}
