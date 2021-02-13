/// \file
/// \ingroup tutorial_geom
/// This macro generates
/// with 2 views of the NA49 detector using the old obsolete geometry package.
///
/// \macro_image
/// \macro_code
///
/// \author Rene Brun

void na49view() {
   TCanvas *c1 = new TCanvas("c1","The NA49 canvas",200,10,700,780);

   gBenchmark->Start("na49view");

   TPad *all = new TPad("all","A Global view of NA49",0.02,0.02,0.48,0.82,28);
   TPad *tof = new TPad("tof","One Time Of Flight element",0.52,0.02,0.98,0.82,28);
   all->Draw();
   tof->Draw();
   TPaveLabel *na49title = new TPaveLabel(0.04,0.86,0.96,0.98,"Two views of the NA49 detector");
   na49title->SetFillColor(32);
   na49title->Draw();
   //
   TFile *nageom = new TFile("na49.root");
   if (!nageom || nageom->IsZombie()) return;
   TGeometry *n49 =(TGeometry*)gROOT->FindObject("na49");
   n49->SetBomb(1.2);
   n49->cd();     //Set current geometry
   all->cd();     //Set current pad
   n49->Draw();
   c1->Update();
   tof->cd();
   TNode *TOFR1 = n49->GetNode("TOFR1");
   TOFR1->Draw();
   c1->Update();

   gBenchmark->Show("na49view");

   // To have a better and dynamic view of any of these pads,
   // you can click with the middle button of your mouse to select it.
   // Then select "View with x3d" in the VIEW menu of the Canvas.
   // Once in x3d, you are in wireframe mode by default.
   // You can switch to:
   //   - Hidden Line mode by typing E
   //   - Solid mode by typing R
   //   - Wireframe mode by typing W
   //   - Stereo mode by clicking S (and you need special glasses)
   //   - To leave x3d type Q
}
