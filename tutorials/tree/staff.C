/// \file
/// \ingroup tutorial_tree
/// \notebook -nodraw
/// Create a plot of the data in `cernstaff.root`
/// To create `cernstaff.root`, execute tutorial `$ROOTSYS/tutorials/tree/cernbuild.C`
///
/// \macro_image
/// \macro_code
///
/// \author Rene Brun

void staff() {
   auto f = TFile::Open("cernstaff.root");
   TTree *T = nullptr;
   f->GetObject("T",T);
   T->Draw("Grade:Age:Cost:Division:Nation","","gl5d");
   if (gPad) gPad->Print("staff.C.png");
}
