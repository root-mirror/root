/// \file
/// \ingroup tutorial_graphics
/// This macro plays a recorded ROOT session showing how to perform various
/// interactive graphical editing operations. The initial graphics setup
/// was created using the following root commands:
/// ~~~{.cpp}
///     TRecorder *t = new TRecorder();
///     t->Start("graphedit_playback.root");
///     gStyle->SetPalette(1);
///     TCanvas *c2 = new TCanvas("c2","c2",0,0,700,500);
///     TH2F* h2 = new TH2F("h2","Random 2D Gaussian",40,-4,4,40,-4,4);
///     h2->SetDirectory(0);
///     TRandom r;
///     for (int i=0;i<50000;i++) h2->Fill(r.Gaus(),r.Gaus());
///     h2->Draw();
///     gPad->Update();
///     TCanvas *c1 = new TCanvas("c1","c1",0,0,700,500);
///     TH1F* h1 = new TH1F("h1","Random 1D Gaussian",100,-4,4);
///     h1->SetDirectory(0);
///     h1->FillRandom("gaus",10000);
///     h1->Draw();
///     gPad->Update();
///
///     // Here the following "sketch" was done.
///
///     t->Stop();
/// ~~~
/// Note: The previous commands should be copy/pasted into a ROOT session, not
/// executed as a macro.
///
/// ### The interactive editing shows
///     - Object editing using object editors
///     - Direct editing on the graphics canvas
///     - Saving PS and bitmap files.
///     - Saving as a .C file: C++ code corresponding to the modifications is saved.
///
/// ### The sketch of the recorded actions is
///
///    #### On the canvas c1
///      - Open View/Editor
///      - Select histogram
///      - Change fill style
///      - Change fill color
///      - Move stat box
///      - Change fill color
///      - Move title
///      - Change fill color using wheel color
///      - Select Y axis
///      - Change axis title
///      - Select X axis
///      - Change axis title
///      - Select histogram
///      - Go in binning
///      - Change range
///      - Move range
///      - On the canvas menu set grid Y
///      - On the canvas menu set grid X
///      - On the canvas menu set log Y
///      - Increase the range
///      - Close View/Editor
///      - Open the Tool Bar
///      - Create a text "Comment"
///      - Create an arrow
///      - Change the arrow size
///      - Close the Tool Bar
///      - Save as PS file
///      - Save as C file
///      - Close c1
///    #### On the canvas c2
///      - Open View/Editor
///      - Select histogram
///      - Select COL
///      - Select Palette
///      - Move Stats
///      - Select Overflows
///      - Select histogram
///      - Select 3D
///      - Select SURF1
///      - Rotate Surface
///      - Go in binning
///      - Change X range
///      - Change Y range
///      - Close View/Editor
///      - Save as GIF file
///      - Save as C file
///      - Close c2
/// \author Olivier Couet

Int_t file_size(const char *filename)
{
   FileStat_t fs;
   gSystem->GetPathInfo(filename, fs);
   return (Int_t)fs.fSize;
}

void graph_edit_playback()
{
   TRecorder *r = new TRecorder();
   r->Replay("http://root.cern/files/graphedit_playback.root");

   // wait for the recorder to finish the replay
   while (r->GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
   }

   Int_t c1_ps_Ref  = 11592 , c1_ps_Err  = 600;
   Int_t c1_C_Ref   =  4729 , c1_C_Err   = 200;
   Int_t c2_gif_Ref = 21184 , c2_gif_Err = 500;
   Int_t c2_C_Ref   = 35471 , c2_C_Err   = 1500;

   Int_t c1_ps  = file_size("c1.ps");
   Int_t c1_C   = file_size("c1.C");
   Int_t c2_gif = file_size("c2.gif");
   Int_t c2_C   = file_size("c2.C");

   cout << "**********************************************************************" <<endl;
   cout << "*  Report of graph_edit_playback.C                                   *" <<endl;
   cout << "**********************************************************************" <<endl;

   if (TMath::Abs(c1_ps_Ref-c1_ps) <= c1_ps_Err) {
      cout << "Canvas c1: PS output............................................... OK" <<endl;
   } else {
      cout << "Canvas c1: PS output........................................... FAILED" <<endl;
   }
   if (TMath::Abs(c1_C_Ref-c1_C) <= c1_C_Err) {
      cout << "           C output................................................ OK" <<endl;
   } else {
      cout << "           C output............................................ FAILED" <<endl;
   }
   if (TMath::Abs(c2_gif_Ref-c2_gif) <= c2_gif_Err) {
      cout << "Canvas c2: GIF output.............................................. OK" <<endl;
   } else {
      cout << "Canvas c2: GIF output.......................................... FAILED" <<endl;
   }
   if (TMath::Abs(c2_C_Ref-c2_C) <= c2_C_Err) {
      cout << "           C output................................................ OK" <<endl;
   } else {
      cout << "           C output............................................ FAILED" <<endl;
   }
   cout << "**********************************************************************" <<endl;

}
