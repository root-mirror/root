/// \file
/// \ingroup tutorial_fit
/// This file will test all the transient frames (aka Dialog windows)
/// displayed in the fitpanel, as the rest of the functionality is
/// tried automatically with the UnitTest.C unit.
///
/// This implies trying the Set Parameters dialog and the Advanced one.
///
/// At every operation, a png file will be saved. These files will be
/// later on compared with some references values, to have an estimation
/// of the goodness of the test.
///
/// \macro_code
///
/// \author David Gonzalez Maline

#include "TSystem.h"
#include "TSystemFile.h"
#include "TRecorder.h"
#include "Riostream.h"

Int_t file_size(char *filename)
{
   FileStat_t fs;
   gSystem->GetPathInfo(filename, fs);
   return (Int_t)fs.fSize;
}

void fitpanel_playback()
{
   r = new TRecorder();
   r->Replay("http://root.cern/files/fitpanel_playback.root");

   // wait for the recorder to finish the replay
   while (r->GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
   }

   Int_t Step_Err    = 100;
   Int_t Step1_Ref   = 15691;
   Int_t Step2_Ref   = 15691;
   Int_t Step3_Ref   = 17632;
   Int_t Step4_Ref   = 12305;
   Int_t Step5_Ref   = 11668;

   Int_t Step1_Size   =  file_size("Step1.png");
   Int_t Step2_Size   =  file_size("Step2.png");
   Int_t Step3_Size   =  file_size("Step3.png");
   Int_t Step4_Size   =  file_size("Step4.png");
   Int_t Step5_Size   =  file_size("Step5.png");


   cout << "**********************************************************************" <<endl;
   cout << "*  Report of fitpanel_playback.C                                     *" <<endl;
   cout << "**********************************************************************" <<endl;

   if (TMath::Abs(Step1_Ref-Step1_Size) <= Step_Err) {
      cout << "Step1: ............................................................ OK" <<endl;
   } else {
      cout << "Step1: ........................................................ FAILED" <<endl;
   }

   if (TMath::Abs(Step2_Ref-Step2_Size) <= Step_Err) {
      cout << "Step2: ............................................................ OK" <<endl;
   } else {
      cout << "Step2: ........................................................ FAILED" <<endl;
   }

   if (TMath::Abs(Step3_Ref-Step3_Size) <= Step_Err) {
      cout << "Step3: ............................................................ OK" <<endl;
   } else {
      cout << "Step3: ........................................................ FAILED" <<endl;
   }

   if (TMath::Abs(Step4_Ref-Step4_Size) <= Step_Err) {
      cout << "Step4: ............................................................ OK" <<endl;
   } else {
      cout << "Step4: ........................................................ FAILED" <<endl;
   }

   if (TMath::Abs(Step5_Ref-Step5_Size) <= Step_Err) {
      cout << "Step5: ............................................................ OK" <<endl;
   } else {
      cout << "Step5: ........................................................ FAILED" <<endl;
   }
   cout << "**********************************************************************" <<endl;

}
