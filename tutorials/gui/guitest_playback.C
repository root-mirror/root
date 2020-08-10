/// \file
/// \ingroup tutorial_gui
/// This macro plays a recorded ROOT session showing how to perform various interactive GUI operations with the guitest.C macro
/// While replaying the session, several temporary macros (guitest0xx.C) macros will be saved.
/// These files will be later on compared with some reference values to verify the validity of different parts of the test.
///
/// \macro_code
///
/// \author Bertrand Bellenot

#include "TSystem.h"
#include "TSystemFile.h"
#include "TRecorder.h"
#include "Riostream.h"

// count characters in the file, skipping cr/lf
Int_t file_size(const char *filename)
{
   FILE *lunin;
   Int_t c, wc = 0;

   lunin = fopen(filename, "rb");
   if (lunin == 0) return -1;
   while (!feof(lunin)) {
      c = fgetc(lunin);
      if (c != 0x0d && c != 0x0a)
         wc++;
   }
   fclose(lunin);
   return wc;
}

// main function
void guitest_playback()
{
   Int_t i;
   Int_t guitest_ref[11], guitest_err[11], guitest_size[11];

   gBenchmark->Start("guitest_playback");

   // first delete old files, if any
   for (i=0;i<11;++i) {
      gSystem->Unlink(TString::Format("guitest%03d.C", i+1));
   }

   TRecorder r("http://root.cern/files/guitest_playback.root");

   // wait for the recorder to finish the replay
   while (r.GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
   }

   for (i=0;i<11;++i) {
      guitest_ref[i] = 0;
      guitest_err[i] = 100;
      guitest_size[i] = file_size(TString::Format("guitest%03d.C", i+1));
   }

   guitest_ref[0]  = 23319;
   guitest_ref[1]  =  5633;
   guitest_ref[2]  = 14939;
   guitest_ref[3]  =  9459;
   guitest_ref[4]  =  5351;
   guitest_ref[5]  = 22982;
   guitest_ref[6]  = 23812;
   guitest_ref[7]  = 23869;
   guitest_ref[8]  = 23918;
   guitest_ref[9]  = 24067;
   guitest_ref[10] = 65517;

   printf("**********************************************************************\n");
   printf("*  Results of guitest_playback.C                                     *\n");
   printf("**********************************************************************\n");

   for (i=0;i<11;++i) {
      printf("guitest %02d: output............................................", i+1);
      if (TMath::Abs(guitest_ref[i] - guitest_size[i]) <= guitest_err[i]) {
         printf("..... OK\n");
         // delete successful tests, keep only the failing ones (for verification)
         gSystem->Unlink(TString::Format("guitest%03d.C", i+1));
      }
      else {
         printf(". FAILED\n");
      }
   }
   printf("**********************************************************************\n");
   gBenchmark->Show("guitest_playback");
}
