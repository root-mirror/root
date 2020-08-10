/// \file
/// \ingroup tutorial_http
///  This program demonstrates simple application control via THttpServer
///  Two histogram are filled within endless loop.
///  Via published commands one can enable/disable histograms filling
///  There are also command to clear histograms content
///
///  After macro started, open in browser with url
/// ~~~
///      http://localhost:8080
/// ~~~
///
///  Histograms will be automatically displayed and
///  monitoring with interval 2000 ms started
///
/// \macro_code
///
/// \author Sergey Linev

#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "THttpServer.h"

Bool_t bFillHist = kTRUE;

void httpcontrol()
{
   // create histograms
   TH1D *hpx = new TH1D("hpx","This is the px distribution",100,-4,4);
   hpx->SetFillColor(48);
   hpx->SetDirectory(0);
   TH2D *hpxpy = new TH2D("hpxpy","py vs px",40,-4,4,40,-4,4);
   hpxpy->SetDirectory(0);

   // start http server
   THttpServer* serv = new THttpServer("http:8080");

   // One could specify location of newer version of JSROOT
   // serv->SetJSROOT("https://root.cern/js/latest/");
   // serv->SetJSROOT("http://jsroot.gsi.de/latest/");

   // register histograms
   serv->Register("/", hpx);
   serv->Register("/", hpxpy);

   // enable monitoring and
   // specify items to draw when page is opened
   serv->SetItemField("/","_monitoring","5000");
   serv->SetItemField("/","_layout","grid2x2");
   serv->SetItemField("/","_drawitem","[hpxpy,hpx,Debug]");
   serv->SetItemField("/","_drawopt","col");

   // register simple start/stop commands
   serv->RegisterCommand("/Start", "bFillHist=kTRUE;", "button;rootsys/icons/ed_execute.png");
   serv->RegisterCommand("/Stop",  "bFillHist=kFALSE;", "button;rootsys/icons/ed_interrupt.png");

   // one could hide commands and let them appear only as buttons
   serv->Hide("/Start");
   serv->Hide("/Stop");

   // register commands, invoking object methods
   serv->RegisterCommand("/ResetHPX","/hpx/->Reset()", "button;rootsys/icons/ed_delete.png");
   serv->RegisterCommand("/ResetHPXPY","/hpxpy/->Reset()", "button;rootsys/icons/bld_delete.png");


   // create debug text element, use MathJax - works only on Firefox
   serv->CreateItem("/Debug","debug output");
   serv->SetItemField("/Debug", "_kind", "Text");
   serv->SetItemField("/Debug", "value","\\(\\displaystyle{x+1\\over y-1}\\)");
   serv->SetItemField("/Debug", "mathjax", "true");

   // Fill histograms randomly
   TRandom3 random;
   Float_t px, py;
   const Long_t kUPDATE = 1000;
   Long_t cnt = 0;
   while (kTRUE) {
      if (bFillHist) {
         random.Rannor(px,py);
         hpx->Fill(px);
         hpxpy->Fill(px,py);
         cnt++;
      } else {
         gSystem->Sleep(10); // sleep minimal time
      }

      if ((cnt % kUPDATE==0) || !bFillHist) {
         // IMPORTANT: one should regularly call ProcessEvents
         // to let http server process requests

         serv->SetItemField("/Debug", "value", Form("\\(\\displaystyle{x+1\\over y-1}\\) Loop:%ld", cnt/kUPDATE));

         if (gSystem->ProcessEvents()) break;
      }
   }
}
