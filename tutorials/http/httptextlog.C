/// \file
/// \ingroup tutorial_http
///  This macro demonstrates custom access and custom drawing for TMsgList class
///  Main motivation for this example - demonstrate how traffic between server and
///  client can be minimized and how one could build custom UI in the browser
///
///  TMsgList in this tutorial keep last N messages, numbering each with unique identifier
///  There is TMsgList::Select() method which selects messages from the list
///  If one specifies identifier, only messages newer than this identifier are selected
///  In the selection list (TList object of TObjString) first item always identifier for
///  the latest message in the list
///
///  In JavaScript code (httptextlog.js) one uses Select() method to receive latest
///  messages, which not yet been seen in the browser and display them as text
///  At maximum, 1000 elements are preserved in the browser.
///
///  After macro started, one could open in browser address
/// ~~~
///    http://localhost:8080?item=log
/// ~~~
///  One could either click item again or enable monitoring to always receive latest messages
///  Or one could open only this output and nothing else:
/// ~~~
///     http://localhost:8080/log/draw.htm?monitoring=2000
/// ~~~
///  In last case it could be used in iframe, also it requires less code to load on the page
///
/// \macro_code
///
/// \author Sergey Linev


#include <stdio.h>
#include <string.h>

#include "TNamed.h"
#include "TList.h"
#include "TObjString.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "THttpServer.h"
#include "TRootSniffer.h"
#include "TDatime.h"
#include "TClass.h"

Bool_t bRun = kTRUE;

class TMsgList : public TNamed {

   protected:

      TList      fMsgs;       //  list messages, stored as TObjString
      Int_t      fLimit;      //  max number of stored messages
      Long64_t   fCounter;    //  current message id
      TList      fSelect;     //! temporary list used for selection
      TObjString fStrCounter; //! current id stored in the string

   public:

      TMsgList(const char* name = "log", Int_t limit = 1000) :
         TNamed(name,"list of log messages"),
         fMsgs(),
         fLimit(limit),
         fCounter(0),
         fSelect(),
         fStrCounter()
      {
         fMsgs.SetOwner(kTRUE);

         // counter initialized from current time
         // if application restarted, id will be bigger and request from browser
         // will not lead to messages lost. Of course, if more than 1000 messages
         // per second are generated, one could have mismatch

         fCounter = ((Long64_t) TDatime().Get()) * 1000;
      }

      virtual ~TMsgList() { fMsgs.Clear(); }

      void AddMsg(const char* msg)
      {
         // add message to the list
         // if number of stored messages bigger than configured, old messages will be removed
         // zero (msg==0) messages will not be add to the list

         while (fMsgs.GetSize() >= fLimit) {
            TObject* last = fMsgs.Last();
            fMsgs.RemoveLast();
            delete last;
         }
         if (msg==0) return;

         fMsgs.AddFirst(new TObjString(msg));
         fCounter++;
      }

      TList* Select(Int_t max = 0, Long64_t id = 0)
      {
         // Central method to select new messages
         // Current id stored as first item and used on the client to request new portion
         // One could limit number of returned messages

         TIter iter(&fMsgs);
         TObject* obj = 0;
         Long64_t curr = fCounter;
         fSelect.Clear();

         if (max == 0) max = fMsgs.GetLast()+1;

         // add current id as first string in the list
         fStrCounter.SetString(TString::LLtoa(fCounter, 10));
         fSelect.Add(&fStrCounter);

         while (((obj = iter()) != 0) && (--curr >= id) && (--max>=0)) fSelect.Add(obj);

         return &fSelect;
      }

   ClassDef(TMsgList, 1); // Custom messages list
};

void httptextlog()
{
   // create logging instance
   TMsgList* log = new TMsgList("log", 200);

   if ((TMsgList::Class()->GetMethodAllAny("Select") == 0) || (strcmp(log->ClassName(), "TMsgList")!=0)) {
      printf("Most probably, macro runs in interpreter mode\n");
      printf("To access new methods from TMsgList class,\n");
      printf("one should run macro with ACLiC like:\n");
      printf("   shell> root -b httpextlog.C+\n");
      return;
    }

   if (gSystem->AccessPathName("httptextlog.js")!=0) {
      printf("Please start macro from directory where httptextlog.js is available\n");
      printf("Only in this case web interface can work\n");
      return;
   }

   // create histograms, just for fun
   TH1D *hpx = new TH1D("hpx","This is the px distribution",100,-4,4);
   hpx->SetFillColor(48);
   hpx->SetDirectory(0);
   TH2F *hpxpy = new TH2F("hpxpy","py vs px",40,-4,4,40,-4,4);
   hpxpy->SetDirectory(0);

   // start http server
   THttpServer* serv = new THttpServer("http:8080");

   // One could specify location of newer version of JSROOT
   // serv->SetJSROOT("https://root.cern.ch/js/latest/");
   // serv->SetJSROOT("http://jsroot.gsi.de/latest/");

   // let always load httptextlog.js script in the browser
   serv->GetSniffer()->SetAutoLoad("currentdir/httptextlog.js");

   // register histograms
   serv->Register("/", hpx);
   serv->Register("/", hpxpy);

   // register log instance
   serv->Register("/", log);

   // while server runs in read-only mode, we should allow methods execution
   serv->Restrict("/log", "allow_method=Select,GetTitle");

   // register exit command
   serv->RegisterCommand("/Stop","bRun=kFALSE;", "rootsys/icons/ed_delete.png");
   serv->RegisterCommand("/ExitRoot","gSystem->Exit(1);", "rootsys/icons/ed_delete.png");

   // Fill histograms randomly
   TRandom3 random;
   Float_t px, py;
   const Long_t kUPDATE = 1000;
   Long_t cnt = 0;
   while (bRun) {
      random.Rannor(px,py);
      hpx->Fill(px);
      hpxpy->Fill(px,py);

      // IMPORTANT: one should regularly call ProcessEvents
      if (cnt++ % kUPDATE == 0) {
         if (gSystem->ProcessEvents()) break;

         Long_t loop = cnt / kUPDATE;

         // make messages not very often
         if (loop % 1000 == 0) {
            loop = loop/1000;
            int shift = loop % 40;
            // make a 'stairs' with spaces
            log->AddMsg(TString::Format("%*s Message %ld", shift, "", loop));
         }
      }
   }

   delete serv; // delete http server
}
