// Author: Sergey Linev   7/12/2016

/*************************************************************************
 * Copyright (C) 2016, Sergey Linev                                      *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebCanvas.h"

#include "TWebSnapshot.h"
#include "TWebPadPainter.h"
#include "TWebVirtualX.h"
#include "TWebMenuItem.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TClass.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TList.h"
#include "TH1.h"
#include "TGraph.h"
#include "TBufferJSON.h"
#include "Riostream.h"

#include <ROOT/TWebWindowsManager.hxx>

#include <stdio.h>
#include <string.h>


ClassImp(TWebCanvas);

TWebCanvas::TWebCanvas() :
   TCanvasImp(),
   fWebConn(),
   fHasSpecials(kFALSE),
   fCanvVersion(1)
{
}

TWebCanvas::TWebCanvas(TCanvas *c, const char *name, Int_t x, Int_t y, UInt_t width, UInt_t height) :
   TCanvasImp(c, name, x, y, width, height),
   fWebConn(),
   fHasSpecials(kFALSE),
   fCanvVersion(1)
{
}

TWebCanvas::~TWebCanvas()
{
}

Int_t TWebCanvas::InitWindow()
{
   TWebVirtualX *vx = dynamic_cast<TWebVirtualX *> (gVirtualX);
   if (vx) vx->SetWebCanvasSize(Canvas()->GetWw(), Canvas()->GetWh());

   // at this place canvas is not yet register to the list of canvases - we cannot start browser
   return 777111777; // magic number, should be catch by TWebVirtualX
}

TVirtualPadPainter* TWebCanvas::CreatePadPainter()
{
   return new TWebPadPainter();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE when object is fully supported on JSROOT side
/// In ROOT7 Paint function will just return appropriate flag that object can be displayed on JSROOT side

Bool_t TWebCanvas::IsJSSupportedClass(TObject* obj)
{
   if (!obj) return kTRUE;

   static const struct {
      const char *name;
      bool with_derived;
   } supported_classes[] = {
       { "TH1", true },
       { "TF1", true },
       { "TGraph", true },
       { "TFrame", false },
       { "THStack", false },
       { "TMultiGraph", false },
       { "TGraphPolargram", true },
       { "TPave", true },
       { "TGaxis", false },
       { "TPave", true },
       { "TArrow", false },
       { "TBox", false },      // in principle, can be handled via TWebPainter
       { "TWbox", false },     // some extra calls which cannout be handled via TWebPainter
       { "TLine", false },     // also can be handler via TWebPainter
       { "TText", false },
       { "TLatex", false },
       { "TMathText", false },
       { "TPolyMarker3D", false },
       { "TGraph2D", false },
       { 0, false }
   };

   // fast check of class name
   for (int i = 0; supported_classes[i].name != 0; ++i)
      if (strcmp(supported_classes[i].name, obj->ClassName()) == 0) return kTRUE;

   // now check inheritance only for configured classes
   for (int i = 0; supported_classes[i].name != 0; ++i)
      if (supported_classes[i].with_derived)
         if (obj->InheritsFrom(supported_classes[i].name)) return kTRUE;

   printf("Unsupported class %s\n", obj->ClassName());

   return kFALSE;
}

TObject* TWebCanvas::FindPrimitive(const char *sid, TPad *pad)
{
   // search of object with given id in list of primitives

   if (!pad) pad = Canvas();

   const char *kind = "";
   const char *separ = strchr(sid, '#');
   UInt_t id = 0;

   if (separ == 0) {
      id = (UInt_t) TString(sid).Atoll();
   } else {
      kind = separ + 1;
      id = (UInt_t) TString(sid, separ-sid).Atoll();
   }

   if (TString::Hash(&pad, sizeof(pad)) == id) return pad;

   TIter iter(pad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = iter()) != 0) {
      TH1 *h1 = obj->InheritsFrom(TH1::Class()) ? (TH1 *)obj : 0;
      if (TString::Hash(&obj, sizeof(obj)) == id) {
         if (h1 && (*kind == 'x')) return h1->GetXaxis();
         if (h1 && (*kind == 'y')) return h1->GetYaxis();
         if (h1 && (*kind == 'z')) return h1->GetZaxis();
         return obj;
      }
      if (h1 != 0) {
         TIter fiter(h1->GetListOfFunctions());
         TObject *fobj = 0;
         while ((fobj = fiter()) != 0)
            if (TString::Hash(&fobj, sizeof(fobj)) == id) return fobj;
      } else if (obj->InheritsFrom(TPad::Class())) {
         obj = FindPrimitive(sid, (TPad*) obj);
         if (obj) return obj;
      }
   }

   return 0;
}


TWebSnapshot *TWebCanvas::CreateObjectSnapshot(TObject *obj, const char *opt)
{
   TWebSnapshot *sub = new TWebSnapshot();
   sub->SetObjectIDAsPtr(obj);
   sub->SetOption(opt);
   TWebPainting *p = 0;

   if (!IsJSSupportedClass(obj)) {
      TWebPadPainter *painter = dynamic_cast<TWebPadPainter *> (Canvas()->GetCanvasPainter());
      if (painter) {
         painter->ResetPainting(); // ensure painter is created
         painter->SetWebCanvasSize(Canvas()->GetWw(), Canvas()->GetWh()); // provide canvas dimension
      }

      TWebVirtualX *vx = dynamic_cast<TWebVirtualX *> (gVirtualX);
      if (vx) {
         vx->SetWebCanvasSize(Canvas()->GetWw(), Canvas()->GetWh());
         vx->SetWebPainter(painter); // redirect virtualx back to pad painter
      }

      // calling Paint function for the object
      obj->Paint(opt);

      if (vx) vx->SetWebPainter(0);

      if (painter) p = painter->TakePainting();
      fHasSpecials = kTRUE;
   }

   // when paint method was used and resultative,

   if (p) {
      p->FixSize();
      sub->SetSnapshot(TWebSnapshot::kSVG, p, kTRUE);
   } else {
      sub->SetSnapshot(TWebSnapshot::kObject, obj);
   }

   return sub;
}

Bool_t TWebCanvas::AddCanvasSpecials(TPadWebSnapshot *master)
{
   // if (!TColor::DefinedColors()) return 0;
   TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();

   if (!colors) return kFALSE;
   Int_t cnt = 0;
   for (Int_t n=0;n<=colors->GetLast();++n)
      if (colors->At(n) != 0) cnt++;
   if (cnt <= 598) return kFALSE; // normally there are 598 colors defined

   TWebSnapshot *sub = new TWebSnapshot();
   sub->SetSnapshot(TWebSnapshot::kSpecial, colors);
   master->Add(sub);

   if (gDebug>1) Info("AddCanvasSpecials" ,"ADD COLORS TABLES %d", cnt);

   //save the current palette
   TArrayI pal = TColor::GetPalette();
   Int_t palsize = pal.GetSize();
   TObjArray *CurrentColorPalette = new TObjArray();
   CurrentColorPalette->SetName("CurrentColorPalette");
   for (Int_t i=0; i<palsize; i++) CurrentColorPalette->Add(gROOT->GetColor(pal[i]));

   sub = new TWebSnapshot();
   sub->SetSnapshot(TWebSnapshot::kSpecial, CurrentColorPalette, kTRUE);
   master->Add(sub);

   return kTRUE;
}


TString TWebCanvas::CreateSnapshot(TPad* pad, TPadWebSnapshot *master, TList *primitives_lst)
{
   TList master_lst; // main list of TList object which are primitives or functions
   if (!master && !primitives_lst) primitives_lst = &master_lst;

   TPadWebSnapshot *curr = new TPadWebSnapshot();
   if (master) {
      curr->SetObjectIDAsPtr(pad);
      master->Add(curr);
   }

   TWebSnapshot* padshot = new TWebSnapshot();
   padshot->SetObjectIDAsPtr(pad);
   padshot->SetSnapshot(TWebSnapshot::kObject, pad);
   curr->Add(padshot);

   if (primitives_lst == &master_lst) AddCanvasSpecials(curr);

   TList *primitives = pad->GetListOfPrimitives();

   primitives_lst->Add(primitives); // add list of primitives

   TIter iter(primitives);
   TObject* obj = 0;
   while ((obj = iter()) != 0) {
      if (obj->InheritsFrom(TPad::Class())) {
         CreateSnapshot((TPad*) obj, curr, primitives_lst);
      } else if (obj->InheritsFrom(TH1::Class())) {
         TWebSnapshot *sub = new TWebSnapshot();
         TH1 *hist = (TH1 *) obj;
         sub->SetObjectIDAsPtr(hist);
         sub->SetOption(iter.GetOption());
         sub->SetSnapshot(TWebSnapshot::kObject, obj);
         curr->Add(sub);

         TIter fiter(hist->GetListOfFunctions());
         TObject *fobj = 0;
         while ((fobj = fiter()) != 0)
            if (!fobj->InheritsFrom("TPaveStats") && !fobj->InheritsFrom("TPaletteAxis"))
               curr->Add(CreateObjectSnapshot(fobj, fiter.GetOption()));

         primitives_lst->Add(hist->GetListOfFunctions());
      } else if (obj->InheritsFrom(TGraph::Class())) {
         TWebSnapshot *sub = new TWebSnapshot();
         TGraph *gr = (TGraph *) obj;
         sub->SetObjectIDAsPtr(gr);
         sub->SetOption(iter.GetOption());
         sub->SetSnapshot(TWebSnapshot::kObject, obj);
         curr->Add(sub);

         TIter fiter(gr->GetListOfFunctions());
         TObject *fobj = 0;
         while ((fobj = fiter()) != 0)
            if (!fobj->InheritsFrom("TPaveStats")) // stats should be created on the client side
               curr->Add(CreateObjectSnapshot(fobj, fiter.GetOption()));

         primitives_lst->Add(gr->GetListOfFunctions());
      } else {
         curr->Add(CreateObjectSnapshot(obj, iter.GetOption()));
      }
   }

   if (primitives_lst != &master_lst) return "";

   // now move all primitives and functions into separate list to perform I/O

   // TBufferJSON::ExportToFile("canvas.json", pad);

   TList save_lst;
   TIter diter(&master_lst);
   TList *dlst = 0;
   while ((dlst = (TList*) diter()) != 0) {
      TIter fiter(dlst);
      while ((obj = fiter()) != 0)
         save_lst.Add(obj, fiter.GetOption());
      save_lst.Add(dlst); // add list itslef to have marker
      dlst->Clear("nodelete");
   }

   // TBufferJSON::ExportToFile("canvas_empty.json", pad);

   //gDebug = 4;

   // Info("CreateSnapshot","In canvas primitives are %d", pad->GetListOfPrimitives()->GetSize());

   TString res = TBufferJSON::ConvertToJSON(curr, 23);
   // gDebug = 0;

   // TODO: this is only for debugging, remove it later
   TBufferJSON::ExportToFile("snapshot.json", curr);

   delete curr; // destroy created snapshot

   TIter siter(&save_lst);
   diter.Reset();
   while ((dlst = (TList*) diter()) != 0) {
      while ((obj = siter()) != 0) {
         if (obj == dlst) break;
         dlst->Add(obj, siter.GetOption());
      }
   }

   save_lst.Clear("nodelete");

   master_lst.Clear("nodelete");

   return res;
}


void TWebCanvas::CheckDataToSend()
{
   if (!Canvas()) return;

   for (WebConnList::iterator citer = fWebConn.begin(); citer != fWebConn.end(); ++citer) {
      WebConn& conn = *citer;

      // check if direct data sending is possible
      if (!fWindow->CanSend(conn.fConnId, true))
         continue;

      TString buf;

      if (conn.fGetMenu.Length()>0) {

         TObject* obj = FindPrimitive(conn.fGetMenu.Data());
         if (!obj) obj = Canvas();

         TWebMenuItems items;
         items.PopulateObjectMenu(obj, obj->IsA());
         buf = "MENU:";
         buf.Append(conn.fGetMenu);
         buf.Append(":");
         buf += items.ProduceJSON();

         conn.fGetMenu.Clear();
      } else if (conn.fDrawVersion < fCanvVersion) {
         buf = "SNAP6:";
         buf.Append(TString::LLtoa(fCanvVersion, 10));
         buf.Append(":");
         buf += CreateSnapshot(Canvas());

         // printf("Snapshot created %d\n", buf.Length());
         //if (buf.Length() < 10000) printf("Snapshot %s\n", buf.Data());
      } else if (conn.fSend.Length() > 0) {
         buf = conn.fSend;
         conn.fSend.Clear();
      }

      if (buf.Length() > 0) {
         // sending of data can be moved into separate thread - not to block user code
         fWindow->Send(buf.Data(), conn.fConnId);
      }
   }
}

void TWebCanvas::Close()
{
   printf("Call TWebCanvas::Close\n");
}

void TWebCanvas::Show()
{
   const char *swhere = gSystem->Getenv("WEBGUI_WHERE"); // let configure place like with ROOT7
   std::string where = swhere ? swhere : "browser";

   if (!fWindow) {
      fWindow = ROOT::Experimental::TWebWindowsManager::Instance()->CreateWindow(gROOT->IsBatch());

      fWindow->SetConnLimit(0); // allow any number of connections

      fWindow->SetDefaultPage("file:$jsrootsys/files/canvas6.htm");

      fWindow->SetDataCallBack([this](unsigned connid, const std::string &arg) { ProcessData(connid, arg); });

      // fWindow->SetGeometry(500,300);
   }

   fWindow->Show(where);
}

void TWebCanvas::ShowCmd(const char *arg, Bool_t show)
{
   // command used to toggle showing of menu, toolbar, editors, ...
   for (WebConnList::iterator citer = fWebConn.begin(); citer != fWebConn.end(); ++citer) {
      WebConn& conn = *citer;

      if (!conn.fConnId) continue;

      conn.fSend = "SHOW:";
      conn.fSend.Append(arg);
      conn.fSend.Append(show ? ":1" : ":0");
   }

   CheckDataToSend();
}


Bool_t TWebCanvas::DecodePadRanges(TPad *pad, const char *arg)
{
   if (!pad || !arg || !*arg) return kFALSE;

   Double_t ux1,ux2,uy1,uy2,px1,px2,py1,py2;
   Int_t cnt = sscanf(arg, "%lf:%lf:%lf:%lf:%lf:%lf:%lf:%lf",&ux1,&ux2,&uy1,&uy2,&px1,&px2,&py1,&py2);
   if (cnt!=8) return kFALSE;

   Double_t ux1_,ux2_,uy1_,uy2_,px1_,px2_,py1_,py2_;

   pad->GetRange(px1_,py1_,px2_,py2_);
   pad->GetRangeAxis(ux1_,uy1_,ux2_,uy2_);

   if ((ux1==ux1_) && (ux2==ux2_) && (uy1==uy1_) && (uy2==uy2_) && (px1==px1_) && (px2==px2_) && (py1==py1_) && (py2==py2_)) {
      //Info("DecodePadRanges","Ranges not changed");
      return kFALSE;
   }

   pad->Range(px1,py1,px2,py2);
   pad->RangeAxis(ux1,uy1,ux2,uy2);

   if (gDebug > 0)
      Info("DecodePadRanges", "Apply new ranges %s for pad %s", arg, pad->GetName());

   // without special objects no need for explicit update of the canvas
   if (!fHasSpecials) return kFALSE;

   pad->Modified(kTRUE);
   return kTRUE;
}

Bool_t TWebCanvas::DecodeAllRanges(const char *arg)
{
   if (!arg || !*arg) return kFALSE;
   //Bool_t isany = kFALSE;

   const char *curr = arg, *pos = 0;

   while ((pos = strstr(curr, "id=")) != 0) {
      curr = pos + 3;
      const char *next = strstr(curr,":");
      if (!next) break;

      TString sid(curr, next-curr);
      TPad *pad = dynamic_cast<TPad *>(FindPrimitive(sid.Data()));

      curr = next+1;
      DecodePadRanges(pad, curr);
   }

   return kTRUE;
}

void TWebCanvas::ProcessData(unsigned connid, const std::string &arg)
{
   if (arg.empty()) return;

   if (arg == "CONN_READY") {

      WebConn newconn;
      newconn.fConnId = connid;

      fWebConn.push_back(newconn);

      CheckDataToSend();

      return;
   }

   // try to identify connection for given WS request
   WebConn* conn(nullptr);
   WebConnList::iterator iter = fWebConn.begin();
   while (iter != fWebConn.end()) {
      if (iter->fConnId == connid) {
         conn = &(*iter); break;
      }
      ++iter;
   }

   if (!conn) {
      printf("Get data without not existing connection %u\n", connid);
      return;
   }

   const char *cdata = arg.c_str();

   if (arg == "CONN_CLOSED") {
      fWebConn.erase(iter);
   } else if (strncmp(cdata,"READY",5)==0) {
      CheckDataToSend();
   } else if (strncmp(cdata, "RREADY:", 7)==0) {
      cdata += 7;

      const char *separ = strchr(cdata, ':');
      if (!separ) {
         conn->fDrawVersion = TString(cdata).Atoll();
      } else {
         conn->fDrawVersion = TString(cdata, separ-cdata).Atoll();
         cdata = separ+1;
         if (gDebug>1) Info("ProcessData", "RANGES %s", cdata);
         if (connid==0) DecodeAllRanges(cdata); // only first connection get ranges
      }
      CheckDataToSend();
   } else if (strncmp(cdata,"GETMENU:",8)==0) {
      conn->fGetMenu = cdata+8;
      CheckDataToSend();
   } else if (strncmp(cdata,"OBJEXEC:",8)==0) {
      TString buf(cdata+8);
      Int_t pos = buf.First(':');

      if (pos>0) {
         TString sid(buf, pos);
         buf.Remove(0, pos+1);

         TObject* obj = FindPrimitive(sid.Data());
         if (obj && (buf.Length()>0)) {
            TString exec;
            exec.Form("((%s*) %p)->%s;", obj->ClassName(), obj, buf.Data());
            Info("ProcessWS", "Obj %s Execute %s", obj->GetName(), exec.Data());
            gROOT->ProcessLine(exec);

            // PerformUpdate(); // check that canvas was changed
            if (IsAnyPadModified(Canvas())) fCanvVersion++;
            CheckDataToSend();
         }
      }
   } else if (strncmp(cdata,"EXECANDSEND:",12)==0) {
      TString buf(cdata+12), reply;
      TObject *obj = 0;

      Int_t pos = buf.First(':');

      if (pos>0) {
         reply.Append(buf, pos);
         buf.Remove(0, pos+1);
         pos = buf.First(':');
         if (pos>0) {
            TString sid(buf, pos);
            buf.Remove(0, pos+1);
            obj = FindPrimitive(sid.Data());
         }
      }

      if (obj && (buf.Length()>0) && (reply.Length()>0)) {
         TString exec;
         exec.Form("((%s*) %p)->%s;", obj->ClassName(), obj, buf.Data());
         if (gDebug > 1) Info("ProcessData", "Obj %s Exec %s", obj->GetName(), exec.Data());

         Long_t res = gROOT->ProcessLine(exec);
         TObject *resobj = (TObject *) res;
         if (resobj) {
            conn->fSend = reply;
            conn->fSend.Append(":");
            conn->fSend.Append(TBufferJSON::ConvertToJSON(resobj,23));
            if (reply[0]=='D') delete resobj; // delete object if first symbol in reply is D
         }

         CheckDataToSend(); // check if data should be send
      }
   } else if (strncmp(cdata,"QUIT",4)==0) {
      // use window manager to correctly terminate http server
      ROOT::Experimental::TWebWindowsManager::Instance()->Terminate();
   } else if (strncmp(cdata,"RELOAD",6)==0) {
      conn->fDrawVersion = 0;
      CheckDataToSend();
   } else if (strncmp(cdata,"GETIMG:",7)==0) {
      const char* img = cdata+7;

      const char* separ = strchr(img,':');
      if (separ) {
         TString filename(img, separ-img);
         img = separ+1;
         filename.Append(".svg"); // temporary - JSROOT returns SVG

         std::ofstream ofs(filename);
         ofs << "<?xml version=\"1.0\" standalone=\"no\"?>";
         ofs << img;
         ofs.close();

         Info("ProcessWS", "SVG file %s has been created", filename.Data());
      }
      CheckDataToSend();
   } else if (arg == "KEEPALIVE") {
   } else {
      Error("ProcessWS", "GET unknown request %d %30s", (int) arg.length(), cdata);
   }
}

Bool_t TWebCanvas::IsAnyPadModified(TPad *pad)
{
   // returns true when any pad or sub pad modified
   // reset modified flags

   Bool_t res = kFALSE;

   if (pad->IsModified()) {
      pad->Modified(kFALSE);
      res = kTRUE;
   }

   TIter iter(pad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = iter()) != 0) {
      if (obj->InheritsFrom(TPad::Class()) && IsAnyPadModified((TPad*) obj)) res = kTRUE;
   }

   return res;
}

UInt_t TWebCanvas::GetWindowGeometry(Int_t &x, Int_t &y, UInt_t &w, UInt_t &h)
{
   // reset dimension in gVirtualX  - it will be requested immediately
   TWebVirtualX *vx = dynamic_cast<TWebVirtualX *> (gVirtualX);
   if (vx) vx->SetWebCanvasSize(Canvas()->GetWw(), Canvas()->GetWh());

   x = 0; y = 0;
   w = Canvas()->GetWw() + 4;
   h = Canvas()->GetWh() + 28;
   return 0;
}

Bool_t TWebCanvas::PerformUpdate()
{
   // check if canvas modified. If true and communication allowed,
   // It scan all primitives in the TCanvas and subpads and convert them into
   // the structure which will be delivered to JSROOT client

   if (IsAnyPadModified(Canvas())) fCanvVersion++;

   CheckDataToSend();

   // block in canvas update, can it be optional
   WaitWhenCanvasPainted(fCanvVersion);

   return kTRUE;
}


Bool_t TWebCanvas::WaitWhenCanvasPainted(Long64_t ver)
{
   // simple polling loop until specified version delivered to the clients

   long cnt = 0;
   bool had_connection = false;

   if (gDebug>2) Info("WaitWhenCanvasPainted", "version %ld", (long) ver);

   while (cnt++ < 1000) {
      if (fWebConn.size() > 0) had_connection = true;

      if ((fWebConn.size() == 0) && (had_connection || (cnt > 800))) {
         if (gDebug>2) Info("WaitWhenCanvasPainted", "no connections - abort");
         return kFALSE; // wait ~1 min if no new connection established
      }

      if ((fWebConn.size() > 0) && (fWebConn.front().fDrawVersion >= ver)) {
         if (gDebug>2) Info("WaitWhenCanvasPainted", "ver %ld got painted", (long) ver);
         return kTRUE;
      }

      gSystem->ProcessEvents();

      gSystem->Sleep((cnt < 500) ? 1 : 100); // increase sleep interval when do very often
   }

   if (gDebug>2) Info("WaitWhenCanvasPainted", "timeout");

   return kFALSE;
}

