// Author: Sergey Linev   7/12/2016

/*************************************************************************
 * Copyright (C) 2016, Sergey Linev                                      *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebCanvas.h"

#include "THttpCallArg.h"
#include "THttpEngine.h"
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
#include "TBufferJSON.h"
#include "TApplication.h"
#include "Riostream.h"

#include <stdio.h>
#include <string.h>


ClassImp(TWebCanvas)

TWebCanvas::TWebCanvas() :
   THttpWSHandler("", ""),
   TCanvasImp(),
   fWebConn(),
   fAddr(),
   fServer(0),
   fHasSpecials(kFALSE)
{
}

TWebCanvas::TWebCanvas(TCanvas *c, const char *name, Int_t x, Int_t y, UInt_t width, UInt_t height, TString addr, void *serv) :
   THttpWSHandler(name, "title"),
   TCanvasImp(c, name, x, y, width, height),
   fWebConn(),
   fAddr(addr),
   fServer(serv),
   fHasSpecials(kFALSE)
{
   UInt_t hash = TString::Hash(&c, sizeof(c));
   SetName(Form("0x%u", hash)); // make name very screwed
}

TWebCanvas::~TWebCanvas()
{
   for (WebConnList::iterator iter = fWebConn.begin(); iter != fWebConn.end(); ++iter) {
      if (iter->fHandle) {
         iter->fHandle->ClearHandle();
         delete iter->fHandle;
         iter->fHandle = 0;
      }
   }
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

   if (obj->InheritsFrom("TH1") || obj->InheritsFrom("TGraph") || obj->InheritsFrom("TF1") ||
       obj->InheritsFrom("TFrame") || obj->InheritsFrom("TBox") ||
       obj->InheritsFrom("TPave") || obj->InheritsFrom("TArrow")) return kTRUE;

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
      if (painter) painter->ResetPainting(); // ensure painter is created

      // calling Paint function for the object
      obj->Paint(opt);

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

   printf("ADD COLORS TABLES %d\n", cnt);

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


TString TWebCanvas::CreateSnapshot(TPad* pad, TPadWebSnapshot *master, TList *tempbuf)
{
   TList main_buf;
   if (!master && !tempbuf) tempbuf = &main_buf;

   TPadWebSnapshot *curr = new TPadWebSnapshot();
   if (master) {
      curr->SetObjectIDAsPtr(pad);
      master->Add(curr);
   }

   TWebSnapshot* padshot = new TWebSnapshot();
   padshot->SetObjectIDAsPtr(pad);
   padshot->SetSnapshot(TWebSnapshot::kObject, pad);
   curr->Add(padshot);

   if (tempbuf == &main_buf) AddCanvasSpecials(curr);

   TList *primitives = pad->GetListOfPrimitives();
   TList hlist; // list of histograms, required for functions handling

   TIter iter(primitives);
   TObject* obj = 0;
   while ((obj = iter()) != 0) {
      if (obj->InheritsFrom(TPad::Class())) {
         CreateSnapshot((TPad*) obj, curr, tempbuf);
      } else if (obj->InheritsFrom(TH1::Class())) {
         TWebSnapshot *sub = new TWebSnapshot();
         TH1 *hist = (TH1*) obj;
         sub->SetObjectIDAsPtr(hist);
         sub->SetOption(iter.GetOption());
         sub->SetSnapshot(TWebSnapshot::kObject, obj);
         curr->Add(sub);

         TIter fiter(hist->GetListOfFunctions());
         TObject *fobj = 0;
         while ((fobj = fiter()) != 0)
            if (!fobj->InheritsFrom("TPaveStats") && !fobj->InheritsFrom("TPaletteAxis"))
               curr->Add(CreateObjectSnapshot(fobj, fiter.GetOption()));

         hlist.Add(hist);
      } else {
         curr->Add(CreateObjectSnapshot(obj, iter.GetOption()));
      }
   }

   const char *pad_marker = "!!!pad!!!";
   const char *hist_marker = "!!!hist!!!";

   // remove primitives and keep them in extra list
   tempbuf->Add(pad, pad_marker); // special marker for pad
   iter.Reset();
   while ((obj = iter()) != 0)
      tempbuf->Add(obj, iter.GetOption());
   primitives->Clear("nodelete");

   // remove functions from all histograms and also add them to extra list
   TIter hiter(&hlist);
   TH1 *h1 = 0;
   while ((h1 = (TH1*) hiter()) != 0) {
      tempbuf->Add(h1, hist_marker); // special marker for pad
      TIter fiter(h1->GetListOfFunctions());
      while ((obj = fiter()) != 0)
         tempbuf->Add(obj, fiter.GetOption());
      h1->GetListOfFunctions()->Clear("nodelete");
   }

   hlist.Clear("nodelete");

   if (tempbuf != &main_buf) return "";

   TString res = TBufferJSON::ConvertToJSON(curr, 23);

   // TBufferJSON::ExportToFile("debug.json", curr);

   delete curr; // destroy created snapshot

   TPad *rpad = 0;
   h1 = 0;
   TIter miter(&main_buf);
   while ((obj = miter()) != 0) {
      TString opt = miter.GetOption();

      if (opt == pad_marker) {
         rpad = (TPad*) obj; h1 = 0;
      } else if (opt == hist_marker) {
         rpad = 0; h1 = (TH1*) obj;
      } else if (rpad != 0) {
         rpad->GetListOfPrimitives()->Add(obj, opt);
      } else if (h1 != 0) {
         h1->GetListOfFunctions()->Add(obj, opt);
      }
   }

   main_buf.Clear("nodelete");

   return res;
}


void TWebCanvas::CheckModifiedFlag()
{
   if (!Canvas()) return;

   for (WebConnList::iterator citer = fWebConn.begin(); citer != fWebConn.end(); ++citer) {
      WebConn& conn = *citer;

      if (!conn.fReady || !conn.fHandle) continue;

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
      } else if (conn.fModified) {
         buf = "SNAP6:";
         buf += CreateSnapshot(Canvas());

         // printf("Snapshot created %d\n", buf.Length());
         //if (buf.Length() < 10000) printf("Snapshot %s\n", buf.Data());
         conn.fModified = kFALSE;
      } else if (conn.fSend.Length() > 0) {
         buf = conn.fSend;
         conn.fSend.Clear();
      }

      if (buf.Length() > 0) {
         // sending of data can be moved into separate thread - not to block user code
         conn.fReady = kFALSE;
         conn.fHandle->SendCharStar(buf.Data());
      }
   }
}

void TWebCanvas::Close()
{
   printf("Call TWebCanvas::Close\n");
}

void TWebCanvas::Show()
{
   TString addr;

   Func_t symbol_qt5 = gSystem->DynFindSymbol("*", "webgui_start_browser_in_qt5");
   if (symbol_qt5) {
      typedef void (*FunctionQt5)(const char *, void *, bool);

      addr.Form("://dummy:8080/web6gui/%s/draw.htm?longpollcanvas&no_root_json%s&qt5", GetName(),
                (gROOT->IsBatch() ? "&batch_mode" : ""));
      // addr.Form("example://localhost:8080/Canvases/%s/draw.htm", Canvas()->GetName());

      Info("NewDisplay", "Show canvas in Qt5 window:  %s", addr.Data());

      FunctionQt5 func = (FunctionQt5)symbol_qt5;
      func(addr.Data(), fServer, gROOT->IsBatch());
      return;
   }

   // TODO: one should try to load CEF libraries only when really needed
   // probably, one should create separate DLL with CEF-related code
   Func_t symbol_cef = gSystem->DynFindSymbol("*", "webgui_start_browser_in_cef3");
   const char *cef_path = gSystem->Getenv("CEF_PATH");
   const char *rootsys = gSystem->Getenv("ROOTSYS");
   if (symbol_cef && cef_path && !gSystem->AccessPathName(cef_path) && rootsys) {
      typedef void (*FunctionCef3)(const char *, void *, bool, const char *, const char *);

      addr.Form("/web6gui/%s/draw.htm?cef_canvas&no_root_json%s", GetName(), (gROOT->IsBatch() ? "&batch_mode" : ""));

      Info("NewDisplay", "Show canvas in CEF window:  %s", addr.Data());

      FunctionCef3 func = (FunctionCef3)symbol_cef;
      func(addr.Data(), fServer, gROOT->IsBatch(), rootsys, cef_path);

      return;
   }

   addr.Form("%s/web6gui/%s/draw.htm?webcanvas", fAddr.Data(), GetName());

   Info("Show", "Call TWebCanvas::Show:  %s", addr.Data());

   // exec.Form("setsid chromium --app=http://localhost:8080/Canvases/%s/draw.htm?websocket </dev/null >/dev/null 2>/dev/null &", Canvas()->GetName());


   const char *swhere = gSystem->Getenv("WEBGUI_WHERE"); // let configure place like with ROOT7
   std::string where = swhere ? swhere : "browser";

   TString exec;

   if (where != "browser") {
      if (where.find("$url") != std::string::npos) {
         exec = where.c_str();
         exec.ReplaceAll("$url", addr);
      } else {
         exec.Form("%s %s", where.c_str(), addr.Data());
      }
   } else
   if (gSystem->InheritsFrom("TMacOSXSystem"))
      exec.Form("open %s", addr.Data());
   else
      exec.Form("xdg-open %s &", addr.Data());
   printf("Exec %s\n", exec.Data());

   gSystem->Exec(exec);
}

void TWebCanvas::ShowCmd(const char *arg, Bool_t show)
{
   // command used to toggle showing of menu, toolbar, editors, ...
   for (WebConnList::iterator citer = fWebConn.begin(); citer != fWebConn.end(); ++citer) {
      WebConn& conn = *citer;

      if (!conn.fHandle) continue;

      conn.fSend = "SHOW:";
      conn.fSend.Append(arg);
      conn.fSend.Append(show ? ":1" : ":0");
   }

   CheckModifiedFlag();
}



Bool_t TWebCanvas::DecodePadRanges(TPad *pad, const char *arg)
{
   if (!arg || !*arg) return kFALSE;

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
  Bool_t isany = kFALSE;

  const char *curr = arg, *pos = 0;

  while ((pos = strstr(curr, "id=")) != 0) {
     curr = pos + 3;
     const char *next = strstr(curr,":");
     if (!next) break;

     TString sid(curr, next-curr);
     TPad *pad = dynamic_cast<TPad *>(FindPrimitive(sid.Data()));

     curr = next+1;
     if (pad && DecodePadRanges(pad, curr)) isany = kTRUE;
  }

  if (isany) PerformUpdate();
  return kTRUE;
}


Bool_t TWebCanvas::ProcessWS(THttpCallArg *arg)
{
   if (!arg) return kTRUE;

   // try to identify connection for given WS request
   WebConn* conn = 0;
   WebConnList::iterator iter = fWebConn.begin();
   while (iter != fWebConn.end()) {
      if (iter->fHandle && (iter->fHandle->GetId() == arg->GetWSId()) && arg->GetWSId()) {
         conn = &(*iter); break;
      }
      ++iter;
   }

   if (strcmp(arg->GetMethod(),"WS_CONNECT")==0) {

      // accept all requests, in future one could limit number of connections
      // arg->Set404(); // refuse connection

   } else
   if (strcmp(arg->GetMethod(),"WS_READY")==0) {
      THttpWSEngine* wshandle = dynamic_cast<THttpWSEngine*> (arg->TakeWSHandle());

      if (conn != 0) Error("ProcessWS","WSHandle with given websocket id exists!!!");

      WebConn newconn;
      newconn.fHandle = wshandle;
      newconn.fModified = kTRUE;

      fWebConn.push_back(newconn);

      // if (gDebug>0) Info("ProcessRequest","Set WebSocket handle %p", wshandle);

      // connection is established
   } else
   if (strcmp(arg->GetMethod(),"WS_DATA")==0) {
      // process received data

      if (!conn) {
         Error("ProcessWS","Get websocket data without valid connection - ignore!!!");
         return kFALSE;
      }

      if (conn->fHandle->PreviewData(arg)) return kTRUE;

      const char* cdata = (arg->GetPostDataLength()<=0) ? "" : (const char*) arg->GetPostData();

      if (strncmp(cdata,"READY",5)==0) {
         conn->fReady = kTRUE;
         CheckModifiedFlag();
      } else
      if (strncmp(cdata, "RREADY:", 7)==0) {
         conn->fReady = kTRUE;
         if (!DecodeAllRanges(cdata+7)) CheckModifiedFlag();
      } else
      if (strncmp(cdata,"GETMENU:",8)==0) {
         conn->fReady = kTRUE;
         conn->fGetMenu = cdata+8;
         CheckModifiedFlag();
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

               PerformUpdate(); // check that canvas was changed
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
            if (gDebug > 1) Info("ProcessWS", "Obj %s Exec %s", obj->GetName(), exec.Data());

            Long_t res = gROOT->ProcessLine(exec);
            TObject *resobj = (TObject *) res;
            if (resobj) {
               conn->fSend = reply;
               conn->fSend.Append(":");
               conn->fSend.Append(TBufferJSON::ConvertToJSON(resobj,23));
               if (reply[0]=='D') delete resobj; // delete object if first symbol in reply is D
            }

            // PerformUpdate(); // check that canvas was changed

            CheckModifiedFlag(); // check if data should be send
         }
      } else if (strncmp(cdata,"QUIT",4)==0) {
         if (gApplication)
            gApplication->Terminate(0);
      } else if (strncmp(cdata,"RELOAD",6)==0) {
         conn->fModified = kTRUE;
         CheckModifiedFlag();
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
         conn->fReady = kTRUE;
         CheckModifiedFlag();
      } else if (strncmp(cdata,"KEEPALIVE",9)==0) {
      } else {
         Error("ProcessWS", "GET unknown request %d %s", (int) strlen(cdata), cdata);
      }

   } else
   if (strcmp(arg->GetMethod(),"WS_CLOSE")==0) {
      // connection is closed, one can remove handle

      if (conn && conn->fHandle) {
         conn->fHandle->ClearHandle();
         delete conn->fHandle;
         conn->fHandle = 0;
      }

      if (conn) fWebConn.erase(iter);
   }

   return kTRUE;
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

   if (IsAnyPadModified(Canvas())) {
      for (WebConnList::iterator iter = fWebConn.begin(); iter != fWebConn.end(); ++iter)
         iter->fModified = kTRUE;
   }

   CheckModifiedFlag();

   return kTRUE;
}
