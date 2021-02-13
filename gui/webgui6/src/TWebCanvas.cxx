// Author: Sergey Linev, GSI   7/12/2016

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebCanvas.h"

#include "TWebSnapshot.h"
#include "TWebPadPainter.h"
#include "TWebPS.h"
#include "TWebMenuItem.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TROOT.h"
#include "TClass.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TList.h"
#include "TH1.h"
#include "TEnv.h"
#include "TError.h"
#include "TGraph.h"
#include "TBufferJSON.h"
#include "TBase64.h"
#include "TAtt3D.h"
#include "TView.h"

#include <ROOT/RMakeUnique.hxx>

#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>

/** \class TWebCanvas
\ingroup webgui6

Basic TCanvasImp ABI implementation for Web-based GUI
Provides painting of main ROOT6 classes in web browsers
Major interactive features implemented in TWebCanvasFull class.

*/

using namespace std::string_literals;

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TWebCanvas::TWebCanvas(TCanvas *c, const char *name, Int_t x, Int_t y, UInt_t width, UInt_t height, Bool_t readonly)
   : TCanvasImp(c, name, x, y, width, height)
{
   fReadOnly = readonly;
   fStyleDelivery = gEnv->GetValue("WebGui.StyleDelivery", 0);
   fPaletteDelivery = gEnv->GetValue("WebGui.PaletteDelivery", 1);
   fPrimitivesMerge = gEnv->GetValue("WebGui.PrimitivesMerge", 100);
   fJsonComp = gEnv->GetValue("WebGui.JsonComp", TBufferJSON::kSameSuppression + TBufferJSON::kNoSpaces);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize window for the web canvas
/// At this place canvas is not yet register to the list of canvases - one cannot call RWebWindow::Show()

Int_t TWebCanvas::InitWindow()
{
   return 111222333; // should not be used at all
}

////////////////////////////////////////////////////////////////////////////////
/// Creates web-based pad painter

TVirtualPadPainter *TWebCanvas::CreatePadPainter()
{
   return new TWebPadPainter();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE when object is fully supported on JSROOT side
/// In ROOT7 Paint function will just return appropriate flag that object can be displayed on JSROOT side

Bool_t TWebCanvas::IsJSSupportedClass(TObject *obj)
{
   if (!obj)
      return kTRUE;

   static const struct {
      const char *name;
      bool with_derived;
   } supported_classes[] = {{"TH1", true},
                            {"TF1", true},
                            {"TGraph", true},
                            {"TFrame", false},
                            {"THStack", false},
                            {"TMultiGraph", false},
                            {"TGraphPolargram", true},
                            {"TPave", true},
                            {"TGaxis", false},
                            {"TPave", true},
                            {"TArrow", false},
                            {"TBox", false},  // in principle, can be handled via TWebPainter
                            {"TWbox", false}, // some extra calls which cannot be handled via TWebPainter
                            {"TLine", false}, // also can be handler via TWebPainter
                            {"TText", false},
                            {"TLatex", false},
                            {"TMathText", false},
                            {"TMarker", false},
                            {"TPolyMarker", false},
                            {"TPolyMarker3D", false},
                            {"TPolyLine3D", false},
                            {"TGraph2D", false},
                            {"TGraph2DErrors", false},
                            {"TASImage", false},
                            {nullptr, false}};

   // fast check of class name
   for (int i = 0; supported_classes[i].name != nullptr; ++i)
      if (strcmp(supported_classes[i].name, obj->ClassName()) == 0)
         return kTRUE;

   // now check inheritance only for configured classes
   for (int i = 0; supported_classes[i].name != nullptr; ++i)
      if (supported_classes[i].with_derived)
         if (obj->InheritsFrom(supported_classes[i].name))
            return kTRUE;

   return IsCustomClass(obj->IsA());
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Configures custom script for canvas.
/// If started from "load:" or "assert:" prefix will be loaded with JSROOT.AssertPrerequisites function
/// Script should implement custom user classes, which transferred as is to client
/// In the script draw handler for appropriate classes would be assigned

void TWebCanvas::SetCustomScripts(const std::string &src)
{
   fCustomScripts = src;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Assign custom class

void TWebCanvas::AddCustomClass(const std::string &clname, bool with_derived)
{
   if (with_derived)
      fCustomClasses.emplace_back("+"s + clname);
   else
      fCustomClasses.emplace_back(clname);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Checks if class belongs to custom

bool TWebCanvas::IsCustomClass(const TClass *cl) const
{
   for (auto &name : fCustomClasses) {
      if (name[0] == '+') {
         if (cl->InheritsFrom(name.substr(1).c_str()))
            return true;
      } else if (name.compare(cl->GetName()) == 0) {
         return true;
      }
   }
   return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates representation of the object for painting in web browser

void TWebCanvas::CreateObjectSnapshot(TPadWebSnapshot &master, TPad *pad, TObject *obj, const char *opt, TWebPS *masterps)
{
   if (IsJSSupportedClass(obj)) {
      master.NewPrimitive(obj, opt).SetSnapshot(TWebSnapshot::kObject, obj);
      return;
   }

   // painter is not necessary for batch canvas, but keep configuring it for a while
   auto *painter = dynamic_cast<TWebPadPainter *>(Canvas()->GetCanvasPainter());

   fHasSpecials = kTRUE;

   TView *view = nullptr;
   TVirtualPad *savepad = gPad;

   pad->cd();

   if (obj->InheritsFrom(TAtt3D::Class()) && !pad->GetView()) {
      pad->GetViewer3D("pad");
      view = TView::CreateView(1, 0, 0); // Cartesian view by default
      pad->SetView(view);

      // Set view to perform first auto-range (scaling) pass
      view->SetAutoRange(kTRUE);
   }

   TVirtualPS *saveps = gVirtualPS;

   TWebPS ps;
   gVirtualPS = masterps ? masterps : &ps;
   if (painter)
      painter->SetPainting(ps.GetPainting());

   // calling Paint function for the object
   obj->Paint(opt);

   if (view) {
      view->SetAutoRange(kFALSE);
      // call 3D paint once again to make real drawing
      obj->Paint(opt);
      pad->SetView(nullptr);
   }

   if (painter)
      painter->SetPainting(nullptr);

   gVirtualPS = saveps;
   if (savepad)
      savepad->cd();

   // if there are master PS, do not create separate entries
   if (!masterps && !ps.IsEmptyPainting())
      master.NewPrimitive(obj, opt).SetSnapshot(TWebSnapshot::kSVG, ps.TakePainting(), kTRUE);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Add special canvas objects like colors list at selected palette

void TWebCanvas::AddColorsPalette(TPadWebSnapshot &master)
{
   TObjArray *colors = (TObjArray *)gROOT->GetListOfColors();

   if (!colors)
      return;

   Int_t cnt = 0;
   for (Int_t n = 0; n <= colors->GetLast(); ++n)
      if (colors->At(n))
         cnt++;

   if (cnt <= 598)
      return; // normally there are 598 colors defined

   TArrayI pal = TColor::GetPalette();

   auto *listofcols = new TWebPainting;
   for (Int_t n = 0; n <= colors->GetLast(); ++n)
      if (colors->At(n))
         listofcols->AddColor(n, (TColor *)colors->At(n));

   // store palette in the buffer
   auto *tgt = listofcols->Reserve(pal.GetSize());
   for (Int_t i = 0; i < pal.GetSize(); i++)
      tgt[i] = pal[i];
   listofcols->FixSize();

   master.NewSpecials().SetSnapshot(TWebSnapshot::kColors, listofcols, kTRUE);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Create snapshot for pad and all primitives
/// Callback function is used to create JSON in the middle of data processing -
/// when all misc objects removed from canvas list of primitives or histogram list of functions
/// After that objects are moved back to their places

void TWebCanvas::CreatePadSnapshot(TPadWebSnapshot &paddata, TPad *pad, Long64_t version, PadPaintingReady_t resfunc)
{
   paddata.SetActive(pad == gPad);
   paddata.SetObjectIDAsPtr(pad);
   paddata.SetSnapshot(TWebSnapshot::kSubPad, pad); // add ref to the pad

   if (resfunc && (GetStyleDelivery() > (version > 0 ? 1 : 0)))
      paddata.NewPrimitive().SetSnapshot(TWebSnapshot::kStyle, gStyle);

   TList *primitives = pad->GetListOfPrimitives();

   if (primitives) fPrimitivesLists.Add(primitives); // add list of primitives

   TWebPS masterps;
   bool usemaster = primitives ? (primitives->GetSize() > fPrimitivesMerge) : false;

   TIter iter(primitives);
   TObject *obj = nullptr;
   TFrame *frame = nullptr;
   TPaveText *title = nullptr;
   bool need_frame = false;
   std::string need_title;

   while ((obj = iter()) != nullptr) {
      if (obj->InheritsFrom(TFrame::Class())) {
         frame = static_cast<TFrame *>(obj);
      } else if (obj->InheritsFrom(TH1::Class())) {
         need_frame = true;
         if (!obj->TestBit(TH1::kNoTitle) && (strlen(obj->GetTitle())>0)) need_title = obj->GetTitle();
      } else if (obj->InheritsFrom(TGraph::Class())) {
         need_frame = true;
         if (strlen(obj->GetTitle())>0) need_title = obj->GetTitle();
      } else if (obj->InheritsFrom(TPaveText::Class())) {
         if (strcmp(obj->GetName(),"title") == 0)
            title = static_cast<TPaveText *>(obj);
      }
   }

   if (need_frame && !frame && CanCreateObject("TFrame")) {
      frame = pad->GetFrame();
      primitives->AddFirst(frame);
   }

   if (!need_title.empty()) {
      if (title) {
         auto line0 = title->GetLine(0);
         if (line0 && !IsReadOnly()) line0->SetTitle(need_title.c_str());
      } else if (CanCreateObject("TPaveText")) {
         title = new TPaveText(0, 0, 0, 0, "blNDC");
         title->SetFillColor(gStyle->GetTitleFillColor());
         title->SetFillStyle(gStyle->GetTitleStyle());
         title->SetName("title");
         title->SetBorderSize(gStyle->GetTitleBorderSize());
         title->SetTextColor(gStyle->GetTitleTextColor());
         title->SetTextFont(gStyle->GetTitleFont(""));
         if (gStyle->GetTitleFont("") % 10 > 2)
            title->SetTextSize(gStyle->GetTitleFontSize());
         title->AddText(need_title.c_str());
         title->SetBit(kCanDelete);
         primitives->Add(title);
      }
   }

   auto flush_master = [&]() {
      if (!usemaster || masterps.IsEmptyPainting()) return;

      paddata.NewPrimitive(pad).SetSnapshot(TWebSnapshot::kSVG, masterps.TakePainting(), kTRUE);
      masterps.CreatePainting(); // create for next operations
   };

   iter.Reset();

   bool first_obj = true;

   while ((obj = iter()) != nullptr) {
      if (obj->InheritsFrom(TPad::Class())) {
         flush_master();
         CreatePadSnapshot(paddata.NewSubPad(), (TPad *)obj, version, nullptr);
      } else if (obj->InheritsFrom(TH1::Class())) {
         flush_master();

         TH1 *hist = (TH1 *)obj;
         TIter fiter(hist->GetListOfFunctions());
         TObject *fobj = nullptr;
         TPaveStats *stats = nullptr;
         TObject *palette = nullptr;

         while ((fobj = fiter()) != nullptr) {
           if (fobj->InheritsFrom(TPaveStats::Class()))
               stats = dynamic_cast<TPaveStats *> (fobj);
           else if (fobj->InheritsFrom("TPaletteAxis"))
              palette = fobj;
         }

         if (!stats && first_obj && CanCreateObject("TPaveStats")) {
            stats  = new TPaveStats(
                           gStyle->GetStatX() - gStyle->GetStatW(),
                           gStyle->GetStatY() - gStyle->GetStatH(),
                           gStyle->GetStatX(),
                           gStyle->GetStatY(), "brNDC");

             stats->SetParent(hist);
             stats->SetOptFit(gStyle->GetOptFit());
             stats->SetOptStat(gStyle->GetOptStat());
             stats->SetFillColor(gStyle->GetStatColor());
             stats->SetFillStyle(gStyle->GetStatStyle());
             stats->SetBorderSize(gStyle->GetStatBorderSize());
             stats->SetTextFont(gStyle->GetStatFont());
             if (gStyle->GetStatFont()%10 > 2)
                stats->SetTextSize(gStyle->GetStatFontSize());
             stats->SetFitFormat(gStyle->GetFitFormat());
             stats->SetStatFormat(gStyle->GetStatFormat());
             stats->SetName("stats");

             stats->SetTextColor(gStyle->GetStatTextColor());
             stats->SetTextAlign(12);
             stats->SetBit(kCanDelete);
             stats->SetBit(kMustCleanup);

             hist->GetListOfFunctions()->Add(stats);
         }

         TString hopt = iter.GetOption();

         if (!palette && CanCreateObject("TPaletteAxis") && (hist->GetDimension() > 1) &&
             (hopt.Index("colz", 0, TString::kIgnoreCase) != kNPOS)) {
            std::stringstream exec;
            exec << "new TPaletteAxis(0,0,0,0, (TH1*)" << std::hex << std::showbase << (size_t)hist << ");";
            palette = (TObject *)gROOT->ProcessLine(exec.str().c_str());
            if (palette)
               hist->GetListOfFunctions()->AddFirst(palette);
         }

         if (title && first_obj) hopt.Append(";;use_pad_title");
         if (stats) hopt.Append(";;use_pad_stats");
         if (palette) hopt.Append(";;use_pad_palette");

         paddata.NewPrimitive(obj, hopt.Data()).SetSnapshot(TWebSnapshot::kObject, obj);

         fiter.Reset();
         while ((fobj = fiter()) != nullptr)
            CreateObjectSnapshot(paddata, pad, fobj, fiter.GetOption());

         fPrimitivesLists.Add(hist->GetListOfFunctions());
         first_obj = false;
      } else if (obj->InheritsFrom(TGraph::Class())) {
         flush_master();

         TGraph *gr = (TGraph *)obj;

         TIter fiter(gr->GetListOfFunctions());
         TObject *fobj = nullptr;
         TPaveStats *stats = nullptr;

         while ((fobj = fiter()) != nullptr) {
           if (fobj->InheritsFrom(TPaveStats::Class()))
               stats = dynamic_cast<TPaveStats *> (fobj);
         }

         // ensure histogram exists on server to draw it properly on clients side
         if (!IsReadOnly() && first_obj)
            gr->GetHistogram();

         TString gropt = iter.GetOption();
         if (title && first_obj) gropt.Append(";;use_pad_title");
         if (stats) gropt.Append(";;use_pad_stats");

         paddata.NewPrimitive(obj, gropt.Data()).SetSnapshot(TWebSnapshot::kObject, obj);

         fiter.Reset();
         while ((fobj = fiter()) != nullptr)
            CreateObjectSnapshot(paddata, pad, fobj, fiter.GetOption());

         fPrimitivesLists.Add(gr->GetListOfFunctions());
         first_obj = false;
      } else if (IsJSSupportedClass(obj)) {
         flush_master();
         paddata.NewPrimitive(obj, iter.GetOption()).SetSnapshot(TWebSnapshot::kObject, obj);
      } else {
         CreateObjectSnapshot(paddata, pad, obj, iter.GetOption(), usemaster ? &masterps : nullptr);
      }
   }

   flush_master();

   bool provide_colors = GetPaletteDelivery() > 0;
   if (GetPaletteDelivery() == 1)
      provide_colors = !!resfunc && (version <= 0);
   else if (GetPaletteDelivery() == 2)
      provide_colors = !!resfunc;

   // add specials after painting is performed - new colors may be generated only during painting
   if (provide_colors)
      AddColorsPalette(paddata);

   if (!resfunc)
      return;

   // now move all primitives and functions into separate list to perform I/O

   TList save_lst;
   TIter diter(&fPrimitivesLists);
   TList *dlst = nullptr;
   while ((dlst = (TList *)diter()) != nullptr) {
      TIter fiter(dlst);
      while ((obj = fiter()) != nullptr)
         save_lst.Add(obj, fiter.GetOption());
      save_lst.Add(dlst); // add list itself to have marker
      dlst->Clear("nodelete");
   }

   // execute function to prevent storing of colors with custom TCanvas streamer
   // TODO: Olivier - we need to change logic here!
   TColor::DefinedColors();

   // invoke callback for master painting
   resfunc(&paddata);

   TIter siter(&save_lst);
   diter.Reset();
   while ((dlst = (TList *)diter()) != nullptr) {
      while ((obj = siter()) != nullptr) {
         if (obj == dlst)
            break;
         dlst->Add(obj, siter.GetOption());
      }
   }

   save_lst.Clear("nodelete");

   fPrimitivesLists.Clear("nodelete");
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Add message to send queue for specified connection
/// If connid == 0, message will be add to all connections
/// Return kFALSE if queue is full or connection is not exists

Bool_t TWebCanvas::AddToSendQueue(unsigned connid, const std::string &msg)
{
   Bool_t res = false;
   for (auto &conn : fWebConn) {
      if ((conn.fConnId == connid) || (connid == 0)) {
         conn.fSend.emplace(msg);
         res = kTRUE;
      }
   }
   return res;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/// Check if any data should be send to client
/// If connid != 0, only selected connection will be checked

void TWebCanvas::CheckDataToSend(unsigned connid)
{
   if (!Canvas())
      return;

   for (auto &conn : fWebConn) {
      if (connid && (conn.fConnId != connid))
         continue;

      // check if direct data sending is possible
      if (!fWindow->CanSend(conn.fConnId, true))
         continue;

      std::string buf;

      if ((conn.fSendVersion < fCanvVersion) && (conn.fSendVersion == conn.fDrawVersion)) {

         buf = "SNAP6:";

         TCanvasWebSnapshot holder(IsReadOnly(), fCanvVersion);

         // scripts send only when canvas drawn for the first time
         if (!conn.fSendVersion)
            holder.SetScripts(fCustomScripts);

         CreatePadSnapshot(holder, Canvas(), conn.fSendVersion, [&buf,this](TPadWebSnapshot *snap) {
            buf.append(TBufferJSON::ToJSON(snap, fJsonComp).Data());
         });

         conn.fSendVersion = fCanvVersion;

      } else if (!conn.fSend.empty()) {

         std::swap(buf, conn.fSend.front());
         conn.fSend.pop();

      }

      if (!buf.empty())
         fWindow->Send(conn.fConnId, buf);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Close web canvas - not implemented

void TWebCanvas::Close()
{
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Show canvas in specified place.
/// If parameter args not specified, default ROOT web display will be used

void TWebCanvas::ShowWebWindow(const ROOT::Experimental::RWebDisplayArgs &args)
{
   if (!fWindow) {
      fWindow = ROOT::Experimental::RWebWindow::Create();

      fWindow->SetConnLimit(0); // configure connections limit

      fWindow->SetDefaultPage("file:rootui5sys/canv/canvas6.html");

      fWindow->SetCallBacks(
         // connection
         [this](unsigned connid) {
            fWebConn.emplace_back(connid);
            CheckDataToSend(connid);
         },
         // data
         [this](unsigned connid, const std::string &arg) {
            ProcessData(connid, arg);
            CheckDataToSend(connid);
         },
         // disconnect
         [this](unsigned connid) {
            unsigned indx = 0;
            for (auto &c : fWebConn) {
               if (c.fConnId == connid) {
                  fWebConn.erase(fWebConn.begin() + indx);
                  break;
               }
               indx++;
            }
         });
   }

   auto w = Canvas()->GetWw(), h = Canvas()->GetWh();

   if ((w > 10) && (w < 50000) && (h > 10) && (h < 30000))
      fWindow->SetGeometry(w + 6, h + 22);

   fWindow->Show(args);
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Show canvas in browser window

void TWebCanvas::Show()
{
   ShowWebWindow();
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Function used to send command to browser to toggle menu, toolbar, editors, ...

void TWebCanvas::ShowCmd(const std::string &arg, Bool_t show)
{
   if (AddToSendQueue(0, "SHOW:"s + arg + (show ? ":1"s : ":0"s)))
      CheckDataToSend();
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Activate object in editor in web browser

void TWebCanvas::ActivateInEditor(TPad *pad, TObject *obj)
{
   if (!pad || !obj) return;

   UInt_t hash = TString::Hash(&obj, sizeof(obj));

   if (AddToSendQueue(0, "EDIT:"s + std::to_string(hash)))
      CheckDataToSend();
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if web canvas has graphical editor

Bool_t TWebCanvas::HasEditor() const
{
   return (fClientBits & TCanvas::kShowEditor) != 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if web canvas has menu bar

Bool_t TWebCanvas::HasMenuBar() const
{
   return (fClientBits & TCanvas::kMenuBar) != 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if web canvas has status bar

Bool_t TWebCanvas::HasStatusBar() const
{
   return (fClientBits & TCanvas::kShowEventStatus) != 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if tooltips are activated in web canvas

Bool_t TWebCanvas::HasToolTips() const
{
   return (fClientBits & TCanvas::kShowToolTips) != 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Assign clients bits

void TWebCanvas::AssignStatusBits(UInt_t bits)
{
   fClientBits = bits;
   Canvas()->SetBit(TCanvas::kShowEventStatus, bits & TCanvas::kShowEventStatus);
   Canvas()->SetBit(TCanvas::kShowEditor, bits & TCanvas::kShowEditor);
   Canvas()->SetBit(TCanvas::kShowToolTips, bits & TCanvas::kShowToolTips);
   Canvas()->SetBit(TCanvas::kMenuBar, bits & TCanvas::kMenuBar);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Decode all pad options, which includes ranges plus objects options

Bool_t TWebCanvas::DecodePadOptions(const std::string &msg)
{
   if (IsReadOnly() || msg.empty())
      return kFALSE;

   auto arr = TBufferJSON::FromJSON<std::vector<TWebPadOptions>>(msg);

   if (!arr)
      return kFALSE;

   for (unsigned n = 0; n < arr->size(); ++n) {
      auto &r = arr->at(n);
      TPad *pad = dynamic_cast<TPad *>(FindPrimitive(r.snapid));

      if (!pad)
         continue;

      if (pad == Canvas()) AssignStatusBits(r.bits);

      if (r.active && (pad != gPad)) gPad = pad;

      pad->SetTicks(r.tickx, r.ticky);
      pad->SetGrid(r.gridx, r.gridy);
      if (r.logx != pad->GetLogx())
         pad->SetLogx(r.logx);
      if (r.logy != pad->GetLogy())
         pad->SetLogy(r.logy);
      if (r.logz != pad->GetLogz())
         pad->SetLogz(r.logz);

      pad->SetLeftMargin(r.mleft);
      pad->SetRightMargin(r.mright);
      pad->SetTopMargin(r.mtop);
      pad->SetBottomMargin(r.mbottom);

      if (r.ranges) {

         Double_t ux1_, ux2_, uy1_, uy2_, px1_, px2_, py1_, py2_;

         pad->GetRange(px1_, py1_, px2_, py2_);
         pad->GetRangeAxis(ux1_, uy1_, ux2_, uy2_);

         bool same_range = (r.ux1 == ux1_) && (r.ux2 == ux2_) && (r.uy1 == uy1_) && (r.uy2 == uy2_) &&
                           (r.px1 == px1_) && (r.px2 == px2_) && (r.py1 == py1_) && (r.py2 == py2_);

         if (!same_range) {
            pad->RangeAxis(r.ux1, r.uy1, r.ux2, r.uy2);

            pad->Range(r.px1, r.py1, r.px2, r.py2);

            if (gDebug > 1)
               Info("DecodeAllRanges", "Change ranges for pad %s", pad->GetName());
         }
      }

      pad->SetPad(r.mleft, r.mbottom, 1-r.mright, 1-r.mtop);

      TH1 *hist = static_cast<TH1 *>(FindPrimitive("histogram", pad));

      if (hist) {
         double hmin = 0, hmax = 0;

         if (r.zx1 == r.zx2)
            hist->GetXaxis()->SetRange(0,0);
         else
            hist->GetXaxis()->SetRangeUser(r.zx1, r.zx2);

         if (hist->GetDimension() == 1) {
            hmin = r.zy1;
            hmax = r.zy2;
         } else if (r.zy1 == r.zy2) {
            hist->GetYaxis()->SetRange(0,0);
         } else {
            hist->GetYaxis()->SetRangeUser(r.zy1, r.zy2);
         }

         if (hist->GetDimension() == 2) {
            hmin = r.zz1;
            hmax = r.zz2;
         } else if (hist->GetDimension() == 3) {
            if (r.zz1 == r.zz2) {
               hist->GetZaxis()->SetRange(0,0);
            } else {
              hist->GetZaxis()->SetRangeUser(r.zz1, r.zz2);
            }
         }

         if (hmin == hmax) { hist->SetMinimum(); hist->SetMaximum(); }
                      else { hist->SetMinimum(hmin); hist->SetMaximum(hmax); }

      }

      for (auto &item : r.primitives)
         ProcessObjectOptions(item, pad);

      // without special objects no need for explicit update of the pad
      if (fHasSpecials)
         pad->Modified(kTRUE);
   }

   if (fUpdatedSignal) fUpdatedSignal(); // invoke signal

   return kTRUE;
}


//////////////////////////////////////////////////////////////////////////////////////////
/// Handle data from web browser
/// Returns kFALSE if message was not processed

Bool_t TWebCanvas::ProcessData(unsigned connid, const std::string &arg)
{
   if (arg.empty())
      return kTRUE;

   // try to identify connection for given WS request
   unsigned indx = 0;
   for (auto &c : fWebConn) {
      if (c.fConnId == connid) break;
      indx++;
   }
   if (indx >= fWebConn.size())
      return kTRUE;

   const char *cdata = arg.c_str();

   if (arg == "KEEPALIVE") {
      // do nothing

   } else if (arg == "QUIT") {

      // use window manager to correctly terminate http server
      fWindow->TerminateROOT();

   } else if (arg.compare(0, 7, "READY6:") == 0) {

      // this is reply on drawing of ROOT6 snapshot
      // it confirms when drawing of specific canvas version is completed

      cdata += 7;

      const char *separ = strchr(cdata, ':');
      if (!separ) {
         fWebConn[indx].fDrawVersion = std::stoll(cdata);
      } else {
         fWebConn[indx].fDrawVersion = std::stoll(std::string(cdata, separ - cdata));
         if ((indx == 0) && !IsReadOnly())
            DecodePadOptions(separ+1);
      }

   } else if (arg == "RELOAD") {

      // trigger reload of canvas data
      fWebConn[indx].fSendVersion = fWebConn[indx].fDrawVersion = 0;

   } else if (arg.compare(0, 5, "SAVE:") == 0) {

      // save image produced by the client side - like png or svg
      const char *img = cdata + 5;

      const char *separ = strchr(img, ':');
      if (separ) {
         TString filename(img, separ - img);
         img = separ + 1;

         std::ofstream ofs(filename.Data());

         if (filename.Index(".svg") != kNPOS) {
            // ofs << "<?xml version=\"1.0\" standalone=\"no\"?>";
            ofs << img;
         } else {
            TString binary = TBase64::Decode(img);
            ofs.write(binary.Data(), binary.Length());
         }
         ofs.close();

         Info("ProcessData", "File %s has been created", filename.Data());
      }

   } else if (arg.compare(0, 8, "PRODUCE:") == 0) {

      // create ROOT, PDF, ... files using native ROOT functionality
      Canvas()->Print(arg.c_str() + 8);

   } else if (arg.compare(0, 9, "OPTIONS6:") == 0) {

      if ((indx == 0) && !IsReadOnly())
         DecodePadOptions(arg.substr(9));

   } else if (arg.compare(0, 11, "STATUSBITS:") == 0) {

      if (indx == 0) {
         AssignStatusBits(std::stoul(arg.substr(11)));
         if (fUpdatedSignal) fUpdatedSignal(); // invoke signal
      }

   } else if (IsReadOnly()) {

      // all following messages are not allowed in readonly mode
      return kFALSE;

   } else if (arg.compare(0, 8, "GETMENU:") == 0) {

      TObject *obj = FindPrimitive(arg.substr(8));
      if (!obj)
         obj = Canvas();

      TWebMenuItems items(arg.c_str() + 8);
      items.PopulateObjectMenu(obj, obj->IsA());
      std::string buf = "MENU:";
      buf.append(TBufferJSON::ToJSON(&items, 103).Data());

      AddToSendQueue(connid, buf);

   } else if (arg.compare(0, 8, "PRIMIT6:") == 0) {

      if (IsFirstConn(connid) && !IsReadOnly()) { // only first connection can modify object

         auto opt = TBufferJSON::FromJSON<TWebObjectOptions>(arg.c_str() + 8);

         if (opt) {
            TPad *modpad = ProcessObjectOptions(*opt, nullptr);

            // indicate that pad was modified
            if (modpad)
               modpad->Modified();
         }
      }

   } else if (arg.compare(0, 11, "PADCLICKED:") == 0) {

      auto click = TBufferJSON::FromJSON<TWebPadClick>(arg.c_str() + 11);

      if (click && IsFirstConn(connid) && !IsReadOnly()) {

         TPad *pad = dynamic_cast<TPad *>(FindPrimitive(click->padid));
         if (pad && (pad != gPad)) {
            Info("ProcessData", "Activate pad %s", pad->GetName());
            gPad = pad;
            Canvas()->SetClickSelectedPad(pad);
            if (fActivePadChangedSignal)
               fActivePadChangedSignal(pad);
         }

         if (!click->objid.empty()) {
            TObject *selobj = FindPrimitive(click->objid);
            Canvas()->SetClickSelected(selobj);
            if (pad && selobj && fObjSelectSignal)
               fObjSelectSignal(pad, selobj);
         }

         if ((click->x >= 0) && (click->y >= 0)) {
            if (click->dbl && fPadDblClickedSignal)
               fPadDblClickedSignal(pad, click->x, click->y);
            else if (fPadClickedSignal)
               fPadClickedSignal(pad, click->x, click->y);
         }
      }

   } else if (arg.compare(0, 8, "OBJEXEC:") == 0) {

     auto buf = arg.substr(8);
     auto pos = buf.find(":");

     if ((pos > 0) && IsFirstConn(connid) && !IsReadOnly()) { // only first client can execute commands
        auto sid = buf.substr(0, pos);
        buf.erase(0, pos + 1);

        TObject *obj = FindPrimitive(sid);
        if (obj && !buf.empty()) {

           while (!buf.empty()) {
              std::string sub = buf;
              pos = buf.find(";;");
              if (pos == std::string::npos) {
                 sub = buf;
                 buf.clear();
              } else {
                 sub = buf.substr(0,pos);
                 buf = buf.substr(pos+2);
              }
              if (sub.empty()) continue;

              std::stringstream exec;
              exec << "((" << obj->ClassName() << " *) " << std::hex << std::showbase << (size_t)obj << ")->" << sub << ";";
              Info("ProcessData", "Obj %s Execute %s", obj->GetName(), exec.str().c_str());
              gROOT->ProcessLine(exec.str().c_str());
           }

           CheckPadModified(Canvas());
        }
     }

   } else if (arg.compare(0, 12, "EXECANDSEND:") == 0) {

      // execute method and send data, used by drawing projections

      std::string buf = arg.substr(12);
      std::string reply;
      TObject *obj = nullptr;

      auto pos = buf.find(":");

      if ((pos > 0) && IsFirstConn(connid) && !IsReadOnly()) {
         // only first client can execute commands
         reply = buf.substr(0, pos);
         buf.erase(0, pos + 1);
         pos = buf.find(":");
         if (pos > 0) {
            auto sid = buf.substr(0, pos);
            buf.erase(0, pos + 1);
            obj = FindPrimitive(sid);
         }
      }

      if (obj && !buf.empty() && !reply.empty()) {
         std::stringstream exec;
         exec << "((" << obj->ClassName() << " *) " << std::hex << std::showbase << (size_t)obj
              << ")->" << buf << ";";
         if (gDebug > 1)
            Info("ProcessData", "Obj %s Exec %s", obj->GetName(), exec.str().c_str());

         Long_t res = gROOT->ProcessLine(exec.str().c_str());
         TObject *resobj = (TObject *)res;
         if (resobj) {
            std::string send = reply;
            send.append(":");
            send.append(TBufferJSON::ToJSON(resobj, 23).Data());
            AddToSendQueue(connid, send);
            if (reply[0] == 'D')
               delete resobj; // delete object if first symbol in reply is D
         }
      }

   }  else {

      // unknown message, probably should be processed by other implementation
      return kFALSE;
   }

   return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns true if any pad in the canvas were modified
/// Reset modified flags, increment canvas version (if inc_version is true)

Bool_t TWebCanvas::CheckPadModified(TPad *pad, Bool_t inc_version)
{
   Bool_t modified = kFALSE;

   if (pad->IsModified()) {
      pad->Modified(kFALSE);
      modified = kTRUE;
   }

   TIter iter(pad->GetListOfPrimitives());
   TObject *obj = nullptr;
   while ((obj = iter()) != nullptr) {
      if (obj->InheritsFrom(TPad::Class()) && CheckPadModified(static_cast<TPad *>(obj), kFALSE))
         modified = kTRUE;
   }

   if (inc_version && modified)
      fCanvVersion++;

   return modified;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns window geometry including borders and menus

UInt_t TWebCanvas::GetWindowGeometry(Int_t &x, Int_t &y, UInt_t &w, UInt_t &h)
{
   x = 0;
   y = 0;
   w = Canvas()->GetWw() + 4;
   h = Canvas()->GetWh() + 28;
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// if canvas or any subpad was modified,
/// scan all primitives in the TCanvas and subpads and convert them into
/// the structure which will be delivered to JSROOT client

Bool_t TWebCanvas::PerformUpdate()
{
   CheckPadModified(Canvas());

   CheckDataToSend();

   // block in canvas update, can it be optional?
   WaitWhenCanvasPainted(fCanvVersion);

   return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Increment canvas version and force sending data to client - do not wit for reply

void TWebCanvas::ForceUpdate()
{
   fCanvVersion++;

   CheckDataToSend();
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Wait when specified version of canvas was painted and confirmed by browser

Bool_t TWebCanvas::WaitWhenCanvasPainted(Long64_t ver)
{
   // simple polling loop until specified version delivered to the clients

   long cnt = 0;

   if (gDebug > 2)
      Info("WaitWhenCanvasPainted", "version %ld", (long)ver);

   while (cnt++ < 1000) {

      if (!fWindow->HasConnection(0, false)) {
         if (gDebug > 2)
            Info("WaitWhenCanvasPainted", "no connections - abort");
         return kFALSE; // wait ~1 min if no new connection established
      }

      if ((fWebConn.size() > 0) && (fWebConn.front().fDrawVersion >= ver)) {
         if (gDebug > 2)
            Info("WaitWhenCanvasPainted", "ver %ld got painted", (long)ver);
         return kTRUE;
      }

      gSystem->ProcessEvents();

      gSystem->Sleep((cnt < 500) ? 1 : 100); // increase sleep interval when do very often
   }

   if (gDebug > 2)
      Info("WaitWhenCanvasPainted", "timeout");

   return kFALSE;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Create JSON painting output for given canvas
/// Produce JSON can be used for offline drawing with JSROOT

TString TWebCanvas::CreateCanvasJSON(TCanvas *c, Int_t json_compression)
{
   TString res;

   if (!c)
      return res;

   Bool_t isbatch = c->IsBatch();
   c->SetBatch(kTRUE);

   {
      auto imp = std::make_unique<TWebCanvas>(c, c->GetName(), 0, 0, 1000, 500);

      TCanvasWebSnapshot holder(true, 1); // always readonly

      imp->CreatePadSnapshot(holder, c, 0, [&res, json_compression](TPadWebSnapshot *snap) {
         res = TBufferJSON::ToJSON(snap, json_compression);
      });
   }

   c->SetBatch(isbatch);
   return res;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Create JSON painting output for given canvas and store into the file
/// See TBufferJSON::ExportToFile() method for more details

Int_t TWebCanvas::StoreCanvasJSON(TCanvas *c, const char *filename, const char *option)
{
   Int_t res{0};

   if (!c)
      return res;

   Bool_t isbatch = c->IsBatch();
   c->SetBatch(kTRUE);

   {
      auto imp = std::make_unique<TWebCanvas>(c, c->GetName(), 0, 0, 1000, 500);

      TCanvasWebSnapshot holder(true, 1); // always readonly

      imp->CreatePadSnapshot(holder, c, 0, [&res, filename, option](TPadWebSnapshot *snap) {
         res = TBufferJSON::ExportToFile(filename, snap, option);
      });
   }

   c->SetBatch(isbatch);
   return res;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Create image using batch (headless) capability of Chrome browser
/// Supported png, jpeg, svg, pdf formats

bool TWebCanvas::ProduceImage(TCanvas *c, const char *fileName, Int_t width, Int_t height)
{
   if (!c)
      return false;

   auto json = TWebCanvas::CreateCanvasJSON(c, TBufferJSON::kNoSpaces + TBufferJSON::kSameSuppression);
   if (!json.Length())
      return false;

   return ROOT::Experimental::RWebDisplayHandle::ProduceImage(fileName, json.Data(), width ? width : c->GetWw(), height ? height : c->GetWh());
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Process data for single primitive
/// Returns object pad if object was modified

TPad *TWebCanvas::ProcessObjectOptions(TWebObjectOptions &item, TPad *pad)
{
   TObjLink *lnk = nullptr;
   TPad *objpad = nullptr;
   TObject *obj = FindPrimitive(item.snapid, pad, &lnk, &objpad);

   if (item.fcust.compare("exec") == 0) {
      auto pos = item.opt.find("(");
      if (obj && (pos != std::string::npos) && obj->IsA()->GetMethodAllAny(item.opt.substr(0,pos).c_str())) {
         std::stringstream exec;
         exec << "((" << obj->ClassName() << " *) " << std::hex << std::showbase
                      << (size_t)obj << ")->" << item.opt << ";";
         Info("ProcessObjectOptions", "Obj %s Execute %s", obj->GetName(), exec.str().c_str());
         gROOT->ProcessLine(exec.str().c_str());
      } else {
         Error("ProcessObjectOptions", "Fail to execute %s for object %p %s", item.opt.c_str(), obj, obj ? obj->ClassName() : "---");
         objpad = nullptr;
      }
      return objpad;
   }

   bool modified = false;

   if (obj && lnk) {
      auto pos = item.opt.find(";;use_"); // special coding of extra options
      if (pos != std::string::npos) item.opt.resize(pos);

      if (gDebug > 1)
         Info("ProcessObjectOptions", "Set draw option %s for object %s %s", item.opt.c_str(),
               obj->ClassName(), obj->GetName());

      lnk->SetOption(item.opt.c_str());

      modified = true;
   }

   if (item.fcust.compare("frame") == 0) {
      if (obj && obj->InheritsFrom(TFrame::Class())) {
         TFrame *frame = static_cast<TFrame *>(obj);
         if (item.fopt.size() >= 4) {
            frame->SetX1(item.fopt[0]);
            frame->SetY1(item.fopt[1]);
            frame->SetX2(item.fopt[2]);
            frame->SetY2(item.fopt[3]);
            modified = true;
         }
      }
   } else if (item.fcust.compare("pave") == 0) {
      if (obj && obj->InheritsFrom(TPave::Class())) {
         TPave *pave = static_cast<TPave *>(obj);
         if ((item.fopt.size() >= 4) && objpad) {
            auto save = gPad;
            gPad = objpad;

            // first time need to overcome init problem
            pave->ConvertNDCtoPad();

            pave->SetX1NDC(item.fopt[0]);
            pave->SetY1NDC(item.fopt[1]);
            pave->SetX2NDC(item.fopt[2]);
            pave->SetY2NDC(item.fopt[3]);
            modified = true;

            pave->ConvertNDCtoPad();
            gPad = save;
         }
      }
   }

   return modified ? objpad : nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Search of object with given id in list of primitives
/// One could specify pad where search could be start
/// Also if object is in list of primitives, one could ask for entry link for such object,
/// This can allow to change draw option

TObject *TWebCanvas::FindPrimitive(const std::string &sid, TPad *pad, TObjLink **padlnk, TPad **objpad)
{

   if (!pad)
      pad = Canvas();

   std::string kind;
   auto separ = sid.find("#");
   long unsigned id = 0;
   bool search_hist = false;

   if (sid == "histogram") {
      search_hist = true;
   } else if (separ == std::string::npos) {
      id = std::stoul(sid);
   } else {
      kind = sid.substr(separ + 1);
      id = std::stoul(sid.substr(0, separ));
   }

   if (!search_hist && TString::Hash(&pad, sizeof(pad)) == id)
      return pad;

   TObjLink *lnk = pad->GetListOfPrimitives()->FirstLink();
   while (lnk) {
      TObject *obj = lnk->GetObject();
      if (!obj) {
         lnk = lnk->Next();
         continue;
      }
      TH1 *h1 = obj->InheritsFrom(TH1::Class()) ? static_cast<TH1 *>(obj) : nullptr;
      TGraph *gr = obj->InheritsFrom(TGraph::Class()) ? static_cast<TGraph *>(obj) : nullptr;

      if (search_hist) {
         if (h1) return h1;
         if (gr) {
            auto offset = TGraph::Class()->GetDataMemberOffset("fHistogram");
            if (offset > 0) {
               return *((TH1 **)((char*) gr + offset));
            } else {
               printf("ERROR: Cannot access fHistogram data member in TGraph\n");
               return nullptr;
            }
         }
         lnk = lnk->Next();
         continue;
      }

      if (TString::Hash(&obj, sizeof(obj)) == id) {
         if (objpad)
            *objpad = pad;

         if (gr && (kind.find("hist")==0)) {
            // access to graph histogram
            obj = h1 = gr->GetHistogram();
            kind.erase(0,4);
            if (!kind.empty() && (kind[0]=='#')) kind.erase(0,1);
            padlnk = nullptr;
         }

         if (h1 && (kind == "x"))
            return h1->GetXaxis();
         if (h1 && (kind == "y"))
            return h1->GetYaxis();
         if (h1 && (kind == "z"))
            return h1->GetZaxis();

         if (!kind.empty() && (kind.compare(0,7,"member_") == 0)) {
            auto member = kind.substr(7);
            auto offset = obj->IsA()->GetDataMemberOffset(member.c_str());
            if (offset > 0) {
               TObject **mobj = (TObject **)((char*) obj + offset);
               return *mobj;
            }
            return nullptr;
         }

         if (padlnk)
            *padlnk = lnk;
         return obj;
      }
      if (h1 || gr) {
         TIter fiter(h1 ? h1->GetListOfFunctions() : gr->GetListOfFunctions());
         TObject *fobj = nullptr;
         while ((fobj = fiter()) != nullptr)
            if (TString::Hash(&fobj, sizeof(fobj)) == id) {
               if (objpad)
                  *objpad = pad;
               return fobj;
            }
      } else if (obj->InheritsFrom(TPad::Class())) {
         obj = FindPrimitive(sid, (TPad *)obj, padlnk, objpad);
         if (objpad && !*objpad)
            *objpad = pad;
         if (obj)
            return obj;
      }
      lnk = lnk->Next();
   }

   return nullptr;
}


