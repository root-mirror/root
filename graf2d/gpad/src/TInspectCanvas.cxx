// @(#)root/gpad:$Id$
// Author: Rene Brun   08/01/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TDatime.h"
#include "TGuiFactory.h"
#include "TInspectCanvas.h"
#include "TButton.h"
#include "TClass.h"
#include "TLine.h"
#include "TLink.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TRealData.h"
#include "strlcpy.h"
#include "snprintf.h"

ClassImp(TInspectCanvas);


/** \class TInspectorObject
\ingroup gpad

This class is designed to wrap a Foreign object in order to
inject it into the Browse sub-system.
*/

class TInspectorObject : public TObject
{
public:

   TInspectorObject(void *obj, TClass *cl) : fObj(obj),fClass(cl) {};
   ~TInspectorObject(){;}

   void   *GetObject() const { return fObj; };
   void    Inspect() const {
      gGuiFactory->CreateInspectorImp(this, 400, 200);
   };
   TClass *IsA() const { return fClass; }

private:
   void     *fObj;   //! pointer to the foreign object
   TClass   *fClass; //! pointer to class of the foreign object

};


/** \class TInspectCanvas
\ingroup gpad

A TInspectCanvas is a canvas specialized to inspect Root objects.
*/

////////////////////////////////////////////////////////////////////////////////
/// InspectCanvas default constructor.

TInspectCanvas::TInspectCanvas() : TCanvas()
{
   fBackward   = 0;
   fForward    = 0;
   fCurObject  = 0;
   fObjects    = 0;
   fLogx       = kFALSE;
   fLogy       = kFALSE;
   SetFillColor(0);
}

////////////////////////////////////////////////////////////////////////////////
/// InspectCanvas constructor.

TInspectCanvas::TInspectCanvas(UInt_t ww, UInt_t wh)
            : TCanvas("inspect","ROOT Object Inspector",ww,wh)
{
   fBackward   = 0;
   fForward    = 0;
   fCurObject  = 0;
   fObjects    = new TList;
   fLogx       = kFALSE;
   fLogy       = kFALSE;
   SetFillColor(0);
}

////////////////////////////////////////////////////////////////////////////////
/// InspectCanvas default destructor.

TInspectCanvas::~TInspectCanvas()
{
   if (fObjects) {
      fObjects->Clear("nodelete");
      delete fObjects;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Dump contents of obj in a graphics canvas.
/// Same action as TObject::Dump but in a graphical form.
/// In addition pointers to other objects can be followed.
///
/// The following picture is the Inspect of a histogram object:
/// \image html gpad_inspect.png

void TInspectCanvas::InspectObject(TObject *obj)
{
   Int_t cdate = 0;
   Int_t ctime = 0;
   UInt_t *cdatime = 0;
   Bool_t isdate = kFALSE;
   Bool_t isbits = kFALSE;
   const Int_t kname  = 1;
   const Int_t kvalue = 25;
   const Int_t ktitle = 37;
   const Int_t kline  = 1024;
   char line[kline];
   char *pname;

   TClass *cl = obj->IsA();
   if (cl == 0) return;
   TInspectorObject *proxy=0;
   if (!cl->IsTObject()) {
      // This is possible only if obj is actually a TInspectorObject
      // wrapping a non-TObject.
      proxy = (TInspectorObject*)obj;
      obj = (TObject*)proxy->GetObject();
   }

   if (!cl->GetListOfRealData()) cl->BuildRealData(obj);

   // Count number of data members in order to resize the canvas
   TRealData *rd;
   TIter      next(cl->GetListOfRealData());
   Int_t nreal = cl->GetListOfRealData()->GetSize();
   if (nreal == 0) return;

   Int_t nrows = 33;
   if (nreal+7 > nrows) nrows = nreal+7;
   Int_t nh = nrows*15;
   Int_t nw = 700;
   TVirtualPad *canvas = GetVirtCanvas();
   if (canvas) {
      canvas->Clear();                // remove primitives from canvas
      canvas->SetCanvasSize(nw, nh);  // set new size of drawing area
      canvas->Range(0,-3,20,nreal+4);
   }

   Float_t xvalue = 5;
   Float_t xtitle = 8;
   Float_t dy     = 1;
   Float_t ytext  = Float_t(nreal) - 1.5;
   Float_t tsize  = 0.99/ytext;
   if (tsize < 0.02) tsize = 0.02;
   if (tsize > 0.03) tsize = 0.03;

   // Create text objects
   TText tname, tvalue, ttitle;
   TText *tval;
   tname.SetTextFont(61);
   tname.SetTextAngle(0);
   tname.SetTextAlign(12);
   tname.SetTextColor(1);
   tname.SetTextSize(tsize);
   tvalue.SetTextFont(61);
   tvalue.SetTextAngle(0);
   tvalue.SetTextAlign(12);
   tvalue.SetTextColor(1);
   tvalue.SetTextSize(tsize);
   ttitle.SetTextFont(62);
   ttitle.SetTextAngle(0);
   ttitle.SetTextAlign(12);
   ttitle.SetTextColor(1);
   ttitle.SetTextSize(tsize);

   Float_t x1 = 0.2;
   Float_t x2 = 19.8;
   Float_t y1 = -0.5;
   Float_t y2 = Float_t(nreal) - 0.5;
   Float_t y3 = y2 + 1;
   Float_t y4 = y3 + 1.5;
   Float_t db = 25./GetWh();
   Float_t btop = 0.999;

   // Draw buttons
   fBackward = new TButton("backward","TInspectCanvas::GoBackward();",.01,btop-db,.15,btop);
   fBackward->Draw();
   fBackward->SetToolTipText("Inspect previous object");
   fForward  = new TButton("forward", "TInspectCanvas::GoForward();", .21,btop-db,.35,btop);
   fForward->Draw();
   fForward->SetToolTipText("Inspect next object");

   // Draw surrounding box and title areas
   TLine frame;
   frame.SetLineColor(1);
   frame.SetLineStyle(1);
   frame.SetLineWidth(1);
   frame.DrawLine(x1, y1, x2, y1);
   frame.DrawLine(x2, y1, x2, y4);
   frame.DrawLine(x2, y4, x1, y4);
   frame.DrawLine(x1, y4, x1, y1);
   frame.DrawLine(x1, y2, x2, y2);
   frame.DrawLine(x1, y3, x2, y3);
   frame.DrawLine(xvalue, y1, xvalue, y3);
   frame.DrawLine(xtitle, y1, xtitle, y3);
   ttitle.SetTextSize(0.8*tsize);
   ttitle.SetTextAlign(21);
   ttitle.DrawText(0.5*(x1+xvalue), y2+0.1, "Member Name");
   ttitle.DrawText(0.5*(xvalue+xtitle), y2+0.1, "Value");
   ttitle.DrawText(0.5*(xtitle+x2), y2+0.1, "Title");
   ttitle.SetTextSize(1.2*tsize);
   ttitle.SetTextColor(2);
   ttitle.SetTextAlign(11);
   ttitle.DrawText(x1+0.2, y3+0.1, cl->GetName());
   if (proxy==0) {
      ttitle.SetTextColor(4);
      strlcpy(line,obj->GetName(),kline);
      ttitle.DrawText(xvalue+0.2, y3+0.1, line);
      ttitle.SetTextColor(6);
      ttitle.DrawText(xtitle+2, y3+0.1, obj->GetTitle());
   } else {
      ttitle.SetTextColor(4);
      snprintf(line,1023,"%s:%d","Foreign object",0);
      ttitle.DrawText(xvalue+0.2, y3+0.1, line);
      ttitle.SetTextColor(6);
      ttitle.DrawText(xtitle+2, y3+0.1, "no title given");
   }
   ttitle.SetTextSize(tsize);
   ttitle.SetTextColor(1);
   ttitle.SetTextFont(11);
   ttitle.SetTextAlign(12);

   //---Now loop on data members-----------------------
   // We make 3 passes. Faster than one single pass because changing
   // font parameters is time consuming
   for (Int_t pass = 0; pass < 3; pass++) {
      ytext  = y2 - 0.5;
      next.Reset();
      while ((rd = (TRealData*) next())) {
         TDataMember *member = rd->GetDataMember();
         if (!member) continue;
         TDataType *membertype = member->GetDataType();
         isdate = kFALSE;
         if (strcmp(member->GetName(),"fDatime") == 0 && membertype && membertype->GetType() == kUInt_t) {
            isdate = kTRUE;
         }
         isbits = kFALSE;
         if (strcmp(member->GetName(),"fBits") == 0 && membertype && membertype->GetType() == kUInt_t) {
            isbits = kTRUE;
         }

         // Encode data member name
         pname = &line[kname];
         for (Int_t i=0;i<kline;i++) line[i] = ' ';
         line[kline-1] = 0;
         strlcpy(pname,rd->GetName(),kline-kname);
         if (strstr(member->GetFullTypeName(),"**")) strlcat(pname,"**",kline-kname);

         // Encode data value or pointer value
         tval = &tvalue;
         Int_t offset = rd->GetThisOffset();
         char *pointer = (char*)obj + offset;
         char **ppointer = (char**)(pointer);
         TLink *tlink = 0;

         TClass *clm=0;
         if (!membertype) {
            clm = member->GetClass();
         }

         if (member->IsaPointer()) {
            char **p3pointer = (char**)(*ppointer);
            if (clm && !clm->IsStartingWithTObject() ) {
               //NOTE: memory leak!
               p3pointer = (char**)new TInspectorObject(p3pointer,clm);
            }

            if (!p3pointer) {
               snprintf(&line[kvalue],kline-kvalue,"->0");
            } else if (!member->IsBasic()) {
               if (pass == 1) {
                  tlink = new TLink(xvalue+0.1, ytext, p3pointer);
               }
            } else if (membertype) {
               if (!strcmp(membertype->GetTypeName(), "char"))
                  strlcpy(&line[kvalue], *ppointer,kline-kvalue);
               else
                  strlcpy(&line[kvalue], membertype->AsString(p3pointer),kline-kvalue);
            } else if (!strcmp(member->GetFullTypeName(), "char*") ||
                     !strcmp(member->GetFullTypeName(), "const char*")) {
               strlcpy(&line[kvalue], *ppointer,kline-kvalue);
            } else {
               if (pass == 1) tlink = new TLink(xvalue+0.1, ytext, p3pointer);
            }
         } else if (membertype)
            if (isdate) {
               cdatime = (UInt_t*)pointer;
               TDatime::GetDateTime(cdatime[0],cdate,ctime);
               snprintf(&line[kvalue],kline-kvalue,"%d/%d",cdate,ctime);
            } else if (isbits) {
               snprintf(&line[kvalue],kline-kvalue,"0x%08x", *(UInt_t*)pointer);
            } else {
               strlcpy(&line[kvalue], membertype->AsString(pointer),kline-kvalue);
            }
         else
            snprintf(&line[kvalue],kline-kvalue,"->%lx ", (Long_t)pointer);

         // Encode data member title
         Int_t ltit = 0;
         if (isdate == kFALSE && strcmp(member->GetFullTypeName(), "char*") &&
             strcmp(member->GetFullTypeName(), "const char*")) {
            Int_t lentit = strlen(member->GetTitle());
            if (lentit >= kline-ktitle) lentit = kline-ktitle-1;
            strlcpy(&line[ktitle],member->GetTitle(),kline-ktitle);
            line[ktitle+lentit] = 0;
            ltit = ktitle;
         }

         // Ready to draw the name, value and title columns
         if (pass == 0)tname.DrawText( x1+0.1,  ytext, &line[kname]);
         if (pass == 1) {
            if (tlink) {
               tlink->SetTextFont(61);
               tlink->SetTextAngle(0);
               tlink->SetTextAlign(12);
               tlink->SetTextColor(2);
               tlink->SetTextSize(tsize);
               tlink->SetBit(kCanDelete);
               tlink->Draw();
               if (strstr(member->GetFullTypeName(),"**")) tlink->SetBit(TLink::kIsStarStar);
               tlink->SetName(member->GetTypeName());
            } else {
               tval->DrawText(xvalue+0.1, ytext, &line[kvalue]);
            }
         }
         if (pass == 2 && ltit) ttitle.DrawText(xtitle+0.3, ytext, &line[ltit]);
         ytext -= dy;
      }
   }
   Update();
   fCurObject = obj;
}

////////////////////////////////////////////////////////////////////////////////
/// static function , inspect previous object

void TInspectCanvas::GoBackward()
{
   TInspectCanvas *inspect = (TInspectCanvas*)(gROOT->GetListOfCanvases())->FindObject("inspect");
   if (!inspect) return;
   TObject *cur = inspect->GetCurObject();
   TObject *obj = inspect->GetObjects()->Before(cur);
   if (obj)       inspect->InspectObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// static function , inspect next object

void TInspectCanvas::GoForward()
{
   TInspectCanvas *inspect = (TInspectCanvas*)(gROOT->GetListOfCanvases())->FindObject("inspect");
   if (!inspect) return;
   TObject *cur = inspect->GetCurObject();
   TObject *obj = inspect->GetObjects()->After(cur);
   if (obj)       inspect->InspectObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// static function , interface to InspectObject.
/// Create the InspectCanvas if it does not exist yet.

void TInspectCanvas::Inspector(TObject *obj)
{
   TVirtualPad *padsav = gPad;
   TInspectCanvas *inspect = (TInspectCanvas*)(gROOT->GetListOfCanvases())->FindObject("inspect");
   if (!inspect) inspect = new TInspectCanvas(700,600);
   else          inspect->cd();

   inspect->InspectObject(obj);
   inspect->GetObjects()->Add(obj);
   //obj->SetBit(kMustCleanup);

   if (padsav) padsav->cd();
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from the list of objects.

void TInspectCanvas::RecursiveRemove(TObject *obj)
{
   fObjects->Remove(obj);
   TPad::RecursiveRemove(obj);
}
