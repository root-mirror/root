// @(#)root/gl:$Id$
// Author:  Matevz and Alja Tadel  20/02/2009

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLAnnotation.h"

#include "TGLIncludes.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGLUtil.h"
#include "TGLCamera.h"
#include "TGLRnrCtx.h"
#include "TGLSelectRecord.h"
#include "TGLViewerBase.h"
#include "TObjString.h"
#include "TGTextEdit.h"
#include "TGButton.h"
#include "TGLViewer.h"

#include "TMath.h"

#include <KeySymbols.h>

/** \class TGLAnnotation
\ingroup opengl
GL-overlay annotation.
*/

ClassImp(TGLAnnotation);

Color_t  TGLAnnotation::fgBackColor = kAzure + 10;
Color_t  TGLAnnotation::fgTextColor = kOrange;

////////////////////////////////////////////////////////////////////////////////

TGLAnnotation::TGLAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy) :
   TGLOverlayElement(TGLOverlayElement::kAnnotation),

   fPosX(posx), fPosY(posy),
   fMouseX(0),  fMouseY(0),
   fDrag(kNone),
   fDrawW(0), fDrawH(0), fTextSizeDrag(0),
   fActive(kFALSE),
   fMainFrame(0), fTextEdit(0),

   fParent(0),

   fText(text),
   fTextSize(0.03),
   fTextAlign(TGLFont::kLeft),
   fBackColor(fgBackColor),
   fTextColor(fgTextColor),
   fTransparency(100),
   fDrawRefLine(kFALSE),
   fUseColorSet(kTRUE),
   fAllowClose(kTRUE)
{
   // Constructor.
   // Create annotation as plain text

   parent->AddOverlayElement(this);
   fParent = (TGLViewer*)parent;
}

////////////////////////////////////////////////////////////////////////////////

TGLAnnotation::TGLAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy, TGLVector3 ref) :
   TGLOverlayElement(TGLOverlayElement::kAnnotation),
   fPosX(posx), fPosY(posy),
   fMouseX(0),  fMouseY(0),
   fDrag(kNone),
   fDrawW(0), fDrawH(0), fTextSizeDrag(0),
   fActive(kFALSE),
   fMainFrame(0), fTextEdit(0),

   fParent(0),

   fText(text),
   fTextSize(0.03),
   fTextAlign(TGLFont::kLeft),
   fBackColor(fgBackColor),
   fTextColor(fgTextColor),
   fTransparency(40),
   fDrawRefLine(kTRUE),
   fUseColorSet(kTRUE),
   fAllowClose(kTRUE)
{
   // Constructor.
   // Create annotation by picking an object.

   fPointer = ref;
   parent->AddOverlayElement(this);
   fParent = (TGLViewer*)parent;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLAnnotation::~TGLAnnotation()
{
   fParent->RemoveOverlayElement(this);
   delete fMainFrame;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle overlay event.
/// Return TRUE if event was handled.

Bool_t TGLAnnotation::Handle(TGLRnrCtx&          rnrCtx,
                             TGLOvlSelectRecord& selRec,
                             Event_t*            event)
{
   if (selRec.GetN() < 2) return kFALSE;
   Int_t recID = selRec.GetItem(1);
   switch (event->fType)
   {
      case kButtonPress:
      {
         fMouseX = event->fX;
         fMouseY = event->fY;
         fDrag = (recID == kResizeID) ? kResize : kMove;
         fTextSizeDrag = fTextSize;
         return kTRUE;
      }
      case kButtonRelease:
      {
         fDrag = kNone;
         if (recID == kDeleteID)
         {
            TGLViewer *v = fParent;
            delete this;
            v->RequestDraw(rnrCtx.ViewerLOD());
         }
         else if (recID == kEditID)
         {
            MakeEditor();
         }
         return kTRUE;
      }
      case kMotionNotify:
      {
         const TGLRect& vp = rnrCtx.RefCamera().RefViewport();
         if (vp.Width() == 0 || vp.Height() == 0) return kFALSE;

         if (fDrag == kMove)
         {
            fPosX += (Float_t)(event->fX - fMouseX) / vp.Width();
            fPosY -= (Float_t)(event->fY - fMouseY) / vp.Height();
            fMouseX = event->fX;
            fMouseY = event->fY;
            // Make sure we don't go offscreen (use fDraw variables set in draw)
            if (fPosX < 0)
               fPosX = 0;
            else if (fPosX + fDrawW > 1.0f)
               fPosX = 1.0f - fDrawW;
            if (fPosY < fDrawH)
               fPosY = fDrawH;
            else if (fPosY > 1.0f)
               fPosY = 1.0f;
         }
         else if (fDrag == kResize)
         {
            using namespace TMath;
            Float_t oovpw = 1.0f / vp.Width(), oovph = 1.0f / vp.Height();

            Float_t xw = oovpw * Min(Max(0, event->fX), vp.Width());
            Float_t yw = oovph * Min(Max(0, vp.Height() - event->fY), vp.Height());

            Float_t rx = Max((xw - fPosX) / (oovpw * fMouseX - fPosX), 0.0f);
            Float_t ry = Max((yw - fPosY) / (oovph*(vp.Height() - fMouseY) - fPosY), 0.0f);

            fTextSize  = Max(fTextSizeDrag * Min(rx, ry), 0.01f);
         }
         return kTRUE;
      }
      default:
      {
         return kFALSE;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse has entered overlay area.

Bool_t TGLAnnotation::MouseEnter(TGLOvlSelectRecord& /*rec*/)
{
   fActive = kTRUE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse has left overlay area.

void TGLAnnotation::MouseLeave()
{
   fActive = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Render the annotation.

void TGLAnnotation::Render(TGLRnrCtx& rnrCtx)
{
   const TGLRect& vp = rnrCtx.RefCamera().RefViewport();
   if (vp.Width() == 0 && vp.Height() == 0)
      return;

   Float_t old_depth_range[2];
   glGetFloatv(GL_DEPTH_RANGE, old_depth_range);
   glDepthRange(0, 0.001);


   glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
   TGLCapabilitySwitch lights_off(GL_LIGHTING, kFALSE);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glDisable(GL_CULL_FACE);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   // prepare colors
   Color_t bgCol = fBackColor;
   Color_t fgCol = fTextColor;

   if (fUseColorSet)
   {
      fgCol = rnrCtx.ColorSet().Markup().GetColorIndex();

      TColor* c1 = gROOT->GetColor(rnrCtx.ColorSet().Markup().GetColorIndex());
      TColor* c2 = gROOT->GetColor(rnrCtx.ColorSet().Background().GetColorIndex());

      if (c1 && c2) {
         Float_t f1 = 0.5, f2 = 0.5;
         bgCol = TColor::GetColor(c1->GetRed()  *f1  + c2->GetRed()  *f2,
                                  c1->GetGreen()*f1  + c2->GetGreen()*f2,
                                  c1->GetBlue() *f1  + c2->GetBlue() *f2);
      }
   }

   // reset matrix
   rnrCtx.ProjectionMatrixPushIdentity();

   glPushMatrix();
   // set ortho camera to [0,1] [0.1]
   glLoadIdentity();
   glTranslatef(-1.0f, -1.0f, 0.0f);
   glScalef(2.0f, 2.0f, 1.0f);

   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(0.1f, 1.0f);

   glPushMatrix();

   TGLUtil::LineWidth(1.0f);

   // move to pos
   glTranslatef(fPosX, fPosY, 0.0f);

   TObjArray  *lines = fText.Tokenize("\n");
   TIter       line_iter(lines);
   TObjString *osl;

   Float_t widthTxt, heightTxt, sx, sy, descent, line_height;
   {
      // get unscaled text size
      Int_t fs = TGLFontManager::GetFontSize(TMath::Nint(vp.Height()*fTextSize), 12, 64);
      rnrCtx.RegisterFontNoScale(fs, "arial", TGLFont::kTexture, fFont);
      descent     = fFont.GetDescent();
      line_height = fFont.GetLineHeight();

      Float_t llx, lly, llz, urx, ury, urz;
      widthTxt = heightTxt = 0;
      while ((osl = (TObjString*) line_iter()) != 0)
      {
         fFont.BBox(osl->GetString().Data(), llx, lly, llz, urx, ury, urz);
         widthTxt   = TMath::Max(widthTxt, urx);
         heightTxt += line_height;
      }
      widthTxt  += 2.0f * descent;
      heightTxt += 2.0f * descent;

      // keep proportions
      sy = fTextSize / (line_height + descent);
      sx = sy / vp.Aspect();
      fDrawW = sx*widthTxt;
      fDrawH = sy*heightTxt;
   }
   glScalef(sx, sy, 1.0f);

   glPushName(kMoveID);

   Float_t x1, x2, y1, y2;
   Float_t z3 =  0.0f;  // main background
   Float_t z2 = -0.01f; // outlines and text
   Float_t z1 = -0.02f; // button on top of text
   Float_t z0 = -0.03f; // button on top of text

   // main background
   glLoadName(kMoveID);
   x1 =  0.0f;
   x2 =  widthTxt;
   y1 = -heightTxt;
   y2 =  0.0f;
   TGLUtil::ColorTransparency(bgCol, fTransparency);
   glBegin(GL_QUADS);
   glVertex3f(x1, y1, z3);
   glVertex3f(x2, y1, z3);
   glVertex3f(x2, y2, z3);
   glVertex3f(x1, y2, z3);
   glEnd();
   // main polygon outline
   TGLUtil::ColorTransparency(fgCol, GetLineTransparency());
   glBegin(GL_LINE_LOOP);
   glVertex3f(x1, y1, z2);
   glVertex3f(x2, y1, z2);
   glVertex3f(x2, y2, z2);
   glVertex3f(x1, y2, z2);
   glEnd();

   // annotation text
   TGLUtil::Color(fgCol);
   fFont.PreRender();
   glPushMatrix();
   Float_t tx = 0;
   line_iter.Reset();
   while ((osl = (TObjString*) line_iter()) != 0)
   {
      if (fTextAlign == TGLFont::kLeft) {
         tx = 0;
      }
      else if  (fTextAlign == TGLFont::kCenterH) {
         tx = 0.5f * widthTxt - descent ;
      }
      else {
         tx = widthTxt - 2.0f * descent;
      }
      glTranslatef(0.0f, -line_height, 0.0f);
      fFont.Render(osl->GetString(), tx+descent, 0, z2, fTextAlign, TGLFont::kTop) ;
   }
   glPopMatrix();
   fFont.PostRender();

   delete lines;

   // buttons
   if (fActive)
   {
      Float_t bbox[6];
      fFont.PreRender();
      fFont.BBox("X", bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
      glLoadName(kEditID);
      fFont.Render("E", descent, descent, z2, fTextAlign, TGLFont::kTop);
      x2 = bbox[3] + 2.0f * descent;
      if (fAllowClose)
      {
         glLoadName(kDeleteID);
         fFont.Render("X", x2 + descent, descent, z2, fTextAlign, TGLFont::kTop);
      }
      fFont.PostRender();

      x1 = 0.0f;
      y1 = 0.0f;
      y2 = line_height + descent;
      {
         // edit button
         glLoadName(kEditID);
         // polygon
         TGLUtil::ColorTransparency(bgCol, fTransparency);
         glBegin(GL_QUADS);
         glVertex3f(x1, y1, z3);
         glVertex3f(x2, y1, z3);
         glVertex3f(x2, y2, z3);
         glVertex3f(x1, y2, z3);
         glEnd();
         //  outline
         TGLUtil::ColorTransparency(fgCol, GetLineTransparency());
         glBegin(GL_LINE_LOOP);
         glVertex3f(x1, y1, z0);
         glVertex3f(x2, y1, z0);
         glVertex3f(x2, y2, z0);
         glVertex3f(x1, y2, z0);
         glEnd();
      }
      x1 += x2;
      x2 += x2;
      if (fAllowClose)
      {
         // close button
         glLoadName(kDeleteID);
         // polygon
         TGLUtil::ColorTransparency(bgCol, fTransparency);
         glBegin(GL_QUADS);
         glVertex3f(x1, y1, z3);
         glVertex3f(x2, y1, z3);
         glVertex3f(x2, y2, z3);
         glVertex3f(x1, y2, z3);
         glEnd();
         //  outline
         TGLUtil::ColorTransparency(fgCol, GetLineTransparency());
         glBegin(GL_LINE_LOOP);
         glVertex3f(x1, y1, z0);
         glVertex3f(x2, y1, z0);
         glVertex3f(x2, y2, z0);
         glVertex3f(x1, y2, z0);
         glEnd();
      }
      {
         // resize button
         glLoadName(kResizeID);
         // polygon
         x1 =  widthTxt - line_height;
         x2 =  widthTxt;
         y1 = -heightTxt;
         y2 = -heightTxt + line_height;
         TGLUtil::ColorTransparency(bgCol, fTransparency);
         glBegin(GL_QUADS);
         glVertex3f(x1, y1, z1);
         glVertex3f(x2, y1, z1);
         glVertex3f(x2, y2, z1);
         glVertex3f(x1, y2, z1);
         glEnd();
         // draw resize corner lines
         TGLUtil::ColorTransparency(fgCol, GetLineTransparency());
         glBegin(GL_LINES);
         Float_t aOff = 0.25*line_height;
         glVertex3f(x1+aOff, y1+aOff, z0);
         glVertex3f(x2-aOff, y1+aOff, z0);
         glVertex3f(x2-aOff, y1+aOff, z0);
         glVertex3f(x2-aOff, y2-aOff, z0);
         glEnd();
      }
   }

   glPopName();

   glPopMatrix();

   if (fDrawRefLine)
   {
      TGLVertex3 op = rnrCtx.RefCamera().WorldToViewport(fPointer);
      op[0] /= vp.Width();  op[1] /= vp.Height();

      Float_t fx = op[0] < fPosX ? 0.0f : (op[0] > fPosX + fDrawW ? 1.0f : 0.5f);
      Float_t fy = op[1] < fPosY-fDrawH ? 1.0f : (op[1] > fPosY ? 0.0f : 0.5f);

      if (fx != 0.5f || fy != 0.5f)
      {
         TGLUtil::ColorTransparency(bgCol, fTransparency);
         TGLUtil::LineWidth(2);
         glBegin(GL_LINES);
         glVertex3f(fPosX + fx*fDrawW, fPosY - fy*fDrawH, z3);
         glVertex3f(op[0], op[1], z3);
         glEnd();
      }
   }

   glPopMatrix();
   rnrCtx.ProjectionMatrixPop();

   glDepthRange(old_depth_range[0], old_depth_range[1]);
   glPopAttrib();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns transparency of annotation outline.
/// If annotation is selected enforce visibility of outline.

Char_t TGLAnnotation::GetLineTransparency() const
{
   if (fActive)
      return TMath::Min(70, fTransparency);
   else
      return fTransparency;
}

////////////////////////////////////////////////////////////////////////////////
/// Show the annotation editor.

void TGLAnnotation::MakeEditor()
{
   if (fMainFrame == 0)
   {
      fMainFrame = new TGMainFrame(gClient->GetRoot(), 1000, 1000);
      fMainFrame->SetWindowName("Annotation Editor");

      TGVerticalFrame* vf = new TGVerticalFrame(fMainFrame);

      fTextEdit = new TGTextEdit(vf,  1000, 1000, kSunkenFrame);
      vf->AddFrame(fTextEdit,  new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));

      TGHorizontalFrame* hf = new TGHorizontalFrame(vf);

      TGTextButton* btt1 = new TGTextButton(hf, "OK");
      hf->AddFrame(btt1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

      TGTextButton* btt2 = new TGTextButton(hf, "Cancel");
      hf->AddFrame(btt2, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

      btt1->Connect("Clicked()", "TGLAnnotation", this, "UpdateText()");
      btt2->Connect("Clicked()", "TGLAnnotation", this, "CloseEditor()");

      vf->AddFrame(hf, new TGLayoutHints(kLHintsBottom | kLHintsRight | kLHintsExpandX, 2, 2, 5, 1));

      fMainFrame->AddFrame(vf,  new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));
      fMainFrame->SetCleanup(kDeepCleanup);
      fMainFrame->MapSubwindows();
   }

   TGText *tgt = new TGText();
   tgt->LoadBuffer(fText.Data());
   fTextEdit->SetText(tgt);

   Int_t nrow = tgt->RowCount();
   Int_t h = nrow*20;
   Int_t w = fTextEdit->ReturnLongestLineWidth();
   fMainFrame->Resize(TMath::Max(100, w+30), TMath::Max(100, h+40));

   fMainFrame->Layout();
   fMainFrame->MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Close the annotation editor.

void TGLAnnotation::CloseEditor()
{
   fMainFrame->UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Modify the annotation text from the text-edit widget.

void TGLAnnotation::UpdateText()
{
   fText = fTextEdit->GetText()->AsString();
   fMainFrame->UnmapWindow();
   fParent->RequestDraw();
}
