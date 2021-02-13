// @(#)root/gl:$Id$
// Author: Matevz Tadel   25/09/2006

#include <cstring>

#include "TGLPShapeObjEditor.h"
#include "TGLPShapeObj.h"
#include "TGedEditor.h"

#include "TG3DLine.h"
#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TString.h"
#include "TGLabel.h"
#include "TGSlider.h"
#include "TGNumberEntry.h"
#include "TROOT.h"
#include "TVirtualMutex.h"

#include "TVirtualX.h"
#include "TGLViewer.h"
#include "TGLUtil.h"
#include "TGLPhysicalShape.h"
#include "TGLWidget.h"
#include "TGLIncludes.h"

#include "Buttons.h"

/** \class TGLPShapeObjEditor
\ingroup opengl
GUI editor for TGLPShapeObj.
*/

ClassImp(TGLPShapeObjEditor);

enum EGeometry {
      kCenterX,
      kCenterY,
      kCenterZ,
      kScaleX,
      kScaleY,
      kScaleZ,
      kTot
};

enum EApplyButtonIds {
      kTBcp,
      kTBcpm,
      kTBda,
      kTBa,
      kTBaf,
      kTBEndOfList
};

enum EGLEditorIdent {
      kCPa = kTBEndOfList + 1,
      kCPd, kCPs, kCPe,
      kHSr, kHSg, kHSb,
      kHSa, kHSs, kHSe,
      kNExc, kNEyc, kNEzc,
      kNExs, kNEys, kNEzs,
      kNExp, kNEyp, kNEzp,
      kNEat
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor of TGLPhysicalShape editor GUI.

TGLPShapeObjEditor::TGLPShapeObjEditor(const TGWindow *p,  Int_t width, Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p,  width, height, options | kVerticalFrame, back),
     fLb(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 2, 3, 3), //button
     fLe(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 0, 0, 3, 3), //entries
     fLl(kLHintsLeft, 0, 8, 6, 0), // labels
     fLs(kLHintsTop | kLHintsCenterX, 2, 2, 0, 0),  ///sliders
     fGeoFrame(0),fGeoApplyButton(0),
     fColorFrame(0),
     fRedSlider(0), fGreenSlider(0), fBlueSlider(0), fAlphaSlider(0), fShineSlider(0),
     fColorApplyButton(0), fColorApplyFamily(0),
     fRGBA(),
     fPShapeObj(0)
{
   fRGBA[12] = fRGBA[13] = fRGBA[14] = 0.0f;
   fRGBA[15] = 1.0f;
   fRGBA[16] = 60.0f;

   CreateColorControls();
   CreateGeoControls();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy color editor GUI component.
/// Done automatically.

TGLPShapeObjEditor::~TGLPShapeObjEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Shape has changed.
/// Check if set to zero and make sure we're no longer in editor.

void TGLPShapeObjEditor::SetPShape(TGLPhysicalShape * shape)
{
   TGLPShapeRef::SetPShape(shape);
   if (shape == 0 && fGedEditor->GetModel() == fPShapeObj)
      fGedEditor->SetModel(fGedEditor->GetPad(), fPShapeObj->fViewer, kButton1Down);
}

////////////////////////////////////////////////////////////////////////////////
/// Shape has been modified.
/// Update editor if we're still shown. Otherwise unref.

void TGLPShapeObjEditor::PShapeModified()
{
   if (fGedEditor->GetModel() == fPShapeObj)
      fGedEditor->SetModel(fGedEditor->GetPad(), fPShapeObj, kButton1Down);
   else
      SetPShape(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets model or disables/hides viewer.

void TGLPShapeObjEditor::SetModel(TObject* obj)
{
   fPShapeObj = 0;

   fPShapeObj = static_cast<TGLPShapeObj *>(obj);
   SetPShape(fPShapeObj->fPShape);

   SetRGBA(fPShapeObj->fPShape->Color());
   SetCenter(fPShapeObj->fPShape->GetTranslation().CArr());
   SetScale(fPShapeObj->fPShape->GetScale().CArr());
   fGeoApplyButton->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// Set internal center data from 3 component 'c'.

void TGLPShapeObjEditor::SetCenter(const Double_t *c)
{
   fGeomData[kCenterX]->SetNumber(c[0]);
   fGeomData[kCenterY]->SetNumber(c[1]);
   fGeomData[kCenterZ]->SetNumber(c[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Set internal scale data from 3 component 'c'.

void TGLPShapeObjEditor::SetScale(const Double_t *s)
{
   fGeomData[kScaleX]->SetNumber(s[0]);
   fGeomData[kScaleY]->SetNumber(s[1]);
   fGeomData[kScaleZ]->SetNumber(s[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Process 'Apply' - update the viewer object from GUI.

void TGLPShapeObjEditor::DoGeoButton()
{
   TGLVertex3 trans;
   TGLVector3 scale;
   GetObjectData(trans.Arr(), scale.Arr());
   if (fPShape) {
      fPShape->SetTranslation(trans);
      fPShape->Scale(scale);
   }
   fPShapeObj->fViewer->RequestDraw();
   fGeoApplyButton->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// Extract the GUI object data, return center in 3 component 'center'
/// scale in 3 component 'scale'.

void TGLPShapeObjEditor::GetObjectData(Double_t *center, Double_t *scale)
{
   center[0] = fGeomData[kCenterX]->GetNumber();
   center[1] = fGeomData[kCenterY]->GetNumber();
   center[2] = fGeomData[kCenterZ]->GetNumber();
   scale[0] = fGeomData[kScaleX]->GetNumber();
   scale[1] = fGeomData[kScaleY]->GetNumber();
   scale[2] = fGeomData[kScaleZ]->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Process setting of value in edit box - activate 'Apply' button.

void TGLPShapeObjEditor::GeoValueSet(Long_t)
{
   if (fGeoApplyButton->GetState() != kButtonUp)
       fGeoApplyButton->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Create GUI for setting scale and position.

void TGLPShapeObjEditor::CreateGeoControls()
{
   fGeoFrame = CreateEditorTabSubFrame("Geometry");

   TGLabel *label=0;

   // Postion container
   TGGroupFrame* container = new TGGroupFrame(fGeoFrame, "Object position:");
   container->SetTitlePos(TGGroupFrame::kLeft);
   fGeoFrame->AddFrame(container, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 8, 8, 3, 3));//-
   TGLayoutHints lh =  TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 0, 0, 0, 0);

   TGHorizontalFrame* hf;

   hf = new TGHorizontalFrame(container);
   label = new TGLabel(hf, "X:");
   hf->AddFrame(label, new TGLayoutHints(fLl));
   fGeomData[kCenterX] = new TGNumberEntry(hf, 0.0, 8, kNExc);
   hf->AddFrame(fGeomData[kCenterX], new TGLayoutHints(fLe));
   fGeomData[kCenterX]->Connect("ValueSet(Long_t)", "TGLPShapeObjEditor",
                                this, "GeoValueSet(Long_t)");
   container->AddFrame(hf, new TGLayoutHints(lh));

   hf = new TGHorizontalFrame(container);
   label = new TGLabel(hf, "Y:");
   hf->AddFrame(label, new TGLayoutHints(fLl));
   fGeomData[kCenterY] = new TGNumberEntry(hf, 0.0, 8, kNEyc);
   hf->AddFrame(fGeomData[kCenterY], new TGLayoutHints(fLe));
   fGeomData[kCenterY]->Connect("ValueSet(Long_t)", "TGLPShapeObjEditor",
                                this, "GeoValueSet(Long_t)");
   container->AddFrame(hf, new TGLayoutHints(lh));

   hf = new TGHorizontalFrame(container);
   hf->AddFrame(new TGLabel(hf, "Z:"), new TGLayoutHints(fLl));
   fGeomData[kCenterZ] = new TGNumberEntry(hf, 1.0, 8, kNEzc);
   hf->AddFrame(fGeomData[kCenterZ], new TGLayoutHints(fLe));
   fGeomData[kCenterZ]->Connect("ValueSet(Long_t)", "TGLPShapeObjEditor",
                                this, "GeoValueSet(Long_t)");
   container->AddFrame(hf, new TGLayoutHints(lh));

   // Scale container
   TGGroupFrame* osf = new TGGroupFrame(fGeoFrame, "Object scale:", kLHintsTop | kLHintsCenterX);
   osf->SetTitlePos(TGGroupFrame::kLeft);
   fGeoFrame->AddFrame(osf, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 8, 8, 3, 3));

   hf = new TGHorizontalFrame(osf);
   hf->AddFrame(new TGLabel(hf, "X:"),new TGLayoutHints(fLl));
   fGeomData[kScaleX] = new TGNumberEntry(hf, 1.0, 5, kNExs);
   hf->AddFrame(fGeomData[kScaleX], new TGLayoutHints(fLe));
   fGeomData[kScaleX]->Connect("ValueSet(Long_t)", "TGLPShapeObjEditor",
                               this, "GeoValueSet(Long_t)");
   osf->AddFrame(hf, new TGLayoutHints(lh));

   hf = new TGHorizontalFrame(osf);
   hf->AddFrame(new TGLabel(hf, "Y:"),new TGLayoutHints(fLl));
   fGeomData[kScaleY] = new TGNumberEntry(hf, 1.0, 5, kNEys);
   hf->AddFrame(fGeomData[kScaleY], new TGLayoutHints(fLe));
   fGeomData[kScaleY]->Connect("ValueSet(Long_t)", "TGLPShapeObjEditor",
                               this, "GeoValueSet(Long_t)");
   osf->AddFrame(hf, new TGLayoutHints(lh));

   hf = new TGHorizontalFrame(osf);
   hf->AddFrame(new TGLabel(hf, "Z:"),new TGLayoutHints(fLl));
   fGeomData[kScaleZ] = new TGNumberEntry(hf, 1.0, 5, kNEzs);
   hf->AddFrame(fGeomData[kScaleZ], new TGLayoutHints(fLe));
   fGeomData[kScaleZ]->Connect("ValueSet(Long_t)", "TGLPShapeObjEditor",
                               this, "GeoValueSet(Long_t)");
   osf->AddFrame(hf, new TGLayoutHints(lh));

   hf = new TGHorizontalFrame(osf);
   fGeomData[kScaleX]->SetLimits(TGNumberFormat::kNELLimitMin, 0.1);
   fGeomData[kScaleY]->SetLimits(TGNumberFormat::kNELLimitMin, 0.1);
   fGeomData[kScaleZ]->SetLimits(TGNumberFormat::kNELLimitMin, 0.1);
   osf->AddFrame(hf, new TGLayoutHints(lh));

   // Modify button
   fGeoApplyButton = new TGTextButton(fGeoFrame, "Modify object");
   fGeoFrame->AddFrame(fGeoApplyButton, new TGLayoutHints(fLb));
   fGeoApplyButton->SetState(kButtonDisabled);
   fGeoApplyButton->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoGeoButton()");
}

////////////////////////////////////////////////////////////////////////////////
/// Set color sliders from 17 component 'rgba'.

void TGLPShapeObjEditor::SetRGBA(const Float_t *rgba)
{
   fColorApplyButton->SetState(kButtonDisabled);
   fColorApplyFamily->SetState(kButtonDisabled);

   for (Int_t i = 0; i < 17; ++i) fRGBA[i] = rgba[i];

   fRedSlider->SetPosition(Int_t(fRGBA[fLMode * 4] * 100));
   fGreenSlider->SetPosition(Int_t(fRGBA[fLMode * 4 + 1] * 100));
   fBlueSlider->SetPosition(Int_t(fRGBA[fLMode * 4 + 2] * 100));
   fShineSlider->SetPosition(Int_t(fRGBA[16]));

   DrawSphere();
}

////////////////////////////////////////////////////////////////////////////////
/// Process slider movement.

void TGLPShapeObjEditor::DoColorSlider(Int_t val)
{
   TGSlider *frm = (TGSlider *)gTQSender;

   if (frm) {
      Int_t wid = frm->WidgetId();

      switch (wid) {
      case kHSr:
         fRGBA[fLMode * 4] = val / 100.f;
         break;
      case kHSg:
         fRGBA[fLMode * 4 + 1] = val / 100.f;
         break;
      case kHSb:
         fRGBA[fLMode * 4 + 2] = val / 100.f;
         break;
      case kHSa:
         fRGBA[fLMode * 4 + 3] = val / 100.f;
         break;
      case kHSs:
         fRGBA[16] = val;
         break;
      }

      fColorApplyButton->SetState(kButtonUp);
      fColorApplyFamily->SetState(kButtonUp);
      DrawSphere();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process button action.

void TGLPShapeObjEditor::DoColorButton()
{
   TGButton *btn = (TGButton *) gTQSender;
   Int_t id = btn->WidgetId();

   switch (id) {
   case kCPd:
      fLightTypes[fLMode]->SetState(kButtonUp);
      fLMode = kDiffuse;
      SetColorSlidersPos();
      break;
   case kCPa:
      fLightTypes[fLMode]->SetState(kButtonUp);
      fLMode = kAmbient;
      SetColorSlidersPos();
      break;
   case kCPs:
      fLightTypes[fLMode]->SetState(kButtonUp);
      fLMode = kSpecular;
      SetColorSlidersPos();
      break;
   case kCPe:
      fLightTypes[fLMode]->SetState(kButtonUp);
      fLMode = kEmission;
      SetColorSlidersPos();
      break;
   case kTBa:
      fColorApplyButton->SetState(kButtonDisabled);
      fColorApplyFamily->SetState(kButtonDisabled);
      if (fPShape) {
         fPShape->SetColor(GetRGBA());
      }
      fPShapeObj->fViewer->RequestDraw();
      break;
   case kTBaf:
      fColorApplyButton->SetState(kButtonDisabled);
      fColorApplyFamily->SetState(kButtonDisabled);
      if (fPShape) {
         fPShape->SetColorOnFamily(GetRGBA());
      }
      fPShapeObj->fViewer->RequestDraw();
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create Diffuse/Ambient/Specular/Emissive radio buttons and sub-frames.

void TGLPShapeObjEditor::CreateColorRadioButtons()
{
   TGGroupFrame *partFrame = new TGGroupFrame(fColorFrame, "Color components:", kLHintsTop | kLHintsCenterX);
   fColorFrame->AddFrame(partFrame, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 0, 2, 2));

   partFrame->SetTitlePos(TGGroupFrame::kLeft);
   TGMatrixLayout *ml = new TGMatrixLayout(partFrame, 0, 1, 10);
   partFrame->SetLayoutManager(ml);

   // partFrame will delete the layout manager ml for us so don't add to fTrash
   fLightTypes[kDiffuse] = new TGRadioButton(partFrame, "Diffuse", kCPd);
   fLightTypes[kDiffuse]->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoColorButton()");
   fLightTypes[kDiffuse]->SetToolTipText("Diffuse component of color");
   partFrame->AddFrame(fLightTypes[kDiffuse]);

   fLightTypes[kAmbient] = new TGRadioButton(partFrame, "Ambient", kCPa);
   fLightTypes[kAmbient]->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoColorButton()");
   fLightTypes[kAmbient]->SetToolTipText("Ambient component of color");
   partFrame->AddFrame(fLightTypes[kAmbient]);

   fLightTypes[kSpecular] = new TGRadioButton(partFrame, "Specular", kCPs);
   fLightTypes[kSpecular]->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoColorButton()");
   fLightTypes[kSpecular]->SetToolTipText("Specular component of color");
   partFrame->AddFrame(fLightTypes[kSpecular]);

   fLightTypes[kEmission] = new TGRadioButton(partFrame, "Emissive", kCPe);
   fLightTypes[kEmission]->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoColorButton()");
   fLightTypes[kEmission]->SetToolTipText("Emissive component of color");
   partFrame->AddFrame(fLightTypes[kEmission]);

   fLMode = kDiffuse;
   fLightTypes[fLMode]->SetState(kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Create GUI for setting light color.

void TGLPShapeObjEditor::CreateColorSliders()
{
   UInt_t sw = 120; //fColorFrame->GetDefalutWidth();,

   // Create Red/Green/BlueAlpha/Shine sliders
   fColorFrame->AddFrame(new TGLabel(fColorFrame, "Red :"), new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 0, 0, 0));
   fRedSlider = new TGHSlider(fColorFrame, sw, kSlider1 | kScaleBoth, kHSr);
   fRedSlider->Connect("PositionChanged(Int_t)", "TGLPShapeObjEditor", this, "DoColorSlider(Int_t)");
   fRedSlider->SetRange(0, 100);
   fRedSlider->SetPosition(Int_t(fRGBA[0] * 100));
   fColorFrame->AddFrame(fRedSlider, new TGLayoutHints(fLs));


   fColorFrame->AddFrame(new TGLabel(fColorFrame, "Green :"), new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 0, 0, 0));
   fGreenSlider = new TGHSlider(fColorFrame, sw, kSlider1 | kScaleBoth, kHSg);
   fGreenSlider->Connect("PositionChanged(Int_t)", "TGLPShapeObjEditor", this, "DoColorSlider(Int_t)");
   fGreenSlider->SetRange(0, 100);
   fGreenSlider->SetPosition(Int_t(fRGBA[1] * 100));
   fColorFrame->AddFrame(fGreenSlider, new TGLayoutHints(fLs));


   fColorFrame->AddFrame(new TGLabel(fColorFrame, "Blue :"), new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 0, 0, 0));
   fBlueSlider = new TGHSlider(fColorFrame, sw, kSlider1 | kScaleBoth, kHSb);
   fBlueSlider->Connect("PositionChanged(Int_t)", "TGLPShapeObjEditor", this, "DoColorSlider(Int_t)");
   fBlueSlider->SetRange(0, 100);
   fBlueSlider->SetPosition(Int_t(fRGBA[2] * 100));
   fColorFrame->AddFrame(fBlueSlider, new TGLayoutHints(fLs));

   fColorFrame->AddFrame(new TGLabel(fColorFrame, "Shine :"), new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 0, 0, 0));
   fShineSlider = new TGHSlider(fColorFrame, sw, kSlider1 | kScaleBoth, kHSs);
   fShineSlider->Connect("PositionChanged(Int_t)", "TGLPShapeObjEditor", this, "DoColorSlider(Int_t)");
   fShineSlider->SetRange(0, 128);
   fColorFrame->AddFrame(fShineSlider, new TGLayoutHints(fLs));
}

////////////////////////////////////////////////////////////////////////////////
/// Update GUI sliders from internal data.

void TGLPShapeObjEditor::SetColorSlidersPos()
{
   fRedSlider->SetPosition(Int_t(fRGBA[fLMode * 4] * 100));
   fGreenSlider->SetPosition(Int_t(fRGBA[fLMode * 4 + 1] * 100));
   fBlueSlider->SetPosition(Int_t(fRGBA[fLMode * 4 + 2] * 100));
   //   fAlphaSlider->SetPosition(Int_t(fRGBA[fLMode * 4 + 3] * 100));

   if (fRGBA[16] >= 0.f)
      fShineSlider->SetPosition(Int_t(fRGBA[16]));
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw widget. Render sphere and pass to base-class.

void TGLPShapeObjEditor::DoRedraw()
{
   DrawSphere();
   TGedFrame::DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////

namespace {
   GLUquadric *GetQuadric()
   {
      // GLU quadric.

      static struct Init {
         Init()
         {
            fQuad = gluNewQuadric();
            if (!fQuad) {
               Error("GetQuadric::Init", "could not create quadric object");
            } else {
               gluQuadricOrientation(fQuad, (GLenum)GLU_OUTSIDE);
               gluQuadricDrawStyle(fQuad,   (GLenum)GLU_FILL);
               gluQuadricNormals(fQuad,     (GLenum)GLU_FLAT);
            }
         }
         ~Init()
         {
            if(fQuad)
               gluDeleteQuadric(fQuad);
         }
         GLUquadric *fQuad;
      }singleton;

      return singleton.fQuad;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Draw local sphere reflecting current color options.

void TGLPShapeObjEditor::DrawSphere()const
{
   if (!gVirtualX->IsCmdThread()) {
      gROOT->ProcessLineFast(Form("((TGLPShapeObjEditor *)0x%lx)->DrawSphere()", (ULong_t)this));
      return;
   }

   R__LOCKGUARD(gROOTMutex);

   fMatView->MakeCurrent();
   glViewport(0, 0, fMatView->GetWidth(), fMatView->GetHeight());
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-0.5, 0.5, -0.5, 0.5, 1., 10.);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   Float_t ligPos[] = {0.f, 0.f, 0.f, 1.f};
   glLightfv(GL_LIGHT0, GL_POSITION, ligPos);
   glTranslated(0., 0., -3.);

   const Float_t whiteColor[] = {1.f, 1.f, 1.f, 1.f};
   const Float_t nullColor[] = {0.f, 0.f, 0.f, 1.f};

   if (fRGBA[16] < 0.f) {
      glLightfv(GL_LIGHT0, GL_DIFFUSE, fRGBA);
      glLightfv(GL_LIGHT0, GL_AMBIENT, fRGBA + 4);
      glLightfv(GL_LIGHT0, GL_SPECULAR, fRGBA + 8);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, whiteColor);
      glMaterialfv(GL_FRONT, GL_AMBIENT, nullColor);
      glMaterialfv(GL_FRONT, GL_SPECULAR, whiteColor);
      glMaterialfv(GL_FRONT, GL_EMISSION, nullColor);
      glMaterialf(GL_FRONT, GL_SHININESS, 60.f);
   } else {
      glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteColor);
      glLightfv(GL_LIGHT0, GL_AMBIENT, nullColor);
      glLightfv(GL_LIGHT0, GL_SPECULAR, whiteColor);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, fRGBA);
      glMaterialfv(GL_FRONT, GL_AMBIENT, fRGBA + 4);
      glMaterialfv(GL_FRONT, GL_SPECULAR, fRGBA + 8);
      glMaterialfv(GL_FRONT, GL_EMISSION, fRGBA + 12);
      glMaterialf(GL_FRONT, GL_SHININESS, fRGBA[16]);
   }

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   GLUquadric * quad = GetQuadric();
   if (quad) {
      glRotated(-90., 1., 0., 0.);
      gluSphere(quad, 1., 100, 100);
   }
   glDisable(GL_BLEND);

   fMatView->SwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////
/// Create widgets to chose colors component and its RGBA values on fGedEditor
/// model or family it belongs to.

void TGLPShapeObjEditor::CreateColorControls()
{
   fColorFrame = this;

   fMatView = TGLWidget::Create(fColorFrame, kFALSE, kTRUE, 0, 120, 120);
   fColorFrame->AddFrame(fMatView, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 0, 2, 2));

   CreateColorRadioButtons();

   CreateColorSliders();

   //apply button creation
   fColorApplyButton = new TGTextButton(fColorFrame, "Apply", kTBa);
   fColorFrame->AddFrame(fColorApplyButton, new TGLayoutHints(fLb));
   fColorApplyButton->SetState(kButtonDisabled);
   fColorApplyButton->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoColorButton()");
   //apply to family button creation
   fColorApplyFamily = new TGTextButton(fColorFrame, "Apply to family", kTBaf);
   fColorFrame->AddFrame(fColorApplyFamily, new TGLayoutHints(fLb));
   fColorApplyFamily->SetState(kButtonDisabled);
   fColorApplyFamily->Connect("Pressed()", "TGLPShapeObjEditor", this, "DoColorButton()");
}
