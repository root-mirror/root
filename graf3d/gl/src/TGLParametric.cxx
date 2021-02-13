// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  26/01/2007

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include <cctype>

#ifdef WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include "TVirtualX.h"
#include "TString.h"
#include "TROOT.h"
#include "TH3.h"

#include "TGLPlotCamera.h"
#include "TGLParametric.h"
#include "TGLIncludes.h"
#include "KeySymbols.h"
#include "Buttons.h"
#include "TMath.h"

namespace
{

   /////////////////////////////////////////////////////////////////////////////
   ///User defines equations using names 'u' and 'v' for
   ///parameters. But TF2 works with 'x' and 'y'. So,
   ///find 'u' and 'v' (which are not parts of other names)
   ///and replace them with 'x' and 'y' correspondingly.

   void ReplaceUVNames(TString &equation)
   {
      using namespace std;
      const Ssiz_t len = equation.Length();
      //TF2 requires 'y' in formula.
      //'v' <=> 'y', so if none 'v' was found, I'll append "+0*y" to the equation.
      Int_t vFound = 0;

      for (Ssiz_t i = 0; i < len;) {
         const char c = equation[i];
         if (!isalpha(c)) {
            ++i;
            continue;
         } else{
            ++i;
            if (c == 'u' || c == 'v') {
               //1. This 'u' or 'v' is the last symbol in a string or
               //2. After this 'u' or 'v' symbol, which cannot be part of longer name.
               if (i == len || (!isalpha(equation[i]) && !isdigit(equation[i]) && equation[i] != '_')) {
                  //Replace 'u' with 'x' or 'v' with 'y'.
                  equation[i - 1] = c == 'u' ? 'x' : (++vFound, 'y');
               } else {
                  //This 'u' or 'v' is the beginning of some longer name.
                  //Skip the remaining part of this name.
                  while (i < len && (isalpha(equation[i]) || isdigit(equation[i]) || equation[i] == '_'))
                     ++i;
               }
            } else {
               while (i < len && (isalpha(equation[i]) || isdigit(equation[i]) || equation[i] == '_'))
                  ++i;
            }
         }
      }

      if (!vFound)
         equation += "+0*y";
   }

}

/** \class TGLParametricEquation
\ingroup opengl
A parametric surface is a surface defined by a parametric equation, involving
two parameters (u, v):

~~~ {.cpp}
S(u, v) = (x(u, v), y(u, v), z(u, v)).
~~~

For example, "limpet torus" surface can be defined as:

~~~ {.cpp}
    x = cos(u) / (sqrt(2) + sin(v))
    y = sin(u) / (sqrt(2) + sin(v))
    z = 1 / (sqrt(2) + cos(v)),
~~~

where -pi <= u <= pi, -pi <= v <= pi.

~~~ {.cpp}
 TGLParametricEquation * eq =
    new TGLParametricEquation("Limpet_torus", "cos(u) / (sqrt(2.) + sin(v))",
                              "sin(u) / (sqrt(2.) + sin(v))",
                              "1 / (sqrt(2) + cos(v))");
~~~

 `$ROOTSYS/tutorials/gl/glparametric.C` contains more examples.

 Parametric equations can be specified:
  - 1. by string expressions, as with TF2, but with 'u' instead of 'x' and
       'v' instead of 'y'.
  - 2. by function - see ParametricEquation_t declaration.

*/

ClassImp(TGLParametricEquation);

////////////////////////////////////////////////////////////////////////////////
///Surface is defined by three strings.
///ROOT does not use exceptions in ctors,
///so, I have to use MakeZombie to let
///external user know about errors.

TGLParametricEquation::TGLParametricEquation(const TString &name, const TString &xFun, const TString &yFun,
                             const TString &zFun, Double_t uMin, Double_t uMax,
                             Double_t vMin, Double_t vMax)
                  : TNamed(name, name),
                    fEquation(0),
                    fURange(uMin, uMax),
                    fVRange(vMin, vMax),
                    fConstrained(kFALSE),
                    fModified(kFALSE)
{
   if (!xFun.Length() || !yFun.Length() || !zFun.Length()) {
      Error("TGLParametricEquation", "One of string expressions is empty");
      MakeZombie();
      return;
   }

   TString equation(xFun);
   equation.ToLower();
   ReplaceUVNames(equation);
   fXEquation.reset(new TF2(name + "xEquation", equation.Data(), uMin, uMax, vMin, vMax));
   //Formula was incorrect.
   if (fXEquation->IsZombie()) {
      MakeZombie();
      return;
   }

   equation = yFun;
   equation.ToLower();
   ReplaceUVNames(equation);
   fYEquation.reset(new TF2(name + "yEquation", equation.Data(), uMin, uMax, vMin, vMax));
   //Formula was incorrect.
   if (fYEquation->IsZombie()) {
      MakeZombie();
      return;
   }

   equation = zFun;
   equation.ToLower();
   ReplaceUVNames(equation);
   fZEquation.reset(new TF2(name + "zEquation", equation.Data(), uMin, uMax, vMin, vMax));
   //Formula was incorrect.
   if (fZEquation->IsZombie())
      MakeZombie();
}

////////////////////////////////////////////////////////////////////////////////
///Surface defined by user's function (see ParametricEquation_t declaration in TGLParametricEquation.h)

TGLParametricEquation::TGLParametricEquation(const TString &name, ParametricEquation_t equation,
                             Double_t uMin, Double_t uMax, Double_t vMin, Double_t vMax)
                  : TNamed(name, name),
                    fEquation(equation),
                    fURange(uMin, uMax),
                    fVRange(vMin, vMax),
                    fConstrained(kFALSE),
                    fModified(kFALSE)
{
   if (!fEquation) {
      Error("TGLParametricEquation", "Function ptr is null");
      MakeZombie();
   }
}

////////////////////////////////////////////////////////////////////////////////
///[uMin, uMax]

Rgl::Range_t TGLParametricEquation::GetURange()const
{
   return fURange;
}

////////////////////////////////////////////////////////////////////////////////
///[vMin, vMax]

Rgl::Range_t TGLParametricEquation::GetVRange()const
{
   return fVRange;
}

////////////////////////////////////////////////////////////////////////////////
///Check is constrained.

Bool_t TGLParametricEquation::IsConstrained()const
{
   return fConstrained;
}

////////////////////////////////////////////////////////////////////////////////
///Set constrained.

void TGLParametricEquation::SetConstrained(Bool_t c)
{
   fConstrained = c;
}

////////////////////////////////////////////////////////////////////////////////
///Something was changed in parametric equation (or constrained option was changed).

Bool_t TGLParametricEquation::IsModified()const
{
   return fModified;
}

////////////////////////////////////////////////////////////////////////////////
///Set modified.

void TGLParametricEquation::SetModified(Bool_t m)
{
   fModified = m;
}

////////////////////////////////////////////////////////////////////////////////
///Calculate vertex.

void TGLParametricEquation::EvalVertex(TGLVertex3 &newVertex, Double_t u, Double_t v)const
{
   if (fEquation)
      return fEquation(newVertex, u, v);

   if (IsZombie())
      return;

   newVertex.X() = fXEquation->Eval(u, v);
   newVertex.Y() = fYEquation->Eval(u, v);
   newVertex.Z() = fZEquation->Eval(u, v);
}

////////////////////////////////////////////////////////////////////////////////
///Check, if parametric surface is under cursor.

Int_t TGLParametricEquation::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (fPainter.get())
      return fPainter->DistancetoPrimitive(px, py);
   return 9999;
}

////////////////////////////////////////////////////////////////////////////////
///Pass event to painter.

void TGLParametricEquation::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (fPainter.get())
      return fPainter->ExecuteEvent(event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
///No object info yet.

char *TGLParametricEquation::GetObjectInfo(Int_t /*px*/, Int_t /*py*/) const
{
   static char mess[] = { "parametric surface" };
   return mess;
}

////////////////////////////////////////////////////////////////////////////////
///Delegate paint.

void TGLParametricEquation::Paint(Option_t * /*option*/)
{
   if (!fPainter.get())
      fPainter.reset(new TGLHistPainter(this));
   fPainter->Paint("dummyoption");
}

/** \class TGLParametricPlot
\ingroup opengl
*/

ClassImp(TGLParametricPlot);

////////////////////////////////////////////////////////////////////////////////
///Constructor.

TGLParametricPlot::TGLParametricPlot(TGLParametricEquation *eq,
                                     TGLPlotCamera *camera)
                      : TGLPlotPainter(camera),
                        fMeshSize(90),
                        fShowMesh(kFALSE),
                        fColorScheme(4),
                        fEquation(eq)
{
   fXAxis = &fCartesianXAxis;
   fYAxis = &fCartesianYAxis;
   fZAxis = &fCartesianZAxis;

   fCoord = &fCartesianCoord;

   InitGeometry();
   InitColors();
}

////////////////////////////////////////////////////////////////////////////////
///Build mesh. The surface is 'immutable':
///the only reason to rebuild it - the change in size or
///if one of equations contain reference to TF2 function, whose
///parameters were changed.

Bool_t TGLParametricPlot::InitGeometry()
{
   // const Bool_t constrained = kTRUE;//fEquation->IsConstrained();

   if (fMeshSize * fMeshSize != (Int_t)fMesh.size() || fEquation->IsModified()) {
      if (fEquation->IsZombie())
         return kFALSE;

      fEquation->SetModified(kFALSE);

      fMesh.resize(fMeshSize * fMeshSize);
      fMesh.SetRowLen(fMeshSize);

      const Rgl::Range_t uRange(fEquation->GetURange());
      const Rgl::Range_t vRange(fEquation->GetVRange());

      const Double_t dU = (uRange.second - uRange.first) / (fMeshSize - 1);
      const Double_t dV = (vRange.second - vRange.first) / (fMeshSize - 1);
      const Double_t dd = 0.001;
      Double_t u = uRange.first;

      TGLVertex3 min;
      fEquation->EvalVertex(min, uRange.first, vRange.first);
      TGLVertex3 max(min), newVert, v1, v2;
      using namespace TMath;

      for (Int_t i = 0; i < fMeshSize; ++i) {
         Double_t v = vRange.first;
         for (Int_t j = 0; j < fMeshSize; ++j) {
            fEquation->EvalVertex(newVert, u, v);
            min.X() = Min(min.X(), newVert.X());
            max.X() = Max(max.X(), newVert.X());
            min.Y() = Min(min.Y(), newVert.Y());
            max.Y() = Max(max.Y(), newVert.Y());
            min.Z() = Min(min.Z(), newVert.Z());
            max.Z() = Max(max.Z(), newVert.Z());

            fMesh[i][j].fPos = newVert;

            v += dV;
         }
         u += dU;
      }

      TH3F hist("tmp", "tmp", 2, -1., 1., 2, -1., 1., 2, -1., 1.);
      hist.SetDirectory(0);
      //TAxis has a lot of attributes, defaults, set by ctor,
      //are not enough to be correctly painted by TGaxis object.
      //To simplify their initialization - I use temporary histogram.
      hist.GetXaxis()->Copy(fCartesianXAxis);
      hist.GetYaxis()->Copy(fCartesianYAxis);
      hist.GetZaxis()->Copy(fCartesianZAxis);

      fCartesianXAxis.Set(fMeshSize, min.X(), max.X());
      fCartesianXAxis.SetTitle("x");//it's lost when copying from temp. hist.
      fCartesianYAxis.Set(fMeshSize, min.Y(), max.Y());
      fCartesianYAxis.SetTitle("y");
      fCartesianZAxis.Set(fMeshSize, min.Z(), max.Z());
      fCartesianZAxis.SetTitle("z");

      if (!fCoord->SetRanges(&fCartesianXAxis, &fCartesianYAxis, &fCartesianZAxis))
         return kFALSE;

      for (Int_t i = 0; i < fMeshSize; ++i) {
         for (Int_t j = 0; j < fMeshSize; ++j) {
            TGLVertex3 &ver = fMesh[i][j].fPos;
            ver.X() *= fCoord->GetXScale(), ver.Y() *= fCoord->GetYScale(), ver.Z() *= fCoord->GetZScale();
         }
      }

      u = uRange.first;
      for (Int_t i = 0; i < fMeshSize; ++i) {
         Double_t v = vRange.first;
         for (Int_t j = 0; j < fMeshSize; ++j) {
            TGLVertex3 &ver = fMesh[i][j].fPos;
            fEquation->EvalVertex(v1, u + dd, v);
            fEquation->EvalVertex(v2, u, v + dd);
            v1.X() *= fCoord->GetXScale(), v1.Y() *= fCoord->GetYScale(), v1.Z() *= fCoord->GetZScale();
            v2.X() *= fCoord->GetXScale(), v2.Y() *= fCoord->GetYScale(), v2.Z() *= fCoord->GetZScale();
            Normal2Plane(ver.CArr(), v1.CArr(), v2.CArr(), fMesh[i][j].fNormal.Arr());
            v += dV;
         }
         u += dU;
      }

      fBackBox.SetPlotBox(fCoord->GetXRangeScaled(),
                          fCoord->GetYRangeScaled(),
                          fCoord->GetZRangeScaled());
      if (fCamera) fCamera->SetViewVolume(fBackBox.Get3DBox());
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///User clicks right mouse button (in a pad).

void TGLParametricPlot::StartPan(Int_t px, Int_t py)
{
   fMousePosition.fX = px;
   fMousePosition.fY = fCamera->GetHeight() - py;
   fCamera->StartPan(px, py);
   fBoxCut.StartMovement(px, fCamera->GetHeight() - py);
}

////////////////////////////////////////////////////////////////////////////////
///User's moving mouse cursor, with middle mouse button pressed (for pad).
///Calculate 3d shift related to 2d mouse movement.

void TGLParametricPlot::Pan(Int_t px, Int_t py)
{
   if (fSelectedPart) {
      SaveModelviewMatrix();
      SaveProjectionMatrix();

      fCamera->SetCamera();
      fCamera->Apply(fPadPhi, fPadTheta);

      if (fBoxCut.IsActive() && (fSelectedPart >= kXAxis && fSelectedPart <= kZAxis))
         fBoxCut.MoveBox(px, fCamera->GetHeight() - py, fSelectedPart);
      else
         fCamera->Pan(px, py);

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   }

   fUpdateSelection = kTRUE;//
}

////////////////////////////////////////////////////////////////////////////////
///No object info yet.

char *TGLParametricPlot::GetPlotInfo(Int_t /*px*/, Int_t /*py*/)
{
   static char mess[] = { "parametric surface" };
   return mess;
}

////////////////////////////////////////////////////////////////////////////////
///No additional options for parametric surfaces.

void TGLParametricPlot::AddOption(const TString &/*option*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///Change color/mesh size or switch on/off mesh/box cut.
///Left double click - remove box cut.

void TGLParametricPlot::ProcessEvent(Int_t event, Int_t /*px*/, Int_t py)
{
   if (event == kButton1Double && fBoxCut.IsActive()) {
      fBoxCut.TurnOnOff();
      if (!gVirtualX->IsCmdThread())
         gROOT->ProcessLineFast(Form("((TGLPlotPainter *)0x%lx)->Paint()", (ULong_t)this));
      else
         Paint();
   } else if (event == kKeyPress) {
      if (py == kKey_c || py == kKey_C) {
         if (fHighColor)
            Info("ProcessEvent", "Switch to true color to use box cut");
         else {
            fBoxCut.TurnOnOff();
            fUpdateSelection = kTRUE;
         }
      } else if (py == kKey_s || py == kKey_S) {
         fColorScheme == 20 ? fColorScheme = -1 : ++fColorScheme;
         InitColors();//color scheme was changed! recalculate vertices colors.
      } else if (py == kKey_w || py == kKey_W) {
         fShowMesh = !fShowMesh;
      } else if (py == kKey_l || py == kKey_L) {
         fMeshSize == kHigh ? fMeshSize = kLow : fMeshSize += 15;
         InitGeometry();
         InitColors();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///Initialize gl state.

void TGLParametricPlot::InitGL()const
{
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Initialize gl state.

void TGLParametricPlot::DeInitGL()const
{
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
}

////////////////////////////////////////////////////////////////////////////////
///Draw parametric surface.

void TGLParametricPlot::DrawPlot()const
{
   //Shift plot to point of origin.
   const Rgl::PlotTranslation trGuard(this);

   if (!fSelectionPass) {
      SetSurfaceColor();
      if (fShowMesh) {
         glEnable(GL_POLYGON_OFFSET_FILL);
         glPolygonOffset(1.f, 1.f);
      }
   } else {
      Rgl::ObjectIDToColor(fSelectionBase, fHighColor);
   }

   glBegin(GL_TRIANGLES);

   for (Int_t i = 0; i < fMeshSize - 1; ++i) {
      for (Int_t j = 0; j < fMeshSize - 1; ++j) {
         if (fBoxCut.IsActive()) {
            using TMath::Min;
            using TMath::Max;
            const Double_t xMin = Min(Min(fMesh[i][j].fPos.X(), fMesh[i + 1][j].fPos.X()), Min(fMesh[i][j + 1].fPos.X(), fMesh[i + 1][j + 1].fPos.X()));
            const Double_t xMax = Max(Max(fMesh[i][j].fPos.X(), fMesh[i + 1][j].fPos.X()), Max(fMesh[i][j + 1].fPos.X(), fMesh[i + 1][j + 1].fPos.X()));
            const Double_t yMin = Min(Min(fMesh[i][j].fPos.Y(), fMesh[i + 1][j].fPos.Y()), Min(fMesh[i][j + 1].fPos.Y(), fMesh[i + 1][j + 1].fPos.Y()));
            const Double_t yMax = Max(Max(fMesh[i][j].fPos.Y(), fMesh[i + 1][j].fPos.Y()), Max(fMesh[i][j + 1].fPos.Y(), fMesh[i + 1][j + 1].fPos.Y()));
            const Double_t zMin = Min(Min(fMesh[i][j].fPos.Z(), fMesh[i + 1][j].fPos.Z()), Min(fMesh[i][j + 1].fPos.Z(), fMesh[i + 1][j + 1].fPos.Z()));
            const Double_t zMax = Max(Max(fMesh[i][j].fPos.Z(), fMesh[i + 1][j].fPos.Z()), Max(fMesh[i][j + 1].fPos.Z(), fMesh[i + 1][j + 1].fPos.Z()));

            if (fBoxCut.IsInCut(xMin, xMax, yMin, yMax, zMin, zMax))
               continue;
         }

         glNormal3dv(fMesh[i + 1][j + 1].fNormal.CArr());
         if(fColorScheme != -1)
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fMesh[i + 1][j + 1].fRGBA);
         glVertex3dv(fMesh[i + 1][j + 1].fPos.CArr());

         glNormal3dv(fMesh[i][j + 1].fNormal.CArr());
         if(fColorScheme != -1)
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fMesh[i][j + 1].fRGBA);
         glVertex3dv(fMesh[i][j + 1].fPos.CArr());

         glNormal3dv(fMesh[i][j].fNormal.CArr());
         if(fColorScheme != -1)
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fMesh[i][j].fRGBA);
         glVertex3dv(fMesh[i][j].fPos.CArr());

         glNormal3dv(fMesh[i + 1][j].fNormal.CArr());
         if(fColorScheme != -1)
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fMesh[i + 1][j].fRGBA);
         glVertex3dv(fMesh[i + 1][j].fPos.CArr());

         glNormal3dv(fMesh[i + 1][j + 1].fNormal.CArr());
         if(fColorScheme != -1)
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fMesh[i + 1][j + 1].fRGBA);
         glVertex3dv(fMesh[i + 1][j + 1].fPos.CArr());

         glNormal3dv(fMesh[i][j].fNormal.CArr());
         if(fColorScheme != -1)
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fMesh[i][j].fRGBA);
         glVertex3dv(fMesh[i][j].fPos.CArr());
      }
   }

   glEnd();

   if (!fSelectionPass && fShowMesh) {
      glDisable(GL_POLYGON_OFFSET_FILL);
      const TGLDisableGuard lightGuard(GL_LIGHTING);
      const TGLEnableGuard blendGuard(GL_BLEND);
      const TGLEnableGuard smoothGuard(GL_LINE_SMOOTH);

      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4d(0., 0., 0., 0.5);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

      for (Int_t i = 0; i < fMeshSize - 1; ++i) {
         for (Int_t j = 0; j < fMeshSize - 1; ++j) {
            if (fBoxCut.IsActive()) {
               using TMath::Min;
               using TMath::Max;
               const Double_t xMin = Min(Min(fMesh[i][j].fPos.X(), fMesh[i + 1][j].fPos.X()), Min(fMesh[i][j + 1].fPos.X(), fMesh[i + 1][j + 1].fPos.X()));
               const Double_t xMax = Max(Max(fMesh[i][j].fPos.X(), fMesh[i + 1][j].fPos.X()), Max(fMesh[i][j + 1].fPos.X(), fMesh[i + 1][j + 1].fPos.X()));
               const Double_t yMin = Min(Min(fMesh[i][j].fPos.Y(), fMesh[i + 1][j].fPos.Y()), Min(fMesh[i][j + 1].fPos.Y(), fMesh[i + 1][j + 1].fPos.Y()));
               const Double_t yMax = Max(Max(fMesh[i][j].fPos.Y(), fMesh[i + 1][j].fPos.Y()), Max(fMesh[i][j + 1].fPos.Y(), fMesh[i + 1][j + 1].fPos.Y()));
               const Double_t zMin = Min(Min(fMesh[i][j].fPos.Z(), fMesh[i + 1][j].fPos.Z()), Min(fMesh[i][j + 1].fPos.Z(), fMesh[i + 1][j + 1].fPos.Z()));
               const Double_t zMax = Max(Max(fMesh[i][j].fPos.Z(), fMesh[i + 1][j].fPos.Z()), Max(fMesh[i][j + 1].fPos.Z(), fMesh[i + 1][j + 1].fPos.Z()));

               if (fBoxCut.IsInCut(xMin, xMax, yMin, yMax, zMin, zMax))
                  continue;
            }
            glBegin(GL_POLYGON);
            glVertex3dv(fMesh[i][j].fPos.CArr());
            glVertex3dv(fMesh[i][j + 1].fPos.CArr());
            glVertex3dv(fMesh[i + 1][j + 1].fPos.CArr());
            glVertex3dv(fMesh[i + 1][j].fPos.CArr());
            glEnd();
         }
      }

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   }

   if (fBoxCut.IsActive())
      fBoxCut.DrawBox(fSelectionPass, fSelectedPart);
}

////////////////////////////////////////////////////////////////////////////////
///Calculate colors for vertices,
///using one of 20 color themes.
///-1 simple 'metal' surface.

void TGLParametricPlot::InitColors()
{
   if (fColorScheme == -1)
      return;

   const Rgl::Range_t uRange(fEquation->GetURange());

   const Float_t dU = Float_t((uRange.second - uRange.first) / (fMeshSize - 1));
   Float_t u = Float_t(uRange.first);

   for (Int_t i = 0; i < fMeshSize; ++i) {
      for (Int_t j = 0; j < fMeshSize; ++j)
         Rgl::GetColor(u, uRange.first, uRange.second, fColorScheme, fMesh[i][j].fRGBA);
      u += dU;
   }
}

////////////////////////////////////////////////////////////////////////////////
///No such sections.

void TGLParametricPlot::DrawSectionXOZ()const
{
}

////////////////////////////////////////////////////////////////////////////////
///No such sections.

void TGLParametricPlot::DrawSectionYOZ()const
{
}

////////////////////////////////////////////////////////////////////////////////
///No such sections.

void TGLParametricPlot::DrawSectionXOY()const
{
}

////////////////////////////////////////////////////////////////////////////////
///Set material properties.

void TGLParametricPlot::SetSurfaceColor()const
{
   const Float_t specular[] = {1.f, 1.f, 1.f, 1.f};
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.f);

   if (fColorScheme == -1) {
      const Float_t outerDiff[] = {0.5f, 0.42f, 0.f, 1.f};
      glMaterialfv(GL_FRONT, GL_DIFFUSE, outerDiff);
      const Float_t innerDiff[] = {0.5f, 0.2f,  0.f, 1.f};
      glMaterialfv(GL_BACK,  GL_DIFFUSE, innerDiff);
   }
}
