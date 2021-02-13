// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  17/11/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <stdexcept>
#include <cstring>

#include "TVirtualGL.h"
#include "KeySymbols.h"
#include "Buttons.h"
#include "TH2Poly.h"
#include "TROOT.h"
#include "TGL5D.h"
#include "TVirtualPad.h"
#include "TH3.h"
#include "TF3.h"

#include "TGLSurfacePainter.h"
#include "TGLTH3Composition.h"
#include "TGLH2PolyPainter.h"
#include "TGLVoxelPainter.h"
#include "TGLHistPainter.h"
#include "TGLLegoPainter.h"
#include "TGLBoxPainter.h"
#include "TGLTF3Painter.h"
#include "TGLParametric.h"
#include "TGL5DPainter.h"
#include "TGLUtil.h"

ClassImp(TGLHistPainter);

/** \class TGLHistPainter
\ingroup opengl
The histogram painter class using OpenGL.

Histograms are, by default, drawn via the `THistPainter` class.
`TGLHistPainter` allows to paint them using the OpenGL 3D graphics
library. The plotting options provided by `TGLHistPainter` start with
`GL` keyword.

### General information: plot types and supported options

The following types of plots are provided:

#### Lego - (`TGLLegoPainter`)
   The supported options are:

  - `"GLLEGO"  :` Draw a lego plot.
  - `"GLLEGO2" :` Bins with color levels.
  - `"GLLEGO3" :` Cylindrical bars.

  Lego painter in cartesian supports logarithmic scales for X, Y, Z.
  In polar only Z axis can be logarithmic, in cylindrical only Y (if you see
  what it means).


#### Surfaces (`TF2` and `TH2` with `"GLSURF"` options) - (`TGLSurfacePainter`)
   The supported options are:

  - `"GLSURF"  :` Draw a surface.
  - `"GLSURF1" :` Surface with color levels
  - `"GLSURF2" :` The same as `"GLSURF1"` but without polygon outlines.
  - `"GLSURF3" :` Color level projection on top of plot (works only in cartesian coordinate system).
  - `"GLSURF4" :` Same as `"GLSURF"` but without polygon outlines.


  The surface painting in cartesian coordinates supports logarithmic scales along X, Y, Z axis.
  In polar coordinates only the Z axis can be logarithmic, in cylindrical coordinates only the Y axis.

#### Additional options to `SURF` and `LEGO` - Coordinate systems:
   The supported options are:

  - `" "   :` Default, cartesian coordinates system.
  - `"POL" :` Polar coordinates system.
  - `"CYL" :` Cylindrical coordinates system.
  - `"SPH" :` Spherical coordinates system.


#### `TH3` as boxes (spheres) - (`TGLBoxPainter`)
   The supported options are:

  - `"GLBOX" :` TH3 as a set of boxes, size of box is proportional to bin content.
  - `"GLBOX1":` the same as "glbox", but spheres are drawn instead of boxes.


#### `TH3` as iso-surface(s) - (`TGLIsoPainter`)
   The supported option is:

  - `"GLISO" :` TH3 is drawn using iso-surfaces.


#### `TH3` as color boxes - (`TGLVoxelPainter`)
   The supported option is:

  - `"GLCOL" :` TH3 is drawn using semi-transparent colored boxes.
  See `$ROOTSYS/tutorials/gl/glvox1.C`.


#### `TF3` (implicit function) - (`TGLTF3Painter`)
   The supported option is:

  - `"GLTF3" :` Draw a `TF3`.


#### Parametric surfaces - (`TGLParametricPlot`)
  `$ROOTSYS/tutorials/gl/glparametric.C` shows how to create parametric equations and
  visualize the surface.


### Interaction with the plots


#### General information.

  All the interactions are implemented via standard methods `DistancetoPrimitive` and
  `ExecuteEvent`. That's why all the interactions with the OpenGL plots are possible i
  only when the mouse cursor is in the plot's area (the plot's area is the part of a the pad
  occupied by gl-produced picture). If the mouse cursor is not above gl-picture,
  the standard pad interaction is performed.

#### Selectable parts.

  Different parts of the plot can be selected:

  - *xoz, yoz, xoy back planes*:
     When such a plane selected, it's highlighted in green if the dynamic slicing
     by this plane is supported, and it's highlighted in red, if the dynamic slicing
     is not supported.
  -*The plot itself*:
     On surfaces, the selected surface is outlined in red. (TF3 and ISO are not
     outlined). On lego plots, the selected bin is highlihted. The bin number and content are displayed in pad's status
     bar. In box plots, the box or sphere is highlighted and the bin info is displayed in pad's status bar.

#### Rotation and zooming.

  - *Rotation*:

  When the plot is selected, it can be rotated by pressing and holding the left mouse button and move the cursor.
  - *Zoom/Unzoom*:

  Mouse wheel or `'j'`, `'J'`, `'k'`, `'K'` keys.


#### Panning.

  The selected plot can be moved in a pad's area by
  pressing and holding the left mouse button and the shift key.

### Box cut
  Surface, iso, box, TF3 and parametric painters support box cut by pressing the `'c'` or
  `'C'` key when the mouse cursor is in a plot's area. That will display a transparent box,
  cutting away part of the surface (or boxes) in order to show internal part of plot.
  This box can be moved inside the plot's area (the full size of the box is equal to the plot's
  surrounding box) by selecting one of the box cut axes and pressing the left mouse button to move it.

### Plot specific interactions (dynamic slicing etc.)
  Currently, all gl-plots support some form of slicing.
  When back plane is selected (and if it's highlighted in green)
  you can press and hold left mouse button and shift key
  and move this back plane inside plot's area, creating the slice.
  During this "slicing" plot becomes semi-transparent. To remove all slices (and projected curves for surfaces)
  - double click with left mouse button in a plot's area.

  #### Surface with option `"GLSURF"`

  The surface profile is displayed on the slicing plane.
  The profile projection is drawn on the back plane
  by pressing `'p'` or `'P'` key.

  #### TF3

  The contour plot is drawn on the slicing plane.
  For `TF3` the color scheme can be changed by pressing `'s'` or `'S'`.

  #### Box

  The contour plot corresponding to slice plane position is drawn in real time.

  #### Iso

  Slicing is similar to `"GLBOX"` option.

  #### Parametric plot

  No slicing. Additional keys: `'s'` or `'S'` to change color scheme - about 20 color schemes supported
  (`'s'` for "scheme"); `'l'` or `'L'` to increase number of polygons (`'l'` for "level" of details),
  `'w'` or `'W'` to show outlines (`'w'` for "wireframe").
*/

////////////////////////////////////////////////////////////////////////////////
/// ROOT does not use exceptions, so, if default painter's creation failed,
/// fDefaultPainter is 0. In each function, which use it, I have to check the pointer first.

TGLHistPainter::TGLHistPainter(TH1 *hist)
                   : fDefaultPainter(TVirtualHistPainter::HistPainter(hist)),
                     fEq(0),
                     fHist(hist),
                     fF3(0),
                     fStack(0),
                     fPlotType(kGLDefaultPlot)//THistPainter
{
}

////////////////////////////////////////////////////////////////////////////////
///This ctor creates gl-parametric plot's painter.

TGLHistPainter::TGLHistPainter(TGLParametricEquation *equation)
                   : fEq(equation),
                     fHist(0),
                     fF3(0),
                     fStack(0),
                     fPlotType(kGLParametricPlot)//THistPainter
{
   fGLPainter.reset(new TGLParametricPlot(equation, &fCamera));
}

////////////////////////////////////////////////////////////////////////////////
///This ctor creates plot painter for TGL5DDataSet.

TGLHistPainter::TGLHistPainter(TGL5DDataSet *data)
                   : fEq(0),
                     fHist(0),
                     fF3(0),
                     fStack(0),
                     fPlotType(kGL5D)//THistPainter
{
   fGLPainter.reset(new TGL5DPainter(data, &fCamera, &fCoord));
}

////////////////////////////////////////////////////////////////////////////////
///This ctor creates plot painter for TGL5DDataSet.

TGLHistPainter::TGLHistPainter(TGLTH3Composition *data)
                   : fEq(0),
                     fHist(data),
                     fF3(0),
                     fStack(0),
                     fPlotType(kGLTH3Composition)
{
   fGLPainter.reset(new TGLTH3CompositionPainter(data, &fCamera, &fCoord));
}

////////////////////////////////////////////////////////////////////////////////
///Selects plot or axis.
///9999 is the magic number, ROOT's classes use in DistancetoPrimitive.

Int_t TGLHistPainter::DistancetoPrimitive(Int_t px, Int_t py)
{
   //[tp: return statement added.
   //tp]

   if (fPlotType == kGLDefaultPlot)
      return fDefaultPainter.get() ? fDefaultPainter->DistancetoPrimitive(px, py) : 9999;
   else {
      //Adjust px and py - canvas can have several pads inside, so we need to convert
      //the from canvas' system into pad's.

      //Retina-related adjustments must be done inside!!!

      py = gPad->GetWh() - py;

      //One hist can be appended to several pads,
      //the current pad should have valid OpenGL context.
      const Int_t glContext = gPad->GetGLDevice();

      if (glContext != -1) {
         //Add "viewport" extraction here.
         PadToViewport(kTRUE);

         if (!gGLManager->PlotSelected(fGLPainter.get(), px, py))
            gPad->SetSelected(gPad);
      } else {
         Error("DistancetoPrimitive",
               "Attempt to use TGLHistPainter, while the current pad (gPad) does not support gl");
         gPad->SetSelected(gPad);
      }

      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Default implementation is OK
///This function is called from a context menu
///after right click on a plot's area. Opens window
///("panel") with several controls.

void TGLHistPainter::DrawPanel()
{
   if (fDefaultPainter.get())
      fDefaultPainter->DrawPanel();
}

////////////////////////////////////////////////////////////////////////////////
///Execute event.
///Events are: mouse events in a plot's area,
///key presses (while mouse cursor is in plot's area).
///"Event execution" means one of the following actions:
///  1. Rotation.
///  2. Panning.
///  3. Zoom changing.
///  4. Moving dynamic profile.
///  5. Plot specific events - for example, 's' or 'S' key press for TF3.

void TGLHistPainter::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (fPlotType == kGLDefaultPlot) {
      if(fDefaultPainter.get()) {
         fDefaultPainter->ExecuteEvent(event, px, py);
      }
   } else {
      //One hist can be appended to several pads,
      //the current pad should have valid OpenGL context.
      const Int_t glContext = gPad->GetGLDevice();

      if (glContext == -1) {
         Error("ExecuteEvent",
               "Attempt to use TGLHistPainter, while the current pad (gPad) does not support gl");
         return;
      } else {
         //Add viewport extraction here.
         /*fGLDevice.SetGLDevice(glContext);
         fGLPainter->SetGLDevice(&fGLDevice);*/
         PadToViewport();
      }

      if (event != kKeyPress) {
         //Adjust px and py - canvas can have several pads inside, so we need to convert
         //the from canvas' system into pad's. If it was a key press event,
         //px and py ARE NOT coordinates.
         py -= Int_t((1 - gPad->GetHNDC() - gPad->GetYlowNDC()) * gPad->GetWh());
         px -= Int_t(gPad->GetXlowNDC() * gPad->GetWw());

         //We also have to take care of retina displays with a different viewports.
         TGLUtil::InitializeIfNeeded();
         const Float_t scale = TGLUtil::GetScreenScalingFactor();
         if (scale > 1) {
            px *= scale;
            py *= scale;
         }
      }

      switch (event) {
      case kButton1Double:
         //Left double click removes dynamic sections, user created (if plot type supports sections).
         fGLPainter->ProcessEvent(event, px, py);
         break;
      case kButton1Down:
         //Left mouse down in a plot area starts rotation.
         if (!fGLPainter->CutAxisSelected())
            fCamera.StartRotation(px, py);
         else
            fGLPainter->StartPan(px, py);
         //During rotation, usual TCanvas/TPad machinery (CopyPixmap/Flush/UpdateWindow/etc.)
         //is skipped - I use "bit blasting" functions to copy picture directly onto window.
         //gGLManager->MarkForDirectCopy(glContext, kTRUE);
         break;
      case kButton1Motion:
         //Rotation invalidates "selection buffer"
         // - (color-to-object map, previously read from gl-buffer).
         fGLPainter->InvalidateSelection();
         if (fGLPainter->CutAxisSelected())
            gGLManager->PanObject(fGLPainter.get(), px, py);
         else
            fCamera.RotateCamera(px, py);
         //Draw modified scene onto canvas' window.
         //gGLManager->PaintSingleObject(fGLPainter.get());
         gPad->Update();
         break;
      case kButton1Up:
      case kButton2Up:
         gGLManager->MarkForDirectCopy(glContext, kFALSE);
         break;
      case kMouseMotion:
         gPad->SetCursor(kRotate);
         break;
      case 7://kButton1Down + shift modifier
         //The current version of ROOT does not
         //have enumerators for button events + key modifiers,
         //so I use hardcoded literals. :(
         //With left mouse button down and shift pressed
         //we can move plot as the whole or move
         //plot's parts - dynamic sections.
         fGLPainter->StartPan(px, py);
         gGLManager->MarkForDirectCopy(glContext, kTRUE);
         break;
      case 8://kButton1Motion + shift modifier
         gGLManager->PanObject(fGLPainter.get(), px, py);
         //gGLManager->PaintSingleObject(fGLPainter.get());
         gPad->Update();
         break;
      case kKeyPress:
      case 5:
      case 6:
         //5, 6 are mouse wheel events (see comment about literals above).
         //'p'/'P' - specific events processed by TGLSurfacePainter,
         //'s'/'S' - specific events processed by TGLTF3Painter,
         //'c'/'C' - turn on/off box cut.
         gGLManager->MarkForDirectCopy(glContext, kTRUE);
         if (event == 6 || py == kKey_J || py == kKey_j) {
            fCamera.ZoomIn();
            fGLPainter->InvalidateSelection();
            //gGLManager->PaintSingleObject(fGLPainter.get());
            gPad->Update();
         } else if (event == 5 || py == kKey_K || py == kKey_k) {
            fCamera.ZoomOut();
            fGLPainter->InvalidateSelection();
            //gGLManager->PaintSingleObject(fGLPainter.get());
            gPad->Update();
         } else if (py == kKey_p || py == kKey_P || py == kKey_S || py == kKey_s
                    || py == kKey_c || py == kKey_C || py == kKey_x || py == kKey_X
                    || py == kKey_y || py == kKey_Y || py == kKey_z || py == kKey_Z
                    || py == kKey_w || py == kKey_W || py == kKey_l || py == kKey_L
                    /*|| py == kKey_r || py == kKey_R*/)
         {
            fGLPainter->ProcessEvent(event, px, py);
            //gGLManager->PaintSingleObject(fGLPainter.get());
            gPad->Update();
         }
         gGLManager->MarkForDirectCopy(glContext, kFALSE);
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///Get contour list.
///I do not use this function. Contours are implemented in
///a completely different way by gl-painters.

TList *TGLHistPainter::GetContourList(Double_t contour)const
{
   return fDefaultPainter.get() ? fDefaultPainter->GetContourList(contour) : 0;
}

////////////////////////////////////////////////////////////////////////////////
///Overrides TObject::GetObjectInfo.
///For lego info is: bin numbers (i, j), bin content.
///For TF2 info is: x,y,z 3d surface-point for 2d screen-point under cursor
///(this can work incorrectly now, because of wrong code in TF2).
///For TF3 no info now.
///For box info is: bin numbers (i, j, k), bin content.

char *TGLHistPainter::GetObjectInfo(Int_t px, Int_t py)const
{
   static char errMsg[] = { "TGLHistPainter::GetObjectInfo: Error in a hist painter\n" };
   if (fPlotType == kGLDefaultPlot)
      return fDefaultPainter.get() ? fDefaultPainter->GetObjectInfo(px, py)
                                   : errMsg;
   else {
      TGLUtil::InitializeIfNeeded();
      const Float_t scale = TGLUtil::GetScreenScalingFactor();
      if (scale > 1.f) {
         px *= scale;
         py *= scale;
      }

      return gGLManager->GetPlotInfo(fGLPainter.get(), px, py);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get stack.

TList *TGLHistPainter::GetStack()const
{
   return fStack;
}

////////////////////////////////////////////////////////////////////////////////
///Returns kTRUE if the cell ix, iy is inside one of the graphical cuts.
///I do not use this function anywhere, this is a "default implementation".

Bool_t TGLHistPainter::IsInside(Int_t x, Int_t y)
{
   if (fPlotType == kGLDefaultPlot)
      return fDefaultPainter.get() ? fDefaultPainter->IsInside(x, y) : kFALSE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///Returns kTRUE if the cell x, y is inside one of the graphical cuts.
///I do not use this function anywhere, this is a "default implementation".

Bool_t TGLHistPainter::IsInside(Double_t x, Double_t y)
{
   if (fPlotType == kGLDefaultPlot)
      return fDefaultPainter.get() ? fDefaultPainter->IsInside(x, y) : kFALSE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///Paint statistics.
///This does not work on windows.

void TGLHistPainter::PaintStat(Int_t dostat, TF1 *fit)
{
   if (fDefaultPainter.get())
      fDefaultPainter->PaintStat(dostat, fit);
}

////////////////////////////////////////////////////////////////////////////////
/// Process message.

void TGLHistPainter::ProcessMessage(const char *m, const TObject *o)
{
   if (!std::strcmp(m, "SetF3"))
      fF3 = (TF3 *)o;

   if (fDefaultPainter.get())
      fDefaultPainter->ProcessMessage(m, o);
}

////////////////////////////////////////////////////////////////////////////////
/// Set highlight mode

void TGLHistPainter::SetHighlight()
{
   if (fDefaultPainter.get())
      fDefaultPainter->SetHighlight();
}

////////////////////////////////////////////////////////////////////////////////
/// Set histogram.

void TGLHistPainter::SetHistogram(TH1 *h)
{
   fHist = h;

   if (fDefaultPainter.get())
      fDefaultPainter->SetHistogram(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Set stack.

void TGLHistPainter::SetStack(TList *s)
{
   fStack = s;

   if (fDefaultPainter.get())
      fDefaultPainter->SetStack(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Make cuts.

Int_t TGLHistPainter::MakeCuts(char *o)
{
   if (fPlotType == kGLDefaultPlot && fDefaultPainter.get())
      return fDefaultPainter->MakeCuts(o);

   return 0;
}

struct TGLHistPainter::PlotOption_t {
   EGLPlotType  fPlotType;
   EGLCoordType fCoordType;
   Bool_t       fBackBox;
   Bool_t       fFrontBox;
   Bool_t       fDrawAxes;
   Bool_t       fLogX;
   Bool_t       fLogY;
   Bool_t       fLogZ;
};

////////////////////////////////////////////////////////////////////////////////
///Final-overrider for TObject::Paint.

void TGLHistPainter::Paint(Option_t *o)
{
   TString option(o);
   option.ToLower();

   const Ssiz_t glPos = option.Index("gl");
   if (glPos != kNPOS)
      option.Remove(glPos, 2);
   else if (fPlotType != kGLParametricPlot && fPlotType != kGL5D && fPlotType != kGLTH3Composition) {
      gPad->SetCopyGLDevice(kFALSE);
      if (fDefaultPainter.get())
         fDefaultPainter->Paint(o);//option.Data());
      return;
   }

   if (fPlotType != kGLParametricPlot && fPlotType != kGL5D && fPlotType != kGLTH3Composition)
      CreatePainter(ParsePaintOption(option), option);

   if (fPlotType == kGLDefaultPlot) {
      //In case of default plot pad
      //should not copy gl-buffer (it will be simply black)

      //[tp: code was commented.
      //gPad->SetCopyGLDevice(kFALSE);
      //tp]

      if (fDefaultPainter.get())
         fDefaultPainter->Paint(option.Data());
   } else {
      Int_t glContext = gPad->GetGLDevice();

      if (glContext != -1) {
         //With gl-plot, pad should copy
         //gl-buffer into the final pad/canvas pixmap/DIB.
         //fGLDevice.SetGLDevice(glContext);

         //[tp: code commented.
         //gPad->SetCopyGLDevice(kTRUE);
         //tp]
         //fGLPainter->SetGLDevice(&fGLDevice);
         //Add viewport extraction here.
         PadToViewport();
         if (gPad->GetFrameFillColor() != kWhite)
            fGLPainter->SetFrameColor(gROOT->GetColor(gPad->GetFrameFillColor()));
         fGLPainter->SetPadColor(gROOT->GetColor(gPad->GetFillColor()));
         if (fGLPainter->InitGeometry())
            gGLManager->PaintSingleObject(fGLPainter.get());
      }
   }
}

namespace {

Bool_t FindAndRemoveOption(TString &options, const char *toFind)
{
   const UInt_t len = std::strlen(toFind);
   const Ssiz_t index = options.Index(toFind);

   if (index != kNPOS) {
      options.Remove(index, len);
      return kTRUE;
   }

   return kFALSE;
}

}

////////////////////////////////////////////////////////////////////////////////
///In principle, we can have several conflicting options: "lego surf pol sph", surfbb: surf, fb, bb.
///but only one will be selected, which one - depends on parsing order in this function.

TGLHistPainter::PlotOption_t
TGLHistPainter::ParsePaintOption(const TString &o)const
{
   TString options(o);

   PlotOption_t parsedOption = {kGLDefaultPlot, kGLCartesian,
                                kTRUE, kTRUE, kTRUE, //Show back box, show front box, show axes.
                                Bool_t(gPad->GetLogx()), Bool_t(gPad->GetLogy()),
                                Bool_t(gPad->GetLogz())};

   //Check coordinate system type.
   if (FindAndRemoveOption(options, "pol"))
      parsedOption.fCoordType = kGLPolar;
   if (FindAndRemoveOption(options, "cyl"))
      parsedOption.fCoordType = kGLCylindrical;
   if (FindAndRemoveOption(options, "sph"))
      parsedOption.fCoordType = kGLSpherical;

   //Define plot type.
   if (FindAndRemoveOption(options, "lego"))
      fStack ? parsedOption.fPlotType = kGLStackPlot : parsedOption.fPlotType = kGLLegoPlot;
   if (FindAndRemoveOption(options, "surf"))
      parsedOption.fPlotType = kGLSurfacePlot;
   if (FindAndRemoveOption(options, "tf3"))
      parsedOption.fPlotType = kGLTF3Plot;
   if (FindAndRemoveOption(options, "box"))
      parsedOption.fPlotType = kGLBoxPlot;
   if (FindAndRemoveOption(options, "iso"))
      parsedOption.fPlotType = kGLIsoPlot;
   if (FindAndRemoveOption(options, "col"))
      parsedOption.fPlotType = kGLVoxel;

   //Check BB and FB options.
   if (FindAndRemoveOption(options, "bb"))
      parsedOption.fBackBox = kFALSE;
   if (FindAndRemoveOption(options, "fb"))
      parsedOption.fFrontBox = kFALSE;

   //Check A option.
   if (FindAndRemoveOption(options, "a"))
      parsedOption.fDrawAxes = kFALSE;

   return parsedOption;
}

////////////////////////////////////////////////////////////////////////////////
/// Create painter.

void TGLHistPainter::CreatePainter(const PlotOption_t &option, const TString &addOption)
{
   if (option.fPlotType != fPlotType) {
      fCoord.ResetModified();
      fGLPainter.reset(0);
   }

   if (option.fPlotType == kGLLegoPlot) {
      if (!fGLPainter.get()) {
         if (dynamic_cast<TH2Poly*>(fHist))
            fGLPainter.reset(new TGLH2PolyPainter(fHist, &fCamera, &fCoord));
         else
            fGLPainter.reset(new TGLLegoPainter(fHist, &fCamera, &fCoord));
      }
   } else if (option.fPlotType == kGLSurfacePlot) {
      if (!fGLPainter.get())
         fGLPainter.reset(new TGLSurfacePainter(fHist, &fCamera, &fCoord));
   } else if (option.fPlotType == kGLBoxPlot) {
      if (!fGLPainter.get())
         fGLPainter.reset(new TGLBoxPainter(fHist, &fCamera, &fCoord));
   } else if (option.fPlotType == kGLTF3Plot) {
      if (!fGLPainter.get())
         fGLPainter.reset(new TGLTF3Painter(fF3, fHist, &fCamera, &fCoord));
   } else if (option.fPlotType == kGLIsoPlot) {
      if (!fGLPainter.get())
         fGLPainter.reset(new TGLIsoPainter(fHist, &fCamera, &fCoord));
   } else if (option.fPlotType == kGLVoxel) {
      if (!fGLPainter.get())
         fGLPainter.reset(new TGLVoxelPainter(fHist, &fCamera, &fCoord));
   }

   if (fGLPainter.get()) {
      fPlotType = option.fPlotType;
      fCoord.SetXLog(gPad->GetLogx());
      fCoord.SetYLog(gPad->GetLogy());
      fCoord.SetZLog(gPad->GetLogz());
      fCoord.SetCoordType(option.fCoordType);
      fGLPainter->AddOption(addOption);

      fGLPainter->SetDrawFrontBox(option.fFrontBox);
      fGLPainter->SetDrawBackBox(option.fBackBox);
      fGLPainter->SetDrawAxes(option.fDrawAxes);
   } else
      fPlotType = kGLDefaultPlot;
}

////////////////////////////////////////////////////////////////////////////////
/// Set show projection.

void TGLHistPainter::SetShowProjection(const char *option, Int_t nbins)
{
   if (fDefaultPainter.get()) fDefaultPainter->SetShowProjection(option, nbins);
}

////////////////////////////////////////////////////////////////////////////////

void TGLHistPainter::PadToViewport(Bool_t /*selectionPass*/)
{
   if (!fGLPainter.get())
      return;

   TGLRect vp;
   vp.Width()  = Int_t(gPad->GetAbsWNDC() * gPad->GetWw());
   vp.Height() = Int_t(gPad->GetAbsHNDC() * gPad->GetWh());

   vp.X() = Int_t(gPad->XtoAbsPixel(gPad->GetX1()));
   vp.Y() = Int_t((gPad->GetWh() - gPad->YtoAbsPixel(gPad->GetY1())));

   TGLUtil::InitializeIfNeeded();
   const Float_t scale = TGLUtil::GetScreenScalingFactor();

   if (scale > 1.f) {
      vp.X() = Int_t(vp.X() * scale);
      vp.Y() = Int_t(vp.Y() * scale);

      vp.Width() = Int_t(vp.Width() * scale);
      vp.Height() = Int_t(vp.Height() * scale);
   }

   fCamera.SetViewport(vp);
   if (fCamera.ViewportChanged() && fGLPainter.get())
      fGLPainter->InvalidateSelection();
}
