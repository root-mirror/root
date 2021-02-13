// @(#)root/gpad:$Id$
// Author: Nenad Buncic   20/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TControlBar
\ingroup gpad

A Control Bar is a fully user configurable tool which provides fast
access to frequently used operations. The user can choose between
buttons and drawn buttons (let's say icons) and assign to them their
own actions (let's say ROOT or C++ commands).

The macro below shows an example of control bar.
To execute an item, click with the left mouse button.
To see the HELP of a button, click on the right mouse button.

You have access to the last clicked button via the method
GetClicked(). For example, bar->GetClicked()->GetName()
will return the name of the last clicked button.

~~~ {.cpp}
{
   gROOT->Reset("a");
   TControlBar bar("vertical");
   bar.AddButton("Help to run demos",".x demoshelp.C",
                 "Explains how to run the demos");
   bar.AddButton("framework",        ".x framework.C",
                 "An Example of Object Oriented User Interface");
   bar.AddButton("hsimple",          ".x hsimple.C",
                 "An Example Creating Histograms/Ntuples on File");
   bar.AddButton("hsum",             ".x hsum.C",
                 "Filling histograms and some graphics options");
   bar.AddButton("canvas",           ".x canvas.C",
                 "Canvas and Pad Management");
   bar.AddButton("formula1",         ".x formula1.C",
                 "Simple Formula and Functions");
   bar.AddButton("fillrandom",       ".x fillrandom.C",
                 "Histograms with Random Numbers from a Function");
   bar.AddButton("fit1",             ".x fit1.C",
                 "A Simple Fitting Example");
   bar.AddButton("h1draw",           ".x h1draw.C",
                 "Drawing Options for 1D Histograms");
   bar.AddButton("graph",            ".x graph.C",
                 "Examples of a simple graph");
   bar.AddButton("tornado",          ".x tornado.C",
                 "Examples of 3-D PolyMarkers");
   bar.AddButton("shapes",           ".x shapes.C",
                 "The Geometry Shapes");
   bar.AddButton("atlasna49",        ".x atlasna49.C",
                 "Creating and Viewing Geometries");
   bar.AddButton("file_layout",      ".x file.C",
                 "The ROOT file format");
   bar.AddButton("tree_layout",      ".x tree.C",
                 "The Tree Data Structure");
   bar.AddButton("ntuple1",          ".x ntuple1.C",
                 "Ntuples and Selections");
   bar.AddButton("run benchmarks",   ".x benchmarks.C",
                 "Runs all the ROOT benchmarks");
   bar.AddButton("rootmarks",        ".x rootmarks.C",
                 "Prints an estimated ROOTMARKS for your machine");
   bar.AddButton("edit_hsimple",     ".!ved hsimple.C &",
                 "Invokes the text editor on file hsimple.C");
   bar.AddButton("Close Bar",        "gROOT.Reset(\"a\")",
                 "Close ControlBar");
   bar.Show();
   gROOT->SaveContext();
}
~~~
\image html gpad_controlbar.png
*/

#include "TApplication.h"
#include "TControlBar.h"
#include "TControlBarImp.h"
#include "TGuiFactory.h"
#include "TList.h"
#include "TStyle.h"


ClassImp(TControlBar);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TControlBar::TControlBar() : TControlBarButton()
{
   fControlBarImp = 0;
   fOrientation   = 0;
   fButtons       = 0;
   fNoroc         = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TControlBar::TControlBar(const char *orientation, const char *title)
            : TControlBarButton(title, "", "", "button")
{
   SetOrientation( orientation );
   Initialize(-999, -999);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TControlBar::TControlBar(const char *orientation, const char *title, Int_t x, Int_t y)
            : TControlBarButton(title, "", "", "button")
{
   Int_t xs = (Int_t)(x*gStyle->GetScreenFactor());
   Int_t ys = (Int_t)(y*gStyle->GetScreenFactor());
   SetOrientation( orientation );
   Initialize(xs, ys);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TControlBar::~TControlBar()
{
   delete fControlBarImp;

   if( fButtons )
      fButtons->Delete();

   fButtons       = 0;
   fControlBarImp = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add button.

void TControlBar::AddButton(TControlBarButton &button)
{
   AddButton( &button );
}

////////////////////////////////////////////////////////////////////////////////
/// Add button.

void TControlBar::AddButton(TControlBarButton *button)
{
   if( fButtons && button )
      fButtons->Add( button );
}

////////////////////////////////////////////////////////////////////////////////
/// Add button.

void TControlBar::AddButton(const char *label, const char *action, const char *hint, const char *type)
{
   TControlBarButton *button = new TControlBarButton( label, action, hint, type );
   AddButton( button );
}

////////////////////////////////////////////////////////////////////////////////
/// Add control bar.

void TControlBar::AddControlBar(TControlBar &controlBar)
{
   AddControlBar( &controlBar );
}

////////////////////////////////////////////////////////////////////////////////
/// Add control bar.

void TControlBar::AddControlBar(TControlBar *controlBar)
{
   if( fButtons && controlBar )
      fButtons->Add( controlBar );
}

////////////////////////////////////////////////////////////////////////////////
/// Add separator.

void TControlBar::AddSeparator()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create control bar.

void TControlBar::Create()
{
   if( fControlBarImp ) {
      fControlBarImp->Create();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Hide control bar.

void TControlBar::Hide()
{
   if( fControlBarImp ) {
      fControlBarImp->Hide();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize control bar.

void TControlBar::Initialize(Int_t x, Int_t y)
{
   // Load and initialize graphics libraries if
   // TApplication::NeedGraphicsLibs() has been called by a
   // library static initializer.
   if (gApplication)
      gApplication->InitializeGraphics();

   if (x == -999) {
      fControlBarImp = gGuiFactory->CreateControlBarImp( this, GetName() );
   } else {
      fControlBarImp = gGuiFactory->CreateControlBarImp( this, GetName(), x, y );
   }

   fButtons       = new TList();
   fNoroc         = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets new font for control bar buttons, e.g.:
/// ~~~ {.cpp}
/// root > .x tutorials/demos.C
/// root > bar->SetFont("-adobe-helvetica-bold-r-*-*-24-*-*-*-*-*-iso8859-1")
/// ~~~

void TControlBar::SetFont(const char *fontName)
{
   fControlBarImp->SetFont(fontName);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets text color for control bar buttons, e.g.:
/// ~~~ {.cpp}
/// root > .x tutorials/demos.C
/// root > bar->SetTextColor("red")
/// ~~~

void TControlBar::SetTextColor(const char *colorName)
{
   fControlBarImp->SetTextColor(colorName);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets a state for control bar button 'label'; possible states are
/// 0-kButtonUp, 1-kButtonDown, 2-kButtonEngaged, 3-kButtonDisabled,
///
/// e.g.:
/// ~~~ {.cpp}
/// root > .x tutorials/demos.C
/// ~~~
/// to disable the button 'first' do:
/// ~~~ {.cpp}
/// root > bar->SetButtonState("first", 3)
/// ~~~
/// to enable the button 'first' do:
/// ~~~ {.cpp}
/// root > bar->SetButtonState("first", 0)
/// ~~~

void TControlBar::SetButtonState(const char *label, Int_t state)
{
   if (state > 3) {
      Error("SetButtonState", "not valid button state (expecting 0, 1, 2 or 3)");
      return;
   }
   fControlBarImp->SetButtonState(label, state);
}


 ///////////////////////////////////////////////////////////////////////////////
 /// Sets the width in pixels for control bar button.

void TControlBar::SetButtonWidth(UInt_t width)
{
   fControlBarImp->SetButtonWidth(width);
}

////////////////////////////////////////////////////////////////////////////////
/// Set control bar orientation.

void TControlBar::SetOrientation(const char *o)
{
   fOrientation = kVertical;

   if( *o ) {
      if( !strcasecmp( o, "vertical" ) )
         fOrientation = kVertical;
      else if( !strcasecmp( o, "horizontal" ) )
         fOrientation = kHorizontal;
      else
         Error( "SetOrientation", "Unknown orientation: '%s' !\n\t\t(choice of: %s, %s)",
                 o, "vertical", "horizontal" );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set control bar orientation.

void TControlBar::SetOrientation(Int_t o)
{
   fOrientation = kVertical;

   if( ( o == kVertical ) || ( o == kHorizontal ) )
      fOrientation = o;
   else
      Error( "SetOrientation", "Unknown orientation: %d !\n\t\t(choice of: %d, %d)",
              o, kVertical, kHorizontal );
}

////////////////////////////////////////////////////////////////////////////////
/// Show control bar.

void TControlBar::Show()
{
   if( fControlBarImp )
      fControlBarImp->Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the last clicked control bar button;
/// null if no button was clicked yet

TControlBarButton *TControlBar::GetClicked() const
{
   if (!fControlBarImp->GetClicked())
      Printf("None of the control bar buttons is clicked yet");
   return fControlBarImp->GetClicked();
}
