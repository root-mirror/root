/*
 * $Header$
 * $Log$
 *
 * Various routines, and global variables
 */

#define __XSVARIOUS_CXX

#include <TEnv.h>
#include <TGClient.h>
#include <TGWindow.h>
#include <TGWindow.h>
#include <TVirtualX.h>

#include "XSGraph.h"
#include "XSVarious.h"

/* --------------- XSinitialise ---------------- */
void
XSinitialise( )
{
   XSelements = new XSElements(ISOTOPES_DESC_FILE);
   XSReactionDesc = new NdbMTReacDesc(MT_DESC_FILE);
   graphList = new TList();

   // --- Initialise GUI variables ---
   fixedFontStruct = gClient->GetFontByName(
      gEnv->GetValue("Gui.NormalFont","fixed"));

   // Define new graphics context.
   memset(&gval,0,sizeof(gval));
   gval.fMask =    kGCForeground | kGCBackground | kGCFont |
   kGCFillStyle | kGCGraphicsExposures;
   gval.fFillStyle = kFillSolid;
   gval.fGraphicsExposures = kFALSE;
   gval.fFont = gVirtualX->GetFontHandle(fixedFontStruct);
   fixedGC=gVirtualX->CreateGC(gClient->GetRoot()->GetId(), &gval);

   blueFontStruct = gClient->GetFontByName(
         gEnv->GetValue("Gui.NormalFont",
         "-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1"));
   // Define new graphics context.
   gval.fMask = kGCForeground | kGCFont;
   gval.fFont = gVirtualX->GetFontHandle(blueFontStruct);
   gClient->GetColorByName("blue", gval.fForeground);
   blueBoldGC=gVirtualX->CreateGC(gClient->GetRoot()->GetId(), &gval);

} // XSinitialise

/* ---------------- XSfinalise ----------------- */
void
XSfinalise()
{
   delete   XSelements;
   delete   XSReactionDesc;
   delete   graphList;
} // XSfinalise

/* --------------- Add2GraphList --------------- */
void
Add2GraphList( XSGraph *gr)
{
   graphList->Add(gr);

   /* --- Scan all graphs and update the canvas --- */
   // ....
} // Add2GraphList
