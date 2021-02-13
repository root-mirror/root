// @(#)root/gl:$Id$
// Author:  Richard Maunder, Olivier Couet  02/07/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLOutput.h"

#include <iostream>
#include <fstream>
#include "TVirtualPad.h"
#include "TVirtualPS.h"
#include "TGLViewer.h"
#include "TSystem.h"
#include "gl2ps.h"
#include "TError.h"
#include <cassert>

/** \class TGLOutput
\ingroup opengl
Wrapper class for GL capture & output routines
*/

ClassImp(TGLOutput);

////////////////////////////////////////////////////////////////////////////////
/// Capture viewer to file. Arguments are:
///  - 'viewer' - viewer object to capture from
///  - 'format' - output format - only postscript types presently.
///               One of kEPS_SIMPLE, kEPS_BSP, kPDF_SIMPLE or kPDF_BSP
///               See TGLOutput::CapturePostscript() for meanings
///  - 'filePath' - file output name. If null defaults to './viewer.eps' or './viewer.pdf'
///                depending on format requested
///
/// Note : Output files can be large and take considerable time (up to mins)
/// to generate

Bool_t TGLOutput::Capture(TGLViewer & viewer, EFormat format, const char * filePath)
{
   switch(format) {
      case(kEPS_SIMPLE):
      case(kEPS_BSP):
      case(kPDF_SIMPLE):
      case(kPDF_BSP): {
         return CapturePostscript(viewer, format, filePath);
      }
   }

   assert(kFALSE);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Capture viewer to postscript file. Arguments are:
///  - 'viewer' - viewer object to capture from
///  - 'format' - output format
///              -  kEPS_SIMPLE - lower quality EPS
///              -  kEPS_BSP    - higher quality EPS
///              -  kPDF_SIMPLE - lower quality PDF
///              -  kPDF_BSP    - higher quality PDF
///  - 'filePath' - file output name. If null defaults to './viewer.eps' or './viewer.pdf'
/// depending on format requested

Bool_t TGLOutput::CapturePostscript(TGLViewer & viewer, EFormat format, const char * filePath)
{
   if (!filePath || !filePath[0]) {
      if (format == kEPS_SIMPLE || format == kEPS_BSP) {
         filePath = "viewer.eps";
      } else if (format == kPDF_SIMPLE || format == kPDF_BSP) {
         filePath = "viewer.pdf";
      }
   }
   Info("TGLOutput::Postscript()", "Start creating %s.", filePath);
   std::cout << "Please wait.";

   if (FILE *output = fopen (filePath, "w+b"))
   {
      Int_t gl2psFormat;
      Int_t gl2psSort;

      switch(format) {
         case kEPS_SIMPLE:
            gl2psFormat = GL2PS_EPS;
            gl2psSort = GL2PS_SIMPLE_SORT;
            break;
         case kEPS_BSP:
            gl2psFormat = GL2PS_EPS;
            gl2psSort = GL2PS_BSP_SORT;
            break;
         case kPDF_SIMPLE:
            gl2psFormat = GL2PS_PDF;
            gl2psSort = GL2PS_SIMPLE_SORT;
            break;
         case kPDF_BSP:
            gl2psFormat = GL2PS_PDF;
            gl2psSort = GL2PS_BSP_SORT;
            break;
         default:
            assert(kFALSE);
            return kFALSE;
      }
      Int_t buffsize = 0, state = GL2PS_OVERFLOW;
      viewer.DoDraw();
      viewer.fIsPrinting = kTRUE;
      while (state == GL2PS_OVERFLOW) {
         buffsize += 1024*1024;
         gl2psBeginPage ("ROOT Scene Graph", "ROOT", NULL,
         gl2psFormat, gl2psSort, GL2PS_USE_CURRENT_VIEWPORT
         | GL2PS_SILENT | GL2PS_BEST_ROOT | GL2PS_OCCLUSION_CULL | 0,
         GL_RGBA, 0, NULL,0, 0, 0,
         buffsize, output, NULL);
         viewer.DoDraw();
         state = gl2psEndPage();
         std::cout << ".";
      }
      std::cout << std::endl;
      fclose (output);
      viewer.fIsPrinting = kFALSE;
      if (!gSystem->AccessPathName(filePath)) {
         Info("TGLOutput::Postscript", "Finished creating %s.", filePath);
         return kTRUE;
      }
   } else {
      Error("TGLOutput::Postscript", "Failed to create %s. ", filePath);
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// this function used by gl-in-pad

void TGLOutput::StartEmbeddedPS()
{
   Info("TGLOutput::StartEmbeddedPS", "PS output started ...");

   gVirtualPS->PrintStr("@");
   gVirtualPS->PrintStr("% Start gl2ps EPS@");
   gVirtualPS->PrintStr("newpath gsave save@");
   Double_t xx[2] = {0.}, yy[2] = {0.};
   xx[0] = gPad->GetUxmin();
   yy[0] = gPad->GetUymin();
   xx[1] = gPad->GetUxmax();
   yy[1] = gPad->GetUymax();
   gVirtualPS->PrintStr("@");

   GLint vp[4];
   glGetIntegerv(GL_VIEWPORT,vp);
   gVirtualPS->DrawPS(0, xx, yy);
   gVirtualPS->PrintStr(" exch");
   xx[0] = xx[1];
   yy[0] = yy[1];
   gVirtualPS->DrawPS(0, xx, yy);
   gVirtualPS->PrintStr(" 4 1 roll exch sub 3 1 roll sub");
   gVirtualPS->WriteInteger(2*4*gPad->GetBorderSize());
   gVirtualPS->PrintStr(" sub exch");
   gVirtualPS->WriteInteger(2*4*gPad->GetBorderSize());
   gVirtualPS->PrintStr(" sub exch");
   gVirtualPS->WriteInteger((Int_t)(vp[3]));
   gVirtualPS->WriteInteger((Int_t)(vp[2]));
   gVirtualPS->PrintStr(" 4 1 roll div 3 1 roll exch div exch scale@");
   gVirtualPS->PrintStr("@");
   gVirtualPS->PrintStr("countdictstack@");
   gVirtualPS->PrintStr("mark@");
   gVirtualPS->PrintStr("/showpage {} def@");

   // Close the gVirtualPS output stream
   std::ofstream *fs = (std::ofstream*)gVirtualPS->GetStream();
   fs->close();

}

////////////////////////////////////////////////////////////////////////////////
///this function used by gl-in-pad
/// Restore the gVirtualPS output stream

void TGLOutput::CloseEmbeddedPS()
{
   std::ofstream *fs = new std::ofstream(gVirtualPS->GetName(),std::ios::app);
   gVirtualPS->SetStream(fs);
   gVirtualPS->PrintStr("@");
   gVirtualPS->PrintStr("cleartomark@");
   gVirtualPS->PrintStr("countdictstack exch sub { end } repeat@");
   gVirtualPS->PrintStr("restore grestore@");
   gVirtualPS->PrintStr("% End gl2ps EPS@");

   Info("TGLOutput::CloseEmbeddedPS", "PS output finished");
}

////////////////////////////////////////////////////////////////////////////////
///this function used by gl-viewer, embedded into pad

void TGLOutput::Capture(TGLViewer & viewer)
{
   StartEmbeddedPS();

   FILE *output = fopen (gVirtualPS->GetName(), "a");
   if (!output) {
      //Quite stupid fix, since if fopen fails, CloseEmbeddedPS will also
      //fail but still I have to do it.
      Error("TGLOutput::Capture", "can not open file for embedding ps");
      CloseEmbeddedPS();
      return;
   }

   Int_t gl2psFormat = GL2PS_EPS;
   Int_t gl2psSort = GL2PS_BSP_SORT;
   Int_t buffsize = 0, state = GL2PS_OVERFLOW;
   viewer.DoDraw();
   viewer.fIsPrinting = kTRUE;

   while (state == GL2PS_OVERFLOW) {
      buffsize += 1024*1024;
      gl2psBeginPage ("ROOT Scene Graph", "ROOT", NULL,
      gl2psFormat, gl2psSort, GL2PS_USE_CURRENT_VIEWPORT
      | GL2PS_SILENT | GL2PS_BEST_ROOT | GL2PS_OCCLUSION_CULL | 0,
      GL_RGBA, 0, NULL,0, 0, 0,
      buffsize, output, NULL);
      viewer.DoDraw();
      state = gl2psEndPage();
      std::cout << ".";
   }

   std::cout << std::endl;
   fclose (output);
   viewer.fIsPrinting = kFALSE;

   CloseEmbeddedPS();
}
