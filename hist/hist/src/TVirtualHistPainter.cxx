// @(#)root/hist:$Id$
// Author: Rene Brun   30/08/99
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TVirtualHistPainter.h"
#include "TPluginManager.h"

TClass *TVirtualHistPainter::fgPainter = nullptr;


ClassImp(TVirtualHistPainter);

/** \class TVirtualHistPainter
 \ingroup Hist
 Abstract interface to a histogram painter
*/


////////////////////////////////////////////////////////////////////////////////
/// Static function returning a pointer to the current histogram painter.
/// The painter will paint the specified obj. If the histogram painter
/// does not exist a default painter is created.

TVirtualHistPainter *TVirtualHistPainter::HistPainter(TH1 *obj)
{
   // if no painter set yet, create a default painter via the PluginManager
   if (!fgPainter) {
      TPluginHandler *h;
      if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualHistPainter"))) {
         if (h->LoadPlugin() == -1)
            return 0;
         TVirtualHistPainter::SetPainter(h->GetClass());
         if (!fgPainter) return 0;
      } else {
         // fgPainter is still null
         return 0;
      }
   }

   //create an instance of the histogram painter
   TVirtualHistPainter *p = (TVirtualHistPainter*)fgPainter->New();
   if (p) p->SetHistogram(obj);
   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to set an alternative histogram painter.

void TVirtualHistPainter::SetPainter(const char *painter)
{
   fgPainter = TClass::GetClass(painter);
}
