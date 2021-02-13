// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2017-06-29
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "rootwebpage.h"

#include <ROOT/RLogger.hxx>
#include "TString.h"
#include "TEnv.h"

ROOT::Experimental::RLogChannel &QtWebDisplayLog()
{
   static ROOT::Experimental::RLogChannel sChannel("ROOT.QtWebDisplay");
   return sChannel;
}

/** \class RootWebPage
\ingroup qt5webdisplay
*/


RootWebPage::RootWebPage(QObject *parent) : QWebEnginePage(parent)
{
   fConsole = gEnv->GetValue("WebGui.Console", (int)0);
}

void RootWebPage::javaScriptConsoleMessage(JavaScriptConsoleMessageLevel lvl, const QString &message, int lineNumber,
                                           const QString &src)
{
   switch (lvl) {
   case InfoMessageLevel:
      if (fConsole > 0)
         R__LOG_DEBUG(0, QtWebDisplayLog()) << Form("%s:%d: %s", src.toLatin1().constData(), lineNumber,
                                     message.toLatin1().constData());
      break;
   case WarningMessageLevel:
      if (fConsole > -1)
         R__LOG_WARNING(QtWebDisplayLog()) << Form("%s:%d: %s", src.toLatin1().constData(), lineNumber,
                                       message.toLatin1().constData());
      break;
   case ErrorMessageLevel:
      if (fConsole > -2)
         R__LOG_ERROR(QtWebDisplayLog()) << Form("%s:%d: %s", src.toLatin1().constData(), lineNumber,
                                     message.toLatin1().constData());
      break;
   }
}
