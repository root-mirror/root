/// \file rootwebpage.cpp
/// \ingroup WebGui
/// \author Sergey Linev <S.Linev@gsi.de>
/// \date 2017-06-29
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

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
         R__DEBUG_HERE("Qt") << Form("%s:%d: %s", src.toLatin1().constData(), lineNumber,
                                     message.toLatin1().constData());
      break;
   case WarningMessageLevel:
      if (fConsole > -1)
         R__WARNING_HERE("Qt") << Form("%s:%d: %s", src.toLatin1().constData(), lineNumber,
                                       message.toLatin1().constData());
      break;
   case ErrorMessageLevel:
      if (fConsole > -2)
         R__ERROR_HERE("Qt") << Form("%s:%d: %s", src.toLatin1().constData(), lineNumber,
                                     message.toLatin1().constData());
      break;
   }
}
