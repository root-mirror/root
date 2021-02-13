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

#ifndef ROOT_RootWebPage
#define ROOT_RootWebPage

#include <QWebEnginePage>

namespace ROOT {
namespace Experimental {
class RLogChannel;
}
}

ROOT::Experimental::RLogChannel &QtWebDisplayLog();

class RootWebPage : public QWebEnginePage {
   Q_OBJECT
protected:
   int fConsole{0};
   virtual void javaScriptConsoleMessage(QWebEnginePage::JavaScriptConsoleMessageLevel level, const QString &message,
                                         int lineNumber, const QString &sourceID);

public:
   RootWebPage(QObject *parent = nullptr);
   virtual ~RootWebPage() = default;
};

#endif
