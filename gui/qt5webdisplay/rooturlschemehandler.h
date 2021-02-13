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

#ifndef ROOT_UrlSchemeHandler
#define ROOT_UrlSchemeHandler

#include <QWebEngineUrlSchemeHandler>

class THttpServer;

class UrlRequestJobHolder : public QObject {
   Q_OBJECT

   QWebEngineUrlRequestJob *fRequest{nullptr};

public:
   UrlRequestJobHolder(QWebEngineUrlRequestJob *req);

   QWebEngineUrlRequestJob *req() const { return fRequest; }

   void reset();

public slots:

   void onRequestDeleted(QObject *obj);

};

// ===============================================================


class RootUrlSchemeHandler : public QWebEngineUrlSchemeHandler {
   Q_OBJECT
protected:

   THttpServer *fServer{nullptr}; ///< server instance which should handle requests

public:
   QString MakeFullUrl(THttpServer *serv, const QString &url);

   void requestStarted(QWebEngineUrlRequestJob *request) override;
};


#endif
