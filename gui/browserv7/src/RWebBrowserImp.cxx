// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2021-02-11
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RWebBrowserImp.hxx>

using namespace ROOT::Experimental;

////////////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RWebBrowserImp::RWebBrowserImp(TBrowser *b) : TBrowserImp(b)
{
   fWebBrowser = std::make_shared<RBrowser>();
}

////////////////////////////////////////////////////////////////////////////////////////
/// Constructor with width and height parameters

RWebBrowserImp::RWebBrowserImp(TBrowser *b, const char *title, UInt_t width, UInt_t height, Option_t *opt) : TBrowserImp(b,title, width, height, opt)
{
   fWidth = width;
   fHeight = height;
   fWebBrowser = std::make_shared<RBrowser>();
}

////////////////////////////////////////////////////////////////////////////////////////
/// Constructor with x,y, width and height parameters

RWebBrowserImp::RWebBrowserImp(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height, Option_t *opt) : TBrowserImp(b,title, x, y, width, height, opt)
{
   fX = x;
   fY = y;
   fWidth = width;
   fHeight = height;
   fWebBrowser = std::make_shared<RBrowser>();
}

////////////////////////////////////////////////////////////////////////////////////////
/// Constructor with width and height parameters

RWebBrowserImp::~RWebBrowserImp()
{
}

////////////////////////////////////////////////////////////////////////////////////////
/// Iconify browser

void RWebBrowserImp::Iconify()
{
}

////////////////////////////////////////////////////////////////////////////////////////
/// Refresh browser

void RWebBrowserImp::Refresh(Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////////////
/// Show browser

void RWebBrowserImp::Show()
{
   fWebBrowser->Show();
}

////////////////////////////////////////////////////////////////////////////////////////
/// Factory method to create RWebBrowserImp via plugin

TBrowserImp *RWebBrowserImp::NewBrowser(TBrowser *b, const char *title, UInt_t width, UInt_t height, Option_t *opt)
{
   return new RWebBrowserImp(b, title, width, height, opt);
}

////////////////////////////////////////////////////////////////////////////////////////
/// Factory method to create RWebBrowserImp via plugin

TBrowserImp *RWebBrowserImp::NewBrowser(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height, Option_t *opt)
{
   return new RWebBrowserImp(b, title, x, y, width, height, opt);
}
