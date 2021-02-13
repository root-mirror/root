// @(#)root/base:$Id$
// Author: Nenad Buncic   08/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TContextMenu
#define ROOT_TContextMenu


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TContextMenu                                                               //
//                                                                            //
// This class provides an interface to  context sensitive popup menus.        //
// These menus pop up when the user hits  the right mouse button,  and        //
// are destroyed when the menu pops downs.                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#ifdef R__LESS_INCLUDES
class TContextMenuImp;
#else
#include "TContextMenuImp.h"
#endif

class TMethod;
class TFunction;
class TMethodArg;
class TVirtualPad;
class TObjArray;
class TBrowser;
class TToggle;
class TClassMenuItem;


class TContextMenu : public TNamed {

friend class  TContextMenuImp;

private:
   TContextMenu(const TContextMenu&) = delete;            // TContextMenu can not be copied since we do not know the actual type of the TContextMenuImp (and it can not be 'Cloned')
   TContextMenu& operator=(const TContextMenu&) = delete; // TContextMenu can not be copied since we do not know the actual type of the TContextMenuImp (and it can not be 'Cloned')

protected:
   TContextMenuImp *fContextMenuImp;      //!Context menu system specific implementation
   TFunction       *fSelectedMethod;      //selected method
   TObject         *fSelectedObject;      //selected object
   TObject         *fCalledObject;        //object to call
   TClassMenuItem  *fSelectedMenuItem;    //selected class menu item
   TVirtualPad     *fSelectedCanvas;      //selected canvas (if exist)
   TVirtualPad     *fSelectedPad;         //selected pad (if exist)
   TBrowser        *fBrowser;             //selected browser (if exist)

   virtual void DisplayPopUp(Int_t x, Int_t y);

private:
   TContextMenu();

public:
   TContextMenu(const char *name, const char *title = "Context sensitive popup menu");
   virtual ~TContextMenu();

   virtual void Action(TObject *object, TMethod *method);
   virtual void Action(TObject *object, TToggle *toggle);
   virtual void Action(TClassMenuItem *classmenuitem);
   void Action(TMethod *method) { Action(fSelectedObject, method); }
   void Action(TToggle *toggle) { Action(fSelectedObject, toggle); }
   virtual const char *CreateArgumentTitle(TMethodArg *argument);
   virtual const char *CreateDialogTitle(TObject *object, TFunction *method);
   virtual const char *CreatePopupTitle(TObject *object );
   virtual void Execute(const char *method,  const char *params, Int_t *error=nullptr) { TObject::Execute(method, params, error); }
   virtual void Execute(TMethod *method, TObjArray *params, Int_t *error=nullptr) { TObject::Execute(method, params, error); }
   virtual void Execute(TObject *object, TFunction *method, const char *params);
   virtual void Execute(TObject *object, TFunction *method, TObjArray *params);
   void Execute(const char *params) { Execute(fCalledObject, fSelectedMethod, params); }
   void Execute(TObjArray *params) { Execute(fCalledObject, fSelectedMethod, params); }
   virtual TBrowser *GetBrowser() { return fBrowser; }
   virtual TContextMenuImp *GetContextMenuImp() { return fContextMenuImp; }
   virtual TVirtualPad *GetSelectedCanvas() { return fSelectedCanvas; }
   virtual TFunction *GetSelectedMethod() { return fSelectedMethod; }
   virtual TObject *GetSelectedObject() { return fSelectedObject; }
   virtual TObject *GetCalledObject() { return fCalledObject; }
   virtual TClassMenuItem *GetSelectedMenuItem() { return fSelectedMenuItem; }
   virtual TVirtualPad *GetSelectedPad() { return fSelectedPad; }
   virtual void Popup(Int_t x, Int_t y, TObject *obj, TVirtualPad *c=nullptr, TVirtualPad *p=nullptr); // Create menu from canvas
   virtual void Popup(Int_t x, Int_t y, TObject *obj, TBrowser *b);  // Create menu from Browser
   virtual void SetCanvas(TVirtualPad *c) { fSelectedCanvas = c; }
   virtual void SetBrowser(TBrowser *b) { fBrowser = b; }
   virtual void SetMethod(TFunction *m) { fSelectedMethod = m; }
   virtual void SetCalledObject(TObject *o) { fCalledObject = o; }
   virtual void SetSelectedMenuItem(TClassMenuItem *mi) { fSelectedMenuItem = mi; }
   virtual void SetNameTitle(const char *name, const char *title) { TNamed::SetNameTitle(name, title); }
   virtual void SetObject(TObject *o) { fSelectedObject = o; }
   virtual void SetPad(TVirtualPad *p) { fSelectedPad = p; }

   ClassDef(TContextMenu,0)  //Context sensitive popup menu
};

#endif
