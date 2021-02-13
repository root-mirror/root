// Author: Enric Tejedor CERN  07/2020
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPyDispatcher
#define ROOT_TPyDispatcher

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyDispatcher                                                            //
//                                                                          //
// Dispatcher for C++ callbacks into Python code.                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// ROOT
#include "TObject.h"

class TDNDData;
class TEveDigitSet;
class TEveElement;
class TEveTrack;
class TEveWindow;
class TGFrame;
class TGListTreeItem;
class TGMdiFrame;
class TGLPhysicalShape;
class TGShutterItem;
class TGLVEntry;
class TGLViewerBase;
class TGVFileSplitter;
class TList;
class TObject;
class TPad;
class TProofProgressInfo;
class TQCommand;
class TSlave;
class TSocket;
class TVirtualPad;

struct Event_t;

// Python
struct _object;
typedef _object PyObject;

class TPyDispatcher : public TObject {
public:
   TPyDispatcher(PyObject *callable);
   TPyDispatcher(const TPyDispatcher &);
   TPyDispatcher &operator=(const TPyDispatcher &);
   ~TPyDispatcher();

public:
#ifndef __CINT__
   PyObject *DispatchVA(const char *format = 0, ...);
#else
   PyObject *DispatchVA(const char *format, ...);
#endif
   PyObject *DispatchVA1(const char *clname, void *obj, const char *format, ...);

   // pre-defined dispatches, same as per TQObject::Emit(); note that
   // Emit() maps exclusively to this set, so several builtin types (e.g.
   // Int_t, Bool_t, Float_t, etc.) have been omitted here
   PyObject *Dispatch() { return DispatchVA(0); }
   PyObject *Dispatch(const char *param) { return DispatchVA("s", param); }
   PyObject *Dispatch(Double_t param) { return DispatchVA("d", param); }
   PyObject *Dispatch(Long_t param) { return DispatchVA("l", param); }
   PyObject *Dispatch(Long64_t param) { return DispatchVA("L", param); }

   // further selection of pre-defined, existing dispatches
   PyObject *Dispatch(Bool_t param) { return DispatchVA("i", param); }
   PyObject *Dispatch(char *param) { return DispatchVA("s", param); }
   PyObject *Dispatch(const char *text, Int_t len) { return DispatchVA("si", text, len); }
   PyObject *Dispatch(Int_t param) { return DispatchVA("i", param); }
   PyObject *Dispatch(Int_t x, Int_t y) { return DispatchVA("ii", x, y); }
   PyObject *Dispatch(ULong_t param) { return DispatchVA("k", param); }
   // ULong_t also for Handle_t (and Window_t, etc. ... )

   PyObject *Dispatch(Event_t *event) { return DispatchVA1("Event_t", event, 0); }
   PyObject *Dispatch(Event_t *event, ULong_t wid) { return DispatchVA1("Event_t", event, "k", wid); }
   PyObject *Dispatch(TEveDigitSet *qs, Int_t idx) { return DispatchVA1("TEveDigitSet", qs, "i", idx); }
   PyObject *Dispatch(TEveElement *el) { return DispatchVA1("TEveElement", el, 0); }
   PyObject *Dispatch(TEveTrack *et) { return DispatchVA1("TEveTrack", et, 0); }
   PyObject *Dispatch(TEveWindow *window) { return DispatchVA1("TEveWindow", window, 0); }
   PyObject *Dispatch(TGFrame *frame) { return DispatchVA1("TGFrame", frame, 0); }
   PyObject *Dispatch(TGFrame *frame, Int_t btn) { return DispatchVA1("TGFrame", frame, "i", btn); }
   PyObject *Dispatch(TGFrame *frame, Int_t btn, Int_t x, Int_t y)
   {
      return DispatchVA1("TGFrame", frame, "iii", btn, x, y);
   }
   PyObject *Dispatch(TGFrame *frame, UInt_t keysym, UInt_t mask)
   {
      return DispatchVA1("TGFrame", frame, "II", keysym, mask);
   }
   PyObject *Dispatch(TGListTreeItem *entry) { return DispatchVA1("TGListTreeItem", entry, 0); }
   PyObject *Dispatch(TGListTreeItem *entry, UInt_t mask) { return DispatchVA1("TGListTreeItem", entry, "I", mask); }
   PyObject *Dispatch(TGListTreeItem *entry, UInt_t keysym, UInt_t mask)
   {
      return DispatchVA1("TGListTreeItem", entry, "II", keysym, mask);
   }
   PyObject *Dispatch(TGListTreeItem *entry, Int_t btn) { return DispatchVA1("TGListTreeItem", entry, "i", btn); }
   PyObject *Dispatch(TGListTreeItem *entry, Int_t btn, Int_t x, Int_t y)
   {
      return DispatchVA1("TGListTreeItem", entry, "iii", btn, x, y);
   }
   PyObject *Dispatch(TGLVEntry *entry, Int_t btn) { return DispatchVA1("TGLVEntry", entry, "i", btn); }
   PyObject *Dispatch(TGLVEntry *entry, Int_t btn, Int_t x, Int_t y)
   {
      return DispatchVA1("TGLVEntry", entry, "iii", btn, x, y);
   }
   PyObject *Dispatch(TGLViewerBase *viewer) { return DispatchVA1("TGLViewerBase", viewer, 0); }
   PyObject *Dispatch(TGLPhysicalShape *shape) { return DispatchVA1("TGLPhysicalShape", shape, 0); }
   PyObject *Dispatch(TGLPhysicalShape *shape, UInt_t u1, UInt_t u2)
   {
      return DispatchVA1("TGLPhysicalShape", shape, "II", u1, u2);
   }
   PyObject *Dispatch(TGMdiFrame *frame) { return DispatchVA1("TGMdiFrame", frame, 0); }
   PyObject *Dispatch(TGShutterItem *item) { return DispatchVA1("TGShutterItem", item, 0); }
   PyObject *Dispatch(TGVFileSplitter *frame) { return DispatchVA1("TGVFileSplitter", frame, 0); }
   PyObject *Dispatch(TList *objs) { return DispatchVA1("TList", objs, 0); }
   PyObject *Dispatch(TObject *obj) { return DispatchVA1("TObject", obj, 0); }
   PyObject *Dispatch(TObject *obj, Bool_t check) { return DispatchVA1("TObject", obj, "i", check); }
   PyObject *Dispatch(TObject *obj, UInt_t state) { return DispatchVA1("TObject", obj, "I", state); }
   PyObject *Dispatch(TObject *obj, UInt_t button, UInt_t state)
   {
      return DispatchVA1("TObject", obj, "II", button, state);
   }
   PyObject *Dispatch(TSocket *sock) { return DispatchVA1("TSocket", sock, 0); }
   PyObject *Dispatch(TVirtualPad *pad) { return DispatchVA1("TVirtualPad", pad, 0); }

   PyObject *Dispatch(TPad *selpad, TObject *selected, Int_t event);
   PyObject *Dispatch(Int_t event, Int_t x, Int_t y, TObject *selected);
   PyObject *Dispatch(TVirtualPad *pad, TObject *obj, Int_t event);
   PyObject *Dispatch(TGListTreeItem *item, TDNDData *data);
   PyObject *Dispatch(const char *name, const TList *attr);

   // for PROOF
   PyObject *Dispatch(const char *msg, Bool_t all) { return DispatchVA("si", msg, all); }
   PyObject *Dispatch(Long64_t total, Long64_t processed) { return DispatchVA("LL", total, processed); }
   PyObject *Dispatch(Long64_t total, Long64_t processed, Long64_t bytesread, Float_t initTime, Float_t procTime,
                      Float_t evtrti, Float_t mbrti)
   {
      return DispatchVA("LLLffff", total, processed, bytesread, initTime, procTime, evtrti, mbrti);
   }
   PyObject *Dispatch(Long64_t total, Long64_t processed, Long64_t bytesread, Float_t initTime, Float_t procTime,
                      Float_t evtrti, Float_t mbrti, Int_t actw, Int_t tses, Float_t eses)
   {
      return DispatchVA("LLLffffiif", total, processed, bytesread, initTime, procTime, evtrti, mbrti, actw, tses, eses);
   }
   PyObject *Dispatch(const char *sel, Int_t sz, Long64_t fst, Long64_t ent)
   {
      return DispatchVA("siLL", sel, sz, fst, ent);
   }
   PyObject *Dispatch(const char *msg, Bool_t status, Int_t done, Int_t total)
   {
      return DispatchVA("siii", msg, status, done, total);
   }

   PyObject *Dispatch(TSlave *slave, Long64_t total, Long64_t processed)
   {
      return DispatchVA1("TSlave", slave, "LL", total, processed);
   }
   PyObject *Dispatch(TProofProgressInfo *pi) { return DispatchVA1("TProofProgressInfo", pi, 0); }
   PyObject *Dispatch(TSlave *slave) { return DispatchVA("TSlave", slave, 0); }
   PyObject *Dispatch(TSlave *slave, TProofProgressInfo *pi);

private:
   PyObject *fCallable; //! callable object to be dispatched
};

#endif
