// @(#)root/base:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   15/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TQObject
#define ROOT_TQObject

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This is the ROOT implementation of the Qt object communication       //
// mechanism (see also http://www.troll.no/qt/metaobjects.html)         //
//                                                                      //
// Signals and slots are used for communication between objects.        //
// When an object has changed in some way that might be interesting     //
// for the outside world, it emits a signal to tell whoever is          //
// listening. All slots that are connected to this signal will be       //
// activated (called).  It is even possible to connect a signal         //
// directly to  another signal (this will emit the second signal        //
// immediately whenever the first is emitted.) There is no limitation   //
// on the number of slots that can be connected to a signal.            //
// The slots will be activated in the order they were connected         //
// to the signal. This mechanism allows objects to be easily reused,    //
// because the object that emits a signal does not need to know         //
// to what the signals are connected to.                                //
// Together, signals and slots make up a powerfull component            //
// programming mechanism.                                               //
//                                                                      //
// This implementation is provided by                                   //
// Valeriy Onuchin (onuchin@sirius.ihep.su).                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TList.h"
#include "TString.h"
#include "TVirtualQConnection.h"

class TClass;

R__EXTERN void *gTQSender;   // the latest sender object

class TQObject {

protected:
   TList   *fListOfSignals;        //! list of signals from this object
   TList   *fListOfConnections;    //! list of connections to this object
   Bool_t   fSignalsBlocked;       //! flag used for suppression of signals

   static Bool_t fgAllSignalsBlocked;  // flag used for suppression of all signals

   virtual void       *GetSender() { return this; }
   virtual const char *GetSenderClassName() const { return ""; }


   static Bool_t ConnectToClass(TQObject *sender,
                                const char *signal,
                                TClass *receiver_class,
                                void *receiver,
                                const char *slot);

   static Bool_t ConnectToClass(const char *sender_class,
                                const char *signal,
                                TClass *receiver_class,
                                void *receiver,
                                const char *slot);

   static Int_t CheckConnectArgs(TQObject *sender,
                                 TClass *sender_class, const char *signal,
                                 TClass *receiver_class, const char *slot);

   static TString CompressName(const char *method_name);

private:
   TQObject(const TQObject &) = delete;
   TQObject& operator=(const TQObject &) = delete;

public:
   TQObject();
   virtual ~TQObject();

   TList   *GetListOfClassSignals() const;
   TList   *GetListOfSignals() const { return fListOfSignals; }
   TList   *GetListOfConnections() const { return fListOfConnections; }

   Bool_t   AreSignalsBlocked() const { return fSignalsBlocked; }
   Bool_t   BlockSignals(Bool_t b)
            { Bool_t ret = fSignalsBlocked; fSignalsBlocked = b; return ret; }

   void  CollectClassSignalLists(TList& list, TClass* cls);

   ///////////////////////////////////////////////////////////////////////////////
   /// Emit a signal with a varying number of arguments.
   ///
   template <typename... T> void EmitVA(const char *signal_name, Int_t /* nargs */, const T&... params)
   {
      // Activate signal with variable argument list.
      // For internal use and for var arg EmitVA() in RQ_OBJECT.h.

      if (fSignalsBlocked || AreAllSignalsBlocked())
         return;

      TList classSigLists;
      CollectClassSignalLists(classSigLists, IsA());

      if (classSigLists.IsEmpty() && !fListOfSignals)
         return;

      TString signal = CompressName(signal_name);

      TVirtualQConnection *connection = 0;

      // execute class signals
      TList *sigList;
      TIter  nextSigList(&classSigLists);
      while ((sigList = (TList*) nextSigList())) {
         TIter nextcl((TList*) sigList->FindObject(signal));
         while ((connection = static_cast<TVirtualQConnection*>(nextcl()))) {
            gTQSender = GetSender();
            connection->SetArgs(params...);
            connection->SendSignal();
         }
      }
      if (!fListOfSignals)
         return;

      // execute object signals
      TIter next((TList*) fListOfSignals->FindObject(signal));
      while (fListOfSignals && (connection = static_cast<TVirtualQConnection*>(next()))) {
         gTQSender = GetSender();
         connection->SetArgs(params...);
         connection->SendSignal();
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Activate signal with single parameter.
   /// Example:
   /// ~~~ {.cpp}
   ///   theButton->Emit("Progress(Long64_t)",processed)
   /// ~~~
   ///
   /// If we call Emit with an array of the parameters, they should be converted
   /// to long type.
   /// Example:
   /// ~~~ {.cpp}
   ///    TQObject *processor; // data processor
   ///    TH1F     *hist;      // filled with processor results
   ///
   ///    processor->Connect("Evaluated(Float_t,Float_t)",
   ///                       "TH1F",hist,"Fill12(Axis_t,Axis_t)");
   ///
   ///    Long_t args[2];
   ///    args[0] = (Long_t)processor->GetValue(1);
   ///    args[1] = (Long_t)processor->GetValue(2);
   ///
   ///    processor->Emit("Evaluated(Float_t,Float_t)",args);
   /// ~~~
   template <typename T> void Emit(const char *signal, const T& arg) {
      Int_t placeholder = 0;
      EmitVA(signal, placeholder, arg);
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Acitvate signal without args.
   /// Example:
   ///          theButton->Emit("Clicked()");
   void  Emit(const char *signal) { EmitVA(signal, (Int_t) 0); }

   Bool_t Connect(const char *signal,
                  const char *receiver_class,
                  void *receiver,
                  const char *slot);

   Bool_t Disconnect(const char *signal = 0,
                     void *receiver = 0,
                     const char *slot = 0);

   virtual void   HighPriority(const char *signal_name,
                               const char *slot_name = 0);

   virtual void   LowPriority(const char *signal_name,
                              const char *slot_name = 0);

   virtual Bool_t HasConnection(const char *signal_name) const;
   virtual Int_t  NumberOfSignals() const;
   virtual Int_t  NumberOfConnections() const;
   virtual void   Connected(const char * /*signal_name*/) { }
   virtual void   Disconnected(const char * /*signal_name*/) { }

   virtual void   Destroyed()
                  { Emit("Destroyed()"); }                 // *SIGNAL*
   virtual void   ChangedBy(const char *method)
                  { Emit("ChangedBy(char*)", method); }    // *SIGNAL*
   virtual void   Message(const char *msg)
                  { Emit("Message(char*)", msg); }         // *SIGNAL*

   static Bool_t  Connect(TQObject *sender,
                          const char *signal,
                          const char *receiver_class,
                          void *receiver,
                          const char *slot);

   static Bool_t  Connect(const char *sender_class,
                          const char *signal,
                          const char *receiver_class,
                          void *receiver,
                          const char *slot);

   static Bool_t  Disconnect(TQObject *sender,
                             const char *signal = 0,
                             void *receiver = 0,
                             const char *slot = 0);

   static Bool_t  Disconnect(const char *class_name,
                             const char *signal,
                             void *receiver = 0,
                             const char *slot = 0);

   static Bool_t  AreAllSignalsBlocked();
   static Bool_t  BlockAllSignals(Bool_t b);

   ClassDef(TQObject,1) //Base class for object communication mechanism
};


class TQObjSender : public TQObject {

protected:
   void    *fSender;        //delegation object
   TString  fSenderClass;   //class name of delegation object

   virtual void       *GetSender() { return fSender; }
   virtual const char *GetSenderClassName() const { return fSenderClass; }

private:
   TQObjSender(const TQObjSender&);            // not implemented
   TQObjSender& operator=(const TQObjSender&); // not implemented

public:
   TQObjSender() : TQObject(), fSender(0), fSenderClass() { }
   virtual ~TQObjSender() { Disconnect(); }

   virtual void SetSender(void *sender) { fSender = sender; }
   void SetSenderClassName(const char *sclass = "") { fSenderClass = sclass; }

   ClassDef(TQObjSender,0) //Used to "delegate" TQObject functionality
                           //to interpreted classes, see also RQ_OBJECT.h
};



// Global function which simplifies making connections in interpreted
// ROOT session
//
//  ConnectCINT      - connects to interpreter(CINT) command

extern Bool_t ConnectCINT(TQObject *sender, const char *signal,
                          const char *slot);

#ifdef G__DICTIONARY
// This include makes it possible to have a single connection
// from all objects of the same class but is only needed in
// the dictionary.
#include "TQClass.h"
#endif


//---- ClassImpQ macro ----------------------------------------------
//
// This macro used to correspond to the ClassImp macro and should be used
// for classes derived from TQObject instead of the ClassImp macro.
// This macro makes it possible to have a single connection from
// all objects of the same class.
// *** It is now obsolete ***

#define ClassImpQ(name) \
   ClassImp(name)

#endif
