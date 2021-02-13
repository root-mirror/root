// @(#)root/base:$Id: 5d6810ad46b864564f576f88aa9b154789d91d48 $
// Author: Valeriy Onuchin & Fons Rademakers   15/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TQObject
\ingroup Base

This is the ROOT implementation of the Qt object communication
mechanism (see also http://www.troll.no/qt/metaobjects.html)

Signals and slots are used for communication between objects.
When an object has changed in some way that might be interesting
for the outside world, it emits a signal to tell whoever is
listening. All slots that are connected to this signal will be
activated (called). It is even possible to connect a signal
directly to another signal (this will emit the second signal
immediately whenever the first is emitted.) There is no limitation
on the number of slots that can be connected to a signal.
The slots will be activated in the order they were connected
to the signal. This mechanism allows objects to be easily reused,
because the object that emits a signal does not need to know
to which objects the signals are connected.
Together, signals and slots make up a powerfull component
programming mechanism.

### Signals

~~~ {.cpp}
   Destroyed()
~~~
Signal emitted when object is destroyed.
This signal could be connected to some garbage-collector object.

~~~ {.cpp}
   ChangedBy(const char *method_name)
~~~
This signal is emitted when some important data members of
the object were changed. method_name parameter can be used
as an identifier of the modifier method.

~~~ {.cpp}
   Message(const char *msg)
~~~

General purpose message signal
*/

#include "TQObject.h"
#include "TQConnection.h"
#include "THashList.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TMethod.h"
#include "TBaseClass.h"
#include "TDataType.h"
#include "TInterpreter.h"
#include "TQClass.h"
#include "TError.h"
#include <iostream>
#include "RQ_OBJECT.h"
#include "TVirtualMutex.h"
#include "RConfigure.h"

void *gTQSender; // A pointer to the object that sent the last signal.
                 // Getting access to the sender might be practical
                 // when many signals are connected to a single slot.

Bool_t TQObject::fgAllSignalsBlocked = kFALSE;


ClassImpQ(TQObject)
ClassImpQ(TQObjSender)
ClassImpQ(TQClass)

////////////////////////////////////////////////////////////////////////////////
/// Removes "const" words and blanks from full (with prototype)
/// method name and resolve any typedefs in the method signature.
/// If a null or empty string is passed in, an empty string
/// is returned.
///
/// Example:
/// ~~~ {.cpp}
///   CompressName(" Draw(const char *, const char *,
///                  Option_t * , Int_t , Int_t)");
/// ~~~
/// returns the string "Draw(char*,char*,char*,int,int)".

TString TQObject::CompressName(const char *method_name)
{
   TString res(method_name);
   if (res.IsNull())
      return res;

   {
      static TVirtualMutex *  lock = 0;
      R__LOCKGUARD2(lock);

      static TPMERegexp *constRe = 0, *wspaceRe = 0;
      if (constRe == 0) {
         constRe  = new TPMERegexp("(?<=\\(|\\s|,|&|\\*)const(?=\\s|,|\\)|&|\\*)", "go");
         wspaceRe = new TPMERegexp("\\s+(?=([^\"]*\"[^\"]*\")*[^\"]*$)", "go");
      }
      constRe ->Substitute(res, "");
      wspaceRe->Substitute(res, "");
   }

   TStringToken methargs(res, "\\(|\\)", kTRUE);

   methargs.NextToken();
   res = methargs;
   res += "(";

   methargs.NextToken();
   TStringToken arg(methargs, ",");
   while (arg.NextToken())
   {
      Int_t  pri = arg.Length() - 1;
      Char_t prc = 0;
      if (arg[pri] == '*' || arg[pri] == '&') {
         prc = arg[pri];
         arg.Remove(pri);
      }
      TDataType *dt = gROOT->GetType(arg.Data());
      if (dt) {
         res += dt->GetFullTypeName();
      } else {
         res += arg;
      }
      if (prc)          res += prc;
      if (!arg.AtEnd()) res += ",";
   }
   res += ")";
   return res;
}

namespace {

////////////////////////////////////////////////////////////////////////////////
/// Almost the same as TClass::GetMethodWithPrototype().

TMethod *GetMethodWithPrototype(TClass *cl, const char *method,
                                const char *proto, Int_t &nargs)
{
   nargs = 0;

   if (!gInterpreter || cl == 0) return 0;

   TMethod *m = cl->GetMethodWithPrototype(method,proto);
   if (m) nargs = m->GetNargs();
   return m;
}

////////////////////////////////////////////////////////////////////////////////
/// Almost the same as TClass::GetMethod().

static TMethod *GetMethod(TClass *cl, const char *method, const char *params)
{
   if (!gInterpreter || cl == 0) return 0;
   return cl->GetMethod(method,params);
}

}

////////////////////////////////////////////////////////////////////////////////
/// Checking of consitency of sender/receiver methods/arguments.
/// Returns -1 on error, otherwise number or arguments of signal function.
/// Static method.

Int_t TQObject::CheckConnectArgs(TQObject *sender,
                                 TClass *sender_class, const char *signal,
                                 TClass *receiver_class, const char *slot)
{
   char *signal_method = new char[strlen(signal)+1];
   if (signal_method) strcpy(signal_method, signal);

   char *signal_proto;
   char *tmp;

   if ((signal_proto = strchr(signal_method,'('))) {
      // substitute first '(' symbol with '\0'
      *signal_proto++ = '\0';
      // substitute last ')' symbol with '\0'
      if ((tmp = strrchr(signal_proto,')'))) *tmp = '\0';
   }

   if (!signal_proto) signal_proto = (char*)""; // avoid zero strings

   // if delegation object TQObjSender is used get the real sender class
   if (sender && sender_class == TQObjSender::Class()) {
      sender_class = TClass::GetClass(sender->GetSenderClassName());
      if (!sender_class) {
         ::Error("TQObject::CheckConnectArgs", "for signal/slot consistency\n"
                 "checking need to specify class name as argument to "
                 "RQ_OBJECT macro");
         delete [] signal_method;
         return -1;
      }
   }

   Int_t nargs;
   TMethod *signalMethod = GetMethodWithPrototype(sender_class,
                                                  signal_method,
                                                  signal_proto,
                                                  nargs);
   if (!signalMethod) {
      ::Error("TQObject::CheckConnectArgs",  "signal %s::%s(%s) does not exist",
              sender_class->GetName(), signal_method, signal_proto);
      delete [] signal_method;
      return -1;
   }
   Int_t nsigargs = nargs;

#if defined(CHECK_COMMENT_STRING)
   const char *comment = 0;
   if (signalMethod != (TMethod *) -1)   // -1 in case of interpreted class
      comment = signalMethod->GetCommentString();

   if (!comment || !comment[0] || strstr(comment,"*SIGNAL")){
      ::Error("TQObject::CheckConnectArgs",
              "signal %s::%s(%s), to declare signal use comment //*SIGNAL*",
      sender_class->GetName(), signal_method, signal_proto);
      delete [] signal_method;
      return -1;
   }
#endif

   // cleaning
   delete [] signal_method;

   char *slot_method = new char[strlen(slot)+1];
   if (slot_method) strcpy(slot_method, slot);

   char *slot_proto;
   char *slot_params = 0;

   if ((slot_proto = strchr(slot_method,'('))) {

      // substitute first '(' symbol with '\0'
      *slot_proto++ = '\0';

      // substitute last ')' symbol with '\0'
      if ((tmp = strrchr(slot_proto,')'))) *tmp = '\0';
   }

   if (!slot_proto) slot_proto = (char*)"";     // avoid zero strings
   if ((slot_params = strchr(slot_proto,'='))) *slot_params = ' ';

   TFunction *slotMethod = 0;
   if (!receiver_class) {
      // case of slot_method is compiled/intrepreted function
      slotMethod = gROOT->GetGlobalFunction(slot_method,0,kFALSE);
   } else {
      slotMethod  = !slot_params ?
                          GetMethodWithPrototype(receiver_class,
                                                 slot_method,
                                                 slot_proto,
                                                 nargs) :
                          GetMethod(receiver_class,
                                    slot_method, slot_params);
   }

   if (!slotMethod) {
      if (!slot_params) {
         ::Error("TQObject::CheckConnectArgs", "slot %s(%s) does not exist",
                 receiver_class ? Form("%s::%s", receiver_class->GetName(),
                 slot_method) : slot_method, slot_proto);
      } else {
         ::Error("TQObject::CheckConnectArgs", "slot %s(%s) does not exist",
                 receiver_class ? Form("%s::%s", receiver_class->GetName(),
                 slot_method) : slot_method, slot_params);
      }
      delete [] slot_method;
      return -1;
   }

#if defined(CHECK_ARGS_NUMBER)
   if (slotMethod != (TMethod *) -1 && slotMethod->GetNargsOpt() >= 0 &&
       nsigargs < (slotMethod->GetNargs() - slotMethod->GetNargsOpt())) {
      ::Error("TQObject::CheckConnectArgs",
              "inconsistency in numbers of arguments");
      delete [] slot_method;
      return -1;
   }
#endif

   // cleaning
   delete [] slot_method;

   return nsigargs;
}

/** \class TQConnectionList
   TQConnectionList is the named list of connections,
   see also TQConnection class.
*/

class TQConnectionList : public TList {

private:
   Int_t   fSignalArgs;    // number of arguments in signal function

public:
   TQConnectionList(const char *name, Int_t nsigargs) : TList()
      { fName = name; fSignalArgs = nsigargs; }
   virtual ~TQConnectionList();

   Bool_t Disconnect(void *receiver=0, const char *slot_name=0);
   Int_t  GetNargs() const { return fSignalArgs; }
   void   ls(Option_t *option = "") const;
};

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TQConnectionList::~TQConnectionList()
{
   TIter next(this);
   TQConnection *connection;

   while ((connection = (TQConnection*)next())) {
      // remove this from feed back reference list
      connection->Remove(this);
      if (connection->IsEmpty()) delete connection;
   }
   Clear("nodelete");
}

////////////////////////////////////////////////////////////////////////////////
/// Remove connection from the list. For more info see
/// TQObject::Disconnect()

Bool_t TQConnectionList::Disconnect(void *receiver, const char *slot_name)
{
   TQConnection *connection = 0;
   Bool_t return_value = kFALSE;

   TObjLink *lnk = FirstLink();
   TObjLink *savlnk; // savlnk is used when link is deleted

   while (lnk) {
      connection = (TQConnection*)lnk->GetObject();
      const char *name = connection->GetName();
      void *obj = connection->GetReceiver();

      if (!slot_name || !slot_name[0]
                     || !strcmp(name,slot_name)) {

         if (!receiver || (receiver == obj)) {
            return_value = kTRUE;
            savlnk = lnk->Next();   // keep next link ..
            Remove(lnk);
            lnk = savlnk;           // current link == saved ...
            connection->Remove(this);      // remove back reference
            if (connection->IsEmpty()) SafeDelete(connection);
            continue;               // .. continue from saved link
         }
      }
      lnk = lnk->Next();
   }
   return return_value;
}

////////////////////////////////////////////////////////////////////////////////
/// List signal name and list all connections in this signal list.

void TQConnectionList::ls(Option_t *option) const
{
   std::cout <<  "TQConnectionList:" << "\t" << GetName() << std::endl;
   ((TQConnectionList*)this)->R__FOR_EACH(TQConnection,Print)(option);
}


////////////////////////////////////////////////////////////////////////////////
/// TQObject Constructor.
/// Comment:
///  - In order to minimize memory allocation fListOfSignals and
///    fListOfConnections are allocated only if it is neccesary
///  - When fListOfSignals/fListOfConnections are empty they will
///    be deleted

TQObject::TQObject()
{
   fListOfSignals     = 0;
   fListOfConnections = 0;
   fSignalsBlocked    = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// TQObject Destructor.
///    - delete all connections and signal list

TQObject::~TQObject()
{
   if (!gROOT) return;

   Destroyed();   // emit "Destroyed()" signal

   if (fListOfSignals) {
      fListOfSignals->Delete();
      SafeDelete(fListOfSignals);   // delete list of signals
   }

   // loop over all connections and remove references to this object
   if (fListOfConnections) {
      TIter next_connection(fListOfConnections);
      TQConnection *connection;

      while ((connection = (TQConnection*)next_connection())) {
         TIter next_list(connection);
         TQConnectionList *list;
         while ((list = (TQConnectionList*)next_list())) {
            list->Remove(connection);
            if (list->IsEmpty()) SafeDelete(list);
         }
      }
      SafeDelete(fListOfConnections);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pointer to list of signals of this class.

TList *TQObject::GetListOfClassSignals() const
{
   TQClass *qcl = 0;

   qcl = dynamic_cast<TQClass*>(IsA());

   return qcl ? qcl->fListOfSignals : 0; //!!
}

////////////////////////////////////////////////////////////////////////////////
/// Collect class signal lists from class cls and all its
/// base-classes.
///
/// The recursive traversal is not performed for classes not
/// deriving from TQClass.

void TQObject::CollectClassSignalLists(TList& list, TClass* cls)
{
   TQClass *qcl = dynamic_cast<TQClass*>(cls);
   if (qcl)
   {
      if (qcl->fListOfSignals)
         list.Add(qcl->fListOfSignals);

      // Descend into base-classes.
      TIter       next_base_class(cls->GetListOfBases());
      TBaseClass *base;
      while ((base = (TBaseClass*) next_base_class()))
      {
         CollectClassSignalLists(list, base->GetClassPointer());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// 1. If slot_name = 0 => makes signal defined by the signal_name
///    to be the first in the fListOfSignals, this decreases
///    the time for lookup.
/// 2. If slot_name != 0 => makes slot defined by the slot_name
///    to be executed first when signal_name is emitted.
/// Signal name is not compressed.

void TQObject::HighPriority(const char *signal_name, const char *slot_name)
{
   if (!fListOfSignals) return;
   TQConnectionList *clist = (TQConnectionList*)
      fListOfSignals->FindObject(signal_name);

   if (!clist)  return;      // not found
   if (!slot_name)  {        // update list of signal lists
      fListOfSignals->Remove(clist);   // remove and add first
      fListOfSignals->AddFirst(clist);
      return;
   } else {                   // slot_name != 0 , update signal list
      TQConnection *con = (TQConnection*) clist->FindObject(slot_name);
      if (!con) return;       // not found
      clist->Remove(con);     // remove and add as first
      clist->AddFirst(con);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// 1. If slot_name = 0 => makes signal defined by the signal_name
///    to be the last in the fListOfSignals, this increase the time
///    for lookup.
/// 2. If slot_name != 0 => makes slot defined by the slot_name
///    to  be executed last when signal_name is emitted.
/// Signal name is not compressed.

void TQObject::LowPriority(const char *signal_name, const char *slot_name)
{
   if (!fListOfSignals) return;
   TQConnectionList *clist = (TQConnectionList*)
      fListOfSignals->FindObject(signal_name);

   if (!clist)   return;
   if (!slot_name)  {
      fListOfSignals->Remove(clist);   // remove and add first
      fListOfSignals->AddLast(clist);
      return;
   } else  {                  // slot_name != 0 , update signal list
      TQConnection *con = (TQConnection*) clist->FindObject(slot_name);
      if (!con) return;
      clist->Remove(con);     // remove and add as last
      clist->AddLast(con);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if there is any object connected to this signal.
/// Only checks for object signals.

Bool_t TQObject::HasConnection(const char *signal_name) const
{
   if (!fListOfSignals)
      return kFALSE;

   TString signal = CompressName(signal_name);

   return (fListOfSignals->FindObject(signal) != 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of signals for this object.
/// Only checks for object signals.

Int_t TQObject::NumberOfSignals() const
{
   if (fListOfSignals)
      return fListOfSignals->GetSize();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of connections for this object.

Int_t TQObject::NumberOfConnections() const
{
   if (fListOfConnections)
      return fListOfConnections->GetSize();
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Create connection between sender and receiver.
/// Receiver class needs to have a dictionary.

Bool_t TQObject::ConnectToClass(TQObject *sender,
                                const char *signal,
                                TClass *cl,
                                void *receiver,
                                const char *slot)
{
   // sender should be TQObject
   if (!sender->IsA()->InheritsFrom(TQObject::Class()))
      return kFALSE;

   // remove "const" and strip blanks
   TString signal_name = CompressName(signal);
   TString slot_name   = CompressName(slot);

   // check consitency of signal/slot methods/args
   Int_t nsigargs;
   if ((nsigargs = CheckConnectArgs(sender, sender->IsA(), signal_name, cl, slot_name)) == -1)
      return kFALSE;

   if (!sender->fListOfSignals)
      sender->fListOfSignals = new THashList();

   TQConnectionList *clist = (TQConnectionList*)
      sender->fListOfSignals->FindObject(signal_name);

   if (!clist) {
      clist = new TQConnectionList(signal_name, nsigargs);
      sender->fListOfSignals->Add(clist);
   }

   TIter next(clist);
   TQConnection *connection = 0;

   while ((connection = (TQConnection*)next())) {
      if (!strcmp(slot_name,connection->GetName()) &&
          (receiver == connection->GetReceiver())) break;
   }

   if (!connection)
      connection = new TQConnection(cl, receiver, slot_name);

   // check to prevent multiple entries
   if (!clist->FindObject(connection)) {
      clist->Add(connection);
      if (!connection->FindObject(clist)) connection->Add(clist);
      sender->Connected(signal_name);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This method allows to make connection from any object
/// of the same class to the receiver object.
/// Receiver class needs to have a dictionary.

Bool_t TQObject::ConnectToClass(const char *class_name,
                                const char *signal,
                                TClass *cl,
                                void *receiver,
                                const char *slot)
{
   TClass *sender = TClass::GetClass(class_name);

   // sender class should be TQObject (i.e. TQClass)
   if (!sender || !sender->IsA()->InheritsFrom(TQObject::Class()))
      return kFALSE;

   TList *slist = ((TQClass*)sender)->fListOfSignals;
   TString signal_name = CompressName(signal);
   TString slot_name   = CompressName(slot);

   // check consitency of signal/slot methods/args
   Int_t nsigargs;
   if ((nsigargs = CheckConnectArgs(0, sender, signal_name, cl, slot_name)) == -1)
      return kFALSE;

   if (!slist)
      ((TQClass*)sender)->fListOfSignals = slist = new THashList();

   TQConnectionList *clist = (TQConnectionList*) slist->FindObject(signal_name);

   if (!clist) {
      clist = new TQConnectionList(signal_name, nsigargs);
      slist->Add(clist);
   }

   TQConnection *connection = 0;
   TIter next(clist);

   while ((connection = (TQConnection*)next())) {
      if (!strcmp(slot_name,connection->GetName()) &&
          (receiver == connection->GetReceiver())) break;
   }

   if (!connection)
      connection = new TQConnection(cl, receiver, slot_name);

   // check to prevent multiple entries
   if (!clist->FindObject(connection)) {
      clist->Add(connection);
      if (!connection->FindObject(clist)) connection->Add(clist);
      ((TQClass*)sender)->Connected(signal_name);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create connection between sender and receiver.
/// Signal and slot string must have a form:
///    "Draw(char*, Option_t* ,Int_t)"
/// All blanks and "const" words will be removed,
///
/// cl != 0 - class name, it can be class with or
///           without dictionary, e.g interpreted class.
/// Example:
/// ~~~ {.cpp}
///       TGButton *myButton;
///       TH2F     *myHist;
///
///       TQObject::Connect(myButton,"Clicked()",
///                         "TH2F", myHist,"Draw(Option_t*)");
/// ~~~
/// cl == 0 - corresponds to function (interpereted or global)
///           the name of the function is defined by the slot string,
///           parameter receiver should be 0.
/// Example:
/// ~~~ {.cpp}
///       TGButton *myButton;
///       TH2F     *myHist;
///
///       TQObject::Connect(myButton,"Clicked()",
///                         0, 0,"hsimple()");
/// ~~~
/// Warning:
///  If receiver is class not derived from TQObject and going to be
///  deleted, disconnect all connections to this receiver.
///  In case of class derived from TQObject it is done automatically.

Bool_t TQObject::Connect(TQObject *sender,
                         const char *signal,
                         const char *cl,
                         void *receiver,
                         const char *slot)
{
   if (cl) {
      TClass *rcv_cl = TClass::GetClass(cl);
      if (rcv_cl) return ConnectToClass(sender, signal, rcv_cl, receiver, slot);
   }

   // the following is the case of receiver class without dictionary
   // e.g. interpreted class or function.

   // sender should be TQObject
   if (!sender->IsA()->InheritsFrom(TQObject::Class()))
      return kFALSE;

   // remove "const" and strip blanks
   TString signal_name = CompressName(signal);
   TString slot_name   = CompressName(slot);

   // check consitency of signal/slot methods/args
   Int_t nsigargs;
   if ((nsigargs = CheckConnectArgs(sender, sender->IsA(), signal_name, 0, slot_name)) == -1)
      return kFALSE;

   if (!sender->fListOfSignals) sender->fListOfSignals = new THashList();

   TQConnectionList *clist = (TQConnectionList*)
      sender->fListOfSignals->FindObject(signal_name);

   if (!clist) {
      clist = new TQConnectionList(signal_name, nsigargs);
      sender->fListOfSignals->Add(clist);
   }

   TQConnection *connection = 0;
   TIter next(clist);

   while ((connection = (TQConnection*)next())) {
      if (!strcmp(slot_name,connection->GetName()) &&
          (receiver == connection->GetReceiver())) break;
   }

   if (!connection)
      connection = new TQConnection(cl, receiver, slot_name);

   // check to prevent multiple entries
   if (!clist->FindObject(connection)) {
      clist->Add(connection);
      if (!connection->FindObject(clist)) connection->Add(clist);
      sender->Connected(signal_name);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This method allows to make a connection from any object
/// of the same class to a single slot.
/// Signal and slot string must have a form:
///    "Draw(char*, Option_t* ,Int_t)"
/// All blanks and "const" words will be removed,
///
/// cl != 0 - class name, it can be class with or
///           without dictionary, e.g interpreted class.
/// Example:
/// ~~~ {.cpp}
///       TGButton *myButton;
///       TH2F     *myHist;
///
///       TQObject::Connect("TGButton", "Clicked()",
///                         "TH2F", myHist, "Draw(Option_t*)");
/// ~~~
/// cl == 0 - corresponds to function (interpereted or global)
///           the name of the function is defined by the slot string,
///           parameter receiver should be 0.
/// Example:
/// ~~~ {.cpp}
///       TGButton *myButton;
///       TH2F     *myHist;
///
///       TQObject::Connect("TGButton", "Clicked()",
///                         0, 0, "hsimple()");
/// ~~~
/// Warning:
///
///  If receiver class not derived from TQObject and going to be
///  deleted, disconnect all connections to this receiver.
///  In case of class derived from TQObject it is done automatically.

Bool_t TQObject::Connect(const char *class_name,
                         const char *signal,
                         const char *cl,
                         void *receiver,
                         const char *slot)
{
   if (cl) {
      TClass *rcv_cl = TClass::GetClass(cl);
      if (rcv_cl) return ConnectToClass(class_name, signal, rcv_cl, receiver,
                                        slot);
   }

   // the following is case of receiver class without dictionary
   // e.g. interpreted class or function.

   TClass *sender = TClass::GetClass(class_name);

   // sender class should be TQObject (i.e. TQClass)
   if (!sender || !sender->IsA()->InheritsFrom(TQObject::Class()))
      return kFALSE;

   TList *slist = ((TQClass*)sender)->fListOfSignals;

   TString signal_name = CompressName(signal);
   TString slot_name   = CompressName(slot);

   // check consitency of signal/slot methods/args
   Int_t nsigargs;
   if ((nsigargs = CheckConnectArgs(0, sender, signal_name, 0, slot_name)) == -1)
      return kFALSE;

   if (!slist) {
      slist = ((TQClass*)sender)->fListOfSignals = new THashList();
   }

   TQConnectionList *clist = (TQConnectionList*)
      slist->FindObject(signal_name);

   if (!clist) {
      clist = new TQConnectionList(signal_name, nsigargs);
      slist->Add(clist);
   }

   TQConnection *connection = 0;
   TIter next(clist);

   while ((connection = (TQConnection*)next())) {
      if (!strcmp(slot_name,connection->GetName()) &&
          (receiver == connection->GetReceiver())) break;
   }

   if (!connection)
      connection = new TQConnection(cl, receiver, slot_name);

   // check to prevent multiple entries
   if (!clist->FindObject(connection)) {
      clist->Add(connection);
      if (!connection->FindObject(clist)) connection->Add(clist);
      ((TQClass*)sender)->Connected(signal_name);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Non-static method is used to connect from the signal
/// of this object to the receiver slot.
///
/// Warning! No check on consistency of sender/receiver
/// classes/methods.
///
/// This method makes possible to have connection/signals from
/// interpreted class. See also RQ_OBJECT.h.

Bool_t TQObject::Connect(const char *signal,
                         const char *receiver_class,
                         void *receiver,
                         const char *slot)
{
   // remove "const" and strip blanks
   TString signal_name = CompressName(signal);
   TString slot_name   = CompressName(slot);

   // check consitency of signal/slot methods/args
   TClass *cl = 0;
   if (receiver_class)
      cl = TClass::GetClass(receiver_class);
   Int_t nsigargs;
   if ((nsigargs = CheckConnectArgs(this, IsA(), signal_name, cl, slot_name)) == -1)
      return kFALSE;

   if (!fListOfSignals) fListOfSignals = new THashList();

   TQConnectionList *clist = (TQConnectionList*)
      fListOfSignals->FindObject(signal_name);

   if (!clist) {
      clist = new TQConnectionList(signal_name, nsigargs);
      fListOfSignals->Add(clist);
   }

   TIter next(clist);
   TQConnection *connection = 0;

   while ((connection = (TQConnection*)next())) {
      if (!strcmp(slot_name,connection->GetName()) &&
          (receiver == connection->GetReceiver())) break;
   }

   if (!connection)
      connection = new TQConnection(receiver_class, receiver, slot_name);

   // check to prevent multiple entries
   if (!clist->FindObject(connection)) {
      clist->Add(connection);
      if (!connection->FindObject(clist)) connection->Add(clist);
      Connected(signal_name);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnects signal in object sender from slot_method in
/// object receiver. For objects derived from TQObject signal-slot
/// connection is removed when either of the objects involved
/// are destroyed.
///
/// Disconnect() is typically used in three ways, as the following
/// examples shows:
///
///  - Disconnect everything connected to an object's signals:
/// ~~~ {.cpp}
///       Disconnect(myObject);
/// ~~~
///  - Disconnect everything connected to a signal:
/// ~~~ {.cpp}
///       Disconnect(myObject, "mySignal()");
/// ~~~
///  - Disconnect a specific receiver:
/// ~~~ {.cpp}
///       Disconnect(myObject, 0, myReceiver, 0);
/// ~~~
///
/// 0 may be used as a wildcard in three of the four arguments,
/// meaning "any signal", "any receiving object" or
/// "any slot in the receiving object", respectively.
///
/// The sender has no default and may never be 0
/// (you cannot disconnect signals from more than one object).
///
/// If signal is 0, it disconnects receiver and slot_method
/// from any signal. If not, only the specified signal is
/// disconnected.
///
/// If  receiver is 0, it disconnects anything connected to signal.
/// If not, slots in objects other than receiver are not
/// disconnected
///
/// If slot_method is 0, it disconnects anything that is connected
/// to receiver.  If not, only slots named slot_method will be
/// disconnected, and all other slots are left alone.
/// The slot_method must be 0 if receiver is left out, so you
/// cannot disconnect a specifically-named slot on all objects.

Bool_t TQObject::Disconnect(TQObject *sender,
                            const char *signal,
                            void *receiver,
                            const char *slot)
{
   Bool_t return_value = kFALSE;
   Bool_t next_return  = kFALSE;

   if (!sender->GetListOfSignals()) return kFALSE;

   TString signal_name = CompressName(signal);
   TString slot_name   = CompressName(slot);

   TQConnectionList *slist = 0;
   TIter next_signal(sender->GetListOfSignals());

   while ((slist = (TQConnectionList*)next_signal()))   {
      if (!signal || signal_name.IsNull()) { // disconnect all signals
         next_return = slist->Disconnect(receiver,slot_name);
         return_value = return_value || next_return;

         if (slist->IsEmpty()) {
            sender->GetListOfSignals()->Remove(slist);
            SafeDelete(slist);            // delete empty list
         }
      } else if (signal && !strcmp(signal_name,slist->GetName())) {
         next_return = slist->Disconnect(receiver,slot_name);
         return_value = return_value || next_return;

         if (slist->IsEmpty()) {
            sender->GetListOfSignals()->Remove(slist);
            SafeDelete(slist);            // delete empty list
            break;
         }
      }
   }

   if (sender->GetListOfSignals() && sender->GetListOfSignals()->IsEmpty()) {
      SafeDelete(sender->fListOfSignals);
   }

   return return_value;
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnects "class signal". The class is defined by class_name.
/// See also Connect(class_name,signal,receiver,slot).

Bool_t TQObject::Disconnect(const char *class_name,
                            const char *signal,
                            void *receiver,
                            const char *slot)
{
   TClass *sender = TClass::GetClass(class_name);

   // sender should be TQClass (which derives from TQObject)
   if (!sender || !sender->IsA()->InheritsFrom(TQObject::Class()))
      return kFALSE;

   TQClass *qcl = (TQClass*)sender;   // cast TClass to TQClass
   return Disconnect(qcl, signal, receiver, slot);
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnects signal of this object from slot of receiver.
/// Equivalent to Disconnect(this, signal, receiver, slot)

Bool_t TQObject::Disconnect(const char *signal,
                            void *receiver,
                            const char *slot)
{
   return Disconnect(this, signal, receiver, slot);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TQObject.

void TQObject::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      // nothing to read
   } else {
      // nothing to write
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if all signals are blocked.

Bool_t TQObject::AreAllSignalsBlocked()
{
   return fgAllSignalsBlocked;
}

////////////////////////////////////////////////////////////////////////////////
/// Block or unblock all signals. Returns the previous block status.

Bool_t TQObject::BlockAllSignals(Bool_t b)
{
   Bool_t ret = fgAllSignalsBlocked;
   fgAllSignalsBlocked = b;
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Global function which simplifies making connection in interpreted ROOT session
///
///  ConnectCINT      - connects to interpreter(CINT) command

Bool_t ConnectCINT(TQObject *sender, const char *signal, const char *slot)
{
   TString str = "ProcessLine(=";
   str += '"';
   str += slot;
   str += '"';
   str += ")";
   return TQObject::Connect(sender, signal, "TInterpreter",
                            gInterpreter, str.Data());
}
