// @(#)root/base:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   15/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TQConnection
\ingroup Base

TQConnection class is an internal class, used in the object
communication mechanism.

TQConnection:
   -  is a list of signal_lists containing pointers
      to this connection
   -  receiver is the object to which slot-method is applied
*/

#include "TQConnection.h"
#include "TROOT.h"
#include "TRefCnt.h"
#include "TClass.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TDataType.h"
#include "TInterpreter.h"
#include <iostream>
#include "TVirtualMutex.h"
#include "THashTable.h"

ClassImpQ(TQConnection)

char *gTQSlotParams; // used to pass string parameter

/** \class TQSlot
Slightly modified TMethodCall class used in the object communication mechanism.
*/

class TQSlot : public TObject, public TRefCnt {

protected:
   CallFunc_t    *fFunc;      // CINT method invocation environment
   ClassInfo_t   *fClass;     // CINT class for fFunc
   TFunction     *fMethod;    // slot method or global function
   Long_t         fOffset;    // offset added to object pointer
   TString        fName;      // full name of method
   Int_t          fExecuting; // true if one of this slot's ExecuteMethod methods is being called
public:
   TQSlot(TClass *cl, const char *method, const char *funcname);
   TQSlot(const char *class_name, const char *funcname);
   virtual ~TQSlot();

   Bool_t      CheckSlot(Int_t nargs) const;
   Long_t      GetOffset() const { return fOffset; }
   CallFunc_t *StartExecuting();
   CallFunc_t *GetFunc() const { return fFunc; }
   void        EndExecuting();

   const char *GetName() const {
      return fName.Data();
   }

   Int_t GetMethodNargs() { return fMethod->GetNargs(); }

   void ExecuteMethod(void *object, Int_t nargs, va_list ap) = delete;
   void ExecuteMethod(void *object);
   void ExecuteMethod(void *object, Long_t param);
   void ExecuteMethod(void *object, Long64_t param);
   void ExecuteMethod(void *object, Double_t param);
   void ExecuteMethod(void *object, const char *params);
   void ExecuteMethod(void *object, Long_t *paramArr, Int_t nparam = -1);
   void Print(Option_t *opt = "") const;
   void ls(Option_t *opt = "") const {
      Print(opt);
   }

   Bool_t IsExecuting() const {
      return fExecuting > 0;
   }
};


////////////////////////////////////////////////////////////////////////////////
/// Create the method invocation environment. Necessary input
/// information: the class, full method name with prototype
/// string of the form: method(char*,int,float).
/// To initialize class method with default arguments, method
/// string with default parameters should be of the form:
///
/// method(=\"ABC\",1234,3.14) (!! parameter string should
/// consists of '=').
///
/// To execute the method call TQSlot::ExecuteMethod(object,...).

TQSlot::TQSlot(TClass *cl, const char *method_name,
               const char *funcname) : TObject(), TRefCnt()
{
   fFunc      = 0;
   fClass     = 0;
   fOffset    = 0;
   fMethod    = 0;
   fName      = "";
   fExecuting = 0;

   // cl==0, is the case of interpreted function.

   fName = method_name;

   char *method = new char[strlen(method_name) + 1];
   if (method) strcpy(method, method_name);

   char *proto;
   char *tmp;
   char *params = 0;

   // separate method and prototype strings

   if ((proto = strchr(method, '('))) {

      // substitute first '(' symbol with '\0'
      *proto++ = '\0';

      // last ')' symbol with '\0'
      if ((tmp = strrchr(proto, ')'))) * tmp = '\0';
      if ((params = strchr(proto, '='))) * params = ' ';
   }

   R__LOCKGUARD(gInterpreterMutex);
   fFunc = gCling->CallFunc_Factory();

   // initiate class method (function) with proto
   // or with default params

   if (cl) {
      if (params) {
         gCling->CallFunc_SetFunc(fFunc, cl->GetClassInfo(), method, params, &fOffset);
         fMethod = cl->GetMethod(method, params);
      } else {
         gCling->CallFunc_SetFuncProto(fFunc, cl->GetClassInfo(), method, proto, &fOffset);
         fMethod = cl->GetMethodWithPrototype(method, proto);
      }
   } else {
      fClass = gCling->ClassInfo_Factory();
      if (params) {
         gCling->CallFunc_SetFunc(fFunc, fClass, (char *)funcname, params, &fOffset);
         fMethod = gROOT->GetGlobalFunction(funcname, params, kFALSE);
      } else {
         gCling->CallFunc_SetFuncProto(fFunc, fClass, (char *)funcname, proto, &fOffset);
         fMethod = gROOT->GetGlobalFunctionWithPrototype(funcname, proto, kFALSE);
      }
   }

   // cleaning
   delete [] method;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the method invocation environment. Necessary input
/// information: the name of class (could be interpreted class),
/// full method name with prototype or parameter string
/// of the form: method(char*,int,float).
/// To initialize class method with default arguments, method
/// string with default parameters  should be of the form:
///
/// method(=\"ABC\",1234,3.14) (!! parameter string should
/// consists of '=').
///
/// To execute the method call TQSlot::ExecuteMethod(object,...).

TQSlot::TQSlot(const char *class_name, const char *funcname) :
   TObject(), TRefCnt()
{
   fFunc      = 0;
   fClass     = 0;
   fOffset    = 0;
   fMethod    = 0;
   fName      = funcname;
   fExecuting = 0;

   char *method = new char[strlen(funcname) + 1];
   if (method) strcpy(method, funcname);

   char *proto;
   char *tmp;
   char *params = 0;

   // separate method and prototype strings

   if ((proto =  strchr(method, '('))) {
      *proto++ = '\0';
      if ((tmp = strrchr(proto, ')'))) * tmp  = '\0';
      if ((params = strchr(proto, '='))) * params = ' ';
   }

   R__LOCKGUARD(gInterpreterMutex);
   fFunc = gCling->CallFunc_Factory();
   gCling->CallFunc_IgnoreExtraArgs(fFunc, true);

   fClass = gCling->ClassInfo_Factory();
   TClass *cl = 0;

   if (!class_name)
      ;                       // function
   else {
      gCling->ClassInfo_Init(fClass, class_name);  // class
      cl = TClass::GetClass(class_name);
   }

   if (params) {
      gCling->CallFunc_SetFunc(fFunc, fClass, method, params, &fOffset);
      if (cl)
         fMethod = cl->GetMethod(method, params);
      else
         fMethod = gROOT->GetGlobalFunction(method, params, kTRUE);
   } else {
      gCling->CallFunc_SetFuncProto(fFunc, fClass, method, proto , &fOffset);
      if (cl)
         fMethod = cl->GetMethodWithPrototype(method, proto);
      else
         fMethod = gROOT->GetGlobalFunctionWithPrototype(method, proto, kTRUE);
   }

   delete [] method;
}

////////////////////////////////////////////////////////////////////////////////
/// TQSlot dtor.

TQSlot::~TQSlot()
{
   // don't delete executing environment of a slot that is being executed
   if (!fExecuting) {
      gCling->CallFunc_Delete(fFunc);
      gCling->ClassInfo_Delete(fClass);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ExecuteMethod the method (with preset arguments) for
/// the specified object.

inline void TQSlot::ExecuteMethod(void *object)
{
   ExecuteMethod(object, (Long_t*)nullptr, 0);

}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the method is valid and the number of arguments is
/// acceptable.

inline Bool_t TQSlot::CheckSlot(Int_t nargs) const
{
   if (!fMethod) {
      Error("ExecuteMethod", "method %s not found,"
            "\n(note: interpreted methods are not supported with varargs)",
            fName.Data());
      return kFALSE;
   }

   if (nargs < fMethod->GetNargs() - fMethod->GetNargsOpt() ||
       nargs > fMethod->GetNargs()) {
      Error("ExecuteMethod", "nargs (%d) not consistent with expected number of arguments ([%d-%d])",
            nargs, fMethod->GetNargs() - fMethod->GetNargsOpt(),
            fMethod->GetNargs());
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Mark the slot as executing.

CallFunc_t *TQSlot::StartExecuting() {
   fExecuting++;
   return fFunc;
}

////////////////////////////////////////////////////////////////////////////////
/// Mark the slot as no longer executing and cleanup if need be.

void TQSlot::EndExecuting() {
   fExecuting--;
   if (!TestBit(kNotDeleted) && !fExecuting)
      gCling->CallFunc_Delete(fFunc);
}

////////////////////////////////////////////////////////////////////////////////
/// ExecuteMethod the method for the specified object and
/// with single argument value.

inline void TQSlot::ExecuteMethod(void *object, Long_t param)
{
   ExecuteMethod(object, &param, 1);

}

////////////////////////////////////////////////////////////////////////////////
/// ExecuteMethod the method for the specified object and
/// with single argument value.

inline void TQSlot::ExecuteMethod(void *object, Long64_t param)
{
   Long_t *arg = reinterpret_cast<Long_t *>(&param);
   ExecuteMethod(object, arg, 1);

}

////////////////////////////////////////////////////////////////////////////////
/// ExecuteMethod the method for the specified object and
/// with single argument value.

inline void TQSlot::ExecuteMethod(void *object, Double_t param)
{
   Long_t *arg = reinterpret_cast<Long_t *>(&param);
   ExecuteMethod(object, arg, 1);

}

////////////////////////////////////////////////////////////////////////////////
/// ExecuteMethod the method for the specified object and text param.

inline void TQSlot::ExecuteMethod(void *object, const char *param)
{
   Long_t arg = reinterpret_cast<Long_t>(param);
   ExecuteMethod(object, &arg, 1);

}

////////////////////////////////////////////////////////////////////////////////
/// ExecuteMethod the method for the specified object and with
/// several argument values.
/// ParamArr is an array containing the function argument values.
/// If nparam = -1 then paramArr must contain values for all function
/// arguments, otherwise Nargs-NargsOpt <= nparam <= Nargs, where
/// Nargs is the number of all arguments and NargsOpt is the number
/// of default arguments.

inline void TQSlot::ExecuteMethod(void *object, Long_t *paramArr, Int_t nparam)
{
   void *address = 0;
   R__LOCKGUARD(gInterpreterMutex);
   if (paramArr) gCling->CallFunc_SetArgArray(fFunc, paramArr, nparam);
   if (object) address = (void *)((Long_t)object + fOffset);
   fExecuting++;
   gCling->CallFunc_Exec(fFunc, address);
   fExecuting--;
   if (!TestBit(kNotDeleted) && !fExecuting)
      gCling->CallFunc_Delete(fFunc);
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about slot.

void TQSlot::Print(Option_t *) const
{
   std::cout << IsA()->GetName() << "\t" << GetName() << "\t"
             << "Number of Connections = " << References() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

class TQSlotPool {
private:
   THashTable *fTable;
public:
   TQSlotPool() {
      fTable = new THashTable(50);
   }
   virtual ~TQSlotPool() {
      fTable->Clear("nodelete");
   }

   TQSlot  *New(const char *class_name, const char *funcname);
   TQSlot  *New(TClass *cl, const char *method, const char *func);
   void     Free(TQSlot *slot);
};

////////////////////////////////////////////////////////////////////////////////
/// Create new slot or return already existing one.

TQSlot *TQSlotPool::New(const char *class_name, const char *funcname)
{
   TString name = class_name;
   name += "::";
   name += funcname;

   TQSlot *slot = (TQSlot *)fTable->FindObject(name.Data());

   if (!slot) {
      slot = new TQSlot(class_name, funcname);
      fTable->Add(slot);
   }
   slot->AddReference();
   return slot;
}

////////////////////////////////////////////////////////////////////////////////
/// Create new slot or return already existing one.

TQSlot *TQSlotPool::New(TClass *cl, const char *method, const char *func)
{
   TString name;

   if (cl) {
      name = cl->GetName();
      name += "::";
      name += method;
   } else {
      name = "::";
      name += func;
   }

   TQSlot *slot = (TQSlot *)fTable->FindObject(name.Data());

   if (!slot) {
      slot = new TQSlot(cl, method, func);
      fTable->Add(slot);
   }
   slot->AddReference();
   return slot;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete slot if there is no reference to it.

void TQSlotPool::Free(TQSlot *slot)
{
   slot->RemoveReference();  // decrease references to slot

   if (slot->References() <= 0) {
      fTable->Remove(slot);
      if (!slot->IsExecuting()) SafeDelete(slot);
   }
}

static TQSlotPool gSlotPool;  // global pool of slots

void TQConnection::SetArg(const Long_t *params, Int_t nparam/* = -1*/) {
   if (nparam == -1)
      nparam = fSlot->GetMethodNargs();

   // FIXME: Why TInterpreter needs non-const SetArgArray. TClingCallFunc
   // doesn't modify the value.
   gInterpreter->CallFunc_SetArgArray(fSlot->GetFunc(), const_cast<Long_t*>(params), nparam);
}


////////////////////////////////////////////////////////////////////////////////
/// TQConnection ctor.
///    cl != 0  - connection to object == receiver of class == cl
///               and method == method_name
///    cl == 0  - connection to function with name == method_name

TQConnection::TQConnection(TClass *cl, void *receiver, const char *method_name)
   : TQObject()
{
   const char *funcname = 0;
   fReceiver = receiver;      // fReceiver is pointer to receiver

   if (!cl) {
      Error("SetFCN", "Not used anymore.");
      /*
      funcname = gCling->Getp2f2funcname(fReceiver);
      if (!funcname)
         Warning("TQConnection", "%s cannot be compiled", method_name);
      */
   }

   if (cl) fClassName = cl->GetName();
   fSlot = gSlotPool.New(cl, method_name, funcname);
}

////////////////////////////////////////////////////////////////////////////////
/// TQConnection ctor.
///    Creates connection to method of class specified by name,
///    it could be interpreted class and with method == funcname.

TQConnection::TQConnection(const char *class_name, void *receiver,
                           const char *funcname) : TQObject()
{
   fClassName = class_name;
   fSlot = gSlotPool.New(class_name, funcname);  // new slot-method
   fReceiver = receiver;      // fReceiver is pointer to receiver
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Ignore connections to this TQConnections

TQConnection::TQConnection(const TQConnection &con) : TQObject()
{
   fClassName = con.fClassName;
   fSlot = con.fSlot;
   fSlot->AddReference();
   fReceiver = con.fReceiver;
}

////////////////////////////////////////////////////////////////////////////////
/// TQConnection dtor.
///    - remove this connection from all signal lists
///    - we do not delete fSlot if it has other connections,
///      TQSlot::fCounter > 0 .

TQConnection::~TQConnection()
{
   TIter next(this);
   TList *list;

   while ((list = (TList *)next())) {
      list->Remove(this);
      if (list->IsEmpty()) delete list;   // delete empty list
   }
   Clear("nodelete");

   if (!fSlot) return;
   gSlotPool.Free(fSlot);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of connection (aka name of slot)

const char *TQConnection::GetName() const
{
   return fSlot->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Signal Destroyed tells that connection is destroyed.

void TQConnection::Destroyed()
{
   MakeZombie();
   Emit("Destroyed()");
}

////////////////////////////////////////////////////////////////////////////////
/// List TQConnection full method name and list all signals
/// connected to this connection.

void TQConnection::ls(Option_t *option) const
{
   std::cout << "\t" <<  IsA()->GetName() << "\t" << GetName() << std::endl;
   ((TQConnection *)this)->R__FOR_EACH(TList, ls)(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Print TQConnection full method name and print all
/// signals connected to this connection.

void TQConnection::PrintCollectionHeader(Option_t *) const
{
   TROOT::IndentLevel();
   std::cout << IsA()->GetName() << "\t" << fReceiver << "\t" << GetName() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply slot-method to the fReceiver object without arguments.

void TQConnection::ExecuteMethod()
{
   // This connection might be deleted in result of the method execution
   // (for example in case of a Disconnect).  Hence we do not assume
   // the object is still valid on return.
   TQSlot *s = fSlot;
   fSlot->ExecuteMethod(fReceiver);
   if (s->References() <= 0) delete s;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply slot-method to the fReceiver object with
/// single argument value.

void TQConnection::ExecuteMethod(Long_t param)
{
   // This connection might be deleted in result of the method execution
   // (for example in case of a Disconnect).  Hence we do not assume
   // the object is still valid on return.
   TQSlot *s = fSlot;
   fSlot->ExecuteMethod(fReceiver, param);
   if (s->References() <= 0) delete s;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply slot-method to the fReceiver object with
/// single argument value.

void TQConnection::ExecuteMethod(Long64_t param)
{
   // This connection might be deleted in result of the method execution
   // (for example in case of a Disconnect).  Hence we do not assume
   // the object is still valid on return.
   TQSlot *s = fSlot;
   fSlot->ExecuteMethod(fReceiver, param);
   if (s->References() <= 0) delete s;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply slot-method to the fReceiver object with
/// single argument value.

void TQConnection::ExecuteMethod(Double_t param)
{
   // This connection might be deleted in result of the method execution
   // (for example in case of a Disconnect).  Hence we do not assume
   // the object is still valid on return.
   TQSlot *s = fSlot;
   fSlot->ExecuteMethod(fReceiver, param);
   if (s->References() <= 0) delete s;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply slot-method to the fReceiver object with variable
/// number of argument values.

void TQConnection::ExecuteMethod(Long_t *params, Int_t nparam)
{
   // This connection might be deleted in result of the method execution
   // (for example in case of a Disconnect).  Hence we do not assume
   // the object is still valid on return.
   TQSlot *s = fSlot;
   fSlot->ExecuteMethod(fReceiver, params, nparam);
   if (s->References() <= 0) delete s;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply slot-method to the fReceiver object and
/// with string parameter.

void TQConnection::ExecuteMethod(const char *param)
{
   // This connection might be deleted in result of the method execution
   // (for example in case of a Disconnect).  Hence we do not assume
   // the object is still valid on return.
   TQSlot *s = fSlot;
   fSlot->ExecuteMethod(fReceiver, param);
   if (s->References() <= 0) delete s;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the underlying method is value and the number of argument
/// is compatible.

Bool_t TQConnection::CheckSlot(Int_t nargs) const {
   return fSlot->CheckSlot(nargs);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object address to be passed to the function.

void *TQConnection::GetSlotAddress() const {
   if (fReceiver) return (void *)((Long_t)fReceiver + fSlot->GetOffset());
   else return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Lock the interpreter and mark the slot as executing.

CallFunc_t *TQConnection::LockSlot() const {
   if (gInterpreterMutex) gInterpreterMutex->Lock();
   return fSlot->StartExecuting();
}

////////////////////////////////////////////////////////////////////////////////
/// Unlock the interpreter and mark the slot as no longer executing.

void TQConnection::UnLockSlot(TQSlot *s) const {
   s->EndExecuting();
   if (s->References() <= 0) delete s;
   if (gInterpreterMutex) gInterpreterMutex->UnLock();
}

CallFunc_t* TQConnection::GetSlotCallFunc() const {
   return fSlot->GetFunc();
}
