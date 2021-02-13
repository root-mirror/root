// @(#)root/meta:$Id$
// Author: Fons Rademakers   13/06/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMethodCall
#define ROOT_TMethodCall


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMethodCall                                                          //
//                                                                      //
// Method or function calling interface. Objects of this class contain  //
// the (CINT) environment to call a global function or a method for an  //
// object of a specific class with the desired arguments. This class is //
// espicially useful when a method has to be called more times for      //
// different objects and/or with different arguments. If a function or  //
// method needs to be called only once one better uses                  //
// TInterpreter::Execute().                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include "TInterpreter.h"

class TClass;
class TFunction;

class TMethodCall : public TObject {

public:
   using EReturnType = TInterpreter::EReturnType;

   // For backward compatibility:
   static const EReturnType kLong = TInterpreter::EReturnType::kLong;
   static const EReturnType kDouble = TInterpreter::EReturnType::kDouble;
   static const EReturnType kString = TInterpreter::EReturnType::kString;
   static const EReturnType kOther = TInterpreter::EReturnType::kOther;
   static const EReturnType kNoReturnType = TInterpreter::EReturnType::kNoReturnType;
   // Historical name.
   static const EReturnType kNone = TInterpreter::EReturnType::kNoReturnType;

   // enum EReturnType { kLong, kDouble, kString, kOther, kNone };

private:
   CallFunc_t    *fFunc;      //CINT method invocation environment
   Long_t         fOffset;    //offset added to object pointer before method invocation
   TClass        *fClass;     //pointer to the class info
   TFunction     *fMetPtr;    //pointer to the method or function info
   TString        fMethod;    //method name
   TString        fParams;    //argument string
   TString        fProto;     //prototype string
   Bool_t         fDtorOnly;  //call only dtor and not delete when calling ~xxx
   EReturnType    fRetType;   //method return type

   void Execute(const char *,  const char *, int * /*error*/ = 0) { }    // versions of TObject
   void Execute(TMethod *, TObjArray *, int * /*error*/ = 0) { }

   void InitImplementation(const char *methodname, const char *params, const char *proto, Bool_t objectIsConst, TClass *cl, const ClassInfo_t *cinfo, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);

public:
   TMethodCall();
   TMethodCall(TClass *cl, CallFunc_t *callfunc, Long_t offset = 0);
   TMethodCall(TClass *cl, const char *method, const char *params);
   TMethodCall(const char *function, const char *params);
   TMethodCall(const TFunction *func);
   TMethodCall(const TMethodCall &org);
   TMethodCall& operator=(const TMethodCall &rhs);
   ~TMethodCall();

   void           Init(const TFunction *func);
   void           Init(TClass *cl, CallFunc_t *func, Long_t offset = 0);
   void           Init(TClass *cl, const char *method, const char *params, Bool_t objectIsConst = kFALSE);
   void           Init(const char *function, const char *params);
   void           InitWithPrototype(TClass *cl, const char *method, const char *proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   void           InitWithPrototype(const char *function, const char *proto, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   Bool_t         IsValid() const;
   TObject       *Clone(const char *newname="") const;
   void           CallDtorOnly(Bool_t set = kTRUE) { fDtorOnly = set; }

   TFunction     *GetMethod();
   const char    *GetMethodName() const { return fMethod.Data(); }
   const char    *GetParams() const { return fParams.Data(); }
   const char    *GetProto() const { return fProto.Data(); }
   CallFunc_t    *GetCallFunc() const { return fFunc; }
   EReturnType    ReturnType();

   void     SetParamPtrs(void *paramArr, Int_t nparam = -1);
   void     ResetParam();
   void     SetParam(Long_t l);
   void     SetParam(Float_t f);
   void     SetParam(Double_t d);
   void     SetParam(Long64_t ll);
   void     SetParam(ULong64_t ull);

   template <typename... T> void SetParams(const T&... params) {
      if (!fFunc) return;
      gInterpreter->CallFunc_SetArguments(fFunc,params...);
   }

   void     Execute(void *object);
   void     Execute(void *object, const char *params);
   void     Execute(void *object, Long_t &retLong);
   void     Execute(void *object, const char *params, Long_t &retLong);
   void     Execute(void *object, Double_t &retDouble);
   void     Execute(void *object, const char *params, Double_t &retDouble);

   void     Execute(void *object, char **retText);
   void     Execute(void *object, const char *params, char **retText);

   void     Execute();
   void     Execute(const char *params);
   void     Execute(Long_t &retLong);
   void     Execute(const char *params, Long_t &retLong);
   void     Execute(Double_t &retDouble);
   void     Execute(const char *params, Double_t &retDouble);

   void     Execute(void *objAddress, const void* args[], int nargs, void *ret = 0);

   ClassDef(TMethodCall,0)  //Method calling interface
};

inline void TMethodCall::Execute()
   { Execute((void *)0); }
inline void TMethodCall::Execute(const char *params)
   { Execute((void *)0, params); }
inline void TMethodCall::Execute(Long_t &retLong)
   { Execute((void *)0, retLong); }
inline void TMethodCall::Execute(const char *params, Long_t &retLong)
   { Execute((void *)0, params, retLong); }
inline void TMethodCall::Execute(Double_t &retDouble)
   { Execute((void *)0, retDouble); }
inline void TMethodCall::Execute(const char *params, Double_t &retDouble)
   { Execute((void *)0, params, retDouble); }

#endif
