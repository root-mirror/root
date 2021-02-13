// root/core/meta
// vim: sw=3
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_CallFunc
#define ROOT_CallFunc

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingCallFunc                                                       //
//                                                                      //
// Emulation of the CINT CallFunc class.                                //
//                                                                      //
// The CINT C++ interpreter provides an interface for calling           //
// functions through the generated wrappers in dictionaries with        //
// the CallFunc class. This class provides the same functionality,      //
// using an interface as close as possible to CallFunc but the          //
// function metadata and calling service comes from the Cling           //
// C++ interpreter and the Clang C++ compiler, not CINT.                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingMethodInfo.h"
#include "TClingClassInfo.h"
#include "TClingUtils.h"
#include "TInterpreter.h"
#include <string>

#include "cling/Interpreter/Value.h"

#include <llvm/ADT/SmallVector.h>

namespace clang {
class BuiltinType;
class CXXMethodDecl;
class DeclContext;
class Expr;
class FunctionDecl;
}

namespace cling {
class Interpreter;
}

class TClingClassInfo;
class TInterpreterValue;

typedef void (*tcling_callfunc_Wrapper_t)(void*, int, void**, void*);
typedef void (*tcling_callfunc_ctor_Wrapper_t)(void**, void*, unsigned long);
typedef void (*tcling_callfunc_dtor_Wrapper_t)(void*, unsigned long, int);

class TClingCallFunc {

private:

   /// Cling interpreter, we do *not* own.
   cling::Interpreter* fInterp;
   /// ROOT normalized context for that interpreter
   const ROOT::TMetaUtils::TNormalizedCtxt &fNormCtxt;
   /// Current method, we own.
   std::unique_ptr<TClingMethodInfo> fMethod;
   /// Decl for the method
   const clang::FunctionDecl *fDecl = nullptr;
   /// Number of required arguments
   size_t fMinRequiredArguments = -1;
   /// Pointer to compiled wrapper, we do *not* own.
   tcling_callfunc_Wrapper_t fWrapper;
   /// Stored function arguments, we own.
   mutable llvm::SmallVector<cling::Value, 8> fArgVals;
   /// If true, do not limit number of function arguments to declared number.
   bool fIgnoreExtraArgs : 1;
   bool fReturnIsRecordType : 1;

private:
   enum EReferenceType {
      kNotReference,
      kLValueReference,
      kRValueReference
   };

   using ExecWithRetFunc_t =  std::function<void(void* address, cling::Value &ret)>;

   void* compile_wrapper(const std::string& wrapper_name,
                         const std::string& wrapper,
                         bool withAccessControl = true);

   void collect_type_info(clang::QualType& QT, std::ostringstream& typedefbuf,
                          std::ostringstream& callbuf, std::string& type_name,
                          EReferenceType& refType, bool& isPointer, int indent_level,
                          bool forArgument);

   void make_narg_call(const std::string &return_type, const unsigned N, std::ostringstream &typedefbuf,
                       std::ostringstream &callbuf, const std::string &class_name, int indent_level);

   void make_narg_ctor(const unsigned N, std::ostringstream& typedefbuf,
                       std::ostringstream& callbuf,
                       const std::string& class_name, int indent_level);

   void make_narg_call_with_return(const unsigned N,
                                   const std::string& class_name,
                                   std::ostringstream& buf, int indent_level);

   void make_narg_ctor_with_return(const unsigned N,
                                   const std::string& class_name,
                                   std::ostringstream& buf, int indent_level);

   tcling_callfunc_Wrapper_t make_wrapper();

   tcling_callfunc_ctor_Wrapper_t
   make_ctor_wrapper(const TClingClassInfo *, ROOT::TMetaUtils::EIOCtorCategory, const std::string &);

   tcling_callfunc_dtor_Wrapper_t
   make_dtor_wrapper(const TClingClassInfo* info);

   // Implemented in source file.
   template <typename T>
   void execWithLL(void* address, cling::Value* val);
   template <typename T>
   void execWithULL(void* address, cling::Value* val);
   template <class T>
   ExecWithRetFunc_t InitRetAndExecIntegral(clang::QualType QT, cling::Value &ret);

   ExecWithRetFunc_t InitRetAndExecBuiltin(clang::QualType QT, const clang::BuiltinType *BT, cling::Value &ret);
   ExecWithRetFunc_t InitRetAndExecNoCtor(clang::QualType QT, cling::Value &ret);
   ExecWithRetFunc_t InitRetAndExec(const clang::FunctionDecl *FD, cling::Value &ret);

   void exec(void* address, void* ret);

   void exec_with_valref_return(void* address,
                                cling::Value* ret);
   void EvaluateArgList(const std::string& ArgList);

   size_t CalculateMinRequiredArguments();

   size_t GetMinRequiredArguments() {
      if (fMinRequiredArguments == (size_t)-1)
         fMinRequiredArguments = CalculateMinRequiredArguments();
      return fMinRequiredArguments;
   }

   // Implemented in source file.
   template <typename T>
   T ExecT(void* address);


public:

   ~TClingCallFunc() = default;

   explicit TClingCallFunc(cling::Interpreter *interp, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
      : fInterp(interp), fNormCtxt(normCtxt), fWrapper(0), fIgnoreExtraArgs(false), fReturnIsRecordType(false)
   {
      fMethod = std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(interp));
   }

   explicit TClingCallFunc(const TClingMethodInfo &minfo, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
   : fInterp(minfo.GetInterpreter()), fNormCtxt(normCtxt), fWrapper(0), fIgnoreExtraArgs(false),
     fReturnIsRecordType(false)

   {
      fMethod = std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(minfo));
   }

   TClingCallFunc(const TClingCallFunc &rhs)
      : fInterp(rhs.fInterp), fNormCtxt(rhs.fNormCtxt), fWrapper(rhs.fWrapper), fArgVals(rhs.fArgVals),
        fIgnoreExtraArgs(rhs.fIgnoreExtraArgs), fReturnIsRecordType(rhs.fReturnIsRecordType)
   {
      fMethod = std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(*rhs.fMethod));
   }

   TClingCallFunc &operator=(const TClingCallFunc &rhs) = delete;

   void* ExecDefaultConstructor(const TClingClassInfo* info,
                                ROOT::TMetaUtils::EIOCtorCategory kind,
                                const std::string &type_name,
                                void* address = nullptr, unsigned long nary = 0UL);
   void ExecDestructor(const TClingClassInfo* info, void* address = nullptr,
                       unsigned long nary = 0UL, bool withFree = true);
   void ExecWithReturn(void* address, void *ret = nullptr);
   void ExecWithArgsAndReturn(void* address,
                              const void* args[] = 0,
                              int nargs = 0,
                              void* ret = 0);
   void Exec(void* address, TInterpreterValue* interpVal = 0);
   long ExecInt(void* address);
   long long ExecInt64(void* address);
   double ExecDouble(void* address);
   TClingMethodInfo* FactoryMethod() const;
   void IgnoreExtraArgs(bool ignore) { fIgnoreExtraArgs = ignore; }
   void Init();
   void Init(const TClingMethodInfo&);
   void Init(std::unique_ptr<TClingMethodInfo>);
   void Invoke(cling::Value* result = 0) const;
   void* InterfaceMethod();
   bool IsValid() const;
   TInterpreter::CallFuncIFacePtr_t IFacePtr();
   const clang::DeclContext *GetDeclContext() const;

   int get_wrapper_code(std::string &wrapper_name, std::string &wrapper);

   const clang::FunctionDecl *GetDecl() {
      if (!fDecl)
         fDecl = fMethod->GetTargetFunctionDecl();
      return fDecl;
   }
   const clang::FunctionDecl* GetDecl() const {
      if (fDecl)
         return fDecl;
      return fMethod->GetTargetFunctionDecl();
   }
   const clang::Decl *GetFunctionOrShadowDecl() const {
      return fMethod->GetDecl();
   }
   void ResetArg();
   void SetArg(long arg);
   void SetArg(unsigned long arg);
   void SetArg(float arg);
   void SetArg(double arg);
   void SetArg(long long arg);
   void SetArg(unsigned long long arg);
   void SetArgArray(long* argArr, int narg);
   void SetArgs(const char* args);
   void SetFunc(const TClingClassInfo* info, const char* method,
                const char* arglist, long* poffset);
   void SetFunc(const TClingClassInfo* info, const char* method,
                const char* arglist, bool objectIsConst, long* poffset);
   void SetFunc(const TClingMethodInfo* info);
   void SetFuncProto(const TClingClassInfo* info, const char* method,
                     const char* proto, long* poffset,
                     ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   void SetFuncProto(const TClingClassInfo* info, const char* method,
                     const char* proto, bool objectIsConst, long* poffset,
                     ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   void SetFuncProto(const TClingClassInfo* info, const char* method,
                     const llvm::SmallVectorImpl<clang::QualType>& proto,
                     long* poffset,
                     ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   void SetFuncProto(const TClingClassInfo* info, const char* method,
                     const llvm::SmallVectorImpl<clang::QualType>& proto,
                     bool objectIsConst, long* poffset,
                     ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
};

#endif // ROOT_CallFunc
