// @(#)root/metautils:$Id$
// Author: Axel Naumann, Nov 2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMetaUtils
#define ROOT_TMetaUtils

#include "RConversionRuleParser.h"

#include <functional>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <list>
#include <map>
#include <utility>

#include <cstdlib>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#include "clang/Basic/Module.h"

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace llvm {
   class StringRef;
}

namespace clang {
   class ASTContext;
   class Attr;
   class ClassTemplateDecl;
   class ClassTemplateSpecializationDecl;
   class CompilerInstance;
   class CXXBaseSpecifier;
   class CXXRecordDecl;
   class Decl;
   class DeclContext;
   class DeclaratorDecl;
   class FieldDecl;
   class FunctionDecl;
   class NamedDecl;
   class ParmVarDecl;
   class PresumedLoc;
   class QualType;
   class RecordDecl;
   class SourceLocation;
   class TagDecl;
   class TemplateDecl;
   class TemplateName;
   class TemplateArgument;
   class TemplateArgumentList;
   class TemplateParameterList;
   class Type;
   class TypeDecl;
   class TypedefNameDecl;
   struct PrintingPolicy;
}

namespace cling {
   class Interpreter;
   class LookupHelper;
   namespace utils {
      namespace Transform {
         struct Config;
      }
   }
}

// For ROOT::ESTLType
#include "ESTLType.h"

// for TClassEdit::TInterpreterLookupHelper
#include "TClassEdit.h"

#include "Varargs.h"

namespace ROOT {
namespace TMetaUtils {

///\returns the resolved normalized absolute path possibly resolving symlinks.
std::string GetRealPath(const std::string &path);

// Forward Declarations --------------------------------------------------------
class AnnotatedRecordDecl;

// Constants, typedefs and Enums -----------------------------------------------

// Convention for the ROOT relevant properties
namespace propNames{
   static const std::string separator("@@@");
   static const std::string iotype("iotype");
   static const std::string name("name");
   static const std::string pattern("pattern");
   static const std::string ioname("ioname");
   static const std::string comment("comment");
   static const std::string nArgsToKeep("nArgsToKeep");
   static const std::string persistent("persistent");
   static const std::string transient("transient");
}

// Get the array index information for a data member.
enum DataMemberInfo__ValidArrayIndex_error_code { VALID, NOT_INT, NOT_DEF, IS_PRIVATE, UNKNOWN };

typedef void (*CallWriteStreamer_t)(const AnnotatedRecordDecl &cl,
                                    const cling::Interpreter &interp,
                                    const TNormalizedCtxt &normCtxt,
                                    std::ostream& dictStream,
                                    bool isAutoStreamer);

const int kInfo            =      0;
const int kNote            =    500;
const int kWarning         =   1000;
const int kError           =   2000;
const int kSysError        =   3000;
const int kFatal           =   4000;
const int kMaxLen          =   1024;

// Classes ---------------------------------------------------------------------
class TNormalizedCtxtImpl;

//______________________________________________________________________________
class TNormalizedCtxt {
private:
   TNormalizedCtxtImpl* fImpl;
public:
   using Config_t = cling::utils::Transform::Config;
   using TypesCont_t = std::set<const clang::Type*>;
   using TemplPtrIntMap_t = std::map<const clang::ClassTemplateDecl*, int>;

   TNormalizedCtxt(const cling::LookupHelper &lh);
   TNormalizedCtxt(const TNormalizedCtxt& other);
   ~TNormalizedCtxt();
   const Config_t& GetConfig() const;
   const TypesCont_t &GetTypeWithAlternative() const;

   void AddTemplAndNargsToKeep(const clang::ClassTemplateDecl* templ, unsigned int i);
   int GetNargsToKeep(const clang::ClassTemplateDecl* templ) const;
   const TemplPtrIntMap_t GetTemplNargsToKeepMap() const;
   void keepTypedef(const cling::LookupHelper &lh, const char* name,
                    bool replace = false);
};

//______________________________________________________________________________
class TClingLookupHelper : public TClassEdit::TInterpreterLookupHelper {
public:
   typedef bool (*ExistingTypeCheck_t)(const std::string &tname, std::string &result);
   typedef bool (*AutoParse_t)(const char *name);

private:
   cling::Interpreter *fInterpreter;
   TNormalizedCtxt    *fNormalizedCtxt;
   ExistingTypeCheck_t fExistingTypeCheck;
   AutoParse_t         fAutoParse;
   bool               *fInterpreterIsShuttingDownPtr;
   const int          *fPDebug; // debug flag, might change at runtime thus *
   bool WantDiags() const { return fPDebug && *fPDebug > 5; }

public:
   TClingLookupHelper(cling::Interpreter &interpreter, TNormalizedCtxt &normCtxt,
                      ExistingTypeCheck_t existingTypeCheck,
                      AutoParse_t autoParse,
                      bool *shuttingDownPtr,
                      const int *pgDebug = 0);
   virtual ~TClingLookupHelper() { /* we're not owner */ }

   virtual bool ExistingTypeCheck(const std::string &tname, std::string &result);
   virtual void GetPartiallyDesugaredName(std::string &nameLong);
   virtual bool IsAlreadyPartiallyDesugaredName(const std::string &nondef, const std::string &nameLong);
   virtual bool IsDeclaredScope(const std::string &base, bool &isInlined);
   virtual bool GetPartiallyDesugaredNameWithScopeHandling(const std::string &tname, std::string &result, bool dropstd = true);
   virtual void ShuttingDownSignal();
};

//______________________________________________________________________________
class AnnotatedRecordDecl {
private:
   static std::string BuildDemangledTypeInfo(const clang::RecordDecl *rDecl,
                                             const std::string &normalizedName);
   long fRuleIndex;
   const clang::RecordDecl* fDecl;
   std::string fRequestedName;
   std::string fNormalizedName;
   std::string fDemangledTypeInfo;
   bool fRequestStreamerInfo;
   bool fRequestNoStreamer;
   bool fRequestNoInputOperator;
   bool fRequestOnlyTClass;
   int  fRequestedVersionNumber;

public:
   enum ERootFlag {
      kNoStreamer      = 0x01,
      kNoInputOperator = 0x02,
      kUseByteCount    = 0x04,
      kStreamerInfo    = 0x04,
      kHasVersion      = 0x08
   };

   AnnotatedRecordDecl(long index,
                       const clang::RecordDecl *decl,
                       bool rStreamerInfo,
                       bool rNoStreamer,
                       bool rRequestNoInputOperator,
                       bool rRequestOnlyTClass,
                       int rRequestedVersionNumber,
                       const cling::Interpreter &interpret,
                       const TNormalizedCtxt &normCtxt);

   AnnotatedRecordDecl(long index,
                       const clang::RecordDecl *decl,
                       const char *requestName,
                       bool rStreamerInfo,
                       bool rNoStreamer,
                       bool rRequestNoInputOperator,
                       bool rRequestOnlyTClass,
                       int rRequestedVersionNumber,
                       const cling::Interpreter &interpret,
                       const TNormalizedCtxt &normCtxt);

   AnnotatedRecordDecl(long index,
                       const clang::Type *requestedType,
                       const clang::RecordDecl *decl,
                       const char *requestedName,
                       bool rStreamerInfo,
                       bool rNoStreamer,
                       bool rRequestNoInputOperator,
                       bool rRequestOnlyTClass,
                       int rRequestedVersionNumber,
                       const cling::Interpreter &interpret,
                       const TNormalizedCtxt &normCtxt);

   AnnotatedRecordDecl(long index,
                       const clang::Type *requestedType,
                       const clang::RecordDecl *decl,
                       const char *requestedName,
                       unsigned int nTemplateArgsToSkip,
                       bool rStreamerInfo,
                       bool rNoStreamer,
                       bool rRequestNoInputOperator,
                       bool rRequestOnlyTClass,
                       int rRequestedVersionNumber,
                       const cling::Interpreter &interpret,
                       const TNormalizedCtxt &normCtxt);

   ~AnnotatedRecordDecl() {
      // Nothing to do we do not own the pointer;
   }


   long GetRuleIndex() const { return fRuleIndex; }

   const char *GetRequestedName() const { return fRequestedName.c_str(); }
   const char *GetNormalizedName() const { return fNormalizedName.c_str(); }
   const std::string &GetDemangledTypeInfo() const { return fDemangledTypeInfo; }
   bool HasClassVersion() const { return fRequestedVersionNumber >=0 ; }
   bool RequestStreamerInfo() const {
      // Equivalent to CINT's cl.RootFlag() & G__USEBYTECOUNT
      return fRequestStreamerInfo;
   }
   bool RequestNoInputOperator() const { return fRequestNoInputOperator; }
   bool RequestNoStreamer() const { return fRequestNoStreamer; }
   bool RequestOnlyTClass() const { return fRequestOnlyTClass; }
   int  RequestedVersionNumber() const { return fRequestedVersionNumber; }
   int  RootFlag() const {
      // Return the request (streamerInfo, has_version, etc.) combined in a single
      // int.  See RScanner::AnnotatedRecordDecl::ERootFlag.
      int result = 0;
      if (fRequestNoStreamer) result = kNoStreamer;
      if (fRequestNoInputOperator) result |= kNoInputOperator;
      if (fRequestStreamerInfo) result |= kStreamerInfo;
      if (fRequestedVersionNumber > -1) result |= kHasVersion;
      return result;
   }
   const clang::RecordDecl* GetRecordDecl() const { return fDecl; }

   operator clang::RecordDecl const *() const {
      return fDecl;
   }

   bool operator<(const AnnotatedRecordDecl& right) const
   {
      return fRuleIndex < right.fRuleIndex;
   }

   struct CompareByName {
      bool operator() (const AnnotatedRecordDecl& right, const AnnotatedRecordDecl& left) const
      {
         return left.fNormalizedName < right.fNormalizedName;
      }
   };
};

//______________________________________________________________________________
class RConstructorType {
private:
   const std::string           fArgTypeName;
   const clang::CXXRecordDecl *fArgType;

public:
   RConstructorType(const char *type_of_arg, const cling::Interpreter&);

   const char *GetName() const ;
   const clang::CXXRecordDecl *GetType() const;
};
typedef std::list<RConstructorType> RConstructorTypes;

// Functions -------------------------------------------------------------------

//______________________________________________________________________________
int extractAttrString(clang::Attr* attribute, std::string& attrString);

//______________________________________________________________________________
int extractPropertyNameValFromString(const std::string attributeStr,std::string& attrName, std::string& attrValue);

//______________________________________________________________________________
int extractPropertyNameVal(clang::Attr* attribute, std::string& attrName, std::string& attrValue);

//______________________________________________________________________________
bool ExtractAttrPropertyFromName(const clang::Decl& decl,
                                 const std::string& propName,
                                 std::string& propValue);

//______________________________________________________________________________
bool ExtractAttrIntPropertyFromName(const clang::Decl& decl,
                                    const std::string& propName,
                                    int& propValue);

//______________________________________________________________________________
bool RequireCompleteType(const cling::Interpreter &interp, const clang::CXXRecordDecl *cl);

//______________________________________________________________________________
bool RequireCompleteType(const cling::Interpreter &interp, clang::SourceLocation Loc, clang::QualType Type);

//______________________________________________________________________________
// Add default template parameters.
clang::QualType AddDefaultParameters(clang::QualType instanceType,
                                     const cling::Interpreter &interpret,
                                     const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
llvm::StringRef DataMemberInfo__ValidArrayIndex(const clang::DeclaratorDecl &m, int *errnum = 0, llvm::StringRef  *errstr = 0);

enum class EIOCtorCategory : short { kAbsent, kDefault, kIOPtrType, kIORefType };

//______________________________________________________________________________
EIOCtorCategory CheckConstructor(const clang::CXXRecordDecl*, const RConstructorType&, const cling::Interpreter& interp);

//______________________________________________________________________________
bool CheckDefaultConstructor(const clang::CXXRecordDecl*, const cling::Interpreter& interp);

//______________________________________________________________________________
EIOCtorCategory CheckIOConstructor(const clang::CXXRecordDecl*, const char *, const clang::CXXRecordDecl *, const cling::Interpreter& interp);

//______________________________________________________________________________
const clang::FunctionDecl* ClassInfo__HasMethod(const clang::DeclContext *cl, char const*, const cling::Interpreter& interp);

//______________________________________________________________________________
void CreateNameTypeMap(clang::CXXRecordDecl const&, std::map<std::string, ROOT::Internal::TSchemaType>&);

//______________________________________________________________________________
int ElementStreamer(std::ostream& finalString,
                    const clang::NamedDecl &forcontext,
                    const clang::QualType &qti,
                    const char *t,
                    int rwmode,
                    const cling::Interpreter &interp,
                    const char *tcl=0);

//______________________________________________________________________________
bool IsBase(const clang::CXXRecordDecl *cl, const clang::CXXRecordDecl *base, const clang::CXXRecordDecl *context,const cling::Interpreter &interp);

//______________________________________________________________________________
bool IsBase(const clang::FieldDecl &m, const char* basename, const cling::Interpreter &interp);

//______________________________________________________________________________
bool HasCustomOperatorNewArrayPlacement(clang::RecordDecl const&, const cling::Interpreter &interp);

//______________________________________________________________________________
bool HasCustomOperatorNewPlacement(char const*, clang::RecordDecl const&, const cling::Interpreter&);

//______________________________________________________________________________
bool HasCustomOperatorNewPlacement(clang::RecordDecl const&, const cling::Interpreter&);

//______________________________________________________________________________
bool HasDirectoryAutoAdd(clang::CXXRecordDecl const*, const cling::Interpreter&);

//______________________________________________________________________________
bool HasIOConstructor(clang::CXXRecordDecl const*, std::string&, const RConstructorTypes&, const cling::Interpreter&);

//______________________________________________________________________________
bool HasNewMerge(clang::CXXRecordDecl const*, const cling::Interpreter&);

//______________________________________________________________________________
bool HasOldMerge(clang::CXXRecordDecl const*, const cling::Interpreter&);

//______________________________________________________________________________
bool hasOpaqueTypedef(clang::QualType instanceType, const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
bool hasOpaqueTypedef(const AnnotatedRecordDecl &cl, const cling::Interpreter &interp, const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
bool HasResetAfterMerge(clang::CXXRecordDecl const*, const cling::Interpreter&);

//______________________________________________________________________________
bool NeedDestructor(clang::CXXRecordDecl const*, const cling::Interpreter&);

//______________________________________________________________________________
bool NeedTemplateKeyword(clang::CXXRecordDecl const*);

//______________________________________________________________________________
bool CheckPublicFuncWithProto(clang::CXXRecordDecl const*, char const*, char const*,
                              const cling::Interpreter&, bool diagnose);

//______________________________________________________________________________
long GetLineNumber(clang::Decl const*);

//______________________________________________________________________________
bool GetNameWithinNamespace(std::string&, std::string&, std::string&, clang::CXXRecordDecl const*);

//______________________________________________________________________________
void GetQualifiedName(std::string &qual_name, const clang::QualType &type, const clang::NamedDecl &forcontext);

//----
std::string GetQualifiedName(const clang::QualType &type, const clang::NamedDecl &forcontext);

//______________________________________________________________________________
void GetQualifiedName(std::string &qual_name, const clang::Type &type, const clang::NamedDecl &forcontext);

//----
std::string GetQualifiedName(const clang::Type &type, const clang::NamedDecl &forcontext);

//______________________________________________________________________________
void GetQualifiedName(std::string &qual_name, const clang::NamedDecl &nd);

//----
std::string GetQualifiedName(const clang::NamedDecl &nd);

//______________________________________________________________________________
void GetQualifiedName(std::string &qual_name, const AnnotatedRecordDecl &annotated);

//----
std::string GetQualifiedName(const AnnotatedRecordDecl &annotated);

//______________________________________________________________________________
void GetQualifiedName(std::string &qual_name, const clang::RecordDecl &recordDecl);

//----
std::string GetQualifiedName(const clang::RecordDecl &recordDecl);

//______________________________________________________________________________
int WriteNamespaceHeader(std::ostream&, const clang::RecordDecl *);

//______________________________________________________________________________
int WriteNamespaceHeader(std::ostream&, const clang::DeclContext *);

//______________________________________________________________________________
void WritePointersSTL(const AnnotatedRecordDecl &cl, const cling::Interpreter &interp, const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
int GetClassVersion(const clang::RecordDecl *cl, const cling::Interpreter &interp);

//______________________________________________________________________________
std::pair<bool, int> GetTrivialIntegralReturnValue(const clang::FunctionDecl *funcCV, const cling::Interpreter &interp);

//______________________________________________________________________________
int IsSTLContainer(const AnnotatedRecordDecl &annotated);

//______________________________________________________________________________
ROOT::ESTLType IsSTLContainer(const clang::FieldDecl &m);

//______________________________________________________________________________
int IsSTLContainer(const clang::CXXBaseSpecifier &base);

void foreachHeaderInModule(const clang::Module &module,
                           const std::function<void(const clang::Module::Header &)> &closure,
                           bool includeDirectlyUsedModules = true);

//______________________________________________________________________________
const char *ShortTypeName(const char *typeDesc);

//______________________________________________________________________________
std::string ShortTypeName(const clang::FieldDecl &m);

//______________________________________________________________________________
bool IsStreamableObject(const clang::FieldDecl &m, const cling::Interpreter& interp);

//______________________________________________________________________________
clang::RecordDecl *GetUnderlyingRecordDecl(clang::QualType type);

//______________________________________________________________________________
std::string TrueName(const clang::FieldDecl &m);

//______________________________________________________________________________
const clang::CXXRecordDecl *ScopeSearch(const char *name,
                                        const cling::Interpreter &gInterp,
                                        bool diagnose,
                                        const clang::Type** resultType);

//______________________________________________________________________________
void WriteAuxFunctions(std::ostream& finalString,
                       const AnnotatedRecordDecl &cl,
                       const clang::CXXRecordDecl *decl,
                       const cling::Interpreter &interp,
                       const RConstructorTypes& ctorTypes,
                       const TNormalizedCtxt &normCtxt);


//______________________________________________________________________________
const clang::FunctionDecl *GetFuncWithProto(const clang::Decl* cinfo,
                                            const char *method,
                                            const char *proto,
                                            const cling::Interpreter &gInterp,
                                            bool diagnose);

//______________________________________________________________________________
void WriteClassCode(CallWriteStreamer_t WriteStreamerFunc,
                    const AnnotatedRecordDecl &cl,
                    const cling::Interpreter &interp,
                    const TNormalizedCtxt &normCtxt,
                    std::ostream& finalString,
                    const RConstructorTypes& ctorTypes,
                    bool isGenreflex);

//______________________________________________________________________________
void WriteClassInit(std::ostream& finalString,
                    const AnnotatedRecordDecl &cl,
                    const clang::CXXRecordDecl *decl,
                    const cling::Interpreter &interp,
                    const TNormalizedCtxt &normCtxt,
                    const RConstructorTypes& ctorTypes,
                    bool& needCollectionProxy);

//______________________________________________________________________________
bool HasCustomStreamerMemberFunction(const AnnotatedRecordDecl &cl,
                                     const clang::CXXRecordDecl* clxx,
                                     const cling::Interpreter &interp,
                                     const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
bool HasCustomConvStreamerMemberFunction(const AnnotatedRecordDecl &cl,
                                         const clang::CXXRecordDecl* clxx,
                                         const cling::Interpreter &interp,
                                         const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
// Return the header file to be included to declare the Decl
llvm::StringRef GetFileName(const clang::Decl& decl,
                            const cling::Interpreter& interp);

//______________________________________________________________________________
// Return the dictionary file name for a module
std::string GetModuleFileName(const char* moduleName);

//______________________________________________________________________________
// Return (in the argument 'output') a mangled version of the C++ symbol/type (pass as 'input')
// that can be used in C++ as a variable name.
void GetCppName(std::string &output, const char *input);

//______________________________________________________________________________
// Return the type with all parts fully qualified (most typedefs),
// including template arguments, appended to name.
void GetFullyQualifiedTypeName(std::string &name, const clang::QualType &type, const cling::Interpreter &interpreter);

//______________________________________________________________________________
// Return the type with all parts fully qualified (most typedefs),
// including template arguments, appended to name, without using the interpreter
void GetFullyQualifiedTypeName(std::string &name, const clang::QualType &type, const clang::ASTContext &);

//______________________________________________________________________________
// Return the type normalized for ROOT,
// keeping only the ROOT opaque typedef (Double32_t, etc.) and
// adding default template argument for all types except those explicitly
// requested to be drop by the user.
// Default template for STL collections are not yet removed by this routine.
clang::QualType GetNormalizedType(const clang::QualType &type, const cling::Interpreter &interpreter, const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
// Return the type name normalized for ROOT,
// keeping only the ROOT opaque typedef (Double32_t, etc.) and
// adding default template argument for all types except the STL collections
// where we remove the default template argument if any.
void GetNormalizedName(std::string &norm_name, const clang::QualType &type, const cling::Interpreter &interpreter, const TNormalizedCtxt &normCtxt);

//______________________________________________________________________________
// Alternative signature
void GetNormalizedName(std::string &norm_name,
                       const clang::TypeDecl* typeDecl,
                       const cling::Interpreter &interpreter);

//______________________________________________________________________________
// Analog to GetNameForIO but with types.
// It uses the LookupHelper of Cling to transform the name in type.
clang::QualType GetTypeForIO(const clang::QualType& templateInstanceType,
                             const cling::Interpreter &interpreter,
                             const TNormalizedCtxt &normCtxt,
                             TClassEdit::EModType mode = TClassEdit::kNone);

//______________________________________________________________________________
// Get the name and the type for the IO given a certain type. In some sense the
// combination of GetNameForIO and GetTypeForIO.
std::pair<std::string,clang::QualType> GetNameTypeForIO(const clang::QualType& templateInstanceType,
                                                        const cling::Interpreter &interpreter,
                                                        const TNormalizedCtxt &normCtxt,
                                                        TClassEdit::EModType mode = TClassEdit::kNone);

//______________________________________________________________________________
// Returns comment in a meaningful way
llvm::StringRef GetComment(const clang::Decl &decl, clang::SourceLocation *loc = 0);

//______________________________________________________________________________
// Returns the comment of the ClassDef macro
llvm::StringRef GetClassComment(const clang::CXXRecordDecl &decl, clang::SourceLocation *loc, const cling::Interpreter &interpreter);

//______________________________________________________________________________
// Return the base/underlying type of a chain of array or pointers type.
const clang::Type *GetUnderlyingType(clang::QualType type);

//______________________________________________________________________________
// Scans the redeclaration chain for an annotation.
//
// returns 0 if no annotation was found.
//
template<typename T>
const T* GetAnnotatedRedeclarable(const T* Redecl) {
   if (!Redecl)
      return 0;

   Redecl = Redecl->getMostRecentDecl();
   while (Redecl && !Redecl->hasAttrs())
      Redecl = Redecl->getPreviousDecl();

   return Redecl;
}

//______________________________________________________________________________
// Overload the template for typedefs, because they don't contain
// isThisDeclarationADefinition method. (Use inline to avoid violating ODR)
const clang::TypedefNameDecl* GetAnnotatedRedeclarable(const clang::TypedefNameDecl* TND);

//______________________________________________________________________________
// Overload the template for tags, because we only check definitions.
const clang::TagDecl* GetAnnotatedRedeclarable(const clang::TagDecl* TND);

//______________________________________________________________________________
// Return true if the decl is part of the std namespace.
bool IsStdClass(const clang::RecordDecl &cl);

//______________________________________________________________________________
// Return true, if the decl is part of the std namespace and we want
// its default parameter dropped.
bool IsStdDropDefaultClass(const clang::RecordDecl &cl);

//______________________________________________________________________________
// See if the CXXRecordDecl matches the current of any of the previous CXXRecordDecls
bool MatchWithDeclOrAnyOfPrevious(const clang::CXXRecordDecl &cl, const clang::CXXRecordDecl &currentCl);

//______________________________________________________________________________
// Return true if the decl is of type
bool IsOfType(const clang::CXXRecordDecl &cl, const std::string& type, const cling::LookupHelper& lh);

//______________________________________________________________________________
// Return which kind of STL container the decl is, if any.
ROOT::ESTLType IsSTLCont(const clang::RecordDecl &cl);

//______________________________________________________________________________
// Check if 'input' or any of its template parameter was substituted when
// instantiating the class template instance and replace it with the
// partially sugared type we have from 'instance'.
clang::QualType ReSubstTemplateArg(clang::QualType input, const clang::Type *instance);

//______________________________________________________________________________
// Remove the last n template arguments from the name
int RemoveTemplateArgsFromName(std::string& name, unsigned int);

//______________________________________________________________________________
clang::TemplateName ExtractTemplateNameFromQualType(const clang::QualType& qt);

//______________________________________________________________________________
bool QualType2Template(const clang::QualType& qt,
                       clang::ClassTemplateDecl*& ctd,
                       clang::ClassTemplateSpecializationDecl*& ctsd);

//______________________________________________________________________________
clang::ClassTemplateDecl* QualType2ClassTemplateDecl(const clang::QualType& qt);

//______________________________________________________________________________
// Extract the namespaces enclosing a DeclContext
void ExtractCtxtEnclosingNameSpaces(const clang::DeclContext&,
                                    std::list<std::pair<std::string,bool> >&);
//______________________________________________________________________________
void ExtractEnclosingNameSpaces(const clang::Decl&,
                                std::list<std::pair<std::string,bool> >&);

//______________________________________________________________________________
const clang::RecordDecl* ExtractEnclosingScopes(const clang::Decl& decl,
                                          std::list<std::pair<std::string,unsigned int> >& enclosingSc);
//______________________________________________________________________________
// Kind of stl container
ROOT::ESTLType STLKind(const llvm::StringRef type);

//______________________________________________________________________________
// Set the toolchain and the include paths for relocatability
void SetPathsForRelocatability(std::vector<std::string>& clingArgs);

//______________________________________________________________________________
void ReplaceAll(std::string& str, const std::string& from, const std::string& to, bool recurse=false);

// Functions for the printouts -------------------------------------------------

//______________________________________________________________________________
inline unsigned int &GetNumberOfErrors()
{
   static unsigned int gNumberOfErrors = 0;
   return gNumberOfErrors;
}

//______________________________________________________________________________
// True if printing a warning should increase GetNumberOfErrors
inline bool &GetWarningsAreErrors()
{
   static bool gWarningsAreErrors = false;
   return gWarningsAreErrors;
}

//______________________________________________________________________________
// Inclusive minimum error level a message needs to get handled
inline int &GetErrorIgnoreLevel() {
   static int gErrorIgnoreLevel = ROOT::TMetaUtils::kError;
   return gErrorIgnoreLevel;
}

//______________________________________________________________________________
inline void LevelPrint(bool prefix, int level, const char *location, const char *fmt, va_list ap)
{
   if (level < GetErrorIgnoreLevel())
      return;

   const char *type = 0;

   if (level >= ROOT::TMetaUtils::kInfo)
      type = "Info";
   if (level >= ROOT::TMetaUtils::kNote)
      type = "Note";
   if (level >= ROOT::TMetaUtils::kWarning)
      type = "Warning";
   if (level >= ROOT::TMetaUtils::kError)
      type = "Error";
   if (level >= ROOT::TMetaUtils::kSysError)
      type = "SysError";
   if (level >= ROOT::TMetaUtils::kFatal)
      type = "Fatal";

   if (!location || !location[0]) {
      if (prefix) fprintf(stderr, "%s: ", type);
      vfprintf(stderr, (const char*)va_(fmt), ap);
   } else {
      if (prefix) fprintf(stderr, "%s in <%s>: ", type, location);
      else fprintf(stderr, "In <%s>: ", location);
      vfprintf(stderr, (const char*)va_(fmt), ap);
   }

   fflush(stderr);

   // Keep track of the warnings/errors we printed.
   if (level >= ROOT::TMetaUtils::kError || (level == ROOT::TMetaUtils::kWarning && GetWarningsAreErrors())) {
      ++GetNumberOfErrors();
   }
}

//______________________________________________________________________________
// Use this function in case an error occured.
inline void Error(const char *location, const char *va_(fmt), ...)
{
   va_list ap;
   va_start(ap,va_(fmt));
   LevelPrint(true, ROOT::TMetaUtils::kError, location, va_(fmt), ap);
   va_end(ap);
}

//______________________________________________________________________________
// Use this function in case a system (OS or GUI) related error occured.
inline void SysError(const char *location, const char *va_(fmt), ...)
{
   va_list ap;
   va_start(ap, va_(fmt));
   LevelPrint(true, ROOT::TMetaUtils::kSysError, location, va_(fmt), ap);
   va_end(ap);
}

//______________________________________________________________________________
// Use this function for informational messages.
inline void Info(const char *location, const char *va_(fmt), ...)
{
   va_list ap;
   va_start(ap,va_(fmt));
   LevelPrint(true, ROOT::TMetaUtils::kInfo, location, va_(fmt), ap);
   va_end(ap);
}

//______________________________________________________________________________
// Use this function in warning situations.
inline void Warning(const char *location, const char *va_(fmt), ...)
{
   va_list ap;
   va_start(ap,va_(fmt));
   LevelPrint(true, ROOT::TMetaUtils::kWarning, location, va_(fmt), ap);
   va_end(ap);
}

//______________________________________________________________________________
// Use this function in case of a fatal error. It will abort the program.
inline void Fatal(const char *location, const char *va_(fmt), ...)
{
   va_list ap;
   va_start(ap,va_(fmt));
   LevelPrint(true, ROOT::TMetaUtils::kFatal, location, va_(fmt), ap);
   va_end(ap);
}

//______________________________________________________________________________
const std::string& GetPathSeparator();

//______________________________________________________________________________
bool EndsWith(const std::string &theString, const std::string &theSubstring);

//______________________________________________________________________________
bool BeginsWith(const std::string &theString, const std::string &theSubstring);

//______________________________________________________________________________
bool IsLinkdefFile(const char *filename);

//______________________________________________________________________________
bool IsHeaderName(const std::string &filename);

//______________________________________________________________________________
namespace AST2SourceTools {

//______________________________________________________________________________
const std::string Decls2FwdDecls(const std::vector<const clang::Decl*> &decls,
                                 bool (*ignoreFiles)(const clang::PresumedLoc&) ,
                                 const cling::Interpreter& interp);

//______________________________________________________________________________
int PrepareArgsForFwdDecl(std::string& templateArgs,
                          const clang::TemplateParameterList& tmplParamList,
                          const cling::Interpreter& interpreter);

//______________________________________________________________________________
int EncloseInNamespaces(const clang::Decl& decl, std::string& defString);

//______________________________________________________________________________
const clang::RecordDecl* EncloseInScopes(const clang::Decl& decl, std::string& defString);

//______________________________________________________________________________
int FwdDeclFromRcdDecl(const clang::RecordDecl& recordDecl,
                       const cling::Interpreter& interpreter,
                       std::string& defString,
                       bool acceptStl=false);

//______________________________________________________________________________
int FwdDeclFromTmplDecl(const clang::TemplateDecl& tmplDecl,
                        const cling::Interpreter& interpreter,
                        std::string& defString);
//______________________________________________________________________________
int GetDefArg(const clang::ParmVarDecl& par, std::string& valAsString, const clang::PrintingPolicy& pp);

//______________________________________________________________________________
int FwdDeclFromFcnDecl(const clang::FunctionDecl& fcnDecl,
                       const cling::Interpreter& interpreter,
                       std::string& defString);
//______________________________________________________________________________
int FwdDeclFromTypeDefNameDecl(const clang::TypedefNameDecl& tdnDecl,
                               const cling::Interpreter& interpreter,
                               std::string& fwdDeclString,
                               std::unordered_set<std::string>* fwdDeclSet=nullptr);

} // namespace AST2SourceTools

} // namespace TMetaUtils

} // namespace ROOT

#endif // ROOT_TMetaUtils
