// @(#)root/metautils:$Id$
// Author: Victor Perev   10/04/2003
//         Philippe Canal 05/2004

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClassEdit
#define ROOT_TClassEdit

#include <ROOT/RConfig.hxx>
#include "RConfigure.h"
#include <stdlib.h>
#ifdef R__WIN32
#ifndef UNDNAME_COMPLETE
#define UNDNAME_COMPLETE 0x0000
#endif
extern "C" {
   char *__unDName(char *demangled, const char *mangled, int out_len,
                   void * (* pAlloc )(size_t), void (* pFree )(void *),
                   unsigned short int flags);
}
#else
#include <cxxabi.h>
#endif
#include <string>
#include <vector>
#include <array>

#include "ESTLType.h"

#ifdef R__OLDHPACC
namespace std {
   using ::string;
   using ::vector;
}
#endif

#if defined(__CYGWIN__)
// std::to_string is missing on cygwin with gcc 4.8.2-2 and 4.8.3
#include <sstream>
namespace std {
  template <typename T>
  string to_string(T value) {
    ostringstream os;
    os << value;
    return os.str();
  }
}
#endif

namespace cling {
   class Interpreter;
}
namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}
#include "ROOT/RStringView.hxx"

// TClassEdit is used to manipulate class and type names.
//
// This class does not dependent on any other ROOT facility
// so that it can be used by rootcint.

namespace TClassEdit {

   enum EModType {
      kNone             = 0,
      kDropTrailStar    = 1<<0,
      kDropDefaultAlloc = 1<<1,
      kDropAlloc        = 1<<2,
      kInnerClass       = 1<<3,
      kInnedMostClass   = 1<<4,
      kDropStlDefault   = 1<<5, /* implies kDropDefaultAlloc */
      kDropComparator   = 1<<6, /* if the class has a comparator, drops BOTH the comparator and the Allocator */
      kDropAllDefault   = 1<<7, /* Drop default template parameter even in non STL classes */
      kLong64           = 1<<8, /* replace all 'long long' with Long64_t. */
      kDropStd          = 1<<9, /* Drop any std:: */
      kKeepOuterConst   = 1<<10,/* Make sure to keep the const keyword even outside the template parameters */
      kResolveTypedef   = 1<<11,/* Strip all typedef except Double32_t and co. */
      kDropPredicate    = 1<<12,/* Drop the predicate if applies to the collection */
      kDropHash         = 1<<13 /* Drop the hash if applies to the collection */
   };

   enum ESTLType {
      kNotSTL            = ROOT::kNotSTL,
      kVector            = ROOT::kSTLvector,
      kList              = ROOT::kSTLlist,
      kForwardlist       = ROOT::kSTLforwardlist,
      kDeque             = ROOT::kSTLdeque,
      kMap               = ROOT::kSTLmap,
      kMultiMap          = ROOT::kSTLmultimap,
      kSet               = ROOT::kSTLset,
      kMultiSet          = ROOT::kSTLmultiset,
      kUnorderedSet      = ROOT::kSTLunorderedset,
      kUnorderedMultiSet = ROOT::kSTLunorderedmultiset,
      kUnorderedMap      = ROOT::kSTLunorderedmap,
      kUnorderedMultiMap = ROOT::kSTLunorderedmultimap,
      kBitSet            = ROOT::kSTLbitset,
      kEnd               = ROOT::kSTLend
   };

   enum class EComplexType : short {
      kNone,
      kDouble,
      kFloat,
      kInt,
      kLong
   };

   EComplexType GetComplexType(const char*);

   class TInterpreterLookupHelper {
   public:
      TInterpreterLookupHelper() { }
      virtual ~TInterpreterLookupHelper();

      virtual bool ExistingTypeCheck(const std::string & /*tname*/,
                                     std::string & /*result*/) = 0;
      virtual void GetPartiallyDesugaredName(std::string & /*nameLong*/) = 0;
      virtual bool IsAlreadyPartiallyDesugaredName(const std::string & /*nondef*/,
                                                   const std::string & /*nameLong*/) = 0;
      virtual bool IsDeclaredScope(const std::string & /*base*/, bool & /*isInlined*/) = 0;
      virtual bool GetPartiallyDesugaredNameWithScopeHandling(const std::string & /*tname*/,
                                                              std::string & /*result*/,
                                                              bool /* dropstd */ = true) = 0;
      virtual void ShuttingDownSignal() = 0;
   };

   struct TSplitType {

      const char *fName; // Original spelling of the name.
      std::vector<std::string> fElements;
      int fNestedLocation; // Stores the location of the tail (nested names) in nestedLoc (0 indicates no tail).

      TSplitType(const char *type2split, EModType mode = TClassEdit::kNone);

      int  IsSTLCont(int testAlloc=0) const;
      ROOT::ESTLType  IsInSTL() const;
      void ShortType(std::string &answer, int mode);
      bool IsTemplate();

   private:
      TSplitType(const TSplitType&) = delete;
      TSplitType &operator=(const TSplitType &) = delete;
   };

   void        Init(TClassEdit::TInterpreterLookupHelper *helper);

   std::string CleanType (const char *typeDesc,int mode = 0,const char **tail=0);
   inline bool IsArtificial(std::string_view name) { return name.find('@') != name.npos; }
   inline bool IsArtificial(ROOT::Internal::TStringView name) {return IsArtificial(std::string_view(name)); }
   bool        IsDefAlloc(const char *alloc, const char *classname);
   bool        IsDefAlloc(const char *alloc, const char *keyclassname, const char *valueclassname);
   bool        IsDefComp (const char *comp , const char *classname);
   bool        IsDefPred(const char *predname, const char *classname);
   bool        IsDefHash(const char *hashname, const char *classname);
   bool        IsInterpreterDetail(const char *type);
   bool        IsSTLBitset(const char *type);
   ROOT::ESTLType UnderlyingIsSTLCont(std::string_view type);
   inline ROOT::ESTLType UnderlyingIsSTLCont (ROOT::Internal::TStringView type) { return UnderlyingIsSTLCont(std::string_view(type)); }
   ROOT::ESTLType IsSTLCont (std::string_view type);
   inline ROOT::ESTLType IsSTLCont (ROOT::Internal::TStringView type) { return IsSTLCont(std::string_view(type)); }
   int         IsSTLCont (const char *type,int testAlloc);
   bool        IsStdClass(const char *type);
   bool        IsVectorBool(const char *name);
   void        GetNormalizedName(std::string &norm_name, std::string_view name);
   inline void GetNormalizedName (std::string &norm_name, ROOT::Internal::TStringView name) { return GetNormalizedName(norm_name, std::string_view(name)); }
   std::string GetLong64_Name(const char *original);
   std::string GetLong64_Name(const std::string& original);
   int         GetSplit  (const char *type, std::vector<std::string> &output, int &nestedLoc, EModType mode = TClassEdit::kNone);
   ROOT::ESTLType STLKind(std::string_view type);    //Kind of stl container
   inline ROOT::ESTLType STLKind(ROOT::Internal::TStringView type) { return STLKind(std::string_view(type)); }
   int         STLArgs   (int kind);            //Min number of arguments without allocator
   std::string ResolveTypedef(const char *tname, bool resolveAll = false);
   std::string ShortType (const char *typeDesc, int mode);
   std::string InsertStd(const char *tname);
   const char* GetUnqualifiedName(const char*name);
   inline bool IsUniquePtr(std::string_view name) {return 0 == name.compare(0, 11, "unique_ptr<");}
   inline bool IsUniquePtr(ROOT::Internal::TStringView name) {return IsUniquePtr(std::string_view(name)); }
   inline bool IsStdArray(std::string_view name) {return 0 == name.compare(0, 6, "array<");}
   inline bool IsStdArray(ROOT::Internal::TStringView name) {return IsStdArray(std::string_view(name)); }
   inline bool IsStdPair(std::string_view name)
   {
      return 0 == name.compare(0, 10, "std::pair<") || 0 == name.compare(0, 5, "pair<");
   }
   inline bool IsStdPair(ROOT::Internal::TStringView name) {return IsStdPair(std::string_view(name)); }
   inline bool IsStdPairBase(std::string_view name)
   {
      return 0 == name.compare(0, 17, "std::__pair_base<") || 0 == name.compare(0, 12, "__pair_base<");
   }
   inline bool IsStdPairBase(ROOT::Internal::TStringView name) {return IsStdPair(std::string_view(name)); }
   inline std::string GetUniquePtrType(std::string_view name)
   {
      // Find the first template parameter
      std::vector<std::string> v;
      int i;
      GetSplit(name.data(), v, i);
      return v[1];
   }
   inline std::string GetUniquePtrType(ROOT::Internal::TStringView name) {return GetUniquePtrType(std::string_view(name)); }
   std::string GetNameForIO(const std::string& templateInstanceName,
                           TClassEdit::EModType mode = TClassEdit::kNone,
                           bool* hasChanged = nullptr);
   bool GetStdArrayProperties(const char* typeName,
                              std::string& typeNameBuf,
                              std::array<int, 5>& maxIndices,
                              int& ndim);

   inline char* DemangleName(const char* mangled_name, int& errorCode)
   {
   // Demangle in a portable way the name.
   // IMPORTANT: The caller is responsible for freeing the returned const char*

   errorCode=0;
#ifdef R__WIN32
   char *demangled_name = __unDName(0, mangled_name, 0, malloc, free, UNDNAME_COMPLETE);
   if (!demangled_name) {
      errorCode = -1;
      return nullptr;
   }
   std::string demangledName = demangled_name;
   if (demangledName.compare(0, 6, "class ") == 0)
      demangledName.erase(0, 6);
   else if (demangledName.compare(0, 7, "struct ") == 0)
      demangledName.erase(0, 7);
   strcpy(demangled_name, demangledName.c_str());
#else
   char *demangled_name = abi::__cxa_demangle(mangled_name, 0, 0, &errorCode);
   if (!demangled_name || errorCode) {
      free(demangled_name);
      return nullptr;
   }
#endif
   return demangled_name;
   }
   char* DemangleTypeIdName(const std::type_info& ti, int& errorCode);


   /// Result of splitting a function declaration into
   /// fReturnType fScopeName::fFunctionName<fFunctionTemplateArguments>(fFunctionParameters)
   struct FunctionSplitInfo {
      /// Return type of the function, might be empty if the function declaration string did not provide it.
      std::string fReturnType;

      /// Name of the scope qualification of the function, possibly empty
      std::string fScopeName;

      /// Name of the function
      std::string fFunctionName;

      /// Template arguments of the function template specialization, if any; will contain one element "" for
      /// `function<>()`
      std::vector<std::string> fFunctionTemplateArguments;

      /// Function parameters.
      std::vector<std::string> fFunctionParameters;
   };

   /// Split a function declaration into its different parts.
   bool SplitFunction(std::string_view decl, FunctionSplitInfo &result);
}

#endif
