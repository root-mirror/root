// @(#)root/metautils:$Id$
/// \file TClassEdit.cxx
/// \ingroup Base
/// \author Victor Perev
/// \author Philippe Canal
/// \date 04/10/2003

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "TClassEdit.h"
#include <cctype>
#include "Rstrstream.h"
#include <set>
#include <stack>
// for shared_ptr
#include <memory>
#include "ROOT/RStringView.hxx"
#include <algorithm>

using namespace std;

namespace {
   static TClassEdit::TInterpreterLookupHelper *gInterpreterHelper = nullptr;

   template <typename T>
   struct ShuttingDownSignaler : public T {
      using T::T;

      ~ShuttingDownSignaler()
      {
         if (gInterpreterHelper)
            gInterpreterHelper->ShuttingDownSignal();
      }
   };
}

////////////////////////////////////////////////////////////////////////////////
/// Return the length, if any, taken by std:: and any
/// potential inline namespace (well compiler detail namespace).

static size_t StdLen(const std::string_view name)
{
   size_t len = 0;
   if (name.compare(0,5,"std::")==0) {
      len = 5;

      // TODO: This is likely to induce unwanted autoparsing, those are reduced
      // by the caching of the result.
      if (gInterpreterHelper) {
         for(size_t i = 5; i < name.length(); ++i) {
            if (name[i] == '<') break;
            if (name[i] == ':') {
               bool isInlined;
               std::string scope(name.data(),i);
               std::string scoperesult;
               // We assume that we are called in already serialized code.
               // Note: should we also cache the negative answers?
               static ShuttingDownSignaler<std::set<std::string>> gInlined;

               if (gInlined.find(scope) != gInlined.end()) {
                  len = i;
                  if (i+1<name.length() && name[i+1]==':') {
                     len += 2;
                  }
               }
               if (!gInterpreterHelper->ExistingTypeCheck(scope, scoperesult)
                   && gInterpreterHelper->IsDeclaredScope(scope,isInlined)) {
                  if (isInlined) {
                     gInlined.insert(scope);
                     len = i;
                     if (i+1<name.length() && name[i+1]==':') {
                        len += 2;
                     }
                  }
               }
            }
         }
      }
   }

   return len;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove std:: and any potential inline namespace (well compiler detail
/// namespace.

static void RemoveStd(std::string &name, size_t pos = 0)
{
   size_t len = StdLen({name.data()+pos,name.length()-pos});
   if (len) {
      name.erase(pos,len);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove std:: and any potential inline namespace (well compiler detail
/// namespace.

static void RemoveStd(std::string_view &name)
{
   size_t len = StdLen(name);
   if (len) {
      name.remove_prefix(len);
   }
}

////////////////////////////////////////////////////////////////////////////////

TClassEdit::EComplexType TClassEdit::GetComplexType(const char* clName)
{
   if (0 == strncmp(clName, "complex<", 8)) {
      const char *clNamePlus8 = clName + 8;
      if (0 == strcmp("float>", clNamePlus8)) {
         return EComplexType::kFloat;
         }
      if (0 == strcmp("double>", clNamePlus8)) {
         return EComplexType::kDouble;
      }
      if (0 == strcmp("int>", clNamePlus8)) {
         return EComplexType::kInt;
      }
      if (0 == strcmp("long>", clNamePlus8)) {
         return EComplexType::kLong;
      }
   }
   return EComplexType::kNone;
}

////////////////////////////////////////////////////////////////////////////////
TClassEdit::TInterpreterLookupHelper::~TInterpreterLookupHelper()
{
   // Already too late to call this->ShuttingDownSignal
   // the virtual table has already lost (on some platform) the
   // address of the derived function that we would need to call.
   // But at least forget about this instance!

   if (this == gInterpreterHelper)
      gInterpreterHelper = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void TClassEdit::Init(TClassEdit::TInterpreterLookupHelper *helper)
{
   gInterpreterHelper = helper;
}

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TClassEdit::TSplitType::TSplitType(const char *type2split, EModType mode) : fName(type2split), fNestedLocation(0)
{
   TClassEdit::GetSplit(type2split, fElements, fNestedLocation, mode);
}

////////////////////////////////////////////////////////////////////////////////
///  type     : type name: vector<list<classA,allocator>,allocator>[::iterator]
///  result:    0          : not stl container and not declared inside an stl container.
///             result: code of container that the type or is the scope of the type

ROOT::ESTLType TClassEdit::TSplitType::IsInSTL() const
{
   if (fElements[0].empty()) return ROOT::kNotSTL;
   return STLKind(fElements[0]);
}

////////////////////////////////////////////////////////////////////////////////
///  type     : type name: vector<list<classA,allocator>,allocator>
///  testAlloc: if true, we test allocator, if it is not default result is negative
///  result:    0          : not stl container
///             abs(result): code of container 1=vector,2=list,3=deque,4=map
///                           5=multimap,6=set,7=multiset
///             positive val: we have a vector or list with default allocator to any depth
///                   like vector<list<vector<int>>>
///             negative val: STL container other than vector or list, or non default allocator
///                           For example: vector<deque<int>> has answer -1

int TClassEdit::TSplitType::IsSTLCont(int testAlloc) const
{

   if (fElements[0].empty()) return 0;
   int numb = fElements.size();
   if (!fElements[numb-1].empty() && fElements[numb-1][0]=='*') --numb;

   if ( fNestedLocation ) {
      // The type has been defined inside another namespace and/or class
      // this couldn't possibly be an STL container
      return 0;
   }

   int kind = STLKind(fElements[0]);

   if (kind==ROOT::kSTLvector || kind==ROOT::kSTLlist || kind==ROOT::kSTLforwardlist) {

      int nargs = STLArgs(kind);
      if (testAlloc && (numb-1 > nargs) && !IsDefAlloc(fElements[numb-1].c_str(),fElements[1].c_str())) {

         // We have a non default allocator,
         // let's return a negative value.

         kind = -kind;

      } else {

         // We has a default allocator, let's continue to
         // look inside the argument list.
         int k = TClassEdit::IsSTLCont(fElements[1].c_str(),testAlloc);
         if (k<0) kind = -kind;

      }
   }

   // We return a negative value for anything which is not a vector or a list.
   if(kind>2) kind = - kind;
   return kind;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/// Return the absolute type of typeDesc into the string answ.

void TClassEdit::TSplitType::ShortType(std::string &answ, int mode)
{
   // E.g.: typeDesc = "class const volatile TNamed**", returns "TNamed**".
   // if (mode&1) remove last "*"s                     returns "TNamed"
   // if (mode&2) remove default allocators from STL containers
   // if (mode&4) remove all     allocators from STL containers
   // if (mode&8) return inner class of stl container. list<innerClass>
   // if (mode&16) return deepest class of stl container. vector<list<deepest>>
   // if (mode&kDropAllDefault) remove default template arguments
   /////////////////////////////////////////////////////////////////////////////

   answ.clear();
   int narg = fElements.size();
   int tailLoc = 0;

   if (narg == 0) {
      answ = fName;
      return ;
   }
   //      fprintf(stderr,"calling ShortType %d for %s with narg %d\n",mode,typeDesc,narg);
   //      {for (int i=0;i<narg;i++) fprintf(stderr,"calling ShortType %d for %s with %d %s \n",
   //                                        mode,typeDesc,i,arglist[i].c_str());
   //      }
   if (fElements[narg-1].empty() == false &&
       (fElements[narg-1][0]=='*'
        || fElements[narg-1][0]=='&'
        || fElements[narg-1][0]=='['
        || 0 == fElements[narg-1].compare(0,6,"const*")
        || 0 == fElements[narg-1].compare(0,6,"const&")
        || 0 == fElements[narg-1].compare(0,6,"const[")
        || 0 == fElements[narg-1].compare("const")
        )
       ) {
      if ((mode&1)==0) tailLoc = narg-1;
   }
   else { assert(fElements[narg-1].empty()); };
   narg--;
   mode &= (~1);

   if (fNestedLocation) narg--;

   //    fprintf(stderr,"calling ShortType %d for %s with narg %d tail %d\n",imode,typeDesc,narg,tailLoc);

   //kind of stl container
   const int kind = STLKind(fElements[0]);
   const int iall = STLArgs(kind);

   // Only class is needed
   if (mode&(8|16)) {
      while(narg-1>iall) { fElements.pop_back(); narg--;}
      if (!fElements[0].empty() && tailLoc) {
         tailLoc = 0;
      }
      fElements[0].clear();
      mode&=(~8);
   }

   if (mode & kDropAllDefault) mode |= kDropStlDefault;
   if (mode & kDropStlDefault) mode |= kDropDefaultAlloc;

   if (kind) {
      bool allocRemoved = false;

      if ( mode & (kDropDefaultAlloc|kDropAlloc) ) {
         // remove allocators


         if (narg-1 == iall+1) {
            // has an allocator specified
            bool dropAlloc = false;
            if (mode & kDropAlloc) {

               dropAlloc = true;

            } else if (mode & kDropDefaultAlloc) {
               switch (kind) {
                  case ROOT::kSTLvector:
                  case ROOT::kSTLlist:
                  case ROOT::kSTLforwardlist:
                  case ROOT::kSTLdeque:
                  case ROOT::kSTLset:
                  case ROOT::kSTLmultiset:
                  case ROOT::kSTLunorderedset:
                  case ROOT::kSTLunorderedmultiset:
                     dropAlloc = IsDefAlloc(fElements[iall+1].c_str(),fElements[1].c_str());
                     break;
                  case ROOT::kSTLmap:
                  case ROOT::kSTLmultimap:
                  case ROOT::kSTLunorderedmap:
                  case ROOT::kSTLunorderedmultimap:
                     dropAlloc = IsDefAlloc(fElements[iall+1].c_str(),fElements[1].c_str(),fElements[2].c_str());
                     break;
                  default:
                     dropAlloc = false;
               }

            }
            if (dropAlloc) {
               narg--;
               allocRemoved = true;
            }
         } else {
            // has no allocator specified (hence it is already removed!)
            allocRemoved = true;
         }
      }

      if ( allocRemoved && (mode & kDropStlDefault) && narg-1 == iall) { // remove default comparator
         if ( IsDefComp( fElements[iall].c_str(), fElements[1].c_str() ) ) {
            narg--;
         }
      } else if ( mode & kDropComparator ) {

         switch (kind) {
            case ROOT::kSTLvector:
            case ROOT::kSTLlist:
            case ROOT::kSTLforwardlist:
            case ROOT::kSTLdeque:
               break;
            case ROOT::kSTLset:
            case ROOT::kSTLmultiset:
            case ROOT::kSTLmap:
            case ROOT::kSTLmultimap:
               if (!allocRemoved && narg-1 == iall+1) {
                  narg--;
                  allocRemoved = true;
               }
               if (narg-1 == iall) narg--;
               break;
            default:
               break;
         }
      }

      // Treat now Pred and Hash for unordered set/map containers. Signature is:
      // template < class Key,
      //            class Hash = hash<Key>,
      //            class Pred = equal_to<Key>,
      //            class Alloc = allocator<Key>
      //          > class unordered_{set,multiset}
      // template < class Key,
      //            class Val,
      //            class Hash = hash<Key>,
      //            class Pred = equal_to<Key>,
      //            class Alloc = allocator<Key>
      //          > class unordered_{map,multimap}


      if (kind == ROOT::kSTLunorderedset || kind == ROOT::kSTLunorderedmultiset || kind == ROOT::kSTLunorderedmap || kind == ROOT::kSTLunorderedmultimap){

         bool predRemoved = false;

         if ( allocRemoved && (mode & kDropStlDefault) && narg-1 == iall) { // remove default predicate
            if ( IsDefPred( fElements[iall].c_str(), fElements[1].c_str() ) ) {
               predRemoved=true;
               narg--;
            }
         }

         if ( predRemoved && (mode & kDropStlDefault) && narg == iall) { // remove default hash
            if ( IsDefHash( fElements[iall-1].c_str(), fElements[1].c_str() ) ) {
               narg--;
            }
         }
      }
   } // End of treatment of stl containers
   else {
      if ( (mode & kDropStlDefault) && (narg >= 3)) {
         unsigned int offset = (0==strncmp("const ",fElements[0].c_str(),6)) ? 6 : 0;
         offset += (0==strncmp("std::",fElements[0].c_str()+offset,5)) ? 5 : 0;
         if (0 == strcmp(fElements[0].c_str()+offset,"__shared_ptr"))
         {
#ifdef _CONCURRENCE_H
            static const std::string sharedPtrDef = std::to_string(__gnu_cxx::__default_lock_policy); // to_string is C++11
#else
            static const std::string sharedPtrDef = std::to_string(2); // to_string is C++11
#endif
            if (fElements[2] == sharedPtrDef) {
               narg--;
            }
         }
      }
   }

   //   do the same for all inside
   for (int i=1;i<narg; i++) {
      if (strchr(fElements[i].c_str(),'<')==0) {
         if (mode&kDropStd) {
            unsigned int offset = (0==strncmp("const ",fElements[i].c_str(),6)) ? 6 : 0;
            RemoveStd( fElements[i], offset );
         }
         if (mode&kResolveTypedef) {
            fElements[i] = ResolveTypedef(fElements[i].c_str(),true);
         }
         continue;
      }
      fElements[i] = TClassEdit::ShortType(fElements[i].c_str(),mode | TClassEdit::kKeepOuterConst);
      if (mode&kResolveTypedef) {
         // We 'just' need to check whether the outer type is a typedef or not;
         // this also will add the default template parameter if any needs to
         // be added.
         string typeresult;
         if (gInterpreterHelper &&
             (gInterpreterHelper->ExistingTypeCheck(fElements[i], typeresult)
              || gInterpreterHelper->GetPartiallyDesugaredNameWithScopeHandling(fElements[i], typeresult))) {
            if (!typeresult.empty()) fElements[i] = typeresult;
         }
      }
   }

   unsigned int tailOffset = 0;
   if (tailLoc && fElements[tailLoc].compare(0,5,"const") == 0) {
      if (mode & kKeepOuterConst) answ += "const ";
      tailOffset = 5;
   }
   if (!fElements[0].empty()) {answ += fElements[0]; answ +="<";}

#if 0
   // This code is no longer use, the moral equivalent would be to get
   // the 'fixed' number of argument the user told us to ignore and drop those.
   // However, the name we get here might be (usually) normalized enough that
   // this is not necessary (at the very least nothing break in roottest without
   // the aforementioned new code or this old code).
   if (mode & kDropAllDefault) {
      int nargNonDefault = 0;
      std::string nonDefName = answ;
      // "superlong" because tLong might turn fName into an even longer name
      std::string nameSuperLong = fName;
      if (gInterpreterHelper)
         gInterpreterHelper->GetPartiallyDesugaredName(nameSuperLong);
      while (++nargNonDefault < narg) {
         // If T<a> is a "typedef" (aka default template params)
         // to T<a,b> then we can strip the "b".
         const char* closeTemplate = " >";
         if (nonDefName[nonDefName.length() - 1] != '>')
            ++closeTemplate;
         string nondef = nonDefName + closeTemplate;
         if (gInterpreterHelper &&
             gInterpreterHelper->IsAlreadyPartiallyDesugaredName(nondef, nameSuperLong))
            break;
         if (nargNonDefault>1) nonDefName += ",";
         nonDefName += fElements[nargNonDefault];
      }
      if (nargNonDefault < narg)
         narg = nargNonDefault;
   }
#endif

   { for (int i=1;i<narg-1; i++) { answ += fElements[i]; answ+=",";} }
   if (narg>1) { answ += fElements[narg-1]; }

   if (!fElements[0].empty()) {
      if ( answ.at(answ.size()-1) == '>') {
         answ += " >";
      } else {
         answ += '>';
      }
   }
   if (fNestedLocation) {
      // Treat X pf A<B>::X
      fElements[fNestedLocation] = TClassEdit::ShortType(fElements[fNestedLocation].c_str(),mode);
      answ += fElements[fNestedLocation];
   }
   // tail is not a type name, just [2], &, * etc.
   if (tailLoc) answ += fElements[tailLoc].c_str()+tailOffset;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the type is a template
bool TClassEdit::TSplitType::IsTemplate()
{
   return !fElements[0].empty();
}

////////////////////////////////////////////////////////////////////////////////
/// Converts STL container name to number. vector -> 1, etc..
/// If len is greater than 0, only look at that many characters in the string.

ROOT::ESTLType TClassEdit::STLKind(std::string_view type)
{
   if (type.length() == 0)
      return ROOT::kNotSTL;
   size_t offset = 0;
   if (type.compare(0,6,"const ")==0) { offset += 6; }
   offset += StdLen(type.substr(offset));
   const auto len = type.length() - offset;
   if (len == 0)
      return ROOT::kNotSTL;

   //container names
   static const char *stls[] =
      { "any", "vector", "list", "deque", "map", "multimap", "set", "multiset", "bitset",
         "forward_list", "unordered_set", "unordered_multiset", "unordered_map", "unordered_multimap", 0};
   static const size_t stllen[] =
      { 3, 6, 4, 5, 3, 8, 3, 8, 6,
         12, 13, 18, 13, 18, 0};
   static const ROOT::ESTLType values[] =
      {  ROOT::kNotSTL, ROOT::kSTLvector,
         ROOT::kSTLlist, ROOT::kSTLdeque,
         ROOT::kSTLmap, ROOT::kSTLmultimap,
         ROOT::kSTLset, ROOT::kSTLmultiset,
         ROOT::kSTLbitset,
         // New C++11
         ROOT::kSTLforwardlist,
         ROOT::kSTLunorderedset, ROOT::kSTLunorderedmultiset,
         ROOT::kSTLunorderedmap, ROOT::kSTLunorderedmultimap,
         ROOT::kNotSTL
      };

   // kind of stl container
   // find the correct ESTLType, skipping std::any (because I/O for it is not implemented yet?)
   for (int k = 1; stls[k]; ++k) {
      if (len == stllen[k]) {
         if (type.compare(offset, len, stls[k]) == 0)
            return values[k];
      }
   }
   return ROOT::kNotSTL;
}

////////////////////////////////////////////////////////////////////////////////
///      Return number of arguments for STL container before allocator

int   TClassEdit::STLArgs(int kind)
{
   static const char  stln[] =// min number of container arguments
      //     vector, list, deque, map, multimap, set, multiset, bitset,
      {    1,     1,    1,     1,   3,        3,   2,        2,      1,
      // forward_list, unordered_set, unordered_multiset, unordered_map, unordered_multimap
                    1,             3,                  3,             4,                  4};

   return stln[kind];
}

////////////////////////////////////////////////////////////////////////////////

static size_t findNameEnd(const std::string_view full)
{
   int level = 0;
   for(size_t i = 0; i < full.length(); ++i) {
      switch(full[i]) {
         case '<': { ++level; break; }
         case '>': {
            if (level == 0) return i;
            else --level;
            break;
         }
         case ',': {
            if (level == 0) return i;
            break;
         }
         default: break;
      }
   }
   return full.length();
}

////////////////////////////////////////////////////////////////////////////////

static size_t findNameEnd(const std::string &full, size_t pos)
{
   return pos + findNameEnd( {full.data()+pos,full.length()-pos} );
}

////////////////////////////////////////////////////////////////////////////////
/// return whether or not 'allocname' is the STL default allocator for type
/// 'classname'

bool TClassEdit::IsDefAlloc(const char *allocname, const char *classname)
{
   string_view a( allocname );
   RemoveStd(a);

   if (a=="alloc")                              return true;
   if (a=="__default_alloc_template<true,0>")   return true;
   if (a=="__malloc_alloc_template<0>")         return true;

   const static int alloclen = strlen("allocator<");
   if (a.compare(0,alloclen,"allocator<") != 0) {
      return false;
   }
   a.remove_prefix(alloclen);

   RemoveStd(a);

   string_view k = classname;
   RemoveStd(k);

   if (a.compare(0,k.length(),k) != 0) {
      // Now we need to compare the normalized name.
      size_t end = findNameEnd(a);

      std::string valuepart;
      GetNormalizedName(valuepart,std::string_view(a.data(),end));

      std::string norm_value;
      GetNormalizedName(norm_value,k);

      if (valuepart != norm_value) {
         return false;
      }
      a.remove_prefix(end);
   } else {
      a.remove_prefix(k.length());
   }

   if (a.compare(0,1,">")!=0 && a.compare(0,2," >")!=0) {
      return false;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// return whether or not 'allocname' is the STL default allocator for a key
/// of type 'keyclassname' and a value of type 'valueclassname'

bool TClassEdit::IsDefAlloc(const char *allocname,
                            const char *keyclassname,
                            const char *valueclassname)
{
   if (IsDefAlloc(allocname,keyclassname)) return true;

   string_view a( allocname );
   RemoveStd(a);

   const static int alloclen = strlen("allocator<");
   if (a.compare(0,alloclen,"allocator<") != 0) {
      return false;
   }
   a.remove_prefix(alloclen);

   RemoveStd(a);

   const static int pairlen = strlen("pair<");
   if (a.compare(0,pairlen,"pair<") != 0) {
      return false;
   }
   a.remove_prefix(pairlen);

   const static int constlen = strlen("const");
   if (a.compare(0,constlen+1,"const ") == 0) {
      a.remove_prefix(constlen+1);
   }

   RemoveStd(a);

   string_view k = keyclassname;
   RemoveStd(k);
   if (k.compare(0,constlen+1,"const ") == 0) {
      k.remove_prefix(constlen+1);
   }

   if (a.compare(0,k.length(),k) != 0) {
      // Now we need to compare the normalized name.
      size_t end = findNameEnd(a);

      std::string alloc_keypart;
      GetNormalizedName(alloc_keypart,std::string_view(a.data(),end));

      std::string norm_key;
      GetNormalizedName(norm_key,k);

      if (alloc_keypart != norm_key) {
         if ( norm_key[norm_key.length()-1] == '*' ) {
            // also check with a trailing 'const'.
            norm_key += "const";
         } else {
            norm_key += " const";
         }
         if (alloc_keypart != norm_key) {
           return false;
         }
      }
      a.remove_prefix(end);
   } else {
      size_t end = k.length();
      if ( (a[end-1] == '*') || a[end]==' ' ) {
         size_t skipSpace = (a[end] == ' ');
         if (a.compare(end+skipSpace,constlen,"const") == 0) {
            end += constlen+skipSpace;
         }
      }
      a.remove_prefix(end);
   }

   if (a[0] != ',') {
      return false;
   }
   a.remove_prefix(1);
   RemoveStd(a);

   string_view v = valueclassname;
   RemoveStd(v);

   if (a.compare(0,v.length(),v) != 0) {
      // Now we need to compare the normalized name.
      size_t end = findNameEnd(a);

      std::string valuepart;
      GetNormalizedName(valuepart,std::string_view(a.data(),end));

      std::string norm_value;
      GetNormalizedName(norm_value,v);

      if (valuepart != norm_value) {
         return false;
      }
      a.remove_prefix(end);
   } else {
      a.remove_prefix(v.length());
   }

   if (a.compare(0,1,">")!=0 && a.compare(0,2," >")!=0) {
      return false;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// return whether or not 'elementName' is the STL default Element for type
/// 'classname'

static bool IsDefElement(const char *elementName, const char* defaultElementName, const char *classname)
{
   string c = elementName;

   size_t pos = StdLen(c);

   const int elementlen = strlen(defaultElementName);
   if (c.compare(pos,elementlen,defaultElementName) != 0) {
      return false;
   }
   pos += elementlen;

   string k = classname;
   if (c.compare(pos,k.length(),k) != 0) {
      // Now we need to compare the normalized name.
      size_t end = findNameEnd(c,pos);

      std::string keypart;
      if (pos != end) {  // i.e. elementName != "std::less<>", see ROOT-11000.
         TClassEdit::GetNormalizedName(keypart,std::string_view(c.c_str()+pos,end-pos));
      }

      std::string norm_key;
      TClassEdit::GetNormalizedName(norm_key,k);

      if (keypart != norm_key) {
         return false;
      }
      pos = end;
   } else {
      pos += k.length();
   }

   if (c.compare(pos,1,">")!=0 && c.compare(pos,2," >")!=0) {
      return false;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// return whether or not 'compare' is the STL default comparator for type
/// 'classname'

bool TClassEdit::IsDefComp(const char *compname, const char *classname)
{
   return IsDefElement(compname, "less<", classname);
}

////////////////////////////////////////////////////////////////////////////////
/// return whether or not 'predname' is the STL default predicate for type
/// 'classname'

bool TClassEdit::IsDefPred(const char *predname, const char *classname)
{
   return IsDefElement(predname, "equal_to<", classname);
}

////////////////////////////////////////////////////////////////////////////////
/// return whether or not 'hashname' is the STL default hash for type
/// 'classname'

bool TClassEdit::IsDefHash(const char *hashname, const char *classname)
{
   return IsDefElement(hashname, "hash<", classname);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the normalized name.  See TMetaUtils::GetNormalizedName.
///
/// Return the type name normalized for ROOT,
/// keeping only the ROOT opaque typedef (Double32_t, etc.) and
/// removing the STL collections default parameter if any.
///
/// Compare to TMetaUtils::GetNormalizedName, this routines does not
/// and can not add default template parameters.

void TClassEdit::GetNormalizedName(std::string &norm_name, std::string_view name)
{
   if (name.empty()) {
      norm_name.clear();
      return;
   }

   norm_name = std::string(name); // NOTE: Is that the shortest version?

   if (TClassEdit::IsArtificial(name)) {
      // If there is a @ symbol (followed by a version number) then this is a synthetic class name created
      // from an already normalized name for the purpose of supporting schema evolution.
      return;
   }

   // Remove the std:: and default template argument and insert the Long64_t and change basic_string to string.
   TClassEdit::TSplitType splitname(norm_name.c_str(),(TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd | TClassEdit::kDropStlDefault | TClassEdit::kKeepOuterConst));
   splitname.ShortType(norm_name, TClassEdit::kDropStd | TClassEdit::kDropStlDefault | TClassEdit::kResolveTypedef | TClassEdit::kKeepOuterConst);

   // 4 elements expected: "pair", "first type name", "second type name", "trailing stars"
   if (splitname.fElements.size() == 4 && (splitname.fElements[0] == "std::pair" || splitname.fElements[0] == "pair" || splitname.fElements[0] == "__pair_base")) {
      // We don't want to lookup the std::pair itself.
      std::string first, second;
      GetNormalizedName(first, splitname.fElements[1]);
      GetNormalizedName(second, splitname.fElements[2]);
      norm_name = splitname.fElements[0] + "<" + first + "," + second;
      if (!second.empty() && second.back() == '>')
         norm_name += " >";
      else
         norm_name += ">";
      return;
   }

   // Depending on how the user typed their code, in particular typedef
   // declarations, we may end up with an explicit '::' being
   // part of the result string.  For consistency, we must remove it.
   if (norm_name.length()>2 && norm_name[0]==':' && norm_name[1]==':') {
      norm_name.erase(0,2);
   }

   if (gInterpreterHelper) {
      // See if the expanded name itself is a typedef.
      std::string typeresult;
      if (gInterpreterHelper->ExistingTypeCheck(norm_name, typeresult)
          || gInterpreterHelper->GetPartiallyDesugaredNameWithScopeHandling(norm_name, typeresult)) {

         if (!typeresult.empty()) norm_name = typeresult;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Replace 'long long' and 'unsigned long long' by 'Long64_t' and 'ULong64_t'

string TClassEdit::GetLong64_Name(const char* original)
{
   if (original==0)
      return "";
   else
      return GetLong64_Name(string(original));
}

////////////////////////////////////////////////////////////////////////////////
/// Replace 'long long' and 'unsigned long long' by 'Long64_t' and 'ULong64_t'

string TClassEdit::GetLong64_Name(const string& original)
{
   static const char* longlong_s  = "long long";
   static const char* ulonglong_s = "unsigned long long";
   static const unsigned int longlong_len  = strlen(longlong_s);
   static const unsigned int ulonglong_len = strlen(ulonglong_s);

   string result = original;

   int pos = 0;
   while( (pos = result.find(ulonglong_s,pos) ) >=0 ) {
      result.replace(pos, ulonglong_len, "ULong64_t");
   }
   pos = 0;
   while( (pos = result.find(longlong_s,pos) ) >=0 ) {
      result.replace(pos, longlong_len, "Long64_t");
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the start of the unqualified name include in 'original'.

const char *TClassEdit::GetUnqualifiedName(const char *original)
{
   const char *lastPos = original;
   {
      long depth = 0;
      for(auto cursor = original; *cursor != '\0'; ++cursor) {
         if ( *cursor == '<' || *cursor == '(') ++depth;
         else if ( *cursor == '>' || *cursor == ')' ) --depth;
         else if ( *cursor == ':' ) {
            if (depth==0 && *(cursor+1) == ':' && *(cursor+2) != '\0') {
               lastPos = cursor+2;
            }
         }
      }
   }
   return lastPos;
}

////////////////////////////////////////////////////////////////////////////////

static void R__FindTrailing(std::string &full,  /*modified*/
                            std::string &stars /* the literal output */
                            )
{
   const char *t = full.c_str();
   const unsigned int tlen( full.size() );

   const char *starloc = t + tlen - 1;
   bool hasconst = false;
   if ( (*starloc)=='t'
       && (starloc-t) > 4 && 0 == strncmp((starloc-4),"const",5)
       && ( (*(starloc-5)) == ' ' || (*(starloc-5)) == '*' || (*(starloc-5)) == '&'
           || (*(starloc-5)) == '>' || (*(starloc-5)) == ']') ) {
      // we are ending on a const.
      starloc -= 4;
      if ((*starloc-1)==' ') {
         // Take the space too.
         starloc--;
      }
      hasconst = true;
   }
   if ( hasconst || (*starloc)=='*' || (*starloc)=='&' || (*starloc)==']' ) {
      bool isArray = ( (*starloc)==']' );
      while( t<=(starloc-1) && ((*(starloc-1))=='*' || (*(starloc-1))=='&' || (*(starloc-1))=='t' || isArray)) {
         if (isArray) {
            starloc--;
            isArray = ! ( (*starloc)=='[' );
         } else if ( (*(starloc-1))=='t' ) {
            if ( (starloc-1-t) > 5 && 0 == strncmp((starloc-5),"const",5)
                && ( (*(starloc-6)) == ' ' || (*(starloc-6)) == '*' || (*(starloc-6)) == '&'
                    || (*(starloc-6)) == '>' || (*(starloc-6)) == ']')) {
               // we have a const.
               starloc -= 5;
            } else {
               break;
            }
         } else {
            starloc--;
         }
      }
      stars = starloc;
      if ((*(starloc-1))==' ') {
         // erase the space too.
         starloc--;
      }

      const unsigned int starlen = strlen(starloc);
      full.erase(tlen-starlen,starlen);
   } else if (hasconst) {
      stars = starloc;
      const unsigned int starlen = strlen(starloc);
      full.erase(tlen-starlen,starlen);
   }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///  Stores in output (after emptying it) the split type.
///  Stores the location of the tail (nested names) in nestedLoc (0 indicates no tail).
///  Return the number of elements stored.
///
///  First in list is the template name or is empty
///         "vector<list<int>,alloc>**" to "vector" "list<int>" "alloc" "**"
///   or    "TNamed*" to "" "TNamed" "*"
////////////////////////////////////////////////////////////////////////////

int TClassEdit::GetSplit(const char *type, vector<string>& output, int &nestedLoc, EModType mode)
{
   nestedLoc = 0;
   output.clear();
   if (strlen(type)==0) return 0;

   int cleantypeMode = 1 /* keepInnerConst */;
   if (mode & kKeepOuterConst) {
      cleantypeMode = 0; /* remove only the outer class keyword */
   }
   string full( mode & kLong64 ? TClassEdit::GetLong64_Name( CleanType(type, cleantypeMode) )
               : CleanType(type, cleantypeMode) );

   // We need to replace basic_string with string.
   {
      unsigned int const_offset = (0==strncmp("const ",full.c_str(),6)) ? 6 : 0;
      bool isString = false;
      bool isStdString = false;
      size_t std_offset = const_offset;
      static const char* basic_string_std = "std::basic_string<char";
      static const unsigned int basic_string_std_len = strlen(basic_string_std);

      if (full.compare(const_offset,basic_string_std_len,basic_string_std) == 0
          && full.size() > basic_string_std_len) {
         isString = true;
         isStdString = true;
         std_offset += 5;
      } else if (full.compare(const_offset,basic_string_std_len-5,basic_string_std+5) == 0
                 && full.size() > (basic_string_std_len-5)) {
         // no std.
         isString = true;
      } else if (full.find("basic_string") != std::string::npos) {
         size_t len = StdLen(full.c_str() + const_offset);
         if (len && len != 5 && full.compare(const_offset + len, basic_string_std_len-5, basic_string_std+5) == 0) {
            isString = true;
            isStdString = true;
            std_offset += len;
         }
      }
      if (isString) {
         size_t offset = basic_string_std_len - 5;
         offset += std_offset; // std_offset includs both the size of std prefix and const prefix.
         if ( full[offset] == '>' ) {
            // done.
         } else if (full[offset] == ',') {
            ++offset;
            if (full.compare(offset, 5, "std::") == 0) {
               offset += 5;
            }
            static const char* char_traits_s = "char_traits<char>";
            static const unsigned int char_traits_len = strlen(char_traits_s);
            if (full.compare(offset, char_traits_len, char_traits_s) == 0) {
               offset += char_traits_len;
               if ( full[offset] == '>') {
                  // done.
               } else if (full[offset] == ' ' && full[offset+1] == '>') {
                  ++offset;
                  // done.
               } else if (full[offset] == ',') {
                  ++offset;
                  if (full.compare(offset, 5, "std::") == 0) {
                     offset += 5;
                  }
                  static const char* allocator_s = "allocator<char>";
                  static const unsigned int allocator_len = strlen(allocator_s);
                  if (full.compare(offset, allocator_len, allocator_s) == 0) {
                     offset += allocator_len;
                     if ( full[offset] == '>') {
                        // done.
                     } else if (full[offset] == ' ' && full[offset+1] == '>') {
                        ++offset;
                        // done.
                     } else {
                        // Not std::string
                        isString = false;
                     }
                  }
               } else {
                  // Not std::string
                  isString = false;
               }
            } else {
               // Not std::string.
               isString = false;
            }
         } else {
            // Not std::string.
            isString = false;
         }
         if (isString) {
            output.push_back(string());
            if (const_offset && (mode & kKeepOuterConst)) {
               if (isStdString && !(mode & kDropStd)) {
                  output.push_back("const std::string");
               } else {
                  output.push_back("const string");
               }
            } else {
               if (isStdString && !(mode & kDropStd)) {
                  output.push_back("std::string");
               } else {
                  output.push_back("string");
               }
            }
            if (offset < full.length()) {
               // Copy the trailing text.
               // keep the '>' inside right for R__FindTrailing to work
               string right( full.substr(offset) );
               string stars;
               R__FindTrailing(right, stars);
               output.back().append(right.c_str()+1); // skip the '>'
               output.push_back(stars);
            } else {
               output.push_back("");
            }
            return output.size();
         }
      }
   }

   if ( mode & kDropStd) {
      unsigned int offset = (0==strncmp("const ",full.c_str(),6)) ? 6 : 0;
      RemoveStd( full, offset );
   }

   string stars;
   if ( !full.empty() ) {
      R__FindTrailing(full, stars);
   }

   const char *c = strchr(full.c_str(),'<');
   if (c) {
      //we have 'something<'
      output.push_back(string(full,0,c - full.c_str()));

      const char *cursor;
      int level = 0;
      int parenthesis = 0;
      for(cursor = c + 1; *cursor != '\0' && !(level==0 && *cursor == '>'); ++cursor) {
         if (*cursor == '(') {
            ++parenthesis;
            continue;
         } else if (*cursor == ')') {
            --parenthesis;
            continue;
         }
         if (parenthesis)
            continue;
         switch (*cursor) {
            case '<': ++level; break;
            case '>': --level; break;
            case ',':
               if (level == 0) {
                  output.push_back(std::string(c+1,cursor));
                  c = cursor;
               }
               break;
         }
      }
      if (*cursor=='>') {
         if (*(cursor-1) == ' ') {
            output.push_back(std::string(c+1,cursor-1));
         } else {
            output.push_back(std::string(c+1,cursor));
         }
         // See what's next!
         if (*(cursor+1)==':') {
            // we have a name specified inside the class/namespace
            // For now we keep it in one piece
            nestedLoc = output.size();
            output.push_back((cursor+1));
         }
      } else if (level >= 0) {
         // Unterminated template
         output.push_back(std::string(c+1,cursor));
      }
   } else {
      //empty
      output.push_back(string());
      output.push_back(full);
   }

   if (!output.empty()) output.push_back(stars);
   return output.size();
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///      Cleanup type description, redundant blanks removed
///      and redundant tail ignored
///      return *tail = pointer to last used character
///      if (mode==0) keep keywords
///      if (mode==1) remove keywords outside the template params
///      if (mode>=2) remove the keywords everywhere.
///      if (tail!=0) cut before the trailing *
///
///      The keywords currently are: "const" , "volatile" removed
///
///
///      CleanType(" A<B, C< D, E> > *,F,G>") returns "A<B,C<D,E> >*"
////////////////////////////////////////////////////////////////////////////

string TClassEdit::CleanType(const char *typeDesc, int mode, const char **tail)
{
   static const char* remove[] = {"class","const","volatile",0};
   static bool isinit = false;
   static std::vector<size_t> lengths;
   if (!isinit) {
      for (int k=0; remove[k]; ++k) {
         lengths.push_back(strlen(remove[k]));
      }
      isinit = true;
   }

   string result;
   result.reserve(strlen(typeDesc)*2);
   int lev=0,kbl=1;
   const char* c;

   for(c=typeDesc;*c;c++) {
      if (c[0]==' ') {
         if (kbl)       continue;
         if (!isalnum(c[ 1]) && c[ 1] !='_')    continue;
      }
      if (kbl && (mode>=2 || lev==0)) { //remove "const' etc...
         int done = 0;
         int n = (mode) ? 999 : 1;

         // loop on all the keywords we want to remove
         for (int k=0; k<n && remove[k]; k++) {
            int rlen = lengths[k];

            // Do we have a match
            if (strncmp(remove[k],c,rlen)) continue;

            // make sure that the 'keyword' is not part of a longer indentifier
            if (isalnum(c[rlen]) || c[rlen]=='_' ||  c[rlen]=='$') continue;

            c+=rlen-1; done = 1; break;
         }
         if (done) continue;
      }

      kbl = (!isalnum(c[ 0]) && c[ 0]!='_' && c[ 0]!='$' && c[0]!='[' && c[0]!=']' && c[0]!='-' && c[0]!='@');
      // '@' is special character used only the artifical class name used by ROOT to implement the
      // I/O customization rules that requires caching of the input data.

      if (*c == '<' || *c == '(')   lev++;
      if (lev==0 && !isalnum(*c)) {
         if (!strchr("*&:._$ []-@",*c)) break;
         // '.' is used as a module/namespace separator by PyROOT, see
         // TPyClassGenerator::GetClass.
      }
      if (c[0]=='>' && result.size() && result[result.size()-1]=='>') result+=" ";

      result += c[0];

      if (*c == '>' || *c == ')')    lev--;
   }
   if(tail) *tail=c;
   return result;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/// Return the absolute type of typeDesc.
/// E.g.: typeDesc = "class const volatile TNamed**", returns "TNamed**".
/// if (mode&1) remove last "*"s                     returns "TNamed"
/// if (mode&2) remove default allocators from STL containers
/// if (mode&4) remove all     allocators from STL containers
/// if (mode&8) return inner class of stl container. list<innerClass>
/// if (mode&16) return deapest class of stl container. vector<list<deapest>>
/// if (mode&kDropAllDefault) remove default template arguments
//////////////////////////////////////////////////////////////////////////////

string TClassEdit::ShortType(const char *typeDesc, int mode)
{
   string answer;

   // get list of all arguments
   if (typeDesc) {
      TSplitType arglist(typeDesc, (EModType) mode);
      arglist.ShortType(answer, mode);
   }

   return answer;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the type is one the interpreter details which are
/// only forward declared (ClassInfo_t etc..)

bool TClassEdit::IsInterpreterDetail(const char *type)
{
   size_t len = strlen(type);
   if (len < 2 || strncmp(type+len-2,"_t",2) != 0) return false;

   unsigned char offset = 0;
   if (strncmp(type,"const ",6)==0) { offset += 6; }
   static const char *names[] = { "CallFunc_t","ClassInfo_t","BaseClassInfo_t",
      "DataMemberInfo_t","FuncTempInfo_t","MethodInfo_t","MethodArgInfo_t",
      "TypeInfo_t","TypedefInfo_t",0};

   for(int k=1;names[k];k++) {if (strcmp(type+offset,names[k])==0) return true;}
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true is the name is std::bitset<number> or bitset<number>

bool TClassEdit::IsSTLBitset(const char *classname)
{
   size_t offset = StdLen(classname);
   if ( strncmp(classname+offset,"bitset<",strlen("bitset<"))==0) return true;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the type of STL collection, if any, that is the underlying type
/// of the given type.   Namely return the value of IsSTLCont after stripping
/// pointer, reference and constness from the type.
///    UnderlyingIsSTLCont("vector<int>*") == IsSTLCont("vector<int>")
/// See TClassEdit::IsSTLCont
///
///  type     : type name: vector<list<classA,allocator>,allocator>*
///  result:    0          : not stl container
///             code of container 1=vector,2=list,3=deque,4=map
///                     5=multimap,6=set,7=multiset

ROOT::ESTLType TClassEdit::UnderlyingIsSTLCont(std::string_view type)
{
   if (type.compare(0,6,"const ",6) == 0)
      type.remove_prefix(6);

   while(type[type.length()-1]=='*' ||
         type[type.length()-1]=='&' ||
         type[type.length()-1]==' ') {
      type.remove_suffix(1);
   }
   return IsSTLCont(type);
}

////////////////////////////////////////////////////////////////////////////////
///  type     : type name: vector<list<classA,allocator>,allocator>
///  result:    0          : not stl container
///             code of container 1=vector,2=list,3=deque,4=map
///                     5=multimap,6=set,7=multiset

ROOT::ESTLType TClassEdit::IsSTLCont(std::string_view type)
{
   auto pos = type.find('<');
   if (pos==std::string_view::npos) return ROOT::kNotSTL;

   auto c = pos+1;
   for (decltype(type.length()) level = 1; c < type.length(); ++c) {
      if (type[c] == '<') ++level;
      if (type[c] == '>') --level;
      if (level == 0) break;
   }
   if (c != (type.length()-1) ) {
      return ROOT::kNotSTL;
   }

   return STLKind(type.substr(0,pos));
}

////////////////////////////////////////////////////////////////////////////////
///  type     : type name: vector<list<classA,allocator>,allocator>
///  testAlloc: if true, we test allocator, if it is not default result is negative
///  result:    0          : not stl container
///             abs(result): code of container 1=vector,2=list,3=deque,4=map
///                           5=multimap,6=set,7=multiset
///             positive val: we have a vector or list with default allocator to any depth
///                   like vector<list<vector<int>>>
///             negative val: STL container other than vector or list, or non default allocator
///                           For example: vector<deque<int>> has answer -1

int TClassEdit::IsSTLCont(const char *type, int testAlloc)
{
   if (strchr(type,'<')==0) return 0;

   TSplitType arglist( type );
   return arglist.IsSTLCont(testAlloc);
}

////////////////////////////////////////////////////////////////////////////////
/// return true if the class belongs to the std namespace

bool TClassEdit::IsStdClass(const char *classname)
{
   classname += StdLen( classname );
   if ( strcmp(classname,"string")==0 ) return true;
   if ( strncmp(classname,"bitset<",strlen("bitset<"))==0) return true;
   if ( IsStdPair(classname) ) return true;
   if ( strcmp(classname,"allocator")==0) return true;
   if ( strncmp(classname,"allocator<",strlen("allocator<"))==0) return true;
   if ( strncmp(classname,"greater<",strlen("greater<"))==0) return true;
   if ( strncmp(classname,"less<",strlen("less<"))==0) return true;
   if ( strncmp(classname,"equal_to<",strlen("equal_to<"))==0) return true;
   if ( strncmp(classname,"hash<",strlen("hash<"))==0) return true;
   if ( strncmp(classname,"auto_ptr<",strlen("auto_ptr<"))==0) return true;

   if ( strncmp(classname,"vector<",strlen("vector<"))==0) return true;
   if ( strncmp(classname,"list<",strlen("list<"))==0) return true;
   if ( strncmp(classname,"forward_list<",strlen("forward_list<"))==0) return true;
   if ( strncmp(classname,"deque<",strlen("deque<"))==0) return true;
   if ( strncmp(classname,"map<",strlen("map<"))==0) return true;
   if ( strncmp(classname,"multimap<",strlen("multimap<"))==0) return true;
   if ( strncmp(classname,"set<",strlen("set<"))==0) return true;
   if ( strncmp(classname,"multiset<",strlen("multiset<"))==0) return true;
   if ( strncmp(classname,"unordered_set<",strlen("unordered_set<"))==0) return true;
   if ( strncmp(classname,"unordered_multiset<",strlen("unordered_multiset<"))==0) return true;
   if ( strncmp(classname,"unordered_map<",strlen("unordered_map<"))==0) return true;
   if ( strncmp(classname,"unordered_multimap<",strlen("unordered_multimap<"))==0) return true;
   if ( strncmp(classname,"bitset<",strlen("bitset<"))==0) return true;

   return false;
}


////////////////////////////////////////////////////////////////////////////////

bool TClassEdit::IsVectorBool(const char *name) {
   TSplitType splitname( name );

   return ( TClassEdit::STLKind( splitname.fElements[0] ) == ROOT::kSTLvector)
      && ( splitname.fElements[1] == "bool" || splitname.fElements[1]=="Bool_t");
}

////////////////////////////////////////////////////////////////////////////////

static void ResolveTypedefProcessType(const char *tname,
                                      unsigned int /* len */,
                                      unsigned int cursor,
                                      bool constprefix,
                                      unsigned int start_of_type,
                                      unsigned int end_of_type,
                                      unsigned int mod_start_of_type,
                                      bool &modified,
                                      std::string &result)
{
   std::string type(modified && (mod_start_of_type < result.length()) ?
                    result.substr(mod_start_of_type, string::npos)
                    : string(tname, start_of_type, end_of_type == 0 ? cursor - start_of_type : end_of_type - start_of_type));  // we need to try to avoid this copy
   string typeresult;
   if (gInterpreterHelper->ExistingTypeCheck(type, typeresult)
       || gInterpreterHelper->GetPartiallyDesugaredNameWithScopeHandling(type, typeresult, false)) {
      // it is a known type
      if (!typeresult.empty()) {
         // and it is a typedef, we need to replace it in the output.
         if (modified) {
            result.replace(mod_start_of_type, string::npos,
                           typeresult);
         }
         else {
            modified = true;
            mod_start_of_type = start_of_type;
            result += string(tname,0,start_of_type);
            if (constprefix && typeresult.compare(0,6,"const ",6) == 0) {
               result += typeresult.substr(6,string::npos);
            } else {
               result += typeresult;
            }
         }
      } else if (modified) {
         result.replace(mod_start_of_type, string::npos,
                        type);
      }
      if (modified) {
         if (end_of_type != 0 && end_of_type!=cursor) {
            result += std::string(tname,end_of_type,cursor-end_of_type);
         }
      }
   } else {
      // no change needed.
      if (modified) {
         // result += type;
         if (end_of_type != 0 && end_of_type!=cursor) {
            result += std::string(tname,end_of_type,cursor-end_of_type);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

static void ResolveTypedefImpl(const char *tname,
                               unsigned int len,
                               unsigned int &cursor,
                               bool &modified,
                               std::string &result)
{
   // Need to parse and deal with
   // A::B::C< D, E::F, G::H<I,J>::K::L >::M
   // where E might be replace by N<O,P>
   // and G::H<I,J>::K or G might be a typedef.

   bool constprefix = false;

   if (tname[cursor]==' ') {
      if (!modified) {
         modified = true;
         result += string(tname,0,cursor);
      }
      while (tname[cursor]==' ') ++cursor;
   }

   if (tname[cursor]=='c' && (cursor+6<len)) {
      if (strncmp(tname+cursor,"const ",6) == 0) {
         cursor += 6;
         if (modified) result += "const ";
      }
      constprefix = true;

   }

   if (len > 2 && strncmp(tname+cursor,"::",2) == 0) {
      cursor += 2;
   }

   unsigned int start_of_type = cursor;
   unsigned int end_of_type = 0;
   unsigned int mod_start_of_type = result.length();
   unsigned int prevScope = cursor;
   for ( ; cursor<len; ++cursor) {
      switch (tname[cursor]) {
         case ':': {
            if ((cursor+1)>=len || tname[cursor+1]!=':') {
               // we expected another ':', malformed, give up.
               if (modified) result += (tname+prevScope);
               return;
            }
            string scope;
            if (modified) {
               scope = result.substr(mod_start_of_type, string::npos);
               scope += std::string(tname+prevScope,cursor-prevScope);
            } else {
               scope = std::string(tname, start_of_type, cursor - start_of_type); // we need to try to avoid this copy
            }
            std::string scoperesult;
            bool isInlined = false;
            if (gInterpreterHelper->ExistingTypeCheck(scope, scoperesult)
                ||gInterpreterHelper->GetPartiallyDesugaredNameWithScopeHandling(scope, scoperesult)) {
               // it is a known type
               if (!scoperesult.empty()) {
                  // and it is a typedef
                  if (modified) {
                     if (constprefix && scoperesult.compare(0,6,"const ",6) != 0) mod_start_of_type -= 6;
                     result.replace(mod_start_of_type, string::npos,
                                    scoperesult);
                     result += "::";
                  } else {
                     modified = true;
                     mod_start_of_type = start_of_type;
                     result += string(tname,0,start_of_type);
                     //if (constprefix) result += "const ";
                     result += scoperesult;
                     result += "::";
                  }
               } else if (modified) {
                  result += std::string(tname+prevScope,cursor+2-prevScope);
               }
            } else if (!gInterpreterHelper->IsDeclaredScope(scope,isInlined)) {
               // the nesting namespace is not declared, just ignore it and move on
               if (modified) result += std::string(tname+prevScope,cursor+2-prevScope);
            } else if (isInlined) {
               // humm ... just skip it.
               if (!modified) {
                  modified = true;
                  mod_start_of_type = start_of_type;
                  result += string(tname,0,start_of_type);
                  //if (constprefix) result += "const ";
                  result += string(tname,start_of_type,prevScope - start_of_type);
               }
            } else if (modified) {
               result += std::string(tname+prevScope,cursor+2-prevScope);
            }
            // Consume the 1st semi colon, the 2nd will be consume by the for loop.
            ++cursor;
            prevScope = cursor+1;
            break;
         }
         case '<': {
            // push information on stack
            if (modified) {
               result += std::string(tname+prevScope,cursor+1-prevScope);
               // above includes the '<' .... result += '<';
            }
            do {
               ++cursor;
               ResolveTypedefImpl(tname,len,cursor,modified,result);
            } while( cursor<len && tname[cursor] == ',' );

            while (cursor<len && tname[cursor+1]==' ') ++cursor;

            // Since we already checked the type, skip the next section
            // (respective the scope section and final type processing section)
            // as they would re-do the same job.
            if (cursor+2<len && tname[cursor+1]==':' && tname[cursor+2]==':') {
               if (modified) result += "::";
               cursor += 2;
               prevScope = cursor+1;
            }
            if ( (cursor+1)<len && tname[cursor+1] == ',') {
               ++cursor;
               if (modified) result += ',';
               return;
            }
            if ( (cursor+1)<len && tname[cursor+1] == '>') {
               ++cursor;
               if (modified) result += " >";
               return;
            }
            if ( (cursor+1) >= len) {
               return;
            }
            if (tname[cursor] != ' ') break;
            if (modified) prevScope = cursor+1;
            // If the 'current' character is a space we need to treat it,
            // since this the next case statement, we can just fall through,
            // otherwise we should need to do:
            // --cursor; break;
         }
         case ' ': {
            end_of_type = cursor;
            // let's see if we have 'long long' or 'unsigned int' or 'signed char' or what not.
            while ((cursor+1)<len && tname[cursor+1] == ' ') ++cursor;

            auto next = cursor+1;
            if (strncmp(tname+next,"const",5) == 0 && ((next+5)==len || tname[next+5] == ' ' || tname[next+5] == '*' || tname[next+5] == '&' || tname[next+5] == ',' || tname[next+5] == '>' || tname[next+5] == ']'))
            {
               // A first const after the type needs to be move in the front.
               if (!modified) {
                  modified = true;
                  result += string(tname,0,start_of_type);
                  result += "const ";
                  mod_start_of_type = start_of_type + 6;
                  result += string(tname,start_of_type,end_of_type-start_of_type);
               } else if (mod_start_of_type < result.length()) {
                  result.insert(mod_start_of_type,"const ");
                  mod_start_of_type += 6;
               } else {
                  result += "const ";
                  mod_start_of_type += 6;
                  result += string(tname,start_of_type,end_of_type-start_of_type);
               }
               cursor += 5;
               end_of_type = cursor+1;
               prevScope = end_of_type;
               if ((next+5)==len || tname[next+5] == ',' || tname[next+5] == '>' || tname[next+5] == '[') {
                  break;
               }
            } else if (next!=len && tname[next] != '*' && tname[next] != '&') {
               // the type is not ended yet.
               end_of_type = 0;
               break;
            }
            ++cursor;
            // Intentional fall through;
         }
         case '*':
         case '&': {
            if (tname[cursor] != ' ') end_of_type = cursor;
            // check and skip const (followed by *,&, ,) ... what about followed by ':','['?
            auto next = cursor+1;
            if (strncmp(tname+next,"const",5) == 0) {
               if ((next+5)==len || tname[next+5] == ' ' || tname[next+5] == '*' || tname[next+5] == '&' || tname[next+5] == ',' || tname[next+5] == '>' || tname[next+5] == '[') {
                  next += 5;
               }
            }
            while (next<len &&
                   (tname[next] == ' ' || tname[next] == '*' || tname[next] == '&')) {
               ++next;
               // check and skip const (followed by *,&, ,) ... what about followed by ':','['?
               if (strncmp(tname+next,"const",5) == 0) {
                  if ((next+5)==len || tname[next+5] == ' ' || tname[next+5] == '*' || tname[next+5] == '&' || tname[next+5] == ',' || tname[next+5] == '>' || tname[next+5] == '[') {
                     next += 5;
                  }
               }
            }
            cursor = next-1;
//            if (modified && mod_start_of_type < result.length()) {
//               result += string(tname,end_of_type,cursor-end_of_type);
//            }
            break;
         }
         case ',': {
            if (modified && prevScope) {
               result += std::string(tname+prevScope,(end_of_type == 0 ? cursor : end_of_type)-prevScope);
            }
            ResolveTypedefProcessType(tname,len,cursor,constprefix,start_of_type,end_of_type,mod_start_of_type,
                                      modified, result);
            if (modified) result += ',';
            return;
         }
         case '>': {
            if (modified && prevScope) {
               result += std::string(tname+prevScope,(end_of_type == 0 ? cursor : end_of_type)-prevScope);
            }
            ResolveTypedefProcessType(tname,len,cursor,constprefix,start_of_type,end_of_type,mod_start_of_type,
                                      modified, result);
            if (modified) result += '>';
            return;
         }
         default:
            end_of_type = 0;
      }
   }

   if (prevScope && modified) result += std::string(tname+prevScope,(end_of_type == 0 ? cursor : end_of_type)-prevScope);

   ResolveTypedefProcessType(tname,len,cursor,constprefix,start_of_type,end_of_type,mod_start_of_type,
                             modified, result);
}


////////////////////////////////////////////////////////////////////////////////

string TClassEdit::ResolveTypedef(const char *tname, bool /* resolveAll */)
{
   // Return the name of type 'tname' with all its typedef components replaced
   // by the actual type its points to
   // For example for "typedef MyObj MyObjTypedef;"
   //    vector<MyObjTypedef> return vector<MyObj>
   //

   if (tname == 0 || tname[0] == 0)
      return "";
   if (!gInterpreterHelper)
      return tname;

   std::string result;

   // Check if we already know it is a normalized typename or a registered
   // typedef (i.e. known to gROOT).
   if (gInterpreterHelper->ExistingTypeCheck(tname, result))
   {
      if (result.empty()) return tname;
      else return result;
   }

   unsigned int len = strlen(tname);

   unsigned int cursor = 0;
   bool modified = false;
   ResolveTypedefImpl(tname,len,cursor,modified,result);

   if (!modified) return tname;
   else return result;
}


////////////////////////////////////////////////////////////////////////////////

string TClassEdit::InsertStd(const char *tname)
{
   // Return the name of type 'tname' with all STL classes prepended by "std::".
   // For example for "vector<set<auto_ptr<int*> > >" it returns
   //    "std::vector<std::set<std::auto_ptr<int*> > >"
   //

   static const char* sSTLtypes[] = {
      "allocator",
      "auto_ptr",
      "bad_alloc",
      "bad_cast",
      "bad_exception",
      "bad_typeid",
      "basic_filebuf",
      "basic_fstream",
      "basic_ifstream",
      "basic_ios",
      "basic_iostream",
      "basic_istream",
      "basic_istringstream",
      "basic_ofstream",
      "basic_ostream",
      "basic_ostringstream",
      "basic_streambuf",
      "basic_string",
      "basic_stringbuf",
      "basic_stringstream",
      "binary_function",
      "binary_negate",
      "bitset",
      "char_traits",
      "codecvt_byname",
      "codecvt",
      "collate",
      "collate_byname",
      "compare",
      "complex",
      "ctype_byname",
      "ctype",
      "deque",
      "divides",
      "domain_error",
      "equal_to",
      "exception",
      "forward_list",
      "fpos",
      "greater_equal",
      "greater",
      "gslice_array",
      "gslice",
      "hash",
      "indirect_array",
      "invalid_argument",
      "ios_base",
      "istream_iterator",
      "istreambuf_iterator",
      "istrstream",
      "iterator_traits",
      "iterator",
      "length_error",
      "less_equal",
      "less",
      "list",
      "locale",
      "localedef utility",
      "locale utility",
      "logic_error",
      "logical_and",
      "logical_not",
      "logical_or",
      "map",
      "mask_array",
      "mem_fun",
      "mem_fun_ref",
      "messages",
      "messages_byname",
      "minus",
      "modulus",
      "money_get",
      "money_put",
      "moneypunct",
      "moneypunct_byname",
      "multimap",
      "multiplies",
      "multiset",
      "negate",
      "not_equal_to",
      "num_get",
      "num_put",
      "numeric_limits",
      "numpunct",
      "numpunct_byname",
      "ostream_iterator",
      "ostreambuf_iterator",
      "ostrstream",
      "out_of_range",
      "overflow_error",
      "pair",
      "plus",
      "pointer_to_binary_function",
      "pointer_to_unary_function",
      "priority_queue",
      "queue",
      "range_error",
      "raw_storage_iterator",
      "reverse_iterator",
      "runtime_error",
      "set",
      "slice_array",
      "slice",
      "stack",
      "string",
      "strstream",
      "strstreambuf",
      "time_get_byname",
      "time_get",
      "time_put_byname",
      "time_put",
      "unary_function",
      "unary_negate",
      "unique_pointer",
      "underflow_error",
      "unordered_map",
      "unordered_multimap",
      "unordered_multiset",
      "unordered_set",
      "valarray",
      "vector",
      "wstring"
   };
   static ShuttingDownSignaler<set<string>> sSetSTLtypes;

   if (tname==0 || tname[0]==0) return "";

   if (sSetSTLtypes.empty()) {
      // set up static set
      const size_t nSTLtypes = sizeof(sSTLtypes) / sizeof(const char*);
      for (size_t i = 0; i < nSTLtypes; ++i)
         sSetSTLtypes.insert(sSTLtypes[i]);
   }

   size_t b = 0;
   size_t len = strlen(tname);
   string ret;
   ret.reserve(len + 20); // expect up to 4 extra "std::" to insert
   string id;
   while (b < len) {
      // find beginning of next identifier
      bool precScope = false; // whether the identifier was preceded by "::"
      while (!(isalnum(tname[b]) || tname[b] == '_') && b < len) {
         precScope = (b < len - 2) && (tname[b] == ':') && (tname[b + 1] == ':');
         if (precScope) {
            ret += "::";
            b += 2;
         } else
            ret += tname[b++];
      }

      // now b is at the beginning of an identifier or len
      size_t e = b;
      // find end of identifier
      id.clear();
      while (e < len && (isalnum(tname[e]) || tname[e] == '_'))
         id += tname[e++];
      if (!id.empty()) {
         if (!precScope) {
            set<string>::const_iterator iSTLtype = sSetSTLtypes.find(id);
            if (iSTLtype != sSetSTLtypes.end())
               ret += "std::";
         }

         ret += id;
         b = e;
      }
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// An helper class to dismount the name and remount it changed whenever
/// necessary

class NameCleanerForIO {
   std::string fName;
   std::vector<std::unique_ptr<NameCleanerForIO>> fArgumentNodes = {};
   NameCleanerForIO* fMother;
   bool fHasChanged = false;
   bool AreAncestorsSTLContOrArray()
   {
      NameCleanerForIO* mother = fMother;
      if (!mother) return false;
      bool isSTLContOrArray = true;
      while (nullptr != mother){
         auto stlType = TClassEdit::IsSTLCont(mother->fName+"<>");
         isSTLContOrArray &= ROOT::kNotSTL != stlType || TClassEdit::IsStdArray(mother->fName+"<");
         mother = mother->fMother;
      }

      return isSTLContOrArray;
   }

public:
   NameCleanerForIO(const std::string& clName = "",
                    TClassEdit::EModType mode = TClassEdit::kNone,
                    NameCleanerForIO* mother = nullptr):fMother(mother)
   {
      if (clName.back() != '>') {
         fName = clName;
         return;
      }

      std::vector<std::string> v;
      int dummy=0;
      TClassEdit::GetSplit(clName.c_str(), v, dummy, mode);

      // We could be in presence of templates such as A1<T1>::A2<T2>::A3<T3>
      auto argsEnd = v.end();
      auto argsBeginPlusOne = ++v.begin();
      auto argPos = std::find_if(argsBeginPlusOne, argsEnd,
           [](std::string& arg){return (!arg.empty() && arg.front() == ':');});
      if (argPos != argsEnd) {
         const int lenght = clName.size();
         int wedgeBalance = 0;
         int lastOpenWedge = 0;
         for (int i=lenght-1;i>-1;i--) {
            auto& c = clName.at(i);
            if (c == '<') {
               wedgeBalance++;
               lastOpenWedge = i;
            } else if (c == '>') {
               wedgeBalance--;
            } else if (c == ':' && 0 == wedgeBalance) {
               // This would be A1<T1>::A2<T2>
               auto nameToClean = clName.substr(0,i-1);
               NameCleanerForIO node(nameToClean, mode);
               auto cleanName = node.ToString();
               fHasChanged = node.HasChanged();
               // We got A1<T1>::A2<T2> cleaned

               // We build the changed A1<T1>::A2<T2>::A3
               cleanName += "::";
               // Now we get A3 and append it
               cleanName += clName.substr(i+1,lastOpenWedge-i-1);

               // We now get the args of what in our case is A1<T1>::A2<T2>::A3
               auto lastTemplate = &clName.data()[i+1];

               // We split it
               TClassEdit::GetSplit(lastTemplate, v, dummy, mode);
               // We now replace the name of the template
               v[0] = cleanName;
               break;
            }
         }
      }

      fName = v.front();
      unsigned int nargs = v.size() - 2;
      for (unsigned int i=0;i<nargs;++i) {
         fArgumentNodes.emplace_back(new NameCleanerForIO(v[i+1],mode,this));
      }
   }

   bool HasChanged() const {return fHasChanged;}

   std::string ToString()
   {
      std::string name(fName);

      if (fArgumentNodes.empty()) return name;

      // We have in hands a case like unique_ptr< ... >
      // Perhaps we could treat atomics as well like this?
      if (!fMother && TClassEdit::IsUniquePtr(fName+"<")) {
         name = fArgumentNodes.front()->ToString();
         // ROOT-9933: we remove const if present.
         TClassEdit::TSplitType tst(name.c_str());
         tst.ShortType(name, 1);
         fHasChanged = true;
         return name;
      }

      // Now we treat the case of the collections of unique ptrs
      auto stlContType = AreAncestorsSTLContOrArray();
      if (stlContType != ROOT::kNotSTL && TClassEdit::IsUniquePtr(fName+"<")) {
         name = fArgumentNodes.front()->ToString();
         name += "*";
         fHasChanged = true;
         return name;
      }

      name += "<";
      for (auto& node : fArgumentNodes) {
         name += node->ToString() + ",";
         fHasChanged |= node->HasChanged();
      }
      name.pop_back(); // Remove the last comma.
      name += name.back() == '>' ? " >" : ">"; // Respect name normalisation
      return name;
   }

   const std::string& GetName() {return fName;}
   const std::vector<std::unique_ptr<NameCleanerForIO>>* GetChildNodes() const {return &fArgumentNodes;}
};

////////////////////////////////////////////////////////////////////////////////

std::string TClassEdit::GetNameForIO(const std::string& templateInstanceName,
                                     TClassEdit::EModType mode,
                                     bool* hasChanged)
{
   // Decompose template name into pieces and remount it applying the necessary
   // transformations necessary for the ROOT IO subsystem, namely:
   // - Transform std::unique_ptr<T> into T (for selections) (also nested)
   // - Transform std::unique_ptr<const T> into T (for selections) (also nested)
   // - Transform std::COLL<std::unique_ptr<T>> into std::COLL<T*> (also nested)
   // Name normalisation is respected (e.g. spaces).
   // The implementation uses an internal class defined in the cxx file.
   NameCleanerForIO node(templateInstanceName, mode);
   auto nameForIO = node.ToString();
   if (hasChanged) {
      *hasChanged = node.HasChanged();
      }
   return nameForIO;
}

////////////////////////////////////////////////////////////////////////////////
// We could introduce a tuple as return type, but we be consistent with the rest
// of the code.
bool TClassEdit::GetStdArrayProperties(const char* typeName,
                                       std::string& typeNameBuf,
                                       std::array<int, 5>& maxIndices,
                                       int& ndim)
{
   if (!IsStdArray(typeName)) return false;

   // We have an array, it's worth continuing
   NameCleanerForIO node(typeName);

   // We now recurse updating the data according to what we find
   auto childNodes = node.GetChildNodes();
   for (ndim = 1;ndim <=5 ; ndim++) {
      maxIndices[ndim-1] = std::atoi(childNodes->back()->GetName().c_str());
      auto& frontNode = childNodes->front();
      typeNameBuf = frontNode->GetName();
      if (! IsStdArray(typeNameBuf+"<")) {
         typeNameBuf = frontNode->ToString();
         return true;
      }
      childNodes = frontNode->GetChildNodes();
   }

   return true;
}


////////////////////////////////////////////////////////////////////////////////
/// Demangle in a portable way the type id name.
/// IMPORTANT: The caller is responsible for freeing the returned const char*

char* TClassEdit::DemangleTypeIdName(const std::type_info& ti, int& errorCode)
{
   const char* mangled_name = ti.name();
   return DemangleName(mangled_name, errorCode);
}
/*
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
*/

namespace {
   /// Find the first occurrence of any of needle's characters in haystack that
   /// is not nested in a <>, () or [] pair.
   std::size_t FindNonNestedNeedles(std::string_view haystack, string_view needles)
   {
      std::stack<char> expected;
      for (std::size_t pos = 0, end = haystack.length(); pos < end; ++pos) {
         char c = haystack[pos];
         if (expected.empty()) {
            if (needles.find(c) != std::string_view::npos)
               return pos;
         } else {
            if (c == expected.top()) {
               expected.pop();
               continue;
            }
         }
         switch (c) {
            case '<': expected.emplace('>'); break;
            case '(': expected.emplace(')'); break;
            case '[': expected.emplace(']'); break;
         }
      }
      return std::string_view::npos;
   }

   /// Find the first occurrence of `::` that is not nested in a <>, () or [] pair.
   std::size_t FindNonNestedDoubleColons(std::string_view haystack)
   {
      std::size_t lenHaystack = haystack.length();
      std::size_t prevAfterColumn = 0;
      while (true) {
         std::size_t posColumn = FindNonNestedNeedles(haystack.substr(prevAfterColumn), ":");
         if (posColumn == std::string_view::npos)
            return std::string_view::npos;
         prevAfterColumn += posColumn;
         // prevAfterColumn must have "::", i.e. two characters:
         if (prevAfterColumn + 1 >= lenHaystack)
            return std::string_view::npos;

         ++prevAfterColumn; // done with first (or only) ':'
         if (haystack[prevAfterColumn] == ':')
            return prevAfterColumn - 1;
         ++prevAfterColumn; // That was not a ':'.
      }

      return std::string_view::npos;
   }

   std::string_view StripSurroundingSpace(std::string_view str)
   {
      while (!str.empty() && std::isspace(str[0]))
         str.remove_prefix(1);
      while (!str.empty() && std::isspace(str.back()))
         str.remove_suffix(1);
      return str;
   }

   std::string ToString(std::string_view sv)
   {
      // ROOT's string_view backport doesn't add the new std::string contructor and assignment;
      // convert to std::string instead and assign that.
      return std::string(sv.data(), sv.length());
   }
} // unnamed namespace

/// Split a function declaration into its different parts.
bool TClassEdit::SplitFunction(std::string_view decl, TClassEdit::FunctionSplitInfo &result)
{
   // General structure:
   // `...` last-space `...` (`...`)
   // The first `...` is the return type.
   // The second `...` is the (possibly scoped) function name.
   // The third `...` are the parameters.
   // The function name can be of the form `...`<`...`>
   std::size_t posArgs = FindNonNestedNeedles(decl, "(");
   std::string_view declNoArgs = decl.substr(0, posArgs);

   std::size_t prevAfterWhiteSpace = 0;
   static const char whitespace[] = " \t\n";
   while (declNoArgs.length() > prevAfterWhiteSpace) {
      std::size_t posWS = FindNonNestedNeedles(declNoArgs.substr(prevAfterWhiteSpace), whitespace);
      if (posWS == std::string_view::npos)
         break;
      prevAfterWhiteSpace += posWS + 1;
      while (declNoArgs.length() > prevAfterWhiteSpace
             && strchr(whitespace, declNoArgs[prevAfterWhiteSpace]))
         ++prevAfterWhiteSpace;
   }

   /// Include any '&*' in the return type:
   std::size_t endReturn = prevAfterWhiteSpace;
   while (declNoArgs.length() > endReturn
          && strchr("&* \t \n", declNoArgs[endReturn]))
          ++endReturn;

   result.fReturnType = ToString(StripSurroundingSpace(declNoArgs.substr(0, endReturn)));

   /// scope::anotherscope::functionName<tmplt>:
   std::string_view scopeFunctionTmplt = declNoArgs.substr(endReturn);
   std::size_t prevAtScope = FindNonNestedDoubleColons(scopeFunctionTmplt);
   while (prevAtScope != std::string_view::npos
          && scopeFunctionTmplt.length() > prevAtScope + 2) {
      std::size_t posScope = FindNonNestedDoubleColons(scopeFunctionTmplt.substr(prevAtScope + 2));
      if (posScope == std::string_view::npos)
         break;
      prevAtScope += posScope + 2;
   }

   std::size_t afterScope = prevAtScope + 2;
   if (prevAtScope == std::string_view::npos) {
      afterScope = 0;
      prevAtScope = 0;
   }

   result.fScopeName = ToString(StripSurroundingSpace(scopeFunctionTmplt.substr(0, prevAtScope)));
   std::string_view funcNameTmplArgs = scopeFunctionTmplt.substr(afterScope);

   result.fFunctionTemplateArguments.clear();
   std::size_t posTmpltOpen = FindNonNestedNeedles(funcNameTmplArgs, "<");
   if (posTmpltOpen != std::string_view::npos) {
      result.fFunctionName = ToString(StripSurroundingSpace(funcNameTmplArgs.substr(0, posTmpltOpen)));

      // Parse template parameters:
      std::string_view tmpltArgs = funcNameTmplArgs.substr(posTmpltOpen + 1);
      std::size_t posTmpltClose = FindNonNestedNeedles(tmpltArgs, ">");
      if (posTmpltClose != std::string_view::npos) {
         tmpltArgs = tmpltArgs.substr(0, posTmpltClose);
         std::size_t prevAfterArg = 0;
         while (tmpltArgs.length() > prevAfterArg) {
            std::size_t posComma = FindNonNestedNeedles(tmpltArgs.substr(prevAfterArg), ",");
            if (posComma == std::string_view::npos) {
               break;
            }
            result.fFunctionTemplateArguments.emplace_back(ToString(StripSurroundingSpace(tmpltArgs.substr(prevAfterArg, posComma))));
            prevAfterArg += posComma + 1;
         }
         // Add the trailing arg.
         result.fFunctionTemplateArguments.emplace_back(ToString(StripSurroundingSpace(tmpltArgs.substr(prevAfterArg))));
      }
   } else {
      result.fFunctionName = ToString(StripSurroundingSpace(funcNameTmplArgs));
   }

   result.fFunctionParameters.clear();
   if (posArgs != std::string_view::npos) {
      /// (params)
      std::string_view params = decl.substr(posArgs + 1);
      std::size_t posEndArgs = FindNonNestedNeedles(params, ")");
      if (posEndArgs != std::string_view::npos) {
         params = params.substr(0, posEndArgs);
         std::size_t prevAfterArg = 0;
         while (params.length() > prevAfterArg) {
            std::size_t posComma = FindNonNestedNeedles(params.substr(prevAfterArg), ",");
            if (posComma == std::string_view::npos) {
               result.fFunctionParameters.emplace_back(ToString(StripSurroundingSpace(params.substr(prevAfterArg))));
               break;
            }
            result.fFunctionParameters.emplace_back(ToString(StripSurroundingSpace(params.substr(prevAfterArg, posComma))));
            prevAfterArg += posComma + 1; // skip ','
         }
      }
   }

   return true;
}
