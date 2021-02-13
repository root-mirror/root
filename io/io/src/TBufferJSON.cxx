//
// Author: Sergey Linev  4.03.2014

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TBufferJSON
\ingroup IO

Class for serializing object to and from JavaScript Object Notation (JSON) format.
It creates such object representation, which can be directly
used in JavaScript ROOT (JSROOT) for drawing.

TBufferJSON implements TBuffer interface, therefore most of
ROOT and user classes can be converted into JSON.
There are certain limitations for classes with custom streamers,
which should be equipped specially for this purposes (see TCanvas::Streamer()
as example).

To perform conversion into JSON, one should use TBufferJSON::ToJSON method:
~~~{.cpp}
   TH1 *h1 = new TH1I("h1", "title", 100, 0, 10);
   h1->FillRandom("gaus",10000);
   TString json = TBufferJSON::ToJSON(h1);
~~~

To reconstruct object from the JSON string, one should do:
~~~{.cpp}
   TH1 *hnew = nullptr;
   TBufferJSON::FromJSON(hnew, json);
   if (hnew) hnew->Draw("hist");
~~~
JSON does not include stored class version, therefore schema evolution
(reading of older class versions) is not supported. JSON should not be used as
persistent storage for object data - only for live applications.

All STL containers by default converted into JSON Array. Vector of integers:
~~~{.cpp}
   std::vector<int> vect = {1,4,7};
   auto json = TBufferJSON::ToJSON(&vect);
~~~
Will produce JSON code "[1, 4, 7]".

There are special handling for map classes like `map` and `multimap`.
They will create Array of pair objects with "first" and "second" as data members. Code:
~~~{.cpp}
   std::map<int,string> m;
   m[1] = "number 1";
   m[2] = "number 2";
   auto json = TBufferJSON::ToJSON(&m);
~~~
Will generate json string:
~~~{.json}
[
  {"$pair" : "pair<int,string>", "first" : 1, "second" : "number 1"},
  {"$pair" : "pair<int,string>", "first" : 2, "second" : "number 2"}
]
~~~
In special cases map container can be converted into JSON object. For that key parameter
must be `std::string` and compact parameter should be 5.
Like in example:
~~~{.cpp}
std::map<std::string,int> data;
data["name1"] = 11;
data["name2"] = 22;

auto json = TBufferJSON::ToJSON(&data, TBufferJSON::kMapAsObject);
~~~
Will produce JSON output:
~~~
{
  "_typename": "map<string,int>",
  "name1": 11,
  "name2": 22
}
~~~
Another possibility to enforce such conversion - add "JSON_object" into comment line of correspondent
data member like:
~~~{.cpp}
class Container {
   std::map<std::string,int> data;  ///<  JSON_object
};
~~~

*/

#include "TBufferJSON.h"

#include <typeinfo>
#include <string>
#include <cstring>
#include <locale.h>
#include <cmath>
#include <memory>
#include <cstdlib>
#include <fstream>

#include <ROOT/RMakeUnique.hxx>

#include "Compression.h"

#include "TArrayI.h"
#include "TError.h"
#include "TBase64.h"
#include "TROOT.h"
#include "TList.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TClassEdit.h"
#include "TDataType.h"
#include "TRealData.h"
#include "TDataMember.h"
#include "TMap.h"
#include "TRef.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TMemberStreamer.h"
#include "TStreamer.h"
#include "RZip.h"
#include "TClonesArray.h"
#include "TVirtualMutex.h"
#include "TInterpreter.h"
#include "TEmulatedCollectionProxy.h"
#include "snprintf.h"

#include <nlohmann/json.hpp>

ClassImp(TBufferJSON);

enum { json_TArray = 100, json_TCollection = -130, json_TString = 110, json_stdstring = 120 };

///////////////////////////////////////////////////////////////
// TArrayIndexProducer is used to correctly create
/// JSON array separators for multi-dimensional JSON arrays
/// It fully reproduces array dimensions as in original ROOT classes
/// Contrary to binary I/O, which always writes flat arrays

class TArrayIndexProducer {
protected:
   Int_t fTotalLen{0};
   Int_t fCnt{-1};
   const char *fSepar{nullptr};
   TArrayI fIndicies;
   TArrayI fMaxIndex;
   TString fRes;
   Bool_t fIsArray{kFALSE};

public:
   TArrayIndexProducer(TStreamerElement *elem, Int_t arraylen, const char *separ) : fSepar(separ)
   {
      Bool_t usearrayindx = elem && (elem->GetArrayDim() > 0);
      Bool_t isloop = elem && ((elem->GetType() == TStreamerInfo::kStreamLoop) ||
                               (elem->GetType() == TStreamerInfo::kOffsetL + TStreamerInfo::kStreamLoop));
      Bool_t usearraylen = (arraylen > (isloop ? 0 : 1));

      if (usearrayindx && (arraylen > 0)) {
         if (isloop) {
            usearrayindx = kFALSE;
            usearraylen = kTRUE;
         } else if (arraylen != elem->GetArrayLength()) {
            ::Error("TArrayIndexProducer", "Problem with JSON coding of element %s type %d", elem->GetName(),
                    elem->GetType());
         }
      }

      if (usearrayindx) {
         fTotalLen = elem->GetArrayLength();
         fMaxIndex.Set(elem->GetArrayDim());
         for (int dim = 0; dim < elem->GetArrayDim(); dim++)
            fMaxIndex[dim] = elem->GetMaxIndex(dim);
         fIsArray = fTotalLen > 1;
      } else if (usearraylen) {
         fTotalLen = arraylen;
         fMaxIndex.Set(1);
         fMaxIndex[0] = arraylen;
         fIsArray = kTRUE;
      }

      if (fMaxIndex.GetSize() > 0) {
         fIndicies.Set(fMaxIndex.GetSize());
         fIndicies.Reset(0);
      }
   }

   TArrayIndexProducer(TDataMember *member, Int_t extradim, const char *separ) : fSepar(separ)
   {
      Int_t ndim = member->GetArrayDim();
      if (extradim > 0)
         ndim++;

      if (ndim > 0) {
         fIndicies.Set(ndim);
         fIndicies.Reset(0);
         fMaxIndex.Set(ndim);
         fTotalLen = 1;
         for (int dim = 0; dim < member->GetArrayDim(); dim++) {
            fMaxIndex[dim] = member->GetMaxIndex(dim);
            fTotalLen *= member->GetMaxIndex(dim);
         }

         if (extradim > 0) {
            fMaxIndex[ndim - 1] = extradim;
            fTotalLen *= extradim;
         }
      }
      fIsArray = fTotalLen > 1;
   }

   /// returns number of array dimensions
   Int_t NumDimensions() const { return fIndicies.GetSize(); }

   /// return array with current index
   TArrayI &GetIndices() { return fIndicies; };

   /// returns total number of elements in array
   Int_t TotalLength() const { return fTotalLen; }

   Int_t ReduceDimension()
   {
      // reduce one dimension of the array
      // return size of reduced dimension
      if (fMaxIndex.GetSize() == 0)
         return 0;
      Int_t ndim = fMaxIndex.GetSize() - 1;
      Int_t len = fMaxIndex[ndim];
      fMaxIndex.Set(ndim);
      fIndicies.Set(ndim);
      fTotalLen = fTotalLen / len;
      fIsArray = fTotalLen > 1;
      return len;
   }

   Bool_t IsArray() const { return fIsArray; }

   Bool_t IsDone() const
   {
      // return true when iteration over all arrays indexes are done
      return !IsArray() || (fCnt >= fTotalLen);
   }

   const char *GetBegin()
   {
      ++fCnt;
      // return starting separator
      fRes.Clear();
      for (Int_t n = 0; n < fIndicies.GetSize(); ++n)
         fRes.Append("[");
      return fRes.Data();
   }

   const char *GetEnd()
   {
      // return ending separator
      fRes.Clear();
      for (Int_t n = 0; n < fIndicies.GetSize(); ++n)
         fRes.Append("]");
      return fRes.Data();
   }

   /// increment indexes and returns intermediate or last separator
   const char *NextSeparator()
   {
      if (++fCnt >= fTotalLen)
         return GetEnd();

      Int_t cnt = fIndicies.GetSize() - 1;
      fIndicies[cnt]++;

      fRes.Clear();

      while ((cnt >= 0) && (cnt < fIndicies.GetSize())) {
         if (fIndicies[cnt] >= fMaxIndex[cnt]) {
            fRes.Append("]");
            fIndicies[cnt--] = 0;
            if (cnt >= 0)
               fIndicies[cnt]++;
            continue;
         }
         fRes.Append(fIndicies[cnt] == 0 ? "[" : fSepar);
         cnt++;
      }
      return fRes.Data();
   }

   nlohmann::json *ExtractNode(nlohmann::json *topnode, bool next = true)
   {
      if (!IsArray())
         return topnode;
      nlohmann::json *subnode = &((*((nlohmann::json *)topnode))[fIndicies[0]]);
      for (int k = 1; k < fIndicies.GetSize(); ++k)
         subnode = &((*subnode)[fIndicies[k]]);
      if (next)
         NextSeparator();
      return subnode;
   }
};

// TJSONStackObj is used to keep stack of object hierarchy,
// stored in TBuffer. For instance, data for parent class(es)
// stored in subnodes, but initial object node will be kept.

class TJSONStackObj : public TObject {
   struct StlRead {
      Int_t fIndx{0};                   //! index of object in STL container
      Int_t fMap{0};                    //! special iterator over STL map::key members
      Bool_t fFirst{kTRUE};             //! is first or second element is used in the pair
      nlohmann::json::iterator fIter;   //! iterator for std::map stored as JSON object
      const char *fTypeTag{nullptr};    //! type tag used for std::map stored as JSON object
      nlohmann::json fValue;            //! temporary value reading std::map as JSON
      nlohmann::json *GetStlNode(nlohmann::json *prnt)
      {
         if (fMap <= 0)
            return &(prnt->at(fIndx++));

         if (fMap == 1) {
            nlohmann::json *json = &(prnt->at(fIndx));
            if (!fFirst) fIndx++;
            json = &(json->at(fFirst ? "first" : "second"));
            fFirst = !fFirst;
            return  json;
         }

         if (fIndx == 0) {
            // skip _typename if appears
            if (fTypeTag && (fIter.key().compare(fTypeTag) == 0))
               ++fIter;
            fValue = fIter.key();
            fIndx++;
         } else {
            fValue = fIter.value();
             ++fIter;
            fIndx = 0;
         }
         return &fValue;
      }
   };

public:
   TStreamerInfo *fInfo{nullptr};       //!
   TStreamerElement *fElem{nullptr};    //! element in streamer info
   Bool_t fIsStreamerInfo{kFALSE};      //!
   Bool_t fIsElemOwner{kFALSE};         //!
   Bool_t fIsPostProcessed{kFALSE};     //! indicate that value is written
   Bool_t fIsObjStarted{kFALSE};        //! indicate that object writing started, should be closed in postprocess
   Bool_t fAccObjects{kFALSE};          //! if true, accumulate whole objects in values
   Bool_t fBase64{kFALSE};              //! enable base64 coding when writing array
   std::vector<std::string> fValues;    //! raw values
   int fMemberCnt{1};                   //! count number of object members, normally _typename is first member
   int *fMemberPtr{nullptr};            //! pointer on members counter, can be inherit from parent stack objects
   Int_t fLevel{0};                     //! indent level
   std::unique_ptr<TArrayIndexProducer> fIndx; //! producer of ndim indexes
   nlohmann::json *fNode{nullptr};      //! JSON node, used for reading
   std::unique_ptr<StlRead> fStlRead;   //! custom structure for stl container reading
   Version_t fClVersion{0};             //! keep actual class version, workaround for ReadVersion in custom streamer

   TJSONStackObj() = default;

   ~TJSONStackObj()
   {
      if (fIsElemOwner)
         delete fElem;
   }

   Bool_t IsStreamerInfo() const { return fIsStreamerInfo; }

   Bool_t IsStreamerElement() const { return !fIsStreamerInfo && fElem; }

   void PushValue(TString &v)
   {
      fValues.emplace_back(v.Data());
      v.Clear();
   }

   void PushIntValue(Int_t v) { fValues.emplace_back(std::to_string(v)); }

   ////////////////////////////////////////////////////////////////////////
   /// returns separator for data members
   const char *NextMemberSeparator()
   {
      return (!fMemberPtr || ((*fMemberPtr)++ > 0)) ? ","  : "";
   }

   Bool_t IsJsonString() { return fNode && fNode->is_string(); }

   ////////////////////////////////////////////////////////////////////////
   /// checks if specified JSON node is array (compressed or not compressed)
   /// returns length of array (or -1 if failure)
   Int_t IsJsonArray(nlohmann::json *json = nullptr, const char *map_convert_type = nullptr)
   {
      if (!json)
         json = fNode;

      if (map_convert_type) {
         if (!json->is_object()) return -1;
         int sz = 0;
         // count size of object, excluding _typename tag
         for (auto it = json->begin(); it != json->end(); ++it) {
            if ((strlen(map_convert_type)==0) || (it.key().compare(map_convert_type) != 0)) sz++;
         }
         return sz;
      }

      // normal uncompressed array
      if (json->is_array())
         return json->size();

      // compressed array, full array length in "len" attribute, only ReadFastArray
      if (json->is_object() && (json->count("$arr") == 1))
         return json->at("len").get<int>();

      return -1;
   }

   Int_t PopIntValue()
   {
      auto res = std::stoi(fValues.back());
      fValues.pop_back();
      return res;
   }

   std::unique_ptr<TArrayIndexProducer> MakeReadIndexes()
   {
      if (!fElem || (fElem->GetType() <= TStreamerInfo::kOffsetL) ||
          (fElem->GetType() >= TStreamerInfo::kOffsetL + 20) || (fElem->GetArrayDim() < 2))
         return nullptr;

      auto indx = std::make_unique<TArrayIndexProducer>(fElem, -1, "");

      // no need for single dimension - it can be handled directly
      if (!indx->IsArray() || (indx->NumDimensions() < 2))
         return nullptr;

      return indx;
   }

   Bool_t IsStl() const { return fStlRead.get() != nullptr; }

   Bool_t AssignStl(TClass *cl, Int_t map_convert, const char *typename_tag)
   {
      fStlRead = std::make_unique<StlRead>();
      fStlRead->fMap = map_convert;
      if (map_convert == 2) {
         if (!fNode->is_object()) {
            ::Error("TJSONStackObj::AssignStl", "when reading %s expecting JSON object", cl->GetName());
            return kFALSE;
         }
         fStlRead->fIter = fNode->begin();
         fStlRead->fTypeTag = typename_tag && (strlen(typename_tag) > 0) ? typename_tag : nullptr;
      } else {
         if (!fNode->is_array() && !(fNode->is_object() && (fNode->count("$arr") == 1))) {
            ::Error("TJSONStackObj::AssignStl", "when reading %s expecting JSON array", cl->GetName());
            return kFALSE;
         }
      }
      return kTRUE;
   }

   nlohmann::json *GetStlNode()
   {
      return fStlRead ? fStlRead->GetStlNode(fNode) : fNode;
   }

   void ClearStl()
   {
      fStlRead.reset(nullptr);
   }
};

////////////////////////////////////////////////////////////////////////////////
/// Creates buffer object to serialize data into json.

TBufferJSON::TBufferJSON(TBuffer::EMode mode)
   : TBufferText(mode), fOutBuffer(), fOutput(nullptr), fValue(), fStack(), fSemicolon(" : "), fArraySepar(", "),
     fNumericLocale(), fTypeNameTag("_typename")
{
   fOutBuffer.Capacity(10000);
   fValue.Capacity(1000);
   fOutput = &fOutBuffer;

   // checks if setlocale(LC_NUMERIC) returns others than "C"
   // in this case locale will be changed and restored at the end of object conversion

   char *loc = setlocale(LC_NUMERIC, nullptr);
   if (loc && (strcmp(loc, "C") != 0)) {
      fNumericLocale = loc;
      setlocale(LC_NUMERIC, "C");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// destroy buffer

TBufferJSON::~TBufferJSON()
{
   while (fStack.size() > 0)
      PopStack();

   if (fNumericLocale.Length() > 0)
      setlocale(LC_NUMERIC, fNumericLocale.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Converts object, inherited from TObject class, to JSON string
/// Lower digit of compact parameter define formatting rules
///  - 0 - no any compression, human-readable form
///  - 1 - exclude spaces in the begin
///  - 2 - remove newlines
///  - 3 - exclude spaces as much as possible
///
/// Second digit of compact parameter defines algorithm for arrays compression
///  - 0 - no compression, standard JSON array
///  - 1 - exclude leading and trailing zeros
///  - 2 - check values repetition and empty gaps
///
/// Maximal compression achieved when compact parameter equal to 23
/// When member_name specified, converts only this data member

TString TBufferJSON::ConvertToJSON(const TObject *obj, Int_t compact, const char *member_name)
{
   TClass *clActual = nullptr;
   void *ptr = (void *)obj;

   if (obj) {
      clActual = TObject::Class()->GetActualClass(obj);
      if (!clActual)
         clActual = TObject::Class();
      else if (clActual != TObject::Class())
         ptr = (void *)((Long_t)obj - clActual->GetBaseClassOffset(TObject::Class()));
   }

   return ConvertToJSON(ptr, clActual, compact, member_name);
}

////////////////////////////////////////////////////////////////////////////////
/// Set level of space/newline/array compression
/// Lower digit of compact parameter define formatting rules
///  - kNoCompress = 0  - no any compression, human-readable form
///  - kNoIndent = 1    - remove indentation spaces in the begin of each line
///  - kNoNewLine = 2   - remove also newlines
///  - kNoSpaces = 3    - exclude all spaces and new lines
///
/// Second digit of compact parameter defines algorithm for arrays compression
///  - 0 - no compression, standard JSON array
///  - kZeroSuppression = 10  - exclude leading and trailing zeros
///  - kSameSuppression = 20 - check values repetition and empty gaps
///
/// Third digit defines usage of typeinfo
///  - kSkipTypeInfo = 100   - "_typename" field will be skipped, reading by ROOT or JSROOT may be impossible

void TBufferJSON::SetCompact(int level)
{
   if (level < 0)
      level = 0;
   fCompact = level % 10;
   if (fCompact >= kMapAsObject) {
      fMapAsObject = kTRUE;
      fCompact = fCompact % kMapAsObject;
   }
   fSemicolon = (fCompact >= kNoSpaces) ? ":" : " : ";
   fArraySepar = (fCompact >= kNoSpaces) ? "," : ", ";
   fArrayCompact = ((level / 10) % 10) * 10;
   if ((((level / 100) % 10) * 100) == kSkipTypeInfo)
      fTypeNameTag.Clear();
   else if (fTypeNameTag.Length() == 0)
      fTypeNameTag = "_typename";
}

////////////////////////////////////////////////////////////////////////////////
/// Configures _typename tag in JSON structures
/// By default "_typename" field in JSON structures used to store class information
/// One can specify alternative tag like "$typename" or "xy", but such JSON can not be correctly used in JSROOT
/// If empty string is provided, class information will not be stored

void TBufferJSON::SetTypenameTag(const char *tag)
{
   if (!tag)
      fTypeNameTag.Clear();
   else
      fTypeNameTag = tag;
}

////////////////////////////////////////////////////////////////////////////////
/// Configures _typeversion tag in JSON
/// One can specify name of the JSON tag like "_typeversion" or "$tv" which will be used to store class version
/// Such tag can be used to correctly recover objects from JSON
/// If empty string is provided (default), class version will not be stored

void TBufferJSON::SetTypeversionTag(const char *tag)
{
   if (!tag)
      fTypeVersionTag.Clear();
   else
      fTypeVersionTag = tag;
}

////////////////////////////////////////////////////////////////////////////////
/// Specify class which typename will not be stored in JSON
/// Several classes can be configured
/// To exclude typeinfo for all classes, call TBufferJSON::SetTypenameTag("")

void TBufferJSON::SetSkipClassInfo(const TClass *cl)
{
   if (cl && (std::find(fSkipClasses.begin(), fSkipClasses.end(), cl) == fSkipClasses.end()))
      fSkipClasses.emplace_back(cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if class info will be skipped from JSON

Bool_t TBufferJSON::IsSkipClassInfo(const TClass *cl) const
{
   return cl && (std::find(fSkipClasses.begin(), fSkipClasses.end(), cl) != fSkipClasses.end());
}

////////////////////////////////////////////////////////////////////////////////
/// Converts any type of object to JSON string
/// One should provide pointer on object and its class name
/// Lower digit of compact parameter define formatting rules
///  - TBufferJSON::kNoCompress (0) - no any compression, human-readable form
///  - TBufferJSON::kNoIndent (1) - exclude spaces in the begin
///  - TBufferJSON::kNoNewLine (2) - no indent and no newlines
///  - TBufferJSON::kNoSpaces (3) - exclude spaces as much as possible
/// Second digit of compact parameter defines algorithm for arrays compression
///  - 0 - no compression, standard JSON array
///  - TBufferJSON::kZeroSuppression (10) - exclude leading and trailing zeros
///  - TBufferJSON::kSameSuppression (20) - check values repetition and empty gaps
///  - TBufferJSON::kBase64 (30) - arrays will be coded with base64 coding
/// Third digit of compact parameter defines typeinfo storage:
///  - TBufferJSON::kSkipTypeInfo (100) - "_typename" will be skipped, not always can be read back
/// Maximal none-destructive compression can be achieved when
/// compact parameter equal to TBufferJSON::kNoSpaces + TBufferJSON::kSameSuppression
/// When member_name specified, converts only this data member

TString TBufferJSON::ConvertToJSON(const void *obj, const TClass *cl, Int_t compact, const char *member_name)
{
   TClass *clActual = obj ? cl->GetActualClass(obj) : nullptr;
   const void *actualStart = obj;
   if (clActual && (clActual != cl)) {
      actualStart = (char *)obj - clActual->GetBaseClassOffset(cl);
   } else {
      // We could not determine the real type of this object,
      // let's assume it is the one given by the caller.
      clActual = const_cast<TClass *>(cl);
   }

   if (member_name && actualStart) {
      TRealData *rdata = clActual->GetRealData(member_name);
      TDataMember *member = rdata ? rdata->GetDataMember() : nullptr;
      if (!member) {
         TIter iter(clActual->GetListOfRealData());
         while ((rdata = dynamic_cast<TRealData *>(iter())) != nullptr) {
            member = rdata->GetDataMember();
            if (member && strcmp(member->GetName(), member_name) == 0)
               break;
         }
      }
      if (!member)
         return TString();

      Int_t arraylen = -1;
      if (member->GetArrayIndex() != 0) {
         TRealData *idata = clActual->GetRealData(member->GetArrayIndex());
         TDataMember *imember = idata ? idata->GetDataMember() : nullptr;
         if (imember && (strcmp(imember->GetTrueTypeName(), "int") == 0)) {
            arraylen = *((int *)((char *)actualStart + idata->GetThisOffset()));
         }
      }

      void *ptr = (char *)actualStart + rdata->GetThisOffset();
      if (member->IsaPointer())
         ptr = *((char **)ptr);

      return TBufferJSON::ConvertToJSON(ptr, member, compact, arraylen);
   }

   TBufferJSON buf;

   buf.SetCompact(compact);

   return buf.StoreObject(actualStart, clActual);
}

////////////////////////////////////////////////////////////////////////////////
/// Store provided object as JSON structure
/// Allows to configure different TBufferJSON properties before converting object into JSON
/// Actual object class must be specified here
/// Method can be safely called once - after that TBufferJSON instance must be destroyed
/// Code should look like:
///
///   auto obj = new UserClass();
///   TBufferJSON buf;
///   buf.SetCompact(TBufferJSON::kNoSpaces); // change any other settings in TBufferJSON
///   auto json = buf.StoreObject(obj, TClass::GetClass<UserClass>());
///

TString TBufferJSON::StoreObject(const void *obj, const TClass *cl)
{
   if (IsWriting()) {

      InitMap();

      PushStack(); // dummy stack entry to avoid extra checks in the beginning

      JsonWriteObject(obj, cl);

      PopStack();
   } else {
      Error("StoreObject", "Can not store object into TBuffer for reading");
   }

   return fOutBuffer.Length() ? fOutBuffer : fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Converts selected data member into json
/// Parameter ptr specifies address in memory, where data member is located
/// compact parameter defines compactness of produced JSON (from 0 to 3)
/// arraylen (when specified) is array length for this data member,  //[fN] case

TString TBufferJSON::ConvertToJSON(const void *ptr, TDataMember *member, Int_t compact, Int_t arraylen)
{
   if (!ptr || !member)
      return TString("null");

   Bool_t stlstring = !strcmp(member->GetTrueTypeName(), "string");

   Int_t isstl = member->IsSTLContainer();

   TClass *mcl = member->IsBasic() ? nullptr : gROOT->GetClass(member->GetTypeName());

   if (mcl && (mcl != TString::Class()) && !stlstring && !isstl && (mcl->GetBaseClassOffset(TArray::Class()) != 0) &&
       (arraylen <= 0) && (member->GetArrayDim() == 0))
      return TBufferJSON::ConvertToJSON(ptr, mcl, compact);

   TBufferJSON buf;

   buf.SetCompact(compact);

   return buf.JsonWriteMember(ptr, member, mcl, arraylen);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert object into JSON and store in text file
/// Returns size of the produce file
/// Used in TObject::SaveAs()

Int_t TBufferJSON::ExportToFile(const char *filename, const TObject *obj, const char *option)
{
   if (!obj || !filename || (*filename == 0))
      return 0;

   Int_t compact = strstr(filename, ".json.gz") ? 3 : 0;
   if (option && (*option >= '0') && (*option <= '3'))
      compact = TString(option).Atoi();

   TString json = TBufferJSON::ConvertToJSON(obj, compact);

   std::ofstream ofs(filename);

   if (strstr(filename, ".json.gz")) {
      const char *objbuf = json.Data();
      Long_t objlen = json.Length();

      unsigned long objcrc = R__crc32(0, NULL, 0);
      objcrc = R__crc32(objcrc, (const unsigned char *)objbuf, objlen);

      // 10 bytes (ZIP header), compressed data, 8 bytes (CRC and original length)
      Int_t buflen = 10 + objlen + 8;
      if (buflen < 512)
         buflen = 512;

      char *buffer = (char *)malloc(buflen);
      if (!buffer)
         return 0; // failure

      char *bufcur = buffer;

      *bufcur++ = 0x1f; // first byte of ZIP identifier
      *bufcur++ = 0x8b; // second byte of ZIP identifier
      *bufcur++ = 0x08; // compression method
      *bufcur++ = 0x00; // FLAG - empty, no any file names
      *bufcur++ = 0;    // empty timestamp
      *bufcur++ = 0;    //
      *bufcur++ = 0;    //
      *bufcur++ = 0;    //
      *bufcur++ = 0;    // XFL (eXtra FLags)
      *bufcur++ = 3;    // OS   3 means Unix
      // strcpy(bufcur, "item.json");
      // bufcur += strlen("item.json")+1;

      char dummy[8];
      memcpy(dummy, bufcur - 6, 6);

      // R__memcompress fills first 6 bytes with own header, therefore just overwrite them
      unsigned long ziplen = R__memcompress(bufcur - 6, objlen + 6, (char *)objbuf, objlen);

      memcpy(bufcur - 6, dummy, 6);

      bufcur += (ziplen - 6); // jump over compressed data (6 byte is extra ROOT header)

      *bufcur++ = objcrc & 0xff; // CRC32
      *bufcur++ = (objcrc >> 8) & 0xff;
      *bufcur++ = (objcrc >> 16) & 0xff;
      *bufcur++ = (objcrc >> 24) & 0xff;

      *bufcur++ = objlen & 0xff;         // original data length
      *bufcur++ = (objlen >> 8) & 0xff;  // original data length
      *bufcur++ = (objlen >> 16) & 0xff; // original data length
      *bufcur++ = (objlen >> 24) & 0xff; // original data length

      ofs.write(buffer, bufcur - buffer);

      free(buffer);
   } else {
      ofs << json.Data();
   }

   ofs.close();

   return json.Length();
}

////////////////////////////////////////////////////////////////////////////////
/// Convert object into JSON and store in text file
/// Returns size of the produce file

Int_t TBufferJSON::ExportToFile(const char *filename, const void *obj, const TClass *cl, const char *option)
{
   if (!obj || !cl || !filename || (*filename == 0))
      return 0;

   Int_t compact = strstr(filename, ".json.gz") ? 3 : 0;
   if (option && (*option >= '0') && (*option <= '3'))
      compact = TString(option).Atoi();

   TString json = TBufferJSON::ConvertToJSON(obj, cl, compact);

   std::ofstream ofs(filename);

   if (strstr(filename, ".json.gz")) {
      const char *objbuf = json.Data();
      Long_t objlen = json.Length();

      unsigned long objcrc = R__crc32(0, NULL, 0);
      objcrc = R__crc32(objcrc, (const unsigned char *)objbuf, objlen);

      // 10 bytes (ZIP header), compressed data, 8 bytes (CRC and original length)
      Int_t buflen = 10 + objlen + 8;
      if (buflen < 512)
         buflen = 512;

      char *buffer = (char *)malloc(buflen);
      if (!buffer)
         return 0; // failure

      char *bufcur = buffer;

      *bufcur++ = 0x1f; // first byte of ZIP identifier
      *bufcur++ = 0x8b; // second byte of ZIP identifier
      *bufcur++ = 0x08; // compression method
      *bufcur++ = 0x00; // FLAG - empty, no any file names
      *bufcur++ = 0;    // empty timestamp
      *bufcur++ = 0;    //
      *bufcur++ = 0;    //
      *bufcur++ = 0;    //
      *bufcur++ = 0;    // XFL (eXtra FLags)
      *bufcur++ = 3;    // OS   3 means Unix
      // strcpy(bufcur, "item.json");
      // bufcur += strlen("item.json")+1;

      char dummy[8];
      memcpy(dummy, bufcur - 6, 6);

      // R__memcompress fills first 6 bytes with own header, therefore just overwrite them
      unsigned long ziplen = R__memcompress(bufcur - 6, objlen + 6, (char *)objbuf, objlen);

      memcpy(bufcur - 6, dummy, 6);

      bufcur += (ziplen - 6); // jump over compressed data (6 byte is extra ROOT header)

      *bufcur++ = objcrc & 0xff; // CRC32
      *bufcur++ = (objcrc >> 8) & 0xff;
      *bufcur++ = (objcrc >> 16) & 0xff;
      *bufcur++ = (objcrc >> 24) & 0xff;

      *bufcur++ = objlen & 0xff;         // original data length
      *bufcur++ = (objlen >> 8) & 0xff;  // original data length
      *bufcur++ = (objlen >> 16) & 0xff; // original data length
      *bufcur++ = (objlen >> 24) & 0xff; // original data length

      ofs.write(buffer, bufcur - buffer);

      free(buffer);
   } else {
      ofs << json.Data();
   }

   ofs.close();

   return json.Length();
}

////////////////////////////////////////////////////////////////////////////////
/// Read TObject-based class from JSON, produced by ConvertToJSON() method.
/// If object does not inherit from TObject class, return 0.

TObject *TBufferJSON::ConvertFromJSON(const char *str)
{
   TClass *cl = nullptr;
   void *obj = ConvertFromJSONAny(str, &cl);

   if (!cl || !obj)
      return nullptr;

   Int_t delta = cl->GetBaseClassOffset(TObject::Class());

   if (delta < 0) {
      cl->Destructor(obj);
      return nullptr;
   }

   return (TObject *)(((char *)obj) + delta);
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from JSON
/// In class pointer (if specified) read class is returned
/// One must specify expected object class, if it is TArray or STL container

void *TBufferJSON::ConvertFromJSONAny(const char *str, TClass **cl)
{
   TBufferJSON buf(TBuffer::kRead);

   return buf.RestoreObject(str, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from JSON
/// In class pointer (if specified) read class is returned
/// One must specify expected object class, if it is TArray or STL container

void *TBufferJSON::RestoreObject(const char *json_str, TClass **cl)
{
   if (!IsReading())
      return nullptr;

   nlohmann::json docu = nlohmann::json::parse(json_str);

   if (docu.is_null() || (!docu.is_object() && !docu.is_array()))
      return nullptr;

   TClass *objClass = nullptr;

   if (cl) {
      objClass = *cl; // this is class which suppose to created when reading JSON
      *cl = nullptr;
   }

   InitMap();

   PushStack(0, &docu);

   void *obj = JsonReadObject(nullptr, objClass, cl);

   PopStack();

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Read objects from JSON, one can reuse existing object

void *TBufferJSON::ConvertFromJSONChecked(const char *str, const TClass *expectedClass)
{
   if (!expectedClass)
      return nullptr;

   TClass *resClass = const_cast<TClass *>(expectedClass);

   void *res = ConvertFromJSONAny(str, &resClass);

   if (!res || !resClass)
      return nullptr;

   if (resClass == expectedClass)
      return res;

   Int_t offset = resClass->GetBaseClassOffset(expectedClass);
   if (offset < 0) {
      ::Error("TBufferJSON::ConvertFromJSONChecked", "expected class %s is not base for read class %s",
              expectedClass->GetName(), resClass->GetName());
      resClass->Destructor(res);
      return nullptr;
   }

   return (char *)res - offset;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert single data member to JSON structures
/// Returns string with converted member

TString TBufferJSON::JsonWriteMember(const void *ptr, TDataMember *member, TClass *memberClass, Int_t arraylen)
{
   if (!member)
      return "null";

   if (gDebug > 2)
      Info("JsonWriteMember", "Write member %s type %s ndim %d", member->GetName(), member->GetTrueTypeName(),
           member->GetArrayDim());

   Int_t tid = member->GetDataType() ? member->GetDataType()->GetType() : kNoType_t;
   if (strcmp(member->GetTrueTypeName(), "const char*") == 0)
      tid = kCharStar;
   else if (!member->IsBasic() || (tid == kOther_t) || (tid == kVoid_t))
      tid = kNoType_t;

   if (!ptr)
      return (tid == kCharStar) ? "\"\"" : "null";

   PushStack(0);
   fValue.Clear();

   if (tid != kNoType_t) {

      TArrayIndexProducer indx(member, arraylen, fArraySepar.Data());

      Int_t shift = 1;

      if (indx.IsArray() && (tid == kChar_t))
         shift = indx.ReduceDimension();

      char *ppp = (char *)ptr;

      if (indx.IsArray())
         fOutBuffer.Append(indx.GetBegin());

      do {
         fValue.Clear();

         switch (tid) {
         case kChar_t:
            if (shift > 1)
               JsonWriteConstChar((Char_t *)ppp, shift);
            else
               JsonWriteBasic(*((Char_t *)ppp));
            break;
         case kShort_t: JsonWriteBasic(*((Short_t *)ppp)); break;
         case kInt_t: JsonWriteBasic(*((Int_t *)ppp)); break;
         case kLong_t: JsonWriteBasic(*((Long_t *)ppp)); break;
         case kFloat_t: JsonWriteBasic(*((Float_t *)ppp)); break;
         case kCounter: JsonWriteBasic(*((Int_t *)ppp)); break;
         case kCharStar: JsonWriteConstChar((Char_t *)ppp); break;
         case kDouble_t: JsonWriteBasic(*((Double_t *)ppp)); break;
         case kDouble32_t: JsonWriteBasic(*((Double_t *)ppp)); break;
         case kchar: JsonWriteBasic(*((char *)ppp)); break;
         case kUChar_t: JsonWriteBasic(*((UChar_t *)ppp)); break;
         case kUShort_t: JsonWriteBasic(*((UShort_t *)ppp)); break;
         case kUInt_t: JsonWriteBasic(*((UInt_t *)ppp)); break;
         case kULong_t: JsonWriteBasic(*((ULong_t *)ppp)); break;
         case kBits: JsonWriteBasic(*((UInt_t *)ppp)); break;
         case kLong64_t: JsonWriteBasic(*((Long64_t *)ppp)); break;
         case kULong64_t: JsonWriteBasic(*((ULong64_t *)ppp)); break;
         case kBool_t: JsonWriteBasic(*((Bool_t *)ppp)); break;
         case kFloat16_t: JsonWriteBasic(*((Float_t *)ppp)); break;
         case kOther_t:
         case kVoid_t: break;
         }

         fOutBuffer.Append(fValue);
         if (indx.IsArray())
            fOutBuffer.Append(indx.NextSeparator());

         ppp += shift * member->GetUnitSize();

      } while (!indx.IsDone());

      fValue = fOutBuffer;

   } else if (memberClass == TString::Class()) {
      TString *str = (TString *)ptr;
      JsonWriteConstChar(str ? str->Data() : nullptr);
   } else if ((member->IsSTLContainer() == ROOT::kSTLvector) || (member->IsSTLContainer() == ROOT::kSTLlist) ||
              (member->IsSTLContainer() == ROOT::kSTLforwardlist)) {

      if (memberClass)
         memberClass->Streamer((void *)ptr, *this);
      else
         fValue = "[]";

      if (fValue == "0")
         fValue = "[]";

   } else if (memberClass && memberClass->GetBaseClassOffset(TArray::Class()) == 0) {
      TArray *arr = (TArray *)ptr;
      if (arr && (arr->GetSize() > 0)) {
         arr->Streamer(*this);
         // WriteFastArray(arr->GetArray(), arr->GetSize());
         if (Stack()->fValues.size() > 1) {
            Warning("TBufferJSON", "When streaming TArray, more than 1 object in the stack, use second item");
            fValue = Stack()->fValues[1].c_str();
         }
      } else
         fValue = "[]";
   } else if (memberClass && !strcmp(memberClass->GetName(), "string")) {
      // here value contains quotes, stack can be ignored
      memberClass->Streamer((void *)ptr, *this);
   }
   PopStack();

   if (fValue.Length())
      return fValue;

   if (!memberClass || (member->GetArrayDim() > 0) || (arraylen > 0))
      return "<not supported>";

   return TBufferJSON::ConvertToJSON(ptr, memberClass);
}

////////////////////////////////////////////////////////////////////////////////
/// add new level to the structures stack

TJSONStackObj *TBufferJSON::PushStack(Int_t inclevel, void *readnode)
{
   auto next = new TJSONStackObj();
   next->fLevel = inclevel;
   if (IsReading()) {
      next->fNode = (nlohmann::json *)readnode;
   } else if (fStack.size() > 0) {
      auto prev = Stack();
      next->fLevel += prev->fLevel;
      next->fMemberPtr = prev->fMemberPtr;
   }
   fStack.emplace_back(next);
   return next;
}

////////////////////////////////////////////////////////////////////////////////
/// remove one level from stack

TJSONStackObj *TBufferJSON::PopStack()
{
   if (fStack.size() > 0)
      fStack.pop_back();

   return fStack.size() > 0 ? fStack.back().get() : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Append two string to the output JSON, normally separate by line break

void TBufferJSON::AppendOutput(const char *line0, const char *line1)
{
   if (line0)
      fOutput->Append(line0);

   if (line1) {
      if (fCompact < 2)
         fOutput->Append("\n");

      if (strlen(line1) > 0) {
         if (fCompact < 1) {
            if (Stack()->fLevel > 0)
               fOutput->Append(' ', Stack()->fLevel);
         }
         fOutput->Append(line1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Start object element with typeinfo

TJSONStackObj *TBufferJSON::JsonStartObjectWrite(const TClass *obj_class, TStreamerInfo *info)
{
   auto stack = PushStack(2);

   // new object started - assign own member counter
   stack->fMemberPtr = &stack->fMemberCnt;

   if ((fTypeNameTag.Length() > 0) && !IsSkipClassInfo(obj_class)) {
      // stack->fMemberCnt = 1; // default value, comment out here
      AppendOutput("{", "\"");
      AppendOutput(fTypeNameTag.Data());
      AppendOutput("\"");
      AppendOutput(fSemicolon.Data());
      AppendOutput("\"");
      AppendOutput(obj_class->GetName());
      AppendOutput("\"");
      if (fTypeVersionTag.Length() > 0) {
         AppendOutput(stack->NextMemberSeparator(), "\"");
         AppendOutput(fTypeVersionTag.Data());
         AppendOutput("\"");
         AppendOutput(fSemicolon.Data());
         AppendOutput(Form("%d", (int)(info ? info->GetClassVersion() : obj_class->GetClassVersion())));
      }
   } else {
      stack->fMemberCnt = 0; // exclude typename
      AppendOutput("{");
   }

   return stack;
}

////////////////////////////////////////////////////////////////////////////////
/// Start new class member in JSON structures

void TBufferJSON::JsonStartElement(const TStreamerElement *elem, const TClass *base_class)
{
   const char *elem_name = nullptr;
   Int_t special_kind = JsonSpecialClass(base_class);

   switch (special_kind) {
   case 0:
      if (!base_class)
         elem_name = elem->GetName();
      break;
   case TClassEdit::kVector: elem_name = "fVector"; break;
   case TClassEdit::kList: elem_name = "fList"; break;
   case TClassEdit::kForwardlist: elem_name = "fForwardlist"; break;
   case TClassEdit::kDeque: elem_name = "fDeque"; break;
   case TClassEdit::kMap: elem_name = "fMap"; break;
   case TClassEdit::kMultiMap: elem_name = "fMultiMap"; break;
   case TClassEdit::kSet: elem_name = "fSet"; break;
   case TClassEdit::kMultiSet: elem_name = "fMultiSet"; break;
   case TClassEdit::kUnorderedSet: elem_name = "fUnorderedSet"; break;
   case TClassEdit::kUnorderedMultiSet: elem_name = "fUnorderedMultiSet"; break;
   case TClassEdit::kUnorderedMap: elem_name = "fUnorderedMap"; break;
   case TClassEdit::kUnorderedMultiMap: elem_name = "fUnorderedMultiMap"; break;
   case TClassEdit::kBitSet: elem_name = "fBitSet"; break;
   case json_TArray: elem_name = "fArray"; break;
   case json_TString:
   case json_stdstring: elem_name = "fString"; break;
   }

   if (!elem_name)
      return;

   if (IsReading()) {
      nlohmann::json *json = Stack()->fNode;

      if (json->count(elem_name) != 1) {
         Error("JsonStartElement", "Missing JSON structure for element %s", elem_name);
      } else {
         Stack()->fNode = &((*json)[elem_name]);
         if (special_kind == json_TArray) {
            Int_t len = Stack()->IsJsonArray();
            Stack()->PushIntValue(len > 0 ? len : 0);
            if (len < 0)
               Error("JsonStartElement", "Missing array when reading TArray class for element %s", elem->GetName());
         }
         if ((gDebug > 1) && base_class)
            Info("JsonStartElement", "Reading baseclass %s from element %s", base_class->GetName(), elem_name);
      }

   } else {
      AppendOutput(Stack()->NextMemberSeparator(), "\"");
      AppendOutput(elem_name);
      AppendOutput("\"");
      AppendOutput(fSemicolon.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// disable post-processing of the code
void TBufferJSON::JsonDisablePostprocessing()
{
   Stack()->fIsPostProcessed = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// return non-zero value when class has special handling in JSON
/// it is TCollection (-130), TArray (100), TString (110), std::string (120) and STL containers (1..6)

Int_t TBufferJSON::JsonSpecialClass(const TClass *cl) const
{
   if (!cl)
      return 0;

   Bool_t isarray = strncmp("TArray", cl->GetName(), 6) == 0;
   if (isarray)
      isarray = (const_cast<TClass *>(cl))->GetBaseClassOffset(TArray::Class()) == 0;
   if (isarray)
      return json_TArray;

   // negative value used to indicate that collection stored as object
   if ((const_cast<TClass *>(cl))->GetBaseClassOffset(TCollection::Class()) == 0)
      return json_TCollection;

   // special case for TString - it is saved as string in JSON
   if (cl == TString::Class())
      return json_TString;

   bool isstd = TClassEdit::IsStdClass(cl->GetName());
   int isstlcont(ROOT::kNotSTL);
   if (isstd)
      isstlcont = cl->GetCollectionType();
   if (isstlcont > 0)
      return isstlcont;

   // also special handling for STL string, which handled similar to TString
   if (isstd && !strcmp(cl->GetName(), "string"))
      return json_stdstring;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to buffer
/// If object was written before, only pointer will be stored
/// If check_map==kFALSE, object will be stored in any case and pointer will not be registered in the map

void TBufferJSON::JsonWriteObject(const void *obj, const TClass *cl, Bool_t check_map)
{
   if (!cl)
      obj = nullptr;

   if (gDebug > 0)
      Info("JsonWriteObject", "Object %p class %s check_map %s", obj, cl ? cl->GetName() : "null",
           check_map ? "true" : "false");

   Int_t special_kind = JsonSpecialClass(cl), map_convert{0};

   TString fObjectOutput, *fPrevOutput{nullptr};

   TJSONStackObj *stack = Stack();

   if (stack && stack->fAccObjects && ((fValue.Length() > 0) || (stack->fValues.size() > 0))) {
      // accumulate data of super-object in stack

      if (fValue.Length() > 0)
         stack->PushValue(fValue);

      // redirect output to local buffer, use it later as value
      fPrevOutput = fOutput;
      fOutput = &fObjectOutput;
   } else if ((special_kind <= 0) || (special_kind > json_TArray)) {
      // FIXME: later post processing should be active for all special classes, while they all keep output in the value
      JsonDisablePostprocessing();
   } else if ((special_kind == TClassEdit::kMap) || (special_kind == TClassEdit::kMultiMap) ||
              (special_kind == TClassEdit::kUnorderedMap) || (special_kind == TClassEdit::kUnorderedMultiMap)) {

      if ((fMapAsObject && (fStack.size()==1)) || (stack && stack->fElem && strstr(stack->fElem->GetTitle(), "JSON_object")))
         map_convert = 2; // mapped into normal object
      else
         map_convert = 1;
   }

   if (!obj) {
      AppendOutput("null");
      goto post_process;
   }

   if (special_kind <= 0) {
      // add element name which should correspond to the object
      if (check_map) {
         Long64_t refid = GetObjectTag(obj);
         if (refid > 0) {
            // old-style refs, coded into string like "$ref12"
            // AppendOutput(Form("\"$ref:%u\"", iter->second));
            // new-style refs, coded into extra object {"$ref":12}, auto-detected by JSROOT 4.8 and higher
            AppendOutput(Form("{\"$ref\":%u}", (unsigned)(refid - 1)));
            goto post_process;
         }
         MapObject(obj, cl, fJsonrCnt + 1); // +1 used
      }

      fJsonrCnt++; // object counts required in dereferencing part

      stack = JsonStartObjectWrite(cl);

   } else if (map_convert == 2) {
      // special handling of map - it is object, but stored in the fValue

      if (check_map) {
         Long64_t refid = GetObjectTag(obj);
         if (refid > 0) {
            fValue.Form("{\"$ref\":%u}", (unsigned)(refid - 1));
            goto post_process;
         }
         MapObject(obj, cl, fJsonrCnt + 1); // +1 used
      }

      fJsonrCnt++; // object counts required in dereferencing part
      stack = PushStack(0);

   } else {

      bool base64 = ((special_kind == TClassEdit::kVector) && stack && stack->fElem && strstr(stack->fElem->GetTitle(), "JSON_base64"));

      // for array, string and STL collections different handling -
      // they not recognized at the end as objects in JSON
      stack = PushStack(0);

      stack->fBase64 = base64;
   }

   if (gDebug > 3)
      Info("JsonWriteObject", "Starting object %p write for class: %s", obj, cl->GetName());

   stack->fAccObjects = special_kind < ROOT::kSTLend;

   if (special_kind == json_TCollection)
      JsonWriteCollection((TCollection *)obj, cl);
   else
      (const_cast<TClass *>(cl))->Streamer((void *)obj, *this);

   if (gDebug > 3)
      Info("JsonWriteObject", "Done object %p write for class: %s", obj, cl->GetName());

   if (special_kind == json_TArray) {
      if (stack->fValues.size() != 1)
         Error("JsonWriteObject", "Problem when writing array");
      stack->fValues.clear();
   } else if ((special_kind == json_TString) || (special_kind == json_stdstring)) {
      if (stack->fValues.size() > 2)
         Error("JsonWriteObject", "Problem when writing TString or std::string");
      stack->fValues.clear();
      AppendOutput(fValue.Data());
      fValue.Clear();
   } else if ((special_kind > 0) && (special_kind < ROOT::kSTLend)) {
      // here make STL container processing

      if (map_convert == 2) {
         // converting map into object

         if (!stack->fValues.empty() && (fValue.Length() > 0))
            stack->PushValue(fValue);

         const char *separ = (fCompact < 2) ? ", " : ",";
         const char *semi = (fCompact < 2) ? ": " : ":";
         bool first = true;

         fValue = "{";
         if ((fTypeNameTag.Length() > 0) && !IsSkipClassInfo(cl)) {
            fValue.Append("\"");
            fValue.Append(fTypeNameTag);
            fValue.Append("\"");
            fValue.Append(semi);
            fValue.Append("\"");
            fValue.Append(cl->GetName());
            fValue.Append("\"");
            first = false;
         }
         for (Int_t k = 1; k < (int)stack->fValues.size() - 1; k += 2) {
            if (!first)
               fValue.Append(separ);
            first = false;
            fValue.Append(stack->fValues[k].c_str());
            fValue.Append(semi);
            fValue.Append(stack->fValues[k + 1].c_str());
         }
         fValue.Append("}");
         stack->fValues.clear();
      } else if (stack->fValues.empty()) {
         // empty container
         if (fValue != "0")
            Error("JsonWriteObject", "With empty stack fValue!=0");
         fValue = "[]";
      } else {

         auto size = std::stoi(stack->fValues[0]);

         bool trivial_format = false;

         if ((stack->fValues.size() == 1) && ((size > 1) || ((fValue.Length() > 1) && (fValue[0]=='[')))) {
            // prevent case of vector<vector<value_class>>
            const auto proxy = cl->GetCollectionProxy();
            TClass *value_class = proxy ? proxy->GetValueClass() : nullptr;
            if (value_class && TClassEdit::IsStdClass(value_class->GetName()) && (value_class->GetCollectionType() != ROOT::kNotSTL))
               trivial_format = false;
            else
               trivial_format = true;
         }

         if (trivial_format) {
            // case of simple vector, array already in the value
            stack->fValues.clear();
            if (fValue.Length() == 0) {
               Error("JsonWriteObject", "Empty value when it should contain something");
               fValue = "[]";
            }

         } else {
            const char *separ = "[";

            if (fValue.Length() > 0)
               stack->PushValue(fValue);

            if ((size * 2 == (int) stack->fValues.size() - 1) && (map_convert > 0)) {
               // special handling for std::map.
               // Create entries like { '$pair': 'typename' , 'first' : key, 'second' : value }
               TString pairtype = cl->GetName();
               if (pairtype.Index("unordered_map<") == 0)
                  pairtype.Replace(0, 14, "pair<");
               else if (pairtype.Index("unordered_multimap<") == 0)
                  pairtype.Replace(0, 19, "pair<");
               else if (pairtype.Index("multimap<") == 0)
                  pairtype.Replace(0, 9, "pair<");
               else if (pairtype.Index("map<") == 0)
                  pairtype.Replace(0, 4, "pair<");
               else
                  pairtype = "TPair";
               if (fTypeNameTag.Length() == 0)
                  pairtype = "1";
               else
                  pairtype = TString("\"") + pairtype + TString("\"");
               for (Int_t k = 1; k < (int) stack->fValues.size() - 1; k += 2) {
                  fValue.Append(separ);
                  separ = fArraySepar.Data();
                  // fJsonrCnt++; // do not add entry in the map, can conflict with objects inside values
                  fValue.Append("{");
                  fValue.Append("\"$pair\"");
                  fValue.Append(fSemicolon);
                  fValue.Append(pairtype.Data());
                  fValue.Append(fArraySepar);
                  fValue.Append("\"first\"");
                  fValue.Append(fSemicolon);
                  fValue.Append(stack->fValues[k].c_str());
                  fValue.Append(fArraySepar);
                  fValue.Append("\"second\"");
                  fValue.Append(fSemicolon);
                  fValue.Append(stack->fValues[k + 1].c_str());
                  fValue.Append("}");
               }
            } else {
               // for most stl containers write just like blob, but skipping first element with size
               for (Int_t k = 1; k < (int) stack->fValues.size(); k++) {
                  fValue.Append(separ);
                  separ = fArraySepar.Data();
                  fValue.Append(stack->fValues[k].c_str());
               }
            }

            fValue.Append("]");
            stack->fValues.clear();
         }
      }
   }

   // reuse post-processing code for TObject or TRef
   PerformPostProcessing(stack, cl);

   if ((special_kind == 0) && (!stack->fValues.empty() || (fValue.Length() > 0))) {
      if (gDebug > 0)
         Info("JsonWriteObject", "Create blob value for class %s", cl->GetName());

      AppendOutput(fArraySepar.Data(), "\"_blob\"");
      AppendOutput(fSemicolon.Data());

      const char *separ = "[";

      for (auto &elem: stack->fValues) {
         AppendOutput(separ);
         separ = fArraySepar.Data();
         AppendOutput(elem.c_str());
      }

      if (fValue.Length() > 0) {
         AppendOutput(separ);
         AppendOutput(fValue.Data());
      }

      AppendOutput("]");

      fValue.Clear();
      stack->fValues.clear();
   }

   PopStack();

   if ((special_kind <= 0))
      AppendOutput(nullptr, "}");

post_process:

   if (fPrevOutput) {
      fOutput = fPrevOutput;
      // for STL containers and TArray object in fValue itself
      if ((special_kind <= 0) || (special_kind > json_TArray))
         fValue = fObjectOutput;
      else if (fObjectOutput.Length() != 0)
         Error("JsonWriteObject", "Non-empty object output for special class %s", cl->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// store content of ROOT collection

void TBufferJSON::JsonWriteCollection(TCollection *col, const TClass *)
{
   AppendOutput(Stack()->NextMemberSeparator(), "\"name\"");
   AppendOutput(fSemicolon.Data());
   AppendOutput("\"");
   AppendOutput(col->GetName());
   AppendOutput("\"");
   AppendOutput(Stack()->NextMemberSeparator(), "\"arr\"");
   AppendOutput(fSemicolon.Data());

   // collection treated as JS Array
   AppendOutput("[");

   bool islist = col->InheritsFrom(TList::Class());
   TMap *map = nullptr;
   if (col->InheritsFrom(TMap::Class()))
      map = dynamic_cast<TMap *>(col);

   TString sopt;
   if (islist) {
      sopt.Capacity(500);
      sopt = "[";
   }

   TIter iter(col);
   TObject *obj;
   Bool_t first = kTRUE;
   while ((obj = iter()) != nullptr) {
      if (!first)
         AppendOutput(fArraySepar.Data());

      if (map) {
         // fJsonrCnt++; // do not account map pair as JSON object
         AppendOutput("{", "\"$pair\"");
         AppendOutput(fSemicolon.Data());
         AppendOutput("\"TPair\"");
         AppendOutput(fArraySepar.Data(), "\"first\"");
         AppendOutput(fSemicolon.Data());
      }

      WriteObjectAny(obj, TObject::Class());

      if (map) {
         AppendOutput(fArraySepar.Data(), "\"second\"");
         AppendOutput(fSemicolon.Data());
         WriteObjectAny(map->GetValue(obj), TObject::Class());
         AppendOutput("", "}");
      }

      if (islist) {
         if (!first)
            sopt.Append(fArraySepar.Data());
         sopt.Append("\"");
         sopt.Append(iter.GetOption());
         sopt.Append("\"");
      }

      first = kFALSE;
   }

   AppendOutput("]");

   if (islist) {
      sopt.Append("]");
      AppendOutput(Stack()->NextMemberSeparator(), "\"opt\"");
      AppendOutput(fSemicolon.Data());
      AppendOutput(sopt.Data());
   }
   fValue.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// read content of ROOT collection

void TBufferJSON::JsonReadCollection(TCollection *col, const TClass *)
{
   if (!col)
      return;

   TList *lst = nullptr;
   TMap *map = nullptr;
   TClonesArray *clones = nullptr;
   if (col->InheritsFrom(TList::Class()))
      lst = dynamic_cast<TList *>(col);
   else if (col->InheritsFrom(TMap::Class()))
      map = dynamic_cast<TMap *>(col);
   else if (col->InheritsFrom(TClonesArray::Class()))
      clones = dynamic_cast<TClonesArray *>(col);

   nlohmann::json *json = Stack()->fNode;

   std::string name = json->at("name");
   col->SetName(name.c_str());

   nlohmann::json &arr = json->at("arr");
   int size = arr.size();

   for (int n = 0; n < size; ++n) {
      nlohmann::json *subelem = &arr.at(n);

      if (map)
         subelem = &subelem->at("first");

      PushStack(0, subelem);

      TClass *readClass = nullptr, *objClass = nullptr;
      void *subobj = nullptr;

      if (clones) {
         if (n == 0) {
            if (!clones->GetClass() || (clones->GetSize() == 0)) {
               if (fTypeNameTag.Length() > 0) {
                  clones->SetClass(subelem->at(fTypeNameTag.Data()).get<std::string>().c_str(), size);
               } else {
                  Error("JsonReadCollection",
                        "Cannot detect class name for TClonesArray - typename tag not configured");
                  return;
               }
            } else if (size > clones->GetSize()) {
               Error("JsonReadCollection", "TClonesArray size %d smaller than required %d", clones->GetSize(), size);
               return;
            }
         }
         objClass = clones->GetClass();
         subobj = clones->ConstructedAt(n);
      }

      subobj = JsonReadObject(subobj, objClass, &readClass);

      PopStack();

      if (clones)
         continue;

      if (!subobj || !readClass) {
         subobj = nullptr;
      } else if (readClass->GetBaseClassOffset(TObject::Class()) != 0) {
         Error("JsonReadCollection", "Try to add object %s not derived from TObject", readClass->GetName());
         subobj = nullptr;
      }

      TObject *tobj = static_cast<TObject *>(subobj);

      if (map) {
         PushStack(0, &arr.at(n).at("second"));

         readClass = nullptr;
         void *subobj2 = JsonReadObject(nullptr, nullptr, &readClass);

         PopStack();

         if (!subobj2 || !readClass) {
            subobj2 = nullptr;
         } else if (readClass->GetBaseClassOffset(TObject::Class()) != 0) {
            Error("JsonReadCollection", "Try to add object %s not derived from TObject", readClass->GetName());
            subobj2 = nullptr;
         }

         map->Add(tobj, static_cast<TObject *>(subobj2));
      } else if (lst) {
         std::string opt = json->at("opt").at(n).get<std::string>();
         lst->Add(tobj, opt.c_str());
      } else {
         // generic method, all kinds of TCollection should work
         col->Add(tobj);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from current JSON node

void *TBufferJSON::JsonReadObject(void *obj, const TClass *objClass, TClass **readClass)
{
   if (readClass)
      *readClass = nullptr;

   TJSONStackObj *stack = Stack();

   Bool_t process_stl = stack->IsStl();
   nlohmann::json *json = stack->GetStlNode();

   // check if null pointer
   if (json->is_null())
      return nullptr;

   Int_t special_kind = JsonSpecialClass(objClass);

   // Extract pointer
   if (json->is_object() && (json->size() == 1) && (json->find("$ref") != json->end())) {
      unsigned refid = json->at("$ref").get<unsigned>();

      void *ref_obj = nullptr;
      TClass *ref_cl = nullptr;

      GetMappedObject(refid + 1, ref_obj, ref_cl);

      if (!ref_obj || !ref_cl) {
         Error("JsonReadObject", "Fail to find object for reference %u", refid);
         return nullptr;
      }

      if (readClass)
         *readClass = ref_cl;

      if (gDebug > 2)
         Info("JsonReadObject", "Extract object reference %u %p cl:%s expects:%s", refid, ref_obj, ref_cl->GetName(),
              (objClass ? objClass->GetName() : "---"));

      return ref_obj;
   }

   // special case of strings - they do not create JSON object, but just string
   if ((special_kind == json_stdstring) || (special_kind == json_TString)) {
      if (!obj)
         obj = objClass->New();

      if (gDebug > 2)
         Info("JsonReadObject", "Read string from %s", json->dump().c_str());

      if (special_kind == json_stdstring)
         *((std::string *)obj) = json->get<std::string>();
      else
         *((TString *)obj) = json->get<std::string>().c_str();

      if (readClass)
         *readClass = const_cast<TClass *>(objClass);

      return obj;
   }

   Bool_t isBase = (stack->fElem && objClass) ? stack->fElem->IsBase() : kFALSE; // base class

   if (isBase && (!obj || !objClass)) {
      Error("JsonReadObject", "No object when reading base class");
      return obj;
   }

   Int_t map_convert = 0;
   if ((special_kind == TClassEdit::kMap) || (special_kind == TClassEdit::kMultiMap) ||
       (special_kind == TClassEdit::kUnorderedMap) || (special_kind == TClassEdit::kUnorderedMultiMap)) {
      map_convert = json->is_object() ? 2 : 1; // check if map was written as array or as object
   }

   // from now all operations performed with sub-element,
   // stack should be repaired at the end
   if (process_stl)
      stack = PushStack(0, json);

   TClass *jsonClass = nullptr;
   Int_t jsonClassVersion = 0;

   if ((special_kind == json_TArray) || ((special_kind > 0) && (special_kind < ROOT::kSTLend))) {

      jsonClass = const_cast<TClass *>(objClass);

      if (!obj)
         obj = jsonClass->New();

      Int_t len = stack->IsJsonArray(json, map_convert == 2 ? fTypeNameTag.Data() : nullptr);

      stack->PushIntValue(len > 0 ? len : 0);

      if (len < 0) // should never happens
         Error("JsonReadObject", "Not array when expecting such %s", json->dump().c_str());

      if (gDebug > 1)
         Info("JsonReadObject", "Reading special kind %d %s ptr %p", special_kind, objClass->GetName(), obj);

   } else if (isBase) {
      // base class has special handling - no additional level and no extra refid

      jsonClass = const_cast<TClass *>(objClass);

      if (gDebug > 1)
         Info("JsonReadObject", "Reading baseclass %s ptr %p", objClass->GetName(), obj);
   } else {

      if ((fTypeNameTag.Length() > 0) && (json->count(fTypeNameTag.Data()) > 0)) {
         std::string clname = json->at(fTypeNameTag.Data()).get<std::string>();
         jsonClass = TClass::GetClass(clname.c_str());
         if (!jsonClass)
            Error("JsonReadObject", "Cannot find class %s", clname.c_str());
      } else {
         // try to use class which is assigned by streamers - better than nothing
         jsonClass = const_cast<TClass *>(objClass);
      }

      if (!jsonClass) {
         if (process_stl)
            PopStack();
         return obj;
      }

      if ((fTypeVersionTag.Length() > 0) && (json->count(fTypeVersionTag.Data()) > 0))
         jsonClassVersion = json->at(fTypeVersionTag.Data()).get<int>();

      if (objClass && (jsonClass != objClass)) {
         if (obj || (jsonClass->GetBaseClassOffset(objClass) != 0)) {
            if (jsonClass->GetBaseClassOffset(objClass) < 0)
               Error("JsonReadObject", "Not possible to read %s and casting to %s pointer as the two classes are unrelated",
                     jsonClass->GetName(), objClass->GetName());
            else
               Error("JsonReadObject", "Reading %s and casting to %s pointer is currently not supported",
                     jsonClass->GetName(), objClass->GetName());
            if (process_stl)
               PopStack();
            return obj;
         }
      }

      if (!obj)
         obj = jsonClass->New();

      if (gDebug > 1)
         Info("JsonReadObject", "Reading object of class %s refid %u ptr %p", jsonClass->GetName(), fJsonrCnt, obj);

      if (!special_kind)
         special_kind = JsonSpecialClass(jsonClass);

      // add new element to the reading map
      MapObject(obj, jsonClass, ++fJsonrCnt);
   }

   // there are two ways to handle custom streamers
   // either prepare data before streamer and tweak basic function which are reading values like UInt32_t
   // or try re-implement custom streamer here

   if ((jsonClass == TObject::Class()) || (jsonClass == TRef::Class())) {
      // for TObject we re-implement custom streamer - it is much easier

      JsonReadTObjectMembers((TObject *)obj, json);

   } else if (special_kind == json_TCollection) {

      JsonReadCollection((TCollection *)obj, jsonClass);

   } else {

      Bool_t do_read = kTRUE;

      // special handling of STL which coded into arrays
      if ((special_kind > 0) && (special_kind < ROOT::kSTLend))
         do_read = stack->AssignStl(jsonClass, map_convert, fTypeNameTag.Data());

      // if provided - use class version from JSON
      stack->fClVersion = jsonClassVersion ? jsonClassVersion : jsonClass->GetClassVersion();

      if (gDebug > 3)
         Info("JsonReadObject", "Calling streamer of class %s", jsonClass->GetName());

      if (isBase && (special_kind == 0))
         Error("JsonReadObject", "Should not be used for reading of base class %s", jsonClass->GetName());

      if (do_read)
         jsonClass->Streamer((void *)obj, *this);

      stack->fClVersion = 0;

      stack->ClearStl(); // reset STL index for itself to prevent looping
   }

   // return back stack position
   if (process_stl)
      PopStack();

   if (gDebug > 1)
      Info("JsonReadObject", "Reading object of class %s done", jsonClass->GetName());

   if (readClass)
      *readClass = jsonClass;

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Read TObject data members from JSON.
/// Do not call TObject::Streamer() to avoid special tweaking of TBufferJSON interface

void TBufferJSON::JsonReadTObjectMembers(TObject *tobj, void *node)
{
   nlohmann::json *json = node ? (nlohmann::json *)node : Stack()->fNode;

   UInt_t uid = json->at("fUniqueID").get<unsigned>();
   UInt_t bits = json->at("fBits").get<unsigned>();
   // UInt32_t pid = json->at("fPID").get<unsigned>(); // ignore PID for the moment

   tobj->SetUniqueID(uid);

   static auto tobj_fbits_offset = TObject::Class()->GetDataMemberOffset("fBits");

   // there is no method to set all bits directly - do it differently
   if (tobj_fbits_offset > 0) {
      UInt_t *fbits = (UInt_t *) ((char* ) tobj + tobj_fbits_offset);
      *fbits = (*fbits & (TObject::kIsOnHeap | TObject::kNotDeleted)) | bits;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and ReadBuffer functions
/// and indent new level in json structure.
/// This call indicates, that TStreamerInfo functions starts streaming
/// object data of correspondent class

void TBufferJSON::IncrementLevel(TVirtualStreamerInfo *info)
{
   if (gDebug > 2)
      Info("IncrementLevel", "Class: %s", (info ? info->GetClass()->GetName() : "custom"));

   WorkWithClass((TStreamerInfo *)info);
}

////////////////////////////////////////////////////////////////////////////////
/// Prepares buffer to stream data of specified class

void TBufferJSON::WorkWithClass(TStreamerInfo *sinfo, const TClass *cl)
{
   if (sinfo)
      cl = sinfo->GetClass();

   if (!cl)
      return;

   if (gDebug > 3)
      Info("WorkWithClass", "Class: %s", cl->GetName());

   TJSONStackObj *stack = Stack();

   if (IsReading()) {
      stack = PushStack(0, stack->fNode);
   } else if (stack && stack->IsStreamerElement() && !stack->fIsObjStarted &&
              ((stack->fElem->GetType() == TStreamerInfo::kObject) ||
               (stack->fElem->GetType() == TStreamerInfo::kAny))) {

      stack->fIsObjStarted = kTRUE;

      fJsonrCnt++; // count object, but do not keep reference

      stack = JsonStartObjectWrite(cl, sinfo);
   } else {
      stack = PushStack(0);
   }

   stack->fInfo = sinfo;
   stack->fIsStreamerInfo = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and ReadBuffer functions
/// and decrease level in json structure.

void TBufferJSON::DecrementLevel(TVirtualStreamerInfo *info)
{
   if (gDebug > 2)
      Info("DecrementLevel", "Class: %s", (info ? info->GetClass()->GetName() : "custom"));

   TJSONStackObj *stack = Stack();

   if (stack->IsStreamerElement()) {

      if (IsWriting()) {
         if (gDebug > 3)
            Info("DecrementLevel", "    Perform post-processing elem: %s", stack->fElem->GetName());

         PerformPostProcessing(stack);
      }

      stack = PopStack(); // remove stack of last element
   }

   if (stack->fInfo != (TStreamerInfo *)info)
      Error("DecrementLevel", "    Mismatch of streamer info");

   PopStack(); // back from data of stack info

   if (gDebug > 3)
      Info("DecrementLevel", "Class: %s done", (info ? info->GetClass()->GetName() : "custom"));
}

////////////////////////////////////////////////////////////////////////////////
/// Return current streamer info element

TVirtualStreamerInfo *TBufferJSON::GetInfo()
{
   return Stack()->fInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and ReadBuffer functions
/// and add/verify next element of json structure
/// This calls allows separate data, correspondent to one class member, from another

void TBufferJSON::SetStreamerElementNumber(TStreamerElement *elem, Int_t comp_type)
{
   if (gDebug > 3)
      Info("SetStreamerElementNumber", "Element name %s", elem->GetName());

   WorkWithElement(elem, comp_type);
}

////////////////////////////////////////////////////////////////////////////////
/// This is call-back from streamer which indicates
/// that class member will be streamed
/// Name of element used in JSON

void TBufferJSON::WorkWithElement(TStreamerElement *elem, Int_t)
{
   TJSONStackObj *stack = Stack();
   if (!stack) {
      Error("WorkWithElement", "stack is empty");
      return;
   }

   if (gDebug > 0)
      Info("WorkWithElement", "    Start element %s type %d typename %s", elem ? elem->GetName() : "---",
           elem ? elem->GetType() : -1, elem ? elem->GetTypeName() : "---");

   if (stack->IsStreamerElement()) {
      // this is post processing

      if (IsWriting()) {
         if (gDebug > 3)
            Info("WorkWithElement", "    Perform post-processing elem: %s", stack->fElem->GetName());
         PerformPostProcessing(stack);
      }

      stack = PopStack(); // go level back
   }

   fValue.Clear();

   if (!stack) {
      Error("WorkWithElement", "Lost of stack");
      return;
   }

   TStreamerInfo *info = stack->fInfo;
   if (!stack->IsStreamerInfo()) {
      Error("WorkWithElement", "Problem in Inc/Dec level");
      return;
   }

   Int_t number = info ? info->GetElements()->IndexOf(elem) : -1;

   if (!elem) {
      Error("WorkWithElement", "streamer info returns elem = nullptr");
      return;
   }

   TClass *base_class = elem->IsBase() ? elem->GetClassPointer() : nullptr;

   stack = PushStack(0, stack->fNode);
   stack->fElem = elem;
   stack->fIsElemOwner = (number < 0);

   JsonStartElement(elem, base_class);

   if (base_class && IsReading())
      stack->fClVersion = base_class->GetClassVersion();

   if ((elem->GetType() == TStreamerInfo::kOffsetL + TStreamerInfo::kStreamLoop) && (elem->GetArrayDim() > 0)) {
      // array of array, start handling here
      stack->fIndx = std::make_unique<TArrayIndexProducer>(elem, -1, fArraySepar.Data());
      if (IsWriting())
         AppendOutput(stack->fIndx->GetBegin());
   }

   if (IsReading() && (elem->GetType() > TStreamerInfo::kOffsetP) && (elem->GetType() < TStreamerInfo::kOffsetP + 20)) {
      // reading of such array begins with reading of single Char_t value
      // it indicates if array should be read or not
      stack->PushIntValue(stack->IsJsonString() || (stack->IsJsonArray() > 0) ? 1 : 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Should be called in the beginning of custom class streamer.
/// Informs buffer data about class which will be streamed now.
///
/// ClassBegin(), ClassEnd() and ClassMember() should be used in
/// custom class streamers to specify which kind of data are
/// now streamed. Such information is used to correctly
/// convert class data to JSON. Without that functions calls
/// classes with custom streamers cannot be used with TBufferJSON

void TBufferJSON::ClassBegin(const TClass *cl, Version_t)
{
   WorkWithClass(nullptr, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Should be called at the end of custom streamer
/// See TBufferJSON::ClassBegin for more details

void TBufferJSON::ClassEnd(const TClass *)
{
   DecrementLevel(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Method indicates name and typename of class member,
/// which should be now streamed in custom streamer
/// Following combinations are supported:
/// 1. name = "ClassName", typeName = 0 or typename==ClassName
///    This is a case, when data of parent class "ClassName" should be streamed.
///     For instance, if class directly inherited from TObject, custom
///     streamer should include following code:
/// ~~~{.cpp}
///       b.ClassMember("TObject");
///       TObject::Streamer(b);
/// ~~~
/// 2. Basic data type
/// ~~~{.cpp}
///      b.ClassMember("fInt","Int_t");
///      b >> fInt;
/// ~~~
/// 3. Array of basic data types
/// ~~~{.cpp}
///      b.ClassMember("fArr","Int_t", 5);
///      b.ReadFastArray(fArr, 5);
/// ~~~
/// 4. Object as data member
/// ~~~{.cpp}
///      b.ClassMember("fName","TString");
///      fName.Streamer(b);
/// ~~~
/// 5. Pointer on object as data member
/// ~~~{.cpp}
///      b.ClassMember("fObj","TObject*");
///      b.StreamObject(fObj);
/// ~~~
///
/// arrsize1 and arrsize2 arguments (when specified) indicate first and
/// second dimension of array. Can be used for array of basic types.
/// See ClassBegin() method for more details.

void TBufferJSON::ClassMember(const char *name, const char *typeName, Int_t arrsize1, Int_t arrsize2)
{
   if (!typeName)
      typeName = name;

   if (!name || (strlen(name) == 0)) {
      Error("ClassMember", "Invalid member name");
      return;
   }

   TString tname = typeName;

   Int_t typ_id = -1;

   if (strcmp(typeName, "raw:data") == 0)
      typ_id = TStreamerInfo::kMissing;

   if (typ_id < 0) {
      TDataType *dt = gROOT->GetType(typeName);
      if (dt && (dt->GetType() > 0) && (dt->GetType() < 20))
         typ_id = dt->GetType();
   }

   if (typ_id < 0)
      if (strcmp(name, typeName) == 0) {
         TClass *cl = TClass::GetClass(tname.Data());
         if (cl)
            typ_id = TStreamerInfo::kBase;
      }

   if (typ_id < 0) {
      Bool_t isptr = kFALSE;
      if (tname[tname.Length() - 1] == '*') {
         tname.Resize(tname.Length() - 1);
         isptr = kTRUE;
      }
      TClass *cl = TClass::GetClass(tname.Data());
      if (!cl) {
         Error("ClassMember", "Invalid class specifier %s", typeName);
         return;
      }

      if (cl->IsTObject())
         typ_id = isptr ? TStreamerInfo::kObjectp : TStreamerInfo::kObject;
      else
         typ_id = isptr ? TStreamerInfo::kAnyp : TStreamerInfo::kAny;

      if ((cl == TString::Class()) && !isptr)
         typ_id = TStreamerInfo::kTString;
   }

   TStreamerElement *elem = nullptr;

   if (typ_id == TStreamerInfo::kMissing) {
      elem = new TStreamerElement(name, "title", 0, typ_id, "raw:data");
   } else if (typ_id == TStreamerInfo::kBase) {
      TClass *cl = TClass::GetClass(tname.Data());
      if (cl) {
         TStreamerBase *b = new TStreamerBase(tname.Data(), "title", 0);
         b->SetBaseVersion(cl->GetClassVersion());
         elem = b;
      }
   } else if ((typ_id > 0) && (typ_id < 20)) {
      elem = new TStreamerBasicType(name, "title", 0, typ_id, typeName);
   } else if ((typ_id == TStreamerInfo::kObject) || (typ_id == TStreamerInfo::kTObject) ||
              (typ_id == TStreamerInfo::kTNamed)) {
      elem = new TStreamerObject(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kObjectp) {
      elem = new TStreamerObjectPointer(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kAny) {
      elem = new TStreamerObjectAny(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kAnyp) {
      elem = new TStreamerObjectAnyPointer(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kTString) {
      elem = new TStreamerString(name, "title", 0);
   }

   if (!elem) {
      Error("ClassMember", "Invalid combination name = %s type = %s", name, typeName);
      return;
   }

   if (arrsize1 > 0) {
      elem->SetArrayDim(arrsize2 > 0 ? 2 : 1);
      elem->SetMaxIndex(0, arrsize1);
      if (arrsize2 > 0)
         elem->SetMaxIndex(1, arrsize2);
   }

   // we indicate that there is no streamerinfo
   WorkWithElement(elem, -1);
}

////////////////////////////////////////////////////////////////////////////////
/// Function is converts TObject and TString structures to more compact representation

void TBufferJSON::PerformPostProcessing(TJSONStackObj *stack, const TClass *obj_cl)
{
   if (stack->fIsPostProcessed)
      return;

   const TStreamerElement *elem = stack->fElem;

   if (!elem && !obj_cl)
      return;

   stack->fIsPostProcessed = kTRUE;

   // when element was written as separate object, close only braces and exit
   if (stack->fIsObjStarted) {
      AppendOutput("", "}");
      return;
   }

   Bool_t isTObject(kFALSE), isTRef(kFALSE), isTString(kFALSE), isSTLstring(kFALSE), isOffsetPArray(kFALSE),
      isTArray(kFALSE);

   if (obj_cl) {
      if (obj_cl == TObject::Class())
         isTObject = kTRUE;
      else if (obj_cl == TRef::Class())
         isTRef = kTRUE;
      else
         return;
   } else {
      const char *typname = elem->IsBase() ? elem->GetName() : elem->GetTypeName();
      isTObject = (elem->GetType() == TStreamerInfo::kTObject) || (strcmp("TObject", typname) == 0);
      isTString = elem->GetType() == TStreamerInfo::kTString;
      isSTLstring = elem->GetType() == TStreamerInfo::kSTLstring;
      isOffsetPArray = (elem->GetType() > TStreamerInfo::kOffsetP) && (elem->GetType() < TStreamerInfo::kOffsetP + 20);
      isTArray = (strncmp("TArray", typname, 6) == 0);
   }

   if (isTString || isSTLstring) {
      // just remove all kind of string length information

      if (gDebug > 3)
         Info("PerformPostProcessing", "reformat string value = '%s'", fValue.Data());

      stack->fValues.clear();
   } else if (isOffsetPArray) {
      // basic array with [fN] comment

      if (stack->fValues.empty() && (fValue == "0")) {
         fValue = "[]";
      } else if ((stack->fValues.size() == 1) && (stack->fValues[0] == "1")) {
         stack->fValues.clear();
      } else {
         Error("PerformPostProcessing", "Wrong values for kOffsetP element %s", (elem ? elem->GetName() : "---"));
         stack->fValues.clear();
         fValue = "[]";
      }
   } else if (isTObject || isTRef) {
      // complex workaround for TObject/TRef streamer
      // would be nice if other solution can be found
      // Here is not supported TRef on TRef (double reference)

      Int_t cnt = stack->fValues.size();
      if (fValue.Length() > 0)
         cnt++;

      if (cnt < 2 || cnt > 3) {
         if (gDebug > 0)
            Error("PerformPostProcessing", "When storing TObject/TRef, strange number of items %d", cnt);
         AppendOutput(stack->NextMemberSeparator(), "\"dummy\"");
         AppendOutput(fSemicolon.Data());
      } else {
         AppendOutput(stack->NextMemberSeparator(), "\"fUniqueID\"");
         AppendOutput(fSemicolon.Data());
         AppendOutput(stack->fValues[0].c_str());
         AppendOutput(stack->NextMemberSeparator(), "\"fBits\"");
         AppendOutput(fSemicolon.Data());
         auto tbits = std::atol((stack->fValues.size() > 1) ? stack->fValues[1].c_str() : fValue.Data());
         AppendOutput(std::to_string(tbits & ~TObject::kNotDeleted & ~TObject::kIsOnHeap).c_str());
         if (cnt == 3) {
            AppendOutput(stack->NextMemberSeparator(), "\"fPID\"");
            AppendOutput(fSemicolon.Data());
            AppendOutput((stack->fValues.size() > 2) ? stack->fValues[2].c_str() : fValue.Data());
         }

         stack->fValues.clear();
         fValue.Clear();
         return;
      }

   } else if (isTArray) {
      // for TArray one deletes complete stack
      stack->fValues.clear();
   }

   if (elem && elem->IsBase() && (fValue.Length() == 0)) {
      // here base class data already completely stored
      return;
   }

   if (!stack->fValues.empty()) {
      // append element blob data just as abstract array, user is responsible to decode it
      AppendOutput("[");
      for (auto &blob: stack->fValues) {
         AppendOutput(blob.c_str());
         AppendOutput(fArraySepar.Data());
      }
   }

   if (fValue.Length() == 0) {
      AppendOutput("null");
   } else {
      AppendOutput(fValue.Data());
      fValue.Clear();
   }

   if (!stack->fValues.empty())
      AppendOutput("]");
}

////////////////////////////////////////////////////////////////////////////////
/// suppressed function of TBuffer

TClass *TBufferJSON::ReadClass(const TClass *, UInt_t *)
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// suppressed function of TBuffer

void TBufferJSON::WriteClass(const TClass *) {}

////////////////////////////////////////////////////////////////////////////////
/// read version value from buffer

Version_t TBufferJSON::ReadVersion(UInt_t *start, UInt_t *bcnt, const TClass *cl)
{
   Version_t res = cl ? cl->GetClassVersion() : 0;

   if (start)
      *start = 0;
   if (bcnt)
      *bcnt = 0;

   if (!cl && Stack()->fClVersion) {
      res = Stack()->fClVersion;
      Stack()->fClVersion = 0;
   }

   if (gDebug > 3)
      Info("ReadVersion", "Result: %d Class: %s", res, (cl ? cl->GetName() : "---"));

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Ignored in TBufferJSON

UInt_t TBufferJSON::WriteVersion(const TClass * /*cl*/, Bool_t /* useBcnt */)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from buffer. Only used from TBuffer

void *TBufferJSON::ReadObjectAny(const TClass *expectedClass)
{
   if (gDebug > 2)
      Info("ReadObjectAny", "From current JSON node");
   void *res = JsonReadObject(nullptr, expectedClass);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Skip any kind of object from buffer

void TBufferJSON::SkipObjectAny() {}

////////////////////////////////////////////////////////////////////////////////
/// Write object to buffer. Only used from TBuffer

void TBufferJSON::WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse)
{
   if (gDebug > 3)
      Info("WriteObjectClass", "Class %s", (actualClass ? actualClass->GetName() : " null"));

   JsonWriteObject(actualObjStart, actualClass, cacheReuse);
}

////////////////////////////////////////////////////////////////////////////////
/// If value exists, push in the current stack for post-processing

void TBufferJSON::JsonPushValue()
{
   if (fValue.Length() > 0)
      Stack()->PushValue(fValue);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

Int_t TBufferJSON::ReadArray(Bool_t *&b)
{
   return JsonReadArray(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer

Int_t TBufferJSON::ReadArray(Char_t *&c)
{
   return JsonReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

Int_t TBufferJSON::ReadArray(UChar_t *&c)
{
   return JsonReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

Int_t TBufferJSON::ReadArray(Short_t *&h)
{
   return JsonReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

Int_t TBufferJSON::ReadArray(UShort_t *&h)
{
   return JsonReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

Int_t TBufferJSON::ReadArray(Int_t *&i)
{
   return JsonReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

Int_t TBufferJSON::ReadArray(UInt_t *&i)
{
   return JsonReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

Int_t TBufferJSON::ReadArray(Long_t *&l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

Int_t TBufferJSON::ReadArray(ULong_t *&l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

Int_t TBufferJSON::ReadArray(Long64_t *&l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

Int_t TBufferJSON::ReadArray(ULong64_t *&l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

Int_t TBufferJSON::ReadArray(Float_t *&f)
{
   return JsonReadArray(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

Int_t TBufferJSON::ReadArray(Double_t *&d)
{
   return JsonReadArray(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Read static array from JSON - not used

template <typename T>
R__ALWAYS_INLINE Int_t TBufferJSON::JsonReadArray(T *value)
{
   Info("ReadArray", "Not implemented");
   return value ? 1 : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

Int_t TBufferJSON::ReadStaticArray(Bool_t *b)
{
   return JsonReadArray(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer

Int_t TBufferJSON::ReadStaticArray(Char_t *c)
{
   return JsonReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

Int_t TBufferJSON::ReadStaticArray(UChar_t *c)
{
   return JsonReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

Int_t TBufferJSON::ReadStaticArray(Short_t *h)
{
   return JsonReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

Int_t TBufferJSON::ReadStaticArray(UShort_t *h)
{
   return JsonReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

Int_t TBufferJSON::ReadStaticArray(Int_t *i)
{
   return JsonReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

Int_t TBufferJSON::ReadStaticArray(UInt_t *i)
{
   return JsonReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

Int_t TBufferJSON::ReadStaticArray(Long_t *l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

Int_t TBufferJSON::ReadStaticArray(ULong_t *l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

Int_t TBufferJSON::ReadStaticArray(Long64_t *l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

Int_t TBufferJSON::ReadStaticArray(ULong64_t *l)
{
   return JsonReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

Int_t TBufferJSON::ReadStaticArray(Float_t *f)
{
   return JsonReadArray(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

Int_t TBufferJSON::ReadStaticArray(Double_t *d)
{
   return JsonReadArray(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to read array from the JSON

template <typename T>
R__ALWAYS_INLINE void TBufferJSON::JsonReadFastArray(T *arr, Int_t arrsize, bool asstring)
{
   if (!arr || (arrsize <= 0))
      return;
   nlohmann::json *json = Stack()->fNode;
   if (gDebug > 2)
      Info("ReadFastArray", "Reading array sz %d from JSON %s", arrsize, json->dump().substr(0, 30).c_str());
   auto indexes = Stack()->MakeReadIndexes();
   if (indexes) { /* at least two dims */
      TArrayI &indx = indexes->GetIndices();
      Int_t lastdim = indx.GetSize() - 1;
      if (indexes->TotalLength() != arrsize)
         Error("ReadFastArray", "Mismatch %d-dim array sizes %d %d", lastdim + 1, arrsize, (int)indexes->TotalLength());
      for (int cnt = 0; cnt < arrsize; ++cnt) {
         nlohmann::json *elem = &(json->at(indx[0]));
         for (int k = 1; k < lastdim; ++k)
            elem = &((*elem)[indx[k]]);
         arr[cnt] = (asstring && elem->is_string()) ? elem->get<std::string>()[indx[lastdim]] : (*elem)[indx[lastdim]].get<T>();
         indexes->NextSeparator();
      }
   } else if (asstring && json->is_string()) {
      std::string str = json->get<std::string>();
      for (int cnt = 0; cnt < arrsize; ++cnt)
         arr[cnt] = (cnt < (int)str.length()) ? str[cnt] : 0;
   } else if (json->is_object() && (json->count("$arr") == 1)) {
      if (json->at("len").get<int>() != arrsize)
         Error("ReadFastArray", "Mismatch compressed array size %d %d", arrsize, json->at("len").get<int>());

      for (int cnt = 0; cnt < arrsize; ++cnt)
         arr[cnt] = 0;

      if (json->count("b") == 1) {
         auto base64 = json->at("b").get<std::string>();

         int offset = (json->count("o") == 1) ? json->at("o").get<int>() : 0;

         // TODO: provide TBase64::Decode with direct write into target buffer
         auto decode = TBase64::Decode(base64.c_str());

         if (arrsize * (long) sizeof(T) < (offset + decode.Length())) {
            Error("ReadFastArray", "Base64 data %ld larger than target array size %ld", (long) decode.Length() + offset, (long) (arrsize*sizeof(T)));
         } else if ((sizeof(T) > 1) && (decode.Length() % sizeof(T) != 0)) {
            Error("ReadFastArray", "Base64 data size %ld not matches with element size %ld", (long) decode.Length(), (long) sizeof(T));
         } else {
            memcpy((char *) arr + offset, decode.Data(), decode.Length());
         }
         return;
      }

      int p = 0, id = 0;
      std::string idname = "", pname, vname, nname;
      while (p < arrsize) {
         pname = std::string("p") + idname;
         if (json->count(pname) == 1)
            p = json->at(pname).get<int>();
         vname = std::string("v") + idname;
         if (json->count(vname) != 1)
            break;
         nlohmann::json &v = json->at(vname);
         if (v.is_array()) {
            for (unsigned sub = 0; sub < v.size(); ++sub)
               arr[p++] = v[sub].get<T>();
         } else {
            nname = std::string("n") + idname;
            unsigned ncopy = (json->count(nname) == 1) ? json->at(nname).get<unsigned>() : 1;
            for (unsigned sub = 0; sub < ncopy; ++sub)
               arr[p++] = v.get<T>();
         }
         idname = std::to_string(++id);
      }
   } else {
      if ((int)json->size() != arrsize)
         Error("ReadFastArray", "Mismatch array sizes %d %d", arrsize, (int)json->size());
      for (int cnt = 0; cnt < arrsize; ++cnt)
         arr[cnt] = json->at(cnt).get<T>();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Bool_t from buffer

void TBufferJSON::ReadFastArray(Bool_t *b, Int_t n)
{
   JsonReadFastArray(b, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Char_t from buffer

void TBufferJSON::ReadFastArray(Char_t *c, Int_t n)
{
   JsonReadFastArray(c, n, true);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Char_t from buffer

void TBufferJSON::ReadFastArrayString(Char_t *c, Int_t n)
{
   JsonReadFastArray(c, n, true);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of UChar_t from buffer

void TBufferJSON::ReadFastArray(UChar_t *c, Int_t n)
{
   JsonReadFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Short_t from buffer

void TBufferJSON::ReadFastArray(Short_t *h, Int_t n)
{
   JsonReadFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of UShort_t from buffer

void TBufferJSON::ReadFastArray(UShort_t *h, Int_t n)
{
   JsonReadFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Int_t from buffer

void TBufferJSON::ReadFastArray(Int_t *i, Int_t n)
{
   JsonReadFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of UInt_t from buffer

void TBufferJSON::ReadFastArray(UInt_t *i, Int_t n)
{
   JsonReadFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Long_t from buffer

void TBufferJSON::ReadFastArray(Long_t *l, Int_t n)
{
   JsonReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of ULong_t from buffer

void TBufferJSON::ReadFastArray(ULong_t *l, Int_t n)
{
   JsonReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Long64_t from buffer

void TBufferJSON::ReadFastArray(Long64_t *l, Int_t n)
{
   JsonReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of ULong64_t from buffer

void TBufferJSON::ReadFastArray(ULong64_t *l, Int_t n)
{
   JsonReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Float_t from buffer

void TBufferJSON::ReadFastArray(Float_t *f, Int_t n)
{
   JsonReadFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// read array of Double_t from buffer

void TBufferJSON::ReadFastArray(Double_t *d, Int_t n)
{
   JsonReadFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read an array of 'n' objects from the I/O buffer.
/// Stores the objects read starting at the address 'start'.
/// The objects in the array are assume to be of class 'cl'.
/// Copied code from TBufferFile

void TBufferJSON::ReadFastArray(void *start, const TClass *cl, Int_t n, TMemberStreamer * /* streamer */,
                                const TClass * /* onFileClass */)
{
   if (gDebug > 1)
      Info("ReadFastArray", "void* n:%d cl:%s", n, cl->GetName());

   //   if (streamer) {
   //      Info("ReadFastArray", "(void*) Calling streamer - not handled correctly");
   //      streamer->SetOnFileClass(onFileClass);
   //      (*streamer)(*this, start, 0);
   //      return;
   //   }

   int objectSize = cl->Size();
   char *obj = (char *)start;

   TJSONStackObj *stack = Stack();
   nlohmann::json *topnode = stack->fNode, *subnode = topnode;
   if (stack->fIndx)
      subnode = stack->fIndx->ExtractNode(topnode);

   TArrayIndexProducer indexes(stack->fElem, n, "");

   if (gDebug > 1)
      Info("ReadFastArray", "Indexes ndim:%d totallen:%d", indexes.NumDimensions(), indexes.TotalLength());

   for (Int_t j = 0; j < n; j++, obj += objectSize) {

      stack->fNode = indexes.ExtractNode(subnode);

      JsonReadObject(obj, cl);
   }

   // restore top node - show we use stack here?
   stack->fNode = topnode;
}

////////////////////////////////////////////////////////////////////////////////
/// redefined here to avoid warning message from gcc

void TBufferJSON::ReadFastArray(void **start, const TClass *cl, Int_t n, Bool_t isPreAlloc,
                                TMemberStreamer * /* streamer */, const TClass * /* onFileClass */)
{
   if (gDebug > 1)
      Info("ReadFastArray", "void** n:%d cl:%s prealloc:%s", n, cl->GetName(), (isPreAlloc ? "true" : "false"));

   //   if (streamer) {
   //      Info("ReadFastArray", "(void**) Calling streamer - not handled correctly");
   //      if (isPreAlloc) {
   //         for (Int_t j = 0; j < n; j++) {
   //            if (!start[j])
   //               start[j] = cl->New();
   //         }
   //      }
   //      streamer->SetOnFileClass(onFileClass);
   //      (*streamer)(*this, (void *)start, 0);
   //      return;
   //   }

   TJSONStackObj *stack = Stack();
   nlohmann::json *topnode = stack->fNode, *subnode = topnode;
   if (stack->fIndx)
      subnode = stack->fIndx->ExtractNode(topnode);

   TArrayIndexProducer indexes(stack->fElem, n, "");

   for (Int_t j = 0; j < n; j++) {

      stack->fNode = indexes.ExtractNode(subnode);

      if (!isPreAlloc) {
         void *old = start[j];
         start[j] = JsonReadObject(nullptr, cl);
         if (old && old != start[j] && TStreamerInfo::CanDelete())
            (const_cast<TClass *>(cl))->Destructor(old, kFALSE); // call delete and destruct
      } else {
         if (!start[j])
            start[j] = (const_cast<TClass *>(cl))->New();
         JsonReadObject(start[j], cl);
      }
   }

   stack->fNode = topnode;
}

template <typename T>
R__ALWAYS_INLINE void TBufferJSON::JsonWriteArrayCompress(const T *vname, Int_t arrsize, const char *typname)
{
   bool is_base64 = Stack()->fBase64 || (fArrayCompact == kBase64);

   if (!is_base64 && ((fArrayCompact == 0) || (arrsize < 6))) {
      fValue.Append("[");
      for (Int_t indx = 0; indx < arrsize; indx++) {
         if (indx > 0)
            fValue.Append(fArraySepar.Data());
         JsonWriteBasic(vname[indx]);
      }
      fValue.Append("]");
   } else if (is_base64 && !arrsize) {
      fValue.Append("[]");
   } else {
      fValue.Append("{");
      fValue.Append(TString::Format("\"$arr\":\"%s\"%s\"len\":%d", typname, fArraySepar.Data(), arrsize));
      Int_t aindx(0), bindx(arrsize);
      while ((aindx < arrsize) && (vname[aindx] == 0))
         aindx++;
      while ((aindx < bindx) && (vname[bindx - 1] == 0))
         bindx--;

      if (is_base64) {
         // small initial offset makes no sense - JSON code is large then size gain
         if ((aindx * sizeof(T) < 5) && (aindx < bindx))
            aindx = 0;

         if ((aindx > 0) && (aindx < bindx))
            fValue.Append(TString::Format("%s\"o\":%ld", fArraySepar.Data(), (long) (aindx * (int) sizeof(T))));

         fValue.Append(fArraySepar);
         fValue.Append("\"b\":\"");

         if (aindx < bindx)
            fValue.Append(TBase64::Encode((const char *) (vname + aindx), (bindx - aindx) * sizeof(T)));

         fValue.Append("\"");
      } else if (aindx < bindx) {
         TString suffix("");
         Int_t p(aindx), suffixcnt(-1), lastp(0);
         while (p < bindx) {
            if (vname[p] == 0) {
               p++;
               continue;
            }
            Int_t p0(p++), pp(0), nsame(1);
            if (fArrayCompact != kSameSuppression) {
               pp = bindx;
               p = bindx + 1;
               nsame = 0;
            }
            for (; p <= bindx; ++p) {
               if ((p < bindx) && (vname[p] == vname[p - 1])) {
                  nsame++;
                  continue;
               }
               if (vname[p - 1] == 0) {
                  if (nsame > 9) {
                     nsame = 0;
                     break;
                  }
               } else if (nsame > 5) {
                  if (pp) {
                     p = pp;
                     nsame = 0;
                  } else
                     pp = p;
                  break;
               }
               pp = p;
               nsame = 1;
            }
            if (pp <= p0)
               continue;
            if (++suffixcnt > 0)
               suffix.Form("%d", suffixcnt);
            if (p0 != lastp)
               fValue.Append(TString::Format("%s\"p%s\":%d", fArraySepar.Data(), suffix.Data(), p0));
            lastp = pp; /* remember cursor, it may be the same */
            fValue.Append(TString::Format("%s\"v%s\":", fArraySepar.Data(), suffix.Data()));
            if ((nsame > 1) || (pp - p0 == 1)) {
               JsonWriteBasic(vname[p0]);
               if (nsame > 1)
                  fValue.Append(TString::Format("%s\"n%s\":%d", fArraySepar.Data(), suffix.Data(), nsame));
            } else {
               fValue.Append("[");
               for (Int_t indx = p0; indx < pp; indx++) {
                  if (indx > p0)
                     fValue.Append(fArraySepar.Data());
                  JsonWriteBasic(vname[indx]);
               }
               fValue.Append("]");
            }
         }
      }
      fValue.Append("}");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Bool_t to buffer

void TBufferJSON::WriteArray(const Bool_t *b, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(b, n, "Bool");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer

void TBufferJSON::WriteArray(const Char_t *c, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(c, n, "Int8");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UChar_t to buffer

void TBufferJSON::WriteArray(const UChar_t *c, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(c, n, "Uint8");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Short_t to buffer

void TBufferJSON::WriteArray(const Short_t *h, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(h, n, "Int16");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UShort_t to buffer

void TBufferJSON::WriteArray(const UShort_t *h, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(h, n, "Uint16");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Int_ to buffer

void TBufferJSON::WriteArray(const Int_t *i, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(i, n, "Int32");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UInt_t to buffer

void TBufferJSON::WriteArray(const UInt_t *i, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(i, n, "Uint32");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long_t to buffer

void TBufferJSON::WriteArray(const Long_t *l, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(l, n, "Int64");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong_t to buffer

void TBufferJSON::WriteArray(const ULong_t *l, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(l, n, "Uint64");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long64_t to buffer

void TBufferJSON::WriteArray(const Long64_t *l, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(l, n, "Int64");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong64_t to buffer

void TBufferJSON::WriteArray(const ULong64_t *l, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(l, n, "Uint64");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float_t to buffer

void TBufferJSON::WriteArray(const Float_t *f, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(f, n, "Float32");
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double_t to buffer

void TBufferJSON::WriteArray(const Double_t *d, Int_t n)
{
   JsonPushValue();
   JsonWriteArrayCompress(d, n, "Float64");
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to write array of arbitrary dimensions
/// Different methods can be used for store last array dimension -
/// either JsonWriteArrayCompress<T>() or JsonWriteConstChar()

template <typename T>
R__ALWAYS_INLINE void TBufferJSON::JsonWriteFastArray(const T *arr, Int_t arrsize, const char *typname,
                                                      void (TBufferJSON::*method)(const T *, Int_t, const char *))
{
   JsonPushValue();
   if (arrsize <= 0) { /*fJsonrCnt++;*/
      fValue.Append("[]");
      return;
   }

   TStreamerElement *elem = Stack()->fElem;
   if (elem && (elem->GetArrayDim() > 1) && (elem->GetArrayLength() == arrsize)) {
      TArrayI indexes(elem->GetArrayDim() - 1);
      indexes.Reset(0);
      Int_t cnt = 0, shift = 0, len = elem->GetMaxIndex(indexes.GetSize());
      while (cnt >= 0) {
         if (indexes[cnt] >= elem->GetMaxIndex(cnt)) {
            fValue.Append("]");
            indexes[cnt--] = 0;
            if (cnt >= 0)
               indexes[cnt]++;
            continue;
         }
         fValue.Append(indexes[cnt] == 0 ? "[" : fArraySepar.Data());
         if (++cnt == indexes.GetSize()) {
            (*this.*method)((arr + shift), len, typname);
            indexes[--cnt]++;
            shift += len;
         }
      }
   } else {
      (*this.*method)(arr, arrsize, typname);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Bool_t to buffer

void TBufferJSON::WriteFastArray(const Bool_t *b, Int_t n)
{
   JsonWriteFastArray(b, n, "Bool", &TBufferJSON::JsonWriteArrayCompress<Bool_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer
///
/// Normally written as JSON string, but if string includes \0 in the middle
/// or some special characters, uses regular array. From array size 1000 it
/// will be automatically converted into base64 coding

void TBufferJSON::WriteFastArray(const Char_t *c, Int_t n)
{
   Bool_t need_blob = false;
   Bool_t has_zero = false;
   for (int i=0;i<n;++i) {
      if (!c[i]) {
         has_zero = true; // might be terminal '\0'
      } else if (has_zero || !isprint(c[i])) {
         need_blob = true;
         break;
      }
   }

   if (need_blob && (n >= 1000) && (!Stack()->fElem || (Stack()->fElem->GetArrayDim() < 2)))
      Stack()->fBase64 = true;

   JsonWriteFastArray(c, n, "Int8", need_blob ? &TBufferJSON::JsonWriteArrayCompress<Char_t> : &TBufferJSON::JsonWriteConstChar);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer

void TBufferJSON::WriteFastArrayString(const Char_t *c, Int_t n)
{
   JsonWriteFastArray(c, n, "Int8", &TBufferJSON::JsonWriteConstChar);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UChar_t to buffer

void TBufferJSON::WriteFastArray(const UChar_t *c, Int_t n)
{
   JsonWriteFastArray(c, n, "Uint8", &TBufferJSON::JsonWriteArrayCompress<UChar_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Short_t to buffer

void TBufferJSON::WriteFastArray(const Short_t *h, Int_t n)
{
   JsonWriteFastArray(h, n, "Int16", &TBufferJSON::JsonWriteArrayCompress<Short_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UShort_t to buffer

void TBufferJSON::WriteFastArray(const UShort_t *h, Int_t n)
{
   JsonWriteFastArray(h, n, "Uint16", &TBufferJSON::JsonWriteArrayCompress<UShort_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Int_t to buffer

void TBufferJSON::WriteFastArray(const Int_t *i, Int_t n)
{
   JsonWriteFastArray(i, n, "Int32", &TBufferJSON::JsonWriteArrayCompress<Int_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UInt_t to buffer

void TBufferJSON::WriteFastArray(const UInt_t *i, Int_t n)
{
   JsonWriteFastArray(i, n, "Uint32", &TBufferJSON::JsonWriteArrayCompress<UInt_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long_t to buffer

void TBufferJSON::WriteFastArray(const Long_t *l, Int_t n)
{
   JsonWriteFastArray(l, n, "Int64", &TBufferJSON::JsonWriteArrayCompress<Long_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong_t to buffer

void TBufferJSON::WriteFastArray(const ULong_t *l, Int_t n)
{
   JsonWriteFastArray(l, n, "Uint64", &TBufferJSON::JsonWriteArrayCompress<ULong_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long64_t to buffer

void TBufferJSON::WriteFastArray(const Long64_t *l, Int_t n)
{
   JsonWriteFastArray(l, n, "Int64", &TBufferJSON::JsonWriteArrayCompress<Long64_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong64_t to buffer

void TBufferJSON::WriteFastArray(const ULong64_t *l, Int_t n)
{
   JsonWriteFastArray(l, n, "Uint64", &TBufferJSON::JsonWriteArrayCompress<ULong64_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float_t to buffer

void TBufferJSON::WriteFastArray(const Float_t *f, Int_t n)
{
   JsonWriteFastArray(f, n, "Float32", &TBufferJSON::JsonWriteArrayCompress<Float_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double_t to buffer

void TBufferJSON::WriteFastArray(const Double_t *d, Int_t n)
{
   JsonWriteFastArray(d, n, "Float64", &TBufferJSON::JsonWriteArrayCompress<Double_t>);
}

////////////////////////////////////////////////////////////////////////////////
/// Recall TBuffer function to avoid gcc warning message

void TBufferJSON::WriteFastArray(void *start, const TClass *cl, Int_t n, TMemberStreamer * /* streamer */)
{
   if (gDebug > 2)
      Info("WriteFastArray", "void *start cl:%s n:%d", cl ? cl->GetName() : "---", n);

   //   if (streamer) {
   //      JsonDisablePostprocessing();
   //      (*streamer)(*this, start, 0);
   //      return;
   //   }

   if (n < 0) {
      // special handling of empty StreamLoop
      AppendOutput("null");
      JsonDisablePostprocessing();
   } else {

      char *obj = (char *)start;
      if (!n)
         n = 1;
      int size = cl->Size();

      TArrayIndexProducer indexes(Stack()->fElem, n, fArraySepar.Data());

      if (indexes.IsArray()) {
         JsonDisablePostprocessing();
         AppendOutput(indexes.GetBegin());
      }

      for (Int_t j = 0; j < n; j++, obj += size) {

         if (j > 0)
            AppendOutput(indexes.NextSeparator());

         JsonWriteObject(obj, cl, kFALSE);

         if (indexes.IsArray() && (fValue.Length() > 0)) {
            AppendOutput(fValue.Data());
            fValue.Clear();
         }
      }

      if (indexes.IsArray())
         AppendOutput(indexes.GetEnd());
   }

   if (Stack()->fIndx)
      AppendOutput(Stack()->fIndx->NextSeparator());
}

////////////////////////////////////////////////////////////////////////////////
/// Recall TBuffer function to avoid gcc warning message

Int_t TBufferJSON::WriteFastArray(void **start, const TClass *cl, Int_t n, Bool_t isPreAlloc,
                                  TMemberStreamer * /* streamer */)
{
   if (gDebug > 2)
      Info("WriteFastArray", "void **startp cl:%s n:%d", cl->GetName(), n);

   //   if (streamer) {
   //      JsonDisablePostprocessing();
   //      (*streamer)(*this, (void *)start, 0);
   //      return 0;
   //   }

   if (n <= 0)
      return 0;

   Int_t res = 0;

   TArrayIndexProducer indexes(Stack()->fElem, n, fArraySepar.Data());

   if (indexes.IsArray()) {
      JsonDisablePostprocessing();
      AppendOutput(indexes.GetBegin());
   }

   for (Int_t j = 0; j < n; j++) {

      if (j > 0)
         AppendOutput(indexes.NextSeparator());

      if (!isPreAlloc) {
         res |= WriteObjectAny(start[j], cl);
      } else {
         if (!start[j])
            start[j] = (const_cast<TClass *>(cl))->New();
         // ((TClass*)cl)->Streamer(start[j],*this);
         JsonWriteObject(start[j], cl, kFALSE);
      }

      if (indexes.IsArray() && (fValue.Length() > 0)) {
         AppendOutput(fValue.Data());
         fValue.Clear();
      }
   }

   if (indexes.IsArray())
      AppendOutput(indexes.GetEnd());

   if (Stack()->fIndx)
      AppendOutput(Stack()->fIndx->NextSeparator());

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// stream object to/from buffer

void TBufferJSON::StreamObject(void *obj, const TClass *cl, const TClass * /* onfileClass */)
{
   if (gDebug > 3)
      Info("StreamObject", "Class: %s", (cl ? cl->GetName() : "none"));

   if (IsWriting())
      JsonWriteObject(obj, cl);
   else
      JsonReadObject(obj, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Template function to read basic value from JSON

template <typename T>
R__ALWAYS_INLINE void TBufferJSON::JsonReadBasic(T &value)
{
   value = Stack()->GetStlNode()->get<T>();
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Bool_t value from buffer

void TBufferJSON::ReadBool(Bool_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Char_t value from buffer

void TBufferJSON::ReadChar(Char_t &val)
{
   if (!Stack()->fValues.empty())
      val = (Char_t)Stack()->PopIntValue();
   else
      val = Stack()->GetStlNode()->get<Char_t>();
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UChar_t value from buffer

void TBufferJSON::ReadUChar(UChar_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Short_t value from buffer

void TBufferJSON::ReadShort(Short_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UShort_t value from buffer

void TBufferJSON::ReadUShort(UShort_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Int_t value from buffer

void TBufferJSON::ReadInt(Int_t &val)
{
   if (!Stack()->fValues.empty())
      val = Stack()->PopIntValue();
   else
      JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UInt_t value from buffer

void TBufferJSON::ReadUInt(UInt_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Long_t value from buffer

void TBufferJSON::ReadLong(Long_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads ULong_t value from buffer

void TBufferJSON::ReadULong(ULong_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Long64_t value from buffer

void TBufferJSON::ReadLong64(Long64_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads ULong64_t value from buffer

void TBufferJSON::ReadULong64(ULong64_t &val)
{
   JsonReadBasic(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Float_t value from buffer

void TBufferJSON::ReadFloat(Float_t &val)
{
   nlohmann::json *json = Stack()->GetStlNode();
   if (json->is_null())
      val = std::numeric_limits<Float_t>::quiet_NaN();
   else
      val = json->get<Float_t>();
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Double_t value from buffer

void TBufferJSON::ReadDouble(Double_t &val)
{
   nlohmann::json *json = Stack()->GetStlNode();
   if (json->is_null())
      val = std::numeric_limits<Double_t>::quiet_NaN();
   else
      val = json->get<Double_t>();
}

////////////////////////////////////////////////////////////////////////////////
/// Reads array of characters from buffer

void TBufferJSON::ReadCharP(Char_t *)
{
   Error("ReadCharP", "Not implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Reads a TString

void TBufferJSON::ReadTString(TString &val)
{
   std::string str;
   JsonReadBasic(str);
   val = str.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Reads a std::string

void TBufferJSON::ReadStdString(std::string *val)
{
   JsonReadBasic(*val);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads a char* string

void TBufferJSON::ReadCharStar(char *&s)
{
   std::string str;
   JsonReadBasic(str);

   if (s) {
      delete[] s;
      s = nullptr;
   }

   std::size_t nch = str.length();
   if (nch > 0) {
      s = new char[nch + 1];
      memcpy(s, str.c_str(), nch);
      s[nch] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Bool_t value to buffer

void TBufferJSON::WriteBool(Bool_t b)
{
   JsonPushValue();
   JsonWriteBasic(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Char_t value to buffer

void TBufferJSON::WriteChar(Char_t c)
{
   JsonPushValue();
   JsonWriteBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UChar_t value to buffer

void TBufferJSON::WriteUChar(UChar_t c)
{
   JsonPushValue();
   JsonWriteBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Short_t value to buffer

void TBufferJSON::WriteShort(Short_t h)
{
   JsonPushValue();
   JsonWriteBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UShort_t value to buffer

void TBufferJSON::WriteUShort(UShort_t h)
{
   JsonPushValue();
   JsonWriteBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Int_t value to buffer

void TBufferJSON::WriteInt(Int_t i)
{
   JsonPushValue();
   JsonWriteBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UInt_t value to buffer

void TBufferJSON::WriteUInt(UInt_t i)
{
   JsonPushValue();
   JsonWriteBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Long_t value to buffer

void TBufferJSON::WriteLong(Long_t l)
{
   JsonPushValue();
   JsonWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes ULong_t value to buffer

void TBufferJSON::WriteULong(ULong_t l)
{
   JsonPushValue();
   JsonWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Long64_t value to buffer

void TBufferJSON::WriteLong64(Long64_t l)
{
   JsonPushValue();
   JsonWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes ULong64_t value to buffer

void TBufferJSON::WriteULong64(ULong64_t l)
{
   JsonPushValue();
   JsonWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Float_t value to buffer

void TBufferJSON::WriteFloat(Float_t f)
{
   JsonPushValue();
   JsonWriteBasic(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Double_t value to buffer

void TBufferJSON::WriteDouble(Double_t d)
{
   JsonPushValue();
   JsonWriteBasic(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes array of characters to buffer

void TBufferJSON::WriteCharP(const Char_t *c)
{
   JsonPushValue();

   JsonWriteConstChar(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes a TString

void TBufferJSON::WriteTString(const TString &s)
{
   JsonPushValue();

   JsonWriteConstChar(s.Data(), s.Length());
}

////////////////////////////////////////////////////////////////////////////////
/// Writes a std::string

void TBufferJSON::WriteStdString(const std::string *s)
{
   JsonPushValue();

   if (s)
      JsonWriteConstChar(s->c_str(), s->length());
   else
      JsonWriteConstChar("", 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes a char*

void TBufferJSON::WriteCharStar(char *s)
{
   JsonPushValue();

   JsonWriteConstChar(s);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Char_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Char_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%d", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Short_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Short_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%hd", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Int_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Int_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%d", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Long_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Long_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%ld", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Long64_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Long64_t value)
{
   fValue.Append(std::to_string(value).c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// converts Float_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Float_t value)
{
   if (std::isinf(value)) {
      fValue.Append((value < 0.) ? "-2e308" : "2e308"); // Number.MAX_VALUE is approx 1.79e308
   } else if (std::isnan(value)) {
      fValue.Append("null");
   } else {
      char buf[200];
      ConvertFloat(value, buf, sizeof(buf));
      fValue.Append(buf);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// converts Double_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Double_t value)
{
   if (std::isinf(value)) {
      fValue.Append((value < 0.) ? "-2e308" : "2e308"); // Number.MAX_VALUE is approx 1.79e308
   } else if (std::isnan(value)) {
      fValue.Append("null");
   } else {
      char buf[200];
      ConvertDouble(value, buf, sizeof(buf));
      fValue.Append(buf);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// converts Bool_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(Bool_t value)
{
   fValue.Append(value ? "true" : "false");
}

////////////////////////////////////////////////////////////////////////////////
/// converts UChar_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(UChar_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%u", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts UShort_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(UShort_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%hu", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts UInt_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(UInt_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%u", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts ULong_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(ULong_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%lu", value);
   fValue.Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// converts ULong64_t to string and add to json value buffer

void TBufferJSON::JsonWriteBasic(ULong64_t value)
{
   fValue.Append(std::to_string(value).c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// writes string value, processing all kind of special characters

void TBufferJSON::JsonWriteConstChar(const char *value, Int_t len, const char * /* typname */)
{
   if (!value) {

      fValue.Append("\"\"");

   } else {

      fValue.Append("\"");

      if (len < 0)
         len = strlen(value);

      for (Int_t n = 0; n < len; n++) {
         char c = value[n];
         if (c == 0)
            break;
         switch (c) {
         case '\n': fValue.Append("\\n"); break;
         case '\t': fValue.Append("\\t"); break;
         case '\"': fValue.Append("\\\""); break;
         case '\\': fValue.Append("\\\\"); break;
         case '\b': fValue.Append("\\b"); break;
         case '\f': fValue.Append("\\f"); break;
         case '\r': fValue.Append("\\r"); break;
         case '/': fValue.Append("\\/"); break;
         default:
            if ((c > 31) && (c < 127))
               fValue.Append(c);
            else
               fValue.Append(TString::Format("\\u%04x", (unsigned)c));
         }
      }

      fValue.Append("\"");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read data of base class.

void TBufferJSON::ReadBaseClass(void *start, TStreamerBase *elem)
{
   if (elem->GetClassPointer() == TObject::Class()) {
      JsonReadTObjectMembers((TObject *)start);
   } else {
      TBufferText::ReadBaseClass(start, elem);
   }
}
