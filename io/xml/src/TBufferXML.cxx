// @(#)root/:$Id: 5400e36954e1dc109fcfc306242c30234beb7312 $
// Author: Sergey Linev, Rene Brun  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TBufferXML
\ingroup IO

Class for serializing/deserializing object to/from xml.

It redefines most of TBuffer class function to convert simple types,
array of simple types and objects to/from xml.
Instead of writing a binary data it creates a set of xml structures as
nodes and attributes
TBufferXML class uses streaming mechanism, provided by ROOT system,
therefore most of ROOT and user classes can be stored to xml. There are
limitations for complex objects like TTree, which can not be yet converted to xml.
*/

#include "TBufferXML.h"
#include "Compression.h"
#include "TXMLFile.h"

#include "TObjArray.h"
#include "TROOT.h"
#include "TError.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TDataType.h"
#include "TExMap.h"
#include "TMethodCall.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TProcessID.h"
#include "TFile.h"
#include "TMemberStreamer.h"
#include "TStreamer.h"
#include "TStreamerInfoActions.h"
#include "RZip.h"

#ifdef R__VISUAL_CPLUSPLUS
#define FLong64 "%I64d"
#define FULong64 "%I64u"
#else
#define FLong64 "%lld"
#define FULong64 "%llu"
#endif

ClassImp(TBufferXML);

std::string TBufferXML::fgFloatFmt = "%e";

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TBufferXML::TBufferXML()
   : TBuffer(), TXMLSetup(), fXML(nullptr), fStack(), fVersionBuf(-111), fObjMap(nullptr), fIdArray(nullptr),
     fErrorFlag(0), fCanUseCompact(kFALSE), fExpectedChain(kFALSE), fExpectedBaseClass(0), fCompressLevel(0),
     fIOVersion(3)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates buffer object to serialize/deserialize data to/from xml.
/// Mode should be either TBuffer::kRead or TBuffer::kWrite.

TBufferXML::TBufferXML(TBuffer::EMode mode)
   : TBuffer(mode), TXMLSetup(), fXML(nullptr), fStack(), fVersionBuf(-111), fObjMap(nullptr), fIdArray(nullptr),
     fErrorFlag(0), fCanUseCompact(kFALSE), fExpectedChain(kFALSE), fExpectedBaseClass(0), fCompressLevel(0),
     fIOVersion(3)
{
   fBufSize = 1000000000;

   SetParent(nullptr);
   SetBit(kCannotHandleMemberWiseStreaming);
   // SetBit(kTextBasedStreaming);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates buffer object to serialize/deserialize data to/from xml.
/// This constructor should be used, if data from buffer supposed to be stored in file.
/// Mode should be either TBuffer::kRead or TBuffer::kWrite.

TBufferXML::TBufferXML(TBuffer::EMode mode, TXMLFile *file)
   : TBuffer(mode), TXMLSetup(*file), fXML(nullptr), fStack(), fVersionBuf(-111), fObjMap(nullptr), fIdArray(nullptr),
     fErrorFlag(0), fCanUseCompact(kFALSE), fExpectedChain(kFALSE), fExpectedBaseClass(0), fCompressLevel(0),
     fIOVersion(3)
{
   // this is for the case when StreamerInfo reads elements from
   // buffer as ReadFastArray. When it checks if size of buffer is
   // too small and skip reading. Actually, more improved method should
   // be used here.
   fBufSize = 1000000000;

   SetParent(file);
   SetBit(kCannotHandleMemberWiseStreaming);
   // SetBit(kTextBasedStreaming);
   if (XmlFile()) {
      SetXML(XmlFile()->XML());
      SetCompressionSettings(XmlFile()->GetCompressionSettings());
      SetIOVersion(XmlFile()->GetIOVersion());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy xml buffer.

TBufferXML::~TBufferXML()
{
   if (fObjMap)
      delete fObjMap;
   if (fIdArray)
      delete fIdArray;
   fStack.Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pointer to TXMLFile object.
/// Access to file is necessary to produce unique identifier for object references.

TXMLFile *TBufferXML::XmlFile()
{
   return dynamic_cast<TXMLFile *>(GetParent());
}

////////////////////////////////////////////////////////////////////////////////
/// Converts object, inherited from TObject class, to XML string
/// GenericLayout defines layout choice for XML file
/// UseNamespaces allow XML namespaces.
/// See TXMLSetup class for details

TString TBufferXML::ConvertToXML(const TObject *obj, Bool_t GenericLayout, Bool_t UseNamespaces)
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

   return ConvertToXML(ptr, clActual, GenericLayout, UseNamespaces);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts any type of object to XML string.
/// GenericLayout defines layout choice for XML file
/// UseNamespaces allow XML namespaces.
/// See TXMLSetup class for details

TString TBufferXML::ConvertToXML(const void *obj, const TClass *cl, Bool_t GenericLayout, Bool_t UseNamespaces)
{
   TXMLEngine xml;

   TBufferXML buf(TBuffer::kWrite);
   buf.SetXML(&xml);

   buf.SetXmlLayout(GenericLayout ? TXMLSetup::kGeneralized : TXMLSetup::kSpecialized);
   buf.SetUseNamespaces(UseNamespaces);

   XMLNodePointer_t xmlnode = buf.XmlWriteAny(obj, cl);

   TString res;

   xml.SaveSingleNode(xmlnode, &res);

   xml.FreeNode(xmlnode);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from XML, produced by ConvertToXML() method.
/// If object does not inherit from TObject class, return 0.
/// GenericLayout and UseNamespaces should be the same as in ConvertToXML()

TObject *TBufferXML::ConvertFromXML(const char *str, Bool_t GenericLayout, Bool_t UseNamespaces)
{
   TClass *cl = nullptr;
   void *obj = ConvertFromXMLAny(str, &cl, GenericLayout, UseNamespaces);

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
/// Read object of any class from XML, produced by ConvertToXML() method.
/// If cl!=0, return actual class of object.
/// GenericLayout and UseNamespaces should be the same as in ConvertToXML()

void *TBufferXML::ConvertFromXMLAny(const char *str, TClass **cl, Bool_t GenericLayout, Bool_t UseNamespaces)
{
   TXMLEngine xml;
   TBufferXML buf(TBuffer::kRead);

   buf.SetXML(&xml);

   buf.SetXmlLayout(GenericLayout ? TXMLSetup::kGeneralized : TXMLSetup::kSpecialized);
   buf.SetUseNamespaces(UseNamespaces);

   XMLNodePointer_t xmlnode = xml.ReadSingleNode(str);

   void *obj = buf.XmlReadAny(xmlnode, 0, cl);

   xml.FreeNode(xmlnode);

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert from XML and check if object derived from specified class
/// When possible, cast to given class

void *TBufferXML::ConvertFromXMLChecked(const char *xml, const TClass *expectedClass, Bool_t GenericLayout,
                                        Bool_t UseNamespaces)
{
   TClass *objClass = nullptr;
   void *res = ConvertFromXMLAny(xml, &objClass, GenericLayout, UseNamespaces);

   if (!res || !objClass)
      return nullptr;

   if (objClass == expectedClass)
      return res;

   Int_t offset = objClass->GetBaseClassOffset(expectedClass);
   if (offset < 0) {
      ::Error("TBufferXML::ConvertFromXMLChecked", "expected class %s is not base for read class %s",
              expectedClass->GetName(), objClass->GetName());
      objClass->Destructor(res);
      return nullptr;
   }

   return (char *)res - offset;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert object of any class to xml structures
/// Return pointer on top xml element

XMLNodePointer_t TBufferXML::XmlWriteAny(const void *obj, const TClass *cl)
{
   fErrorFlag = 0;

   if (!fXML)
      return nullptr;

   XMLNodePointer_t res = XmlWriteObject(obj, cl, kTRUE);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Recreate object from xml structure.
/// Return pointer to read object.
/// if (cl!=0) returns pointer to class of object

void *TBufferXML::XmlReadAny(XMLNodePointer_t node, void *obj, TClass **cl)
{
   if (!node)
      return nullptr;

   if (cl)
      *cl = nullptr;

   fErrorFlag = 0;

   if (!fXML)
      return nullptr;

   PushStack(node, kTRUE);

   void *res = XmlReadObject(obj, cl);

   PopStack();

   return res;
}

// TXMLStackObj is used to keep stack of object hierarchy,
// stored in TBuffer. For example, data for parent class(es)
// stored in subnodes, but initial object node will be kept.

class TXMLStackObj : public TObject {
public:
   TXMLStackObj(XMLNodePointer_t node)
      : TObject(), fNode(node), fInfo(0), fElem(0), fElemNumber(0), fCompressedClassNode(kFALSE), fClassNs(0),
        fIsStreamerInfo(kFALSE), fIsElemOwner(kFALSE)
   {
   }

   virtual ~TXMLStackObj()
   {
      if (fIsElemOwner)
         delete fElem;
   }

   Bool_t IsStreamerInfo() const { return fIsStreamerInfo; }

   XMLNodePointer_t fNode;
   TStreamerInfo *fInfo;
   TStreamerElement *fElem;
   Int_t fElemNumber;
   Bool_t fCompressedClassNode;
   XMLNsPointer_t fClassNs;
   Bool_t fIsStreamerInfo;
   Bool_t fIsElemOwner;
};

////////////////////////////////////////////////////////////////////////////////
/// Add new level to xml stack.

TXMLStackObj *TBufferXML::PushStack(XMLNodePointer_t current, Bool_t simple)
{
   if (IsReading() && !simple) {
      current = fXML->GetChild(current);
      fXML->SkipEmpty(current);
   }

   TXMLStackObj *stack = new TXMLStackObj(current);
   fStack.Add(stack);
   return stack;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove one level from xml stack.

TXMLStackObj *TBufferXML::PopStack()
{
   TObject *last = fStack.Last();
   if (last) {
      fStack.Remove(last);
      delete last;
      fStack.Compress();
   }
   return static_cast<TXMLStackObj *>(fStack.Last());
}

////////////////////////////////////////////////////////////////////////////////
/// Return xml stack object of specified depth.

TXMLStackObj *TBufferXML::Stack(Int_t depth)
{
   return (depth <= fStack.GetLast()) ? static_cast<TXMLStackObj *>(fStack.At(fStack.GetLast() - depth)) : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer on current xml node.

XMLNodePointer_t TBufferXML::StackNode()
{
   TXMLStackObj *stack = static_cast<TXMLStackObj *>(fStack.Last());
   return stack ? stack->fNode : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Shift stack node to next.

void TBufferXML::ShiftStack(const char *errinfo)
{
   TXMLStackObj *stack = static_cast<TXMLStackObj *>(fStack.Last());
   if (stack) {
      fXML->ShiftToNext(stack->fNode);
      if (gDebug > 4)
         Info("ShiftStack", "%s to node %s", errinfo, fXML->GetNodeName(stack->fNode));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// See comments for function SetCompressionSettings.

void TBufferXML::SetCompressionAlgorithm(Int_t algorithm)
{
   if (algorithm < 0 || algorithm >= ROOT::kUndefinedCompressionAlgorithm)
      algorithm = 0;
   if (fCompressLevel < 0) {
      // if the level is not defined yet use 1 as a default
      fCompressLevel = 100 * algorithm + 1;
   } else {
      int level = fCompressLevel % 100;
      fCompressLevel = 100 * algorithm + level;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// See comments for function SetCompressionSettings.

void TBufferXML::SetCompressionLevel(Int_t level)
{
   if (level < 0)
      level = 0;
   if (level > 99)
      level = 99;
   if (fCompressLevel < 0) {
      // if the algorithm is not defined yet use 0 as a default
      fCompressLevel = level;
   } else {
      int algorithm = fCompressLevel / 100;
      if (algorithm >= ROOT::kUndefinedCompressionAlgorithm)
         algorithm = 0;
      fCompressLevel = 100 * algorithm + level;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Used to specify the compression level and algorithm.
///
/// See TFile constructor for the details.

void TBufferXML::SetCompressionSettings(Int_t settings)
{
   fCompressLevel = settings;
}

////////////////////////////////////////////////////////////////////////////////
/// Write binary data block from buffer to xml.
/// This data can be produced only by direct call of TBuffer::WriteBuf() functions.

void TBufferXML::XmlWriteBlock(XMLNodePointer_t node)
{
   if (!node || (Length() == 0))
      return;

   const char *src = Buffer();
   int srcSize = Length();

   char *fZipBuffer = 0;

   Int_t compressionLevel = GetCompressionLevel();
   ROOT::ECompressionAlgorithm compressionAlgorithm =
      static_cast<ROOT::ECompressionAlgorithm>(GetCompressionAlgorithm());

   if ((Length() > 512) && (compressionLevel > 0)) {
      int zipBufferSize = Length();
      fZipBuffer = new char[zipBufferSize + 9];
      int dataSize = Length();
      int compressedSize = 0;
      R__zipMultipleAlgorithm(compressionLevel, &dataSize, Buffer(), &zipBufferSize, fZipBuffer, &compressedSize,
                              compressionAlgorithm);
      if (compressedSize > 0) {
         src = fZipBuffer;
         srcSize = compressedSize;
      } else {
         delete[] fZipBuffer;
         fZipBuffer = 0;
      }
   }

   TString res;
   char sbuf[500];
   int block = 0;
   char *tgt = sbuf;
   int srcCnt = 0;

   while (srcCnt++ < srcSize) {
      tgt += sprintf(tgt, " %02x", (unsigned char)*src);
      src++;
      if (block++ == 100) {
         res += sbuf;
         block = 0;
         tgt = sbuf;
      }
   }

   if (block > 0)
      res += sbuf;

   XMLNodePointer_t blocknode = fXML->NewChild(node, 0, xmlio::XmlBlock, res);
   fXML->NewIntAttr(blocknode, xmlio::Size, Length());

   if (fZipBuffer) {
      fXML->NewIntAttr(blocknode, xmlio::Zip, srcSize);
      delete[] fZipBuffer;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read binary block of data from xml.

void TBufferXML::XmlReadBlock(XMLNodePointer_t blocknode)
{
   if (!blocknode)
      return;

   Int_t blockSize = fXML->GetIntAttr(blocknode, xmlio::Size);
   Bool_t blockCompressed = fXML->HasAttr(blocknode, xmlio::Zip);
   char *fUnzipBuffer = 0;

   if (gDebug > 2)
      Info("XmlReadBlock", "Block size = %d, Length = %d, Compressed = %d", blockSize, Length(), blockCompressed);

   if (blockSize > BufferSize())
      Expand(blockSize);

   char *tgt = Buffer();
   Int_t readSize = blockSize;

   TString content = fXML->GetNodeContent(blocknode);

   if (blockCompressed) {
      Int_t zipSize = fXML->GetIntAttr(blocknode, xmlio::Zip);
      fUnzipBuffer = new char[zipSize];

      tgt = fUnzipBuffer;
      readSize = zipSize;
   }

   char *ptr = (char *)content.Data();

   if (gDebug > 3)
      Info("XmlReadBlock", "Content %s", ptr);

   for (int i = 0; i < readSize; i++) {
      while ((*ptr < 48) || ((*ptr > 57) && (*ptr < 97)) || (*ptr > 102))
         ptr++;

      int b_hi = (*ptr > 57) ? *ptr - 87 : *ptr - 48;
      ptr++;
      int b_lo = (*ptr > 57) ? *ptr - 87 : *ptr - 48;
      ptr++;

      *tgt = b_hi * 16 + b_lo;
      tgt++;

      if (gDebug > 4)
         Info("XmlReadBlock", "    Buf[%d] = %d", i, b_hi * 16 + b_lo);
   }

   if (fUnzipBuffer) {

      int srcsize(0), tgtsize(0), unzipRes(0);
      int status = R__unzip_header(&srcsize, (UChar_t *)fUnzipBuffer, &tgtsize);

      if (status == 0)
         R__unzip(&readSize, (unsigned char *)fUnzipBuffer, &blockSize, (unsigned char *)Buffer(), &unzipRes);

      if (status != 0 || unzipRes != blockSize)
         Error("XmlReadBlock", "Decompression error %d", unzipRes);
      else if (gDebug > 2)
         Info("XmlReadBlock", "Unzip ok");

      delete[] fUnzipBuffer;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add "ptr" attribute to node, if ptr is null or
/// if ptr is pointer on object, which is already saved in buffer
/// Automatically add "ref" attribute to node, where referenced object is stored

Bool_t TBufferXML::ProcessPointer(const void *ptr, XMLNodePointer_t node)
{
   if (!node)
      return kFALSE;

   TString refvalue;

   if (!ptr) {
      refvalue = xmlio::Null; // null
   } else {
      if (!fObjMap)
         return kFALSE;

      ULong_t hash = TString::Hash(&ptr, sizeof(void *));

      XMLNodePointer_t refnode = (XMLNodePointer_t)(Long_t)fObjMap->GetValue(hash, (Long_t)ptr);
      if (!refnode)
         return kFALSE;

      if (fXML->HasAttr(refnode, xmlio::Ref)) {
         refvalue = fXML->GetAttr(refnode, xmlio::Ref);
      } else {
         refvalue = xmlio::IdBase;
         if (XmlFile())
            refvalue += XmlFile()->GetNextRefCounter();
         else
            refvalue += GetNextRefCounter();
         fXML->NewAttr(refnode, 0, xmlio::Ref, refvalue.Data());
      }
   }
   if (refvalue.Length() > 0) {
      fXML->NewAttr(node, 0, xmlio::Ptr, refvalue.Data());
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Register pair of object pointer and node, where this object is saved,
/// in object map

void TBufferXML::RegisterPointer(const void *ptr, XMLNodePointer_t node)
{
   if (!node || !ptr)
      return;

   ULong_t hash = TString::Hash(&ptr, sizeof(void *));

   if (!fObjMap)
      fObjMap = new TExMap();

   if (fObjMap->GetValue(hash, (Long_t)ptr) == 0)
      fObjMap->Add(hash, (Long_t)ptr, (Long_t)node);
}

////////////////////////////////////////////////////////////////////////////////
/// Searches for "ptr" attribute and returns pointer to object and class,
/// if "ptr" attribute reference to read object

Bool_t TBufferXML::ExtractPointer(XMLNodePointer_t node, void *&ptr, TClass *&cl)
{
   cl = nullptr;

   if (!fXML->HasAttr(node, xmlio::Ptr))
      return kFALSE;

   const char *ptrid = fXML->GetAttr(node, xmlio::Ptr);

   if (!ptrid)
      return kFALSE;

   // null
   if (strcmp(ptrid, xmlio::Null) == 0) {
      ptr = nullptr;
      return kTRUE;
   }

   if (!fIdArray || !fObjMap)
      return kFALSE;

   TNamed *obj = static_cast<TNamed *>(fIdArray->FindObject(ptrid));
   if (obj) {
      ptr = (void *)(Long_t)fObjMap->GetValue((Long_t)fIdArray->IndexOf(obj));
      cl = TClass::GetClass(obj->GetTitle());
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Analyze if node has "ref" attribute and register it to object map

void TBufferXML::ExtractReference(XMLNodePointer_t node, const void *ptr, const TClass *cl)
{
   if (!node || !ptr)
      return;

   const char *refid = fXML->GetAttr(node, xmlio::Ref);

   if (!refid)
      return;

   if (!fIdArray) {
      fIdArray = new TObjArray;
      fIdArray->SetOwner(kTRUE);
   }
   TNamed *nid = new TNamed(refid, cl->GetName());
   fIdArray->Add(nid);

   if (!fObjMap)
      fObjMap = new TExMap();

   fObjMap->Add((Long_t)fIdArray->IndexOf(nid), (Long_t)ptr);

   if (gDebug > 2)
      Info("ExtractReference", "Find reference %s for object %p", refid, ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if node has specified name

Bool_t TBufferXML::VerifyNode(XMLNodePointer_t node, const char *name, const char *errinfo)
{
   if (!name || !node)
      return kFALSE;

   if (strcmp(fXML->GetNodeName(node), name) != 0) {
      if (errinfo) {
         Error("VerifyNode", "Reading XML file (%s). Get: %s, expects: %s", errinfo, fXML->GetNodeName(node), name);
         fErrorFlag = 1;
      }
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check, if stack node has specified name

Bool_t TBufferXML::VerifyStackNode(const char *name, const char *errinfo)
{
   return VerifyNode(StackNode(), name, errinfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks, that attribute of specified name exists and has specified value

Bool_t TBufferXML::VerifyAttr(XMLNodePointer_t node, const char *name, const char *value, const char *errinfo)
{
   if (!node || !name || !value)
      return kFALSE;

   const char *cont = fXML->GetAttr(node, name);
   if ((!cont || (strcmp(cont, value) != 0))) {
      if (errinfo) {
         Error("VerifyAttr", "%s : attr %s = %s, expected: %s", errinfo, name, cont, value);
         fErrorFlag = 1;
      }
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks stack attribute

Bool_t TBufferXML::VerifyStackAttr(const char *name, const char *value, const char *errinfo)
{
   return VerifyAttr(StackNode(), name, value, errinfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Create item node of specified name

XMLNodePointer_t TBufferXML::CreateItemNode(const char *name)
{
   XMLNodePointer_t node = 0;
   if (GetXmlLayout() == kGeneralized) {
      node = fXML->NewChild(StackNode(), 0, xmlio::Item, 0);
      fXML->NewAttr(node, 0, xmlio::Name, name);
   } else
      node = fXML->NewChild(StackNode(), 0, name, 0);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks, if stack node is item and has specified name

Bool_t TBufferXML::VerifyItemNode(const char *name, const char *errinfo)
{
   Bool_t res = kTRUE;
   if (GetXmlLayout() == kGeneralized)
      res = VerifyStackNode(xmlio::Item, errinfo) && VerifyStackAttr(xmlio::Name, name, errinfo);
   else
      res = VerifyStackNode(name, errinfo);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Create xml node correspondent to TStreamerElement object

void TBufferXML::CreateElemNode(const TStreamerElement *elem)
{
   XMLNodePointer_t elemnode = 0;

   const char *elemxmlname = XmlGetElementName(elem);

   if (GetXmlLayout() == kGeneralized) {
      elemnode = fXML->NewChild(StackNode(), 0, xmlio::Member, 0);
      fXML->NewAttr(elemnode, 0, xmlio::Name, elemxmlname);
   } else {
      // take namesapce for element only if it is not a base class or class name
      XMLNsPointer_t ns = Stack()->fClassNs;
      if ((elem->GetType() == TStreamerInfo::kBase) ||
          ((elem->GetType() == TStreamerInfo::kTNamed) && !strcmp(elem->GetName(), TNamed::Class()->GetName())) ||
          ((elem->GetType() == TStreamerInfo::kTObject) && !strcmp(elem->GetName(), TObject::Class()->GetName())) ||
          ((elem->GetType() == TStreamerInfo::kTString) && !strcmp(elem->GetName(), TString::Class()->GetName())))
         ns = 0;

      elemnode = fXML->NewChild(StackNode(), ns, elemxmlname, 0);
   }

   TXMLStackObj *curr = PushStack(elemnode);
   curr->fElem = (TStreamerElement *)elem;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if stack node correspond to TStreamerElement object

Bool_t TBufferXML::VerifyElemNode(const TStreamerElement *elem)
{
   const char *elemxmlname = XmlGetElementName(elem);

   if (GetXmlLayout() == kGeneralized) {
      if (!VerifyStackNode(xmlio::Member))
         return kFALSE;
      if (!VerifyStackAttr(xmlio::Name, elemxmlname))
         return kFALSE;
   } else {
      if (!VerifyStackNode(elemxmlname))
         return kFALSE;
   }

   PerformPreProcessing(elem, StackNode());

   TXMLStackObj *curr = PushStack(StackNode()); // set pointer to first data inside element
   curr->fElem = (TStreamerElement *)elem;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to buffer
/// If object was written before, only pointer will be stored
/// Return pointer to top xml node, representing object

XMLNodePointer_t TBufferXML::XmlWriteObject(const void *obj, const TClass *cl, Bool_t cacheReuse)
{
   XMLNodePointer_t objnode = fXML->NewChild(StackNode(), 0, xmlio::Object, 0);

   if (!cl)
      obj = 0;
   if (ProcessPointer(obj, objnode))
      return objnode;

   TString clname = XmlConvertClassName(cl->GetName());

   fXML->NewAttr(objnode, 0, xmlio::ObjClass, clname);

   if (cacheReuse)
      RegisterPointer(obj, objnode);

   PushStack(objnode);

   ((TClass *)cl)->Streamer((void *)obj, *this);

   PopStack();

   if (gDebug > 1)
      Info("XmlWriteObject", "Done write for class: %s", cl ? cl->GetName() : "null");

   return objnode;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from the buffer

void *TBufferXML::XmlReadObject(void *obj, TClass **cl)
{
   if (cl)
      *cl = nullptr;

   XMLNodePointer_t objnode = StackNode();

   if (fErrorFlag > 0)
      return obj;

   if (!objnode)
      return obj;

   if (!VerifyNode(objnode, xmlio::Object, "XmlReadObjectNew"))
      return obj;

   TClass *objClass = nullptr;

   if (ExtractPointer(objnode, obj, objClass)) {
      ShiftStack("readobjptr");
      if (cl)
         *cl = objClass;
      return obj;
   }

   TString clname = fXML->GetAttr(objnode, xmlio::ObjClass);
   objClass = XmlDefineClass(clname);
   if (objClass == TDirectory::Class())
      objClass = TDirectoryFile::Class();

   if (!objClass) {
      Error("XmlReadObject", "Cannot find class %s", clname.Data());
      ShiftStack("readobjerr");
      return obj;
   }

   if (gDebug > 1)
      Info("XmlReadObject", "Reading object of class %s", clname.Data());

   if (!obj)
      obj = objClass->New();

   ExtractReference(objnode, obj, objClass);

   PushStack(objnode);

   objClass->Streamer((void *)obj, *this);

   PopStack();

   ShiftStack("readobj");

   if (gDebug > 1)
      Info("XmlReadObject", "Reading object of class %s done", clname.Data());

   if (cl)
      *cl = objClass;

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and ReadBuffer functions
/// and indent new level in xml structure.
/// This call indicates, that TStreamerInfo functions starts streaming
/// object data of correspondent class

void TBufferXML::IncrementLevel(TVirtualStreamerInfo *info)
{
   WorkWithClass((TStreamerInfo *)info);
}

////////////////////////////////////////////////////////////////////////////////
/// Prepares buffer to stream data of specified class.

void TBufferXML::WorkWithClass(TStreamerInfo *sinfo, const TClass *cl)
{
   fCanUseCompact = kFALSE;
   fExpectedChain = kFALSE;

   if (sinfo)
      cl = sinfo->GetClass();

   if (!cl)
      return;

   TString clname = XmlConvertClassName(cl->GetName());

   if (gDebug > 2)
      Info("IncrementLevel", "Class: %s", clname.Data());

   Bool_t compressClassNode = (fExpectedBaseClass == cl);
   fExpectedBaseClass = nullptr;

   TXMLStackObj *stack = Stack();

   if (IsWriting()) {

      XMLNodePointer_t classnode = nullptr;
      if (compressClassNode) {
         classnode = StackNode();
      } else {
         if (GetXmlLayout() == kGeneralized) {
            classnode = fXML->NewChild(StackNode(), 0, xmlio::Class, 0);
            fXML->NewAttr(classnode, 0, "name", clname);
         } else
            classnode = fXML->NewChild(StackNode(), 0, clname, 0);
         stack = PushStack(classnode);
      }

      if (fVersionBuf >= -1) {
         if (fVersionBuf == -1)
            fVersionBuf = 1;
         fXML->NewIntAttr(classnode, xmlio::ClassVersion, fVersionBuf);
         fVersionBuf = -111;
      }

      if (IsUseNamespaces() && (GetXmlLayout() != kGeneralized))
         stack->fClassNs = fXML->NewNS(classnode, XmlClassNameSpaceRef(cl), clname);

   } else {
      if (!compressClassNode) {
         if (GetXmlLayout() == kGeneralized) {
            if (!VerifyStackNode(xmlio::Class, "StartInfo"))
               return;
            if (!VerifyStackAttr("name", clname, "StartInfo"))
               return;
         } else if (!VerifyStackNode(clname, "StartInfo"))
            return;
         stack = PushStack(StackNode());
      }
   }

   stack->fCompressedClassNode = compressClassNode;
   stack->fInfo = sinfo;
   stack->fIsStreamerInfo = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and ReadBuffer functions
/// and decrease level in xml structure.

void TBufferXML::DecrementLevel(TVirtualStreamerInfo *info)
{
   CheckVersionBuf();

   fCanUseCompact = kFALSE;
   fExpectedChain = kFALSE;

   if (gDebug > 2)
      Info("DecrementLevel", "Class: %s", (info ? info->GetClass()->GetName() : "custom"));

   TXMLStackObj *stack = Stack();

   if (!stack->IsStreamerInfo()) {
      PerformPostProcessing();
      stack = PopStack(); // remove stack of last element
   }

   if (stack->fCompressedClassNode) {
      stack->fInfo = 0;
      stack->fIsStreamerInfo = kFALSE;
      stack->fCompressedClassNode = kFALSE;
   } else {
      PopStack(); // back from data of stack info
      if (IsReading())
         ShiftStack("declevel"); // shift to next element after streamer info
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and ReadBuffer functions
/// and add/verify next element of xml structure
/// This calls allows separate data, correspondent to one class member, from another

void TBufferXML::SetStreamerElementNumber(TStreamerElement *elem, Int_t comptype)
{
   WorkWithElement(elem, comptype);
}

////////////////////////////////////////////////////////////////////////////////
/// This function is a part of SetStreamerElementNumber method.
/// It is introduced for reading of data for specified data member of class.
/// Used also in ReadFastArray methods to resolve problem of compressed data,
/// when several data members of the same basic type streamed with single ...FastArray call

void TBufferXML::WorkWithElement(TStreamerElement *elem, Int_t comp_type)
{
   CheckVersionBuf();

   fExpectedChain = kFALSE;
   fCanUseCompact = kFALSE;
   fExpectedBaseClass = nullptr;

   TXMLStackObj *stack = Stack();
   if (!stack) {
      Error("SetStreamerElementNumber", "stack is empty");
      return;
   }

   if (!stack->IsStreamerInfo()) { // this is not a first element
      PerformPostProcessing();
      PopStack(); // go level back
      if (IsReading())
         ShiftStack("startelem"); // shift to next element, only for reading
      stack = static_cast<TXMLStackObj *>(fStack.Last());
   }

   if (!stack) {
      Error("SetStreamerElementNumber", "Lost of stack");
      return;
   }

   if (!elem) {
      Error("SetStreamerElementNumber", "Problem in Inc/Dec level");
      return;
   }

   TStreamerInfo *info = stack->fInfo;

   if (!stack->IsStreamerInfo()) {
      Error("SetStreamerElementNumber", "Problem in Inc/Dec level");
      return;
   }
   Int_t number = info ? info->GetElements()->IndexOf(elem) : -1;

   if (gDebug > 4)
      Info("SetStreamerElementNumber", "    Next element %s", elem->GetName());

   Bool_t isBasicType = (elem->GetType() > 0) && (elem->GetType() < 20);

   fExpectedChain = isBasicType && (comp_type - elem->GetType() == TStreamerInfo::kOffsetL);

   if (fExpectedChain && (gDebug > 3)) {
      Info("SetStreamerElementNumber", "    Expects chain for elem %s number %d", elem->GetName(), number);
   }

   fCanUseCompact =
      isBasicType && ((elem->GetType() == comp_type) || (elem->GetType() == comp_type - TStreamerInfo::kConv) ||
                      (elem->GetType() == comp_type - TStreamerInfo::kSkip));

   if ((elem->GetType() == TStreamerInfo::kBase) ||
       ((elem->GetType() == TStreamerInfo::kTNamed) && !strcmp(elem->GetName(), TNamed::Class()->GetName())))
      fExpectedBaseClass = elem->GetClassPointer();

   if (fExpectedBaseClass && (gDebug > 3))
      Info("SetStreamerElementNumber", "   Expects base class %s with standard streamer",
           fExpectedBaseClass->GetName());

   if (IsWriting()) {
      CreateElemNode(elem);
   } else {
      if (!VerifyElemNode(elem))
         return;
   }

   stack = Stack();
   stack->fElemNumber = number;
   stack->fIsElemOwner = (number < 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Should be called at the beginning of custom class streamer.
///
/// Informs buffer data about class which will be streamed now.
/// ClassBegin(), ClassEnd() and ClassMemeber() should be used in
/// custom class streamers to specify which kind of data are
/// now streamed. Such information is used to correctly
/// convert class data to XML. Without that functions calls
/// classes with custom streamers cannot be used with TBufferXML

void TBufferXML::ClassBegin(const TClass *cl, Version_t)
{
   WorkWithClass(nullptr, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Should be called at the end of custom streamer
/// See TBufferXML::ClassBegin for more details

void TBufferXML::ClassEnd(const TClass *)
{
   DecrementLevel(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Method indicates name and typename of class member,
/// which should be now streamed in custom streamer
///
/// Following combinations are supported:
///   -# name = "ClassName", typeName = 0 or typename==ClassName.
/// This is a case, when data of parent class "ClassName" should be streamed.
/// For instance, if class directly inherited from TObject, custom streamer
/// should include following code:
/// ~~~{.cpp}
/// b.ClassMember("TObject");
/// TObject::Streamer(b);
/// ~~~
///   -# Basic data type
/// ~~~{.cpp}
/// b.ClassMember("fInt","Int_t");
/// b >> fInt;
/// ~~~
///   -# Array of basic data types
/// ~~~{.cpp}
/// b.ClassMember("fArr","Int_t", 5);
/// b.ReadFastArray(fArr, 5);
/// ~~~
///   -# Object as data member
/// ~~~{.cpp}
/// b.ClassMemeber("fName","TString");
/// fName.Streamer(b);
/// ~~~
///   -# Pointer on object as data member
/// ~~~{.cpp}
/// b.ClassMemeber("fObj","TObject*");
/// b.StreamObject(fObj);
/// ~~~
///
/// Arrsize1 and arrsize2 arguments (when specified) indicate first and
/// second dimension of array. Can be used for array of basic types.
/// See ClassBegin() method for more details.

void TBufferXML::ClassMember(const char *name, const char *typeName, Int_t arrsize1, Int_t arrsize2)
{
   if (!typeName)
      typeName = name;

   if (!name || (strlen(name) == 0)) {
      Error("ClassMember", "Invalid member name");
      fErrorFlag = 1;
      return;
   }

   TString tname = typeName;

   Int_t typ_id(-1), comp_type(-1);

   if (strcmp(typeName, "raw:data") == 0)
      typ_id = TStreamerInfo::kMissing;

   if (typ_id < 0) {
      TDataType *dt = gROOT->GetType(typeName);
      if (dt)
         if ((dt->GetType() > 0) && (dt->GetType() < 20))
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
         fErrorFlag = 1;
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
      comp_type = typ_id;
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
      fErrorFlag = 1;
      return;
   }

   if (arrsize1 > 0) {
      elem->SetArrayDim(arrsize2 > 0 ? 2 : 1);
      elem->SetMaxIndex(0, arrsize1);
      if (arrsize2 > 0)
         elem->SetMaxIndex(1, arrsize2);
   }

   // we indicate that there is no streamerinfo
   WorkWithElement(elem, comp_type);
}

////////////////////////////////////////////////////////////////////////////////
/// Function is converts TObject and TString structures to more compact representation

void TBufferXML::PerformPostProcessing()
{
   if (GetXmlLayout() == kGeneralized)
      return;

   const TStreamerElement *elem = Stack()->fElem;
   XMLNodePointer_t elemnode = IsWriting() ? Stack()->fNode : Stack(1)->fNode;

   if (!elem || !elemnode)
      return;

   if (elem->GetType() == TStreamerInfo::kTString) {

      XMLNodePointer_t node = fXML->GetChild(elemnode);
      fXML->SkipEmpty(node);

      XMLNodePointer_t nodecharstar(nullptr), nodeuchar(nullptr), nodeint(nullptr), nodestring(nullptr);

      while (node) {
         const char *name = fXML->GetNodeName(node);
         if (strcmp(name, xmlio::String) == 0) {
            if (nodestring)
               return;
            nodestring = node;
         } else if (strcmp(name, xmlio::UChar) == 0) {
            if (nodeuchar)
               return;
            nodeuchar = node;
         } else if (strcmp(name, xmlio::Int) == 0) {
            if (nodeint)
               return;
            nodeint = node;
         } else if (strcmp(name, xmlio::CharStar) == 0) {
            if (nodecharstar)
               return;
            nodecharstar = node;
         } else
            return; // can not be something else
         fXML->ShiftToNext(node);
      }

      TString str;

      if (GetIOVersion() < 3) {
         if (!nodeuchar)
            return;
         if (nodecharstar)
            str = fXML->GetAttr(nodecharstar, xmlio::v);
         fXML->UnlinkFreeNode(nodeuchar);
         fXML->UnlinkFreeNode(nodeint);
         fXML->UnlinkFreeNode(nodecharstar);
      } else {
         if (nodestring)
            str = fXML->GetAttr(nodestring, xmlio::v);
         fXML->UnlinkFreeNode(nodestring);
      }

      fXML->NewAttr(elemnode, 0, "str", str);
   } else if (elem->GetType() == TStreamerInfo::kTObject) {
      XMLNodePointer_t node = fXML->GetChild(elemnode);
      fXML->SkipEmpty(node);

      XMLNodePointer_t vnode = nullptr, idnode = nullptr, bitsnode = nullptr, prnode = nullptr;

      while (node) {
         const char *name = fXML->GetNodeName(node);

         if (strcmp(name, xmlio::OnlyVersion) == 0) {
            if (vnode)
               return;
            vnode = node;
         } else if (strcmp(name, xmlio::UInt) == 0) {
            if (!idnode)
               idnode = node;
            else if (!bitsnode)
               bitsnode = node;
            else
               return;
         } else if (strcmp(name, xmlio::UShort) == 0) {
            if (prnode)
               return;
            prnode = node;
         } else
            return;
         fXML->ShiftToNext(node);
      }

      if (!vnode || !idnode || !bitsnode)
         return;

      TString str = fXML->GetAttr(idnode, xmlio::v);
      fXML->NewAttr(elemnode, 0, "fUniqueID", str);

      str = fXML->GetAttr(bitsnode, xmlio::v);
      UInt_t bits;
      sscanf(str.Data(), "%u", &bits);

      char sbuf[20];
      snprintf(sbuf, sizeof(sbuf), "%x", bits);
      fXML->NewAttr(elemnode, 0, "fBits", sbuf);

      if (prnode) {
         str = fXML->GetAttr(prnode, xmlio::v);
         fXML->NewAttr(elemnode, 0, "fProcessID", str);
      }

      fXML->UnlinkFreeNode(vnode);
      fXML->UnlinkFreeNode(idnode);
      fXML->UnlinkFreeNode(bitsnode);
      fXML->UnlinkFreeNode(prnode);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Function is unpack TObject and TString structures to be able read
/// them from custom streamers of this objects

void TBufferXML::PerformPreProcessing(const TStreamerElement *elem, XMLNodePointer_t elemnode)
{
   if (GetXmlLayout() == kGeneralized)
      return;
   if (!elem || !elemnode)
      return;

   if (elem->GetType() == TStreamerInfo::kTString) {

      if (!fXML->HasAttr(elemnode, "str"))
         return;
      TString str = fXML->GetAttr(elemnode, "str");
      fXML->FreeAttr(elemnode, "str");

      if (GetIOVersion() < 3) {
         Int_t len = str.Length();
         XMLNodePointer_t ucharnode = fXML->NewChild(elemnode, 0, xmlio::UChar, 0);
         char sbuf[20];
         snprintf(sbuf, sizeof(sbuf), "%d", len);
         if (len < 255) {
            fXML->NewAttr(ucharnode, 0, xmlio::v, sbuf);
         } else {
            fXML->NewAttr(ucharnode, 0, xmlio::v, "255");
            XMLNodePointer_t intnode = fXML->NewChild(elemnode, 0, xmlio::Int, 0);
            fXML->NewAttr(intnode, 0, xmlio::v, sbuf);
         }
         if (len > 0) {
            XMLNodePointer_t node = fXML->NewChild(elemnode, 0, xmlio::CharStar, 0);
            fXML->NewAttr(node, 0, xmlio::v, str);
         }
      } else {
         XMLNodePointer_t node = fXML->NewChild(elemnode, 0, xmlio::String, 0);
         fXML->NewAttr(node, 0, xmlio::v, str);
      }
   } else if (elem->GetType() == TStreamerInfo::kTObject) {
      if (!fXML->HasAttr(elemnode, "fUniqueID"))
         return;
      if (!fXML->HasAttr(elemnode, "fBits"))
         return;

      TString idstr = fXML->GetAttr(elemnode, "fUniqueID");
      TString bitsstr = fXML->GetAttr(elemnode, "fBits");
      TString prstr = fXML->GetAttr(elemnode, "fProcessID");

      fXML->FreeAttr(elemnode, "fUniqueID");
      fXML->FreeAttr(elemnode, "fBits");
      fXML->FreeAttr(elemnode, "fProcessID");

      XMLNodePointer_t node = fXML->NewChild(elemnode, 0, xmlio::OnlyVersion, 0);
      fXML->NewAttr(node, 0, xmlio::v, "1");

      node = fXML->NewChild(elemnode, 0, xmlio::UInt, 0);
      fXML->NewAttr(node, 0, xmlio::v, idstr);

      UInt_t bits;
      sscanf(bitsstr.Data(), "%x", &bits);
      char sbuf[20];
      snprintf(sbuf, sizeof(sbuf), "%u", bits);

      node = fXML->NewChild(elemnode, 0, xmlio::UInt, 0);
      fXML->NewAttr(node, 0, xmlio::v, sbuf);

      if (prstr.Length() > 0) {
         node = fXML->NewChild(elemnode, 0, xmlio::UShort, 0);
         fXML->NewAttr(node, 0, xmlio::v, prstr.Data());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called before any IO operation of TBuffer
/// Now is used to store version value if no proper calls are discovered

void TBufferXML::BeforeIOoperation()
{
   CheckVersionBuf();
}

////////////////////////////////////////////////////////////////////////////////
/// Function to read class from buffer, used in old-style streamers

TClass *TBufferXML::ReadClass(const TClass *, UInt_t *)
{
   const char *clname = nullptr;

   if (VerifyItemNode(xmlio::Class)) {
      clname = XmlReadValue(xmlio::Class);
   }

   if (gDebug > 2)
      Info("ReadClass", "Try to read class %s", clname ? clname : "---");

   return clname ? gROOT->GetClass(clname) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Function to write class into buffer, used in old-style streamers

void TBufferXML::WriteClass(const TClass *cl)
{
   if (gDebug > 2)
      Info("WriteClass", "Try to write class %s", cl->GetName());

   XmlWriteValue(cl->GetName(), xmlio::Class);
}

////////////////////////////////////////////////////////////////////////////////
/// Suppressed function of TBuffer

Int_t TBufferXML::CheckByteCount(UInt_t /*r_s */, UInt_t /*r_c*/, const TClass * /*cl*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Suppressed function of TBuffer

Int_t TBufferXML::CheckByteCount(UInt_t, UInt_t, const char *)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Suppressed function of TBuffer

void TBufferXML::SetByteCount(UInt_t, Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Skip class version from I/O buffer.

void TBufferXML::SkipVersion(const TClass *cl)
{
   ReadVersion(0, 0, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Read version value from buffer

Version_t TBufferXML::ReadVersion(UInt_t *start, UInt_t *bcnt, const TClass * /*cl*/)
{
   BeforeIOoperation();

   Version_t res = 0;

   if (start)
      *start = 0;
   if (bcnt)
      *bcnt = 0;

   if (VerifyItemNode(xmlio::OnlyVersion)) {
      res = AtoI(XmlReadValue(xmlio::OnlyVersion));
   } else if (fExpectedBaseClass && (fXML->HasAttr(Stack(1)->fNode, xmlio::ClassVersion))) {
      res = fXML->GetIntAttr(Stack(1)->fNode, xmlio::ClassVersion);
   } else if (fXML->HasAttr(StackNode(), xmlio::ClassVersion)) {
      res = fXML->GetIntAttr(StackNode(), xmlio::ClassVersion);
   } else {
      Error("ReadVersion", "No correspondent tags to read version");
      fErrorFlag = 1;
   }

   if (gDebug > 2)
      Info("ReadVersion", "Version = %d", res);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Read class version from I/O buffer, when the caller knows for sure that
/// there is no checksum written/involved.

Version_t TBufferXML::ReadVersionNoCheckSum(UInt_t *, UInt_t *)
{
   Info("ReadVersionNoCheckSum", "No idea if required");

   return ReadVersion();
}

////////////////////////////////////////////////////////////////////////////////
/// Checks buffer, filled by WriteVersion
/// if next data is arriving, version should be stored in buffer

void TBufferXML::CheckVersionBuf()
{
   if (IsWriting() && (fVersionBuf >= -100)) {
      char sbuf[20];
      snprintf(sbuf, sizeof(sbuf), "%d", fVersionBuf);
      XmlWriteValue(sbuf, xmlio::OnlyVersion);
      fVersionBuf = -111;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copies class version to buffer, but not writes it to xml
/// Version will be written with next I/O operation or
/// will be added as attribute of class tag, created by IncrementLevel call

UInt_t TBufferXML::WriteVersion(const TClass *cl, Bool_t /* useBcnt */)
{
   BeforeIOoperation();

   if (fExpectedBaseClass != cl)
      fExpectedBaseClass = 0;

   fVersionBuf = cl->GetClassVersion();

   if (gDebug > 2)
      Info("WriteVersion", "Class: %s, version = %d", cl->GetName(), fVersionBuf);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from buffer. Only used from TBuffer

void *TBufferXML::ReadObjectAny(const TClass *)
{
   BeforeIOoperation();
   if (gDebug > 2)
      Info("ReadObjectAny", "From node %s", fXML->GetNodeName(StackNode()));
   void *res = XmlReadObject(nullptr);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Skip any kind of object from buffer
/// Actually skip only one node on current level of xml structure

void TBufferXML::SkipObjectAny()
{
   ShiftStack("skipobjectany");
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to buffer. Only used from TBuffer

void TBufferXML::WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse)
{
   BeforeIOoperation();
   if (gDebug > 2)
      Info("WriteObject", "Class %s", (actualClass ? actualClass->GetName() : " null"));
   XmlWriteObject(actualObjStart, actualClass, cacheReuse);
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to read array content

template <typename T>
R__ALWAYS_INLINE void TBufferXML::XmlReadArrayContent(T *arr, Int_t arrsize)
{
   Int_t indx = 0, cnt, curr;
   while (indx < arrsize) {
      cnt = 1;
      if (fXML->HasAttr(StackNode(), xmlio::cnt))
         cnt = fXML->GetIntAttr(StackNode(), xmlio::cnt);
      XmlReadBasic(arr[indx]);
      curr = indx++;
      while (cnt-- > 1)
         arr[indx++] = arr[curr];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to read array with size attribute
/// If necessary, array is created

template <typename T>
R__ALWAYS_INLINE Int_t TBufferXML::XmlReadArray(T *&arr, bool is_static)
{
   BeforeIOoperation();
   if (!VerifyItemNode(xmlio::Array, is_static ? "ReadStaticArray" : "ReadArray"))
      return 0;
   Int_t n = fXML->GetIntAttr(StackNode(), xmlio::Size);
   if (n <= 0)
      return 0;
   if (!arr) {
      if (is_static)
         return 0;
      arr = new T[n];
   }
   PushStack(StackNode());
   XmlReadArrayContent(arr, n);
   PopStack();
   ShiftStack(is_static ? "readstatarr" : "readarr");
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Float16_t from the buffer

void TBufferXML::ReadFloat16(Float_t *f, TStreamerElement * /*ele*/)
{
   BeforeIOoperation();
   XmlReadBasic(*f);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer

void TBufferXML::ReadDouble32(Double_t *d, TStreamerElement * /*ele*/)
{
   BeforeIOoperation();
   XmlReadBasic(*d);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer when the factor and minimun value have been specified
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().
/// Currently TBufferXML does not optimize space in this case.

void TBufferXML::ReadWithFactor(Float_t *ptr, Double_t /* factor */, Double_t /* minvalue */)
{
   BeforeIOoperation();
   XmlReadBasic(*ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Float16_t from the buffer when the number of bits is specified (explicitly or not)
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16().
/// Currently TBufferXML does not optimize space in this case.

void TBufferXML::ReadWithNbits(Float_t *ptr, Int_t /* nbits */)
{
   BeforeIOoperation();
   XmlReadBasic(*ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer when the factor and minimum value have been specified
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().
/// Currently TBufferXML does not optimize space in this case.

void TBufferXML::ReadWithFactor(Double_t *ptr, Double_t /* factor */, Double_t /* minvalue */)
{
   BeforeIOoperation();
   XmlReadBasic(*ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer when the number of bits is specified (explicitly or not)
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().
/// Currently TBufferXML does not optimize space in this case.

void TBufferXML::ReadWithNbits(Double_t *ptr, Int_t /* nbits */)
{
   BeforeIOoperation();
   XmlReadBasic(*ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Write a Float16_t to the buffer

void TBufferXML::WriteFloat16(Float_t *f, TStreamerElement * /*ele*/)
{
   BeforeIOoperation();
   XmlWriteBasic(*f);
}

////////////////////////////////////////////////////////////////////////////////
/// Write a Double32_t to the buffer

void TBufferXML::WriteDouble32(Double_t *d, TStreamerElement * /*ele*/)
{
   BeforeIOoperation();
   XmlWriteBasic(*d);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

Int_t TBufferXML::ReadArray(Bool_t *&b)
{
   return XmlReadArray(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer

Int_t TBufferXML::ReadArray(Char_t *&c)
{
   return XmlReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

Int_t TBufferXML::ReadArray(UChar_t *&c)
{
   return XmlReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

Int_t TBufferXML::ReadArray(Short_t *&h)
{
   return XmlReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

Int_t TBufferXML::ReadArray(UShort_t *&h)
{
   return XmlReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

Int_t TBufferXML::ReadArray(Int_t *&i)
{
   return XmlReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

Int_t TBufferXML::ReadArray(UInt_t *&i)
{
   return XmlReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

Int_t TBufferXML::ReadArray(Long_t *&l)
{
   return XmlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

Int_t TBufferXML::ReadArray(ULong_t *&l)
{
   return XmlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

Int_t TBufferXML::ReadArray(Long64_t *&l)
{
   return XmlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

Int_t TBufferXML::ReadArray(ULong64_t *&l)
{
   return XmlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

Int_t TBufferXML::ReadArray(Float_t *&f)
{
   return XmlReadArray(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

Int_t TBufferXML::ReadArray(Double_t *&d)
{
   return XmlReadArray(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

Int_t TBufferXML::ReadArrayFloat16(Float_t *&f, TStreamerElement * /*ele*/)
{
   return XmlReadArray(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

Int_t TBufferXML::ReadArrayDouble32(Double_t *&d, TStreamerElement * /*ele*/)
{
   return XmlReadArray(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

Int_t TBufferXML::ReadStaticArray(Bool_t *b)
{
   return XmlReadArray(b, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer

Int_t TBufferXML::ReadStaticArray(Char_t *c)
{
   return XmlReadArray(c, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

Int_t TBufferXML::ReadStaticArray(UChar_t *c)
{
   return XmlReadArray(c, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

Int_t TBufferXML::ReadStaticArray(Short_t *h)
{
   return XmlReadArray(h, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

Int_t TBufferXML::ReadStaticArray(UShort_t *h)
{
   return XmlReadArray(h, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

Int_t TBufferXML::ReadStaticArray(Int_t *i)
{
   return XmlReadArray(i, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

Int_t TBufferXML::ReadStaticArray(UInt_t *i)
{
   return XmlReadArray(i, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

Int_t TBufferXML::ReadStaticArray(Long_t *l)
{
   return XmlReadArray(l, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

Int_t TBufferXML::ReadStaticArray(ULong_t *l)
{
   return XmlReadArray(l, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

Int_t TBufferXML::ReadStaticArray(Long64_t *l)
{
   return XmlReadArray(l, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

Int_t TBufferXML::ReadStaticArray(ULong64_t *l)
{
   return XmlReadArray(l, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

Int_t TBufferXML::ReadStaticArray(Float_t *f)
{
   return XmlReadArray(f, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

Int_t TBufferXML::ReadStaticArray(Double_t *d)
{
   return XmlReadArray(d, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

Int_t TBufferXML::ReadStaticArrayFloat16(Float_t *f, TStreamerElement * /*ele*/)
{
   return XmlReadArray(f, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

Int_t TBufferXML::ReadStaticArrayDouble32(Double_t *d, TStreamerElement * /*ele*/)
{
   return XmlReadArray(d, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to read content of array, which not include size of array
/// Also treated situation, when instead of one single array chain
/// of several elements should be produced

template <typename T>
R__ALWAYS_INLINE void TBufferXML::XmlReadFastArray(T *arr, Int_t n)
{
   BeforeIOoperation();
   if (n <= 0)
      return;
   TStreamerElement *elem = Stack(0)->fElem;
   if (elem && (elem->GetType() > TStreamerInfo::kOffsetL) && (elem->GetType() < TStreamerInfo::kOffsetP) &&
       (elem->GetArrayLength() != n))
      fExpectedChain = kTRUE;
   if (fExpectedChain) {
      fExpectedChain = kFALSE;
      Int_t startnumber = Stack(0)->fElemNumber;
      TStreamerInfo *info = Stack(1)->fInfo;
      Int_t index = 0;
      while (index < n) {
         elem = (TStreamerElement *)info->GetElements()->At(startnumber++);
         if (elem->GetType() < TStreamerInfo::kOffsetL) {
            if (index > 0) {
               PopStack();
               ShiftStack("chainreader");
               VerifyElemNode(elem);
            }
            fCanUseCompact = kTRUE;
            XmlReadBasic(arr[index]);
            index++;
         } else {
            if (!VerifyItemNode(xmlio::Array, "ReadFastArray"))
               return;
            PushStack(StackNode());
            Int_t elemlen = elem->GetArrayLength();
            XmlReadArrayContent((arr + index), elemlen);
            PopStack();
            ShiftStack("readfastarr");
            index += elemlen;
         }
      }
   } else {
      if (!VerifyItemNode(xmlio::Array, "ReadFastArray"))
         return;
      PushStack(StackNode());
      XmlReadArrayContent(arr, n);
      PopStack();
      ShiftStack("readfastarr");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

void TBufferXML::ReadFastArray(Bool_t *b, Int_t n)
{
   XmlReadFastArray(b, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer
/// if nodename==CharStar, read all array as string

void TBufferXML::ReadFastArray(Char_t *c, Int_t n)
{
   if ((n > 0) && VerifyItemNode(xmlio::CharStar)) {
      const char *buf;
      if ((buf = XmlReadValue(xmlio::CharStar))) {
         Int_t size = strlen(buf);
         if (size < n)
            size = n;
         memcpy(c, buf, size);
      }
   } else {
      XmlReadFastArray(c, n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n characters from the I/O buffer.
/// Used only from TLeafC, dummy implementation here

void TBufferXML::ReadFastArrayString(Char_t *c, Int_t n)
{
   ReadFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

void TBufferXML::ReadFastArray(UChar_t *c, Int_t n)
{
   XmlReadFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

void TBufferXML::ReadFastArray(Short_t *h, Int_t n)
{
   XmlReadFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

void TBufferXML::ReadFastArray(UShort_t *h, Int_t n)
{
   XmlReadFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

void TBufferXML::ReadFastArray(Int_t *i, Int_t n)
{
   XmlReadFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

void TBufferXML::ReadFastArray(UInt_t *i, Int_t n)
{
   XmlReadFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

void TBufferXML::ReadFastArray(Long_t *l, Int_t n)
{
   XmlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

void TBufferXML::ReadFastArray(ULong_t *l, Int_t n)
{
   XmlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

void TBufferXML::ReadFastArray(Long64_t *l, Int_t n)
{
   XmlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

void TBufferXML::ReadFastArray(ULong64_t *l, Int_t n)
{
   XmlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

void TBufferXML::ReadFastArray(Float_t *f, Int_t n)
{
   XmlReadFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

void TBufferXML::ReadFastArray(Double_t *d, Int_t n)
{
   XmlReadFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

void TBufferXML::ReadFastArrayFloat16(Float_t *f, Int_t n, TStreamerElement * /*ele*/)
{
   XmlReadFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

void TBufferXML::ReadFastArrayWithFactor(Float_t *f, Int_t n, Double_t /* factor */, Double_t /* minvalue */)
{
   XmlReadFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

void TBufferXML::ReadFastArrayWithNbits(Float_t *f, Int_t n, Int_t /*nbits*/)
{
   XmlReadFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

void TBufferXML::ReadFastArrayDouble32(Double_t *d, Int_t n, TStreamerElement * /*ele*/)
{
   XmlReadFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

void TBufferXML::ReadFastArrayWithFactor(Double_t *d, Int_t n, Double_t /* factor */, Double_t /* minvalue */)
{
   XmlReadFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

void TBufferXML::ReadFastArrayWithNbits(Double_t *d, Int_t n, Int_t /*nbits*/)
{
   XmlReadFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read an array of 'n' objects from the I/O buffer.
/// Stores the objects read starting at the address 'start'.
/// The objects in the array are assume to be of class 'cl'.

void TBufferXML::ReadFastArray(void *start, const TClass *cl, Int_t n, TMemberStreamer *streamer,
                               const TClass *onFileClass)
{
   if (streamer) {
      streamer->SetOnFileClass(onFileClass);
      (*streamer)(*this, start, 0);
      return;
   }

   int objectSize = cl->Size();
   char *obj = (char *)start;
   char *end = obj + n * objectSize;

   for (; obj < end; obj += objectSize)
      ((TClass *)cl)->Streamer(obj, *this, onFileClass);
}

////////////////////////////////////////////////////////////////////////////////
/// Read an array of 'n' objects from the I/O buffer.
///
/// The objects read are stored starting at the address '*start'
/// The objects in the array are assumed to be of class 'cl' or a derived class.
/// 'mode' indicates whether the data member is marked with '->'

void TBufferXML::ReadFastArray(void **start, const TClass *cl, Int_t n, Bool_t isPreAlloc, TMemberStreamer *streamer,
                               const TClass *onFileClass)
{
   if (streamer) {
      if (isPreAlloc) {
         for (Int_t j = 0; j < n; j++) {
            if (!start[j])
               start[j] = cl->New();
         }
      }
      streamer->SetOnFileClass(onFileClass);
      (*streamer)(*this, (void *)start, 0);
      return;
   }

   if (!isPreAlloc) {

      for (Int_t j = 0; j < n; j++) {
         // delete the object or collection
         void *old = start[j];
         start[j] = ReadObjectAny(cl);
         if (old && old != start[j] && TStreamerInfo::CanDelete()
             // There are some cases where the user may set up a pointer in the (default)
             // constructor but not mark this pointer as transient.  Sometime the value
             // of this pointer is the address of one of the object with just created
             // and the following delete would result in the deletion (possibly of the
             // top level object we are goint to return!).
             // Eventhough this is a user error, we could prevent the crash by simply
             // adding:
             // && !CheckObject(start[j],cl)
             // However this can increase the read time significantly (10% in the case
             // of one TLine pointer in the test/Track and run ./Event 200 0 0 20 30000
             //
             // If ReadObjectAny returned the same value as we previous had, this means
             // that when writing this object (start[j] had already been written and
             // is indeed pointing to the same object as the object the user set up
             // in the default constructor).
             ) {
            ((TClass *)cl)->Destructor(old, kFALSE); // call delete and desctructor
         }
      }

   } else {
      // case //-> in comment

      for (Int_t j = 0; j < n; j++) {
         if (!start[j])
            start[j] = ((TClass *)cl)->New();
         ((TClass *)cl)->Streamer(start[j], *this, onFileClass);
      }
   }
}

template <typename T>
R__ALWAYS_INLINE void TBufferXML::XmlWriteArrayContent(const T *arr, Int_t arrsize)
{
   if (fCompressLevel > 0) {
      Int_t indx = 0;
      while (indx < arrsize) {
         XMLNodePointer_t elemnode = XmlWriteBasic(arr[indx]);
         Int_t curr = indx++;
         while ((indx < arrsize) && (arr[indx] == arr[curr]))
            indx++;
         if (indx - curr > 1)
            fXML->NewIntAttr(elemnode, xmlio::cnt, indx - curr);
      }
   } else {
      for (Int_t indx = 0; indx < arrsize; indx++)
         XmlWriteBasic(arr[indx]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array, including it size
/// Content may be compressed

template <typename T>
R__ALWAYS_INLINE void TBufferXML::XmlWriteArray(const T *arr, Int_t arrsize)
{
   BeforeIOoperation();
   XMLNodePointer_t arrnode = CreateItemNode(xmlio::Array);
   fXML->NewIntAttr(arrnode, xmlio::Size, arrsize);
   PushStack(arrnode);
   XmlWriteArrayContent(arr, arrsize);
   PopStack();
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Bool_t to buffer

void TBufferXML::WriteArray(const Bool_t *b, Int_t n)
{
   XmlWriteArray(b, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer

void TBufferXML::WriteArray(const Char_t *c, Int_t n)
{
   XmlWriteArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UChar_t to buffer

void TBufferXML::WriteArray(const UChar_t *c, Int_t n)
{
   XmlWriteArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Short_t to buffer

void TBufferXML::WriteArray(const Short_t *h, Int_t n)
{
   XmlWriteArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UShort_t to buffer

void TBufferXML::WriteArray(const UShort_t *h, Int_t n)
{
   XmlWriteArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Int_ to buffer

void TBufferXML::WriteArray(const Int_t *i, Int_t n)
{
   XmlWriteArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UInt_t to buffer

void TBufferXML::WriteArray(const UInt_t *i, Int_t n)
{
   XmlWriteArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long_t to buffer

void TBufferXML::WriteArray(const Long_t *l, Int_t n)
{
   XmlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong_t to buffer

void TBufferXML::WriteArray(const ULong_t *l, Int_t n)
{
   XmlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long64_t to buffer

void TBufferXML::WriteArray(const Long64_t *l, Int_t n)
{
   XmlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong64_t to buffer

void TBufferXML::WriteArray(const ULong64_t *l, Int_t n)
{
   XmlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float_t to buffer

void TBufferXML::WriteArray(const Float_t *f, Int_t n)
{
   XmlWriteArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double_t to buffer

void TBufferXML::WriteArray(const Double_t *d, Int_t n)
{
   XmlWriteArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float16_t to buffer

void TBufferXML::WriteArrayFloat16(const Float_t *f, Int_t n, TStreamerElement * /*ele*/)
{
   XmlWriteArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double32_t to buffer

void TBufferXML::WriteArrayDouble32(const Double_t *d, Int_t n, TStreamerElement * /*ele*/)
{
   XmlWriteArray(d, n);
}

/////////////////////////////////////////////////////////////////////////////////
/// Write array without size attribute
/// Also treat situation, when instead of one single array
/// chain of several elements should be produced

template <typename T>
R__ALWAYS_INLINE void TBufferXML::XmlWriteFastArray(const T *arr, Int_t n)
{
   BeforeIOoperation();
   if (n <= 0)
      return;
   TStreamerElement *elem = Stack(0)->fElem;
   if (elem && (elem->GetType() > TStreamerInfo::kOffsetL) && (elem->GetType() < TStreamerInfo::kOffsetP) &&
       (elem->GetArrayLength() != n))
      fExpectedChain = kTRUE;
   if (fExpectedChain) {
      TStreamerInfo *info = Stack(1)->fInfo;
      Int_t startnumber = Stack(0)->fElemNumber;
      fExpectedChain = kFALSE;
      Int_t index = 0;
      while (index < n) {
         elem = (TStreamerElement *)info->GetElements()->At(startnumber++);
         if (elem->GetType() < TStreamerInfo::kOffsetL) {
            if (index > 0) {
               PopStack();
               CreateElemNode(elem);
            }
            fCanUseCompact = kTRUE;
            XmlWriteBasic(arr[index]);
            index++;
         } else {
            XMLNodePointer_t arrnode = CreateItemNode(xmlio::Array);
            Int_t elemlen = elem->GetArrayLength();
            PushStack(arrnode);
            XmlWriteArrayContent((arr + index), elemlen);
            index += elemlen;
            PopStack();
         }
      }
   } else {
      XMLNodePointer_t arrnode = CreateItemNode(xmlio::Array);
      PushStack(arrnode);
      XmlWriteArrayContent(arr, n);
      PopStack();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Bool_t to buffer

void TBufferXML::WriteFastArray(const Bool_t *b, Int_t n)
{
   XmlWriteFastArray(b, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer
/// If array does not include any special characters,
/// it will be reproduced as CharStar node with string as attribute

void TBufferXML::WriteFastArray(const Char_t *c, Int_t n)
{
   Bool_t usedefault = (n == 0) || fExpectedChain;
   const Char_t *buf = c;
   if (!usedefault)
      for (int i = 0; i < n; i++) {
         if (*buf < 27) {
            usedefault = kTRUE;
            break;
         }
         buf++;
      }
   if (usedefault) {
      XmlWriteFastArray(c, n);
   } else {
      Char_t *buf2 = new Char_t[n + 1];
      memcpy(buf2, c, n);
      buf2[n] = 0;
      XmlWriteValue(buf2, xmlio::CharStar);
      delete[] buf2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UChar_t to buffer

void TBufferXML::WriteFastArray(const UChar_t *c, Int_t n)
{
   XmlWriteFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Short_t to buffer

void TBufferXML::WriteFastArray(const Short_t *h, Int_t n)
{
   XmlWriteFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UShort_t to buffer

void TBufferXML::WriteFastArray(const UShort_t *h, Int_t n)
{
   XmlWriteFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Int_t to buffer

void TBufferXML::WriteFastArray(const Int_t *i, Int_t n)
{
   XmlWriteFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UInt_t to buffer

void TBufferXML::WriteFastArray(const UInt_t *i, Int_t n)
{
   XmlWriteFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long_t to buffer

void TBufferXML::WriteFastArray(const Long_t *l, Int_t n)
{
   XmlWriteFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong_t to buffer

void TBufferXML::WriteFastArray(const ULong_t *l, Int_t n)
{
   XmlWriteFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long64_t to buffer

void TBufferXML::WriteFastArray(const Long64_t *l, Int_t n)
{
   XmlWriteFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong64_t to buffer

void TBufferXML::WriteFastArray(const ULong64_t *l, Int_t n)
{
   XmlWriteFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float_t to buffer

void TBufferXML::WriteFastArray(const Float_t *f, Int_t n)
{
   XmlWriteFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double_t to buffer

void TBufferXML::WriteFastArray(const Double_t *d, Int_t n)
{
   XmlWriteFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float16_t to buffer

void TBufferXML::WriteFastArrayFloat16(const Float_t *f, Int_t n, TStreamerElement * /*ele*/)
{
   XmlWriteFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double32_t to buffer

void TBufferXML::WriteFastArrayDouble32(const Double_t *d, Int_t n, TStreamerElement * /*ele*/)
{
   XmlWriteFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n characters into the I/O buffer.
/// Used only by TLeafC, just dummy implementation here

void TBufferXML::WriteFastArrayString(const Char_t *c, Int_t n)
{
   WriteFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write an array of object starting at the address 'start' and of length 'n'
/// the objects in the array are assumed to be of class 'cl'

void TBufferXML::WriteFastArray(void *start, const TClass *cl, Int_t n, TMemberStreamer *streamer)
{
   if (streamer) {
      (*streamer)(*this, start, 0);
      return;
   }

   char *obj = (char *)start;
   if (!n)
      n = 1;
   int size = cl->Size();

   for (Int_t j = 0; j < n; j++, obj += size) {
      ((TClass *)cl)->Streamer(obj, *this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write an array of object starting at the address '*start' and of length 'n'
/// the objects in the array are of class 'cl'
/// 'isPreAlloc' indicates whether the data member is marked with '->'
/// Return:
///   - 0: success
///   - 2: truncated success (i.e actual class is missing. Only ptrClass saved.)

Int_t TBufferXML::WriteFastArray(void **start, const TClass *cl, Int_t n, Bool_t isPreAlloc, TMemberStreamer *streamer)
{
   // if isPreAlloc is true (data member has a ->) we can assume that the pointer
   // is never 0.

   if (streamer) {
      (*streamer)(*this, (void *)start, 0);
      return 0;
   }

   int strInfo = 0;

   Int_t res = 0;

   if (!isPreAlloc) {

      for (Int_t j = 0; j < n; j++) {
         // must write StreamerInfo if pointer is null
         if (!strInfo && !start[j]) {
            if (cl->Property() & kIsAbstract) {
               // Do not try to generate the StreamerInfo for an abstract class
            } else {
               TStreamerInfo *info = (TStreamerInfo *)((TClass *)cl)->GetStreamerInfo();
               ForceWriteInfo(info, kFALSE);
            }
         }
         strInfo = 2003;
         res |= WriteObjectAny(start[j], cl);
      }

   } else {
      // case //-> in comment

      for (Int_t j = 0; j < n; j++) {
         if (!start[j])
            start[j] = ((TClass *)cl)->New();
         ((TClass *)cl)->Streamer(start[j], *this);
      }
   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// steram object to/from buffer

void TBufferXML::StreamObject(void *obj, const std::type_info &typeinfo, const TClass * /* onFileClass */)
{
   StreamObject(obj, TClass::GetClass(typeinfo));
}

////////////////////////////////////////////////////////////////////////////////
/// steram object to/from buffer

void TBufferXML::StreamObject(void *obj, const char *className, const TClass * /* onFileClass */)
{
   StreamObject(obj, TClass::GetClass(className));
}

void TBufferXML::StreamObject(TObject *obj)
{
   // steram object to/from buffer

   StreamObject(obj, obj ? obj->IsA() : TObject::Class());
}

////////////////////////////////////////////////////////////////////////////////
/// Stream object to/from buffer

void TBufferXML::StreamObject(void *obj, const TClass *cl, const TClass * /* onfileClass */)
{
   BeforeIOoperation();
   if (gDebug > 1)
      Info("StreamObject", "Class: %s", (cl ? cl->GetName() : "none"));
   if (IsReading())
      XmlReadObject(obj);
   else
      XmlWriteObject(obj, cl, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Bool_t value from buffer

void TBufferXML::ReadBool(Bool_t &b)
{
   BeforeIOoperation();
   XmlReadBasic(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Char_t value from buffer

void TBufferXML::ReadChar(Char_t &c)
{
   BeforeIOoperation();
   XmlReadBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UChar_t value from buffer

void TBufferXML::ReadUChar(UChar_t &c)
{
   BeforeIOoperation();
   XmlReadBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Short_t value from buffer

void TBufferXML::ReadShort(Short_t &h)
{
   BeforeIOoperation();
   XmlReadBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UShort_t value from buffer

void TBufferXML::ReadUShort(UShort_t &h)
{
   BeforeIOoperation();
   XmlReadBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Int_t value from buffer

void TBufferXML::ReadInt(Int_t &i)
{
   BeforeIOoperation();
   XmlReadBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UInt_t value from buffer

void TBufferXML::ReadUInt(UInt_t &i)
{
   BeforeIOoperation();
   XmlReadBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Long_t value from buffer

void TBufferXML::ReadLong(Long_t &l)
{
   BeforeIOoperation();
   XmlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads ULong_t value from buffer

void TBufferXML::ReadULong(ULong_t &l)
{
   BeforeIOoperation();
   XmlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Long64_t value from buffer

void TBufferXML::ReadLong64(Long64_t &l)
{
   BeforeIOoperation();
   XmlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads ULong64_t value from buffer

void TBufferXML::ReadULong64(ULong64_t &l)
{
   BeforeIOoperation();
   XmlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Float_t value from buffer

void TBufferXML::ReadFloat(Float_t &f)
{
   BeforeIOoperation();
   XmlReadBasic(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Double_t value from buffer

void TBufferXML::ReadDouble(Double_t &d)
{
   BeforeIOoperation();
   XmlReadBasic(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads array of characters from buffer

void TBufferXML::ReadCharP(Char_t *c)
{
   BeforeIOoperation();
   const char *buf;
   if ((buf = XmlReadValue(xmlio::CharStar)))
      strcpy(c, buf);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads a TString

void TBufferXML::ReadTString(TString &s)
{
   if (GetIOVersion() < 3) {
      // original TBufferFile method can not be used, while used TString methods are private
      // try to reimplement close to the original
      Int_t nbig;
      UChar_t nwh;
      *this >> nwh;
      if (nwh == 0) {
         s.Resize(0);
      } else {
         if (nwh == 255)
            *this >> nbig;
         else
            nbig = nwh;

         char *data = new char[nbig];
         data[nbig] = 0;
         ReadFastArray(data, nbig);
         s = data;
         delete[] data;
      }
   } else {
      BeforeIOoperation();
      const char *buf = XmlReadValue(xmlio::String);
      if (buf)
         s = buf;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reads a std::string

void TBufferXML::ReadStdString(std::string *obj)
{
   if (GetIOVersion() < 3) {
      if (!obj) {
         Error("ReadStdString", "The std::string address is nullptr but should not");
         return;
      }
      Int_t nbig;
      UChar_t nwh;
      *this >> nwh;
      if (nwh == 0) {
         obj->clear();
      } else {
         if (obj->size()) {
            // Insure that the underlying data storage is not shared
            (*obj)[0] = '\0';
         }
         if (nwh == 255) {
            *this >> nbig;
            obj->resize(nbig, '\0');
            ReadFastArray((char *)obj->data(), nbig);
         } else {
            obj->resize(nwh, '\0');
            ReadFastArray((char *)obj->data(), nwh);
         }
      }
   } else {
      BeforeIOoperation();
      const char *buf = XmlReadValue(xmlio::String);
      if (buf && obj)
         *obj = buf;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a char* string

void TBufferXML::ReadCharStar(char *&s)
{
   delete[] s;
   s = nullptr;

   Int_t nch;
   *this >> nch;
   if (nch > 0) {
      s = new char[nch + 1];
      ReadFastArray(s, nch);
      s[nch] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Bool_t value to buffer

void TBufferXML::WriteBool(Bool_t b)
{
   BeforeIOoperation();
   XmlWriteBasic(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Char_t value to buffer

void TBufferXML::WriteChar(Char_t c)
{
   BeforeIOoperation();
   XmlWriteBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UChar_t value to buffer

void TBufferXML::WriteUChar(UChar_t c)
{
   BeforeIOoperation();
   XmlWriteBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Short_t value to buffer

void TBufferXML::WriteShort(Short_t h)
{
   BeforeIOoperation();
   XmlWriteBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UShort_t value to buffer

void TBufferXML::WriteUShort(UShort_t h)
{
   BeforeIOoperation();
   XmlWriteBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Int_t value to buffer

void TBufferXML::WriteInt(Int_t i)
{
   BeforeIOoperation();
   XmlWriteBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UInt_t value to buffer

void TBufferXML::WriteUInt(UInt_t i)
{
   BeforeIOoperation();
   XmlWriteBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Long_t value to buffer

void TBufferXML::WriteLong(Long_t l)
{
   BeforeIOoperation();
   XmlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes ULong_t value to buffer

void TBufferXML::WriteULong(ULong_t l)
{
   BeforeIOoperation();
   XmlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Long64_t value to buffer

void TBufferXML::WriteLong64(Long64_t l)
{
   BeforeIOoperation();
   XmlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes ULong64_t value to buffer

void TBufferXML::WriteULong64(ULong64_t l)
{
   BeforeIOoperation();
   XmlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Float_t value to buffer

void TBufferXML::WriteFloat(Float_t f)
{
   BeforeIOoperation();
   XmlWriteBasic(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Double_t value to buffer

void TBufferXML::WriteDouble(Double_t d)
{
   BeforeIOoperation();
   XmlWriteBasic(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes array of characters to buffer

void TBufferXML::WriteCharP(const Char_t *c)
{
   BeforeIOoperation();
   XmlWriteValue(c, xmlio::CharStar);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes a TString

void TBufferXML::WriteTString(const TString &s)
{
   if (GetIOVersion() < 3) {
      // original TBufferFile method, keep for compatibility
      Int_t nbig = s.Length();
      UChar_t nwh;
      if (nbig > 254) {
         nwh = 255;
         *this << nwh;
         *this << nbig;
      } else {
         nwh = UChar_t(nbig);
         *this << nwh;
      }
      const char *data = s.Data();
      WriteFastArray(data, nbig);
   } else {
      BeforeIOoperation();
      XmlWriteValue(s.Data(), xmlio::String);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Writes a std::string

void TBufferXML::WriteStdString(const std::string *obj)
{
   if (GetIOVersion() < 3) {
      if (!obj) {
         *this << (UChar_t)0;
         WriteFastArray("", 0);
         return;
      }

      UChar_t nwh;
      Int_t nbig = obj->length();
      if (nbig > 254) {
         nwh = 255;
         *this << nwh;
         *this << nbig;
      } else {
         nwh = UChar_t(nbig);
         *this << nwh;
      }
      WriteFastArray(obj->data(), nbig);
   } else {
      BeforeIOoperation();
      XmlWriteValue(obj ? obj->c_str() : "", xmlio::String);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a char* string

void TBufferXML::WriteCharStar(char *s)
{
   Int_t nch = 0;
   if (s) {
      nch = strlen(s);
      *this << nch;
      WriteFastArray(s, nch);
   } else {
      *this << nch;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Char_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Char_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%d", value);
   return XmlWriteValue(buf, xmlio::Char);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Short_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Short_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%hd", value);
   return XmlWriteValue(buf, xmlio::Short);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Int_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Int_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%d", value);
   return XmlWriteValue(buf, xmlio::Int);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Long_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Long_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%ld", value);
   return XmlWriteValue(buf, xmlio::Long);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Long64_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Long64_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), FLong64, value);
   return XmlWriteValue(buf, xmlio::Long64);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Float_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Float_t value)
{
   char buf[200];
   snprintf(buf, sizeof(buf), fgFloatFmt.c_str(), value);
   return XmlWriteValue(buf, xmlio::Float);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Double_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Double_t value)
{
   char buf[1000];
   snprintf(buf, sizeof(buf), fgFloatFmt.c_str(), value);
   return XmlWriteValue(buf, xmlio::Double);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Bool_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(Bool_t value)
{
   return XmlWriteValue(value ? xmlio::True : xmlio::False, xmlio::Bool);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts UChar_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(UChar_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%u", value);
   return XmlWriteValue(buf, xmlio::UChar);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts UShort_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(UShort_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%hu", value);
   return XmlWriteValue(buf, xmlio::UShort);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts UInt_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(UInt_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%u", value);
   return XmlWriteValue(buf, xmlio::UInt);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts ULong_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(ULong_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%lu", value);
   return XmlWriteValue(buf, xmlio::ULong);
}

////////////////////////////////////////////////////////////////////////////////
/// Converts ULong64_t to string and add xml node to buffer

XMLNodePointer_t TBufferXML::XmlWriteBasic(ULong64_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), FULong64, value);
   return XmlWriteValue(buf, xmlio::ULong64);
}

////////////////////////////////////////////////////////////////////////////////
/// Create xml node with specified name and adds it to stack node

XMLNodePointer_t TBufferXML::XmlWriteValue(const char *value, const char *name)
{
   XMLNodePointer_t node = 0;

   if (fCanUseCompact)
      node = StackNode();
   else
      node = CreateItemNode(name);

   fXML->NewAttr(node, 0, xmlio::v, value);

   fCanUseCompact = kFALSE;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Char_t value

void TBufferXML::XmlReadBasic(Char_t &value)
{
   const char *res = XmlReadValue(xmlio::Char);
   if (res) {
      int n;
      sscanf(res, "%d", &n);
      value = n;
   } else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Short_t value

void TBufferXML::XmlReadBasic(Short_t &value)
{
   const char *res = XmlReadValue(xmlio::Short);
   if (res)
      sscanf(res, "%hd", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Int_t value

void TBufferXML::XmlReadBasic(Int_t &value)
{
   const char *res = XmlReadValue(xmlio::Int);
   if (res)
      sscanf(res, "%d", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Long_t value

void TBufferXML::XmlReadBasic(Long_t &value)
{
   const char *res = XmlReadValue(xmlio::Long);
   if (res)
      sscanf(res, "%ld", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Long64_t value

void TBufferXML::XmlReadBasic(Long64_t &value)
{
   const char *res = XmlReadValue(xmlio::Long64);
   if (res)
      sscanf(res, FLong64, &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Float_t value

void TBufferXML::XmlReadBasic(Float_t &value)
{
   const char *res = XmlReadValue(xmlio::Float);
   if (res)
      sscanf(res, "%f", &value);
   else
      value = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Double_t value

void TBufferXML::XmlReadBasic(Double_t &value)
{
   const char *res = XmlReadValue(xmlio::Double);
   if (res)
      sscanf(res, "%lf", &value);
   else
      value = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to Bool_t value

void TBufferXML::XmlReadBasic(Bool_t &value)
{
   const char *res = XmlReadValue(xmlio::Bool);
   if (res)
      value = (strcmp(res, xmlio::True) == 0);
   else
      value = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to UChar_t value

void TBufferXML::XmlReadBasic(UChar_t &value)
{
   const char *res = XmlReadValue(xmlio::UChar);
   if (res) {
      unsigned int n;
      sscanf(res, "%ud", &n);
      value = n;
   } else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to UShort_t value

void TBufferXML::XmlReadBasic(UShort_t &value)
{
   const char *res = XmlReadValue(xmlio::UShort);
   if (res)
      sscanf(res, "%hud", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to UInt_t value

void TBufferXML::XmlReadBasic(UInt_t &value)
{
   const char *res = XmlReadValue(xmlio::UInt);
   if (res)
      sscanf(res, "%u", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to ULong_t value

void TBufferXML::XmlReadBasic(ULong_t &value)
{
   const char *res = XmlReadValue(xmlio::ULong);
   if (res)
      sscanf(res, "%lu", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reads string from current xml node and convert it to ULong64_t value

void TBufferXML::XmlReadBasic(ULong64_t &value)
{
   const char *res = XmlReadValue(xmlio::ULong64);
   if (res)
      sscanf(res, FULong64, &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// read string value from current stack node

const char *TBufferXML::XmlReadValue(const char *name)
{
   if (fErrorFlag > 0)
      return 0;

   Bool_t trysimple = fCanUseCompact;
   fCanUseCompact = kFALSE;

   if (trysimple) {
      if (fXML->HasAttr(Stack(1)->fNode, xmlio::v))
         fValueBuf = fXML->GetAttr(Stack(1)->fNode, xmlio::v);
      else
         trysimple = kFALSE;
   }

   if (!trysimple) {
      if (!VerifyItemNode(name, "XmlReadValue"))
         return 0;
      fValueBuf = fXML->GetAttr(StackNode(), xmlio::v);
   }

   if (gDebug > 4)
      Info("XmlReadValue", "     Name = %s value = %s", name, fValueBuf.Data());

   if (!trysimple)
      ShiftStack("readvalue");

   return fValueBuf.Data();
}

void TBufferXML::SetFloatFormat(const char *fmt)
{
   // Set printf format for float/double members, default "%e"
   // This method is not thread-safe as it changes a global state.

   if (!fmt)
      fgFloatFmt = "%e";
   fgFloatFmt = fmt;
}

const char *TBufferXML::GetFloatFormat()
{
   // return current printf format for float/double members, default "%e"

   return fgFloatFmt.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Read one collection of objects from the buffer using the StreamerInfoLoopAction.
/// The collection needs to be a split TClonesArray or a split vector of pointers.

Int_t TBufferXML::ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *obj)
{
   TVirtualStreamerInfo *info = sequence.fStreamerInfo;
   IncrementLevel(info);

   if (gDebug) {
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for (TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin(); iter != end;
           ++iter) {
         // Idea: Try to remove this function call as it is really needed only for XML streaming.
         SetStreamerElementNumber((*iter).fConfiguration->fCompInfo->fElem, (*iter).fConfiguration->fCompInfo->fType);
         (*iter).PrintDebug(*this, obj);
         (*iter)(*this, obj);
      }

   } else {
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for (TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin(); iter != end;
           ++iter) {
         // Idea: Try to remove this function call as it is really needed only for XML streaming.
         SetStreamerElementNumber((*iter).fConfiguration->fCompInfo->fElem, (*iter).fConfiguration->fCompInfo->fType);
         (*iter)(*this, obj);
      }
   }

   DecrementLevel(info);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read one collection of objects from the buffer using the StreamerInfoLoopAction.
/// The collection needs to be a split TClonesArray or a split vector of pointers.

Int_t TBufferXML::ApplySequenceVecPtr(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection,
                                      void *end_collection)
{
   TVirtualStreamerInfo *info = sequence.fStreamerInfo;
   IncrementLevel(info);

   if (gDebug) {
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for (TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin(); iter != end;
           ++iter) {
         // Idea: Try to remove this function call as it is really needed only for XML streaming.
         SetStreamerElementNumber((*iter).fConfiguration->fCompInfo->fElem, (*iter).fConfiguration->fCompInfo->fType);
         (*iter).PrintDebug(
            *this, *(char **)start_collection); // Warning: This limits us to TClonesArray and vector of pointers.
         (*iter)(*this, start_collection, end_collection);
      }

   } else {
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for (TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin(); iter != end;
           ++iter) {
         // Idea: Try to remove this function call as it is really needed only for XML streaming.
         SetStreamerElementNumber((*iter).fConfiguration->fCompInfo->fElem, (*iter).fConfiguration->fCompInfo->fType);
         (*iter)(*this, start_collection, end_collection);
      }
   }

   DecrementLevel(info);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read one collection of objects from the buffer using the StreamerInfoLoopAction.

Int_t TBufferXML::ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection,
                                void *end_collection)
{
   TVirtualStreamerInfo *info = sequence.fStreamerInfo;
   IncrementLevel(info);

   TStreamerInfoActions::TLoopConfiguration *loopconfig = sequence.fLoopConfig;
   if (gDebug) {

      // Get the address of the first item for the PrintDebug.
      // (Performance is not essential here since we are going to print to
      // the screen anyway).
      void *arr0 = loopconfig->GetFirstAddress(start_collection, end_collection);
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for (TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin(); iter != end;
           ++iter) {
         // Idea: Try to remove this function call as it is really needed only for XML streaming.
         SetStreamerElementNumber((*iter).fConfiguration->fCompInfo->fElem, (*iter).fConfiguration->fCompInfo->fType);
         (*iter).PrintDebug(*this, arr0);
         (*iter)(*this, start_collection, end_collection, loopconfig);
      }

   } else {
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for (TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin(); iter != end;
           ++iter) {
         // Idea: Try to remove this function call as it is really needed only for XML streaming.
         SetStreamerElementNumber((*iter).fConfiguration->fCompInfo->fElem, (*iter).fConfiguration->fCompInfo->fType);
         (*iter)(*this, start_collection, end_collection, loopconfig);
      }
   }

   DecrementLevel(info);
   return 0;
}

// abstract TBuffer methods, probably dedicated to have them in the TBufferText

////////////////////////////////////////////////////////////////////////////////
/// Check if the specified object of the specified class is already in
/// the buffer. Returns kTRUE if object already in the buffer,
/// kFALSE otherwise (also if obj is 0 ).

Bool_t TBufferXML::CheckObject(const TObject *obj)
{
   return CheckObject(obj, TObject::Class());
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the specified object of the specified class is already in
/// the buffer. Returns kTRUE if object already in the buffer,
/// kFALSE otherwise (also if obj is 0 ).

Bool_t TBufferXML::CheckObject(const void *obj, const TClass *ptrClass)
{
   if (!obj || !ptrClass || !fObjMap)
      return kFALSE;

   TClass *clActual = ptrClass->GetActualClass(obj);

   const char *temp = (const char *)obj;

   if (clActual && (ptrClass != clActual))
      temp -= clActual->GetBaseClassOffset(ptrClass);

   return fObjMap->GetValue(TString::Hash(&temp, sizeof(void *)), (Long_t)temp) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current streamer info element

TVirtualStreamerInfo *TBufferXML::GetInfo()
{
   return Stack()->fInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Deserialize information from a buffer into an object.
///
/// Note: This function is called by the xxx::Streamer() functions in
/// rootcint-generated dictionaries.
/// This function assumes that the class version and the byte count
/// information have been read.
///
/// \param[in] version The version number of the class
/// \param[in] start   The starting position in the buffer b
/// \param[in] count   The number of bytes for this object in the buffer
///

Int_t TBufferXML::ReadClassBuffer(const TClass *cl, void *pointer, Int_t version, UInt_t start, UInt_t count,
                                  const TClass *onFileClass)
{

   //---------------------------------------------------------------------------
   // The ondisk class has been specified so get foreign streamer info
   /////////////////////////////////////////////////////////////////////////////

   TStreamerInfo *sinfo = 0;
   if (onFileClass) {
      sinfo = (TStreamerInfo *)cl->GetConversionStreamerInfo(onFileClass, version);
      if (!sinfo) {
         Error("ReadClassBuffer",
               "Could not find the right streamer info to convert %s version %d into a %s, object skipped at offset %d",
               onFileClass->GetName(), version, cl->GetName(), Length());
         CheckByteCount(start, count, onFileClass);
         return 0;
      }
   }
   //---------------------------------------------------------------------------
   // Get local streamer info
   /////////////////////////////////////////////////////////////////////////////
   /// The StreamerInfo should exist at this point.

   else {
      R__LOCKGUARD(gInterpreterMutex);
      auto infos = cl->GetStreamerInfos();
      auto ninfos = infos->GetSize();
      if (version < -1 || version >= ninfos) {
         Error("ReadBuffer1", "class: %s, attempting to access a wrong version: %d, object skipped at offset %d",
               cl->GetName(), version, Length());
         CheckByteCount(start, count, cl);
         return 0;
      }
      sinfo = (TStreamerInfo *)infos->At(version);
      if (sinfo == 0) {
         // Unless the data is coming via a socket connection from with schema evolution
         // (tracking) was not enabled.  So let's create the StreamerInfo if it is the
         // one for the current version, otherwise let's complain ...
         // We could also get here if there old class version was '1' and the new class version is higher than 1
         // AND the checksum is the same.
         if (version == cl->GetClassVersion() || version == 1) {
            const_cast<TClass *>(cl)->BuildRealData(pointer);
            // This creation is alright since we just checked within the
            // current 'locked' section.
            sinfo = new TStreamerInfo(const_cast<TClass *>(cl));
            const_cast<TClass *>(cl)->RegisterStreamerInfo(sinfo);
            if (gDebug > 0)
               printf("Creating StreamerInfo for class: %s, version: %d\n", cl->GetName(), version);
            sinfo->Build();
         } else if (version == 0) {
            // When the object was written the class was version zero, so
            // there is no StreamerInfo to be found.
            // Check that the buffer position corresponds to the byte count.
            CheckByteCount(start, count, cl);
            return 0;
         } else {
            Error("ReadClassBuffer",
                  "Could not find the StreamerInfo for version %d of the class %s, object skipped at offset %d",
                  version, cl->GetName(), Length());
            CheckByteCount(start, count, cl);
            return 0;
         }
      } else if (!sinfo->IsCompiled()) { // Note this read is protected by the above lock.
         // Streamer info has not been compiled, but exists.
         // Therefore it was read in from a file and we have to do schema evolution.
         const_cast<TClass *>(cl)->BuildRealData(pointer);
         sinfo->BuildOld();
      }
   }

   // Deserialize the object.
   ApplySequence(*(sinfo->GetReadObjectWiseActions()), (char *)pointer);
   if (sinfo->IsRecovered())
      count = 0;

   // Check that the buffer position corresponds to the byte count.
   CheckByteCount(start, count, cl);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Deserialize information from a buffer into an object.
///
/// Note: This function is called by the xxx::Streamer()
/// functions in rootcint-generated dictionaries.
///

Int_t TBufferXML::ReadClassBuffer(const TClass *cl, void *pointer, const TClass *onFileClass)
{
   // Read the class version from the buffer.
   UInt_t R__s = 0; // Start of object.
   UInt_t R__c = 0; // Count of bytes.
   Version_t version;

   if (onFileClass)
      version = ReadVersion(&R__s, &R__c, onFileClass);
   else
      version = ReadVersion(&R__s, &R__c, cl);

   Bool_t v2file = kFALSE;
   TFile *file = (TFile *)GetParent();
   if (file && file->GetVersion() < 30000) {
      version = -1; // This is old file
      v2file = kTRUE;
   }

   //---------------------------------------------------------------------------
   // The ondisk class has been specified so get foreign streamer info
   /////////////////////////////////////////////////////////////////////////////

   TStreamerInfo *sinfo = 0;
   if (onFileClass) {
      sinfo = (TStreamerInfo *)cl->GetConversionStreamerInfo(onFileClass, version);
      if (!sinfo) {
         Error("ReadClassBuffer",
               "Could not find the right streamer info to convert %s version %d into a %s, object skipped at offset %d",
               onFileClass->GetName(), version, cl->GetName(), Length());
         CheckByteCount(R__s, R__c, onFileClass);
         return 0;
      }
   }
   //---------------------------------------------------------------------------
   // Get local streamer info
   /////////////////////////////////////////////////////////////////////////////
   /// The StreamerInfo should exist at this point.

   else {
      TStreamerInfo *guess = (TStreamerInfo *)cl->GetLastReadInfo();
      if (guess && guess->GetClassVersion() == version) {
         sinfo = guess;
      } else {
         // The last one is not the one we are looking for.
         {
            R__LOCKGUARD(gInterpreterMutex);

            const TObjArray *infos = cl->GetStreamerInfos();
            Int_t infocapacity = infos->Capacity();
            if (infocapacity) {
               if (version < -1 || version >= infocapacity) {
                  Error("ReadClassBuffer",
                        "class: %s, attempting to access a wrong version: %d, object skipped at offset %d",
                        cl->GetName(), version, Length());
                  CheckByteCount(R__s, R__c, cl);
                  return 0;
               }
               sinfo = (TStreamerInfo *)infos->UncheckedAt(version);
               if (sinfo) {
                  if (!sinfo->IsCompiled()) {
                     // Streamer info has not been compiled, but exists.
                     // Therefore it was read in from a file and we have to do schema evolution?
                     R__LOCKGUARD(gInterpreterMutex);
                     const_cast<TClass *>(cl)->BuildRealData(pointer);
                     sinfo->BuildOld();
                  }
                  // If the compilation succeeded, remember this StreamerInfo.
                  // const_cast okay because of the lock on gInterpreterMutex.
                  if (sinfo->IsCompiled())
                     const_cast<TClass *>(cl)->SetLastReadInfo(sinfo);
               }
            }
         }

         if (sinfo == 0) {
            // Unless the data is coming via a socket connection from with schema evolution
            // (tracking) was not enabled.  So let's create the StreamerInfo if it is the
            // one for the current version, otherwise let's complain ...
            // We could also get here when reading a file prior to the introduction of StreamerInfo.
            // We could also get here if there old class version was '1' and the new class version is higher than 1
            // AND the checksum is the same.
            if (v2file || version == cl->GetClassVersion() || version == 1) {
               R__LOCKGUARD(gInterpreterMutex);

               // We need to check if another thread did not get here first
               // and did the StreamerInfo creation already.
               auto infos = cl->GetStreamerInfos();
               auto ninfos = infos->GetSize();
               if (!(version < -1 || version >= ninfos)) {
                  sinfo = (TStreamerInfo *)infos->At(version);
               }
               if (!sinfo) {
                  const_cast<TClass *>(cl)->BuildRealData(pointer);
                  sinfo = new TStreamerInfo(const_cast<TClass *>(cl));
                  sinfo->SetClassVersion(version);
                  const_cast<TClass *>(cl)->RegisterStreamerInfo(sinfo);
                  if (gDebug > 0)
                     printf("Creating StreamerInfo for class: %s, version: %d\n", cl->GetName(), version);
                  if (v2file) {
                     sinfo->Build();             // Get the elements.
                     sinfo->Clear("build");      // Undo compilation.
                     sinfo->BuildEmulated(file); // Fix the types and redo compilation.
                  } else {
                     sinfo->Build();
                  }
               }
            } else if (version == 0) {
               // When the object was written the class was version zero, so
               // there is no StreamerInfo to be found.
               // Check that the buffer position corresponds to the byte count.
               CheckByteCount(R__s, R__c, cl);
               return 0;
            } else {
               Error("ReadClassBuffer",
                     "Could not find the StreamerInfo for version %d of the class %s, object skipped at offset %d",
                     version, cl->GetName(), Length());
               CheckByteCount(R__s, R__c, cl);
               return 0;
            }
         }
      }
   }

   // deserialize the object
   ApplySequence(*(sinfo->GetReadObjectWiseActions()), (char *)pointer);
   if (sinfo->TStreamerInfo::IsRecovered())
      R__c = 0; // 'TStreamerInfo::' avoids going via a virtual function.

   // Check that the buffer position corresponds to the byte count.
   CheckByteCount(R__s, R__c, cl);

   if (gDebug > 2)
      printf(" ReadBuffer for class: %s has read %d bytes\n", cl->GetName(), R__c);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Function called by the Streamer functions to serialize object at p
/// to buffer b. The optional argument info may be specified to give an
/// alternative StreamerInfo instead of using the default StreamerInfo
/// automatically built from the class definition.
/// For more information, see class TStreamerInfo.

Int_t TBufferXML::WriteClassBuffer(const TClass *cl, void *pointer)
{
   // build the StreamerInfo if first time for the class
   TStreamerInfo *sinfo = (TStreamerInfo *)const_cast<TClass *>(cl)->GetCurrentStreamerInfo();
   if (sinfo == 0) {
      // Have to be sure between the check and the taking of the lock if the current streamer has changed
      R__LOCKGUARD(gInterpreterMutex);
      sinfo = (TStreamerInfo *)const_cast<TClass *>(cl)->GetCurrentStreamerInfo();
      if (sinfo == 0) {
         const_cast<TClass *>(cl)->BuildRealData(pointer);
         sinfo = new TStreamerInfo(const_cast<TClass *>(cl));
         const_cast<TClass *>(cl)->SetCurrentStreamerInfo(sinfo);
         const_cast<TClass *>(cl)->RegisterStreamerInfo(sinfo);
         if (gDebug > 0)
            printf("Creating StreamerInfo for class: %s, version: %d\n", cl->GetName(), cl->GetClassVersion());
         sinfo->Build();
      }
   } else if (!sinfo->IsCompiled()) {
      R__LOCKGUARD(gInterpreterMutex);
      // Redo the test in case we have been victim of a data race on fIsCompiled.
      if (!sinfo->IsCompiled()) {
         const_cast<TClass *>(cl)->BuildRealData(pointer);
         sinfo->BuildOld();
      }
   }

   // write the class version number and reserve space for the byte count
   UInt_t R__c = WriteVersion(cl, kTRUE);

   // NOTE: In the future Philippe wants this to happen via a custom action
   TagStreamerInfo(sinfo);
   ApplySequence(*(sinfo->GetWriteObjectWiseActions()), (char *)pointer);

   // write the byte count at the start of the buffer
   SetByteCount(R__c, kTRUE);

   if (gDebug > 2)
      printf(" WriteBuffer for class: %s version %d has written %d bytes\n", cl->GetName(), cl->GetClassVersion(),
             UInt_t(fBufCur - fBuffer) - R__c - (UInt_t)sizeof(UInt_t));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Mark the classindex of the current file as using this TStreamerInfo

void TBufferXML::TagStreamerInfo(TVirtualStreamerInfo *info)
{
   TFile *file = (TFile *)GetParent();
   if (file) {
      TArrayC *cindex = file->GetClassIndex();
      Int_t nindex = cindex->GetSize();
      Int_t number = info->GetNumber();
      if (number < 0 || number >= nindex) {
         Error("TagStreamerInfo", "StreamerInfo: %s number: %d out of range[0,%d] in file: %s", info->GetName(), number,
               nindex, file->GetName());
         return;
      }
      if (cindex->fArray[number] == 0) {
         cindex->fArray[0] = 1;
         cindex->fArray[number] = 1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert object into json structures.
/// !!! Should be used only by TBufferJSON itself.
/// Use ConvertToJSON() methods to convert object to json

void TBufferXML::WriteObject(const TObject *obj, Bool_t cacheReuse /* = kTRUE */)
{
   WriteObjectAny(obj, TObject::Class(), cacheReuse);
}

namespace {
struct DynamicType {
   // Helper class to enable typeid on any address
   // Used in code similar to:
   //    typeid( * (DynamicType*) void_ptr );
   virtual ~DynamicType() {}
};
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to I/O buffer.
/// This function assumes that the value in 'obj' is the value stored in
/// a pointer to a "ptrClass". The actual type of the object pointed to
/// can be any class derived from "ptrClass".
/// Return:
///  - 0: failure
///  - 1: success
///  - 2: truncated success (i.e actual class is missing. Only ptrClass saved.)
///
/// If 'cacheReuse' is true (default) upon seeing an object address a second time,
/// we record the offset where its was written the first time rather than streaming
/// the object a second time.
/// If 'cacheReuse' is false, we always stream the object.  This allows the (re)use
/// of temporary object to store different data in the same buffer.

Int_t TBufferXML::WriteObjectAny(const void *obj, const TClass *ptrClass, Bool_t cacheReuse /* = kTRUE */)
{
   if (!obj) {
      WriteObjectClass(0, 0, kTRUE);
      return 1;
   }

   if (!ptrClass) {
      Error("WriteObjectAny", "ptrClass argument may not be 0");
      return 0;
   }

   TClass *clActual = ptrClass->GetActualClass(obj);

   if (!clActual) {
      // The ptrClass is a class with a virtual table and we have no
      // TClass with the actual type_info in memory.

      DynamicType *d_ptr = (DynamicType *)obj;
      Warning("WriteObjectAny", "An object of type %s (from type_info) passed through a %s pointer was truncated (due "
                                "a missing dictionary)!!!",
              typeid(*d_ptr).name(), ptrClass->GetName());
      WriteObjectClass(obj, ptrClass, cacheReuse);
      return 2;
   } else if (clActual && (clActual != ptrClass)) {
      const char *temp = (const char *)obj;
      temp -= clActual->GetBaseClassOffset(ptrClass);
      WriteObjectClass(temp, clActual, cacheReuse);
      return 1;
   } else {
      WriteObjectClass(obj, ptrClass, cacheReuse);
      return 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// force writing the TStreamerInfo to the file

void TBufferXML::ForceWriteInfo(TVirtualStreamerInfo *info, Bool_t force)
{
   if (info)
      info->ForceWriteInfo((TFile *)GetParent(), force);
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure TStreamerInfo is not optimized, otherwise it will not be
/// possible to support schema evolution in read mode.
/// In case the StreamerInfo has already been computed and optimized,
/// one must disable the option BypassStreamer.

void TBufferXML::ForceWriteInfoClones(TClonesArray *a)
{
   TStreamerInfo *sinfo = (TStreamerInfo *)a->GetClass()->GetStreamerInfo();
   ForceWriteInfo(sinfo, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to TStreamerInfo::ReadBufferClones.

Int_t TBufferXML::ReadClones(TClonesArray *a, Int_t nobjects, Version_t objvers)
{
   char **arr = (char **)a->GetObjectRef(0);
   char **end = arr + nobjects;
   // a->GetClass()->GetStreamerInfo()->ReadBufferClones(*this,a,nobjects,-1,0);
   TStreamerInfo *info = (TStreamerInfo *)a->GetClass()->GetStreamerInfo(objvers);
   // return info->ReadBuffer(*this,arr,-1,nobjects,0,1);
   return ApplySequenceVecPtr(*(info->GetReadMemberWiseActions(kTRUE)), arr, end);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to TStreamerInfo::WriteBufferClones.

Int_t TBufferXML::WriteClones(TClonesArray *a, Int_t nobjects)
{
   char **arr = reinterpret_cast<char **>(a->GetObjectRef(0));
   // a->GetClass()->GetStreamerInfo()->WriteBufferClones(*this,(TClonesArray*)a,nobjects,-1,0);
   TStreamerInfo *info = (TStreamerInfo *)a->GetClass()->GetStreamerInfo();
   // return info->WriteBufferAux(*this,arr,-1,nobjects,0,1);
   char **end = arr + nobjects;
   // No need to tell call ForceWriteInfo as it by ForceWriteInfoClones.
   return ApplySequenceVecPtr(*(info->GetWriteMemberWiseActions(kTRUE)), arr, end);
}
