// $Id$
// Author: Sergey Linev   22/12/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TRootSnifferFull.h"

#include "TH1.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TList.h"
#include "TMemFile.h"
#include "TBufferFile.h"
#include "TBufferJSON.h"
#include "TBufferXML.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TClass.h"
#include "TMethod.h"
#include "TFunction.h"
#include "TMethodArg.h"
#include "TMethodCall.h"
#include "TUrl.h"
#include "TImage.h"
#include "TVirtualMutex.h"
#include "TRootSnifferStore.h"
#include "THttpCallArg.h"

#include <stdlib.h>
#include <string.h>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootSnifferFull                                                     //
//                                                                      //
// Subclass of TRootSniffer, which provides access to different         //
// ROOT collections and containers like TTree, TCanvas, ...             //
//////////////////////////////////////////////////////////////////////////

ClassImp(TRootSnifferFull);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TRootSnifferFull::TRootSnifferFull(const char *name, const char *objpath) : TRootSniffer(name, objpath)
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TRootSnifferFull::~TRootSnifferFull()
{
   delete fSinfo;

   delete fMemFile;
}

////////////////////////////////////////////////////////////////////////////////
/// return true if given class can be drawn in JSROOT

Bool_t TRootSnifferFull::IsDrawableClass(TClass *cl)
{
   if (!cl)
      return kFALSE;
   if (cl->InheritsFrom(TH1::Class()))
      return kTRUE;
   if (cl->InheritsFrom(TGraph::Class()))
      return kTRUE;
   if (cl->InheritsFrom(TCanvas::Class()))
      return kTRUE;
   if (cl->InheritsFrom(TProfile::Class()))
      return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// scans object properties
/// here such fields as _autoload or _icon properties depending on class or object name could be assigned
/// By default properties, coded in the Class title are scanned. Example:
///   ClassDef(UserClassName, 1) //  class comments *SNIFF*  _field1=value _field2="string value"
/// Here *SNIFF* mark is important. After it all expressions like field=value are parsed
/// One could use double quotes to code string values with spaces.
/// Fields separated from each other with spaces

void TRootSnifferFull::ScanObjectProperties(TRootSnifferScanRec &rec, TObject *obj)
{
   if (obj && obj->InheritsFrom(TLeaf::Class())) {
      rec.SetField("_more", "false", kFALSE);
      rec.SetField("_can_draw", "false", kFALSE);
      rec.SetField("_player", "JSROOT.drawLeafPlayer");
      rec.SetField("_prereq", "jq2d");
      return;
   }

   TRootSniffer::ScanObjectProperties(rec, obj);
}

////////////////////////////////////////////////////////////////////////////////
/// scans key properties
/// in special cases load objects from the file

void TRootSnifferFull::ScanKeyProperties(TRootSnifferScanRec &rec, TKey *key, TObject *&obj, TClass *&obj_class)
{
   if (strcmp(key->GetClassName(), "TDirectoryFile") == 0) {
      TRootSniffer::ScanKeyProperties(rec, key, obj, obj_class);
   } else {
      obj_class = TClass::GetClass(key->GetClassName());
      if (obj_class && obj_class->InheritsFrom(TTree::Class())) {
         if (rec.CanExpandItem()) {
            // it is requested to expand tree element - read it
            obj = key->ReadObj();
            if (obj)
               obj_class = obj->IsA();
         } else {
            rec.SetField("_ttree", "true", kFALSE); // indicate ROOT TTree
            rec.SetField("_player", "JSROOT.drawTreePlayerKey");
            rec.SetField("_prereq", "jq2d");
            // rec.SetField("_more", "true", kFALSE); // one could allow to extend
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// scans object childs (if any)
/// here one scans collection, branches, trees and so on

void TRootSnifferFull::ScanObjectChilds(TRootSnifferScanRec &rec, TObject *obj)
{
   if (obj->InheritsFrom(TTree::Class())) {
      if (!rec.IsReadOnly(fReadOnly)) {
         rec.SetField("_ttree", "true", kFALSE); // indicate ROOT TTree
         rec.SetField("_player", "JSROOT.drawTreePlayer");
         rec.SetField("_prereq", "jq2d");
      }
      ScanCollection(rec, ((TTree *)obj)->GetListOfLeaves());
   } else if (obj->InheritsFrom(TBranch::Class())) {
      ScanCollection(rec, ((TBranch *)obj)->GetListOfLeaves());
   } else {
      TRootSniffer::ScanObjectChilds(rec, obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns hash value for streamer infos
/// At the moment - just number of items in streamer infos list.

ULong_t TRootSnifferFull::GetStreamerInfoHash()
{
   return fSinfo ? fSinfo->GetSize() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if it is streamer info item name

Bool_t TRootSnifferFull::IsStreamerInfoItem(const char *itemname)
{
   if (!itemname || (*itemname == 0))
      return kFALSE;

   return (strcmp(itemname, "StreamerInfo") == 0) || (strcmp(itemname, "StreamerInfo/") == 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Get hash function for specified item
/// used to detect any changes in the specified object

ULong_t TRootSnifferFull::GetItemHash(const char *itemname)
{
   if (IsStreamerInfoItem(itemname))
      return GetStreamerInfoHash();

   return TRootSniffer::GetItemHash(itemname);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates TMemFile instance, which used for objects streaming
/// One could not use TBufferFile directly,
/// while one also require streamer infos list

void TRootSnifferFull::CreateMemFile()
{
   if (fMemFile)
      return;

   TDirectory *olddir = gDirectory;
   gDirectory = nullptr;
   TFile *oldfile = gFile;
   gFile = nullptr;

   fMemFile = new TMemFile("dummy.file", "RECREATE");
   gROOT->GetListOfFiles()->Remove(fMemFile);

   TH1F *d = new TH1F("d", "d", 10, 0, 10);
   fMemFile->WriteObject(d, "h1");
   delete d;

   TGraph *gr = new TGraph(10);
   gr->SetName("abc");
   //      // gr->SetDrawOptions("AC*");
   fMemFile->WriteObject(gr, "gr1");
   delete gr;

   fMemFile->WriteStreamerInfo();

   // make primary list of streamer infos
   TList *l = new TList();

   l->Add(gROOT->GetListOfStreamerInfo()->FindObject("TGraph"));
   l->Add(gROOT->GetListOfStreamerInfo()->FindObject("TH1F"));
   l->Add(gROOT->GetListOfStreamerInfo()->FindObject("TH1"));
   l->Add(gROOT->GetListOfStreamerInfo()->FindObject("TNamed"));
   l->Add(gROOT->GetListOfStreamerInfo()->FindObject("TObject"));

   fMemFile->WriteObject(l, "ll");
   delete l;

   fMemFile->WriteStreamerInfo();

   fSinfo = fMemFile->GetStreamerInfoList();

   gDirectory = olddir;
   gFile = oldfile;
}

////////////////////////////////////////////////////////////////////////////////
/// Search element with specified path
/// Returns pointer on element
/// Optionally one could obtain element class, member description
/// and number of childs. When chld!=0, not only element is searched,
/// but also number of childs are counted. When member!=0, any object
/// will be scanned for its data members (disregard of extra options)

void *TRootSnifferFull::FindInHierarchy(const char *path, TClass **cl, TDataMember **member, Int_t *chld)
{
   if (IsStreamerInfoItem(path)) {
      // special handling for streamer info
      CreateMemFile();
      if (cl && fSinfo)
         *cl = fSinfo->IsA();
      return fSinfo;
   }

   return TRootSniffer::FindInHierarchy(path, cl, member, chld);
}

////////////////////////////////////////////////////////////////////////////////
/// produce binary data for specified item
/// if "zipped" option specified in query, buffer will be compressed

Bool_t TRootSnifferFull::ProduceBinary(const std::string &path, const std::string & /*query*/, std::string &res)
{
   if (path.empty())
      return kFALSE;

   const char *path_ = path.c_str();
   if (*path_ == '/')
      path_++;

   TClass *obj_cl = nullptr;
   void *obj_ptr = FindInHierarchy(path_, &obj_cl);
   if (!obj_ptr || !obj_cl)
      return kFALSE;

   if (obj_cl->GetBaseClassOffset(TObject::Class()) != 0) {
      Info("ProduceBinary", "Non-TObject class not supported");
      return kFALSE;
   }

   // ensure that memfile exists
   CreateMemFile();

   TDirectory *olddir = gDirectory;
   gDirectory = nullptr;
   TFile *oldfile = gFile;
   gFile = nullptr;

   TObject *obj = (TObject *)obj_ptr;

   TBufferFile *sbuf = new TBufferFile(TBuffer::kWrite, 100000);
   sbuf->SetParent(fMemFile);
   sbuf->MapObject(obj);
   obj->Streamer(*sbuf);
   if (fCurrentArg)
      fCurrentArg->SetExtraHeader("RootClassName", obj_cl->GetName());

   // produce actual version of streamer info
   delete fSinfo;
   fMemFile->WriteStreamerInfo();
   fSinfo = fMemFile->GetStreamerInfoList();

   gDirectory = olddir;
   gFile = oldfile;

   res.resize(sbuf->Length());
   std::copy((const char *)sbuf->Buffer(), (const char *)sbuf->Buffer() + sbuf->Length(), res.begin());

   delete sbuf;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Method to produce image from specified object
///
/// Parameters:
///    kind - image kind TImage::kPng, TImage::kJpeg, TImage::kGif
///    path - path to object
///    options - extra options
///
/// By default, image 300x200 is produced
/// In options string one could provide following parameters:
///    w - image width
///    h - image height
///    opt - draw options
///  For instance:
///     http://localhost:8080/Files/hsimple.root/hpx/get.png?w=500&h=500&opt=lego1
///
///  Return is std::string with binary data

Bool_t TRootSnifferFull::ProduceImage(Int_t kind, const std::string &path, const std::string &options, std::string &res)
{
   if (path.empty())
      return kFALSE;

   const char *path_ = path.c_str();
   if (*path_ == '/')
      path_++;

   TClass *obj_cl(nullptr);
   void *obj_ptr = FindInHierarchy(path_, &obj_cl);
   if (!obj_ptr || !obj_cl)
      return kFALSE;

   if (obj_cl->GetBaseClassOffset(TObject::Class()) != 0) {
      Error("TRootSniffer", "Only derived from TObject classes can be drawn");
      return kFALSE;
   }

   TObject *obj = (TObject *)obj_ptr;

   TImage *img = TImage::Create();
   if (!img)
      return kFALSE;

   if (obj->InheritsFrom(TPad::Class())) {

      if (gDebug > 1)
         Info("TRootSniffer", "Crate IMAGE directly from pad");
      img->FromPad((TPad *)obj);
   } else if (CanDrawClass(obj->IsA())) {

      if (gDebug > 1)
         Info("TRootSniffer", "Crate IMAGE from object %s", obj->GetName());

      Int_t width(300), height(200);
      TString drawopt;

      if (!options.empty()) {
         TUrl url;
         url.SetOptions(options.c_str());
         url.ParseOptions();
         Int_t w = url.GetIntValueFromOptions("w");
         if (w > 10)
            width = w;
         Int_t h = url.GetIntValueFromOptions("h");
         if (h > 10)
            height = h;
         const char *opt = url.GetValueFromOptions("opt");
         if (opt)
            drawopt = opt;
      }

      Bool_t isbatch = gROOT->IsBatch();
      TVirtualPad *save_gPad = gPad;

      if (!isbatch)
         gROOT->SetBatch(kTRUE);

      TCanvas *c1 = new TCanvas("__online_draw_canvas__", "title", width, height);
      obj->Draw(drawopt.Data());
      img->FromPad(c1);
      delete c1;

      if (!isbatch)
         gROOT->SetBatch(kFALSE);
      gPad = save_gPad;

   } else {
      delete img;
      return kFALSE;
   }

   TImage *im = TImage::Create();
   im->Append(img);

   char *png_buffer = nullptr;
   int size(0);

   im->GetImageBuffer(&png_buffer, &size, (TImage::EImageFileTypes)kind);

   if (png_buffer && (size > 0)) {
      res.resize(size);
      memcpy((void *)res.data(), png_buffer, size);
   }

   delete[] png_buffer;
   delete im;

   return !res.empty();
}

////////////////////////////////////////////////////////////////////////////////
/// produce XML data for specified item
/// For object conversion TBufferXML is used

Bool_t TRootSnifferFull::ProduceXml(const std::string &path, const std::string & /*options*/, std::string &res)
{
   if (path.empty())
      return kFALSE;
   const char *path_ = path.c_str();
   if (*path_ == '/')
      path_++;

   TClass *obj_cl = nullptr;
   void *obj_ptr = FindInHierarchy(path_, &obj_cl);
   if (!obj_ptr || !obj_cl)
      return kFALSE;

   // TODO: support std::string in TBufferXML
   res = TBufferXML::ConvertToXML(obj_ptr, obj_cl).Data();

   return !res.empty();
}

////////////////////////////////////////////////////////////////////////////////
/// execute command for specified object
/// options include method and extra list of parameters
/// sniffer should be not-readonly to allow execution of the commands
/// reskind defines kind of result 0 - debug, 1 - json, 2 - binary

Bool_t TRootSnifferFull::ProduceExe(const std::string &path, const std::string &options, Int_t reskind, std::string &res_str)
{
   std::string *debug = (reskind == 0) ? &res_str : nullptr;

   if (path.empty()) {
      if (debug)
         debug->append("Item name not specified\n");
      return debug != nullptr;
   }

   const char *path_ = path.c_str();
   if (*path_ == '/')
      path_++;

   TClass *obj_cl = nullptr;
   void *obj_ptr = FindInHierarchy(path_, &obj_cl);
   if (debug)
      debug->append(Form("Item:%s found:%s\n", path_, obj_ptr ? "true" : "false"));
   if (!obj_ptr || !obj_cl)
      return debug != nullptr;

   TUrl url;
   url.SetOptions(options.c_str());

   const char *method_name = url.GetValueFromOptions("method");
   TString prototype = DecodeUrlOptionValue(url.GetValueFromOptions("prototype"), kTRUE);
   TString funcname = DecodeUrlOptionValue(url.GetValueFromOptions("func"), kTRUE);
   TMethod *method = nullptr;
   TFunction *func = nullptr;
   if (method_name != nullptr) {
      if (prototype.Length() == 0) {
         if (debug)
            debug->append(Form("Search for any method with name \'%s\'\n", method_name));
         method = obj_cl->GetMethodAllAny(method_name);
      } else {
         if (debug)
            debug->append(Form("Search for method \'%s\' with prototype \'%s\'\n", method_name, prototype.Data()));
         method = obj_cl->GetMethodWithPrototype(method_name, prototype);
      }
   }

   if (method) {
      if (debug)
         debug->append(Form("Method: %s\n", method->GetPrototype()));
   } else {
      if (funcname.Length() > 0) {
         if (prototype.Length() == 0) {
            if (debug)
               debug->append(Form("Search for any function with name \'%s\'\n", funcname.Data()));
            func = gROOT->GetGlobalFunction(funcname);
         } else {
            if (debug)
               debug->append(
                  Form("Search for function \'%s\' with prototype \'%s\'\n", funcname.Data(), prototype.Data()));
            func = gROOT->GetGlobalFunctionWithPrototype(funcname, prototype);
         }
      }

      if (func) {
         if (debug)
            debug->append(Form("Function: %s\n", func->GetPrototype()));
      }
   }

   if (!method && !func) {
      if (debug)
         debug->append("Method not found\n");
      return debug != nullptr;
   }

   if ((fReadOnly && (fCurrentRestrict == 0)) || (fCurrentRestrict == 1)) {
      if ((method != nullptr) && (fCurrentAllowedMethods.Index(method_name) == kNPOS)) {
         if (debug)
            debug->append("Server runs in read-only mode, method cannot be executed\n");
         return debug != nullptr;
      } else if ((func != nullptr) && (fCurrentAllowedMethods.Index(funcname) == kNPOS)) {
         if (debug)
            debug->append("Server runs in read-only mode, function cannot be executed\n");
         return debug != nullptr;
      } else {
         if (debug)
            debug->append("For that special method server allows access even read-only mode is specified\n");
      }
   }

   TList *args = method ? method->GetListOfMethodArgs() : func->GetListOfMethodArgs();

   TList garbage;
   garbage.SetOwner(kTRUE);     // use as garbage collection
   TObject *post_obj = nullptr; // object reconstructed from post request
   TString call_args;

   TIter next(args);
   TMethodArg *arg = nullptr;
   while ((arg = (TMethodArg *)next()) != nullptr) {

      if ((strcmp(arg->GetName(), "rest_url_opt") == 0) && (strcmp(arg->GetFullTypeName(), "const char*") == 0) &&
          (args->GetSize() == 1)) {
         // very special case - function requires list of options after method=argument

         const char *pos = strstr(options.c_str(), "method=");
         if (!pos || (strlen(pos) < strlen(method_name) + 7))
            return debug != nullptr;
         const char *rest_url = pos + strlen(method_name) + 7;
         if (*rest_url == '&') ++rest_url;
         call_args.Form("\"%s\"", rest_url);
         break;
      }

      TString sval;
      const char *val = url.GetValueFromOptions(arg->GetName());
      if (val) {
         sval = DecodeUrlOptionValue(val, kFALSE);
         val = sval.Data();
      }

      if ((val != nullptr) && (strcmp(val, "_this_") == 0)) {
         // special case - object itself is used as argument
         sval.Form("(%s*)0x%lx", obj_cl->GetName(), (long unsigned)obj_ptr);
         val = sval.Data();
      } else if ((val != nullptr) && (fCurrentArg != nullptr) && (fCurrentArg->GetPostData() != nullptr)) {
         // process several arguments which are specific for post requests
         if (strcmp(val, "_post_object_xml_") == 0) {
            // post data has extra 0 at the end and can be used as null-terminated string
            post_obj = TBufferXML::ConvertFromXML((const char *)fCurrentArg->GetPostData());
            if (!post_obj) {
               sval = "0";
            } else {
               sval.Form("(%s*)0x%lx", post_obj->ClassName(), (long unsigned)post_obj);
               if (url.HasOption("_destroy_post_"))
                  garbage.Add(post_obj);
            }
            val = sval.Data();
         } else if (strcmp(val, "_post_object_json_") == 0) {
            // post data has extra 0 at the end and can be used as null-terminated string
            post_obj = TBufferJSON::ConvertFromJSON((const char *)fCurrentArg->GetPostData());
            if (!post_obj) {
               sval = "0";
            } else {
               sval.Form("(%s*)0x%lx", post_obj->ClassName(), (long unsigned)post_obj);
               if (url.HasOption("_destroy_post_"))
                  garbage.Add(post_obj);
            }
            val = sval.Data();
         } else if ((strcmp(val, "_post_object_") == 0) && url.HasOption("_post_class_")) {
            TString clname = url.GetValueFromOptions("_post_class_");
            TClass *arg_cl = gROOT->GetClass(clname, kTRUE, kTRUE);
            if ((arg_cl != nullptr) && (arg_cl->GetBaseClassOffset(TObject::Class()) == 0) && (post_obj == nullptr)) {
               post_obj = (TObject *)arg_cl->New();
               if (post_obj == nullptr) {
                  if (debug)
                     debug->append(TString::Format("Fail to create object of class %s\n", clname.Data()).Data());
               } else {
                  if (debug)
                     debug->append(TString::Format("Reconstruct object of class %s from POST data\n", clname.Data()).Data());
                  TBufferFile buf(TBuffer::kRead, fCurrentArg->GetPostDataLength(), (void *)fCurrentArg->GetPostData(), kFALSE);
                  buf.MapObject(post_obj, arg_cl);
                  post_obj->Streamer(buf);
                  if (url.HasOption("_destroy_post_"))
                     garbage.Add(post_obj);
               }
            }
            sval.Form("(%s*)0x%lx", clname.Data(), (long unsigned)post_obj);
            val = sval.Data();
         } else if (strcmp(val, "_post_data_") == 0) {
            sval.Form("(void*)0x%lx", (long unsigned)fCurrentArg->GetPostData());
            val = sval.Data();
         } else if (strcmp(val, "_post_length_") == 0) {
            sval.Form("%ld", (long)fCurrentArg->GetPostDataLength());
            val = sval.Data();
         }
      }

      if (!val)
         val = arg->GetDefault();

      if (debug)
         debug->append(Form("  Argument:%s Type:%s Value:%s \n", arg->GetName(), arg->GetFullTypeName(),
                                       val ? val : "<missed>"));
      if (!val)
         return debug != nullptr;

      if (call_args.Length() > 0)
         call_args += ", ";

      if ((strcmp(arg->GetFullTypeName(), "const char*") == 0) || (strcmp(arg->GetFullTypeName(), "Option_t*") == 0)) {
         int len = strlen(val);
         if ((strlen(val) < 2) || (*val != '\"') || (val[len - 1] != '\"'))
            call_args.Append(TString::Format("\"%s\"", val));
         else
            call_args.Append(val);
      } else {
         call_args.Append(val);
      }
   }

   TMethodCall *call = nullptr;

   if (method != nullptr) {
      call = new TMethodCall(obj_cl, method_name, call_args.Data());
      if (debug)
         debug->append(Form("Calling obj->%s(%s);\n", method_name, call_args.Data()));
   } else {
      call = new TMethodCall(funcname.Data(), call_args.Data());
      if (debug)
         debug->append(Form("Calling %s(%s);\n", funcname.Data(), call_args.Data()));
   }

   garbage.Add(call);

   if (!call->IsValid()) {
      if (debug)
         debug->append("Fail: invalid TMethodCall\n");
      return debug != nullptr;
   }

   Int_t compact = 0;
   if (url.GetValueFromOptions("compact"))
      compact = url.GetIntValueFromOptions("compact");

   TString res = "null";
   void *ret_obj = nullptr;
   TClass *ret_cl = nullptr;
   TBufferFile *resbuf = nullptr;
   if (reskind == 2) {
      resbuf = new TBufferFile(TBuffer::kWrite, 10000);
      garbage.Add(resbuf);
   }

   switch (call->ReturnType()) {
   case TMethodCall::kLong: {
      Long_t l(0);
      if (method)
         call->Execute(obj_ptr, l);
      else
         call->Execute(l);
      if (resbuf)
         resbuf->WriteLong(l);
      else
         res.Form("%ld", l);
      break;
   }
   case TMethodCall::kDouble: {
      Double_t d(0.);
      if (method)
         call->Execute(obj_ptr, d);
      else
         call->Execute(d);
      if (resbuf)
         resbuf->WriteDouble(d);
      else
         res.Form(TBufferJSON::GetFloatFormat(), d);
      break;
   }
   case TMethodCall::kString: {
      char *txt = nullptr;
      if (method)
         call->Execute(obj_ptr, &txt);
      else
         call->Execute(0, &txt); // here 0 is artificial, there is no proper signature
      if (txt != nullptr) {
         if (resbuf)
            resbuf->WriteString(txt);
         else
            res.Form("\"%s\"", txt);
      }
      break;
   }
   case TMethodCall::kOther: {
      std::string ret_kind = func ? func->GetReturnTypeNormalizedName() : method->GetReturnTypeNormalizedName();
      if ((ret_kind.length() > 0) && (ret_kind[ret_kind.length() - 1] == '*')) {
         ret_kind.resize(ret_kind.length() - 1);
         ret_cl = gROOT->GetClass(ret_kind.c_str(), kTRUE, kTRUE);
      }

      if (ret_cl != nullptr) {
         Long_t l(0);
         if (method)
            call->Execute(obj_ptr, l);
         else
            call->Execute(l);
         if (l != 0)
            ret_obj = (void *)l;
      } else {
         if (method)
            call->Execute(obj_ptr);
         else
            call->Execute();
      }

      break;
   }
   case TMethodCall::kNone: {
      if (method)
         call->Execute(obj_ptr);
      else
         call->Execute();
      break;
   }
   }

   const char *_ret_object_ = url.GetValueFromOptions("_ret_object_");
   if (_ret_object_ != nullptr) {
      TObject *obj = nullptr;
      if (gDirectory)
         obj = gDirectory->Get(_ret_object_);
      if (debug)
         debug->append(Form("Return object %s found %s\n", _ret_object_, obj ? "true" : "false"));

      if (obj == nullptr) {
         res = "null";
      } else {
         ret_obj = obj;
         ret_cl = obj->IsA();
      }
   }

   if (ret_obj && ret_cl) {
      if ((resbuf != nullptr) && (ret_cl->GetBaseClassOffset(TObject::Class()) == 0)) {
         TObject *obj = (TObject *)ret_obj;
         resbuf->MapObject(obj);
         obj->Streamer(*resbuf);
         if (fCurrentArg)
            fCurrentArg->SetExtraHeader("RootClassName", ret_cl->GetName());
      } else {
         res = TBufferJSON::ConvertToJSON(ret_obj, ret_cl, compact);
      }
   }

   if ((resbuf != nullptr) && (resbuf->Length() > 0)) {
      res_str.resize(resbuf->Length());
      std::copy((const char *)resbuf->Buffer(), (const char *)resbuf->Buffer() + resbuf->Length(), res_str.begin());
   }

   if (debug)
      debug->append(Form("Result = %s\n", res.Data()));

   if (reskind == 1)
      res_str = res.Data();

   if (url.HasOption("_destroy_result_") && ret_obj && ret_cl) {
      ret_cl->Destructor(ret_obj);
      if (debug)
         debug->append("Destroy result object at the end\n");
   }

   // delete all garbage objects, but should be also done with any return
   garbage.Delete();

   return kTRUE;
}
