// @(#)root/xml:$Id$
// Author: Sergey Linev, Rene Brun  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//________________________________________________________________________
//
// Class for xml code generation
// It should be used for generation of xml streamers, which could be used outside root
// environment. This means, that with help of such streamers user can read and write
// objects from/to xml file, which later can be accepted by ROOT.
//
// At the moment supported only classes, which are not inherited from TObject
// and which not contains any TObject members.
//
// To generate xml code:
//
// 1. ROOT library with required classes should be created.
//    In general, without such library non of user objects can be stored and
//    retrieved from any ROOT file
//
// 2. Generate xml streamers by root script like:
//
//    void generate() {
//      gSystem->Load("libRXML.so");   // load ROOT xml library
//      gSystem->Load("libuser.so");   // load user ROOT library
//
//      TList lst;
//      lst.Add(TClass::GetClass("TUserClass1"));
//      lst.Add(TClass::GetClass("TUserClass2"));
//      ...
//      TXMLPlayer player;
//      player.ProduceCode(&lst, "streamers");    // create xml streamers
//    }
//
//  3. Copy "streamers.h", "streamers.cxx", "TXmlFile.h", "TXmlFile.cxx" files
//     to user project and compile them. TXmlFile class implementation can be taken
//     from http://web-docs.gsi.de/~linev/xmlfile.tar.gz
//
// TXMLPlayer class generates one function per class, which called class streamer.
// Name of such function for class TExample will be TExample_streamer.
//
// Following data members for streamed classes are supported:
//  - simple data types (int, double, float)
//  - array of simple types (int[5], double[5][6])
//  - dynamic array of simple types (int* with comment field // [fSize])
//  - const char*
//  - object of any nonROOT class
//  - pointer on object
//  - array of objects
//  - array of pointers on objects
//  - stl string
//  - stl vector, list, deque, set, multiset, map, multimap
//  - allowed arguments for stl containers are: simple data types, string, object, pointer on object
//  Any other data member can not be (yet) read from xml file and write to xml file.
//
// If data member of class is private or protected, it can not be accessed via
// member name. Two alternative way is supported. First, if for class member fValue
// exists function GetValue(), it will be used to get value from the class, and if
// exists SetValue(), it will be used to set appropriate data member. Names of setter
// and getter methods can be specified in comments filed like:
//
//     int  fValue;   // *OPTION={GetMethod="GetV";SetMethod="SetV"}
//
// If getter or setter methods does not available, address to data member will be
// calculated as predefined offset to object start address. In that case generated code
// should be used only on the same platform (OS + compiler), where it was generated.
//
// Generated streamers resolve inheritance tree for given class. This allows to have
// array (or vector) of object pointers on some basic class, while objects of derived
// class(es) are used.
//
// To access data from xml files, user should use TXmlFile class, which is different from
// ROOT TXMLFile, but provides very similar functionality. For example, to read
// object from xml file:
//
//        TXmlFile file("test.xml");             // open xml file
//        file.ls();                             // show list of keys in file
//        TExample* ex1 = (TExample*) file.Get("ex1", TExample_streamer); // get object
//        file.Close();
//
// To write object to file:
//
//        TXmlFile outfile("test2.xml", "recreate");    // create xml file
//        TExample* ex1 = new TExample;
//        outfile.Write(ex1, "ex1", TExample_streamer);   // write object to file
//        outfile.Close();
//
// Complete example for generating and using of external xml streamers can be taken from
// http://www-docs.gsi.de/~linev/xmlreader.tar.gz
//
// Any bug reports and requests for additional functionality are welcome.
//
// Sergey Linev, S.Linev@gsi.de
//
//________________________________________________________________________

#include "TXMLPlayer.h"

#include "TROOT.h"
#include "TList.h"
#include "TClass.h"
#include "TVirtualStreamerInfo.h"
#include "TStreamerElement.h"
#include "TObjArray.h"
#include "TDataMember.h"
#include "TMethod.h"
#include "TDataType.h"
#include "TMethodCall.h"
#include "TVirtualCollectionProxy.h"
#include "TClassEdit.h"
#include "strlcpy.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

const char *tab1 = "   ";
const char *tab2 = "      ";
const char *tab3 = "         ";
const char *tab4 = "            ";

const char *names_xmlfileclass = "TXmlFile";

ClassImp(TXMLPlayer);

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TXMLPlayer::TXMLPlayer() : TObject()
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor of TXMLPlayer object

TXMLPlayer::~TXMLPlayer()
{
}

////////////////////////////////////////////////////////////////////////////////
/// returns streamer function name for given class

TString TXMLPlayer::GetStreamerName(TClass *cl)
{
   if (cl == 0)
      return "";
   TString res = cl->GetName();
   res += "_streamer";
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Produce streamers for provide class list
/// TList should include list of classes, for which code should be generated.
/// filename specify name of file (without extension), where streamers should be
/// created. Function produces two files: header file and source file.
/// For instance, if filename is "streamers", files "streamers.h" and "streamers.cxx"
/// will be created.

Bool_t TXMLPlayer::ProduceCode(TList *cllist, const char *filename)
{
   if ((cllist == 0) || (filename == 0))
      return kFALSE;

   std::ofstream fh(TString(filename) + ".h");
   std::ofstream fs(TString(filename) + ".cxx");

   fh << "// generated header file" << std::endl << std::endl;
   fh << "#ifndef " << filename << "_h" << std::endl;
   fh << "#define " << filename << "_h" << std::endl << std::endl;

   fh << "#include \"" << names_xmlfileclass << ".h\"" << std::endl << std::endl;

   fs << "// generated source file" << std::endl << std::endl;
   fs << "#include \"" << filename << ".h\"" << std::endl << std::endl;

   // produce appropriate include for all classes

   TObjArray inclfiles;
   TIter iter(cllist);
   TClass *cl = 0;
   while ((cl = (TClass *)iter()) != 0) {
      if (inclfiles.FindObject(cl->GetDeclFileName()) == 0) {
         fs << "#include \"" << cl->GetDeclFileName() << "\"" << std::endl;
         inclfiles.Add(new TNamed(cl->GetDeclFileName(), ""));
      }
   }
   inclfiles.Delete();

   fh << std::endl;
   fs << std::endl;

   // produce streamers declarations and implementations

   iter.Reset();

   while ((cl = (TClass *)iter()) != 0) {

      fh << "extern void* " << GetStreamerName(cl) << "(" << names_xmlfileclass
         << " &buf, void* ptr = 0, bool checktypes = true);" << std::endl
         << std::endl;

      ProduceStreamerSource(fs, cl, cllist);
   }

   fh << "#endif" << std::endl << std::endl;
   fs << std::endl << std::endl;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// returns name of simple data type for given data member

TString TXMLPlayer::GetMemberTypeName(TDataMember *member)
{
   if (member == 0)
      return "int";

   if (member->IsBasic())
      switch (member->GetDataType()->GetType()) {
      case kChar_t: return "char";
      case kShort_t: return "short";
      case kInt_t: return "int";
      case kLong_t: return "long";
      case kLong64_t: return "long long";
      case kFloat16_t:
      case kFloat_t: return "float";
      case kDouble32_t:
      case kDouble_t: return "double";
      case kUChar_t: {
         char first = member->GetDataType()->GetTypeName()[0];
         if ((first == 'B') || (first == 'b'))
            return "bool";
         return "unsigned char";
      }
      case kBool_t: return "bool";
      case kUShort_t: return "unsigned short";
      case kUInt_t: return "unsigned int";
      case kULong_t: return "unsigned long";
      case kULong64_t: return "unsigned long long";
      }

   if (member->IsEnum())
      return "int";

   return member->GetTypeName();
}

////////////////////////////////////////////////////////////////////////////////
/// return simple data types for given TStreamerElement object

TString TXMLPlayer::GetBasicTypeName(TStreamerElement *el)
{
   if (el->GetType() == TVirtualStreamerInfo::kCounter)
      return "int";

   switch (el->GetType() % 20) {
   case TVirtualStreamerInfo::kChar: return "char";
   case TVirtualStreamerInfo::kShort: return "short";
   case TVirtualStreamerInfo::kInt: return "int";
   case TVirtualStreamerInfo::kLong: return "long";
   case TVirtualStreamerInfo::kLong64: return "long long";
   case TVirtualStreamerInfo::kFloat16:
   case TVirtualStreamerInfo::kFloat: return "float";
   case TVirtualStreamerInfo::kDouble32:
   case TVirtualStreamerInfo::kDouble: return "double";
   case TVirtualStreamerInfo::kUChar: {
      char first = el->GetTypeNameBasic()[0];
      if ((first == 'B') || (first == 'b'))
         return "bool";
      return "unsigned char";
   }
   case TVirtualStreamerInfo::kBool: return "bool";
   case TVirtualStreamerInfo::kUShort: return "unsigned short";
   case TVirtualStreamerInfo::kUInt: return "unsigned int";
   case TVirtualStreamerInfo::kULong: return "unsigned long";
   case TVirtualStreamerInfo::kULong64: return "unsigned long long";
   }
   return "int";
}

////////////////////////////////////////////////////////////////////////////////
/// return functions name to read simple data type from xml file

TString TXMLPlayer::GetBasicTypeReaderMethodName(Int_t type, const char *realname)
{
   if (type == TVirtualStreamerInfo::kCounter)
      return "ReadInt";

   switch (type % 20) {
   case TVirtualStreamerInfo::kChar: return "ReadChar";
   case TVirtualStreamerInfo::kShort: return "ReadShort";
   case TVirtualStreamerInfo::kInt: return "ReadInt";
   case TVirtualStreamerInfo::kLong: return "ReadLong";
   case TVirtualStreamerInfo::kLong64: return "ReadLong64";
   case TVirtualStreamerInfo::kFloat16:
   case TVirtualStreamerInfo::kFloat: return "ReadFloat";
   case TVirtualStreamerInfo::kDouble32:
   case TVirtualStreamerInfo::kDouble: return "ReadDouble";
   case TVirtualStreamerInfo::kUChar: {
      Bool_t isbool = false;
      if (realname != 0)
         isbool = (TString(realname).Index("bool", 0, TString::kIgnoreCase) >= 0);
      if (isbool)
         return "ReadBool";
      return "ReadUChar";
   }
   case TVirtualStreamerInfo::kBool: return "ReadBool";
   case TVirtualStreamerInfo::kUShort: return "ReadUShort";
   case TVirtualStreamerInfo::kUInt: return "ReadUInt";
   case TVirtualStreamerInfo::kULong: return "ReadULong";
   case TVirtualStreamerInfo::kULong64: return "ReadULong64";
   }
   return "ReadValue";
}

////////////////////////////////////////////////////////////////////////////////
/// produce code to access member of given class.
/// Parameter specials has following meaning:
///    0 - nothing special
///    1 - cast to data type
///    2 - produce pointer on given member
///    3 - skip casting when produce pointer by buf.P() function

const char *TXMLPlayer::ElementGetter(TClass *cl, const char *membername, int specials)
{
   TClass *membercl = cl ? cl->GetBaseDataMember(membername) : 0;
   TDataMember *member = membercl ? membercl->GetDataMember(membername) : 0;
   TMethodCall *mgetter = member ? member->GetterMethod(0) : 0;

   if ((mgetter != 0) && (mgetter->GetMethod()->Property() & kIsPublic)) {
      fGetterName = "obj->";
      fGetterName += mgetter->GetMethodName();
      fGetterName += "()";
   } else if ((member == 0) || ((member->Property() & kIsPublic) != 0)) {
      fGetterName = "obj->";
      fGetterName += membername;
   } else {
      fGetterName = "";
      Bool_t deref = (member->GetArrayDim() == 0) && (specials != 2);
      if (deref)
         fGetterName += "*(";
      if (specials != 3) {
         fGetterName += "(";
         if (member->Property() & kIsConstant)
            fGetterName += "const ";
         fGetterName += GetMemberTypeName(member);
         if (member->IsaPointer())
            fGetterName += "*";
         fGetterName += "*) ";
      }
      fGetterName += "buf.P(obj,";
      fGetterName += member->GetOffset();
      fGetterName += ")";
      if (deref)
         fGetterName += ")";
      specials = 0;
   }

   if ((specials == 1) && (member != 0)) {
      TString cast = "(";
      cast += GetMemberTypeName(member);
      if (member->IsaPointer() || (member->GetArrayDim() > 0))
         cast += "*";
      cast += ") ";
      cast += fGetterName;
      fGetterName = cast;
   }

   if ((specials == 2) && (member != 0)) {
      TString buf = "&(";
      buf += fGetterName;
      buf += ")";
      fGetterName = buf;
   }

   return fGetterName.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Produce code to set value to given data member.
/// endch should be output after value is specified.

const char *TXMLPlayer::ElementSetter(TClass *cl, const char *membername, char *endch)
{
   strcpy(endch, "");

   TClass *membercl = cl ? cl->GetBaseDataMember(membername) : 0;
   TDataMember *member = membercl ? membercl->GetDataMember(membername) : 0;
   TMethodCall *msetter = member ? member->SetterMethod(cl) : 0;

   if ((msetter != 0) && (msetter->GetMethod()->Property() & kIsPublic)) {
      fSetterName = "obj->";
      fSetterName += msetter->GetMethodName();
      fSetterName += "(";
      strcpy(endch, ")");
   } else if ((member == 0) || (member->Property() & kIsPublic) != 0) {
      fSetterName = "obj->";
      fSetterName += membername;
      fSetterName += " = ";
   } else {
      fSetterName = "";
      if (member->GetArrayDim() == 0)
         fSetterName += "*";
      fSetterName += "((";
      if (member->Property() & kIsConstant)
         fSetterName += "const ";
      fSetterName += GetMemberTypeName(member);
      if (member->IsaPointer())
         fSetterName += "*";
      fSetterName += "*) buf.P(obj,";
      fSetterName += member->GetOffset();
      fSetterName += ")) = ";
   }
   return fSetterName.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Produce source code of streamer function for specified class

void TXMLPlayer::ProduceStreamerSource(std::ostream &fs, TClass *cl, TList *cllist)
{
   if (cl == 0)
      return;
   TVirtualStreamerInfo *info = cl->GetStreamerInfo();
   TObjArray *elements = info->GetElements();
   if (elements == 0)
      return;

   fs << "//__________________________________________________________________________" << std::endl;
   fs << "void* " << GetStreamerName(cl) << "(" << names_xmlfileclass << " &buf, void* ptr, bool checktypes)"
      << std::endl;
   fs << "{" << std::endl;
   fs << tab1 << cl->GetName() << " *obj = (" << cl->GetName() << "*) ptr;" << std::endl;

   fs << tab1 << "if (buf.IsReading()) { " << std::endl;

   TIter iter(cllist);
   TClass *c1 = 0;
   Bool_t firstchild = true;

   while ((c1 = (TClass *)iter()) != 0) {
      if (c1 == cl)
         continue;
      if (c1->GetListOfBases()->FindObject(cl->GetName()) == 0)
         continue;
      if (firstchild) {
         fs << tab2 << "if (checktypes) {" << std::endl;
         fs << tab3 << "void* ";
         firstchild = false;
      } else
         fs << tab3;
      fs << "res = " << GetStreamerName(c1) << "(buf, dynamic_cast<" << c1->GetName() << "*>(obj));" << std::endl;
      fs << tab3 << "if (res) return dynamic_cast<" << cl->GetName() << "*>((" << c1->GetName() << " *) res);"
         << std::endl;
   }
   if (!firstchild)
      fs << tab2 << "}" << std::endl;

   fs << tab2 << "if (!buf.CheckClassNode(\"" << cl->GetName() << "\", " << info->GetClassVersion() << ")) return 0;"
      << std::endl;

   fs << tab2 << "if (obj==0) obj = new " << cl->GetName() << ";" << std::endl;

   int n;
   for (n = 0; n <= elements->GetLast(); n++) {

      TStreamerElement *el = dynamic_cast<TStreamerElement *>(elements->At(n));
      if (el == 0)
         continue;

      Int_t typ = el->GetType();

      switch (typ) {
      // basic types
      case TVirtualStreamerInfo::kBool:
      case TVirtualStreamerInfo::kChar:
      case TVirtualStreamerInfo::kShort:
      case TVirtualStreamerInfo::kInt:
      case TVirtualStreamerInfo::kLong:
      case TVirtualStreamerInfo::kLong64:
      case TVirtualStreamerInfo::kFloat:
      case TVirtualStreamerInfo::kFloat16:
      case TVirtualStreamerInfo::kDouble:
      case TVirtualStreamerInfo::kUChar:
      case TVirtualStreamerInfo::kUShort:
      case TVirtualStreamerInfo::kUInt:
      case TVirtualStreamerInfo::kULong:
      case TVirtualStreamerInfo::kULong64:
      case TVirtualStreamerInfo::kDouble32:
      case TVirtualStreamerInfo::kCounter: {
         char endch[5];
         fs << tab2 << ElementSetter(cl, el->GetName(), endch);
         fs << "buf." << GetBasicTypeReaderMethodName(el->GetType(), 0) << "(\"" << el->GetName() << "\")" << endch
            << ";" << std::endl;
         continue;
      }

      // array of basic types like bool[10]
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kBool:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kChar:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kShort:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kInt:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kLong:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kLong64:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kFloat:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kFloat16:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kDouble:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUChar:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUShort:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUInt:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kULong:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kULong64:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kDouble32: {
         fs << tab2 << "buf.ReadArray(" << ElementGetter(cl, el->GetName(), (el->GetArrayDim() > 1) ? 1 : 0);
         fs << ", " << el->GetArrayLength() << ", \"" << el->GetName() << "\");" << std::endl;
         continue;
      }

      // array of basic types like bool[n]
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kBool:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kChar:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kShort:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kInt:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kLong:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kLong64:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kFloat:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kFloat16:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kDouble:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUChar:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUShort:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUInt:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kULong:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kULong64:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kDouble32: {
         TStreamerBasicPointer *elp = dynamic_cast<TStreamerBasicPointer *>(el);
         if (elp == 0) {
            std::cout << "fatal error with TStreamerBasicPointer" << std::endl;
            continue;
         }
         char endch[5];

         fs << tab2 << ElementSetter(cl, el->GetName(), endch);
         fs << "buf.ReadArray(" << ElementGetter(cl, el->GetName());
         fs << ", " << ElementGetter(cl, elp->GetCountName());
         fs << ", \"" << el->GetName() << "\", true)" << endch << ";" << std::endl;
         continue;
      }

      case TVirtualStreamerInfo::kCharStar: {
         char endch[5];
         fs << tab2 << ElementSetter(cl, el->GetName(), endch);
         fs << "buf.ReadCharStar(" << ElementGetter(cl, el->GetName());
         fs << ", \"" << el->GetName() << "\")" << endch << ";" << std::endl;
         continue;
      }

      case TVirtualStreamerInfo::kBase: {
         fs << tab2 << GetStreamerName(el->GetClassPointer()) << "(buf, dynamic_cast<"
            << el->GetClassPointer()->GetName() << "*>(obj), false);" << std::endl;
         continue;
      }

      // Class*   Class not derived from TObject and with comment field //->
      case TVirtualStreamerInfo::kAnyp:
      case TVirtualStreamerInfo::kAnyp + TVirtualStreamerInfo::kOffsetL: {
         if (el->GetArrayLength() > 0) {
            fs << tab2 << "buf.ReadObjectArr(" << ElementGetter(cl, el->GetName());
            fs << ", " << el->GetArrayLength() << ", -1"
               << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         } else {
            fs << tab2 << "buf.ReadObject(" << ElementGetter(cl, el->GetName());
            fs << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         }
         continue;
      }

      // Class*   Class not derived from TObject and no comment
      case TVirtualStreamerInfo::kAnyP:
      case TVirtualStreamerInfo::kAnyP + TVirtualStreamerInfo::kOffsetL: {
         if (el->GetArrayLength() > 0) {
            fs << tab2 << "for (int n=0;n<" << el->GetArrayLength() << ";n++) "
               << "delete (" << ElementGetter(cl, el->GetName()) << ")[n];" << std::endl;
            fs << tab2 << "buf.ReadObjectPtrArr((void**) " << ElementGetter(cl, el->GetName(), 3);
            fs << ", " << el->GetArrayLength() << ", \"" << el->GetName() << "\", "
               << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         } else {
            char endch[5];

            fs << tab2 << "delete " << ElementGetter(cl, el->GetName()) << ";" << std::endl;
            fs << tab2 << ElementSetter(cl, el->GetName(), endch);
            fs << "(" << el->GetClassPointer()->GetName() << "*) buf.ReadObjectPtr(\"" << el->GetName() << "\", "
               << GetStreamerName(el->GetClassPointer()) << ")" << endch << ";" << std::endl;
         }
         continue;
      }

      // Class  NOT derived from TObject
      case TVirtualStreamerInfo::kAny: {
         fs << tab2 << "buf.ReadObject(" << ElementGetter(cl, el->GetName(), 2);
         fs << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         continue;
      }

      // Class  NOT derived from TObject, array
      case TVirtualStreamerInfo::kAny + TVirtualStreamerInfo::kOffsetL: {
         fs << tab2 << "buf.ReadObjectArr(" << ElementGetter(cl, el->GetName());
         fs << ", " << el->GetArrayLength() << ", sizeof(" << el->GetClassPointer()->GetName() << "), \""
            << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         continue;
      }

      // container with no virtual table (stl) and no comment
      case TVirtualStreamerInfo::kSTLp:
      case TVirtualStreamerInfo::kSTL:
      case TVirtualStreamerInfo::kSTLp + TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kSTL + TVirtualStreamerInfo::kOffsetL: {
         TStreamerSTL *elstl = dynamic_cast<TStreamerSTL *>(el);
         if (elstl == 0)
            break; // to make skip

         if (ProduceSTLstreamer(fs, cl, elstl, false))
            continue;

         fs << tab2 << "// STL type = " << elstl->GetSTLtype() << std::endl;
         break;
      }
      }
      fs << tab2 << "buf.SkipMember(\"" << el->GetName() << "\");   // sinfo type " << el->GetType() << " of class "
         << el->GetClassPointer()->GetName() << " not supported" << std::endl;
   }

   fs << tab2 << "buf.EndClassNode();" << std::endl;

   fs << tab1 << "} else {" << std::endl;

   // generation of writing part of class streamer

   fs << tab2 << "if (obj==0) return 0;" << std::endl;

   firstchild = true;
   iter.Reset();
   while ((c1 = (TClass *)iter()) != 0) {
      if (c1 == cl)
         continue;
      if (c1->GetListOfBases()->FindObject(cl->GetName()) == 0)
         continue;
      if (firstchild) {
         firstchild = false;
         fs << tab2 << "if (checktypes) {" << std::endl;
      }
      fs << tab3 << "if (dynamic_cast<" << c1->GetName() << "*>(obj))" << std::endl;
      fs << tab4 << "return " << GetStreamerName(c1) << "(buf, dynamic_cast<" << c1->GetName() << "*>(obj));"
         << std::endl;
   }
   if (!firstchild)
      fs << tab2 << "}" << std::endl;

   fs << tab2 << "buf.StartClassNode(\"" << cl->GetName() << "\", " << info->GetClassVersion() << ");" << std::endl;

   for (n = 0; n <= elements->GetLast(); n++) {

      TStreamerElement *el = dynamic_cast<TStreamerElement *>(elements->At(n));
      if (el == 0)
         continue;

      Int_t typ = el->GetType();

      switch (typ) {
      // write basic types
      case TVirtualStreamerInfo::kBool:
      case TVirtualStreamerInfo::kChar:
      case TVirtualStreamerInfo::kShort:
      case TVirtualStreamerInfo::kInt:
      case TVirtualStreamerInfo::kLong:
      case TVirtualStreamerInfo::kLong64:
      case TVirtualStreamerInfo::kFloat:
      case TVirtualStreamerInfo::kFloat16:
      case TVirtualStreamerInfo::kDouble:
      case TVirtualStreamerInfo::kUChar:
      case TVirtualStreamerInfo::kUShort:
      case TVirtualStreamerInfo::kUInt:
      case TVirtualStreamerInfo::kULong:
      case TVirtualStreamerInfo::kULong64:
      case TVirtualStreamerInfo::kDouble32:
      case TVirtualStreamerInfo::kCounter: {
         fs << tab2 << "buf.WriteValue(";
         if (typ == TVirtualStreamerInfo::kUChar)
            fs << "(unsigned char) " << ElementGetter(cl, el->GetName());
         else
            fs << ElementGetter(cl, el->GetName());
         fs << ", \"" << el->GetName() << "\");" << std::endl;
         continue;
      }

      // array of basic types
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kBool:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kChar:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kShort:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kInt:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kLong:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kLong64:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kFloat:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kFloat16:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kDouble:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUChar:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUShort:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUInt:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kULong:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kULong64:
      case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kDouble32: {
         fs << tab2 << "buf.WriteArray(" << ElementGetter(cl, el->GetName(), (el->GetArrayDim() > 1) ? 1 : 0);
         fs << ", " << el->GetArrayLength() << ", \"" << el->GetName() << "\");" << std::endl;
         continue;
      }

      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kBool:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kChar:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kShort:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kInt:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kLong:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kLong64:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kFloat:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kFloat16:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kDouble:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUChar:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUShort:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUInt:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kULong:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kULong64:
      case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kDouble32: {
         TStreamerBasicPointer *elp = dynamic_cast<TStreamerBasicPointer *>(el);
         if (elp == 0) {
            std::cout << "fatal error with TStreamerBasicPointer" << std::endl;
            continue;
         }
         fs << tab2 << "buf.WriteArray(" << ElementGetter(cl, el->GetName());
         fs << ", " << ElementGetter(cl, elp->GetCountName()) << ", \"" << el->GetName() << "\", true);" << std::endl;
         continue;
      }

      case TVirtualStreamerInfo::kCharStar: {
         fs << tab2 << "buf.WriteCharStar(" << ElementGetter(cl, el->GetName()) << ", \"" << el->GetName() << "\");"
            << std::endl;
         continue;
      }

      case TVirtualStreamerInfo::kBase: {
         fs << tab2 << GetStreamerName(el->GetClassPointer()) << "(buf, dynamic_cast<"
            << el->GetClassPointer()->GetName() << "*>(obj), false);" << std::endl;
         continue;
      }

      // Class*   Class not derived from TObject and with comment field //->
      case TVirtualStreamerInfo::kAnyp:
      case TVirtualStreamerInfo::kAnyp + TVirtualStreamerInfo::kOffsetL: {
         if (el->GetArrayLength() > 0) {
            fs << tab2 << "buf.WriteObjectArr(" << ElementGetter(cl, el->GetName());
            fs << ", " << el->GetArrayLength() << ", -1"
               << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         } else {
            fs << tab2 << "buf.WriteObject(" << ElementGetter(cl, el->GetName());
            fs << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         }
         continue;
      }

      // Class*   Class not derived from TObject and no comment
      case TVirtualStreamerInfo::kAnyP:
      case TVirtualStreamerInfo::kAnyP + TVirtualStreamerInfo::kOffsetL: {
         if (el->GetArrayLength() > 0) {
            fs << tab2 << "buf.WriteObjectPtrArr((void**) " << ElementGetter(cl, el->GetName(), 3);
            fs << ", " << el->GetArrayLength() << ", \"" << el->GetName() << "\", "
               << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         } else {
            fs << tab2 << "buf.WriteObjectPtr(" << ElementGetter(cl, el->GetName());
            fs << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         }
         continue;
      }

      case TVirtualStreamerInfo::kAny: { // Class  NOT derived from TObject
         fs << tab2 << "buf.WriteObject(" << ElementGetter(cl, el->GetName(), 2);
         fs << ", \"" << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         continue;
      }

      case TVirtualStreamerInfo::kAny + TVirtualStreamerInfo::kOffsetL: {
         fs << tab2 << "buf.WriteObjectArr(" << ElementGetter(cl, el->GetName());
         fs << ", " << el->GetArrayLength() << ", sizeof(" << el->GetClassPointer()->GetName() << "), \""
            << el->GetName() << "\", " << GetStreamerName(el->GetClassPointer()) << ");" << std::endl;
         continue;
      }

      // container with no virtual table (stl) and no comment
      case TVirtualStreamerInfo::kSTLp + TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kSTL + TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kSTLp:
      case TVirtualStreamerInfo::kSTL: {
         TStreamerSTL *elstl = dynamic_cast<TStreamerSTL *>(el);
         if (elstl == 0)
            break; // to make skip

         if (ProduceSTLstreamer(fs, cl, elstl, true))
            continue;
         fs << tab2 << "// STL type = " << elstl->GetSTLtype() << std::endl;
         break;
      }
      }
      fs << tab2 << "buf.MakeEmptyMember(\"" << el->GetName() << "\");   // sinfo type " << el->GetType()
         << " of class " << el->GetClassPointer()->GetName() << " not supported" << std::endl;
   }

   fs << tab2 << "buf.EndClassNode();" << std::endl;

   fs << tab1 << "}" << std::endl;
   fs << tab1 << "return obj;" << std::endl;
   fs << "}" << std::endl << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Produce code to read argument of stl container from xml file

void TXMLPlayer::ReadSTLarg(std::ostream &fs, TString &argname, int argtyp, Bool_t isargptr, TClass *argcl,
                            TString &tname, TString &ifcond)
{
   switch (argtyp) {
   case TVirtualStreamerInfo::kBool:
   case TVirtualStreamerInfo::kChar:
   case TVirtualStreamerInfo::kShort:
   case TVirtualStreamerInfo::kInt:
   case TVirtualStreamerInfo::kLong:
   case TVirtualStreamerInfo::kLong64:
   case TVirtualStreamerInfo::kFloat:
   case TVirtualStreamerInfo::kFloat16:
   case TVirtualStreamerInfo::kDouble:
   case TVirtualStreamerInfo::kUChar:
   case TVirtualStreamerInfo::kUShort:
   case TVirtualStreamerInfo::kUInt:
   case TVirtualStreamerInfo::kULong:
   case TVirtualStreamerInfo::kULong64:
   case TVirtualStreamerInfo::kDouble32:
   case TVirtualStreamerInfo::kCounter: {
      fs << tname << " " << argname << " = buf." << GetBasicTypeReaderMethodName(argtyp, tname.Data()) << "(0);"
         << std::endl;
      break;
   }

   case TVirtualStreamerInfo::kObject: {
      fs << tname << (isargptr ? " " : " *") << argname << " = "
         << "(" << argcl->GetName() << "*)"
         << "buf.ReadObjectPtr(0, " << GetStreamerName(argcl) << ");" << std::endl;
      if (!isargptr) {
         if (ifcond.Length() > 0)
            ifcond += " && ";
         ifcond += argname;
         TString buf = "*";
         buf += argname;
         argname = buf;
      }
      break;
   }

   case TVirtualStreamerInfo::kSTLstring: {
      fs << "string *" << argname << " = "
         << "buf.ReadSTLstring();" << std::endl;
      if (!isargptr) {
         if (ifcond.Length() > 0)
            ifcond += " && ";
         ifcond += argname;
         TString buf = "*";
         buf += argname;
         argname = buf;
      }
      break;
   }

   default: fs << "/* argument " << argname << " not supported */";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Produce code to write argument of stl container to xml file

void TXMLPlayer::WriteSTLarg(std::ostream &fs, const char *accname, int argtyp, Bool_t isargptr, TClass *argcl)
{
   switch (argtyp) {
   case TVirtualStreamerInfo::kBool:
   case TVirtualStreamerInfo::kChar:
   case TVirtualStreamerInfo::kShort:
   case TVirtualStreamerInfo::kInt:
   case TVirtualStreamerInfo::kLong:
   case TVirtualStreamerInfo::kLong64:
   case TVirtualStreamerInfo::kFloat:
   case TVirtualStreamerInfo::kFloat16:
   case TVirtualStreamerInfo::kDouble:
   case TVirtualStreamerInfo::kUChar:
   case TVirtualStreamerInfo::kUShort:
   case TVirtualStreamerInfo::kUInt:
   case TVirtualStreamerInfo::kULong:
   case TVirtualStreamerInfo::kULong64:
   case TVirtualStreamerInfo::kDouble32:
   case TVirtualStreamerInfo::kCounter: {
      fs << "buf.WriteValue(" << accname << ", 0);" << std::endl;
      break;
   }

   case TVirtualStreamerInfo::kObject: {
      fs << "buf.WriteObjectPtr(";
      if (isargptr)
         fs << accname;
      else
         fs << "&(" << accname << ")";
      fs << ", 0, " << GetStreamerName(argcl) << ");" << std::endl;
      break;
   }

   case TVirtualStreamerInfo::kSTLstring: {
      fs << "buf.WriteSTLstring(";
      if (isargptr)
         fs << accname;
      else
         fs << "&(" << accname << ")";
      fs << ");" << std::endl;
      break;
   }

   default: fs << "/* argument not supported */" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Produce code of xml streamer for data member of stl type

Bool_t TXMLPlayer::ProduceSTLstreamer(std::ostream &fs, TClass *cl, TStreamerSTL *el, Bool_t isWriting)
{
   if ((cl == 0) || (el == 0))
      return false;

   TClass *contcl = el->GetClassPointer();

   Bool_t isstr = (el->GetSTLtype() == ROOT::kSTLstring);
   Bool_t isptr = el->IsaPointer();
   Bool_t isarr = (el->GetArrayLength() > 0);
   Bool_t isparent = (strcmp(el->GetName(), contcl->GetName()) == 0);

   int stltyp = -1;
   int narg = 0;
   int argtype[2];
   Bool_t isargptr[2];
   TClass *argcl[2];
   TString argtname[2];

   if (!isstr && contcl->GetCollectionType() != ROOT::kNotSTL) {
      int nestedLoc = 0;
      std::vector<std::string> splitName;
      TClassEdit::GetSplit(contcl->GetName(), splitName, nestedLoc);

      stltyp = contcl->GetCollectionType();
      switch (stltyp) {
      case ROOT::kSTLvector: narg = 1; break;
      case ROOT::kSTLlist: narg = 1; break;
      case ROOT::kSTLforwardlist: narg = 1; break;
      case ROOT::kSTLdeque: narg = 1; break;
      case ROOT::kSTLmap: narg = 2; break;
      case ROOT::kSTLmultimap: narg = 2; break;
      case ROOT::kSTLset: narg = 1; break;
      case ROOT::kSTLmultiset: narg = 1; break;
      case ROOT::kSTLunorderedset: narg = 1; break;
      case ROOT::kSTLunorderedmultiset: narg = 1; break;
      case ROOT::kSTLunorderedmap: narg = 2; break;
      case ROOT::kSTLunorderedmultimap: narg = 2; break;

      default: return false;
      }

      for (int n = 0; n < narg; n++) {
         argtype[n] = -1;
         isargptr[n] = false;
         argcl[n] = 0;
         argtname[n] = "";

         TString buf = splitName[n + 1];

         argtname[n] = buf;

         // nested STL containers not yet supported
         if (TClassEdit::IsSTLCont(buf.Data()))
            return false;

         int pstar = buf.Index("*");

         if (pstar > 0) {
            isargptr[n] = true;
            pstar--;
            while ((pstar > 0) && (buf[pstar] == ' '))
               pstar--;
            buf.Remove(pstar + 1);
         } else
            isargptr[n] = false;

         if (buf.Index("const ") == 0) {
            buf.Remove(0, 6);
            while ((buf.Length() > 0) && (buf[0] == ' '))
               buf.Remove(0, 1);
         }

         TDataType *dt = (TDataType *)gROOT->GetListOfTypes()->FindObject(buf);
         if (dt)
            argtype[n] = dt->GetType();
         else if (buf == "string")
            argtype[n] = TVirtualStreamerInfo::kSTLstring;
         else {
            argcl[n] = TClass::GetClass(buf);
            if (argcl[n] != 0)
               argtype[n] = TVirtualStreamerInfo::kObject;
         }
         if (argtype[n] < 0)
            stltyp = -1;
      } // for narg

      if (stltyp < 0)
         return false;
   }

   Bool_t akaarrayaccess = (narg == 1) && (argtype[0] < 20);

   char tabs[30], tabs2[30];

   if (isWriting) {

      fs << tab2 << "if (buf.StartSTLnode(\"" << fXmlSetup.XmlGetElementName(el) << "\")) {" << std::endl;

      fs << tab3 << contcl->GetName() << " ";

      TString accname;
      if (isptr) {
         if (isarr) {
            fs << "**cont";
            accname = "(*cont)->";
         } else {
            fs << "*cont";
            accname = "cont->";
         }
      } else if (isarr) {
         fs << "*cont";
         accname = "cont->";
      } else {
         fs << "&cont";
         accname = "cont.";
      }

      fs << " = ";

      if (isparent)
         fs << "*dynamic_cast<" << contcl->GetName() << "*>(obj);" << std::endl;
      else
         fs << ElementGetter(cl, el->GetName()) << ";" << std::endl;

      if (isarr && el->GetArrayLength()) {
         strlcpy(tabs, tab4, sizeof(tabs));
         fs << tab3 << "for(int n=0;n<" << el->GetArrayLength() << ";n++) {" << std::endl;
      } else
         strlcpy(tabs, tab3, sizeof(tabs));

      strlcpy(tabs2, tabs, sizeof(tabs2));

      if (isptr) {
         strlcat(tabs2, tab1, sizeof(tabs2));
         fs << tabs << "if (" << (isarr ? "*cont" : "cont") << "==0) {" << std::endl;
         fs << tabs2 << "buf.WriteSTLsize(0" << (isstr ? ",true);" : ");") << std::endl;
         fs << tabs << "} else {" << std::endl;
      }

      fs << tabs2 << "buf.WriteSTLsize(" << accname << (isstr ? "length(), true);" : "size());") << std::endl;

      if (isstr) {
         fs << tabs2 << "buf.WriteSTLstringData(" << accname << "c_str());" << std::endl;
      } else {
         if (akaarrayaccess) {
            fs << tabs2 << argtname[0] << "* arr = new " << argtname[0] << "[" << accname << "size()];" << std::endl;
            fs << tabs2 << "int k = 0;" << std::endl;
         }

         fs << tabs2 << contcl->GetName() << "::const_iterator iter;" << std::endl;
         fs << tabs2 << "for (iter = " << accname << "begin(); iter != " << accname << "end(); iter++)";
         if (akaarrayaccess) {
            fs << std::endl << tabs2 << tab1 << "arr[k++] = *iter;" << std::endl;
            fs << tabs2 << "buf.WriteArray(arr, " << accname << "size(), 0, false);" << std::endl;
            fs << tabs2 << "delete[] arr;" << std::endl;
         } else if (narg == 1) {
            fs << std::endl << tabs2 << tab1;
            WriteSTLarg(fs, "*iter", argtype[0], isargptr[0], argcl[0]);
         } else if (narg == 2) {
            fs << " {" << std::endl;
            fs << tabs2 << tab1;
            WriteSTLarg(fs, "iter->first", argtype[0], isargptr[0], argcl[0]);
            fs << tabs2 << tab1;
            WriteSTLarg(fs, "iter->second", argtype[1], isargptr[1], argcl[1]);
            fs << tabs2 << "}" << std::endl;
         }
      } // if (isstr)

      if (isptr)
         fs << tabs << "}" << std::endl;

      if (isarr && el->GetArrayLength()) {
         if (isptr)
            fs << tabs << "cont++;" << std::endl;
         else
            fs << tabs << "(void*) cont = (char*) cont + sizeof(" << contcl->GetName() << ");" << std::endl;
         fs << tab3 << "}" << std::endl;
      }

      fs << tab3 << "buf.EndSTLnode();" << std::endl;
      fs << tab2 << "}" << std::endl;

   } else {

      fs << tab2 << "if (buf.VerifySTLnode(\"" << fXmlSetup.XmlGetElementName(el) << "\")) {" << std::endl;

      fs << tab3 << contcl->GetName() << " ";
      TString accname, accptr;
      if (isptr) {
         if (isarr) {
            fs << "**cont";
            accname = "(*cont)->";
            accptr = "*cont";
         } else {
            fs << "*cont";
            accname = "cont->";
            accptr = "cont";
         }
      } else if (isarr) {
         fs << "*cont";
         accname = "cont->";
      } else {
         fs << "&cont";
         accname = "cont.";
      }

      fs << " = ";

      if (isparent)
         fs << "*dynamic_cast<" << contcl->GetName() << "*>(obj);" << std::endl;
      else
         fs << ElementGetter(cl, el->GetName()) << ";" << std::endl;

      if (isarr && el->GetArrayLength()) {
         strlcpy(tabs, tab4, sizeof(tabs));
         fs << tab3 << "for(int n=0;n<" << el->GetArrayLength() << ";n++) {" << std::endl;
      } else
         strlcpy(tabs, tab3, sizeof(tabs));

      fs << tabs << "int size = buf.ReadSTLsize(" << (isstr ? "true);" : ");") << std::endl;

      if (isptr) {
         fs << tabs << "delete " << accptr << ";" << std::endl;
         fs << tabs << "if (size==0) " << accptr << " = 0;" << std::endl;
         fs << tabs << "        else " << accptr << " = new " << contcl->GetName() << ";" << std::endl;
         if (!isarr) {
            char endch[5];
            fs << tabs << ElementSetter(cl, el->GetName(), endch);
            fs << "cont" << endch << ";" << std::endl;
         }
      } else {
         fs << tabs << accname << (isstr ? "erase();" : "clear();") << std::endl;
      }

      if (isstr) {
         fs << tabs << "if (size>0) " << accname << "assign(buf.ReadSTLstringData(size));" << std::endl;
      } else {
         if (akaarrayaccess) {
            fs << tabs << argtname[0] << "* arr = new " << argtname[0] << "[size];" << std::endl;
            fs << tabs << "buf.ReadArray(arr, size, 0, false);" << std::endl;
         }

         fs << tabs << "for(int k=0;k<size;k++)";

         if (akaarrayaccess) {
            fs << std::endl << tabs << tab1 << accname;
            if ((stltyp == ROOT::kSTLset) || (stltyp == ROOT::kSTLmultiset))
               fs << "insert";
            else
               fs << "push_back";
            fs << "(arr[k]);" << std::endl;
            fs << tabs << "delete[] arr;" << std::endl;
         } else if (narg == 1) {
            TString arg1("arg"), ifcond;
            fs << " {" << std::endl << tabs << tab1;
            ReadSTLarg(fs, arg1, argtype[0], isargptr[0], argcl[0], argtname[0], ifcond);
            fs << tabs << tab1;
            if (ifcond.Length() > 0)
               fs << "if (" << ifcond << ") ";
            fs << accname;
            if ((stltyp == ROOT::kSTLset) || (stltyp == ROOT::kSTLmultiset))
               fs << "insert";
            else
               fs << "push_back";
            fs << "(" << arg1 << ");" << std::endl;
            fs << tabs << "}" << std::endl;
         } else if (narg == 2) {
            TString arg1("arg1"), arg2("arg2"), ifcond;
            fs << " {" << std::endl << tabs << tab1;
            ReadSTLarg(fs, arg1, argtype[0], isargptr[0], argcl[0], argtname[0], ifcond);
            fs << tabs << tab1;
            ReadSTLarg(fs, arg2, argtype[1], isargptr[1], argcl[1], argtname[1], ifcond);
            fs << tabs << tab1;
            if (ifcond.Length() > 0)
               fs << "if (" << ifcond << ") ";
            fs << accname << "insert(make_pair(" << arg1 << ", " << arg2 << "));" << std::endl;
            fs << tabs << "}" << std::endl;
         }
      }

      if (isarr && el->GetArrayLength()) {
         if (isptr)
            fs << tabs << "cont++;" << std::endl;
         else
            fs << tabs << "(void*) cont = (char*) cont + sizeof(" << contcl->GetName() << ");" << std::endl;
         fs << tab3 << "}" << std::endl;
      }

      fs << tab3 << "buf.EndSTLnode();" << std::endl;
      fs << tab2 << "}" << std::endl;
   }
   return true;
}
