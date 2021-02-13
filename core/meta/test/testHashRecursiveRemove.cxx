#include "ROOTUnitTestSupport.h"

#include "TClass.h"
#include "TClassTable.h"
#include "TInterpreter.h"
#include <memory>

#include "gtest/gtest.h"

const char *gCode = R"CODE(

   #include "TROOT.h"
   #include <iostream>

   class FirstOverload : public TObject
   {
   public:
      virtual ULong_t     Hash() const { return 1; }

      ClassDefInline(FirstOverload, 2);
   };

   class SecondOverload : public FirstOverload // Could also have used TNamed.
   {
   public:
      virtual ULong_t     Hash() const { return 2; }

      ClassDefInline(SecondOverload, 2);
   };

   class SecondNoHash : public FirstOverload // Could also have used TNamed.
   {
   public:

      ClassDefInline(SecondNoHash, 2);
   };

   class SecondAbstract : public FirstOverload // Could also have used TNamed.
   {
   public:
      virtual int Get() = 0;

      ClassDef(SecondAbstract, 2);
   };

   class Third : public SecondAbstract
   {
   public:
      int Get() override { return 0; };

      ClassDefInlineOverride(Third, 2);
   };

   class FirstOverloadCorrect : public TObject
   {
   public:
      ~FirstOverloadCorrect() {
         ROOT::CallRecursiveRemoveIfNeeded(*this);
      }
      virtual ULong_t     Hash() const { return 3; }

      ClassDefInline(FirstOverloadCorrect, 2);
   };

   class SecondCorrectAbstract : public FirstOverloadCorrect // Could also have used TNamed.
   {
   public:
      virtual int Get() = 0;

      ClassDef(SecondCorrectAbstract, 2);
   };

   class SecondCorrectAbstractHash : public FirstOverloadCorrect // Could also have used TNamed.
   {
   public:
      ~SecondCorrectAbstractHash() {
         ROOT::CallRecursiveRemoveIfNeeded(*this);
      }

      virtual ULong_t Hash() const { return 4; }
      virtual int     Get() = 0;

      ClassDef(SecondCorrectAbstractHash, 2);
   };

   class ThirdCorrect : public SecondCorrectAbstract
   {
   public:
      int Get() override { return 0; };

      ClassDefInlineOverride(ThirdCorrect, 2);
   };

   class SecondInCorrectAbstract : public FirstOverloadCorrect // Could also have used TNamed.
   {
   public:
      virtual ULong_t Hash() const { return 5; }
      virtual int     Get() = 0;

      ClassDef(SecondInCorrectAbstract, 2);
   };

   class ThirdInCorrect : public SecondInCorrectAbstract
   {
   public:
      int Get() override { return 0; };

      ClassDefInlineOverride(ThirdInCorrect, 2);
   };

   // Just declare this one so Cling will know it, but
   // do not use it to avoid the TClass being stuck in
   // kInterpreted state.
   class WrongSetup : public TObject
   {
   public:
      virtual ULong_t     Hash() const { return 6; }

      ClassDefInline(WrongSetup, 2);
   };

   // This example is valid according to C++11, 9.2/16: In addition, if class T has a user-declared constructor (12.1),
   // every non-static data member of class T shall have a name different from T.
   //
   class Rho: public TObject
   {
      public:
      Float_t Rho; // rho energy density
      Float_t Edges[2]; // pseudorapidity range edges

      ClassDef(Rho, 1)
   };
)CODE";

class WrongSetup : public TObject {
public:
   virtual ULong_t Hash() const { return 6; }

   ClassDefInline(WrongSetup, 2);
};

class InlineCompiledOnly : public TObject {
public:
   virtual ULong_t Hash() const { return 6; }

   ClassDefInline(InlineCompiledOnly, 2);
};

std::string errormessage(const std::string &baseclass, const std::string &seenclass)
{
   std::string message =
      "The class " + baseclass +
      " overrides TObject::Hash but does not call TROOT::RecursiveRemove in its destructor (seen while checking " +
      seenclass + ").";
   return message;
}

const char *errorlocation = "ROOT::Internal::TCheckHashRecursiveRemoveConsistency::CheckRecursiveRemove";

void DeclareFailingClasses()
{
   gInterpreter->Declare(gCode);
}

TEST(HashRecursiveRemove, GetClassClassDefInline)
{
   DeclareFailingClasses();

   auto getcl = TClass::GetClass("FirstOverload");
   EXPECT_NE(nullptr, getcl);

   std::unique_ptr<TObject> obj((TObject *)getcl->New());
   EXPECT_NE(nullptr, obj.get());

   auto isacl = obj->IsA();
   EXPECT_NE(nullptr, isacl);
   EXPECT_EQ(getcl, isacl);

   getcl = TClass::GetClass("WrongSetup");
   EXPECT_NE(nullptr, getcl);
   // EXPECT_NE(nullptr,TClass::GetClass("WrongSetup"));
   WrongSetup ws;
   isacl = ws.IsA();
   EXPECT_NE(nullptr, isacl);
   EXPECT_EQ(getcl, isacl);
   getcl = TClass::GetClass("WrongSetup");
   EXPECT_EQ(getcl, isacl);

   EXPECT_NE(nullptr, TClassTable::GetDict("InlineCompiledOnly"));
   getcl = TClass::GetClass("InlineCompiledOnly");
   EXPECT_NE(nullptr, getcl);
   // EXPECT_NE(nullptr,TClass::GetClass("InlineCompiledOnly"));
   InlineCompiledOnly ico;
   isacl = ico.IsA();
   EXPECT_NE(nullptr, isacl);
   EXPECT_EQ(getcl, isacl);
   getcl = TClass::GetClass("InlineCompiledOnly");
   EXPECT_EQ(getcl, isacl);
}

TEST(HashRecursiveRemove, RootClasses)
{
   EXPECT_TRUE(TClass::GetClass("TObject")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TNamed")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TH1")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TH1F")->HasConsistentHashMember());

   EXPECT_TRUE(TClass::GetClass("TEnvRec")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TDataType")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TObjArray")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TList")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("THashList")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TClass")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TInterpreter")->HasConsistentHashMember());
   // EXPECT_TRUE(TClass::GetClass("TCling")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("TMethod")->HasConsistentHashMember());
   // EXPECT_TRUE(TClass::GetClass("ROOT::Internal::TCheckHashRecursiveRemoveConsistency")->HasConsistentHashMember());
}

TEST(HashRecursiveRemove, FailingClasses)
{
   EXPECT_NE(nullptr, TClass::GetClass("FirstOverload"));
   EXPECT_TRUE(TClass::GetClass("FirstOverloadCorrect")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("SecondCorrectAbstract")->HasConsistentHashMember());
   EXPECT_FALSE(TClass::GetClass("SecondCorrectAbstractHash")->HasConsistentHashMember());
   EXPECT_TRUE(TClass::GetClass("ThirdCorrect")->HasConsistentHashMember());
   EXPECT_FALSE(TClass::GetClass("SecondInCorrectAbstract")->HasConsistentHashMember());

   ROOT_EXPECT_ERROR(EXPECT_FALSE(TClass::GetClass("FirstOverload")->HasConsistentHashMember()), errorlocation,
                     errormessage("FirstOverload", "FirstOverload"));
   ROOT_EXPECT_ERROR(EXPECT_FALSE(TClass::GetClass("SecondOverload")->HasConsistentHashMember()), errorlocation,
                     errormessage("SecondOverload", "SecondOverload"));
   ROOT_EXPECT_ERROR(EXPECT_FALSE(TClass::GetClass("SecondNoHash")->HasConsistentHashMember()), errorlocation,
                     errormessage("FirstOverload", "SecondNoHash"));
   ROOT_EXPECT_ERROR(EXPECT_FALSE(TClass::GetClass("SecondAbstract")->HasConsistentHashMember()), errorlocation,
                     errormessage("FirstOverload", "SecondAbstract"));
   ROOT_EXPECT_ERROR(EXPECT_FALSE(TClass::GetClass("Third")->HasConsistentHashMember()), errorlocation,
                     errormessage("FirstOverload", "Third"));
   ROOT_EXPECT_ERROR(EXPECT_FALSE(TClass::GetClass("ThirdInCorrect")->HasConsistentHashMember()), errorlocation,
                     errormessage("SecondInCorrectAbstract", "ThirdInCorrect"));
   ROOT_EXPECT_ERROR(EXPECT_FALSE(WrongSetup::Class()->HasConsistentHashMember()), errorlocation,
                     errormessage("WrongSetup", "WrongSetup"));
}

bool CallCheckedHash(const char *clname)
{
   auto cl = TClass::GetClass(clname);
   EXPECT_NE(nullptr, cl);

   std::unique_ptr<TObject> obj((TObject *)cl->New());
   EXPECT_NE(nullptr, obj.get());

   obj->CheckedHash();
   return !obj->HasInconsistentHash();
}
TEST(HashRecursiveRemove, CheckedHashRootClasses)
{
   EXPECT_TRUE(CallCheckedHash("TObject"));
   EXPECT_TRUE(CallCheckedHash("TNamed"));
   EXPECT_TRUE(CallCheckedHash("TH1F"));

   EXPECT_TRUE(CallCheckedHash("TEnvRec"));
   EXPECT_TRUE(CallCheckedHash("TDataType"));
   EXPECT_TRUE(CallCheckedHash("TObjArray"));
   EXPECT_TRUE(CallCheckedHash("TList"));
   EXPECT_TRUE(CallCheckedHash("THashList"));
   EXPECT_TRUE(CallCheckedHash("TClass"));
   EXPECT_TRUE(CallCheckedHash("TMethod"));
   // EXPECT_TRUE(CallCheckedHash("ROOT::Internal::TCheckHashRecursiveRemoveConsistency"));

   TObject *h1 = (TObject *)gInterpreter->Calc("new TH1I(\"histo1\",\"histo title\", 100, -10., 10.);");
   EXPECT_FALSE(h1->HasInconsistentHash());
   delete h1;
}

TEST(HashRecursiveRemove, CheckedHashFailingClasses)
{
   // testing::internal::CaptureStderr();

   EXPECT_FALSE(CallCheckedHash("FirstOverload"));
   EXPECT_FALSE(CallCheckedHash("SecondOverload"));
   EXPECT_FALSE(CallCheckedHash("SecondNoHash"));

   EXPECT_FALSE(CallCheckedHash("Third"));
   EXPECT_TRUE(CallCheckedHash("FirstOverloadCorrect"));

   EXPECT_TRUE(CallCheckedHash("ThirdCorrect"));
   EXPECT_FALSE(CallCheckedHash("ThirdInCorrect"));
   EXPECT_FALSE(CallCheckedHash("WrongSetup"));

   // std::string output = testing::internal::GetCapturedStderr();
   // EXPECT_EQ(gErrorOutput,output);
}

#include "THashList.h"
#include "TROOT.h"

constexpr unsigned int kHowMany = 20000;

TEST(HashRecursiveRemove, SimpleDelete)
{
   THashList cont;

   for (unsigned int i = 0; i < kHowMany; ++i) {
      TNamed *n = new TNamed(TString::Format("n%u", i), TString(""));
      n->SetBit(kMustCleanup);
      cont.Add(n);
   }
   cont.Delete();

   EXPECT_EQ(0, cont.GetSize());
}

TEST(HashRecursiveRemove, SimpleDirectRemove)
{
   THashList cont;
   TList todelete;

   for (unsigned int i = 0; i < kHowMany; ++i) {
      TNamed *n = new TNamed(TString::Format("n%u", i), TString(""));
      n->SetBit(kMustCleanup);
      cont.Add(n);
      todelete.Add(n);
   }
   for (auto o : todelete) {
      cont.Remove(o);
      delete o;
   }

   todelete.Clear("nodelete");

   EXPECT_EQ(0, cont.GetSize());
}

TEST(HashRecursiveRemove, DeleteWithRecursiveRemove)
{
   THashList cont;
   TList todelete;

   for (unsigned int i = 0; i < kHowMany; ++i) {
      TNamed *n = new TNamed(TString::Format("n%u", i), TString(""));
      n->SetBit(kMustCleanup);
      cont.Add(n);
      todelete.Add(n);
   }
   gROOT->GetListOfCleanups()->Add(&cont);
   cont.SetBit(kMustCleanup);

   for (auto o : todelete)
      delete o;

   todelete.Clear("nodelete");

   EXPECT_EQ(0, cont.GetSize());
}

TEST(HashRecursiveRemove, DeleteBadHashWithRecursiveRemove)
{
   THashList cont;
   TList todelete;

   for (size_t i = 0; i < kHowMany / 4; ++i) {
      TObject *o;
      if (i % 2)
         o = (TObject *)TClass::GetClass("FirstOverload")->New();
      else
         o = new WrongSetup;
      o->SetBit(kMustCleanup);
      cont.Add(o);
      todelete.Add(o);
   }
   gROOT->GetListOfCleanups()->Add(&cont);
   cont.SetBit(kMustCleanup);

   for (auto o : todelete) {
      delete o;
   }

   EXPECT_EQ(0, cont.GetSize());

   todelete.Clear("nodelete");

   // Avoid spurrious/redundant error messages in case of failure.
   cont.Clear("nodelete");
}
