#include "TClass.h"
#include "TClassRef.h"
#include "TDataMember.h"
#include "TGlobal.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TListOfDataMembers.h"
#include "TROOT.h"

#include "gtest/gtest.h"

TEST(TClingDataMemberInfo, UsingDecls)
{
   gInterpreter->Declare(R"CODE(
namespace UsingDecls {

struct EmptyStruct {
};

class Base {
protected:
  int protectedVar;
  void protectedFunc();
  void protectedFunc(float);

public:
  int func();
  float publicVar;
  static double publicStatic;
  void publicFunc();
};

class Derived: public Base {
private:
  using Base::func;

protected:
  using Base::publicStatic;

public:
  using Base::protectedVar;
  using Base::publicVar;
  using Base::protectedFunc;
};

class Rederived: public Derived {
  using Derived::publicStatic;
};

struct BagOfThings {
  union {
    int a;
    float b;
  };
  struct {
    short var;
  };
  enum {
    kEnumConst1 = 12,
    kEnumConst2 = 12
  };

  struct Nested {
    int doNotFindMe;
  };
};

struct DerivedBag: BagOfThings {
  using BagOfThings::a;
};

namespace EmptyNamespace {}

namespace SingleVarNamespace {
  int SingleVar;
}
namespace SingleUsingDeclNamespace {
  using SingleVarNamespace::SingleVar;
}
namespace WithInlineNamespace {
  inline namespace InlineNamespace {
    float InlineNamespaceVar;
    inline namespace AnotherInline {
      double InlineInlineNamespaceVar;
    }
    namespace Named {
      float InlineNamedNamespaceVar;
    }
    namespace {
      constexpr char InInlineAnon = 'A';
    }
  }
  namespace {
    Float16_t InAnon;
  }
}

namespace Split {
}
namespace Split {
  int One;
}
namespace Split {
  short Two;
}
namespace Split {
}
namespace Split {
  int Three;
}
namespace Split {
  extern short Two;
}
class Outer {
  struct {
    int InnerPrivate;
  };
public:
  struct {
    long InnerPublic;
  };
};
template <class T>
struct Tmplt {
protected:
  T fMember;
};
}
)CODE");

   struct VarInfo_t {
      const char *fName;
      const char *fType;
      long fPropertyToTest;
   };

   auto checkClass = [](const char *className, std::vector<VarInfo_t> varInfos) {
      TClass *cl = TClass::GetClass(className);
      ASSERT_NE(cl, nullptr) << "Cannot find class " << className;
      TListOfDataMembers *vars = (TListOfDataMembers *)cl->GetListOfDataMembers();
      TList allVars;
      allVars.AddAll(vars);
      allVars.AddAll(cl->GetListOfUsingDataMembers());

      EXPECT_EQ(allVars.GetSize(), varInfos.size()) << "for scope " << className << (allVars.ls(), "") << '\n';
      for (auto &&varInfo : varInfos) {
         TDataMember *dm = (TDataMember *)allVars.FindObject(varInfo.fName);
         ASSERT_NE(dm, nullptr);
         EXPECT_STREQ(dm->GetName(), varInfo.fName) << " for var " << varInfo.fName << " in class " << className;
         if (varInfo.fType) {
            EXPECT_STREQ(dm->GetTrueTypeName(), varInfo.fType)
               << " for var " << varInfo.fName << " in class " << className;
         }
         EXPECT_EQ(dm->Property() & varInfo.fPropertyToTest, varInfo.fPropertyToTest)
            << " for var " << varInfo.fName << " in class " << className;
      }
   };

   checkClass("UsingDecls::EmptyStruct", {});

   checkClass("UsingDecls::Base", {{"protectedVar", "int", kIsProtected},
                                   {"publicVar", "float", kIsPublic},
                                   {"publicStatic", "double", kIsPublic | kIsStatic}});

   checkClass("UsingDecls::Derived", {{"publicStatic", "double", kIsProtected | kIsStatic},
                                      {"protectedVar", "int", kIsPublic},
                                      {"publicVar", "float", kIsPublic}});

   checkClass("UsingDecls::Rederived", {{"publicStatic", "double", kIsPrivate}});

   checkClass("UsingDecls::BagOfThings", {{"a", "int", kIsPublic | kIsUnionMember},
                                          {"b", "float", kIsPublic | kIsUnionMember},
                                          {"var", "short", kIsPublic},
                                          {"kEnumConst1", nullptr, kIsEnum | kIsPublic},
                                          {"kEnumConst2", nullptr, kIsEnum | kIsPublic}});

   checkClass("UsingDecls::EmptyNamespace", {});

   checkClass("UsingDecls::SingleVarNamespace", {{"SingleVar", "int", kIsPublic}});

   checkClass("UsingDecls::SingleUsingDeclNamespace", {{"SingleVar", "int", kIsPublic}});

   checkClass("UsingDecls::WithInlineNamespace",
              {{"InlineNamespaceVar", "float", kIsPublic | kIsStatic},
               {"InlineInlineNamespaceVar", "double",
                kIsPublic | kIsStatic}, // LIMITATION: overrides earlier InlineNamespace::InlineInlineNamespaceVar
               {"InInlineAnon", "const char", kIsConstexpr},
               {"InAnon", "Float16_t", 0}});

   checkClass("UsingDecls::Split", {{"One", "int", 0}, {"Two", "short", 0}, {"Three", "int", 0}});

   checkClass("UsingDecls::Outer", {{"InnerPrivate", "int", kIsPrivate}, {"InnerPublic", "long", kIsPublic}});

   checkClass("UsingDecls::Tmplt<Double32_t>", {{"fMember", "Double32_t", kIsProtected}});
}


TEST(TClingDataMemberInfo, Lookup)
{
   gInterpreter->Declare(R"CODE(
namespace DMLookup {
class Outer {
  struct {
    TObject InnerPrivate;
  };
protected:
  struct {
    TNamed InnerProtected;
  };
};
})CODE");
   ASSERT_TRUE(TClass::GetClass("DMLookup::Outer"));
   ASSERT_TRUE(TClass::GetClass("DMLookup::Outer")->GetListOfDataMembers());

   ASSERT_TRUE(TClass::GetClass("DMLookup::Outer")->GetListOfDataMembers()->FindObject("InnerPrivate"));
   auto *dmInnerPrivate = (TDataMember*)TClass::GetClass("DMLookup::Outer")->GetListOfDataMembers()->FindObject("InnerPrivate");
   EXPECT_EQ(dmInnerPrivate->Property(), kIsPrivate | kIsClass);

   ASSERT_TRUE(TClass::GetClass("DMLookup::Outer")->GetListOfDataMembers()->FindObject("InnerProtected"));
   auto *dmInnerProtected = (TDataMember*)TClass::GetClass("DMLookup::Outer")->GetListOfDataMembers()->FindObject("InnerProtected");
   EXPECT_EQ(dmInnerProtected->Property(), kIsProtected | kIsClass);
}
