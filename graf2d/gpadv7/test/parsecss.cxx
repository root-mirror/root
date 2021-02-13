// @(#)root/graf2d:$Id$
// Author: Sergey Linev <s.linev@gsi.de>, 2020-02-17


#include "gtest/gtest.h"

#include <ROOT/RStyle.hxx>

TEST(RStyleTest, ParseCss)
{
   auto style = std::make_shared<ROOT::Experimental::RStyle>();

   bool res = style->ParseString(
         "/* Basic example of CSS file \n"
         "   shows how different blocks can be defined */\n"
         "csstype {\n"
         "  attr1: true;\n"
         "  attr2: 10;\n"
         "  attr3: 20.5;\n"
         "  attr4: pinc;\n"
         "}\n"
         "\n"
         "// this is rule for css class\n"
         ".cssclass {\n"
         "   attr3: value3;\n"
         "   attr4: value4;\n"
         "}\n"
         "\n"
         "// this is rule for some specific identifier\n"
         "#cssid {\n"
         "   attr5: value5;\n"
         "   attr6: value6;\n"
         "}\n"
         );

   ASSERT_EQ(res, true);

   auto field1 = style->Eval("attr1", "csstype");
   ASSERT_NE(field1, nullptr);
   EXPECT_EQ(field1->GetBool(), true);

   auto field2 = style->Eval("attr2", "csstype");
   ASSERT_NE(field2, nullptr);
   EXPECT_EQ(field2->GetInt(), 10);

   auto field3 = style->Eval("attr3", "csstype");
   ASSERT_NE(field3, nullptr);
   EXPECT_DOUBLE_EQ(field3->GetDouble(), 20.5);

   auto field4 = style->Eval("attr4", "csstype");
   ASSERT_NE(field4, nullptr);
   EXPECT_EQ(field4->GetString(), "pinc");

   field4 = style->Eval("attr4", ".cssclass");
   ASSERT_NE(field4, nullptr);
   EXPECT_EQ(field4->GetString(), "value4");

   auto field5 = style->Eval("attr5", "#cssid");

   ASSERT_NE(field5, nullptr);
   EXPECT_EQ(field5->GetString(), "value5");
}


