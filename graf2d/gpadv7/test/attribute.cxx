#include "gtest/gtest.h"

#include "ROOT/RAttrText.hxx"
#include "ROOT/RAttrFill.hxx"
#include "ROOT/RAttrLine.hxx"

using namespace ROOT::Experimental;

class CustomAttrs : public RAttrBase {
   RAttrLine    fAttrLine{this, "line"};    ///<! line attributes
   RAttrFill    fAttrFill{this, "fill"};    ///<! fill attributes
   RAttrText    fAttrText{this, "text"};    ///<! text attributes

protected:
   // required here while dictionary for CustomAttrs not created
   RAttrMap CollectDefaults() const override { return RAttrMap().AddDefaults(fAttrLine).AddDefaults(fAttrFill).AddDefaults(fAttrText); }

   R__ATTR_CLASS(CustomAttrs, "custom");

   const RAttrLine &AttrLine() const { return fAttrLine; }
   RAttrLine &AttrLine() { return fAttrLine; }

   const RAttrFill &AttrFill() const { return fAttrFill; }
   RAttrFill &AttrFill() { return fAttrFill; }

   const RAttrText &AttrText() const { return fAttrText; }
   RAttrText &AttrText() { return fAttrText; }

   double GetDirect(const std::string &name) const { return GetValue<double>(name); }
};


TEST(OptsTest, AttribStrings) {
   CustomAttrs attrs;

   attrs.AttrLine().SetWidth(42.);
   attrs.AttrText().SetSize(1.7);

   {
      auto val = attrs.GetDirect("line_width");
      EXPECT_FLOAT_EQ(val, 42.f);
   }

   {
      auto val = attrs.GetDirect("text_size");
      EXPECT_FLOAT_EQ(val, 1.7f);
   }
}

TEST(OptsTest, AttribVals) {
   CustomAttrs attrs;

   attrs.AttrText().SetColor(RColor::kBlue);
   auto &line = attrs.AttrLine();
   line.SetWidth(42.);

   {
      // Value was set on this attr, not coming from style:
      EXPECT_FLOAT_EQ(attrs.AttrLine().GetWidth(), 42.f);
      EXPECT_FLOAT_EQ(line.GetWidth(), 42.f);
   }

   {
      // Value was set on this attr, not coming from style:
      EXPECT_EQ(attrs.AttrText().GetColor(), RColor::kBlue);
   }

}

TEST(OptsTest, NullAttribCompare) {
   RAttrLine al1;
   RAttrLine al2;
   EXPECT_EQ(al1, al2);
   EXPECT_EQ(al2, al1);
}

TEST(OptsTest, AttribEqual) {
   CustomAttrs attrs;

   auto &al1 = attrs.AttrLine();
   auto &al2 = attrs.AttrLine();
   EXPECT_EQ(al1, al2);
   EXPECT_EQ(al2, al1);

   al1.SetColor(RColor::kRed);

   EXPECT_EQ(al1, al2);
   EXPECT_EQ(al2, al1);
}

TEST(OptsTest, AttribDiffer) {
   CustomAttrs attrs1;
   CustomAttrs attrs2;
   CustomAttrs attrs3;

   attrs1.AttrLine().SetWidth(7.);
   EXPECT_NE(attrs1, attrs2);
   EXPECT_NE(attrs2, attrs1);
   EXPECT_EQ(attrs2, attrs3);
   EXPECT_EQ(attrs3, attrs2);

   attrs2.AttrLine().SetColor(RColor::kRed);
   EXPECT_NE(attrs1, attrs2);
   EXPECT_NE(attrs2, attrs1);
   EXPECT_NE(attrs1, attrs3);
   EXPECT_NE(attrs3, attrs1);
   EXPECT_NE(attrs2, attrs3);
   EXPECT_NE(attrs3, attrs2);
}


TEST(OptsTest, AttribAssign) {
   CustomAttrs attrs1;
   CustomAttrs attrs2;

   // deep copy - independent from origin
   auto attrLine1 = attrs1.AttrLine();
   auto attrLine2 = attrs2.AttrLine();

   EXPECT_EQ(attrLine2, attrLine1);
   EXPECT_EQ(attrLine1, attrLine2);

   attrLine1.SetWidth(42.);
   EXPECT_NE(attrLine2, attrLine1);

   attrLine2 = attrLine1;
   EXPECT_EQ(attrLine2, attrLine1);
   EXPECT_EQ(attrLine1, attrLine2);

   // But original attributes now differ
   EXPECT_NE(attrs1.AttrLine(), attrLine1);
   EXPECT_NE(attrs2.AttrLine(), attrLine2);

   EXPECT_FLOAT_EQ(attrLine1.GetWidth(), 42.);
   EXPECT_FLOAT_EQ(attrLine2.GetWidth(), 42.);
   // default width return 1
   EXPECT_FLOAT_EQ(attrs1.AttrLine().GetWidth(), 1.);
   EXPECT_FLOAT_EQ(attrs2.AttrLine().GetWidth(), 1.);

   // Are the two attributes disconnected?
   attrLine2.SetWidth(3.);
   EXPECT_EQ(attrs1.AttrLine(), attrs2.AttrLine());
   EXPECT_FLOAT_EQ(attrLine1.GetWidth(), 42.);
   EXPECT_FLOAT_EQ(attrLine2.GetWidth(), 3.);
   EXPECT_FLOAT_EQ(attrs1.AttrLine().GetWidth(), 1.);
   EXPECT_FLOAT_EQ(attrs2.AttrLine().GetWidth(), 1.);
}
