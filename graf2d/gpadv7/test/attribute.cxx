#include "gtest/gtest.h"

#include "ROOT/RAttrText.hxx"
#include "ROOT/RAttrBox.hxx"


using namespace ROOT::Experimental;

class CustomAttrs : public RAttrBase {
   RAttrBox fAttrBox{this, "box"};
   RAttrText fAttrText{this, "text"};

protected:
   RAttrMap CollectDefaults() const override { return RAttrMap().AddDefaults(fAttrBox).AddDefaults(fAttrText); }

   R__ATTR_CLASS(CustomAttrs, "custom");

   const RAttrBox &GetAttrBox() const { return fAttrBox; }
   CustomAttrs &SetAttrBox(const RAttrBox &box) { fAttrBox = box; return *this; }
   RAttrBox &AttrBox() { return fAttrBox; }

   const RAttrText &GetAttrText() const { return fAttrText; }
   CustomAttrs &SetAttrText(const RAttrText &txt) { fAttrText = txt; return *this; }
   RAttrText &AttrText() { return fAttrText; }

   double GetDirect(const std::string &name) const { return GetValue<double>(name); }
};


TEST(OptsTest, AttribStrings) {
   CustomAttrs attrs;

   attrs.AttrBox().AttrBorder().SetWidth(42.);
   attrs.AttrText().SetSize(1.7);

   {
      auto val = attrs.GetDirect("box_border_width");
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
   auto &border = attrs.AttrBox().AttrBorder();
   border.SetWidth(42.);

   {
      // Value was set on this attr, not coming from style:
      EXPECT_FLOAT_EQ(attrs.GetAttrBox().GetAttrBorder().GetWidth(), 42.f);
      EXPECT_FLOAT_EQ(border.GetWidth(), 42.f);
   }

   {
      // Value was set on this attr, not coming from style:
      EXPECT_EQ(attrs.GetAttrText().GetColor(), RColor::kBlue);
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

   auto &al1 = attrs.AttrBox().AttrBorder();
   auto &al2 = attrs.AttrBox().AttrBorder();
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

   attrs1.AttrBox().AttrBorder().SetWidth(7.);
   EXPECT_NE(attrs1, attrs2);
   EXPECT_NE(attrs2, attrs1);
   EXPECT_NE(attrs1, attrs3);
   EXPECT_EQ(attrs2, attrs3);
   EXPECT_EQ(attrs3, attrs2);

   attrs2.AttrBox().AttrBorder().SetColor(RColor::kRed);
   EXPECT_NE(attrs1, attrs2);
   EXPECT_NE(attrs2, attrs1);
   EXPECT_NE(attrs1, attrs3);
   EXPECT_NE(attrs2, attrs3);
   EXPECT_NE(attrs3, attrs2);
}


TEST(OptsTest, AttribAssign) {
   CustomAttrs attrs1;
   CustomAttrs attrs2;

   // deep copy - independent from origin
   auto attrBox1 = attrs1.GetAttrBox();
   auto attrBox2 = attrs2.GetAttrBox();

   EXPECT_EQ(attrBox2, attrBox1);
   EXPECT_EQ(attrBox1, attrBox2);

   attrBox1.AttrBorder().SetWidth(42.);
   EXPECT_NE(attrBox2, attrBox1);

   attrBox2 = attrBox1;
   EXPECT_EQ(attrBox2, attrBox1);
   EXPECT_EQ(attrBox1, attrBox2);

   // But original attributes now differ
   EXPECT_NE(attrs1.GetAttrBox(), attrBox1);
   EXPECT_NE(attrs2.GetAttrBox(), attrBox2);

   EXPECT_FLOAT_EQ(attrBox1.GetAttrBorder().GetWidth(), 42.);
   EXPECT_FLOAT_EQ(attrBox2.GetAttrBorder().GetWidth(), 42.);
   // default width return 1
   EXPECT_FLOAT_EQ(attrs1.GetAttrBox().GetAttrBorder().GetWidth(), 1.);
   EXPECT_FLOAT_EQ(attrs2.GetAttrBox().GetAttrBorder().GetWidth(), 1.);

   // Are the two attributes disconnected?
   attrBox2.AttrBorder().SetWidth(3.);
   EXPECT_EQ(attrs1.GetAttrBox().GetAttrBorder(), attrs2.GetAttrBox().GetAttrBorder());
   EXPECT_FLOAT_EQ(attrBox1.GetAttrBorder().GetWidth(), 42.);
   EXPECT_FLOAT_EQ(attrBox2.GetAttrBorder().GetWidth(), 3.);
   EXPECT_FLOAT_EQ(attrs1.GetAttrBox().GetAttrBorder().GetWidth(), 1.);
   EXPECT_FLOAT_EQ(attrs2.GetAttrBox().GetAttrBorder().GetWidth(), 1.);
}
