// @(#)root/graf2d:$Id$
// Author: Sergey Linev <s.linev@gsi.de>, 2019-10-04


#include "gtest/gtest.h"

#include "ROOT/RStyle.hxx"

#include "ROOT/RDrawable.hxx"
#include "ROOT/RAttrText.hxx"
#include "ROOT/RAttrLine.hxx"
#include "ROOT/RAttrBox.hxx"


using namespace ROOT::Experimental;

class CustomDrawable : public RDrawable {
   RAttrLine  fAttrLine{this, "line"};    ///<! line attributes
   RAttrBox   fAttrBox{this, "box"};      ///<! box attributes
   RAttrText  fAttrText{this, "text"};    ///<! text attributes

public:
   CustomDrawable() : RDrawable("custom") {}

   const RAttrLine &GetAttrLine() const { return fAttrLine; }
   CustomDrawable &SetAttrLine(const RAttrLine &attr) { fAttrLine = attr; return *this; }
   RAttrLine &AttrLine() { return fAttrLine; }

   const RAttrBox &GetAttrBox() const { return fAttrBox; }
   CustomDrawable &SetAttrBox(RAttrBox &box) { fAttrBox = box; return *this; }
   RAttrBox &AttrBox() { return fAttrBox; }

   const RAttrText &GetAttrText() const { return fAttrText; }
   CustomDrawable &SetAttrText(const RAttrText &attr) { fAttrText = attr; return *this; }
   RAttrText &AttrText() { return fAttrText; }
};


TEST(RStyleTest, CreateStyle)
{
   auto style = std::make_shared<RStyle>();

   style->AddBlock("custom").AddDouble("line_width", 2.);

   style->AddBlock("#customid").AddInt("box_fill_style", 5);

   style->AddBlock(".custom_class").AddDouble("text_size", 3.);

   CustomDrawable drawable;
   drawable.SetId("customid");
   drawable.SetCssClass("custom_class");

   drawable.UseStyle(style);

   EXPECT_DOUBLE_EQ(drawable.GetAttrLine().GetWidth(), 2.);

   EXPECT_EQ(drawable.AttrBox().GetAttrFill().GetStyle(), 5);

   EXPECT_DOUBLE_EQ(drawable.GetAttrText().GetSize(), 3.);
}


TEST(RStyleTest, CreateCss)
{
   auto style = RStyle::Parse("custom { line_width: 2; } #customid { box_fill_style: 5; } .custom_class { text_size: 3; }");

   ASSERT_NE(style, nullptr);

   CustomDrawable drawable;
   drawable.SetId("customid");
   drawable.SetCssClass("custom_class");

   drawable.UseStyle(style);

   EXPECT_DOUBLE_EQ(drawable.GetAttrLine().GetWidth(), 2.);

   EXPECT_EQ(drawable.AttrBox().GetAttrFill().GetStyle(), 5);

   EXPECT_DOUBLE_EQ(drawable.GetAttrText().GetSize(), 3.);
}


TEST(RStyleTest, LostStyle)
{
   CustomDrawable drawable;

   {
      auto style = std::make_shared<RStyle>();

      style->AddBlock("custom").AddDouble("line_width", 2.);

      // here weak_ptr will be set, therefore after style is deleted drawable will loose it
      drawable.UseStyle(style);

      EXPECT_DOUBLE_EQ(drawable.GetAttrLine().GetWidth(), 2.);
   }

   // here style no longer exists
   EXPECT_DOUBLE_EQ(drawable.GetAttrLine().GetWidth(), 1.);
}

