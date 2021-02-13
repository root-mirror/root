#include "gtest/gtest.h"

#include "ROOT/RBox.hxx"
#include "ROOT/RLine.hxx"
#include "ROOT/RMarker.hxx"
#include "ROOT/RText.hxx"
#include "ROOT/RPaveText.hxx"
#include "ROOT/RLegend.hxx"
#include "ROOT/RCanvas.hxx"

using namespace ROOT::Experimental;

// Test RBox API
TEST(Primitives, RBox)
{
   RCanvas canv;
   auto box = canv.Draw<RBox>(RPadPos(0.1_normal, 0.3_normal), RPadPos(0.3_normal,0.6_normal));

   box->AttrBox().AttrBorder().SetColor(RColor::kRed).SetWidth(5.).SetStyle(7);
   box->AttrBox().AttrFill().SetColor(RColor::kBlue).SetStyle(6);

   EXPECT_EQ(canv.NumPrimitives(), 1u);

   EXPECT_EQ(box->GetAttrBox().GetAttrBorder().GetColor(), RColor::kRed);
   EXPECT_DOUBLE_EQ(box->GetAttrBox().GetAttrBorder().GetWidth(), 5.);
   EXPECT_EQ(box->GetAttrBox().GetAttrBorder().GetStyle(), 7);

   EXPECT_EQ(box->GetAttrBox().GetAttrFill().GetColor(), RColor::kBlue);
   EXPECT_EQ(box->GetAttrBox().GetAttrFill().GetStyle(), 6);
}

// Test RLine API
TEST(Primitives, RLine)
{
   RCanvas canv;
   auto line = canv.Draw<RLine>(RPadPos(0.1_normal, 0.1_normal), RPadPos(0.9_normal,0.9_normal));

   line->AttrLine().SetColor(RColor::kRed).SetWidth(5.).SetStyle(7);

   EXPECT_EQ(canv.NumPrimitives(), 1u);

   EXPECT_EQ(line->GetAttrLine().GetColor(), RColor::kRed);
   EXPECT_DOUBLE_EQ(line->GetAttrLine().GetWidth(), 5.);
   EXPECT_EQ(line->GetAttrLine().GetStyle(), 7);
}

// Test RMarker API
TEST(Primitives, RMarker)
{
   RCanvas canv;
   auto marker = canv.Draw<RMarker>(RPadPos(0.5_normal, 0.5_normal));

   marker->AttrMarker().SetStyle(7).SetSize(2.5).SetColor(RColor::kGreen);

   EXPECT_EQ(canv.NumPrimitives(), 1u);

   EXPECT_EQ(marker->GetAttrMarker().GetColor(), RColor::kGreen);
   EXPECT_DOUBLE_EQ(marker->GetAttrMarker().GetSize(), 2.5);
   EXPECT_EQ(marker->GetAttrMarker().GetStyle(), 7);
}

// Test RText API
TEST(Primitives, RText)
{
   RCanvas canv;

   auto text = canv.Draw<RText>(RPadPos(0.5_normal, 0.5_normal), "Hello World");

   text->AttrText().SetColor(RColor::kBlack).SetSize(12.5).SetAngle(90.).SetAlign(13).SetFontFamily("Arial");

   EXPECT_EQ(canv.NumPrimitives(), 1u);

   EXPECT_EQ(text->GetText(), "Hello World");
   EXPECT_EQ(text->GetAttrText().GetColor(), RColor::kBlack);
   EXPECT_DOUBLE_EQ(text->GetAttrText().GetSize(), 12.5);
   EXPECT_DOUBLE_EQ(text->GetAttrText().GetAngle(), 90.);
   EXPECT_EQ(text->GetAttrText().GetAlign(), 13);
   EXPECT_EQ(text->GetAttrText().GetFontFamily(), "Arial");
}


// Test same color functionality
TEST(Primitives, SameColor)
{
   RCanvas canv;
   auto line1 = canv.Draw<RLine>(RPadPos(0.1_normal, 0.1_normal), RPadPos(0.9_normal,0.9_normal));
   auto line2 = canv.Draw<RLine>(RPadPos(0.1_normal, 0.9_normal), RPadPos(0.9_normal,0.1_normal));
   auto line3 = canv.Draw<RLine>(RPadPos(0.9_normal, 0.1_normal), RPadPos(0.1_normal,0.9_normal));

   line1->AttrLine().AttrColor().SetAuto();
   line2->AttrLine().AttrColor().SetAuto();
   line3->AttrLine().AttrColor().SetAuto();

   canv.AssignAutoColors();

   EXPECT_EQ(canv.NumPrimitives(), 3u);

   EXPECT_EQ(line1->GetAttrLine().GetColor(), RColor::kRed);
   EXPECT_EQ(line2->GetAttrLine().GetColor(), RColor::kGreen);
   EXPECT_EQ(line3->GetAttrLine().GetColor(), RColor::kBlue);
}

// Test RLegend API
TEST(Primitives, RLegend)
{
   RCanvas canv;
   auto line1 = canv.Draw<RLine>(RPadPos(0.1_normal, 0.1_normal), RPadPos(0.9_normal,0.9_normal));
   auto line2 = canv.Draw<RLine>(RPadPos(0.1_normal, 0.9_normal), RPadPos(0.9_normal,0.1_normal));
   auto line3 = canv.Draw<RLine>(RPadPos(0.9_normal, 0.1_normal), RPadPos(0.1_normal,0.9_normal));

   line1->AttrLine().SetColor(RColor::kRed);
   line2->AttrLine().SetColor(RColor::kGreen);
   line3->AttrLine().SetColor(RColor::kBlue);

   auto legend = canv.Draw<RLegend>("Legend title");
   legend->AttrFill().SetStyle(5).SetColor(RColor::kWhite);
   legend->AttrBorder().SetWidth(2).SetColor(RColor::kRed);
   legend->AddEntry(line1, "RLine 1");
   legend->AddEntry(line2, "RLine 2");
   legend->AddEntry(line3, "RLine 3");

   EXPECT_EQ(canv.NumPrimitives(), 4u);

   EXPECT_EQ(legend->NumEntries(), 3u);
   EXPECT_EQ(legend->GetTitle(), "Legend title");
   EXPECT_EQ(legend->GetAttrFill().GetColor(), RColor::kWhite);
}

// Test RPaveText API
TEST(Primitives, RPaveText)
{
   RCanvas canv;

   auto text = canv.Draw<RPaveText>();

   text->AttrText().SetColor(RColor::kBlack).SetSize(12).SetAlign(13).SetFontFamily("Times New Roman");
   text->AttrBorder().SetColor(RColor::kRed).SetWidth(3);
   text->AttrFill().SetColor(RColor::kBlue).SetStyle(3003);

   text->AddLine("First line");
   text->AddLine("Second line");
   text->AddLine("Third line");

   EXPECT_EQ(canv.NumPrimitives(), 1u);

   EXPECT_EQ(text->NumLines(), 3u);
   EXPECT_EQ(text->GetLine(0), "First line");
   EXPECT_EQ(text->GetLine(1), "Second line");
   EXPECT_EQ(text->GetLine(2), "Third line");

   EXPECT_EQ(text->GetAttrText().GetColor(), RColor::kBlack);
   EXPECT_DOUBLE_EQ(text->GetAttrText().GetSize(), 12);
   EXPECT_EQ(text->GetAttrText().GetAlign(), 13);
   EXPECT_EQ(text->GetAttrText().GetFontFamily(), "Times New Roman");

   EXPECT_EQ(text->GetAttrBorder().GetColor(), RColor::kRed);
   EXPECT_EQ(text->GetAttrBorder().GetWidth(), 3);

   EXPECT_EQ(text->GetAttrFill().GetColor(), RColor::kBlue);
   EXPECT_EQ(text->GetAttrFill().GetStyle(), 3003);
}

