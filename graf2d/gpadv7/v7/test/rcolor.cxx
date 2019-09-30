// @(#)root/graf2d:$Id$
// Author: Sergey Linev <s.linev@gsi.de>, 2019-09-30

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

 #include "gtest/gtest.h"

#include "ROOT/RColor.hxx"

using namespace ROOT::Experimental;

// Test usage of empty color
TEST(RColorTest, Empty) {
   RColor col;
   EXPECT_EQ(col.GetHex(), "");
   EXPECT_DOUBLE_EQ(col.GetAlpha(), 1.);
   EXPECT_EQ(col.GetName(), "");
}

// Test usage of empty color
TEST(RColorTest, AsHex) {
   RColor col;
   col.SetRGB(0,0,0);
   EXPECT_EQ(col.GetHex(), "000000");
   col.SetRGB(1,7,15);
   EXPECT_EQ(col.GetHex(), "01070F");
   col.SetRGB(127,127,127);
   EXPECT_EQ(col.GetHex(), "7F7F7F");
   col.SetRGB(255,255,255);
   EXPECT_EQ(col.GetHex(), "FFFFFF");
}

// Test usage of empty color
TEST(RColorTest, Components) {
   RColor col;
   col.SetHex("012345");
   EXPECT_EQ(col.GetRed(), 0x01);
   EXPECT_EQ(col.GetGreen(), 0x23);
   EXPECT_EQ(col.GetBlue(), 0x45);
}

TEST(RColorTest, Alpha) {

   static constexpr double delta = 0.01; // approx precision of alpha storage

   RColor col;
   col.SetAlpha(0);
   EXPECT_DOUBLE_EQ(col.GetAlpha(), 0.);

   col.SetAlpha(1);
   EXPECT_DOUBLE_EQ(col.GetAlpha(), 1.);

   col.SetAlpha(0.1);
   EXPECT_NEAR(0.1,col.GetAlpha(), delta);

   col.SetAlpha(0.5);
   EXPECT_NEAR(0.5,col.GetAlpha(), delta);

   col.SetAlpha(0.8);
   EXPECT_NEAR(0.8,col.GetAlpha(), delta);
}


TEST(RColorTest, Predef) {
   {
      RColor col{RColor::kRed};
      EXPECT_EQ(col.GetHex(), "FF0000");
      EXPECT_EQ(col.HasAlpha(), false);
   }

   {
      RColor col{RColor::kBlue};
      col.SetAlpha(RColor::kTransparent);
      EXPECT_EQ(col.GetHex(), "0000FF");
      EXPECT_FLOAT_EQ(col.GetAlpha(), 0.);
   }
}
