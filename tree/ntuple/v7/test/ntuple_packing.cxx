#include "ntuple_test.hxx"

TEST(Packing, Bitfield)
{
   ROOT::Experimental::Detail::RColumnElement<bool, ROOT::Experimental::EColumnType::kBit> element(nullptr);
   element.Pack(nullptr, nullptr, 0);
   element.Unpack(nullptr, nullptr, 0);

   bool b = true;
   char c = 0;
   element.Pack(&c, &b, 1);
   EXPECT_EQ(1, c);
   bool e = false;
   element.Unpack(&e, &c, 1);
   EXPECT_TRUE(e);

   bool b8[] = {true, false, true, false, false, true, false, true};
   c = 0;
   element.Pack(&c, &b8, 8);
   bool e8[] = {false, false, false, false, false, false, false, false};
   element.Unpack(&e8, &c, 8);
   for (unsigned i = 0; i < 8; ++i) {
      EXPECT_EQ(b8[i], e8[i]);
   }

   bool b9[] = {true, false, true, false, false, true, false, true, true};
   char c2[2];
   element.Pack(&c2, &b9, 9);
   bool e9[] = {false, false, false, false, false, false, false, false, false};
   element.Unpack(&e9, &c2, 9);
   for (unsigned i = 0; i < 9; ++i) {
      EXPECT_EQ(b9[i], e9[i]);
   }
}
