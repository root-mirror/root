// @(#)root/mathcore:$Id$
// Author: Jonas Hahnfeld 11/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "../src/mulmod.h"

#include "gtest/gtest.h"

#include <cstdint>

TEST(multiply9x9, simple)
{
   uint64_t a[9] = {2, 0, 0, 0, 0, 0, 0, 0, 0};
   uint64_t b[9] = {3, 0, 0, 0, 0, 0, 0, 0, 0};
   uint64_t mul[18];

   multiply9x9(a, b, mul);

   EXPECT_EQ(mul[0], 6);
   for (int i = 1; i < 18; i++) {
      EXPECT_EQ(mul[i], 0);
   }
}

TEST(multiply9x9, full)
{
   uint64_t a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
   uint64_t b[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
   uint64_t mul[18];

   multiply9x9(a, b, mul);

   EXPECT_EQ(mul[0], 1);
   EXPECT_EQ(mul[1], 4);
   EXPECT_EQ(mul[2], 10);
   EXPECT_EQ(mul[3], 20);
   EXPECT_EQ(mul[4], 35);
   EXPECT_EQ(mul[5], 56);
   EXPECT_EQ(mul[6], 84);
   EXPECT_EQ(mul[7], 120);
   EXPECT_EQ(mul[8], 165);
   EXPECT_EQ(mul[9], 200);
   EXPECT_EQ(mul[10], 224);
   EXPECT_EQ(mul[11], 236);
   EXPECT_EQ(mul[12], 235);
   EXPECT_EQ(mul[13], 220);
   EXPECT_EQ(mul[14], 190);
   EXPECT_EQ(mul[15], 144);
   EXPECT_EQ(mul[16], 81);
   EXPECT_EQ(mul[17], 0);
}

TEST(multiply9x9, max)
{
   uint64_t max = UINT64_MAX;
   uint64_t a[9] = {max, max, max, max, max, max, max, max, max};
   uint64_t b[9] = {max, max, max, max, max, max, max, max, max};
   uint64_t mul[18];

   multiply9x9(a, b, mul);

   // Expected result verified using the assembly implementation of ranluxppp
   // and the following Python program:
   //     >>> a = b = (1 << 576) - 1
   //     >>> mul = a * b
   //     >>> while mul > 0:
   //     ...     print(hex(mul & 0xffffffffffffffff))
   //     ...     mul = mul >> 64
   //     ...
   EXPECT_EQ(mul[0], 1);
   for (int i = 1; i < 9; i++) {
      EXPECT_EQ(mul[i], 0);
   }
   EXPECT_EQ(mul[9], 0xfffffffffffffffe);
   for (int i = 10; i < 18; i++) {
      EXPECT_EQ(mul[i], 0xffffffffffffffff);
   }
}

// The modulus m = 2 ** 576 - 2 ** 240 + 1.

TEST(mod_m, nop)
{
   // mul = 2 ** 512 < m
   uint64_t mul[18] = {0};
   mul[8] = 1;
   uint64_t mod[9];

   mod_m(mul, mod);

   for (int i = 0; i < 8; i++) {
      EXPECT_EQ(mod[i], 0);
   }
   EXPECT_EQ(mod[8], 1);
}

TEST(mod_m, modulus)
{
   // mul = m = 2 ** 576 - 2 ** 240 + 1
   uint64_t mul[18] = {0};
   mul[0] = 1;
   mul[3] = 0xffff000000000000;
   for (int i = 4; i < 9; i++) {
      mul[i] = 0xffffffffffffffff;
   }
   uint64_t mod[9];

   mod_m(mul, mod);

   // mod_m does not provide the smallest value, but just a value smaller than
   // 2 ** 576 that is congruent to mul; in this case, mul itself.
   for (int i = 0; i < 9; i++) {
      EXPECT_EQ(mod[i], mul[i]);
   }
}

TEST(mod_m, simple)
{
   // mul = 2 ** 576
   uint64_t mul[18] = {0};
   mul[9] = 1;
   uint64_t mod[9];

   mod_m(mul, mod);

   // The result is mod = mul - m = 2 ** 240 - 1, all first 240 bits set.
   for (int i = 0; i < 3; i++) {
      EXPECT_EQ(mod[i], 0xffffffffffffffff);
   }
   EXPECT_EQ(mod[3], 0x0000ffffffffffff);
   for (int i = 4; i < 9; i++) {
      EXPECT_EQ(mod[i], 0);
   }
}

TEST(mod_m, pattern)
{
   uint64_t pattern1 = 0xdeafdeafdeafdeaf;
   uint64_t pattern2 = 0xf00df00df00df00d;
   uint64_t mul[18] = {
      pattern1, pattern2, pattern1, pattern2, pattern1, pattern1, pattern2, pattern2, pattern2,
      pattern1, pattern2, pattern1, pattern2, pattern2, pattern1, pattern1, pattern2, pattern1,
   };
   uint64_t mod[9];

   mod_m(mul, mod);

   // Expected result verified using the assembly implementation of ranluxppp
   // and the following Python program:
   //     >>> pattern1 = 0xdeafdeafdeafdeaf
   //     >>> pattern2 = 0xf00df00df00df00d
   //     >>> mul = 0
   //     >>> for p in [
   //     ...     pattern1, pattern2, pattern1, pattern2,
   //     ...     pattern1, pattern1, pattern2, pattern2,
   //     ...     pattern2, pattern1, pattern2, pattern1,
   //     ...     pattern2, pattern2, pattern1, pattern1
   //     ... ]:
   //     ...     mul = (mul << 64) + p
   //     ...
   //     >>> m = (1 << 576) - (1 << 240) + 1
   //     >>> mod = mul % m
   //     >>> print(hex(mod))
   EXPECT_EQ(mod[8], 0xf00e016c016c016b);
   EXPECT_EQ(mod[7], 0xf00df00ecebdcebd);
   EXPECT_EQ(mod[6], 0xe01bcebde01be01b);
   EXPECT_EQ(mod[5], 0xcebde01bcebdcebd);
   EXPECT_EQ(mod[4], 0xbd5fac01ac01ac01);
   EXPECT_EQ(mod[3], 0xbd5d215021502150);
   EXPECT_EQ(mod[2], 0x21500ff20ff20ff2);
   EXPECT_EQ(mod[1], 0x0ff2215021502150);
   EXPECT_EQ(mod[0], 0x2150215021502151);
}

TEST(mod_m, r_min)
{
   uint64_t mul[18] = {0};
   mul[14] = 0xffffffffffff0000;
   for (int i = 15; i < 18; i++) {
      mul[i] = UINT64_MAX;
   }
   uint64_t mod[9];

   mod_m(mul, mod);

   // Expected result produced by the assembly implementation.
   EXPECT_EQ(mod[0], 2);
   EXPECT_EQ(mod[1], 0);
   EXPECT_EQ(mod[2], 0);
   EXPECT_EQ(mod[3], 0xfffd000000000000);
   EXPECT_EQ(mod[4], 0xffffffffffffffff);
   EXPECT_EQ(mod[5], 0xffff);
   EXPECT_EQ(mod[6], 0);
   EXPECT_EQ(mod[7], 0x0000000100000000);
   EXPECT_EQ(mod[8], 0);
}

TEST(mod_m, r_max)
{
   uint64_t mul[18] = {0};
   for (int i = 0; i < 14; i++) {
      mul[i] = UINT64_MAX;
   }
   mul[14] = 0xffff;
   uint64_t mod[9];

   mod_m(mul, mod);

   // Expected result produced by the assembly implementation.
   for (int i = 0; i < 5; i++) {
      EXPECT_EQ(mod[i], 0xffffffffffffffff);
   }
   EXPECT_EQ(mod[5], 0xfffffffffffeffff);
   for (int i = 6; i < 9; i++) {
      EXPECT_EQ(mod[i], 0xffffffffffffffff);
   }
}

TEST(powermod, simple)
{
   uint64_t base[9] = {0};
   base[0] = 2;
   uint64_t res[9];

   powermod(base, res, 15);

   EXPECT_EQ(res[0], 1 << 15);
   for (int i = 1; i < 9; i++) {
      EXPECT_EQ(res[i], 0);
   }
}

TEST(powermod, mod)
{
   uint64_t base[9] = {0};
   base[0] = 2;
   uint64_t res[9];

   powermod(base, res, 576);

   // The same result as mod_m::simple above:
   // mod = 2 ** 576 - m = 2 ** 240 - 1, all first 240 bits set.
   for (int i = 0; i < 3; i++) {
      EXPECT_EQ(res[i], 0xffffffffffffffff);
   }
   EXPECT_EQ(res[3], 0x0000ffffffffffff);
   for (int i = 4; i < 9; i++) {
      EXPECT_EQ(res[i], 0);
   }
}

TEST(powermod, alias)
{
   uint64_t a[9] = {0};
   a[0] = 2;

   // Can pass in the same pointer for base and res, the values should be
   // copied before writing the output.
   powermod(a, a, 15);

   EXPECT_EQ(a[0], 1 << 15);
   for (int i = 1; i < 9; i++) {
      EXPECT_EQ(a[i], 0);
   }
}
