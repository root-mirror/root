// @(#)root/mathcore:$Id$
// Author: Jonas Hahnfeld 11/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdint>

/// Compute `a + b` and set `overflow` accordingly.
static inline uint64_t add_overflow(uint64_t a, uint64_t b, unsigned &overflow)
{
   uint64_t add = a + b;
   overflow = (add < a);
   return add;
}

/// Compute `a + b` and increment `carry` if there was an overflow
static inline uint64_t add_carry(uint64_t a, uint64_t b, unsigned &carry)
{
   unsigned overflow;
   uint64_t add = add_overflow(a, b, overflow);
   // Do NOT branch on overflow to avoid jumping code, just add 0 if there was
   // no overflow.
   carry += overflow;
   return add;
}

/// Compute `a - b` and set `overflow` accordingly
static inline uint64_t sub_overflow(uint64_t a, uint64_t b, unsigned &overflow)
{
   uint64_t sub = a - b;
   overflow = (sub > a);
   return sub;
}

/// Compute `a - b` and increment `carry` if there was an overflow
static inline uint64_t sub_carry(uint64_t a, uint64_t b, unsigned &carry)
{
   unsigned overflow;
   uint64_t sub = sub_overflow(a, b, overflow);
   // Do NOT branch on overflow to avoid jumping code, just add 0 if there was
   // no overflow.
   carry += overflow;
   return sub;
}

/// Multiply two 576 bit numbers, stored as 9 numbers of 64 bits each
///
/// \param[in] in1 first factor as 9 numbers of 64 bits each
/// \param[in] in2 second factor as 9 numbers of 64 bits each
/// \param[out] out result with 18 numbers of 64 bits each
void multiply9x9(const uint64_t *in1, const uint64_t *in2, uint64_t *out)
{
   uint64_t next = 0;
   unsigned nextCarry = 0;

#if defined(__clang__) || defined(__INTEL_COMPILER) || defined(__CUDACC__)
#pragma unroll
#elif defined(__GNUC__) && __GNUC__ >= 8
// This pragma was introduced in GCC version 8.
#pragma GCC unroll 18
#endif
   for (int i = 0; i < 18; i++) {
      uint64_t current = next;
      unsigned carry = nextCarry;

      next = 0;
      nextCarry = 0;

#if defined(__clang__) || defined(__INTEL_COMPILER) || defined(__CUDACC__)
#pragma unroll
#elif defined(__GNUC__) && __GNUC__ >= 8
// This pragma was introduced in GCC version 8.
#pragma GCC unroll 9
#endif
      for (int j = 0; j < 9; j++) {
         int k = i - j;
         if (k < 0 || k >= 9)
            continue;

         uint64_t fac1 = in1[j];
         uint64_t fac2 = in2[k];
#if defined(__SIZEOF_INT128__) && !defined(ROOT_NO_INT128)
         unsigned __int128 prod = fac1;
         prod = prod * fac2;

         uint64_t upper = prod >> 64;
         uint64_t lower = static_cast<uint64_t>(prod);
#else
         uint64_t upper1 = fac1 >> 32;
         uint64_t lower1 = static_cast<uint32_t>(fac1);

         uint64_t upper2 = fac2 >> 32;
         uint64_t lower2 = static_cast<uint32_t>(fac2);

         // Multiply 32-bit parts, each product has a maximum value of
         // (2 ** 32 - 1) ** 2 = 2 ** 64 - 2 * 2 ** 32 + 1.
         uint64_t upper = upper1 * upper2;
         uint64_t middle1 = upper1 * lower2;
         uint64_t middle2 = lower1 * upper2;
         uint64_t lower = lower1 * lower2;

         // When adding the two products, the maximum value for middle is
         // 2 * 2 ** 64 - 4 * 2 ** 32 + 2, which exceeds a uint64_t.
         unsigned overflow;
         uint64_t middle = add_overflow(middle1, middle2, overflow);
         // Handling the overflow by a multiplication with 0 or 1 is cheaper
         // than branching with an if statement, which the compiler does not
         // optimize to this equivalent code. Note that we could do entirely
         // without this overflow handling when summing up the intermediate
         // products differently as described in the following SO answer:
         //    https://stackoverflow.com/a/51587262
         // However, this approach takes at least the same amount of thinking
         // why a) the code gives the same results without b) overflowing due
         // to the mixture of 32 bit arithmetic. Moreover, my tests show that
         // the scheme implemented here is actually slightly more performant.
         uint64_t overflow_add = overflow * (uint64_t(1) << 32);
         // This addition can never overflow because the maximum value of upper
         // is 2 ** 64 - 2 * 2 ** 32 + 1 (see above). When now adding another
         // 2 ** 32, the result is 2 ** 64 - 2 ** 32 + 1 and still smaller than
         // the maximum 2 ** 64 - 1 that can be stored in a uint64_t.
         upper += overflow_add;

         uint64_t middle_upper = middle >> 32;
         uint64_t middle_lower = middle << 32;

         lower = add_overflow(lower, middle_lower, overflow);
         upper += overflow;

         // This still can't overflow since the maximum of middle_upper is
         //  - 2 ** 32 - 4 if there was an overflow for middle above, bringing
         //    the maximum value of upper to 2 ** 64 - 2.
         //  - otherwise upper still has the initial maximum value given above
         //    and the addition of a value smaller than 2 ** 32 brings it to
         //    a maximum value of 2 ** 64 - 2 ** 32 + 2.
         // (Both cases include the increment to handle the overflow in lower.)
         //
         // All the reasoning makes perfect sense given that the product of two
         // 64 bit numbers is smaller than or equal to
         //     (2 ** 64 - 1) ** 2 = 2 ** 128 - 2 * 2 ** 64 + 1
         // with the upper bits matching the 2 ** 64 - 2 of the first case.
         upper += middle_upper;
#endif

         // Add to current, remember carry.
         current = add_carry(current, lower, carry);

         // Add to next, remember nextCarry.
         next = add_carry(next, upper, nextCarry);
      }

      next = add_carry(next, carry, nextCarry);

      out[i] = current;
   }
}

/// Compute a value congruent to mul modulo m less than 2 ** 576
///
/// \param[in] mul product from multiply9x9 with 18 numbers of 64 bits each
/// \param[out] out result with 9 numbers of 64 bits each
///
/// \f$ m = 2^{576} - 2^{240} + 1 \f$
///
/// Note that this function does *not* return the smallest value congruent to
/// the modulus, it only guarantees a value smaller than \f$ 2^{576} \$!
void mod_m(const uint64_t *mul, uint64_t *out)
{
   uint64_t r[9] = {0};

   // r = t0 - t1 (24 * 24 = 576 bits)
   unsigned carry = 0;
   for (int i = 0; i < 9; i++) {
      uint64_t t0_i = mul[i];
      uint64_t r_i = sub_overflow(t0_i, carry, carry);

      uint64_t t1_i = mul[i + 9];
      r_i = sub_carry(r_i, t1_i, carry);
      r[i] = r_i;
   }
   int64_t c = -((int64_t)carry);

   // r -= t2 (only 240 bits, so need to extend)
   carry = 0;
   for (int i = 0; i < 9; i++) {
      uint64_t r_i = r[i];
      r_i = sub_overflow(r_i, carry, carry);

      uint64_t t2_bits = 0;
      if (i < 4) {
         t2_bits += mul[i + 14] >> 16;
         if (i < 3) {
            t2_bits += mul[i + 15] << 48;
         }
      }
      r_i = sub_carry(r_i, t2_bits, carry);
      r[i] = r_i;
   }
   c -= carry;

   // r += (t3 + t2) * 2 ** 240
   carry = 0;
   {
      uint64_t r_3 = r[3];
      // 16 upper bits
      uint64_t t2_bits = (mul[14] >> 16) << 48;
      uint64_t t3_bits = (mul[9] << 48);

      r_3 = add_carry(r_3, t2_bits, carry);
      r_3 = add_carry(r_3, t3_bits, carry);

      r[3] = r_3;
   }
   for (int i = 0; i < 3; i++) {
      uint64_t r_i = r[i + 4];
      r_i = add_overflow(r_i, carry, carry);

      uint64_t t2_bits = (mul[14 + i] >> 32) + (mul[15 + i] << 32);
      uint64_t t3_bits = (mul[9 + i] >> 16) + (mul[10 + i] << 48);

      r_i = add_carry(r_i, t2_bits, carry);
      r_i = add_carry(r_i, t3_bits, carry);

      r[i + 4] = r_i;
   }
   {
      uint64_t r_7 = r[7];
      r_7 = add_overflow(r_7, carry, carry);

      uint64_t t2_bits = (mul[17] >> 32);
      uint64_t t3_bits = (mul[12] >> 16) + (mul[13] << 48);

      r_7 = add_carry(r_7, t2_bits, carry);
      r_7 = add_carry(r_7, t3_bits, carry);

      r[7] = r_7;
   }
   {
      uint64_t r_8 = r[8];
      r_8 = add_overflow(r_8, carry, carry);

      uint64_t t3_bits = (mul[13] >> 16) + (mul[14] << 48);

      r_8 = add_carry(r_8, t3_bits, carry);

      r[8] = r_8;
   }
   c += carry;

   // c = floor(r / 2 ** 576) has been computed along the way via the carry
   // flags. Now to update r = r - c * m, it suffices to know c * (-2 ** 240 + 1)
   // because the 2 ** 576 will cancel out. Also note that c may be zero, but
   // the operation is still performed to avoid branching.

   // c * (-2 ** 240 + 1) in 576 bits looks as follows, depending on c:
   //  - if c = 0, the number is zero.
   //  - if c = 1: bits 576 to 240 are set,
   //              bits 239 to 1 are zero, and
   //              the last one is set
   //  - if c = -1, which corresponds to all bits set (signed int64_t):
   //              bits 576 to 240 are zero and the rest is set.
   // Note that all bits except the last are exactly complimentary (unless c = 0)
   // and the last byte is conveniently represented by c already.
   // Now construct the three bit patterns from c, their names correspond to the
   // assembly implementation by Alexei Sibidanov.

   // c = 0 -> t0 = 0; c = 1 -> t0 = 0; c = -1 -> all bits set (sign extension)
   // (The assembly implementation shifts by 63, which gives the same result.)
   int64_t t0 = c >> 1;

   // c = 0 -> t2 = 0; c = 1 -> upper 16 bits set; c = -1 -> lower 48 bits set
   int64_t t2 = t0 - (c << 48);

   // c = 0 -> t1 = 0; c = 1 -> all bits set; c = -1 -> t1 = 0
   // (The assembly implementation shifts by 63, which gives the same result.)
   int64_t t1 = t2 >> 48;

   carry = 0;
   {
      uint64_t r_0 = r[0];

      uint64_t out_0 = sub_carry(r_0, c, carry);
      out[0] = out_0;
   }
   for (int i = 1; i < 3; i++) {
      uint64_t r_i = r[i];
      r_i = sub_overflow(r_i, carry, carry);

      uint64_t out_i = sub_carry(r_i, t0, carry);
      out[i] = out_i;
   }
   {
      uint64_t r_3 = r[3];
      r_3 = sub_overflow(r_3, carry, carry);

      uint64_t out_3 = sub_carry(r_3, t2, carry);
      out[3] = out_3;
   }
   for (int i = 4; i < 9; i++) {
      uint64_t r_i = r[i];
      r_i = sub_overflow(r_i, carry, carry);

      uint64_t out_i = sub_carry(r_i, t1, carry);
      out[i] = out_i;
   }
}

/// Combine multiply9x9 and mod_m with internal temporary storage
///
/// \param[in] in1 first factor with 9 numbers of 64 bits each
/// \param[inout] inout second factor and also the output of the same size
void mulmod(const uint64_t *in1, uint64_t *inout)
{
   uint64_t mul[2 * 9] = {0};
   multiply9x9(in1, inout, mul);
   mod_m(mul, inout);
}

/// Compute base to the n modulo m
///
/// \param[in] base with 9 numbers of 64 bits each
/// \param[out] res output with 9 numbers of 64 bits each
/// \param[in] n exponent
///
/// The arguments base and res may point to the same location.
void powermod(const uint64_t *base, uint64_t *res, uint64_t n)
{
   uint64_t fac[9] = {0};
   fac[0] = base[0];
   res[0] = 1;
   for (int i = 1; i < 9; i++) {
      fac[i] = base[i];
      res[i] = 0;
   }

   uint64_t mul[18] = {0};
   while (n) {
      if (n & 1) {
         multiply9x9(res, fac, mul);
         mod_m(mul, res);
      }
      n >>= 1;
      if (!n)
         break;
      multiply9x9(fac, fac, mul);
      mod_m(mul, fac);
   }
}
