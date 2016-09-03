/* Copyright 2013 Google Inc. All Rights Reserved.

   Distributed under MIT license.
   See file LICENSE for detail or copy at https://opensource.org/licenses/MIT
*/

/* Bit reading helpers */

#include "./bit_reader.h"

#include "./port.h"
#include "./types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

void BrotliInitBitReader(BrotliBitReader* const br) {
  br->val_ = 0;
  br->bit_pos_ = sizeof(br->val_) << 3;
}

int BrotliWarmupBitReader(BrotliBitReader* const br) {
  size_t aligned_read_mask = (sizeof(br->val_) >> 1) - 1;
  /* Fixing alignment after unaligned BrotliFillWindow would result accumulator
     overflow. If unalignment is caused by BrotliSafeReadBits, then there is
     enough space in accumulator to fix aligment. */
  if (!BROTLI_ALIGNED_READ) {
    aligned_read_mask = 0;
  }
  if (BrotliGetAvailableBits(br) == 0) {
    if (!BrotliPullByte(br)) {
      return 0;
    }
  }

  while ((((size_t)br->next_in) & aligned_read_mask) != 0) {
    if (!BrotliPullByte(br)) {
      /* If we consumed all the input, we don't care about the alignment. */
      return 1;
    }
  }
  return 1;
}

#if defined(__cplusplus) || defined(c_plusplus)
}  /* extern "C" */
#endif
