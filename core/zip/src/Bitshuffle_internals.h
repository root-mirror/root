/*
 * Bitshuffle - Filter for improving compression of typed binary data.
 *
 * This file is part of Bitshuffle
 * Author: Kiyoshi Masui <kiyo@physics.ubc.ca>
 * Website: http://www.github.com/kiyo-masui/bitshuffle
 * Created: 2014
 *
 * See LICENSE file for details about copyright and rights to use.
 */


#ifndef BITSHUFFLE_INTERNALS_H
#define BITSHUFFLE_INTERNALS_H

// We assume GNU g++ defining `__cplusplus` has stdint.h
#if (defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199900L) || defined(__cplusplus)
#include <stdint.h>
#else
  typedef unsigned char       uint8_t;
  typedef unsigned short      uint16_t;
  typedef unsigned int        uint32_t;
  typedef   signed int        int32_t;
  #if defined __x86_64__ && !defined __ILP32__
  typedef unsigned long       uint64_t;
  typedef long                int64_t;
  #else
  typedef unsigned long long  uint64_t;
  typedef long long           int64_t;
  #endif
#endif

#include <stdlib.h>
#include "Bitshuffle_iochain.h"


// Constants.
#ifndef BSHUF_MIN_RECOMMEND_BLOCK
#define BSHUF_MIN_RECOMMEND_BLOCK 128
#define BSHUF_BLOCKED_MULT 8    // Block sizes must be multiple of this.
#define BSHUF_TARGET_BLOCK_SIZE_B 8192
#endif


// Macros.
#define CHECK_ERR_FREE(count, buf) if (count < 0) { free(buf); return count; }


#ifdef __cplusplus
extern "C" {
#endif

/* ---- Utility functions for internal use only ---- */

int64_t bshuf_trans_bit_elem(const void* in, void* out, const size_t size,
        const size_t elem_size);

/* Read a 32 bit unsigned integer from a buffer big endian order. */
uint32_t bshuf_read_uint32_BE(const void* buf);

/* Write a 32 bit unsigned integer to a buffer in big endian order. */
void bshuf_write_uint32_BE(void* buf, uint32_t num);

int64_t bshuf_untrans_bit_elem(const void* in, void* out, const size_t size,
        const size_t elem_size);

/* Function definition for worker functions that process a single block. */
typedef int64_t (*bshufBlockFunDef)(ioc_chain* C_ptr,
        const size_t size, const size_t elem_size);

/* Wrap a function for processing a single block to process an entire buffer in
 * parallel. */
int64_t bshuf_blocked_wrap_fun(bshufBlockFunDef fun, const void* in, void* out,
        const size_t size, const size_t elem_size, size_t block_size);

#ifdef __cplusplus
} // extern "C"
#endif

#endif  // BITSHUFFLE_INTERNALS_H
