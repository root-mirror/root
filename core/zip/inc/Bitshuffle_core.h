/*
 * Bitshuffle - Filter for improving compression of typed binary data.
 *
 * This file is part of Bitshuffle
 * Author: Kiyoshi Masui <kiyo@physics.ubc.ca>
 * Website: http://www.github.com/kiyo-masui/bitshuffle
 * Created: 2014
 *
 * See LICENSE file for details about copyright and rights to use.
 *
 *
 * Header File
 *
 * Worker routines return an int64_t which is the number of bytes processed
 * if positive or an error code if negative.
 *
 * Error codes:
 *      -1    : Failed to allocate memory.
 *      -11   : Missing SSE.
 *      -12   : Missing AVX.
 *      -13   : Missing Arm Neon.
 *      -80   : Input size not a multiple of 8.
 *      -81   : block_size not multiple of 8.
 *      -91   : Decompression error, wrong number of bytes processed.
 *      -1YYY : Error internal to compression routine with error code -YYY.
 */


#ifndef BITSHUFFLE_CORE_H
#define BITSHUFFLE_CORE_H

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


// These are usually set in the setup.py.
#ifndef BSHUF_VERSION_MAJOR
#define BSHUF_VERSION_MAJOR 0
#define BSHUF_VERSION_MINOR 3
#define BSHUF_VERSION_POINT 5
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* ---- bshuf_default_block_size ----
 *
 * The default block size as function of element size.
 *
 * This is the block size used by the blocked routines (any routine
 * taking a *block_size* argument) when the block_size is not provided
 * (zero is passed).
 *
 * The results of this routine are guaranteed to be stable such that
 * shuffled/compressed data can always be decompressed.
 *
 * Parameters
 * ----------
 *  elem_size : element size of data to be shuffled/compressed.
 *
 */
size_t bshuf_default_block_size(const size_t elem_size);


#ifdef __cplusplus
} // extern "C"
#endif

#endif  // BITSHUFFLE_CORE_H
