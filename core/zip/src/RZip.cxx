/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "Compression.h"
#include "RConfigure.h"
#include "RZip.h"
#include "Bits.h"
#include "ZipLZMA.h"
#include "ZipLZ4.h"
#include "ZipZSTD.h"

#include "zlib.h"

#include <cstdio>
#include <cassert>

// The size of the ROOT block framing headers for compression:
// - 3 bytes to identify the compression algorithm and version.
// - 3 bytes to identify the deflated buffer size.
// - 3 bytes to identify the inflated buffer size.
#define HDRSIZE 9

/**
 * Forward decl's
 */
static void R__zipOld(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgrt, int *irep);
static void R__zipZLIB(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgrt, int *irep);
static void R__unzipZLIB(int *srcsize, unsigned char *src, int *tgtsize, unsigned char *tgt, int *irep);

/* ===========================================================================
   R__ZipMode is used to select the compression algorithm when R__zip is called
   and when R__zipMultipleAlgorithm is called with its last argument set to 0.
   R__ZipMode = 1 : ZLIB compression algorithm is used (default)
   R__ZipMode = 2 : LZMA compression algorithm is used
   R__ZipMode = 4 : LZ4  compression algorithm is used
   R__ZipMode = 0 or 3 : a very old compression algorithm is used
   (the very old algorithm is supported for backward compatibility)
   The LZMA algorithm requires the external XZ package be installed when linking
   is done. LZMA typically has significantly higher compression factors, but takes
   more CPU time and memory resources while compressing.

  The LZ4 algorithm requires the external LZ4 package to be installed when linking
  is done.  LZ4 typically has the worst compression ratios, but much faster decompression
  speeds - sometimes by an order of magnitude.
*/
#ifdef R__HAS_DEFAULT_ZSTD
ROOT::RCompressionSetting::EAlgorithm::EValues R__ZipMode = ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD;
#elif defined(R__HAS_DEFAULT_LZ4)
ROOT::RCompressionSetting::EAlgorithm::EValues R__ZipMode = ROOT::RCompressionSetting::EAlgorithm::EValues::kLZ4;
#else
ROOT::RCompressionSetting::EAlgorithm::EValues R__ZipMode = ROOT::RCompressionSetting::EAlgorithm::EValues::kZLIB;
#endif

/* ===========================================================================
   Function to set the ZipMode
 */
extern "C" void R__SetZipMode(ROOT::RCompressionSetting::EAlgorithm::EValues mode)
{
   R__ZipMode = mode;
}

unsigned long R__crc32(unsigned long crc, const unsigned char* buf, unsigned int len)
{
   return crc32(crc, buf, len);
}

/* int cxlevel;                      compression level */
/* int  *srcsize, *tgtsize, *irep;   source and target sizes, replay */
/* char *tgt, *src;                  source and target buffers */
/* compressionAlgorithm 0 = use global setting */
/*                      1 = zlib */
/*                      2 = lzma */
/*                      3 = old */
void R__zipMultipleAlgorithm(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep, ROOT::RCompressionSetting::EAlgorithm::EValues compressionAlgorithm)
{

  if (*srcsize < 1 + HDRSIZE + 1) {
     *irep = 0;
     return;
  }

  if (cxlevel <= 0) {
    *irep = 0;
    return;
  }

  if (compressionAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kUseGlobal) {
    compressionAlgorithm = R__ZipMode;
  }

  // The LZMA compression algorithm from the XZ package
  if (compressionAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kLZMA) {
     R__zipLZMA(cxlevel, srcsize, src, tgtsize, tgt, irep);
  } else if (compressionAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kLZ4) {
     R__zipLZ4(cxlevel, srcsize, src, tgtsize, tgt, irep);
  } else if (compressionAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kZSTD) {
     R__zipZSTD(cxlevel, srcsize, src, tgtsize, tgt, irep);
  } else if (compressionAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kOldCompressionAlgo || compressionAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kUseGlobal) {
     R__zipOld(cxlevel, srcsize, src, tgtsize, tgt, irep);
  } else {
     // 1 is for ZLIB (which is the default), ZLIB is also used for any illegal
     // algorithm setting.  This was a poor historic choice, as poor code may result in
     // a surprising change in algorithm in a future version of ROOT.
     R__zipZLIB(cxlevel, srcsize, src, tgtsize, tgt, irep);
  }
}

  // The very old algorithm for backward compatibility
  // 0 for selecting with R__ZipMode in a backward compatible way
  // 3 for selecting in other cases
static void R__zipOld(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep)
{
    int method   = Z_DEFLATED;
    bits_internal_state state;
    ush att      = (ush)UNKNOWN;
    ush flags    = 0;
    if (cxlevel > 9) cxlevel = 9;
    gCompressionLevel = cxlevel;

    *irep        = 0;
    /* error_flag   = 0; */
    if (*tgtsize <= 0) {
       R__error("target buffer too small");
       return;
    }
    if (*srcsize > 0xffffff) {
       R__error("source buffer too big");
       return;
    }

#ifdef DYN_ALLOC
    state.R__window = 0;
    state.R__prev = 0;
#endif

    state.in_buf    = src;
    state.in_size   = (unsigned) (*srcsize);
    state.in_offset = 0;

    state.out_buf     = tgt;
    state.out_size    = (unsigned) (*tgtsize);
    state.out_offset  = HDRSIZE;
    state.R__window_size = 0L;

    if (0 != R__bi_init(&state) ) return;       /* initialize bit routines */
    state.t_state = R__get_thread_tree_state();
    if (0 != R__ct_init(state.t_state,&att, &method)) return; /* initialize tree routines */
    if (0 != R__lm_init(&state, gCompressionLevel, &flags)) return; /* initialize compression */
    R__Deflate(&state,&state.error_flag);                  /* compress data */
    if (state.error_flag != 0) return;

    tgt[0] = 'C';               /* Signature 'C'-Chernyaev, 'S'-Smirnov */
    tgt[1] = 'S';
    tgt[2] = (char) method;

    state.out_size  = state.out_offset - HDRSIZE;         /* compressed size */
    tgt[3] = (char)(state.out_size & 0xff);
    tgt[4] = (char)((state.out_size >> 8) & 0xff);
    tgt[5] = (char)((state.out_size >> 16) & 0xff);

    tgt[6] = (char)(state.in_size & 0xff);         /* decompressed size */
    tgt[7] = (char)((state.in_size >> 8) & 0xff);
    tgt[8] = (char)((state.in_size >> 16) & 0xff);

    *irep     = state.out_offset;
    return;
}

/**
 * Compress buffer contents using the venerable zlib algorithm.
 */
static void R__zipZLIB(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep)
{
  int err;
  int method   = Z_DEFLATED;

    z_stream stream;
    //Don't use the globals but want name similar to help see similarities in code
    unsigned l_in_size, l_out_size;
    *irep = 0;

    if (*tgtsize <= 0) {
       R__error("target buffer too small");
       return;
    }
    if (*srcsize > 0xffffff) {
       R__error("source buffer too big");
       return;
    }

    stream.next_in   = (Bytef*)src;
    stream.avail_in  = (uInt)(*srcsize);

    stream.next_out  = (Bytef*)(&tgt[HDRSIZE]);
    stream.avail_out = (uInt)(*tgtsize);

    stream.zalloc    = (alloc_func)0;
    stream.zfree     = (free_func)0;
    stream.opaque    = (voidpf)0;

    if (cxlevel > 9) cxlevel = 9;
    err = deflateInit(&stream, cxlevel);
    if (err != Z_OK) {
       printf("error %d in deflateInit (zlib)\n",err);
       return;
    }

    while ((err = deflate(&stream, Z_FINISH)) != Z_STREAM_END) {
       if (err != Z_OK) {
          deflateEnd(&stream);
          return;
       }
    }

    err = deflateEnd(&stream);

    tgt[0] = 'Z';               /* Signature ZLib */
    tgt[1] = 'L';
    tgt[2] = (char) method;

    l_in_size   = (unsigned) (*srcsize);
    l_out_size  = stream.total_out;             /* compressed size */
    tgt[3] = (char)(l_out_size & 0xff);
    tgt[4] = (char)((l_out_size >> 8) & 0xff);
    tgt[5] = (char)((l_out_size >> 16) & 0xff);

    tgt[6] = (char)(l_in_size & 0xff);         /* decompressed size */
    tgt[7] = (char)((l_in_size >> 8) & 0xff);
    tgt[8] = (char)((l_in_size >> 16) & 0xff);

    *irep = stream.total_out + HDRSIZE;
    return;
}


void R__zip(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep) {
   R__zipMultipleAlgorithm(cxlevel, srcsize, src, tgtsize, tgt, irep,
                           ROOT::RCompressionSetting::EAlgorithm::kUseGlobal);
}

/**
 * Below are the routines for unzipping (inflating) buffers.
 */

static int is_valid_header_zlib(unsigned char *src)
{
   return src[0] == 'Z' && src[1] == 'L' && src[2] == Z_DEFLATED;
}

static int is_valid_header_old(unsigned char *src)
{
   return src[0] == 'C' && src[1] == 'S' && src[2] == Z_DEFLATED;
}

static int is_valid_header_lzma(unsigned char *src)
{
   return src[0] == 'X' && src[1] == 'Z' && src[2] == 0;
}

static int is_valid_header_lz4(unsigned char *src)
{
   return src[0] == 'L' && src[1] == '4';
}

static int is_valid_header_zstd(unsigned char *src)
{
   return src[0] == 'Z' && src[1] == 'S' && src[2] == '\1';
}

static int is_valid_header(unsigned char *src)
{
   return is_valid_header_zlib(src) || is_valid_header_old(src) || is_valid_header_lzma(src) ||
          is_valid_header_lz4(src) || is_valid_header_zstd(src);
}

int R__unzip_header(int *srcsize, uch *src, int *tgtsize)
{
  // Reads header envelope, and determines target size.
  // Returns 0 in case of success.

  *srcsize = 0;
  *tgtsize = 0;

  /*   C H E C K   H E A D E R   */
  if (!is_valid_header(src)) {
     fprintf(stderr, "Error R__unzip_header: error in header.  Values: %x%x\n", src[0], src[1]);
     return 1;
  }

  *srcsize = HDRSIZE + ((long)src[3] | ((long)src[4] << 8) | ((long)src[5] << 16));
  *tgtsize = (long)src[6] | ((long)src[7] << 8) | ((long)src[8] << 16);

  return 0;
}


/***********************************************************************
 *                                                                     *
 * Name: R__unzip                                    Date:    20.01.95 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: In memory ZIP decompression. Can be issued from FORTRAN.  *
 *           Written for DELPHI collaboration (CERN)                   *
 *                                                                     *
 * Input: scrsize - size of input buffer                               *
 *        src     - input buffer                                       *
 *        tgtsize - size of target buffer                              *
 *                                                                     *
 * Output: tgt - target buffer (decompressed)                          *
 *         irep - size of decompressed data                            *
 *                0 - if error                                         *
 *                                                                     *
 ***********************************************************************/
// N.B. (Brian) - I have kept the original note out of complete awe of the
// age of the original code...
void R__unzip(int *srcsize, uch *src, int *tgtsize, uch *tgt, int *irep)
{
   long isize;
   uch *ibufptr, *obufptr;
   long ibufcnt, obufcnt;

   *irep = 0L;

   /*   C H E C K   H E A D E R   */

   if (*srcsize < HDRSIZE) {
      fprintf(stderr, "R__unzip: too small source\n");
      return;
   }

   /*   C H E C K   H E A D E R   */
   if (!is_valid_header(src)) {
      fprintf(stderr, "Error R__unzip: error in header\n");
      return;
   }

   ibufptr = src + HDRSIZE;
   ibufcnt = (long)src[3] | ((long)src[4] << 8) | ((long)src[5] << 16);
   isize = (long)src[6] | ((long)src[7] << 8) | ((long)src[8] << 16);
   obufptr = tgt;
   obufcnt = *tgtsize;

   if (obufcnt < isize) {
      fprintf(stderr, "R__unzip: too small target\n");
      return;
   }

   if (ibufcnt + HDRSIZE != *srcsize) {
      fprintf(stderr, "R__unzip: discrepancy in source length\n");
      return;
   }

   /* ZLIB and other standard compression algorithms */
   if (is_valid_header_zlib(src)) {
      R__unzipZLIB(srcsize, src, tgtsize, tgt, irep);
      return;
   } else if (is_valid_header_lzma(src)) {
      R__unzipLZMA(srcsize, src, tgtsize, tgt, irep);
      return;
   } else if (is_valid_header_lz4(src)) {
      R__unzipLZ4(srcsize, src, tgtsize, tgt, irep);
      return;
   } else if (is_valid_header_zstd(src)) {
      R__unzipZSTD(srcsize, src, tgtsize, tgt, irep);
      return;
   }

   /* Old zlib format */
   if (R__Inflate(&ibufptr, &ibufcnt, &obufptr, &obufcnt)) {
      fprintf(stderr, "R__unzip: error during decompression\n");
      return;
   }

   /* if (obufptr - tgt != isize) {
     There are some rare cases when a few more bytes are required */
   if (obufptr - tgt > *tgtsize) {
      fprintf(stderr, "R__unzip: discrepancy (%ld) with initial size: %ld, tgtsize=%d\n", (long)(obufptr - tgt), isize,
              *tgtsize);
      *irep = obufptr - tgt;
      return;
   }

   *irep = isize;
}

void R__unzipZLIB(int *srcsize, unsigned char *src, int *tgtsize, unsigned char *tgt, int *irep)
{
     z_stream stream; /* decompression stream */
     int err = 0;

     stream.next_in = (Bytef *)(&src[HDRSIZE]);
     stream.avail_in = (uInt)(*srcsize) - HDRSIZE;
     stream.next_out = (Bytef *)tgt;
     stream.avail_out = (uInt)(*tgtsize);
     stream.zalloc = (alloc_func)0;
     stream.zfree = (free_func)0;
     stream.opaque = (voidpf)0;

     err = inflateInit(&stream);
     if (err != Z_OK) {
        fprintf(stderr, "R__unzip: error %d in inflateInit (zlib)\n", err);
        return;
     }

     while ((err = inflate(&stream, Z_FINISH)) != Z_STREAM_END) {
        if (err != Z_OK) {
           inflateEnd(&stream);
           fprintf(stderr, "R__unzip: error %d in inflate (zlib)\n", err);
           return;
        }
     }

     inflateEnd(&stream);

     *irep = stream.total_out;
     return;
}
