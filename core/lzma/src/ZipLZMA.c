// @(#)root/lzma:$Id$
// Author: David Dagenhart   May 2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef _MSC_VER
#define LZMA_API_STATIC
#endif
#include "ZipLZMA.h"
#include "lzma.h"
#include <stdio.h>

static const int kHeaderSize = 9;

void R__zipLZMA(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep)
{
   uint64_t out_size;             /* compressed size */
   unsigned in_size   = (unsigned) (*srcsize);
   uint32_t dict_size_est = in_size/4;
   lzma_stream stream = LZMA_STREAM_INIT;
   lzma_options_lzma opt_lzma2;
   lzma_filter filters[] = {
      { .id = LZMA_FILTER_LZMA2, .options = &opt_lzma2 },
      { .id = LZMA_VLI_UNKNOWN,  .options = NULL },
   };
   lzma_ret returnStatus;

   *irep = 0;

   if (*tgtsize <= 0) {
      return;
   }

   if (*srcsize > 0xffffff || *srcsize < 0) {
      return;
   }

   if (cxlevel > 9) cxlevel = 9;

   if (lzma_lzma_preset(&opt_lzma2, cxlevel)) {
      return;
   }

   if (LZMA_DICT_SIZE_MIN > dict_size_est) {
      dict_size_est = LZMA_DICT_SIZE_MIN;
   }
   if (opt_lzma2.dict_size > dict_size_est) {
      /* reduce the dictionary size if larger than 1/4 the input size, preset
         dictionaries size can be expensively large
       */
      opt_lzma2.dict_size = dict_size_est;
   }

   returnStatus = lzma_stream_encoder(&stream,
                                      filters,
                                      LZMA_CHECK_CRC32);
   if (returnStatus != LZMA_OK) {
      return;
   }

   stream.next_in   = (const uint8_t *)src;
   stream.avail_in  = (size_t)(*srcsize);

   stream.next_out  = (uint8_t *)(&tgt[kHeaderSize]);
   stream.avail_out = (size_t)(*tgtsize);

   returnStatus = lzma_code(&stream, LZMA_FINISH);
   if (returnStatus != LZMA_STREAM_END) {
      /* No need to print an error message. We simply abandon the compression
         the buffer cannot be compressed or compressed buffer would be larger than original buffer
      */
      lzma_end(&stream);
      return;
   }
   lzma_end(&stream);


   tgt[0] = 'X';  /* Signature of LZMA from XZ Utils */
   tgt[1] = 'Z';
   tgt[2] = 0;

   in_size   = (unsigned) (*srcsize);
   out_size  = stream.total_out;             /* compressed size */

   tgt[3] = (char)(out_size & 0xff);
   tgt[4] = (char)((out_size >> 8) & 0xff);
   tgt[5] = (char)((out_size >> 16) & 0xff);

   tgt[6] = (char)(in_size & 0xff);         /* decompressed size */
   tgt[7] = (char)((in_size >> 8) & 0xff);
   tgt[8] = (char)((in_size >> 16) & 0xff);

   *irep = (int)stream.total_out + kHeaderSize;
}

void R__unzipLZMA(int *srcsize, unsigned char *src, int *tgtsize, unsigned char *tgt, int *irep)
{
   lzma_stream stream = LZMA_STREAM_INIT;
   lzma_ret returnStatus;

   *irep = 0;

   returnStatus = lzma_stream_decoder(&stream,
                                      UINT64_MAX,
                                      0U);
   if (returnStatus != LZMA_OK) {
      fprintf(stderr,
              "R__unzipLZMA: error %d in lzma_stream_decoder\n",
              returnStatus);
      return;
   }

   stream.next_in   = (const uint8_t *)(&src[kHeaderSize]);
   stream.avail_in  = (size_t)(*srcsize);
   stream.next_out  = (uint8_t *)tgt;
   stream.avail_out = (size_t)(*tgtsize);

   returnStatus = lzma_code(&stream, LZMA_FINISH);
   if (returnStatus != LZMA_STREAM_END) {
      fprintf(stderr,
              "R__unzipLZMA: error %d in lzma_code\n",
              returnStatus);
      lzma_end(&stream);
      return;
   }
   lzma_end(&stream);

   *irep = (int)stream.total_out;
}
