// Author: Brian Bockelman March 2015

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_ZipLZ4
#define ROOT_ZipLZ4

// NOTE: the ROOT compression libraries aren't consistently written in C++; hence the
// #ifdef's to avoid problems with C code.
#ifdef __cplusplus
extern "C" {
#endif
void R__zipLZ4(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep);
void R__unzipLZ4(int *srcsize, unsigned char *src, int *tgtsize, unsigned char *tgt, int *irep);
#ifdef __cplusplus
}
#endif

#endif
