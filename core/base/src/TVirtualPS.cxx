// @(#)root/base:$Id$
// Author: Rene Brun   05/09/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVirtualPS
\ingroup Base
\ingroup PS

TVirtualPS is an abstract interface to Postscript, PDF, SVG. TeX etc... drivers
*/

#include <fstream>
#include "strlcpy.h"
#include "TVirtualPS.h"

TVirtualPS *gVirtualPS = nullptr;

const Int_t  kMaxBuffer = 250;

ClassImp(TVirtualPS);


////////////////////////////////////////////////////////////////////////////////
/// VirtualPS default constructor.

TVirtualPS::TVirtualPS()
{
   fStream    = nullptr;
   fNByte     = 0;
   fSizBuffer = kMaxBuffer;
   fBuffer    = new char[fSizBuffer+1];
   fLenBuffer = 0;
   fPrinted   = kFALSE;
   fImplicitCREsc = nullptr;
}


////////////////////////////////////////////////////////////////////////////////
/// VirtualPS constructor.

TVirtualPS::TVirtualPS(const char *name, Int_t)
          : TNamed(name,"Postscript interface")
{
   fStream    = nullptr;
   fNByte     = 0;
   fSizBuffer = kMaxBuffer;
   fBuffer    = new char[fSizBuffer+1];
   fLenBuffer = 0;
   fPrinted   = kFALSE;
   fImplicitCREsc = nullptr;
}


////////////////////////////////////////////////////////////////////////////////
/// VirtualPS destructor

TVirtualPS::~TVirtualPS()
{
   if (fBuffer) delete [] fBuffer;
}


////////////////////////////////////////////////////////////////////////////////
/// Output the string str in the output buffer

void TVirtualPS::PrintStr(const char *str)
{
   if (!str || !str[0])
      return;
   Int_t len = strlen(str);
   while (len) {
      if (str[0] == '@') {
         if (fLenBuffer) {
            fStream->write(fBuffer, fLenBuffer);
            fNByte += fLenBuffer;
            fLenBuffer = 0;
            fStream->write("\n", 1);
            fNByte++;
            fPrinted = kTRUE;
         }
         len--;
         str++;
      } else {
         Int_t lenText = len;
         if (str[len-1] == '@') lenText--;
         PrintFast(lenText, str);
         len -= lenText;
         str += lenText;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Fast version of Print

void TVirtualPS::PrintFast(Int_t len, const char *str)
{
   if (!len || !str) return;
   while ((len + fLenBuffer) > kMaxBuffer) {
      Int_t nWrite = kMaxBuffer;
      if (fImplicitCREsc) {
         if (fLenBuffer > 0) nWrite = fLenBuffer;
      } else {
         if ((len + fLenBuffer) > nWrite) {
            // Search for the nearest preceding space to break a line, if there is no instruction to escape the <end-of-line>.
            while ((nWrite >= fLenBuffer) && (str[nWrite - fLenBuffer] != ' ')) nWrite--;
            if (nWrite < fLenBuffer) {
               while ((nWrite >= 0) && (fBuffer[nWrite] != ' ')) nWrite--;
            }
            if (nWrite <= 0) {
               // Cannot find a convenient place to break a line, so we just break at this location.
               nWrite = kMaxBuffer;
            }
         }
      }
      if (nWrite >= fLenBuffer) {
         if (fLenBuffer > 0) {
            fStream->write(fBuffer, fLenBuffer);
            fNByte += fLenBuffer;
            nWrite -= fLenBuffer;
            fLenBuffer = 0;
         }
         if (nWrite > 0) {
            fStream->write(str, nWrite);
            len -= nWrite;
            str += nWrite;
            fNByte += nWrite;
         }
      } else {
         if (nWrite > 0) {
            fStream->write(fBuffer, nWrite);
            fNByte += nWrite;
            memmove(fBuffer, fBuffer + nWrite, fLenBuffer - nWrite); // not strcpy because source and destination overlap
            fBuffer[fLenBuffer - nWrite] = 0; // not sure if this is needed, but just in case
            fLenBuffer -= nWrite;
         }
      }
      if (fImplicitCREsc) {
         // Write escape characters (if any) before an end-of-line is enforced.
         // For example, in PostScript the <new line> character must be escaped inside strings.
         Int_t crlen = strlen(fImplicitCREsc);
         fStream->write(fImplicitCREsc, crlen);
         fNByte += crlen;
      }
      fStream->write("\n",1);
      fNByte++;
   }
   if (len > 0) {
      strlcpy(fBuffer + fLenBuffer, str, len+1);
      fLenBuffer += len;
      fBuffer[fLenBuffer] = 0;
   }
   fPrinted = kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Write one Integer to the file
///
/// n: Integer to be written in the file.
/// space: If TRUE, a space in written before the integer.

void TVirtualPS::WriteInteger(Int_t n, Bool_t space )
{
   char str[15];
   if (space) {
      snprintf(str,15," %d", n);
   } else {
      snprintf(str,15,"%d", n);
   }
   PrintStr(str);
}


////////////////////////////////////////////////////////////////////////////////
/// Write a Real number to the file

void TVirtualPS::WriteReal(Float_t z, Bool_t space)
{
   char str[15];
   if (space) {
      snprintf(str,15," %g", z);
   } else {
      snprintf(str,15,"%g", z);
   }
   PrintStr(str);
}


////////////////////////////////////////////////////////////////////////////////
/// Print a raw

void TVirtualPS::PrintRaw(Int_t len, const char *str)
{
   fNByte += len;
   if ((len + fLenBuffer) > kMaxBuffer - 1) {
      fStream->write(fBuffer, fLenBuffer);
      while(len > kMaxBuffer-1) {
         fStream->write(str,kMaxBuffer);
         len -= kMaxBuffer;
         str += kMaxBuffer;
      }
      memcpy(fBuffer, str, len);
      fLenBuffer = len;
   } else {
      memcpy(fBuffer + fLenBuffer, str, len);
      fLenBuffer += len;
   }
   fPrinted = kTRUE;
}
