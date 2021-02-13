// @(#)root/base:$Id$
// Author: Fons Rademakers   04/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TString Input/Output functions, put here so the linker will include  //
// them only if I/O is done.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cctype>              // Looking for isspace()
#include <iostream>

#include "TString.h"


////////////////////////////////////////////////////////////////////////////////
/// Replace string with the contents of strm, stopping at an EOF.

std::istream& TString::ReadFile(std::istream& strm)
{
   // get file size
   Ssiz_t end, cur = strm.tellg();
   strm.seekg(0, std::ios::end);
   end = strm.tellg();
   strm.seekg(cur);

   // any positive number of reasonable size for a file
   const Ssiz_t incr = 256;

   Clobber(end-cur);

   while(1) {
      Ssiz_t len = Length();
      Ssiz_t cap = Capacity();
      strm.read(GetPointer()+len, cap-len);
      SetSize(len + strm.gcount());

      if (!strm.good())
         break;                    // EOF encountered

      // If we got here, the read must have stopped because
      // the buffer was going to overflow. Resize and keep
      // going.
      cap = AdjustCapacity(cap, cap+incr);
      Capacity(cap);
   }

   GetPointer()[Length()] = '\0';         // Add null terminator

   return strm;
}

////////////////////////////////////////////////////////////////////////////////
/// Read a line from stream upto newline skipping any whitespace.

std::istream& TString::ReadLine(std::istream& strm, Bool_t skipWhite)
{
   if (skipWhite)
      strm >> std::ws;

   return ReadToDelim(strm, '\n');
}

////////////////////////////////////////////////////////////////////////////////
/// Read a line from stream upto \0, including any newline.

std::istream& TString::ReadString(std::istream& strm)
{
   return ReadToDelim(strm, '\0');
}

////////////////////////////////////////////////////////////////////////////////
/// Read up to an EOF, or a delimiting character, whichever comes
/// first.  The delimiter is not stored in the string,
/// but is removed from the input stream.
/// Because we don't know how big a string to expect, we first read
/// as much as we can and then, if the EOF or null hasn't been
/// encountered, do a resize and keep reading.

std::istream& TString::ReadToDelim(std::istream& strm, char delim)
{
   // any positive number of reasonable size for a string
   const Ssiz_t incr = 32;

   Clobber(incr);

   int p = strm.peek();             // Check if we are already at delim
   if (p == delim) {
      strm.get();                    // eat the delimiter, and return \0.
   } else {
      while (1) {
         Ssiz_t len = Length();
         Ssiz_t cap = Capacity();
         strm.get(GetPointer()+len,          // Address of next byte
                  cap-len+1,                 // Space available (+1 for terminator)
                  delim);                    // Delimiter
         SetSize(len + strm.gcount());
         if (!strm.good()) break;            // Check for EOF or stream failure
         p = strm.peek();
         if (p == delim) {                   // Check for delimiter
            strm.get();                      // eat the delimiter.
            break;
         }
         // Delimiter not seen.  Resize and keep going:
         cap = AdjustCapacity(cap, cap+incr);
         Capacity(cap);
      }
   }

   GetPointer()[Length()] = '\0';                // Add null terminator

   return strm;
}

////////////////////////////////////////////////////////////////////////////////
/// Read a token, delimited by whitespace, from the input stream.

std::istream& TString::ReadToken(std::istream& strm)
{
   // any positive number of reasonable size for a token
   const Ssiz_t incr = 16;

   Clobber(incr);

   strm >> std::ws;                                   // Eat whitespace

   UInt_t wid = strm.width(0);
   char c='\0';
   Int_t hitSpace = 0;
   while ((wid == 0 || Length() < (Int_t)wid) &&
          strm.get(c).good() && (hitSpace = isspace((Int_t)c)) == 0) {
      // Check for overflow:
      Ssiz_t len = Length();
      Ssiz_t cap = Capacity();
      if (len == cap) {
         cap = AdjustCapacity(cap, cap+incr);
         Capacity(cap);
      }
      GetPointer()[len] = c;
      len++;
      SetSize(len);
   }
   if (hitSpace)
      strm.putback(c);

   GetPointer()[Length()] = '\0';                       // Add null terminator

   return strm;
}

////////////////////////////////////////////////////////////////////////////////
/// Read string from stream.

std::istream& operator>>(std::istream& strm, TString& s)
{
   return s.ReadToken(strm);
}

////////////////////////////////////////////////////////////////////////////////
/// Write string to stream.

std::ostream& operator<<(std::ostream& os, const TString& s)
{
   if (os.good()) {
      if (os.tie()) os.tie()->flush(); // instead of opfx
      UInt_t len = s.Length();
      UInt_t wid = os.width();
      wid = (len < wid) ? wid - len : 0;
      os.width(wid);
      long flags = os.flags();
      if (wid && !(flags & std::ios::left))
         os << "";  // let the std::ostream fill
      os.write((char*)s.Data(), s.Length());
      if (wid && (flags & std::ios::left))
         os << "";  // let the std::ostream fill
   }
   // instead of os.osfx();
   if (os.flags() & std::ios::unitbuf)
      os.flush();
   return os;
}

// ------------------- C I/O ------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Read one line from the stream, including the \n, or until EOF.
/// Remove the trailing [\r]\n if chop is true. Returns kTRUE if data was read.

Bool_t TString::Gets(FILE *fp, Bool_t chop)
{
   char buf[256];
   Bool_t r = kFALSE;

   Clobber(256);

   do {
      if (fgets(buf, sizeof(buf), fp) == 0) break;
      *this += buf;
      r = kTRUE;
   } while (!ferror(fp) && !feof(fp) && strchr(buf,'\n') == 0);

   if (chop && EndsWith("\n")) {
      Chop();
      if (EndsWith("\r"))
         Chop();
   }

   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Write string to the stream.

void TString::Puts(FILE *fp)
{
   fputs(GetPointer(), fp);
}
