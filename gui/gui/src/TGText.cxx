// @(#)root/gui:$Id: ba5caabd5d69c640536a71daaa6968de966be4a8 $
// Author: Fons Rademakers   26/04/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**************************************************************************

    This source is based on Xclass95, a Win95-looking GUI toolkit.
    Copyright (C) 1996, 1997 David Barth, Ricky Ralston, Hector Peraza.

    Xclass95 is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

**************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGText                                                               //
//                                                                      //
// A TGText is a multi line text buffer. It allows the text to be       //
// loaded from file, saved to file and edited. It is used in the        //
// TGTextEdit widget. Single line text is handled by TGTextBuffer       //
// and the TGTextEntry widget.                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGText.h"
#include "strlcpy.h"
#include <cctype>

const Int_t kMaxLen = 8000;


ClassImp(TGTextLine);

////////////////////////////////////////////////////////////////////////////////
/// Create empty line of text (default ctor).

TGTextLine::TGTextLine()
{
   fLength = 0;
   fString = 0;
   fPrev = fNext = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize line of text with other line of text (not copy ctor).

TGTextLine::TGTextLine(TGTextLine *line)
{
   fLength = line->GetLineLength();
   fString = 0;
   if (fLength > 0)
      fString = line->GetText(0, line->GetLineLength());
   fPrev = fNext = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize line of text with a const char*.

TGTextLine::TGTextLine(const char *string)
{
   if (string) {
      fLength = strlen(string);
      fString = new char[fLength+1];
      strncpy(fString, string, fLength);
      fString[fLength] = 0;
   } else {
      fLength = 0;
      fString = 0;
   }
   fPrev = fNext = 0;
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGTextLine::TGTextLine(const TGTextLine& tl) : fLength(tl.fLength),
   fPrev(tl.fPrev), fNext(tl.fNext)
{
   fString = 0;
   if (tl.fString) {
      fString = new char[fLength+1];
      strncpy(fString, tl.fString, fLength);
      fString[fLength] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGTextLine& TGTextLine::operator=(const TGTextLine& tl)
{
   if (this != &tl) {
      fLength = tl.fLength;
      if (fString) delete [] fString;
      fString = new char[fLength+1];
      strncpy(fString, tl.fString, fLength);
      fString[fLength] = 0;
      fPrev = tl.fPrev;
      fNext = tl.fNext;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a line of text.

TGTextLine::~TGTextLine()
{
   if (fString)
      delete [] fString;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear a line of text.

void TGTextLine::Clear()
{
   if (fString)
      delete [] fString;
   fString = 0;
   fLength = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete length chars from line starting at position pos.

void TGTextLine::DelText(ULong_t pos, ULong_t length)
{
   if (fLength == 0 || pos >= fLength)
      return;
   if (pos+length > fLength)
      length = fLength - pos;

   if (fLength - length <= 0) {
      delete [] fString;
      fLength = 0;
      fString = 0;
      return;
   }
   char *newstring = new char[fLength - length+1];
   strncpy(newstring, fString, (UInt_t)pos);
   strncpy(newstring+pos, fString+pos+length, UInt_t(fLength-pos-length));
   delete [] fString;
   fString = newstring;
   fLength = fLength - length;
   fString[fLength] = '\0';
}

////////////////////////////////////////////////////////////////////////////////
/// Insert text in line starting at position pos.

void TGTextLine::InsText(ULong_t pos, const char *text)
{
   if (pos > fLength || !text)
      return;

   char *newstring = new char[strlen(text)+fLength+1];
   if (fString != 0)
      strncpy(newstring, fString, (UInt_t)pos);
   // coverity[secure_coding]
   strcpy(newstring+pos, text);
   if (fString != 0 && fLength - pos  > 0)
      strncpy(newstring+pos+strlen(text), fString+pos, UInt_t(fLength-pos));
   fLength = fLength + strlen(text);
   delete [] fString;
   fString = newstring;
   fString[fLength] ='\0';
}

////////////////////////////////////////////////////////////////////////////////
/// Get length characters from line starting at pos. Returns 0
/// in case pos and length are out of range. The returned string
/// must be freed by the user.

char *TGTextLine::GetText(ULong_t pos, ULong_t length)
{
   if (pos >= fLength) {
      return 0;
   }

   if (pos + length > (ULong_t)fString) {
      length = fLength - pos;
   }

   char *retstring = new char[length+1];
   retstring[length] = '\0';
   strncpy(retstring, fString+pos, (UInt_t)length);

   return retstring;
}

////////////////////////////////////////////////////////////////////////////////
/// Get word at position. Returned string must be deleted.

char *TGTextLine::GetWord(ULong_t pos)
{
   if (pos >= fLength) {
      return 0;
   }

   Int_t start = (Int_t)pos;
   UInt_t end = (UInt_t)pos;
   UInt_t i = (UInt_t)pos;

   if (fString[i] == ' ' || fString[i] == '\t') {
      while (start >= 0) {
         if (fString[start] == ' ' || fString[start] == '\t') --start;
         else break;
      }
      ++start;
      while (end < fLength) {
         if (fString[end] == ' ' || fString[end] == '\t') ++end;
         else break;
      }
   } else if (isalnum(fString[i])) {
      while (start >= 0) {
         if (isalnum(fString[start])) --start;
         else break;
      }
      ++start;
      while (end < fLength) {
         if (isalnum(fString[end])) ++end;
         else break;
      }
   } else {
      while (start >= 0) {
         if (isalnum(fString[start]) || fString[start] == ' ' || fString[start] == '\t') {
            break;
         } else {
            --start;
         }
      }
      ++start;
      while (end < fLength) {
         if (isalnum(fString[end]) || fString[end] == ' ' || fString[end] == '\t') {
            break;
         } else {
            ++end;
         }
      }
   }

   UInt_t length = UInt_t(end - start);
   char *retstring = new char[length+1];
   retstring[length] = '\0';
   strncpy(retstring, fString+start, length);

   return retstring;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a character from the line.

void TGTextLine::DelChar(ULong_t pos)
{
   char *newstring;
   if ((fLength <= 0) || (pos > fLength))
      return;
   newstring = new char[fLength];
   strncpy(newstring, fString, (UInt_t)pos-1);
   if (pos < fLength)
      strncpy(newstring+pos-1, fString+pos, UInt_t(fLength-pos+1));
   else
      newstring[pos-1] = 0;
   delete [] fString;
   fString = newstring;
   fLength--;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert a character at the specified position.

void TGTextLine::InsChar(ULong_t pos, char character)
{
   char *newstring;
   if (pos > fLength)
      return;
   newstring = new char[fLength+2];
   newstring[fLength+1] = '\0';
   if (fLength > 0)
      strncpy (newstring, fString, (UInt_t)pos);
   newstring[pos] = character;
   if (fLength - pos > 0)
      strncpy(newstring+pos+1, fString+pos, UInt_t(fLength-pos));
   delete [] fString;
   fString = newstring;
   fLength++;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a character at the specified position from the line.
/// Returns -1 if pos is out of range.

char TGTextLine::GetChar(ULong_t pos)
{
   if ((fLength <= 0) || (pos >= fLength))
      return -1;
   return fString[pos];
}


ClassImp(TGText);

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGText::TGText(const TGText& gt) :
  fFilename(gt.fFilename),
  fIsSaved(gt.fIsSaved),
  fFirst(gt.fFirst),
  fCurrent(gt.fCurrent),
  fCurrentRow(gt.fCurrentRow),
  fRowCount(gt.fRowCount),
  fColCount(gt.fColCount),
  fLongestLine(gt.fLongestLine)
{
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGText& TGText::operator=(const TGText& gt)
{
   if(this!=&gt) {
      fFilename=gt.fFilename;
      fIsSaved=gt.fIsSaved;
      fFirst=gt.fFirst;
      fCurrent=gt.fCurrent;
      fCurrentRow=gt.fCurrentRow;
      fRowCount=gt.fRowCount;
      fColCount=gt.fColCount;
      fLongestLine=gt.fLongestLine;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Common initialization method.

void TGText::Init()
{
   fFirst       = new TGTextLine;
   fCurrent     = fFirst;
   fCurrentRow  = 0;
   fColCount    = 0;
   fRowCount    = 1;
   fLongestLine = 0;
   fIsSaved     = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create default (empty) text buffer.

TGText::TGText()
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create text buffer and initialize with other text buffer.

TGText::TGText(TGText *text)
{
   TGLongPosition pos, end;

   pos.fX = pos.fY = 0;
   end.fY = text->RowCount() - 1;
   end.fX = text->GetLineLength(end.fY) - 1;
   Init();
   InsText(pos, text, pos, end);
}

////////////////////////////////////////////////////////////////////////////////
/// Create text buffer and initialize with single line string.

TGText::TGText(const char *string)
{
   TGLongPosition pos;

   pos.fX = pos.fY = 0;
   Init();
   InsText(pos, string);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy text buffer.

TGText::~TGText()
{
   Clear();
   delete fFirst;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear text buffer.

void TGText::Clear()
{
   TGTextLine *travel = fFirst->fNext;
   TGTextLine *toDelete;
   while (travel != 0) {
      toDelete = travel;
      travel = travel->fNext;
      delete toDelete;
   }
   fFirst->Clear();
   fFirst->fNext = 0;
   fCurrent      = fFirst;
   fCurrentRow   = 0;
   fColCount     = 0;
   fRowCount     = 1;
   fLongestLine  = 0;
   fIsSaved      = kTRUE;
   fFilename     = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Load text from file fn. Startpos is the begin from where to
/// load the file and length is the number of characters to read
/// from the file.

Bool_t TGText::Load(const char *fn, Long_t startpos, Long_t length)
{
   Bool_t      isFirst = kTRUE;
   Bool_t      finished = kFALSE;
   Long_t      count, charcount, i, cnt;
   FILE       *fp;
   char       *buf, c, *src, *dst, *buffer, *buf2;
   TGTextLine *travel, *temp;

   travel = fFirst;

   if (!(fp = fopen(fn, "r"))) return kFALSE;
   buf = new char[kMaxLen];
   i = 0;
   fseek(fp, startpos, SEEK_SET);
   charcount = 0;
   while (fgets(buf, kMaxLen, fp)) {
      if ((length != -1) && (charcount+(Int_t)strlen(buf) > length)) {
         count = length - charcount;
         finished = kTRUE;
      } else
         count = kMaxLen;
      charcount += strlen(buf);
      buf2 = new char[count+1];
      buf2[count] = '\0';
      src = buf;
      dst = buf2;
      cnt = 0;
      while ((c = *src++)) {
         // Don't put CR or NL in buffer
         if (c == 0x0D || c == 0x0A)
            break;
         // Expand tabs
         else if (c == 0x09) {
            *dst++ = '\t';
            while (((dst-buf2) & 0x7) && (cnt++ < count-1))
               *dst++ = 16;  //*dst++ = ' ';
         } else
            *dst++ = c;
         if (cnt++ >= count-1) break;
      }
      *dst = '\0';
      temp = new TGTextLine;
      const size_t bufferSize = strlen(buf2)+1;
      buffer = new char[bufferSize];
      strlcpy(buffer, buf2, bufferSize);
      temp->fLength = strlen(buf2);
      temp->fString = buffer;
      temp->fNext = temp->fPrev = 0;
      if (isFirst) {
         delete fFirst;
         fFirst   = temp;
         fCurrent = temp;
         travel   = fFirst;
         isFirst  = kFALSE;
      } else {
         travel->fNext = temp;
         temp->fPrev   = travel;
         travel        = travel->fNext;
      }
      ++i;
      delete [] buf2;
      if (finished)
         break;
   }
   fclose(fp);
   delete [] buf;

   // Remember the number of lines
   fRowCount = i;
   if (fRowCount == 0)
      fRowCount++;
   fIsSaved  = kTRUE;
   fFilename = fn;
   LongestLine();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Load a 0 terminated buffer. Lines will be split at '\n'.

Bool_t TGText::LoadBuffer(const char *txtbuf)
{
   Bool_t      isFirst = kTRUE;
   Bool_t      finished = kFALSE, lastnl = kFALSE;
   Long_t      i, cnt;
   TGTextLine *travel, *temp;
   char       *buf, c, *src, *dst, *buffer, *buf2, *s;
   const char *tbuf = txtbuf;

   travel = fFirst;

   if (!tbuf || !tbuf[0])
      return kFALSE;

   buf = new char[kMaxLen];
   i = 0;
next:
   if ((s = (char*)strchr(tbuf, '\n'))) {
      if (s-tbuf+1 >= kMaxLen-1) {
         strncpy(buf, tbuf, kMaxLen-2);
         buf[kMaxLen-2] = '\n';
         buf[kMaxLen-1] = 0;
      } else {
         strncpy(buf, tbuf, s-tbuf+1);
         buf[s-tbuf+1] = 0;
      }
      tbuf = s+1;
   } else {
      strncpy(buf, tbuf, kMaxLen-1);
      buf[kMaxLen-1] = 0;
      finished = kTRUE;
   }

   buf2 = new char[kMaxLen+1];
   buf2[kMaxLen] = '\0';
   src = buf;
   dst = buf2;
   cnt = 0;
   while ((c = *src++)) {
      // Don't put CR or NL in buffer
      if (c == 0x0D || c == 0x0A)
         break;
      // Expand tabs
      else if (c == 0x09) {
         *dst++ = '\t';
         while (((dst-buf2) & 0x7) && (cnt++ < kMaxLen-1))
            *dst++ = 16;  //*dst++ = ' ';
      } else
         *dst++ = c;
      if (cnt++ >= kMaxLen-1) break;
   }
   *dst = '\0';
   temp = new TGTextLine;
   const size_t bufferSize = strlen(buf2) + 1;
   buffer = new char[bufferSize];
   strlcpy(buffer, buf2, bufferSize);
   temp->fLength = strlen(buf2);
   temp->fString = buffer;
   temp->fNext = temp->fPrev = 0;
   if (isFirst) {
      delete fFirst;
      fFirst   = temp;
      fCurrent = temp;
      travel   = fFirst;
      isFirst  = kFALSE;
   } else {
      travel->fNext = temp;
      temp->fPrev   = travel;
      travel        = travel->fNext;
   }
   ++i;
   delete [] buf2;

   // make sure that \n generates a single empty line in the TGText
   if (!lastnl && !*tbuf && *(tbuf-1) == '\n') {
      tbuf--;
      lastnl = kTRUE;
   }

   if (!finished && strlen(tbuf))
      goto next;

   delete [] buf;
   // Remember the number of lines
   fRowCount = i;
   if (fRowCount == 0)
      fRowCount++;
   fIsSaved  = kTRUE;
   fFilename = "";
   LongestLine();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save text buffer to file fn.

Bool_t TGText::Save(const char *fn)
{
   char *buffer;
   TGTextLine *travel = fFirst;
   FILE *fp;
   if (!(fp = fopen(fn, "w"))) return kFALSE;

   while (travel) {
      ULong_t i = 0;
      buffer = new char[travel->fLength+2];
      strncpy(buffer, travel->fString, (UInt_t)travel->fLength);
      buffer[travel->fLength]   = '\n';
      buffer[travel->fLength+1] = '\0';
      while (buffer[i] != '\0') {
         if (buffer[i] == '\t') {
            ULong_t j = i+1;
            while (buffer[j] == 16)
               j++;
            // coverity[secure_coding]
            strcpy(buffer+i+1, buffer+j);
         }
         i++;
      }
      if (fputs(buffer, fp) == EOF) {
         delete [] buffer;
         fclose(fp);
         return kFALSE;
      }
      delete [] buffer;
      travel = travel->fNext;
   }
   fIsSaved = kTRUE;
   fFilename = fn;
   fclose(fp);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Append buffer to file fn.

Bool_t TGText::Append(const char *fn)
{
   char *buffer;
   TGTextLine *travel = fFirst;
   FILE *fp;
   if (!(fp = fopen(fn, "a"))) return kFALSE;

   while (travel) {
      ULong_t i = 0;
      buffer = new char[travel->fLength+2];
      strncpy(buffer, travel->fString, (UInt_t)travel->fLength);
      buffer[travel->fLength]   = '\n';
      buffer[travel->fLength+1] = '\0';
      while (buffer[i] != '\0') {
         if (buffer[i] == '\t') {
            ULong_t j = i+1;
            while (buffer[j] == 16 && buffer[j] != '\0')
               j++;
            // coverity[secure_coding]
            strcpy(buffer+i+1, buffer+j);
         }
         i++;
      }
      if (fputs(buffer, fp) == EOF) {
         delete [] buffer;
         fclose(fp);
         return kFALSE;
      }
      delete [] buffer;
      travel = travel->fNext;
   }
   fIsSaved = kTRUE;
   fclose(fp);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete character at specified position pos.

Bool_t TGText::DelChar(TGLongPosition pos)
{
   if ((pos.fY >= fRowCount) || (pos.fY < 0))
      return kFALSE;

   if (!SetCurrentRow(pos.fY)) return kFALSE;
   fCurrent->DelChar(pos.fX);

   fIsSaved = kFALSE;
   LongestLine();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert character c at the specified position pos.

Bool_t TGText::InsChar(TGLongPosition pos, char c)
{
   if ((pos.fY >= fRowCount) || (pos.fY < 0) || (pos.fX < 0))
      return kFALSE;

   if (!SetCurrentRow(pos.fY)) return kFALSE;
   fCurrent->InsChar(pos.fX, c);

   fIsSaved = kFALSE;
   LongestLine();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get character a position pos. If charcater not valid return -1.

char TGText::GetChar(TGLongPosition pos)
{
   if (pos.fY >= fRowCount)
      return -1;

   if (!SetCurrentRow(pos.fY)) return -1;
   return fCurrent->GetChar(pos.fX);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete text between start and end positions. Returns false in
/// case of failure (start and end not being within bounds).

Bool_t TGText::DelText(TGLongPosition start, TGLongPosition end)
{
   if ((start.fY < 0) || (start.fY >= fRowCount) ||
       (end.fY < 0)   || (end.fY >= fRowCount)) {
      return kFALSE;
   }

   if ((end.fX < 0) || (end.fX > GetLineLength(end.fY))) {
      return kFALSE;
   }

   char *tempbuffer;

   if (!SetCurrentRow(start.fY)) return kFALSE;

   if (start.fY == end.fY) {
      fCurrent->DelText(start.fX, end.fX-start.fX+1);
      return kTRUE;
   }
   fCurrent->DelText(start.fX, fCurrent->fLength-start.fX);
   SetCurrentRow(fCurrentRow+1);
   for (Long_t i = start.fY+1; i < end.fY; i++) {
      DelLine(fCurrentRow);
   }

   tempbuffer = fCurrent->GetText(end.fX+1, fCurrent->fLength-end.fX-1);
   DelLine(fCurrentRow);
   SetCurrentRow(start.fY);
   if (tempbuffer) {
      fCurrent->InsText(fCurrent->GetLineLength(), tempbuffer);
      delete [] tempbuffer;
   } else {
      if (fCurrent->fNext) {
         fCurrent->InsText(fCurrent->fLength, fCurrent->fNext->fString);
         DelLine(fCurrentRow+1);
         SetCurrentRow(start.fY);
      }
   }

   fIsSaved = kFALSE;
   LongestLine();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert src text from start_src to end_src into text at position ins_pos.
/// Returns false in case of failure (start_src, end_src out of range for
/// src, and ins_pos out of range for this).

Bool_t TGText::InsText(TGLongPosition ins_pos, TGText *src,
                       TGLongPosition start_src, TGLongPosition end_src)
{
   /*
   if ((start_src.fY < 0) || (start_src.fY >= src->RowCount()) ||
       (end_src.fY < 0)   || (end_src.fY >= src->RowCount()))
      return kFALSE;
   if ((start_src.fX < 0) || (start_src.fX > src->GetLineLength(start_src.fY)) ||
       (end_src.fX < 0)   || (end_src.fX > src->GetLineLength(end_src.fY)))
      return kFALSE;
   if ((ins_pos.fY < 0) || (ins_pos.fY > fRowCount))
      return kFALSE;
   if ((ins_pos.fX < 0) || (ins_pos.fX > GetLineLength(ins_pos.fY)))
      return kFALSE;
   */
   if (ins_pos.fY > fRowCount)
      return kFALSE;

   TGLongPosition pos;
   ULong_t len;
   char *lineString;
   char *restString;
   TGTextLine *following;

   if (ins_pos.fY == fRowCount) {  // for appending text
      pos.fY = fRowCount - 1;
      pos.fX = GetLineLength(pos.fY);
      BreakLine(pos);  // current row is set by this
   } else {
      // otherwise going to the desired row
      if (!SetCurrentRow(ins_pos.fY)) return kFALSE;
   }

   // preparing first line to be inserted
   restString = fCurrent->GetText(ins_pos.fX, fCurrent->fLength - ins_pos.fX);
   fCurrent->DelText(ins_pos.fX, fCurrent->fLength - ins_pos.fX);
   following = fCurrent->fNext;
   // inserting first line
   if (start_src.fY == end_src.fY) {
      len = end_src.fX - start_src.fX+1;
   } else {
      len = src->GetLineLength(start_src.fY) - start_src.fX;
   }

   if (len > 0) {
      lineString = src->GetLine(start_src, len);
      fCurrent->InsText(ins_pos.fX, lineString);
      delete [] lineString;
   }
   // [...] inserting possible lines
   pos.fY = start_src.fY+1;
   pos.fX = 0;
   for ( ; pos.fY < end_src.fY; pos.fY++) {
      Int_t llen = src->GetLineLength(pos.fY);
      lineString = src->GetLine(pos, llen > 0 ? llen : 0);
      fCurrent->fNext = new TGTextLine(lineString);
      fCurrent->fNext->fPrev = fCurrent;
      fCurrent = fCurrent->fNext;
      fRowCount++;
      fCurrentRow++;
      delete [] lineString;
   }
   // last line of inserted text is as special as first line
   if (start_src.fY != end_src.fY) {
      pos.fY = end_src.fY;
      pos.fX = 0;
      lineString = src->GetLine(pos, end_src.fX+1);
      fCurrent->fNext = new TGTextLine(lineString);
      fCurrent->fNext->fPrev = fCurrent;
      fCurrent = fCurrent->fNext;
      fRowCount++;
      fCurrentRow++;
      delete [] lineString;
   }
   // ok, now we have to add the rest of the first destination line
   if (restString) {
#if 0
      if (ins_pos.fX == 0) {
         fCurrent->fNext = new TGTextLine(restString);
         fCurrent->fNext->fPrev = fCurrent;
         fCurrent = fCurrent->fNext;
         fRowCount++;
         fCurrentRow++;
      } else
#endif
         fCurrent->InsText(fCurrent->fLength, restString);
      delete [] restString;
   }
   // now re-linking the rest of the origin text
   fCurrent->fNext = following;
   if (fCurrent->fNext) {
      fCurrent->fNext->fPrev = fCurrent;
   }

   LongestLine();
   fIsSaved = kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert single line at specified position. Return false in case position
/// is out of bounds.

Bool_t TGText::InsText(TGLongPosition pos, const char *buffer)
{
   if (pos.fY < 0 || pos.fY > fRowCount) {
      return kFALSE;
   }

   if (pos.fY == fRowCount) {
      SetCurrentRow(fRowCount-1);
      fCurrent->fNext = new TGTextLine(buffer);
      fCurrent->fNext->fPrev = fCurrent;
      fRowCount++;
   } else {
      SetCurrentRow(pos.fY);
      fCurrent->InsText(pos.fX, buffer);
   }
   LongestLine();
   fIsSaved = kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add another text buffer to this buffer.

Bool_t TGText::AddText(TGText *text)
{
   TGLongPosition end, start_src, end_src;

   end.fY = fRowCount;
   end.fX = 0;
   start_src.fX = start_src.fY = 0;
   end_src.fY   = text->RowCount()-1;
   end_src.fX   = text->GetLineLength(end_src.fY)-1;
   fIsSaved     = kFALSE;
   return InsText(end, text, start_src, end_src);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert string before specified position.
/// Returns false if insertion failed.

Bool_t TGText::InsLine(ULong_t pos, const char *string)
{
   TGTextLine *previous, *newline;
   if ((Long_t)pos > fRowCount) {
      return kFALSE;
   }
   if ((Long_t)pos < fRowCount) {
      SetCurrentRow(pos);
   } else {
      SetCurrentRow(fRowCount-1);
   }

   if (!fCurrent) return kFALSE;

   previous = fCurrent->fPrev;
   newline = new TGTextLine(string);
   newline->fPrev = previous;
   if (previous) {
      previous->fNext = newline;
   } else {
      fFirst = newline;
   }

   newline->fNext = fCurrent;
   fCurrent->fPrev = newline;
   fRowCount++;
   fCurrentRow++;

   LongestLine();
   fIsSaved = kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete specified row. Returns false if row does not exist.

Bool_t TGText::DelLine(ULong_t pos)
{
   if (!SetCurrentRow(pos) || (fRowCount == 1)) {
      return kFALSE;
   }

   TGTextLine *travel = fCurrent;
   if (travel == fFirst) {
      fFirst = fFirst->fNext;
      fFirst->fPrev = 0;
   } else {
      travel->fPrev->fNext = travel->fNext;
      if (travel->fNext) {
         travel->fNext->fPrev = travel->fPrev;
         fCurrent = fCurrent->fNext;
      } else {
         fCurrent = fCurrent->fPrev;
         fCurrentRow--;
      }
   }
   delete travel;

   fRowCount--;
   fIsSaved = kFALSE;
   LongestLine();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return string at position pos. Returns 0 in case pos is not valid.
/// The returned string must be deleted by the user.

char *TGText::GetLine(TGLongPosition pos, ULong_t length)
{
   if (SetCurrentRow(pos.fY)) {
      return fCurrent->GetText(pos.fX, length);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Break line at position pos. Returns false if pos is not valid.

Bool_t TGText::BreakLine(TGLongPosition pos)
{
   if (!SetCurrentRow(pos.fY))
      return kFALSE;
   if ((pos.fX < 0) || (pos.fX > (Long_t)fCurrent->fLength))
      return kFALSE;

   TGTextLine *temp;
   char *tempbuffer;
   if (pos.fX < (Long_t)fCurrent->fLength) {
      tempbuffer = fCurrent->GetText(pos.fX, fCurrent->fLength-pos.fX);
      temp = new TGTextLine(tempbuffer);
      fCurrent->DelText(pos.fX, fCurrent->fLength - pos.fX);
      delete [] tempbuffer;
   } else {
      temp = new TGTextLine;
   }
   temp->fPrev = fCurrent;
   temp->fNext = fCurrent->fNext;
   fCurrent->fNext = temp;
   if (temp->fNext) {
      temp->fNext->fPrev = temp;
   }

   fIsSaved = kFALSE;
   fRowCount++;
   fCurrentRow++;
   fCurrent = fCurrent->fNext;
   LongestLine();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get length of specified line. Returns -1 if row does not exist.

Long_t TGText::GetLineLength(Long_t row)
{
   if (!SetCurrentRow(row)) {
      return -1;
   }
   return (Long_t)fCurrent->GetLineLength();
}

////////////////////////////////////////////////////////////////////////////////
/// Make specified row the current row. Returns false if row does not exist.
/// In which case fCurrent is not changed or set to the last valid line.

Bool_t TGText::SetCurrentRow(Long_t row)
{
   Long_t count;
   if ((row < 0) || (row >= fRowCount)) {
      return kFALSE;
   }
   if (row > fCurrentRow) {
      for (count = fCurrentRow; count < row; count++) {
         if (!fCurrent->fNext) {
            fCurrentRow = count;
            return kFALSE;
         }
         fCurrent = fCurrent->fNext;
      }
   } else {
      if (fCurrentRow == row)
         return kTRUE;
      for (count = fCurrentRow; count > row; count--) {
         if (!fCurrent->fPrev) {
            fCurrentRow = count;
            return kFALSE;
         }
         fCurrent = fCurrent->fPrev;
      }
   }
   fCurrentRow = row;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Redo all tabs in a line. Needed after a new tab is inserted.

void TGText::ReTab(Long_t row)
{
   if (!SetCurrentRow(row)) {
      return;
   }

   // first remove all special tab characters (16)
   char *buffer;
   ULong_t i = 0;

   buffer = fCurrent->fString;
   while (buffer[i] != '\0') {
      if (buffer[i] == '\t') {
         ULong_t j = i+1;
         while (buffer[j] == 16 && buffer[j] != '\0') {
            j++;
         }
         // coverity[secure_coding]
         strcpy(buffer+i+1, buffer+j);
      }
      i++;
   }

   char   c, *src, *dst, *buf2;
   Long_t cnt;

   buf2 = new char[kMaxLen+1];
   buf2[kMaxLen] = '\0';
   src = buffer;
   dst = buf2;
   cnt = 0;
   while ((c = *src++)) {
      // Expand tabs
      if (c == 0x09) {
         *dst++ = '\t';
         while (((dst-buf2) & 0x7) && (cnt++ < kMaxLen-1)) {
            *dst++ = 16;
         }
      } else {
         *dst++ = c;
      }
      if (cnt++ >= kMaxLen-1) break;
   }
   *dst = '\0';

   fCurrent->fString = buf2;
   fCurrent->fLength = strlen(buf2);

   delete [] buffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Search for string searchString starting at the specified position going
/// in forward (direction = true) or backward direction. Returns true if
/// found and foundPos is set accordingly.

Bool_t TGText::Search(TGLongPosition *foundPos, TGLongPosition start,
                      const char *searchString,
                      Bool_t direction, Bool_t caseSensitive)
{
   if (!SetCurrentRow(start.fY))
      return kFALSE;

   Ssiz_t x = kNPOS;

   if (direction) {
      while(1) {
         TString s = fCurrent->fString;
         x = s.Index(searchString, (Ssiz_t)start.fX,
                     caseSensitive ? TString::kExact : TString::kIgnoreCase);
         if (x != kNPOS) {
            foundPos->fX = x;
            foundPos->fY = fCurrentRow;
            return kTRUE;
         }
         if (!SetCurrentRow(fCurrentRow+1))
            break;
         start.fX = 0;
      }
   } else {
      while(1) {
         TString s = fCurrent->fString;
         for (int i = (int)start.fX; i >= 0; i--) {
            x = s.Index(searchString, (Ssiz_t)i,
                        caseSensitive ? TString::kExact : TString::kIgnoreCase);
            if (x >= start.fX) {
               x = kNPOS;
               continue;
            }
            if (x != kNPOS) {
               break;
            }
         }
         if (x != kNPOS) {
            foundPos->fX = x;
            foundPos->fY = fCurrentRow;
            return kTRUE;
         }
         if (!SetCurrentRow(fCurrentRow-1)) {
            break;
         }
         start.fX = fCurrent->fLength;
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace oldText by newText. Returns false if nothing replaced.

Bool_t TGText::Replace(TGLongPosition start, const char *oldText, const char *newText,
                       Bool_t direction, Bool_t caseSensitive)
{
   if (!SetCurrentRow(start.fY)) {
      return kFALSE;
   }

   TGLongPosition foundPos;
   if (!Search(&foundPos, start, oldText, direction, caseSensitive)) {
      return kFALSE;
   }

   TGLongPosition delEnd;
   delEnd.fY = foundPos.fY;
   delEnd.fX = foundPos.fX + strlen(oldText) - 1;
   DelText(foundPos, delEnd);
   InsText(foundPos, newText);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set fLongestLine.

void TGText::LongestLine()
{
   Long_t line_count = 0;
   TGTextLine *travel = fFirst;
   fColCount = 0;
   while (travel) {
      if ((Long_t)travel->fLength > fColCount) {
         fColCount = travel->fLength;
         fLongestLine = line_count;
      }
      travel = travel->fNext;
      line_count++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns content as ROOT string

TString TGText::AsString()
{
   TString ret;

   Long_t line_count = 0;
   TGTextLine *travel = fFirst;
   fColCount = 0;

   while (travel) {
      if ((Long_t)travel->fLength > fColCount) {
         fColCount = travel->fLength;
         fLongestLine = line_count;
      }
      char *txt = travel->GetText();
      ret += txt;
      delete [] txt;
      travel = travel->fNext;
      if (travel) ret += '\n';
      line_count++;
   }

   return ret;
}
