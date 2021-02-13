// @(#)root/new:$Id$
// Author: D.Bertini and M.Ivanov   10/08/2000

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//****************************************************************************//
//
//
// MemCheck is used to check the memory in ROOT based applications.
//
// Principe:
//  A memory leak often arises whenever memory allocated through
//  (new, new[]) is never returned via a corresponding (delete, delete[]).
//  Redefining a special version of (operator new, operator delete) will
//  allow the bookkeeping of created pointers to chunks of memory and, at
//  the same time, pointers to the current function (and all of its callers
//  scanning the stack) in order to trace back where the memory is not freed.
//  This specific bookkeeping of pointers will be done by a kind of
//  "memory info" container class TMemHashTable that will be used via
//  the ROOT memory management defined inside the NewDelete.cxx file.
//
//  To activate the memory checker you have to set in the .rootrc file
//  the resource Root.MemCheck to 1 (e.g.: Root.MemCheck: 1) and you
//  have to link with libNew.so (e.g. use root-config --new --libs) or
//  use rootn.exe. When all this is the case you will find at the end
//  of the program execution a file "memcheck.out" in the directory
//  where you started your ROOT program. Alternatively you can set
//  the resource Root.MemCheckFile to the name of a file to which
//  the leak information will be written. The contents of this
//  "memcheck.out" file can be analyzed and transformed into printable
//  text via the memprobe program (in $ROOTSYS/bin).
//
// (c) 2000 : Gesellschaft fuer Schwerionenforschung GmbH
//            Planckstrasse, 1
//            D-64291 Darmstadt
//            Germany
//
// Created 10/08/2000 by: D.Bertini and M.Ivanov.
// Based on ideas from LeakTracer by Erwin Andreasen.
//
// - Updated:
//    Date: 12/02/2001 Adapt script to new GDB 5.0, new glibc2.2.x and gcc 2.96.
//    Date: 23/10/2000 (hash mechanism speeding up the bookkeeping)
//
// - Documentation:
//
//    http://www-hades.gsi.de/~dbertini/mem.html
//
//****************************************************************************//

#ifndef ROOT_MemCheck
#define ROOT_MemCheck

#include "TROOT.h"


class TStackInfo {
public:
   UInt_t      fSize;             //size of the stack
   Int_t       fTotalAllocCount;  //total number of allocation for stack sequence
   Int_t       fTotalAllocSize;   //total size of allocated memory
   Int_t       fAllocCount;       //current number of allocation-deallocation
   Int_t       fAllocSize;        //current allocated size
   TStackInfo *fNextHash;         //index-pointer to the next info for given hash value

public:
   void     Init(Int_t stacksize, void **stackptrs); //initialization
   void     Inc(Int_t memSize);  //increment counters -when memory allocated
   void     Dec(Int_t memSize);  //decrement counters -when memory deallocated
   ULong_t  Hash();
   Int_t    IsEqual(UInt_t size, void **ptr);
   void    *StackAt(UInt_t i);
   TStackInfo *Next();    //index of the next entries

   static ULong_t HashStack(UInt_t size, void **ptr);
};


class TStackTable {
private:
   char         *fTable;      //pointer to the table
   TStackInfo  **fHashTable;  //pointer to the hash table
   Int_t         fSize;       //current size of the table
   Int_t         fHashSize;   //current size of the hash table
   Int_t         fCount;      //number of entries in table
   char         *fNext;       //pointer to the last stack info

   void Expand(Int_t newsize);

public:
   void        Init();
   TStackInfo *AddInfo(Int_t size, void **stackptrs);
   TStackInfo *FindInfo(Int_t size, void **stackptrs);
   Int_t       GetIndex(TStackInfo *info);
   TStackInfo *GetInfo(Int_t index);
   TStackInfo *First() { return (TStackInfo *)fTable; }
};


class TMemInfo {
public:
   void   *fAddress;    //mem address
   size_t  fSize;       //size of the allocated memory
   Int_t   fStackIndex; //index of the stack info
};

class TMemTable {
public:
   Int_t     fAllocCount;    //number of memory allocation blocks
   Int_t     fMemSize;       //total memory allocated size
   Int_t     fTableSize;     //amount of entries in the below array
   Int_t     fFirstFreeSpot; //where is the first free spot in the leaks array?
   TMemInfo *fLeaks;         //leak table
};

class TDeleteTable {
public:
   Int_t     fAllocCount;    //how many memory blocks do we have
   Int_t     fTableSize;     //amount of entries in the below array
   TMemInfo *fLeaks;         //leak table
};


class TMemHashTable {
public:
   static Int_t        fgSize;            //size of hash table
   static TMemTable  **fgLeak;            //pointer to the hash table
   static Int_t        fgAllocCount;      //number of memory allocation blocks
   static TStackTable  fgStackTable;      //table with stack pointers
   static TDeleteTable fgMultDeleteTable; //pointer to the table

   ~TMemHashTable() { if (TROOT::MemCheck()) Dump(); }

   static void  Init();
   static void  RehashLeak(Int_t newSize);             //rehash leak pointers
   static void *AddPointer(size_t size, void *ptr=0);  //add pointer to the table
   static void  FreePointer(void *p);                  //free pointer
   static void  Dump();                                //write leaks to the output file
};



inline void TStackInfo::Inc(Int_t memSize)
{
   fTotalAllocCount += 1;
   fTotalAllocSize  += memSize;
   fAllocCount      += 1;
   fAllocSize       += memSize;
}

inline void TStackInfo::Dec(int memSize)
{
   fAllocCount -= 1;
   fAllocSize  -= memSize;
}

inline ULong_t TStackInfo::Hash()
{
   return HashStack(fSize, (void**)&(this[1]));
}

inline void *TStackInfo::StackAt(UInt_t i)
{
   //return i<fSize ? ((char*)&(this[1]))+i*sizeof(void*):0;
   void **stptr = (void**)&(this[1]);
   return i < fSize ? stptr[i] : 0;
}

inline TStackInfo *TStackInfo::Next()
{
   //return (TStackInfo*)((char*)(&this[1])+fSize*sizeof(void*));
   return (TStackInfo*)((char*)(this)+fSize*sizeof(void*)+sizeof(TStackInfo));
}

#endif
