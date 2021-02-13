//------------------------------------------------------------------------------
// CLING - the C++ LLVM-based InterpreterG :)
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

// RUN: cat %s | %cling 2>&1 | FileCheck %s
extern "C" int printf(const char*,...);

struct TheStruct {
  char m_Storage[100];
};

// CHECK: About to throw "ABC"
printf("About to throw \"ABC\"\n");
TheStruct getStructButThrows() { throw "ABC!"; }
// CHECK: Exception occurred. Recovering...
getStructButThrows()

struct ThrowInConstructor {
  ThrowInConstructor() { throw 1; }
};
// CHECK: About to throw in constructor
printf("About to throw in constructor\n");
// CHECK: Exception occurred. Recovering...
ThrowInConstructor t;
