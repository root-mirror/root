// Author: Stephan Hageboeck, CERN  2 Mar 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

extern "C" {

/// Default options when address sanitizer starts up in ROOT executables.
/// This is relevant when ROOT's build option `asan` is on.
/// These can be overridden / augmented by the ASAN_OPTIONS environment variable.
/// Using ASAN_OPTIONS=help=1 and starting an instrumented ROOT exectuable, available options will be printed.
const char* __asan_default_options() { 

#ifdef ASAN_DETECT_LEAKS
#define DETECT_LEAKS ":detect_leaks=1"
#else
#define DETECT_LEAKS ":detect_leaks=0"
#endif

  return "strict_string_checks=1"
         ":detect_stack_use_after_return=1"
         ":check_initialization_order=1"
         ":detect_container_overflow=1"
         ":alloc_dealloc_mismatch=0"
         DETECT_LEAKS
         ":verify_asan_link_order=0";
}

/// Default options when leak sanitizer starts up in ROOT exectuables.
/// This is relevant when ROOT's build options `asan` is on.
/// These can be overridden / augmented by the LSAN_OPTIONS environment variable.
/// Using LSAN_OPTIONS=help=1 and starting an instrumented ROOT exectuable, available options will be printed.
const char* __lsan_default_options() {
   return "exitcode=0:max_leaks=10:print_suppressions=1";
}

/// Default suppressions for leak sanitizer in ROOT.
/// Since llvm uses allocators that do not give back memory, many leaks would show up.
/// A customisable version can be found in `$ROOTSYS/build/`.
const char* __lsan_default_suppressions() {
   return "leak:llvm::SmallVectorBase::grow_pod \n"
          "leak:llvm::BumpPtrAllocatorImpl \n"
          "leak:llvm::DenseMap*grow \n"
          "leak:llvm::StringMapImpl \n"
          "leak:llvm::TinyPtr \n"
          "leak:llvm::FoldingSetBase \n"
          "leak:llvm::MemoryBuffer \n"
          "leak:llvm::CodeGenDAG \n"
          "leak:llvm::EmitFastISel \n"
          "leak:llvm-tblgen \n"
          "leak:clang::FileManager::getFile \n" // clang
          "leak:clang::LineTableInfo \n"
          "leak:clang::HeaderSearch \n"
          "leak:clang::Diag \n"
          "leak:clang::Preprocessor:: \n"
          "leak:clang::TextDiagnosticPrinter \n"
          "leak:clang-tblgen \n"
          "leak:cling::IncrementalExecutor\n" //cling macro execution
          "leak:bin/rootcint\n"  // direct calls into non-sanitised libstdc++
          "leak:bin/rootcling\n" // direct calls into non-sanitised libstdc++
          "leak:bin/bash\n";     // When python imports ROOT
}

}
