#ifndef CPPYY_CLINGCWRAPPER
#define CPPYY_CLINGCWRAPPER

#include "capi.h"

#ifdef __cplusplus
extern "C" {
#endif // ifdef __cplusplus

    /* misc helpers */
    void* cppyy_load_dictionary(const char* lib_name);

#ifdef __cplusplus
}
#endif // ifdef __cplusplus

// TODO: pick up from llvm-config --cxxflags
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

// Wrapper callback: except this to become available from Cling directly
typedef void (*CPPYY_Cling_Wrapper_t)(void*, int, void**, void*);

#endif // ifndef CPPYY_CLINGCWRAPPER
