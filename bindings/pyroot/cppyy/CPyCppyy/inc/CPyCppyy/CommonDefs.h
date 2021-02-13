#ifndef CPYCPPYY_COMMONDEFS_H
#define CPYCPPYY_COMMONDEFS_H

// export macros for our own API
// import/export (after precommondefs.h from PyPy)
#ifdef _MSC_VER
// Windows requires symbols to be explicitly exported
#define CPYCPPYY_EXPORT extern __declspec(dllexport)
#define CPYCPPYY_IMPORT extern __declspec(dllimport)
#define CPYCPPYY_CLASS_EXPORT __declspec(dllexport)

// CPYCPPYY_EXTERN is dual use in the public API
#ifndef CPYCPPYY_INTERNAL
#define CPYCPPYY_EXTERN extern __declspec(dllexport)
#define CPYCPPYY_CLASS_EXTERN __declspec(dllexport)
#else
#define CPYCPPYY_EXTERN extern __declspec(dllimport)
#define CPYCPPYY_CLASS_EXTERN __declspec(dllimport)
#endif

#define CPYCPPYY_STATIC

#else
// Linux, Mac, etc.
#define CPYCPPYY_EXPORT extern
#define CPYCPPYY_IMPORT extern
#define CPYCPPYY_CLASS_EXPORT
#define CPYCPPYY_EXTERN extern
#define CPYCPPYY_CLASS_EXTERN
#define CPYCPPYY_STATIC static

#endif

#endif // !CPYCPPYY_COMMONDEFS_H
