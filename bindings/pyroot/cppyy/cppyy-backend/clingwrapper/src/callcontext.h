#ifndef CPYCPPYY_CALLCONTEXT_H
#define CPYCPPYY_CALLCONTEXT_H

//Bindings
#include "cpp_cppyy.h"

//ROOT
#include "Rtypes.h"

namespace CPyCppyy {

// general place holder for function parameters
    struct Parameter {
        union Value {
            bool                 fBool;
            int8_t               fInt8;
            uint8_t              fUInt8;
            short                fShort;
            unsigned short       fUShort;
            int                  fInt;
            unsigned int         fUInt;
            long                 fLong;
            intptr_t             fIntPtr;
            unsigned long        fULong;
            long long            fLLong;
            unsigned long long   fULLong;
            int64_t              fInt64;
            uint64_t             fUInt64;
            float                fFloat;
            double               fDouble;
            long double          fLDouble;
            void*                fVoidp;
        } fValue;
        void* fRef;
        char  fTypeCode;
    };

} // namespace CPyCppyy

#endif // !CPYCPPYY_CALLCONTEXT_H
