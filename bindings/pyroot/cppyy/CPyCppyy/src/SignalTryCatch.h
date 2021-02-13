// Partial reproduction of ROOT's TException.h

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef CPYCPPYY_SIGNALTRYCATCH_H
#define CPYCPPYY_SIGNALTRYCATCH_H

#include <setjmp.h>
#include "CPyCppyy/CommonDefs.h"

#ifndef _WIN32
#define NEED_SIGJMP 1
#endif

struct ExceptionContext_t {
#ifdef NEED_SIGJMP
    sigjmp_buf fBuf;
#else
    jmp_buf fBuf;
#endif
};

#ifdef NEED_SIGJMP
# define SETJMP(buf) sigsetjmp(buf,1)
#else
# define SETJMP(buf) setjmp(buf)
#endif

#define RETRY \
    { \
        static ExceptionContext_t R__curr, *R__old = gException; \
        int R__code; \
        gException = &R__curr; \
        R__code = SETJMP(gException->fBuf); if (R__code) { }; {

#define TRY \
    { \
        static ExceptionContext_t R__curr, *R__old = gException; \
        int R__code; \
        gException = &R__curr; \
        if ((R__code = SETJMP(gException->fBuf)) == 0) {

#define CATCH(n) \
            gException = R__old; \
        } else { \
            int n = R__code; \
            gException = R__old;

#define ENDTRY \
        } \
        gException = R__old; \
    }

CPYCPPYY_IMPORT ExceptionContext_t *gException;

#endif
