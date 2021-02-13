#ifndef CPYCPPYY_PYRESULT_H
#define CPYCPPYY_PYRESULT_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyResult                                                                //
//                                                                          //
// Morphing return type from evaluating python expressions.                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// Python
struct _object;
typedef _object PyObject;

// Bindings
#include "CPyCppyy/CommonDefs.h"


namespace CPyCppyy {

class CPYCPPYY_CLASS_EXTERN PyResult {
public:
    PyResult();
    PyResult(PyObject* pyobject);
    PyResult(const PyResult&);
    PyResult& operator=(const PyResult&);
    virtual ~PyResult();

// conversions to standard types, may fail if unconvertible
    operator char*() const;
    operator const char*() const;
    operator char() const;

    operator long() const;
    operator int() const { return (int)operator long(); }
    operator short() const { return (short)operator long(); }

    operator unsigned long() const;
    operator unsigned int() const {
        return (unsigned int)operator unsigned long();
    }
    operator unsigned short() const {
        return (unsigned short)operator unsigned long();
    }

    operator double() const;
    operator float() const { return (float)operator double(); }

// used for both general object type and PyObject conversions
    operator void*() const;

    template<class T>
    operator T*() const { return (T*)(void*)*this; }

// used strictly for PyObject conversions
    operator PyObject*() const;

private:
    PyObject* fPyObject;            //! actual python object
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_PYRESULT_H
