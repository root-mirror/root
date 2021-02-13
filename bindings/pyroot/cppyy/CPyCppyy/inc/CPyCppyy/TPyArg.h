#ifndef CPYCPPYY_TPYARG_H
#define CPYCPPYY_TPYARG_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyArg                                                                   //
//                                                                          //
// Morphing argument type from evaluating python expressions.               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// Python
struct _object;
typedef _object PyObject;

// Standard
#include <vector>

// Bindings
#include "CPyCppyy/CommonDefs.h"


class CPYCPPYY_CLASS_EXTERN TPyArg {
public:
// converting constructors
    TPyArg(PyObject*);
    TPyArg(int);
    TPyArg(long);
    TPyArg(double);
    TPyArg(const char*);

    TPyArg(const TPyArg&);
    TPyArg& operator=(const TPyArg&);
    virtual ~TPyArg();

// "extractor"
    operator PyObject*() const;

// constructor and generic dispatch
    static void CallConstructor(
        PyObject*& pyself, PyObject* pyclass, const std::vector<TPyArg>& args);
    static void CallConstructor(PyObject*& pyself, PyObject* pyclass);   // default ctor
    static PyObject* CallMethod(PyObject* pymeth, const std::vector<TPyArg>& args);
    static void CallDestructor(
        PyObject*& pyself, PyObject* pymeth, const std::vector<TPyArg>& args);
    static void CallDestructor(PyObject*& pyself);

private:
    mutable PyObject* fPyObject;        //! converted C++ value as python object
};

#endif // !CPYCPPYY_TPYARG_H
