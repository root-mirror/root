// Bindings
#include "CPyCppyy.h"
#include "CPPClassMethod.h"


//- public members --------------------------------------------------------------
PyObject* CPyCppyy::CPPClassMethod::Call(
    CPPInstance*&, PyObject* args, PyObject* kwds, CallContext* ctxt)
{
// preliminary check in case keywords are accidently used (they are ignored otherwise)
    if (kwds && PyDict_Size(kwds)) {
        PyErr_SetString(PyExc_TypeError, "keyword arguments are not yet supported");
        return nullptr;
    }

// setup as necessary
    if (!this->Initialize(ctxt))
        return nullptr;

// translate the arguments
    if (!this->ConvertAndSetArgs(args, ctxt))
        return nullptr;

// execute function
    return this->Execute(nullptr, 0, ctxt);
}
