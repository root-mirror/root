// Bindings
#include "CPyCppyy.h"
#include "CPPInstance.h"
#include "CPPScope.h"
#include "CPPOverload.h"
#include "MemoryRegulator.h"
#include "ProxyWrappers.h"
#include "PyStrings.h"
#include "TypeManip.h"
#include "Utility.h"

#include "CPyCppyy/DispatchPtr.h"

// Standard
#include <algorithm>
#include <sstream>


//______________________________________________________________________________
//                          Python-side proxy objects
//                          =========================
//
// C++ objects are represented in Python by CPPInstances, which encapsulate
// them using either a pointer (normal), pointer-to-pointer (kIsReference set),
// or as an owned value (kIsValue set). Objects held as reference are never
// owned, otherwise the object is owned if kIsOwner is set.
//
// In addition to encapsulation, CPPInstance offers rudimentary comparison
// operators (based on pointer value and class comparisons); stubs (with lazy
// lookups) for numeric operators; and a representation that prints the C++
// pointer values, rather than the PyObject* ones as is the default.
//
// Smart pointers have the underlying type as the Python type, but store the
// pointer to the smart pointer. They carry a pointer to the Python-sode smart
// class for dereferencing to get to the actual instance pointer.


//- private helpers ----------------------------------------------------------
namespace {

// Several specific use cases require extra data in a CPPInstance, but can not
// be a new type. E.g. cross-inheritance derived types are by definition added
// a posterio, and caching of datamembers is up to the datamember, not the
// instance type. To not have normal use of CPPInstance take extra memory, this
// extended data can slot in place of fObject for those use cases.

struct ExtendedData {
    ExtendedData() : fObject(nullptr), fSmartClass(nullptr), fTypeSize(0), fLastState(nullptr), fDispatchPtr(nullptr) {}
    ~ExtendedData() {
        for (auto& pc : fDatamemberCache)
            Py_XDECREF(pc.second);
        fDatamemberCache.clear();
    }

// the original object reference it replaces (Note: has to be first data member, see usage
// in GetObjectRaw(), e.g. for ptr-ptr passing)
    void* fObject;

// for smart pointer types
    CPyCppyy::CPPSmartClass* fSmartClass;
    size_t         fTypeSize;
    void*          fLastState;

// for caching expensive-to-create data member representations
    CPyCppyy::CI_DatamemberCache_t fDatamemberCache;

// for back-referencing from Python-derived instances
    CPyCppyy::DispatchPtr* fDispatchPtr;
};

} // unnamed namespace

#define EXT_OBJECT(pyobj)  ((ExtendedData*)((pyobj)->fObject))->fObject
#define SMART_CLS(pyobj)   ((ExtendedData*)((pyobj)->fObject))->fSmartClass
#define SMART_TYPE(pyobj)  SMART_CLS(pyobj)->fCppType
#define DISPATCHPTR(pyobj) ((ExtendedData*)((pyobj)->fObject))->fDispatchPtr
#define DATA_CACHE(pyobj)  ((ExtendedData*)((pyobj)->fObject))->fDatamemberCache

inline void CPyCppyy::CPPInstance::CreateExtension() {
    if (fFlags & kIsExtended)
        return;
    void* obj = fObject;
    fObject = (void*)new ExtendedData{};
    EXT_OBJECT(this) = obj;
    fFlags |= kIsExtended;
}

void* CPyCppyy::CPPInstance::GetExtendedObject()
{
    if (IsSmart()) {
    // We get the raw pointer from the smart pointer each time, in case it has
    // changed or has been freed.
        return Cppyy::CallR(SMART_CLS(this)->fDereferencer, EXT_OBJECT(this), 0, nullptr);
    }
    return EXT_OBJECT(this);
}


//- public methods -----------------------------------------------------------
CPyCppyy::CPPInstance* CPyCppyy::CPPInstance::Copy(void* cppinst)
{
// create a fresh instance; args and kwds are not used by op_new (see below)
    PyObject* self = (PyObject*)this;
    PyTypeObject* pytype = Py_TYPE(self);
    PyObject* newinst = pytype->tp_new(pytype, nullptr, nullptr);

// set the C++ instance as given
    ((CPPInstance*)newinst)->fObject = cppinst;

// look for user-provided __cpp_copy__ (not reusing __copy__ b/c of differences
// in semantics: need to pass in the new instance) ...
    PyObject* cpy = PyObject_GetAttrString(self, (char*)"__cpp_copy__");
    if (cpy && PyCallable_Check(cpy)) {
        PyObject* args = PyTuple_New(1);
        Py_INCREF(newinst);
        PyTuple_SET_ITEM(args, 0, newinst);
        PyObject* res = PyObject_CallObject(cpy, args);
        Py_DECREF(args);
        Py_DECREF(cpy);
        if (res) {
            Py_DECREF(res);
            return (CPPInstance*)newinst;
        }

    // error already set, but need to return nullptr
        Py_DECREF(newinst);
        return nullptr;
    } else if (cpy)
        Py_DECREF(cpy);
    else
        PyErr_Clear();

// ... otherwise, shallow copy any Python-side dictionary items
    PyObject* selfdct = PyObject_GetAttr(self, PyStrings::gDict);
    PyObject* newdct  = PyObject_GetAttr(newinst, PyStrings::gDict);
    bool bMergeOk = PyDict_Merge(newdct, selfdct, 1) == 0;
    Py_DECREF(newdct);
    Py_DECREF(selfdct);

    if (!bMergeOk) {
    // presume error already set
        Py_DECREF(newinst);
        return nullptr;
    }

    MemoryRegulator::RegisterPyObject((CPPInstance*)newinst, cppinst);
    return (CPPInstance*)newinst;
}


//----------------------------------------------------------------------------
void CPyCppyy::CPPInstance::PythonOwns()
{
    fFlags |= kIsOwner;
    if ((fFlags & kIsExtended) && DISPATCHPTR(this))
        DISPATCHPTR(this)->PythonOwns();
}

//----------------------------------------------------------------------------
void CPyCppyy::CPPInstance::CppOwns()
{
    fFlags &= ~kIsOwner;
    if ((fFlags & kIsExtended) && DISPATCHPTR(this))
        DISPATCHPTR(this)->CppOwns();
}

//----------------------------------------------------------------------------
void CPyCppyy::CPPInstance::SetSmart(PyObject* smart_type)
{
    CreateExtension();
    Py_INCREF(smart_type);
    SMART_CLS(this) = (CPPSmartClass*)smart_type;
    fFlags |= kIsSmartPtr;
}

//----------------------------------------------------------------------------
Cppyy::TCppType_t CPyCppyy::CPPInstance::GetSmartIsA() const
{
    if (!IsSmart()) return (Cppyy::TCppType_t)0;
    return SMART_TYPE(this);
}

//----------------------------------------------------------------------------
CPyCppyy::CI_DatamemberCache_t& CPyCppyy::CPPInstance::GetDatamemberCache()
{
// Return the cache for expensive data objects (and make extended as necessary)
    CreateExtension();
    return DATA_CACHE(this);
}

//----------------------------------------------------------------------------
void CPyCppyy::CPPInstance::SetDispatchPtr(void* ptr)
{
// Set up the dispatch pointer for memory management
    CreateExtension();
    DISPATCHPTR(this) = (DispatchPtr*)ptr;
}


//----------------------------------------------------------------------------
void CPyCppyy::op_dealloc_nofree(CPPInstance* pyobj) {
// Destroy the held C++ object, if owned; does not deallocate the proxy.

    Cppyy::TCppType_t klass = pyobj->ObjectIsA(false /* check_smart */);
    void*& cppobj = pyobj->GetObjectRaw();

    if (pyobj->fFlags & CPPInstance::kIsRegulated)
        MemoryRegulator::UnregisterPyObject(pyobj, (PyObject*)Py_TYPE((PyObject*)pyobj));

    if (pyobj->fFlags & CPPInstance::kIsOwner) {
       if (pyobj->fFlags & CPPInstance::kIsValue) {
           Cppyy::CallDestructor(klass, cppobj);
           Cppyy::Deallocate(klass, cppobj);
       } else {
           if (cppobj) Cppyy::Destruct(klass, cppobj);
       }
    }
    cppobj = nullptr;

    if (pyobj->IsExtended()) delete (ExtendedData*)pyobj->fObject;
    pyobj->fFlags = CPPInstance::kNoWrapConv;
}


namespace CPyCppyy {

//= CPyCppyy object proxy null-ness checking =================================
static int op_nonzero(CPPInstance* self)
{
// Null of the proxy is determined by null-ness of the held C++ object.
    if (!self->GetObject())
        return 0;

// If the object is valid, then the normal Python behavior is to allow __len__
// to determine truth. However, that function is defined in typeobject.c and only
// installed if tp_as_number exists w/o the nb_nonzero/nb_bool slot filled in, so
// it can not be called directly. Instead, since we're only ever dealing with
// CPPInstance derived objects, ignore length from sequence or mapping and call
// the __len__ method, if any, directly.

    PyObject* pylen = PyObject_CallMethodObjArgs((PyObject*)self, PyStrings::gLen, NULL);
    if (!pylen) {
        PyErr_Clear();
        return 1;       // since it's still a valid object
    }

    int result = PyObject_IsTrue(pylen);
    Py_DECREF(pylen);
    return result;
}

//= CPyCppyy object explicit destruction =====================================
static PyObject* op_destruct(CPPInstance* self)
{
// User access to force deletion of the object. Needed in case of a true
// garbage collector (like in PyPy), to allow the user control over when
// the C++ destructor is called. This method requires that the C++ object
// is owned (no-op otherwise).
    op_dealloc_nofree(self);
    Py_RETURN_NONE;
}

//= CPyCppyy object dispatch support =========================================
static PyObject* op_dispatch(PyObject* self, PyObject* args, PyObject* /* kdws */)
{
// User-side __dispatch__ method to allow selection of a specific overloaded
// method. The actual selection is in the __overload__() method of CPPOverload.
    PyObject *mname = nullptr, *sigarg = nullptr;
    if (!PyArg_ParseTuple(args, const_cast<char*>("O!O!:__dispatch__"),
            &CPyCppyy_PyText_Type, &mname, &CPyCppyy_PyText_Type, &sigarg))
        return nullptr;

// get the named overload
    PyObject* pymeth = PyObject_GetAttr(self, mname);
    if (!pymeth)
        return nullptr;

// get the '__overload__' method to allow overload selection
    PyObject* pydisp = PyObject_GetAttrString(pymeth, const_cast<char*>("__overload__"));
    if (!pydisp) {
        Py_DECREF(pymeth);
        return nullptr;
    }

// finally, call dispatch to get the specific overload
    PyObject* oload = PyObject_CallFunctionObjArgs(pydisp, sigarg, nullptr);
    Py_DECREF(pydisp);
    Py_DECREF(pymeth);
    return oload;
}

//= CPyCppyy smart pointer support ===========================================
static PyObject* op_get_smartptr(CPPInstance* self)
{
    if (!self->IsSmart()) {
    // TODO: more likely should raise
        Py_RETURN_NONE;
    }

    return CPyCppyy::BindCppObjectNoCast(self->GetSmartObject(), SMART_TYPE(self), CPPInstance::kNoWrapConv);
}


//----------------------------------------------------------------------------
static PyMethodDef op_methods[] = {
    {(char*)"__destruct__", (PyCFunction)op_destruct, METH_NOARGS, nullptr},
    {(char*)"__dispatch__", (PyCFunction)op_dispatch, METH_VARARGS,
      (char*)"dispatch to selected overload"},
    {(char*)"__smartptr__", (PyCFunction)op_get_smartptr, METH_NOARGS,
      (char*)"get associated smart pointer, if any"},
    {(char*)nullptr, nullptr, 0, nullptr}
};


//= CPyCppyy object proxy construction/destruction ===========================
static CPPInstance* op_new(PyTypeObject* subtype, PyObject*, PyObject*)
{
// Create a new object proxy (holder only).
    CPPInstance* pyobj = (CPPInstance*)subtype->tp_alloc(subtype, 0);
    pyobj->fObject = nullptr;
    pyobj->fFlags = CPPInstance::kNoWrapConv;

    return pyobj;
}

//----------------------------------------------------------------------------
static void op_dealloc(CPPInstance* pyobj)
{
// Remove (Python-side) memory held by the object proxy.
    PyObject_GC_UnTrack((PyObject*)pyobj);
    op_dealloc_nofree(pyobj);
    PyObject_GC_Del((PyObject*)pyobj);
}

//----------------------------------------------------------------------------
static int op_clear(CPPInstance* pyobj)
{
// Garbage collector clear of held python member objects; this is a good time
// to safely remove this object from the memory regulator.
    if (pyobj->fFlags & CPPInstance::kIsRegulated)
        MemoryRegulator::UnregisterPyObject(pyobj, (PyObject*)Py_TYPE((PyObject*)pyobj));;

    return 0;
}

//----------------------------------------------------------------------------
static inline PyObject* eqneq_binop(CPPClass* klass, PyObject* self, PyObject* obj, int op)
{
    using namespace Utility;

    if (!klass->fOperators)
        klass->fOperators = new PyOperators{};

    bool flipit = false;
    PyObject* binop = op == Py_EQ ? klass->fOperators->fEq : klass->fOperators->fNe;
    if (!binop) {
        const char* cppop = op == Py_EQ ? "==" : "!=";
        PyCallable* pyfunc = FindBinaryOperator(self, obj, cppop);
        if (pyfunc) binop = (PyObject*)CPPOverload_New(cppop, pyfunc);
        else {
            Py_INCREF(Py_None);
            binop = Py_None;
        }
    // sets the operator to Py_None if not found, indicating that search was done
        if (op == Py_EQ) klass->fOperators->fEq = binop;
        else klass->fOperators->fNe = binop;
    }

    if (binop == Py_None) {  // can try !== or !!= as alternatives
        binop = op == Py_EQ ? klass->fOperators->fNe : klass->fOperators->fEq;
        if (binop && binop != Py_None) flipit = true;
    }

    if (!binop || binop == Py_None) return nullptr;

    PyObject* args = PyTuple_New(1);
    Py_INCREF(obj);  PyTuple_SET_ITEM(args, 0, obj);
// since this overload is "ours", don't have to worry about rebinding
    ((CPPOverload*)binop)->fSelf = (CPPInstance*)self;
    PyObject* result = CPPOverload_Type.tp_call(binop, args, nullptr);
    ((CPPOverload*)binop)->fSelf = nullptr;
    Py_DECREF(args);

    if (!result) {
        PyErr_Clear();
        return nullptr;
    }

// successful result, but may need to reverse the outcome
    if (!flipit) return result;

    int istrue = PyObject_IsTrue(result);
    Py_DECREF(result);
    if (istrue) {
        Py_RETURN_FALSE;
    }
    Py_RETURN_TRUE;
}

static PyObject* op_richcompare(CPPInstance* self, PyObject* other, int op)
{
// Rich set of comparison objects; only equals and not-equals are defined.
    if (op != Py_EQ && op != Py_NE) {
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }

// special case for None to compare True to a null-pointer
    if ((PyObject*)other == Py_None && !self->fObject) {
        if (op == Py_EQ) { Py_RETURN_TRUE; }
        Py_RETURN_FALSE;
    }

// use C++-side operators if available
    PyObject* result = eqneq_binop((CPPClass*)Py_TYPE(self), (PyObject*)self, other, op);
    if (!result && CPPInstance_Check(other))
        result = eqneq_binop((CPPClass*)Py_TYPE(other), other, (PyObject*)self, op);
    if (result) return result;

// default behavior: type + held pointer value defines identity (covers if
// other is not actually an CPPInstance, as ob_type will be unequal)
    bool bIsEq = false;
    if (Py_TYPE(self) == Py_TYPE(other) && \
            self->GetObject() == ((CPPInstance*)other)->GetObject())
        bIsEq = true;

    if ((op == Py_EQ && bIsEq) || (op == Py_NE && !bIsEq)) {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

//----------------------------------------------------------------------------
static PyObject* op_repr(CPPInstance* self)
{
// Build a representation string of the object proxy that shows the address
// of the C++ object that is held, as well as its type.
    PyObject* pyclass = (PyObject*)Py_TYPE(self);
    PyObject* modname = PyObject_GetAttr(pyclass, PyStrings::gModule);

    Cppyy::TCppType_t klass = self->ObjectIsA();
    std::string clName = klass ? Cppyy::GetFinalName(klass) : "<unknown>";
    if (self->fFlags & CPPInstance::kIsReference)
        clName.append("*");

    PyObject* repr = nullptr;
    if (self->IsSmart()) {
        std::string smartPtrName = Cppyy::GetScopedFinalName(SMART_TYPE(self));
        repr = CPyCppyy_PyText_FromFormat(
            const_cast<char*>("<%s.%s object at %p held by %s at %p>"),
            CPyCppyy_PyText_AsString(modname), clName.c_str(),
            self->GetObject(), smartPtrName.c_str(), self->GetObjectRaw());
    } else {
        repr = CPyCppyy_PyText_FromFormat(const_cast<char*>("<%s.%s object at %p>"),
            CPyCppyy_PyText_AsString(modname), clName.c_str(), self->GetObject());
    }

    Py_DECREF(modname);
    return repr;
}

//----------------------------------------------------------------------------
static inline Py_hash_t CPyCppyy_PyLong_AsHash_t(PyObject* obj)
{
// Cannot use PyLong_AsSize_t here, as it cuts of at PY_SSIZE_T_MAX, which is
// only half of the max of std::size_t returned by the hash.
    if (sizeof(unsigned long) >= sizeof(size_t))
        return (Py_hash_t)PyLong_AsUnsignedLong(obj);
    return (Py_hash_t)PyLong_AsUnsignedLongLong(obj);
}

static Py_hash_t op_hash(CPPInstance* self)
{
// Try to locate an std::hash for this type and use that if it exists
    CPPClass* klass = (CPPClass*)Py_TYPE(self);
    if (klass->fOperators && klass->fOperators->fHash) {
        Py_hash_t h = 0;
        PyObject* hashval = PyObject_CallFunctionObjArgs(klass->fOperators->fHash, (PyObject*)self, nullptr);
        if (hashval) {
            h = CPyCppyy_PyLong_AsHash_t(hashval);
            Py_DECREF(hashval);
        }
        return h;
    }

    Cppyy::TCppScope_t stdhash = Cppyy::GetScope("std::hash<"+Cppyy::GetScopedFinalName(self->ObjectIsA())+">");
    if (stdhash) {
        PyObject* hashcls = CreateScopeProxy(stdhash);
        PyObject* dct = PyObject_GetAttr(hashcls, PyStrings::gDict);
        bool isValid = PyMapping_HasKeyString(dct, (char*)"__call__");
        Py_DECREF(dct);
        if (isValid) {
            PyObject* hashobj = PyObject_CallObject(hashcls, nullptr);
            if (!klass->fOperators) klass->fOperators = new Utility::PyOperators{};
            klass->fOperators->fHash = hashobj;
            Py_DECREF(hashcls);

            Py_hash_t h = 0;
            PyObject* hashval = PyObject_CallFunctionObjArgs(hashobj, (PyObject*)self, nullptr);
            if (hashval) {
                h = CPyCppyy_PyLong_AsHash_t(hashval);
                Py_DECREF(hashval);
            }
            return h;
        }
        Py_DECREF(hashcls);
    }

// if not valid, simply reset the hash function so as to not kill performance
    ((PyTypeObject*)Py_TYPE(self))->tp_hash = PyBaseObject_Type.tp_hash;
    return PyBaseObject_Type.tp_hash((PyObject*)self);
}

//----------------------------------------------------------------------------
static PyObject* op_str_internal(PyObject* pyobj, PyObject* lshift, bool isBound)
{
    static Cppyy::TCppScope_t sOStringStreamID = Cppyy::GetScope("std::ostringstream");
    std::ostringstream s;
    PyObject* pys = BindCppObjectNoCast(&s, sOStringStreamID);
    PyObject* res;
    if (isBound) res = PyObject_CallFunctionObjArgs(lshift, pys, NULL);
    else res = PyObject_CallFunctionObjArgs(lshift, pys, pyobj, NULL);
    Py_DECREF(pys);
    Py_DECREF(lshift);
    if (res) {
        Py_DECREF(res);
        return CPyCppyy_PyText_FromString(s.str().c_str());
    }
    PyErr_Clear();
    return nullptr;
}

static PyObject* op_str(CPPInstance* self)
{
#ifndef _WIN64
// Forward to C++ insertion operator if available, otherwise forward to repr.
    PyObject* result = nullptr;
    PyObject* pyobj = (PyObject*)self;
    PyObject* lshift = PyObject_GetAttr(pyobj, PyStrings::gLShift);
    if (lshift) result = op_str_internal(pyobj, lshift, true);

    if (!result) {
        PyErr_Clear();
        PyObject* pyclass = (PyObject*)Py_TYPE(pyobj);
        lshift = PyObject_GetAttr(pyclass, PyStrings::gLShiftC);
        if (!lshift) {
            PyErr_Clear();
        // attempt lazy install of global operator<<(ostream&)
            std::string rcname = Utility::ClassName(pyobj);
            Cppyy::TCppScope_t rnsID = Cppyy::GetScope(TypeManip::extract_namespace(rcname));
            PyCallable* pyfunc = Utility::FindBinaryOperator("std::ostream", rcname, "<<", rnsID);
            if (pyfunc) {
                Utility::AddToClass(pyclass, "__lshiftc__", pyfunc);
                lshift = PyObject_GetAttr(pyclass, PyStrings::gLShiftC);
            } else
                PyType_Type.tp_setattro(pyclass, PyStrings::gLShiftC, Py_None);
        } else if (lshift == Py_None) {
            Py_DECREF(lshift);
            lshift = nullptr;
        }
        if (lshift) result = op_str_internal(pyobj, lshift, false);
    }

    if (result)
        return result;
#endif  //!_WIN64

    return op_repr(self);
}

//-----------------------------------------------------------------------------
static PyObject* op_getownership(CPPInstance* pyobj, void*)
{
    return PyBool_FromLong((long)(pyobj->fFlags & CPPInstance::kIsOwner));
}

//-----------------------------------------------------------------------------
static int op_setownership(CPPInstance* pyobj, PyObject* value, void*)
{
// Set the ownership (True is python-owns) for the given object.
    long shouldown = PyLong_AsLong(value);
    if (shouldown == -1 && PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "__python_owns__ should be either True or False");
        return -1;
    }

    (bool)shouldown ? pyobj->PythonOwns() : pyobj->CppOwns();

    return 0;
}


//-----------------------------------------------------------------------------
static PyGetSetDef op_getset[] = {
    {(char*)"__python_owns__", (getter)op_getownership, (setter)op_setownership,
      (char*)"If true, python manages the life time of this object", nullptr},
 {(char*)nullptr, nullptr, nullptr, nullptr, nullptr}
};


//= CPyCppyy type number stubs to allow dynamic overrides =====================
#define CPYCPPYY_STUB_BODY(name, op)                                          \
    if (!meth) {                                                              \
        PyErr_Clear();                                                        \
        PyCallable* pyfunc = Utility::FindBinaryOperator(left, right, #op);   \
        if (pyfunc) meth = (PyObject*)CPPOverload_New(#name, pyfunc);         \
        else {                                                                \
            PyErr_SetString(PyExc_NotImplementedError, "");                   \
            return nullptr;                                                   \
        }                                                                     \
    }                                                                         \
    PyObject* res = PyObject_CallFunctionObjArgs(meth, cppobj, other, nullptr);\
    if (!res) {                                                               \
    /* try again, in case there is a better overload out there */             \
        PyErr_Clear();                                                        \
        PyCallable* pyfunc = Utility::FindBinaryOperator(left, right, #op);   \
        if (pyfunc) ((CPPOverload*&)meth)->AdoptMethod(pyfunc);               \
        else {                                                                \
            PyErr_SetString(PyExc_NotImplementedError, "");                   \
            return nullptr;                                                   \
        }                                                                     \
    /* use same overload with newly added function */                         \
        res = PyObject_CallFunctionObjArgs(meth, cppobj, other, nullptr);     \
    }                                                                         \
    return res;


#define CPYCPPYY_OPERATOR_STUB(name, op, ometh)                               \
static PyObject* op_##name##_stub(PyObject* left, PyObject* right)            \
{                                                                             \
/* placeholder to lazily install and forward to 'ometh' if available */       \
    CPPClass* klass = (CPPClass*)Py_TYPE(left);                               \
    if (!klass->fOperators) klass->fOperators = new Utility::PyOperators{};   \
    PyObject*& meth = ometh;                                                  \
    PyObject *cppobj = left, *other = right;                                  \
    CPYCPPYY_STUB_BODY(name, op)                                              \
}

#define CPYCPPYY_ASSOCIATIVE_OPERATOR_STUB(name, op, lmeth, rmeth)            \
static PyObject* op_##name##_stub(PyObject* left, PyObject* right)            \
{                                                                             \
/* placeholder to lazily install and forward do '(l/r)meth' if available  */  \
    CPPClass* klass; PyObject** pmeth;                                        \
    PyObject *cppobj, *other;                                                 \
    if (CPPInstance_Check(left)) {                                            \
        klass = (CPPClass*)Py_TYPE(left);                                     \
        if (!klass->fOperators) klass->fOperators = new Utility::PyOperators{};\
        pmeth = &lmeth; cppobj = left; other = right;                         \
    } else if (CPPInstance_Check(right)) {                                    \
        klass = (CPPClass*)Py_TYPE(right);                                    \
        if (!klass->fOperators) klass->fOperators = new Utility::PyOperators{};\
        pmeth = &rmeth; cppobj = right; other = left;                         \
    } else {                                                                  \
        PyErr_SetString(PyExc_NotImplementedError, "");                       \
        return nullptr;                                                       \
    }                                                                         \
    PyObject*& meth = *pmeth;                                                 \
    CPYCPPYY_STUB_BODY(name, op)                                              \
}

#define CPYCPPYY_UNARY_OPERATOR(name, op, label)                              \
static PyObject* op_##name##_stub(PyObject* pyobj)                            \
{                                                                             \
/* placeholder to lazily install unary operators */                           \
    PyCallable* pyfunc = Utility::FindUnaryOperator((PyObject*)Py_TYPE(pyobj), #op);\
    if (pyfunc && Utility::AddToClass((PyObject*)Py_TYPE(pyobj), #label, pyfunc))\
         return PyObject_CallMethod(pyobj, (char*)#label, nullptr);           \
    PyErr_SetString(PyExc_NotImplementedError, "");                           \
    return nullptr;                                                           \
}

CPYCPPYY_ASSOCIATIVE_OPERATOR_STUB(add, +, klass->fOperators->fLAdd, klass->fOperators->fRAdd)
CPYCPPYY_OPERATOR_STUB(            sub, -, klass->fOperators->fSub)
CPYCPPYY_ASSOCIATIVE_OPERATOR_STUB(mul, *, klass->fOperators->fLMul, klass->fOperators->fRMul)
CPYCPPYY_OPERATOR_STUB(            div, /, klass->fOperators->fDiv)
CPYCPPYY_UNARY_OPERATOR(neg,    -, __neg__)
CPYCPPYY_UNARY_OPERATOR(pos,    +, __pos__)
CPYCPPYY_UNARY_OPERATOR(invert, ~, __invert__)

//-----------------------------------------------------------------------------
static PyNumberMethods op_as_number = {
    (binaryfunc)op_add_stub,       // nb_add
    (binaryfunc)op_sub_stub,       // nb_subtract
    (binaryfunc)op_mul_stub,       // nb_multiply
#if PY_VERSION_HEX < 0x03000000
    (binaryfunc)op_div_stub,       // nb_divide
#endif
    0,                             // nb_remainder
    0,                             // nb_divmod
    0,                             // nb_power
    (unaryfunc)op_neg_stub,        // nb_negative
    (unaryfunc)op_pos_stub,        // nb_positive
    0,                             // nb_absolute
    (inquiry)op_nonzero,           // nb_bool (nb_nonzero in p2)
    (unaryfunc)op_invert_stub,     // nb_invert
    0,                             // nb_lshift
    0,                             // nb_rshift
    0,                             // nb_and
    0,                             // nb_xor
    0,                             // nb_or
#if PY_VERSION_HEX < 0x03000000
    0,                             // nb_coerce
#endif
    0,                             // nb_int
    0,                             // nb_long (nb_reserved in p3)
    0,                             // nb_float
#if PY_VERSION_HEX < 0x03000000
    0,                             // nb_oct
    0,                             // nb_hex
#endif
    0,                             // nb_inplace_add
    0,                             // nb_inplace_subtract
    0,                             // nb_inplace_multiply
#if PY_VERSION_HEX < 0x03000000
    0,                             // nb_inplace_divide
#endif
    0,                             // nb_inplace_remainder
    0,                             // nb_inplace_power
    0,                             // nb_inplace_lshift
    0,                             // nb_inplace_rshift
    0,                             // nb_inplace_and
    0,                             // nb_inplace_xor
    0                              // nb_inplace_or
#if PY_VERSION_HEX >= 0x02020000
    , 0                            // nb_floor_divide
#if PY_VERSION_HEX < 0x03000000
    , 0                            // nb_true_divide
#else
    , (binaryfunc)op_div_stub      // nb_true_divide
#endif
    , 0                            // nb_inplace_floor_divide
    , 0                            // nb_inplace_true_divide
#endif
#if PY_VERSION_HEX >= 0x02050000
    , 0                            // nb_index
#endif
#if PY_VERSION_HEX >= 0x03050000
    , 0                            // nb_matrix_multiply
    , 0                            // nb_inplace_matrix_multiply
#endif
};


//= CPyCppyy object proxy type ===============================================
PyTypeObject CPPInstance_Type = {
    PyVarObject_HEAD_INIT(&CPPScope_Type, 0)
    (char*)"cppyy.CPPInstance",    // tp_name
    sizeof(CPPInstance),           // tp_basicsize
    0,                             // tp_itemsize
    (destructor)op_dealloc,        // tp_dealloc
    0,                             // tp_print
    0,                             // tp_getattr
    0,                             // tp_setattr
    0,                             // tp_compare
    (reprfunc)op_repr,             // tp_repr
    &op_as_number,                 // tp_as_number
    0,                             // tp_as_sequence
    0,                             // tp_as_mapping
    (hashfunc)op_hash,             // tp_hash
    0,                             // tp_call
    (reprfunc)op_str,              // tp_str
    0,                             // tp_getattro
    0,                             // tp_setattro
    0,                             // tp_as_buffer
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE |
        Py_TPFLAGS_HAVE_GC |
        Py_TPFLAGS_CHECKTYPES,     // tp_flags
    (char*)"cppyy object proxy (internal)", // tp_doc
    0,                             // tp_traverse
    (inquiry)op_clear,             // tp_clear
    (richcmpfunc)op_richcompare,   // tp_richcompare
    0,                             // tp_weaklistoffset
    0,                             // tp_iter
    0,                             // tp_iternext
    op_methods,                    // tp_methods
    0,                             // tp_members
    op_getset,                     // tp_getset
    0,                             // tp_base
    0,                             // tp_dict
    0,                             // tp_descr_get
    0,                             // tp_descr_set
    0,                             // tp_dictoffset
    0,                             // tp_init
    0,                             // tp_alloc
    (newfunc)op_new,               // tp_new
    0,                             // tp_free
    0,                             // tp_is_gc
    0,                             // tp_bases
    0,                             // tp_mro
    0,                             // tp_cache
    0,                             // tp_subclasses
    0                              // tp_weaklist
#if PY_VERSION_HEX >= 0x02030000
    , 0                            // tp_del
#endif
#if PY_VERSION_HEX >= 0x02060000
    , 0                            // tp_version_tag
#endif
#if PY_VERSION_HEX >= 0x03040000
    , 0                            // tp_finalize
#endif
};

} // namespace CPyCppyy
