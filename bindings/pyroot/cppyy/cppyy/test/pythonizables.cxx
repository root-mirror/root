#include "pythonizables.h"

#include "Python.h"


//===========================================================================
pyzables::NakedBuffers::NakedBuffers(int size, double valx, double valy) : m_size(size) {
    m_Xbuf = new double[size];
    m_Ybuf = new double[size];

    for (int i=0; i<size; ++i) {
        m_Xbuf[i] = valx*i;
        m_Ybuf[i] = valy*i;
    }
}

pyzables::NakedBuffers::~NakedBuffers() {
    delete[] m_Xbuf;
    delete[] m_Ybuf;
}

int pyzables::NakedBuffers::GetN() { return m_size; }

double* pyzables::NakedBuffers::GetX() { return m_Xbuf; }
double* pyzables::NakedBuffers::GetY() { return m_Ybuf; }


//===========================================================================
pyzables::MyBase::~MyBase() {}
pyzables::MyDerived::~MyDerived() {}

pyzables::MyBase* pyzables::GimeDerived() {
   return new MyDerived();
}


//===========================================================================
int pyzables::Countable::sInstances = 0;
pyzables::SharedCountable_t pyzables::mine =
    pyzables::SharedCountable_t(new pyzables::Countable);

void pyzables::renew_mine() { mine = std::shared_ptr<Countable>(new Countable); }

pyzables::SharedCountable_t pyzables::gime_mine() { return mine; }
pyzables::SharedCountable_t* pyzables::gime_mine_ptr() { return &mine; }
pyzables::SharedCountable_t& pyzables::gime_mine_ref() { return mine; }

unsigned int pyzables::pass_mine_sp(std::shared_ptr<Countable> ptr) { return ptr->m_check; }
unsigned int pyzables::pass_mine_sp_ref(std::shared_ptr<Countable>& ptr) { return ptr->m_check; }
unsigned int pyzables::pass_mine_sp_ptr(std::shared_ptr<Countable>* ptr) { return (*ptr)->m_check; }

unsigned int pyzables::pass_mine_rp(Countable c) { return c.m_check; }
unsigned int pyzables::pass_mine_rp_ref(const Countable& c) { return c.m_check; }
unsigned int pyzables::pass_mine_rp_ptr(const Countable* c) { return c->m_check; }

pyzables::Countable* pyzables::gime_naked_countable() { return new Countable{}; }


//===========================================================================
pyzables::WithCallback1::WithCallback1(int i) : m_int(i) {}

int pyzables::WithCallback1::get_int() { return m_int; }
void pyzables::WithCallback1::set_int(int i) { m_int = i; }

static inline void replace_method_name(PyObject* klass, const char* n1, const char* n2) {
    PyObject* meth = PyObject_GetAttrString(klass, n1);
    PyObject_SetAttrString(klass, n2, meth);
    Py_DECREF(meth);
    PyObject_DelAttrString(klass, n1);
}

void pyzables::WithCallback1::WithCallback1::__cppyy_explicit_pythonize__(PyObject* klass, const std::string& name) {
// change methods to camel case
    replace_method_name(klass, "get_int", "GetInt");
    replace_method_name(klass, "set_int", "SetInt");

// store the provided class name
    klass_name = name;
}

std::string pyzables::WithCallback1::klass_name{"not set"};

pyzables::WithCallback2::WithCallback2(int i) : m_int(i) {}

int pyzables::WithCallback2::get_int() { return m_int; }
void pyzables::WithCallback2::set_int(int i) { m_int = i; }

void pyzables::WithCallback2::WithCallback2::__cppyy_pythonize__(PyObject* klass, const std::string& name) {
// change methods to camel case
    replace_method_name(klass, "get_int", "GetInt");
    replace_method_name(klass, "set_int", "SetInt");

// store the provided class name
    klass_name = name;
}

std::string pyzables::WithCallback2::klass_name{"not set"};

int pyzables::WithCallback3::get_int() { return 2*m_int; }
void pyzables::WithCallback3::set_int(int i) { m_int = 2*i; }
