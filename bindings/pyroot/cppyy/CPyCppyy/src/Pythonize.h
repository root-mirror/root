#ifndef CPYCPPYY_PYTHONIZE_H
#define CPYCPPYY_PYTHONIZE_H

// Standard
#include <string>


namespace CPyCppyy {

// make the named C++ class more python-like
bool Pythonize(PyObject* pyclass, const std::string& name);

} // namespace CPyCppyy

#endif // !CPYCPPYY_PYTHONIZE_H
