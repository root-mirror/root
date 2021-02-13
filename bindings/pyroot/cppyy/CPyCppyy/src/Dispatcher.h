#ifndef CPYCPPYY_DISPATCHER_H
#define CPYCPPYY_DISPATCHER_H

// Standard
#include <sstream>

namespace CPyCppyy {

class CPPScope;

// helper that inserts dispatchers for virtual methods
bool InsertDispatcher(CPPScope* klass, PyObject* bases, PyObject* dct, std::ostringstream& err);

} // namespace CPyCppyy

#endif // !CPYCPPYY_DISPATCHER_H
