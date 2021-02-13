// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jul 2004

#ifndef PYROOT_PYTHONIZE_H
#define PYROOT_PYTHONIZE_H

// Standard
#include <string>


namespace PyROOT {

// make the named ROOT class more python-like
   Bool_t Pythonize( PyObject* pyclass, const std::string& name );

} // namespace PyROOT

#endif // !PYROOT_PYTHONIZE_H
