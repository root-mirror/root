# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# FIXME: The tests in core should require only libCore. OTOH, TQObjectTests uses the interpreter to register the class.
# This means that if we run make CoreBaseTests the executable wouldn't be runnable because it requires libCling and
# onepcm targets to be built.

if(NOT MSVC)
  set(dllib dl)
endif()

ROOT_ADD_GTEST(CoreBaseTests
  TNamedTests.cxx
  TQObjectTests.cxx
  TExceptionHandlerTests.cxx
  LIBRARIES Core Cling RIO ${dllib})

ROOT_ADD_GTEST(CoreErrorTests TErrorTests.cxx LIBRARIES Core)
