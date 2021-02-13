# - Define GNU standard installation directories
# Provides install directory variables as defined for GNU software:
#  http://www.gnu.org/prep/standards/html_node/Directory-Variables.html
# Inclusion of this module defines the following variables:
#  CMAKE_INSTALL_<dir>      - destination for files of a given type
#  CMAKE_INSTALL_FULL_<dir> - corresponding absolute path
# where <dir> is one of:
#  BINDIR           - user executables (bin)
#  LIBDIR           - object code libraries (lib or lib64 or lib/<multiarch-tuple> on Debian)
#  INCLUDEDIR       - C/C++ header files (include)
#  SYSCONFDIR       - read-only single-machine data (etc)
#  PYTHONDIR        - python libraries and modules (same as LIBDIR)
#  DATAROOTDIR      - read-only architecture-independent data (share)
#  DATADIR          - read-only architecture-independent data (DATAROOTDIR/root)
#  MANDIR           - man documentation (DATAROOTDIR/man)
#  MACRODIR         - ROOT macros (DATAROOTDIR/macros)
#  CINTINCDIR       - CINT include files (LIBDIR/cint)
#  ICONDIR          - icons (DATAROOTDIR/icons)
#  SRCDIR           - sources (DATAROOTDIR/src)
#  FONTDIR          - fonts (DATAROOTDIR/fonts)
#  DOCDIR           - documentation root (DATAROOTDIR/doc/PROJECT_NAME)
#  TUTDIR           - tutorials (DOCDIR/tutorials)
#  CMAKEDIR         - cmake modules (DATAROOTDIR/cmake)
#  ELISPDIR         - lisp files (DATAROOTDIR/emacs/site-lisp)
#
# Each CMAKE_INSTALL_<dir> value may be passed to the DESTINATION options of
# install() commands for the corresponding file type.  If the includer does
# not define a value the above-shown default will be used and the value will
# appear in the cache for editing by the user.
# Each CMAKE_INSTALL_FULL_<dir> value contains an absolute path constructed
# from the corresponding destination by prepending (if necessary) the value
# of CMAKE_INSTALL_PREFIX.

#=============================================================================
# Copyright 2011 Nikita Krupen'ko <krnekit@gmail.com>
# Copyright 2011 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

if(NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "user executables (bin)")
endif()

if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  if(gnuinstall)
    set(CMAKE_INSTALL_LIBDIR "lib/root" CACHE PATH "object code libraries (lib/root)")
  else()
    set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "object code libraries (lib)")
  endif()
endif()

if(NOT DEFINED CMAKE_INSTALL_PYTHONDIR)
  if(MSVC)
    set(CMAKE_INSTALL_PYTHONDIR "${CMAKE_INSTALL_BINDIR}" CACHE PATH "python libraries and modules (same as BINDIR)")
  else()
    set(CMAKE_INSTALL_PYTHONDIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "python libraries and modules (same as LIBDIR)")
  endif()
endif()

if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
  if(gnuinstall)
    set(CMAKE_INSTALL_INCLUDEDIR "include/root" CACHE PATH "C header files (include/root)")
  else()
    set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "C header files (include)")
  endif()
endif()

if(NOT DEFINED CMAKE_INSTALL_SYSCONFDIR)
  if(gnuinstall)
    set(CMAKE_INSTALL_SYSCONFDIR "etc/root" CACHE PATH "read-only single-machine data (etc/root)")
  else()
    set(CMAKE_INSTALL_SYSCONFDIR "etc" CACHE PATH "read-only single-machine data (etc)")
  endif()
endif()

if(NOT DEFINED CMAKE_INSTALL_DATAROOTDIR)
  if(gnuinstall)
    set(CMAKE_INSTALL_DATAROOTDIR "share" CACHE PATH "root for the data (share)")
  else()
    set(CMAKE_INSTALL_DATAROOTDIR "." CACHE PATH "root for the data ()")
  endif()
endif()

#-----------------------------------------------------------------------------
# Values whose defaults are relative to DATAROOTDIR.  Store empty values in
# the cache and store the defaults in local variables if the cache values are
# not set explicitly.  This auto-updates the defaults as DATAROOTDIR changes.

if(NOT CMAKE_INSTALL_CINTINCDIR)
  if(gnuinstall)
    set(CMAKE_INSTALL_CINTINCDIR "" CACHE PATH "cint includes and libraries libraries  (LIBDIR/cint)")
    set(CMAKE_INSTALL_CINTINCDIR "${CMAKE_INSTALL_LIBDIR}/cint")
  else()
    set(CMAKE_INSTALL_CINTINCDIR "cint" CACHE PATH "cint includes and libraries libraries (cint)")
  endif()
endif()

if(NOT CMAKE_INSTALL_DATADIR)
  set(CMAKE_INSTALL_DATADIR "" CACHE PATH "read-only architecture-independent data (DATAROOTDIR)/root")
  if(gnuinstall)
    set(CMAKE_INSTALL_DATADIR "${CMAKE_INSTALL_DATAROOTDIR}/root")
  else()
    set(CMAKE_INSTALL_DATADIR ".")
  endif()
endif()

if(NOT CMAKE_INSTALL_MANDIR)
  set(CMAKE_INSTALL_MANDIR "" CACHE PATH "man documentation (DATAROOTDIR/man)")
  if(gnuinstall)
    set(CMAKE_INSTALL_MANDIR "${CMAKE_INSTALL_DATAROOTDIR}/man")
  else()
    set(CMAKE_INSTALL_MANDIR "man")
  endif()
endif()

if(NOT CMAKE_INSTALL_MACRODIR)
  set(CMAKE_INSTALL_MACRODIR "" CACHE PATH "macros documentation (DATADIR/macros)")
  if(gnuinstall)
    set(CMAKE_INSTALL_MACRODIR "${CMAKE_INSTALL_DATADIR}/macros")
  else()
    set(CMAKE_INSTALL_MACRODIR "macros")
  endif()
endif()

if(NOT CMAKE_INSTALL_ICONDIR)
  set(CMAKE_INSTALL_ICONDIR "" CACHE PATH "icons (DATADIR/icons)")
  if(gnuinstall)
    set(CMAKE_INSTALL_ICONDIR "${CMAKE_INSTALL_DATADIR}/icons")
  else()
    set(CMAKE_INSTALL_ICONDIR "icons")
  endif()
endif()

if(NOT CMAKE_INSTALL_FONTDIR)
  set(CMAKE_INSTALL_FONTDIR "" CACHE PATH "fonts (DATADIR/fonts)")
  if(gnuinstall)
    set(CMAKE_INSTALL_FONTDIR "${CMAKE_INSTALL_DATADIR}/fonts")
  else()
    set(CMAKE_INSTALL_FONTDIR "fonts")
  endif()
endif()

if(NOT CMAKE_INSTALL_JSROOTDIR)
  set(CMAKE_INSTALL_JSROOTDIR "" CACHE PATH "js (DATADIR/js)")
  if(gnuinstall)
    set(CMAKE_INSTALL_JSROOTDIR "${CMAKE_INSTALL_DATADIR}/js")
  else()
    set(CMAKE_INSTALL_JSROOTDIR "js")
  endif()
endif()

if(NOT CMAKE_INSTALL_OPENUI5DIR)
  set(CMAKE_INSTALL_OPENUI5DIR "" CACHE PATH "ui5 (DATADIR/ui5)")
  if(gnuinstall)
    set(CMAKE_INSTALL_OPENUI5DIR "${CMAKE_INSTALL_DATADIR}/ui5")
  else()
    set(CMAKE_INSTALL_OPENUI5DIR "ui5")
  endif()
endif()

if(NOT CMAKE_INSTALL_SRCDIR)
  set(CMAKE_INSTALL_SRCDIR "" CACHE PATH "sources (DATADIR/src)")
  if(gnuinstall)
    set(CMAKE_INSTALL_SRCDIR "${CMAKE_INSTALL_DATADIR}/src")
  else()
    set(CMAKE_INSTALL_SRCDIR "src")
  endif()
endif()

if(NOT CMAKE_INSTALL_CMAKEDIR)
  set(CMAKE_INSTALL_CMAKEDIR "" CACHE PATH "CMake modules (DATAROOTDIR/cmake)")
  if(gnuinstall)
    set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATADIR}/cmake")
  else()
    set(CMAKE_INSTALL_CMAKEDIR "cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_ELISPDIR)
  set(CMAKE_INSTALL_ELISPDIR "" CACHE PATH "Lisp files (DATAROOTDIR/emacs/site-lisp)")
  if(gnuinstall)
    set(CMAKE_INSTALL_ELISPDIR "${CMAKE_INSTALL_DATAROOTDIR}/emacs/site-lisp")
  else()
    set(CMAKE_INSTALL_ELISPDIR "emacs/site-lisp")
  endif()
endif()

if(NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (DATAROOTDIR/doc/root)")
  if(gnuinstall)
    set(CMAKE_INSTALL_DOCDIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/root")
  else()
    set(CMAKE_INSTALL_DOCDIR ".")
  endif()
endif()

if(NOT CMAKE_INSTALL_TUTDIR)
  set(CMAKE_INSTALL_TUTDIR "" CACHE PATH "root tutorials (DOCDIR/tutorials)")
  if(gnuinstall)
    set(CMAKE_INSTALL_TUTDIR "${CMAKE_INSTALL_DOCDIR}/tutorials")
  else()
    set(CMAKE_INSTALL_TUTDIR "tutorials")
  endif()
endif()


#-----------------------------------------------------------------------------

mark_as_advanced(
  CMAKE_INSTALL_BINDIR
  CMAKE_INSTALL_LIBDIR
  CMAKE_INSTALL_INCLUDEDIR
  CMAKE_INSTALL_SYSCONFDIR
  CMAKE_INSTALL_MANDIR
  CMAKE_INSTALL_DATAROOTDIR
  CMAKE_INSTALL_DATADIR
  CMAKE_INSTALL_MACRODIR
  CMAKE_INSTALL_CINTINCDIR
  CMAKE_INSTALL_ICONDIR
  CMAKE_INSTALL_FONTDIR
  CMAKE_INSTALL_SRCDIR
  CMAKE_INSTALL_DOCDIR
  CMAKE_INSTALL_TUTDIR
  CMAKE_INSTALL_ELISPDIR
  CMAKE_INSTALL_CMAKEDIR
  )

# Result directories
#
foreach(dir BINDIR
            LIBDIR
            PYTHONDIR
            INCLUDEDIR
            SYSCONFDIR
            MANDIR
            DATAROOTDIR
            DATADIR
            MACRODIR
            CINTINCDIR
            ICONDIR
            FONTDIR
            SRCDIR
            DOCDIR
            TUTDIR
            ELISPDIR
            CMAKEDIR )
  if(NOT IS_ABSOLUTE ${CMAKE_INSTALL_${dir}})
    set(CMAKE_INSTALL_FULL_${dir} "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_${dir}}")
  else()
    set(CMAKE_INSTALL_FULL_${dir} "${CMAKE_INSTALL_${dir}}")
  endif()
endforeach()
