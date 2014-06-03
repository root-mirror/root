#---Check for installed packages depending on the build options/components eamnbled -
include(ExternalProject)
include(FindPackageHandleStandardArgs)

#---Check for Cocoa/Quartz graphics backend (MacOS X only)
if(cocoa)
  if(APPLE)
    set(x11 OFF CACHE BOOL "" FORCE)
    set(builtin_freetype ON CACHE BOOL "" FORCE)
  else()
    message(STATUS "Cocoa option can only be enabled on MacOSX platform")
    set(cocoa OFF CACHE BOOL "" FORCE)
  endif()
endif()

#---Check for Zlib ------------------------------------------------------------------
if(NOT builtin_zlib)
  message(STATUS "Looking for ZLib")
  find_Package(ZLIB)
  if(NOT ZLIB_FOUND)
    message(STATUS "Zlib not found. Switching on builtin_zlib option")
    set(builtin_zlib ON CACHE BOOL "" FORCE)
   endif()
endif()
if(builtin_zlib)
  set(ZLIB_LIBRARY "")
endif()

#---Check for Freetype---------------------------------------------------------------
if(NOT builtin_freetype)
  message(STATUS "Looking for Freetype")
  find_package(Freetype)
  if(FREETYPE_FOUND)
    set(FREETYPE_INCLUDE_DIR ${FREETYPE_INCLUDE_DIR_freetype2})
  else()
    message(STATUS "FreeType not found. Switching on builtin_freetype option")
    set(builtin_freetype ON CACHE BOOL "" FORCE) 	
  endif()
endif()
if(builtin_freetype)  
  set(FREETYPE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/graf2d/freetype/freetype-2.3.12/include)
  set(FREETYPE_INCLUDE_DIRS ${FREETYPE_INCLUDE_DIR})
  if(WIN32)
    set(FREETYPE_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/freetype.lib)
  else()
    set(FREETYPE_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libfreetype.a)
  endif()
endif()

#---Check for PCRE-------------------------------------------------------------------
if(NOT builtin_pcre)
  message(STATUS "Looking for PCRE")
  find_package(PCRE)
  if(PCRE_FOUND)
  else()
    message(STATUS "PCRE not found. Switching on builtin_pcre option")
    set(builtin_pcre ON CACHE BOOL "" FORCE) 	
  endif() 
endif()
if(builtin_pcre)
  set(PCRE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/core/pcre/pcre-7.8)
  if(WIN32)
    set(PCRE_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libpcre.lib) 
  else()
    set(PCRE_LIBRARIES "-L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lpcre") 
  endif()
endif()

#---Check for LZMA-------------------------------------------------------------------
if(NOT builtin_lzma)
  message(STATUS "Looking for LZMA")
  find_package(LZMA)
  if(LZMA_FOUND)
  else()
    message(STATUS "LZMA not found. Switching on builtin_lzma option")
    set(builtin_lzma ON CACHE BOOL "" FORCE) 	
  endif() 
endif()
if(builtin_lzma)
  set(lzma_version 5.0.3)
  message(STATUS "Building LZMA version ${lzma_version} included in ROOT itself")
  if(WIN32)
    ExternalProject_Add(
	  LZMA
	  URL ${CMAKE_SOURCE_DIR}/core/lzma/src/xz-${lzma_version}-win32.tar.gz 
	  URL_MD5  65693dc257802b6778c28ed53ecca678
	  PREFIX LZMA
	  INSTALL_DIR ${CMAKE_BINARY_DIR}
      CONFIGURE_COMMAND "" BUILD_COMMAND ""
	  INSTALL_COMMAND cmake -E copy lib/liblzma.dll <INSTALL_DIR>/bin/${CMAKE_CFG_INTDIR}
	  BUILD_IN_SOURCE 1)
    install(FILES ${CMAKE_BINARY_DIR}/LZMA/src/LZMA/lib/liblzma.dll DESTINATION ${CMAKE_INSTALL_BINDIR})
    set(LZMA_LIBRARIES ${CMAKE_BINARY_DIR}/LZMA/src/LZMA/lib/liblzma.lib)
    set(LZMA_INCLUDE_DIR ${CMAKE_BINARY_DIR}/LZMA/src/LZMA/include)
  else() 
    if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
      set(LZMA_CFLAGS "-Wno-format-nonliteral")
    elseif( CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      set(LZMA_CFLAGS "-wd188 -wd181 -wd1292 -wd10006 -wd10156 -wd2259 -wd981 -wd128 -wd3179")
    endif()
    ExternalProject_Add(
      LZMA
      URL ${CMAKE_SOURCE_DIR}/core/lzma/src/xz-${lzma_version}.tar.gz 
      URL_MD5 858405e79590e9b05634c399497f4ba7
      INSTALL_DIR ${CMAKE_BINARY_DIR}
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix <INSTALL_DIR> --with-pic --disable-shared CFLAGS=${LZMA_CFLAGS}
      BUILD_IN_SOURCE 1)
    set(LZMA_LIBRARIES -L${CMAKE_BINARY_DIR}/lib -llzma)
    set(LZMA_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
  endif()
endif()


#---Check for X11 which is mandatory lib on Unix--------------------------------------
if(x11)
  message(STATUS "Looking for X11")
  if(X11_X11_INCLUDE_PATH)
    set(X11_FIND_QUIETLY 1)
  endif()
  find_package(X11 REQUIRED)
  if(X11_FOUND)
    list(REMOVE_DUPLICATES X11_INCLUDE_DIR)
    if(NOT X11_FIND_QUIETLY)
      message(STATUS "X11_INCLUDE_DIR: ${X11_INCLUDE_DIR}")
      message(STATUS "X11_LIBRARIES: ${X11_LIBRARIES}")
    endif()
  else()
    message(FATAL_ERROR "libX11 and X11 headers must be installed.")
  endif()
  if(X11_Xpm_FOUND)
    if(NOT X11_FIND_QUIETLY)
      message(STATUS "X11_Xpm_INCLUDE_PATH: ${X11_Xpm_INCLUDE_PATH}")
      message(STATUS "X11_Xpm_LIB: ${X11_Xpm_LIB}")
    endif()
  else()
    message(FATAL_ERROR "libXpm and Xpm headers must be installed.")
  endif()
  if(X11_Xft_FOUND)
    if(NOT X11_FIND_QUIETLY)
      message(STATUS "X11_Xft_INCLUDE_PATH: ${X11_Xft_INCLUDE_PATH}")
      message(STATUS "X11_Xft_LIB: ${X11_Xft_LIB}")
    endif()
    set(xft ON)
  else()
    message(FATAL_ERROR "libXft and Xft headers must be installed.")
  endif()
  if(X11_Xext_FOUND)
    if(NOT X11_FIND_QUIETLY)
      message(STATUS "X11_Xext_INCLUDE_PATH: ${X11_Xext_INCLUDE_PATH}")
      message(STATUS "X11_Xext_LIB: ${X11_Xext_LIB}")
    endif()
  else()
    message(FATAL_ERROR "libXext and Xext headers must be installed.")
  endif()
else()
  set(xft OFF)
endif()


#---Check for AfterImage---------------------------------------------------------------
if(NOT builtin_afterimage)
  message(STATUS "Looking for AfterImage")
  find_package(AfterImage)
  if(NOT AFTERIMAGE_FOUND)
    message(STATUS "AfterImage not found. Switching on builtin_afterimage option")
    set(builtin_afterimage ON CACHE BOOL "" FORCE) 	
  endif()
endif()

#---Check for all kind of graphics includes needed by libAfterImage--------------------
if(asimage)
  set(ASEXTRA_LIBRARIES)
  find_Package(GIF)
  if(GIF_FOUND)
    set(ASEXTRA_LIBRARIES ${ASEXTRA_LIBRARIES} ${GIF_LIBRARIES})
  endif()
  find_Package(TIFF)
  if(TIFF_FOUND)
    set(ASEXTRA_LIBRARIES ${ASEXTRA_LIBRARIES} ${TIFF_LIBRARIES})
  endif()
  find_Package(PNG)
  if(PNG_FOUND)
    set(ASEXTRA_LIBRARIES ${ASEXTRA_LIBRARIES} ${PNG_LIBRARIES})
  endif()
  find_Package(JPEG)
  if(JPEG_FOUND)
    set(ASEXTRA_LIBRARIES ${ASEXTRA_LIBRARIES} ${JPEG_LIBRARIES})
  endif()
endif()

#---Check for GSL library---------------------------------------------------------------
if(mathmore)
  message(STATUS "Looking for GSL")
  if(NOT builtin_gsl)
    find_package(GSL 1.10)
    if(NOT GSL_FOUND)
      message(STATUS "GSL not found. Set variable GSL_DIR to point to your GSL installation")
      message(STATUS "               Alternatively, you can also enable the option 'builtin_gsl' to build the GSL libraries internally'") 
      message(STATUS "               For the time being switching OFF 'mathmore' option")
      set(mathmore OFF CACHE BOOL "" FORCE)
    endif()
  else()
    set(gsl_version 1.15)
    message(STATUS "Downloading and building GSL version ${gsl_version}") 
    ExternalProject_Add(
      GSL
      URL http://mirror.switch.ch/ftp/mirror/gnu/gsl/gsl-${gsl_version}.tar.gz
      INSTALL_DIR ${CMAKE_BINARY_DIR}
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix <INSTALL_DIR>
    )
    set(GSL_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
    set(GSL_LIBRARIES -L${CMAKE_BINARY_DIR}/lib -lgsl -lgslcblas -lm)
  endif()
endif()


#---Check for Python installation-------------------------------------------------------
if(python)
  message(STATUS "Looking for Python")
  #---First look for the python interpreter and fix the version of it for the libraries--
  find_package(PythonInterp)
  if(PYTHONINTERP_FOUND)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;sys.stdout.write(sys.version[:3])"
                    OUTPUT_VARIABLE PYTHON_VERSION)
    message(STATUS "Found Python interpreter version ${PYTHON_VERSION}")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;sys.stdout.write(sys.prefix)"
                    OUTPUT_VARIABLE PYTHON_PREFIX)
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PYTHON_PREFIX})
  endif()
  find_package(PythonLibs)
  if(NOT PYTHONLIBS_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "PythonLibs package not found and python component required")
    else()
      set(python OFF CACHE BOOL "" FORCE)
      message(STATUS "Python not found. Switching off python option")
    endif()
  else()
  endif()
endif()

#---Check for Ruby installation-------------------------------------------------------
if(ruby)
  message(STATUS "Looking for Ruby")
  find_package(Ruby)
  if(NOT RUBY_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Ruby package not found and ruby component required")
    else()
      set(ruby OFF CACHE BOOL "" FORCE)
      message(STATUS "Ruby not found. Switching off ruby option")
    endif()
  else()
    string(REGEX REPLACE "([0-9]+).*$" "\\1" RUBY_MAJOR_VERSION "${RUBY_VERSION}")
    string(REGEX REPLACE "[0-9]+\\.([0-9]+).*$" "\\1" RUBY_MINOR_VERSION "${RUBY_VERSION}")
    string(REGEX REPLACE "[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" RUBY_PATCH_VERSION "${RUBY_VERSION}")
  endif()
endif()

#---Check for OpenGL installation-------------------------------------------------------
if(opengl)
  message(STATUS "Looking for OpenGL")
  if(APPLE AND NOT cocoa)
    find_path(OPENGL_INCLUDE_DIR GL/gl.h  PATHS /usr/X11R6/include)
    find_library(OPENGL_gl_LIBRARY NAMES GL PATHS /usr/X11R6/lib)
    find_library(OPENGL_glu_LIBRARY NAMES GLU PATHS /usr/X11R6/lib)
    find_package_handle_standard_args(OpenGL REQUIRED_VARS OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY OPENGL_glu_LIBRARY)
    set(OPENGL_LIBRARIES ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
    mark_as_advanced(OPENGL_INCLUDE_DIR OPENGL_glu_LIBRARY OPENGL_gl_LIBRARY)
  else()
    find_package(OpenGL)
  endif()
  if(NOT OPENGL_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "OpenGL package not found and opengl option required")
    else()
      message(STATUS "OpenGL not found. Switching off opengl option")
      set(opengl OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Graphviz installation-------------------------------------------------------
if(gviz)
  message(STATUS "Looking for Graphviz")
  find_package(Graphviz)
  if(NOT GRAPHVIZ_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Graphviz package not found and gviz option required")
    else()
      message(STATUS "Graphviz not found. Switching off gviz option")
      set(gviz OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Qt installation-------------------------------------------------------
if(qt OR qtgsi)
  message(STATUS "Looking for Qt4")
  find_package(Qt4 4.8 COMPONENTS QtCore QtGui)
  if(NOT QT4_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Qt4 package not found and qt/qtgsi component required")
    else()
      message(STATUS "Qt4 not found. Switching off qt/qtgsi option")
      set(qt OFF CACHE BOOL "" FORCE)
      set(qtgsi OFF CACHE BOOL "" FORCE)
    endif()
  else()
    MATH(EXPR QT_VERSION_NUM "${QT_VERSION_MAJOR}*10000 + ${QT_VERSION_MINOR}*100 + ${QT_VERSION_PATCH}")
  endif()
endif()


#---Check for Bonjour installation-------------------------------------------------------
if(bonjour)
  message(STATUS "Looking for Bonjour")
  find_package(Bonjour)
  if(NOT BONJOUR_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Bonjour/Avahi libraries not found and Bonjour component required")
    else()
      message(STATUS "Bonjour not found. Switching off bonjour option")
      set(bonjour OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()


#---Check for krb5 Support-----------------------------------------------------------
if(krb5)
  message(STATUS "Looking for Kerberos 5")
  find_package(Kerberos5)
  if(NOT KRB5_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Kerberos 5 libraries not found and they are required")
    else()
      message(STATUS "Kerberos 5 not found. Switching off krb5 option")
      set(krb5 OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

if(krb5 OR afs)
  find_library(COMERR_LIBRARY com_err)
  if(COMERR_LIBRARY)
    set(COMERR_LIBRARIES ${COMERR_LIBRARY})
  endif()
  mark_as_advanced(COMERR_LIBRARY)
endif()

#---Check for XML Parser Support-----------------------------------------------------------
if(xml)
  message(STATUS "Looking for LibXml2")
  find_package(LibXml2)
  if(NOT LIBXML2_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "LibXml2 libraries not found and they are required (xml option enabled)")
    else()
      message(STATUS "LibXml2 not found. Switching off xml option")
      set(xml OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for OpenSSL------------------------------------------------------------------
if(ssl)
  message(STATUS "Looking for OpenSSL")
  find_package(OpenSSL)
  if(NOT OPENSSL_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "OpenSSL libraries not found and they are required (ssl option enabled)")
    else()
      message(STATUS "OpenSSL not found. Switching off ssl option")
      set(ssl OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Castor-------------------------------------------------------------------
if(castor OR rfio)
  message(STATUS "Looking for Castor")
  find_package(Castor)
  if(NOT CASTOR_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Castor libraries not found and they are required (castor option enabled)")
    else()
      message(STATUS "Castor not found. Switching off castor/rfio option")
      set(castor OFF CACHE BOOL "" FORCE)
      set(rfio OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for MySQL-------------------------------------------------------------------
if(mysql)
  message(STATUS "Looking for MySQL")
  find_package(MySQL)
  if(NOT MYSQL_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "MySQL libraries not found and they are required (mysql option enabled)")
    else()
      message(STATUS "MySQL not found. Switching off mysql option")
      set(mysql OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Oracle-------------------------------------------------------------------
if(oracle)
  message(STATUS "Looking for Oracle")
  find_package(Oracle)
  if(NOT ORACLE_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Oracle libraries not found and they are required (orable option enabled)")
    else()
      message(STATUS "Oracle not found. Switching off oracle option")
      set(oracle OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for ODBC-------------------------------------------------------------------
if(odbc)
  message(STATUS "Looking for ODBC")
  find_package(ODBC)
  if(NOT ODBC_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "ODBC libraries not found and they are required (odbc option enabled)")
    else()
      message(STATUS "ODBC not found. Switching off odbc option")
      set(odbc OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for PostgreSQL-------------------------------------------------------------------
if(pgsql)
  message(STATUS "Looking for PostgreSQL")
  find_package(PostgreSQL)
  if(NOT POSTGRESQL_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "PostgreSQL libraries not found and they are required (pgsql option enabled)")
    else()
      message(STATUS "PostgreSQL not found. Switching off pgsql option")
      set(pgsql OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for SQLite-------------------------------------------------------------------
if(sqlite)
  message(STATUS "Looking for SQLite")
  find_package(Sqlite)
  if(NOT SQLITE_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "SQLite libraries not found and they are required (sqlite option enabled)")
    else()
      message(STATUS "SQLite not found. Switching off sqlite option")
      set(sqlite OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Pythia6-------------------------------------------------------------------
if(pythia6)
  message(STATUS "Looking for Pythia6")
  find_package(Pythia6 QUIET)
  if(NOT PYTHIA6_FOUND AND NOT pythia6_nolink)
    if(fail-on-missing)
      message(FATAL_ERROR "Pythia6 libraries not found and they are required (pythia6 option enabled)")
    else()
      message(STATUS "Pythia6 not found. Switching off pythia6 option")
      set(pythia6 OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Pythia8-------------------------------------------------------------------
if(pythia8)
  message(STATUS "Looking for Pythia8")
  find_package(Pythia8)
  if(NOT PYTHIA8_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Pythia8 libraries not found and they are required (pythia8 option enabled)")
    else()
      message(STATUS "Pythia8 not found. Switching off pythia8 option")
      set(pythia8 OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for FFTW3-------------------------------------------------------------------
if(fftw3)
  message(STATUS "Looking for FFTW3")
  find_package(FFTW)
  if(NOT FFTW_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "FFTW3 libraries not found and they are required (fftw3 option enabled)")
    else()
      message(STATUS "FFTW3 not found. Switching off fftw3 option")
      set(fftw3 OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for fitsio-------------------------------------------------------------------
if(fitsio)
  if(builtin_cfitsio)
    set(cfitsio_version 3.280)
    string(REPLACE "." "" cfitsio_version_no_dots ${cfitsio_version})
    message(STATUS "Downloading and building CFITSIO version ${cfitsio_version}") 
    ExternalProject_Add(
      CFITSIO
      URL ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio${cfitsio_version_no_dots}.tar.gz 
      INSTALL_DIR ${CMAKE_BINARY_DIR}
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix <INSTALL_DIR>
      BUILD_IN_SOURCE 1
    )
    set(CFITSIO_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
    set(CFITSIO_LIBRARIES -L${CMAKE_BINARY_DIR}/lib -lcfitsio)
  else()
    message(STATUS "Looking for CFITSIO")  
    find_package(CFITSIO)
    if(NOT CFITSIO_FOUND)
      message(STATUS "CFITSIO not found. You can enable the option 'builtin_cfitsio' to build the library internally'") 
      message(STATUS "                   For the time being switching off 'fitsio' option")
      set(fitsio OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()


#---Check Shadow password support----------------------------------------------------
if(shadowpw)
  if(NOT EXISTS /etc/shadow)  #---TODO--The test always succeeds because the actual file is protected
    if(NOT CMAKE_SYSTEM_NAME MATCHES Linux)
      message(STATUS "Support Shadow password not found. Switching off shadowpw option")
      set(shadowpw OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Alien support----------------------------------------------------------------
if(alien)
  find_package(Alien)
  if(NOT ALIEN_FOUND)
    message(STATUS "Alien API not found. Set variable ALIEN_DIR to point to your Alien installation")
    message(STATUS "For the time being switching OFF 'alien' option")
    set(alien OFF CACHE BOOL "" FORCE)
  endif()
endif()

#---Monalisa support----------------------------------------------------------------
if(monalisa)
  find_package(Monalisa)
  if(NOT MONALISA_FOUND)
    message(STATUS "Monalisa not found. Set variable MONALISA_DIR to point to your Monalisa installation")
    message(STATUS "For the time being switching OFF 'monalisa' option")
    set(monalisa OFF CACHE BOOL "" FORCE)
  endif()
endif()

#---Check for Xrootd support---------------------------------------------------------
if(xrootd)
  message(STATUS "Looking for XROOTD")
  if(NOT builtin_xrootd)
    find_package(XROOTD)
    if(NOT XROOTD_FOUND)
      message(STATUS "XROOTD not found. Set environment variable XRDSYS to point to your XROOTD installation")
      message(STATUS "                  Alternatively, you can also enable the option 'builtin_xrootd' to build XROOTD  internally'")
      message(STATUS "                  For the time being switching OFF 'xrootd' option")
      set(xrootd OFF CACHE BOOL "" FORCE)
    else()
      set(xrootd_versionnum ${xrdversnum})  # variable used internally
    endif()
  endif()
endif()
if(builtin_xrootd)
  set(xrootd_version 3.3.6)
  set(xrootd_versionnum 300030006)
  message(STATUS "Downloading and building XROOTD version ${xrootd_version}") 
  string(REPLACE "-Wall " "" __cxxflags "${CMAKE_CXX_FLAGS}")                        # Otherwise it produces tones of warnings
  string(REPLACE "-W " "" __cxxflags "${__cxxflags} -Wno-deprecated-declarations -Wno-duplicate-decl-specifier")
  ExternalProject_Add(
    XROOTD
    URL http://xrootd.org/download/v${xrootd_version}/xrootd-${xrootd_version}.tar.gz
    INSTALL_DIR ${CMAKE_BINARY_DIR}
    CMAKE_ARGS -DENABLE_PERL=FALSE
               -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_CXX_FLAGS=${__cxxflags}
  )
  # We cannot call find_package(XROOTD) becuase the package is not yet built. So, we need to emulate what it defines....
  set(_LIBDIR_DEFAULT "lib")
  if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND NOT CMAKE_CROSSCOMPILING AND NOT EXISTS "/etc/debian_version")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
      set(_LIBDIR_DEFAULT "lib64")
    endif()
  endif()
  set(XROOTD_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include/xrootd ${CMAKE_BINARY_DIR}/include/xrootd/private)
  set(XROOTD_LIBRARIES ${CMAKE_BINARY_DIR}/${_LIBDIR_DEFAULT}/libXrdMain${CMAKE_SHARED_LIBRARY_SUFFIX}
                       ${CMAKE_BINARY_DIR}/${_LIBDIR_DEFAULT}/libXrdUtils${CMAKE_SHARED_LIBRARY_SUFFIX}
                       ${CMAKE_BINARY_DIR}/${_LIBDIR_DEFAULT}/libXrdClient${CMAKE_SHARED_LIBRARY_SUFFIX}
                       ${CMAKE_BINARY_DIR}/${_LIBDIR_DEFAULT}/libXrdCl${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(XROOTD_CFLAGS "-DROOTXRDVERS=${xrootd_versionnum}")
  set(xrootd ON CACHE BOOL "" FORCE)
endif()
if(xrootd AND xrootd_versionnum VERSION_GREATER 300030005)
  set(netxng ON)
else()
  set(netxng OFF)
endif()

#---Check for cling and llvm ----------------------------------------------------------------
find_library(CMAKE_TINFO_LIBS NAMES tinfo ncurses)
mark_as_advanced(CMAKE_TINFO_LIBS)
if(NOT CMAKE_TINFO_LIBS)
  set(CMAKE_TINFO_LIBS "")   #  Often if not found is still OK
endif()
if(cling)
  if(NOT DEFINED LLVM_BUILD_TYPE)
    set(LLVM_BUILD_TYPE Release)
  endif()
  if(builtin_llvm)
    set(LLVM_SOURCE_DIR ${CMAKE_SOURCE_DIR}/interpreter/llvm/src)
    set(LLVM_INSTALL_DIR ${CMAKE_BINARY_DIR}/LLVM-install)
    ExternalProject_Add(
      LLVM
      PREFIX LLVM
      SOURCE_DIR ${LLVM_SOURCE_DIR}
      INSTALL_DIR ${CMAKE_BINARY_DIR}/LLVM-install
      CMAKE_ARGS -DLLVM_INCLUDE_TESTS=OFF
                 -DLLVM_INCLUDE_EXAMPLES=OFF
                 -DLLVM_BUILD_TOOLS=OFF
                 -DLLVM_TARGETS_TO_BUILD=X86
                 -DLLVM_FORCE_USE_OLD_TOOLCHAIN=ON
                 -DCMAKE_BUILD_TYPE=${LLVM_BUILD_TYPE}
                 -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                 -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                 -DLLVM_NOCLING=ON
      INSTALL_COMMAND ${CMAKE_COMMAND} -P <BINARY_DIR>/cmake_install.cmake
      LOG_INSTALL 1
    )
    ExternalProject_Get_Property(LLVM STAMP_DIR BINARY_DIR)
    file(REMOVE ${STAMP_DIR}/LLVM-build)     # Deleting stamp file to force a configure

    #---The list of libraires is optatined by runnning 'llvm-config --libs'
    set(LLVM_INCLUDE_DIR ${CMAKE_BINARY_DIR}/LLVM-install/include)
    set(LLVM_LIBRARIES ${BINARY_DIR}/lib/libclangFrontend.a
                       ${BINARY_DIR}/lib/libclangSerialization.a
                       ${BINARY_DIR}/lib/libclangDriver.a
                       ${BINARY_DIR}/lib/libclangCodeGen.a
                       ${BINARY_DIR}/lib/libclangParse.a
                       ${BINARY_DIR}/lib/libclangSema.a
                       ${BINARY_DIR}/lib/libclangAnalysis.a
                       ${BINARY_DIR}/lib/libclangRewriteCore.a
                       ${BINARY_DIR}/lib/libclangAST.a
                       ${BINARY_DIR}/lib/libclangBasic.a
                       ${BINARY_DIR}/lib/libclangEdit.a
                       ${BINARY_DIR}/lib/libclangLex.a
                       ${BINARY_DIR}/lib/libLLVMInstrumentation.a
                       ${BINARY_DIR}/lib/libLLVMIRReader.a
                       ${BINARY_DIR}/lib/libLLVMAsmParser.a
                       ${BINARY_DIR}/lib/libLLVMDebugInfo.a
                       ${BINARY_DIR}/lib/libLLVMOption.a
                       ${BINARY_DIR}/lib/libLLVMLTO.a
                       ${BINARY_DIR}/lib/libLLVMLinker.a
                       ${BINARY_DIR}/lib/libLLVMipo.a
                       ${BINARY_DIR}/lib/libLLVMVectorize.a
                       ${BINARY_DIR}/lib/libLLVMBitWriter.a
                       ${BINARY_DIR}/lib/libLLVMTableGen.a
                       ${BINARY_DIR}/lib/libLLVMX86Disassembler.a
                       ${BINARY_DIR}/lib/libLLVMX86AsmParser.a
                       ${BINARY_DIR}/lib/libLLVMX86CodeGen.a
                       ${BINARY_DIR}/lib/libLLVMX86Desc.a
                       ${BINARY_DIR}/lib/libLLVMX86Info.a
                       ${BINARY_DIR}/lib/libLLVMX86AsmPrinter.a
                       ${BINARY_DIR}/lib/libLLVMX86Utils.a
                       ${BINARY_DIR}/lib/libLLVMSelectionDAG.a
                       ${BINARY_DIR}/lib/libLLVMAsmPrinter.a
                       ${BINARY_DIR}/lib/libLLVMMCDisassembler.a
                       ${BINARY_DIR}/lib/libLLVMMCParser.a
                       ${BINARY_DIR}/lib/libLLVMInterpreter.a
                       ${BINARY_DIR}/lib/libLLVMMCJIT.a
                       ${BINARY_DIR}/lib/libLLVMJIT.a
                       ${BINARY_DIR}/lib/libLLVMCodeGen.a
                       ${BINARY_DIR}/lib/libLLVMObjCARCOpts.a
                       ${BINARY_DIR}/lib/libLLVMScalarOpts.a
                       ${BINARY_DIR}/lib/libLLVMInstCombine.a
                       ${BINARY_DIR}/lib/libLLVMTransformUtils.a
                       ${BINARY_DIR}/lib/libLLVMipa.a
                       ${BINARY_DIR}/lib/libLLVMAnalysis.a
                       ${BINARY_DIR}/lib/libLLVMRuntimeDyld.a
                       ${BINARY_DIR}/lib/libLLVMExecutionEngine.a
                       ${BINARY_DIR}/lib/libLLVMTarget.a
                       ${BINARY_DIR}/lib/libLLVMMC.a
                       ${BINARY_DIR}/lib/libLLVMObject.a
                       ${BINARY_DIR}/lib/libLLVMBitReader.a
                       ${BINARY_DIR}/lib/libLLVMCore.a
                       ${BINARY_DIR}/lib/libLLVMSupport.a )
    file(READ ${LLVM_SOURCE_DIR}/configure _filestr)
    string(REGEX REPLACE ".*PACKAGE_VERSION='([0-9]+[.][0-9]+).*" "\\1" LLVM_VERSION ${_filestr})
  else()
    find_package(LLVM REQUIRED)
  endif()

  ExternalProject_Add(
    CLING
    PREFIX CLING
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/interpreter/cling
    INSTALL_DIR ${CMAKE_BINARY_DIR}/CLING-install
    CMAKE_ARGS -DCLING_PATH_TO_LLVM_SOURCE=${CMAKE_SOURCE_DIR}/interpreter/llvm/src
               -DCLING_PATH_TO_LLVM_BUILD=${CMAKE_BINARY_DIR}/LLVM-install
               -DLLVM_FORCE_USE_OLD_TOOLCHAIN=ON
               -DCMAKE_BUILD_TYPE=${LLVM_BUILD_TYPE}
               -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    INSTALL_COMMAND ${CMAKE_COMMAND} -P <BINARY_DIR>/cmake_install.cmake
    LOG_INSTALL 1
  )
  ExternalProject_Get_Property(CLING STAMP_DIR BINARY_DIR)
  file(REMOVE ${STAMP_DIR}/CLING-build)     # Deleting stamp file to force a configure
  set(CLING_INCLUDE_DIR ${CMAKE_BINARY_DIR}/CLING-install/include)
  set(CLING_LIBRARIES ${BINARY_DIR}/lib/libclingInterpreter.a
                      ${BINARY_DIR}/lib/libclingMetaProcessor.a
                      ${BINARY_DIR}/lib/libclingUtils.a
                      ${LLVM_LIBRARIES})
  #--Additional flags obtained from llvm-config --cxxflags
  set(CLING_CXXFLAGS "-fvisibility-inlines-hidden -fno-strict-aliasing -Wno-unused-parameter -Wwrite-strings -Wno-long-long -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS")
  if (CMAKE_COMPILER_IS_GNUCXX)
    set(CLING_CXXFLAGS "${CLING_CXXFLAGS} -Wno-missing-field-initializers")
  endif()
  add_dependencies(CLING LLVM)
endif()


#---Check for gfal-------------------------------------------------------------------
if(gfal)
  find_package(GFAL)
  if(NOT GFAL_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "Gfal library not found and is required (gfal option enabled)")
    else()
      message(STATUS "GFAL library not found. Set variable GFAL_DIR to point to your gfal installation
                      and the variable SRM_IFCE_DIR to the srm_ifce installation")
      message(STATUS "For the time being switching OFF 'gfal' option")
      set(gfal OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()


#---Check for dCache-------------------------------------------------------------------
if(dcache)
  find_package(DCAP)
  if(NOT DCAP_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "dCap library not found and is required (dcache option enabled)")
    else()
      message(STATUS "dCap library not found. Set variable DCAP_DIR to point to your dCache installation")
      message(STATUS "For the time being switching OFF 'dcache' option")
      set(dcache OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for Ldap--------------------------------------------------------------------
if(ldap)
  find_package(Ldap)
  if(NOT LDAP_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "ldap library not found and is required (ldap option enabled)")
    else()
      message(STATUS "ldap library not found. Set variable LDAP_DIR to point to your ldap installation")
      message(STATUS "For the time being switching OFF 'ldap' option")
      set(ldap OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for globus--------------------------------------------------------------------
if(globus)
  find_package(Globus)
  if(NOT GLOBUS_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "globus libraries not found and is required ('globus' option enabled)")
    else()
      message(STATUS "globus libraries not found. Set environment var GLOBUS_LOCATION or varibale GLOBUS_DIR to point to your globus installation")
      message(STATUS "For the time being switching OFF 'globus' option")
      set(globus OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

#---Check for ftgl if needed----------------------------------------------------------
if(NOT builtin_ftgl)
  find_package(FTGL)
  if(NOT FTGL_FOUND)
    if(fail-on-missing)
      message(FATAL_ERROR "ftgl library not found and is required ('builtin_ftgl' is OFF). Set varible FTGL_ROOT_DIR to installation location")
    else()
      message(STATUS "ftgl library not found. Set variable FTGL_ROOT_DIR to point to your installation")
      message(STATUS "For the time being switching ON 'builtin_ftgl' option")
      set(builtin_ftgl ON CACHE BOOL "" FORCE)
    endif()
  endif()
endif()
if(builtin_ftgl)
  set(FTGL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/graf3d/ftgl/inc)
  set(FTGL_CFLAGS -DBUILTIN_FTGL)
  set(FTGL_LIBRARIES FTGL)
endif()

#---Check for DavIx library-----------------------------------------------------------
if(davix)
  set(DAVIX_VERSION 0.2.7)
  message(STATUS "Downloading and building Davix version ${DAVIX_VERSION}")
    ExternalProject_Add(
      DAVIX
      URL http://grid-deployment.web.cern.ch/grid-deployment/dms/lcgutil/tar/davix/davix-${DAVIX_VERSION}.tar.gz
      INSTALL_DIR ${CMAKE_BINARY_DIR}
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DLIB_SUFFIX=
    )
    set(DAVIX_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include/davix)
    set(DAVIX_LIBRARY -L${CMAKE_BINARY_DIR}/lib -ldavix)
    set(DAVIX_INCLUDE_DIRS ${DAVIX_INCLUDE_DIR})
    set(DAVIX_LIBRARIES ${DAVIX_LIBRARY})
endif()


#---Report non implemented options---------------------------------------------------
foreach(opt afs chirp clarens glite hdfs pch peac sapdb srp geocad)
  if(${opt})
    message(STATUS ">>> Option '${opt}' not implemented yet! Signal your urgency to pere.mato@cern.ch")
    set(${opt} OFF CACHE BOOL "" FORCE)
  endif()
endforeach()

