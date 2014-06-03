#---Define a function to do not polute the top level namespace with unneeded variables-----------------------
function(RootConfigure)
  
#---Define all sort of variables to bridge between the old Module.mk and the new CMake equivalents-----------
set(valueON yes)
set(valueOFF no)

set(ROOT_DICTTYPE cint)
set(ROOT_CONFIGARGS "")
set(top_srcdir ${CMAKE_SOURCE_DIR})
set(top_builddir ${CMAKE_BINARY_DIR})
set(architecture ${ROOT_ARCHITECTURE})
set(platform ${ROOT_PLATFORM})
set(host)
set(useconfig FALSE)
set(major ${ROOT_MAJOR_VERSION})
set(minor ${ROOT_MINOR_VERSION})
set(revis ${ROOT_PATCH_VERSION})
set(mkliboption "-v ${major} ${minor} ${revis} ")
set(cflags ${CMAKE_CXX_FLAGS})
set(ldflags ${CMAKE_CXX_LINK_FLAGS})

set(winrtdebug ${value${winrtdebug}})
set(exceptions ${value${exceptions}})
set(explicitlink ${value${explicitlink}})

if(gnuinstall)
  set(prefix ${CMAKE_INSTALL_PREFIX})
  set(etcdir ${CMAKE_INSTALL_FULL_SYSCONFDIR})
  set(bindir ${CMAKE_INSTALL_FULL_BINDIR})
  set(libdir ${CMAKE_INSTALL_FULL_LIBDIR})
  set(incdir ${CMAKE_INSTALL_FULL_INCLUDEDIR})
  set(mandir ${CMAKE_INSTALL_FULL_MANDIR})
  set(plugindir ${CMAKE_INSTALL_FULL_SYSCONFDIR}/plugins)
  set(datadir ${CMAKE_INSTALL_FULL_DATADIR})
  set(elispdir ${CMAKE_INSTALL_FULL_ELISPDIR})
  set(ttffontdir ${CMAKE_INSTALL_FULL_FONTDIR})
  set(macrodir ${CMAKE_INSTALL_FULL_MACRODIR})
  set(srcdir ${CMAKE_INSTALL_FULL_SRCDIR})
  set(iconpath ${CMAKE_INSTALL_FULL_ICONDIR})
  set(cintincdir ${CMAKE_INSTALL_FULL_CINTINCDIR})
  set(docdir ${CMAKE_INSTALL_FULL_DOCDIR})
  set(testdir ${CMAKE_INSTALL_FULL_TESTDIR})
  set(tutdir ${CMAKE_INSTALL_FULL_TUTDIR})
  set(aclocaldir ${CMAKE_INSTALL_FULL_ACLOCALDIR})
else()
  set(prefix $(ROOTSYS))
  set(etcdir ${prefix}/${CMAKE_INSTALL_SYSCONFDIR})
  set(bindir ${prefix}/${CMAKE_INSTALL_BINDIR})
  set(libdir ${prefix}/${CMAKE_INSTALL_LIBDIR})
  set(incdir ${prefix}/${CMAKE_INSTALL_INCLUDEDIR})
  set(mandir ${prefix}/${CMAKE_INSTALL_MANDIR})
  set(plugindir ${prefix}/${CMAKE_INSTALL_SYSCONFDIR}/plugins)
  set(datadir ${prefix}/${CMAKE_INSTALL_DATADIR})
  set(elispdir ${prefix}/${CMAKE_INSTALL_ELISPDIR})
  set(ttffontdir ${prefix}/${CMAKE_INSTALL_FONTDIR})
  set(macrodir ${prefix}/${CMAKE_INSTALL_MACRODIR})
  set(srcdir ${prefix}/${CMAKE_INSTALL_SRCDIR})
  set(iconpath ${prefix}/${CMAKE_INSTALL_ICONDIR})
  set(cintincdir ${prefix}/${CMAKE_INSTALL_CINTINCDIR})
  set(docdir ${prefix}/${CMAKE_INSTALL_DOCDIR})
  set(testdir ${prefix}/${CMAKE_INSTALL_TESTDIR})
  set(tutdir ${prefix}/${CMAKE_INSTALL_TUTDIR})
  set(aclocaldir ${prefix}/${CMAKE_INSTALL_ACLOCALDIR})
endif()

set(LibSuffix ${SOEXT})

set(buildx11 ${value${x11}})
set(x11libdir -L${X11_LIBRARY_DIR})
set(xpmlibdir -L${X11_LIBRARY_DIR})
set(xpmlib ${X11_Xpm_LIB})
set(enable_xft ${value${xft}})

set(enable_thread yes)
set(threadflag)
set(threadlibdir)
set(threadlib ${CMAKE_THREAD_LIBS_INIT})

set(builtinfreetype ${value${builtin_freetype}})
set(builtinpcre ${value${builtin_pcre}})

set(builtinzlib ${value${builtin_zlib}})
set(zliblibdir ${ZLIB_LIBRARY_DIR})
set(zliblib ${ZLIB_LIBRARY})
set(zlibincdir ${ZLIB_INCLUDE_DIR})

set(buildgl ${value${opengl}})
set(opengllibdir ${OPENGL_LIBRARY_DIR})
set(openglulib ${OPENGL_glu_LIBRARY})
set(opengllib ${OPENGL_gl_LIBRARY})
set(openglincdir ${OPENGL_INCLUDE_DIR})

set(buildldap ${value${ldap}})
set(ldaplibdir ${LDAP_LIBRARY_DIR})
set(ldaplib ${LDAP_LIBRARY})
set(ldapincdir ${LDAP_INCLUDE_DIR})

set(buildmysql ${value${mysql}})
set(mysqllibdir ${MYSQL_LIBRARY_DIR})
set(mysqllib ${MYSQL_LIBRARY})
set(mysqlincdir ${MYSQL_INCLUDE_DIR})

set(buildoracle ${value${oracle}})
set(oraclelibdir ${ORACLE_LIBRARY_DIR})
set(oraclelib ${ORACLE_LIBRARY})
set(oracleincdir ${ORACLE_INCLUDE_DIR})

set(buildpgsql ${value${pgsql}})
set(pgsqllibdir ${PGSQL_LIBRARY_DIR})
set(pgsqllib ${PGSQL_LIBRARY})
set(pgsqlincdir ${PGSQL_INCLUDE_DIR})

set(buildsqlite ${value${sqlite}})
set(sqlitelibdir ${SQLITE_LIBRARY_DIR})
set(sqlitelib ${SQLITE_LIBRARY})
set(sqliteincdir ${SQLITE_INCLUDE_DIR})

set(buildsapdb ${value${sapdb}})
set(sapdblibdir ${SAPDB_LIBRARY_DIR})
set(sapdblib ${SAPDB_LIBRARY})
set(sapdbincdir ${SAPDB_INCLUDE_DIR})

set(buildodbc ${value${odbc}})
set(odbclibdir ${OCDB_LIBRARY_DIR})
set(odbclib ${OCDB_LIBRARY})
set(odbcincdir ${OCDB_INCLUDE_DIR})

set(buildqt ${value${qt}})
set(buildqtpsi ${value${qtgsi}})
set(qtlibdir ${QT_LIBRARY_DIR})
set(qtlib ${QT_QT_LIBRARY})
set(qtincdir ${QT_INCLUDE_DIR})
set(qtvers ${QT_VERSION_MAJOR})
set(qtmocexe ${QT_MOC_EXECUTABLE})

set(buildrfio ${value${rfio}})
set(shiftlibdir ${RFIO_LIBRARY_DIR})
set(shiftlib ${RFIO_LIBRARY})
set(shiftincdir ${RFIO_INCLUDE_DIR})
set(shiftcflags)

set(buildcastor ${value${castor}})
set(castorlibdir ${CASTOR_LIBRARY_DIR})
set(castorlib ${CASTOR_LIBRARY})
set(castorincdir ${CASTOR_INCLUDE_DIR})
set(castorcflags)


set(builddavix ${value${davix}})
set(davixlibdir ${DAVIX_LIBRARY_DIR})
set(davixlib ${DAVIX_LIBRARY})
set(davixincdir ${DAVIX_INCLUDE_DIR})

set(buildnetxng ${value${netxng}})
if(netxng)
  set(useoldnetx no)
else()
  set(useoldnetx yes)
endif()

set(builddcap ${value${dcap}})
set(dcaplibdir ${DCAP_LIBRARY_DIR})
set(dcaplib ${DCAP_LIBRARY})
set(dcapincdir ${DCAP_INCLUDE_DIR})

set(buildftgl ${value${builtin_ftgl}})
set(ftgllibdir ${FTGL_LIBRARY_DIR})
set(ftgllibs ${FTGL_LIBRARIES})
set(ftglincdir ${FTGL_INCLUDE_DIR})

set(buildglew ${value${builtin_glew}})
set(glewlibdir ${GLEW_LIBRARY_DIR})
set(glewlibs ${GLEW_LIBRARIES})
set(glewincdir ${GLEW_INCLUDE_DIR})

set(buildgfal ${value${gfal}})
set(gfallibdir ${GFAL_LIBRARY_DIR})
set(gfallib ${GFAL_LIBRARY})
set(gfalincdir ${GFAL_INCLUDE_DIR})

set(buildglite ${value${glite}})
set(glitelibdir ${GLITE_LIBRARY_DIR})
set(glitelib ${GLITE_LIBRARY})
set(gaw_cppflags)

set(buildmemstat ${value${memstat}})

set(buildbonjour ${value${bonjour}})
set(dnssdlibdir ${BONJOUR_LIBRARY_DIR})
set(dnssdlib ${BONJOUR_LIBRARY})
set(dnsdincdir ${BONJOUR_INCLUDE_DIR})
set(bonjourcppflags)

set(buildchirp ${value${chirp}})
set(chirplibdir ${CHIRP_LIBRARY_DIR})
set(chirplib ${CHIRP_LIBRARY})
set(chirpincdir ${CHIRP_INCLUDE_DIR})

set(buildhdfs ${value${hdfs}})
set(hdfslibdir ${HDFS_LIBRARY_DIR})
set(hdfslib ${HDFS_LIBRARY})
set(hdfsincdir ${HDFS_INCLUDE_DIR})

set(jniincdir ${Java_INCLUDE_DIRS})
set(jvmlib ${Java_LIBRARIES})
set(jvmlibdir ${Java_LIBRARY_DIR})

set(buildalien ${value${alien}})
set(alienlibdir ${ALIEN_LIBRARY_DIR})
set(alienlib ${ALIEN_LIBRARY})
set(alienincdir ${ALIEN_INCLUDE_DIR})

set(buildasimage ${value${asimage}})
set(builtinafterimage ${builtin_afterimage})
set(asextralib ${ASEXTRA_LIBRARIES})
set(asextralibdir)
set(asjpegincdir ${JPEG_INCLUDE_DIR})
set(aspngincdir ${PNG_INCLUDE_DIR})
set(astiffincdir ${TIFF_INCLUDE_DIR})
set(asgifincdir ${GIF_INCLUDE_DIR})
set(asimageincdir)
set(asimagelib)
set(asimagelibdir)

set(buildpythia6 ${value${pythia6}})
set(pythia6libdir ${PYTHIA6_LIBRARY_DIR})
set(pythia6lib ${PYTHIA6_LIBRARY})
set(pythia6cppflags)
set(buildpythia8 ${value${pythia8}})
set(pythia8libdir ${PYTHIA8_LIBRARY_DIR})
set(pythia8lib ${PYTHIA8_LIBRARY})
set(pythia8cppflags)

set(buildfftw3 ${value${fftw3}})
set(fftw3libdir ${FFTW3_LIBRARY_DIR})
set(fftw3lib ${FFTW3_LIBRARY})
set(fftw3incdir ${FFTW3_INCLUDE_DIR})

set(buildfitsio ${value${fitsio}})
set(fitsiolibdir ${FITSIO_LIBRARY_DIR})
set(fitsiolib ${FITSIO_LIBRARY})
set(fitsioincdir ${FITSIO_INCLUDE_DIR})

set(buildgviz ${value${gviz}})
set(gvizlibdir ${GVIZ_LIBRARY_DIR})
set(gvizlib ${GVIZ_LIBRARY})
set(gvizincdir ${GVIZ_INCLUDE_DIR})
set(gvizcflags)

set(buildpython ${value${python}})
set(pythonlibdir ${PYTHON_LIBRARY_DIR})
set(pythonlib ${PYTHON_LIBRARY})
set(pythonincdir ${PYTHON_INCLUDE_DIR})
set(pythonlibflags)

set(buildruby ${value${ruby}})
set(rubylibdir ${RUBY_LIBRARY_DIR})
set(rubylib ${RUBY_LIBRARY})
set(rubyincdir ${RUBY_INCLUDE_DIR})

set(buildxml ${value${xml}})
set(xmllibdir ${LIBXML2_LIBRARY_DIR})
set(xmllib ${LIBXML2_LIBRARIES})
set(xmlincdir ${LIBXML2_INCLUDE_DIR})

set(buildxrd ${value${xrootd}})
set(xrdlibdir )
set(xrdincdir)
set(xrdaddopts)
set(extraxrdflags)
set(xrdversion)

set(srplibdir)
set(srplib)
set(srpincdir)

set(buildsrputil)
set(srputillibdir)
set(srputillib)
set(srputilincdir)

set(afslib ${AFS_LIBRARY})
set(afslibdir ${AFS_LIBRARY_DIR})
set(afsincdir ${AFS_INCLUDE_DIR})
set(afsextracflags)
set(afsshared)

set(alloclib)
set(alloclibdir)

set(buildkrb5 ${value${krb5}})
set(krb5libdir ${KRB5_LIBRARY_DIR})
set(krb5lib ${KRB5_LIBRARY})
set(krb5incdir ${KRB5_INCLUDE_DIR})
set(krb5init ${KRB5_INIT})

set(comerrlib)
set(comerrlibdir)
set(resolvlib)
set(cryptolib ${CRYPTLIBS})
set(cryptolibdir)

set(buildglobus ${value${globus}})
set(globuslibdir ${GLOBUS_LIBRARY_DIR})
set(globuslib ${GLOBUS_LIBRARY})
set(globusincdir ${GLOBUS_INCLUDE_DIR})
set(buildxrdgsi)

set(buildmonalisa ${value${monalisa}})
set(monalisalibdir ${MONALISA_LIBRARY_DIR})
set(monalisalib ${MONALISA_LIBRARY})
set(monalisaincdir ${MONALISA_INCLUDE_DIR})

set(ssllib ${OPENSSL_LIBRARIES})
set(ssllibdir)
set(sslincdir ${OPENSSL_INCLUDE_DIR})
set(sslshared)

set(gsllibs ${GSL_LIBRARIES})
set(gsllibdir)
set(gslincdir ${GSL_INCLUDE_DIR})
set(gslflags)

set(shadowpw ${value${shadowpw}})
set(buildgenvector ${value${genvector}})
set(buildmathmore ${value${mathmore}})
set(buildcling ${value${cling}})
set(buildroofit ${value${roofit}})
set(buildminuit2 ${value${minuit2}})
set(buildunuran ${value${unuran}})
set(buildgdml ${value${gdml}})
set(buildhttp ${value${http}})
set(buildtable ${value${table}})
set(buildtmva ${value${tmva}})

set(buildclarens ${value${clarens}})
set(clarensincdir ${CLARENS_INCLUDE_DIR})
set(clarenslibs ${CLARENS_LIBRARIES})
set(buildpeac ${value${peac}})


set(cursesincdir ${CURSES_INCLUDE_DIR})
set(curseslibdir)
set(curseslib ${CURSES_LIBRARIES})
set(curseshdr ${CURSES_HEADER_FILE})
set(buildeditline ${value${editline}})
set(cppunit)
set(dicttype ${ROOT_DICTTYPE})

#---RConfigure-------------------------------------------------------------------------------------------------
set(hasON define)
set(hasOFF undef)
set(hason define)
set(hasoff undef)
set(has1 define)
set(has0 undef)
set(has undef)
set(setresuid undef)
set(hasmathmore ${has${mathmore}})
set(haspthread ${has${CMAKE_USE_PTHREADS_INIT}})
set(hasxft ${has${xft}})
set(hascling ${has${cling}})
set(haslzmacompression ${has${lzma}})
set(hascocoa ${has${cocoa}})
set(hasvc ${has${vc}})
set(usec++11 ${has${cxx11}})
set(uselibc++ ${has${libcxx}})
set(hasllvm undef)
set(llvmdir /**/)
if(gcctoolchain)
  set(setgcctoolchain define)
else()
  set(setgcctoolchain undef)
endif()

#---root-config----------------------------------------------------------------------------------------------
ROOT_SHOW_OPTIONS(features)
string(REPLACE "c++11" "cxx11" features ${features}) # change the name of the c++11 feature needed for root-config.in
set(configfeatures ${features})
set(configargs ${ROOT_CONFIGARGS})
set(configoptions ${ROOT_CONFIGARGS})
get_filename_component(altcc ${CMAKE_C_COMPILER} NAME)
get_filename_component(altcxx ${CMAKE_CXX_COMPILER} NAME)
get_filename_component(altf77 "${CMAKE_Fortran_COMPILER}" NAME)
get_filename_component(altld ${CMAKE_CXX_COMPILER} NAME)

set(pythonvers ${PYTHON_VERSION})

#---RConfigure.h---------------------------------------------------------------------------------------------
configure_file(${PROJECT_SOURCE_DIR}/config/RConfigure.in include/RConfigure.h)
install(FILES ${CMAKE_BINARY_DIR}/include/RConfigure.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#---Configure and install various files----------------------------------------------------------------------
execute_Process(COMMAND hostname OUTPUT_VARIABLE BuildNodeInfo OUTPUT_STRIP_TRAILING_WHITESPACE )

configure_file(${CMAKE_SOURCE_DIR}/config/rootrc.in ${CMAKE_BINARY_DIR}/etc/system.rootrc @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/RConfigOptions.in include/RConfigOptions.h)
if(ruby)
  file(APPEND ${CMAKE_BINARY_DIR}/include/RConfigOptions.h "\#define R__RUBY_MAJOR ${RUBY_MAJOR_VERSION}\n\#define R__RUBY_MINOR ${RUBY_MINOR_VERSION}\n")
endif()

configure_file(${CMAKE_SOURCE_DIR}/config/Makefile-comp.in config/Makefile.comp)
configure_file(${CMAKE_SOURCE_DIR}/config/Makefile.in config/Makefile.config)
configure_file(${CMAKE_SOURCE_DIR}/config/mimes.unix.in ${CMAKE_BINARY_DIR}/etc/root.mimes)

#---Generate the ROOTConfig files to be used by CMake projects-----------------------------------------------
ROOT_SHOW_OPTIONS(ROOT_ENABLED_OPTIONS)
configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/ROOTConfig-version.cmake.in
               ${CMAKE_BINARY_DIR}/ROOTConfig-version.cmake @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/RootUseFile.cmake.in
               ${CMAKE_BINARY_DIR}/ROOTUseFile.cmake @ONLY)
#---To be used from the binary tree--------------------------------------------------------------------------
get_property(buildtree_include_dirs GLOBAL PROPERTY ROOT_INCLUDE_DIRS)
list(REMOVE_DUPLICATES buildtree_include_dirs)
set(ROOT_INCLUDE_DIR_SETUP "
# ROOT configured for use from the build tree - absolute paths are used.
set(ROOT_INCLUDE_DIRS ${buildtree_include_dirs})
")
set(ROOT_LIBRARY_DIR_SETUP "
# ROOT configured for use from the build tree - absolute paths are used.
set(ROOT_LIBRARY_DIR ${CMAKE_BINARY_DIR}/lib)
")
set(ROOT_BINARY_DIR_SETUP "
# ROOT configured for use from the build tree - absolute paths are used.
set(ROOT_BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
")
set(ROOT_MODULE_PATH_SETUP "
# ROOT configured for use CMake modules from source tree
set(CMAKE_MODULE_PATH \${CMAKE_MODULE_PATH} ${CMAKE_MODULE_PATH})
")

get_property(exported_targets GLOBAL PROPERTY ROOT_EXPORTED_TARGETS)
export(TARGETS ${exported_targets} FILE ${PROJECT_BINARY_DIR}/ROOTConfig-targets.cmake)
configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/ROOTConfig.cmake.in
               ${CMAKE_BINARY_DIR}/ROOTConfig.cmake @ONLY)

#---To be used from the install tree--------------------------------------------------------------------------
set(ROOT_INCLUDE_DIR_SETUP "
# ROOT configured for the install with relative paths, so use these
get_filename_component(ROOT_INCLUDE_DIRS \"\${_thisdir}/../include\" ABSOLUTE)
")
set(ROOT_LIBRARY_DIR_SETUP "
# ROOT configured for the install with relative paths, so use these
get_filename_component(ROOT_LIBRARY_DIR \"\${_thisdir}/../lib\" ABSOLUTE)
")
set(ROOT_BINARY_DIR_SETUP "
# ROOT configured for the install with relative paths, so use these
get_filename_component(ROOT_BINARY_DIR \"\${_thisdir}/../bin\" ABSOLUTE)
")
set(ROOT_MODULE_PATH_SETUP "
# ROOT configured for use CMake modules from installation tree
set(CMAKE_MODULE_PATH \${CMAKE_MODULE_PATH}  \${_thisdir}/modules)
")
configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/ROOTConfig.cmake.in
               ${CMAKE_BINARY_DIR}/installtree/ROOTConfig.cmake @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/ROOTConfig-version.cmake
              ${CMAKE_BINARY_DIR}/ROOTUseFile.cmake
              ${CMAKE_BINARY_DIR}/installtree/ROOTConfig.cmake DESTINATION ${CMAKE_INSTALL_CMAKEDIR})
install(EXPORT ${CMAKE_PROJECT_NAME}Exports FILE ROOTConfig-targets.cmake DESTINATION ${CMAKE_INSTALL_CMAKEDIR})


#---Especial definitions for root-config et al.--------------------------------------------------------------
if(prefix STREQUAL "$(ROOTSYS)")
  set(prefix $ROOTSYS)
  set(bindir $ROOTSYS/bin)
  set(libdir $ROOTSYS/lib)
  set(incdir $ROOTSYS/include)
  set(etcdir $ROOTSYS/etc)
  set(mandir $ROOTSYS/man/man1)
endif()


#---compiledata.h--------------------------------------------------------------------------------------------
if(WIN32)
  # We cannot use the compiledata.sh script for windows
  configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/compiledata.win32.in include/compiledata.h)
else()
  execute_process(COMMAND ${CMAKE_SOURCE_DIR}/build/unix/compiledata.sh include/compiledata.h "${CXX}" ""
       "${CMAKE_CXX_FLAGS_${uppercase_CMAKE_BUILD_TYPE}}"
	     "${CMAKE_CXX_FLAGS}" "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS}" "${CMAKE_EXE_FLAGS}" "${LibSuffix}" "${SYSLIBS}"
	     "${libdir}" "-lCore" "-lRint" "${incdir}" "" "" "${ROOT_ARCHITECTURE}" "" "${explicitlink}" )
endif()

configure_file(${CMAKE_SOURCE_DIR}/config/root-config.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/root-config @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/memprobe.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/memprobe @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/thisroot.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/thisroot.sh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/thisroot.csh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/thisroot.csh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/setxrd.csh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/setxrd.csh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/setxrd.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/setxrd.sh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/proofserv.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/proofserv @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/roots.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/roots @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/config/root-help.el.in root-help.el @ONLY)

if(WIN32)
  set(thisrootbat ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/thisroot.bat)
  configure_file(${CMAKE_SOURCE_DIR}/config/thisroot.bat ${thisrootbat} @ONLY)
endif()

install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/memprobe
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/thisroot.sh
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/thisroot.csh
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/setxrd.csh
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/setxrd.sh
              ${thisrootbat}
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/root-config
              ${CMAKE_SOURCE_DIR}/cmake/scripts/setenvwrap.csh
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/roots
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/proofserv
              PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
                          GROUP_EXECUTE GROUP_READ 
                          WORLD_EXECUTE WORLD_READ 
              DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${CMAKE_BINARY_DIR}/include/RConfigOptions.h
              ${CMAKE_BINARY_DIR}/include/compiledata.h 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/etc/root.mimes 
              ${CMAKE_BINARY_DIR}/etc/system.rootrc
              DESTINATION ${CMAKE_INSTALL_SYSCONFDIR})
              
install(FILES ${CMAKE_BINARY_DIR}/root-help.el DESTINATION ${CMAKE_INSTALL_ELISPDIR})


#install(FILES ${CMAKE_BINARY_DIR}/config/Makefile.comp
#              ${CMAKE_BINARY_DIR}/config/Makefile.config
#              DESTINATION config)


endfunction()
RootConfigure()
