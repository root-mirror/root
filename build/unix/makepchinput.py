#! /usr/bin/env python
#
# Extract the input needed to build a PCH for the main (enabled) ROOT modules.
# Script takes as argument the source directory path, the set of enabled
# modules and extra headers (from cling) to be included in the PCH.
#
# Copyright (c) 2014 Rene Brun and Fons Rademakers
# Author: Axel Naumann <axel@cern.ch>, 2014-10-16
# Translated to Python by Danilo Piparo, 2015-04-22

from __future__ import print_function
import sys
import os
import glob
import shutil

#-------------------------------------------------------------------------------
def removeFiles(filesList):
   existingFilesList = filter(os.path.exists,filesList)
   map(os.unlink, existingFilesList)

#-------------------------------------------------------------------------------
def removeLeftOvers(filesToRemove):
   """
   Remove leftover files from old versions of this script.
   """
   filesToRemove.extend(
                    [os.path.join("include","allHeaders.h.pch"),
                    os.path.join("etc","allDict.cxx"),
                    os.path.join("etc","allDict.cxx.h"),
                    os.path.join("etc","allDict.cxx.pch")]
                    )
   removeFiles(filesToRemove)

#-------------------------------------------------------------------------------
def getParams():
   """
   Extract parameters from the commandline, which looks like
   makePCHInput.py WWW XXX YYY ZZZ -- CXXFLAGS
   """
   argv = sys.argv
   rootSrcDir, modules, legacyPyROOT = argv[1:4]
   clingetpchList = argv[4:]

   return rootSrcDir, modules, legacyPyROOT == 'ON', clingetpchList

#-------------------------------------------------------------------------------
def getGuardedStlInclude(headerName):
   return '#if __has_include("%s")\n' %headerName +\
          '#include <%s>\n' %headerName +\
          '#endif\n'

#-------------------------------------------------------------------------------
def getSTLIncludes():
   """
   Here we include the list of c++11 stl headers
   From http://en.cppreference.com/w/cpp/header
   valarray is removed because it causes lots of compilation at startup.
   """
   stlHeadersList = ("cstdlib",
                     "csignal",
                     "csetjmp",
                     "cstdarg",
                     "typeinfo",
                     "typeindex",
                     "type_traits",
                     "bitset",
                     "functional",
                     "utility",
                     "ctime",
                     "chrono",
                     "cstddef",
                     "initializer_list",
                     "tuple",
                     "new",
                     "memory",
                     "scoped_allocator",
                     "climits",
                     "cfloat",
                     "cstdint",
                     "cinttypes",
                     "limits",
                     "exception",
                     "stdexcept",
                     "cassert",
                     "system_error",
                     "cerrno",
                     "cctype",
                     "cwctype",
                     "cstring",
                     "cwchar",
                     "cuchar",
                     "string",
                     "array",
                     "vector",
                     "deque",
                     "list",
                     "forward_list",
                     "set",
                     "map",
                     "unordered_set",
                     "unordered_map",
                     "stack",
                     "queue",
                     "algorithm",
                     "iterator",
                     "cmath",
                     "complex",
#                     "valarray",
                     "random",
                     "numeric",
                     "ratio",
                     "cfenv",
                     "iosfwd",
                     "ios",
                     "istream",
                     "ostream",
                     "iostream",
                     "fstream",
                     "sstream",
                     "iomanip",
                     "streambuf",
                     "cstdio",
                     "locale",
                     "clocale",
                     "codecvt",
                     "atomic",
                     "thread",
                     "mutex",
                     "future",
                     "condition_variable",
                     "ciso646",
                     "ccomplex",
                     "ctgmath",
                     "regex",
                     "cstdbool")

   allHeadersPartContent = "// STL headers\n"

   for header in stlHeadersList:
      allHeadersPartContent += getGuardedStlInclude(header)

   # Special case for regex
   allHeadersPartContent += '// treat regex separately\n' +\
                            '#if __has_include("regex") && !defined __APPLE__\n' +\
                            '#include <regex>\n' +\
                            '#endif\n'

   # treat this deprecated headers in a special way
   stlDeprecatedHeadersList=["strstream"]
   allHeadersPartContent += '// STL Deprecated headers\n' +\
                            '#define _BACKWARD_BACKWARD_WARNING_H\n' +\
                            "#pragma clang diagnostic push\n" +\
                            '#pragma GCC diagnostic ignored "-Wdeprecated"\n'

   for header in stlDeprecatedHeadersList:
      allHeadersPartContent += getGuardedStlInclude(header)

   allHeadersPartContent += '#pragma clang diagnostic pop\n' +\
                            '#undef _BACKWARD_BACKWARD_WARNING_H\n'
   return allHeadersPartContent

#-------------------------------------------------------------------------------
def getExtraIncludes(headers):
   """
   Add include files according to list
   """
   allHeadersPartContent=""
   for header in headers:
      allHeadersPartContent+='#include "%s"\n' %header
   return allHeadersPartContent

#-------------------------------------------------------------------------------
def getDictNames(theDirName):
   """
   Get a list of all dictionaries in a directory
   """
   #`find $modules -name 'G__*.cxx' 2> /dev/null | grep -v core/metautils/src/G__std_`; do
   wildcards = (os.path.join(theDirName , "*", "*", "G__*.cxx"),
                os.path.join(theDirName , "*", "G__*.cxx"))
   allDictNames = []
   for wildcard in wildcards:
      allDictNames += glob.glob(wildcard)
   stdDictpattern = os.path.join("core","metautils","src","G__std_")
   dictNames = filter (lambda dictName: not (stdDictpattern in dictName or "/roottest/" in dictName),allDictNames )
   return dictNames

#-------------------------------------------------------------------------------
def getDirName(dictName):
   """
   foo/src/G__PIPPO.cxx --> foo/
   """
   # get rid of the drive on win
   dirName = os.path.splitdrive(dictName)

   # foo/src/G__PIPPO.cxx --> foo/src/
   dirName = os.path.split(dictName)[0]

   # foo/src/ --> foo/src
   dirName = dirName[:-1]

   # foo/src --> foo/
   return os.path.normpath(os.path.split(dictName)[0])

#-------------------------------------------------------------------------------
def isAnyPatternInString(patterns,theString):
   """
   Check if any of the patterns is contained in the string
   """
   for pattern in patterns:
      if os.path.normpath(pattern) in theString: return True
   return False

#-------------------------------------------------------------------------------
def isDirForPCH(dirName, legacyPyROOT):
   """
   Check if the directory corresponds to a module whose headers must belong to
   the PCH
   """
   PCHPatternsWhitelist = ["interpreter/",
                           "core/",
                           "io/io",
                           "net/net",
                           "math/",
                           "hist/",
                           "tree/",
                           "graf2d",
                           "graf3d/ftgl",
                           "graf3d/g3d",
                           "graf3d/gl",
                           "gui/gui",
                           "gui/fitpanel",
                           "rootx",
                           "roofit/",
                           "tmva",
                           "main"]
   if legacyPyROOT:
      PCHPatternsWhitelist.append("bindings/pyroot_legacy")
   else:
      PCHPatternsWhitelist.append("bindings/tpython")

   PCHPatternsBlacklist = ["gui/guihtml",
                           "gui/guibuilder",
                           "math/fftw",
                           "math/foam",
                           "math/fumili",
                           "math/mlp",
                           "math/quadp",
                           "math/rtools",
                           "math/splot",
                           "math/unuran",
                           "math/vdt",
                           "tmva/rmva"]

   if (sys.platform != 'win32' and sys.maxsize <= 2**32): # https://docs.python.org/3/library/platform.html#cross-platform
      PCHPatternsBlacklist.append("tree/dataframe")

   accepted = isAnyPatternInString(PCHPatternsWhitelist,dirName) and \
               not isAnyPatternInString(PCHPatternsBlacklist,dirName)

   return accepted

#-------------------------------------------------------------------------------
def getLinesFromDict(marker, dictFileName):
   """
   Search for the line marker
   in the dictionary and save all lines until the line '0'
   Return them as List
   """
   selectedLines = []
   ifile = open(dictFileName)
   lines = ifile.readlines()
   ifile.close()
   recording = False
   for line in lines:
      if marker in line:
         recording = True
         continue

      if recording and "0" == line[0]: break

      if recording:
         selectedLines.append(line[:-1])

   return selectedLines

#-------------------------------------------------------------------------------
def getIncludeLinesFromDict(dictFileName):
   """
   Search for the headers after the line
   'static const char* headers[]'
   Return the code to be added to the allHeaders as string
   """
   allHeadersPartContent=""
   selectedLines = getLinesFromDict('static const char* headers[]', dictFileName)
   allHeadersPartContent += "// %s\n"  % dictFileName
   for selectedLine in selectedLines:
      header = selectedLine[:-1] # remove the ","
      allHeadersPartContent += "#include %s\n" %header
   return allHeadersPartContent

#-------------------------------------------------------------------------------
def getIncludePathsFromDict(dictFileName):
   """
   Search for the include paths after the line
   'static const char* includePaths[]'
   Return them as list
   """
   incPathsPart=[]
   selectedLines = getLinesFromDict('static const char* includePaths[]', dictFileName)
   for selectedLine in selectedLines:
      incPath = selectedLine[1:-2] # remove the "," and the two '"'
      incPathsPart.append(incPath)
   return incPathsPart

#-------------------------------------------------------------------------------
def getDefUndefLines(dirName):
   """
   Add undefines and defines if the directory needs them
   """
   allHeadersPartContent=""
   if "%sqt" %os.sep in dirName:
      allHeadersPartContent += '#ifdef emit\n' +\
                               '# undef emit\n' +\
                               '#endif\n' +\
                               '#ifdef signals\n' +\
                               '# undef signals\n' +\
                               '#endif\n'
   return allHeadersPartContent

#-------------------------------------------------------------------------------
def mkdirIfNotThere(dirName):
   if not os.path.exists(dirName):
      os.makedirs(dirName)

#-------------------------------------------------------------------------------
def copyLinkDefs(rootSrcDir, outdir):
   """
   Extract the linkdef files
   """
   linkDefPartContent = ""
   curDir = os.getcwd()
   os.chdir(rootSrcDir)
   wildcards = (os.path.join("*", "inc", "*LinkDef*.h"),
                os.path.join("*", "*", "inc", "*LinkDef*.h"),
                os.path.join("*", "*", "inc", "*" , "*LinkDef*.h"))
   linkDefNames = []
   for wildcard in wildcards:
      linkDefNames += glob.glob(wildcard)
   os.chdir(curDir)
   for linkDefName in linkDefNames:
      linkDefDirName = os.path.dirname(linkDefName)
      mkdirIfNotThere(os.path.join(outdir,linkDefDirName))
      srcName = os.path.join(rootSrcDir,linkDefName)
      destName = os.path.join(outdir,linkDefName)
      shutil.copyfile(srcName, destName)

#-------------------------------------------------------------------------------
def getLocalLinkDefs(rootSrcDir, outdir , dirName):
   linkDefPartContent = ""
   curDir = os.getcwd()
   os.chdir(rootSrcDir)
   wildcards = (os.path.join(dirName , "*", "*", "*LinkDef*.h"),
                os.path.join(dirName , "*", "*LinkDef*.h"))
   linkDefNames = []
   for wildcard in wildcards:
      linkDefNames += glob.glob(wildcard)

   # now get the ones in the inc directory
   linkDefNames = filter (lambda name: "%sinc%s" %(os.sep,os.sep) in name, linkDefNames)

   for linkDefName in linkDefNames:
      fullLinkDefName = os.path.join(outdir,linkDefName)
      linkDefPartContent += '#include "%s"\n' %fullLinkDefName
   os.chdir(curDir)
   return linkDefPartContent

#-------------------------------------------------------------------------------
def resolveSoftLinks(thePaths):
   return map(os.path.realpath,thePaths)

#-------------------------------------------------------------------------------
def getCppFlags(rootSrcDir,allIncPaths):
   """
   Sort, clean, no duplicates
   cat $cppflags.tmp | sort | uniq | grep -v $srcdir | grep -v `pwd` > $cppflags
   We must resolve softlinks.
   returns a string
   """
   allHeadersPartContent = ""
   filteredIncPaths = sorted(list(set(resolveSoftLinks(allIncPaths))))
   for name in resolveSoftLinks((rootSrcDir,os.getcwd())):
      filteredIncPaths = filter (lambda incPath: not name in incPath,filteredIncPaths)
   for incPath in filteredIncPaths:
      allHeadersPartContent += "-I%s\n" %incPath
   return allHeadersPartContent

#-------------------------------------------------------------------------------
def writeToFile(content, filename):
   ofile = open(filename, "w")
   ofile.write(content)
   ofile.close()

#-------------------------------------------------------------------------------
def writeFiles(contentFileNamePairs):
   for content, filename in contentFileNamePairs:
      writeToFile(content, filename)

#-------------------------------------------------------------------------------
def printModulesMessageOnScreen(selModules):
   modulesList = sorted(list(selModules))
   print ("\nGenerating PCH for %s\n" %" ".join(modulesList))

#-------------------------------------------------------------------------------
def getExtraHeaders():
   """ Get extra headers which do not fall in other special categories
   """
   extraHeaders=["ROOT/TSeq.hxx","ROOT/StringConv.hxx"]
   code = "// Extra headers\n"
   for extraHeader in extraHeaders:
      code += '#include "%s"\n' %extraHeader
   return code

#-------------------------------------------------------------------------------
def removeUnwantedHeaders(allHeadersContent):
   """ remove unwanted headers, e.g. the ones used for dictionaries but not desirable in the pch
   """
   unwantedHeaders = []
   deprecatedHeaders = ['']
   unwantedHeaders.extend(deprecatedHeaders)
   for unwantedHeader in unwantedHeaders:
      allHeadersContent = allHeadersContent.replace('#include "%s"' %unwantedHeader,"")
   return allHeadersContent


#-------------------------------------------------------------------------------
def makePCHInput():
   """
   Create the input for the pch file, i.e. 3 files:
      * etc/dictpch/allLinkDefs.h
      * etc/dictpch/allHeaders.h
      * etc/dictpch/allCppflags.txt
   """
   rootSrcDir, modules, legacyPyROOT, clingetpchList = getParams()

   outdir = os.path.join("etc","dictpch")
   allHeadersFilename = os.path.join(outdir,"allHeaders.h")
   allLinkdefsFilename = os.path.join(outdir,"allLinkDefs.h")
   cppFlagsFilename = os.path.join(outdir, "allCppflags.txt")

   if sys.platform == 'win32':
      outdir.replace("\\","/")
      allHeadersFilename.replace("\\","/")
      allLinkdefsFilename.replace("\\","/")
      cppFlagsFilename.replace("\\","/")

   mkdirIfNotThere(outdir)
   removeLeftOvers([allHeadersFilename, allLinkdefsFilename, cppFlagsFilename])

   allHeadersContent = getSTLIncludes()
   allHeadersContent += getExtraIncludes(clingetpchList)

   allLinkdefsContent = ""

   # Loop over the dictionaries, ROOT modules
   dictNames = getDictNames(modules)
   selModules = set([])
   allIncPathsList = []
   for dictName in dictNames:
      dirName = getDirName(dictName)
      if not isDirForPCH(dirName, legacyPyROOT): continue

      selModules.add(dirName)

      allHeadersContent += getIncludeLinesFromDict(dictName)
      allIncPathsList += getIncludePathsFromDict(dictName)

      allHeadersContent += getDefUndefLines(dictName)

      allLinkdefsContent += getLocalLinkDefs(rootSrcDir, outdir , dirName)

   allHeadersContent += getExtraHeaders()

   allHeadersContent = removeUnwantedHeaders(allHeadersContent)

   copyLinkDefs(rootSrcDir, outdir)

   cppFlagsContent = getCppFlags(rootSrcDir, allIncPathsList) + '\n'

   writeFiles(((allHeadersContent, allHeadersFilename),
               (allLinkdefsContent, allLinkdefsFilename),
               (cppFlagsContent, cppFlagsFilename)))


   printModulesMessageOnScreen(selModules)

if __name__ == "__main__":
   makePCHInput()
