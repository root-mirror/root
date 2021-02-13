# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/meta package
############################################################################

set(Meta_dict_headers
  TBaseClass.h
  TClassGenerator.h
  TClass.h
  TClassRef.h
  TClassStreamer.h
  TDataMember.h
  TDataType.h
  TDictAttributeMap.h
  TDictionary.h
  TEnumConstant.h
  TEnum.h
  TFileMergeInfo.h
  TFunction.h
  TFunctionTemplate.h
  TGenericClassInfo.h
  TGlobal.h
  TInterpreter.h
  TInterpreterValue.h
  TIsAProxy.h
  TListOfDataMembers.h
  TListOfEnums.h
  TListOfEnumsWithLock.h
  TListOfFunctions.h
  TListOfFunctionTemplates.h
  TMemberStreamer.h
  TMethodArg.h
  TMethodCall.h
  TMethod.h
  TProtoClass.h
  TRealData.h
  TSchemaHelper.h
  TSchemaRule.h
  TSchemaRuleSet.h
  TStatusBitsChecker.h
  TStreamerElement.h
  TStreamer.h
  TVirtualIsAProxy.h
  TVirtualRefProxy.h
  TVirtualStreamerInfo.h
  TVirtualArray.h
  TVirtualObject.h
  PARENT_SCOPE
)

ROOT_OBJECT_LIBRARY(Meta
  src/TBaseClass.cxx
  src/TClass.cxx
  src/TClassGenerator.cxx
  src/TClassRef.cxx
  src/TDataMember.cxx
  src/TDataType.cxx
  src/TDictAttributeMap.cxx
  src/TDictionary.cxx
  src/TEnum.cxx
  src/TEnumConstant.cxx
  src/TFunction.cxx
  src/TFunctionTemplate.cxx
  src/TGenericClassInfo.cxx
  src/TGlobal.cxx
  src/TInterpreter.cxx
  src/TIsAProxy.cxx
  src/TListOfDataMembers.cxx
  src/TListOfEnums.cxx
  src/TListOfEnumsWithLock.cxx
  src/TListOfFunctions.cxx
  src/TListOfFunctionTemplates.cxx
  src/TMethod.cxx
  src/TMethodArg.cxx
  src/TMethodCall.cxx
  src/TProtoClass.cxx
  src/TRealData.cxx
  src/TSchemaRule.cxx
  src/TSchemaRuleSet.cxx
  src/TStatusBitsChecker.cxx
  src/TStreamerElement.cxx
  src/TViewPubDataMembers.cxx
  src/TViewPubFunctions.cxx
  src/TVirtualStreamerInfo.cxx
)

target_include_directories(Meta PRIVATE
  ${CMAKE_SOURCE_DIR}/core/meta/res
  ${CMAKE_SOURCE_DIR}/core/foundation/res
  ${CMAKE_SOURCE_DIR}/core/base/inc
  ${CMAKE_SOURCE_DIR}/core/clib/inc
  ${CMAKE_SOURCE_DIR}/core/cont/inc
  ${CMAKE_SOURCE_DIR}/core/thread/inc
  ${CMAKE_SOURCE_DIR}/core/gui/inc
  ${CMAKE_SOURCE_DIR}/core/foundation/inc
  ${CMAKE_BINARY_DIR}/ginclude
)

if(MSVC)
   target_include_directories(Meta PRIVATE
      ${CMAKE_SOURCE_DIR}/core/winnt/inc
   )
endif()


ROOT_INSTALL_HEADERS()
ROOT_ADD_TEST_SUBDIRECTORY(test)
