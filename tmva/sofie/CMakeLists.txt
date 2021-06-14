# Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building TMVA SOFIE package
############################################################################
#Author: Sitong An, Lorenzo Moneta 10/03/2021

#sofie is built only if protobuf is found

protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS "onnx_proto3")
set_source_files_properties(${PROTO_SRCS} ${PROTO_HDRS} PROPERTIES GENERATED TRUE)
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS_Development_Main} ${NUMPY_INCLUDE_DIRS})

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTTMVASofie
  HEADERS
   TMVA/OperatorList.hxx
   TMVA/RModel.hxx
   TMVA/RModelParser_ONNX.hxx
   TMVA/RModelParser_Keras.hxx
   TMVA/ROperator.hxx
   TMVA/ROperator_Gemm.hxx
   TMVA/ROperator_Relu.hxx
   TMVA/ROperator_Transpose.hxx
   TMVA/SOFIE_common.hxx
   ${PROTO_HDRS}
  SOURCES
    src/RModel.cxx
    src/RModelParser_ONNX.cxx
    src/RModelParser_Keras.cxx
    src/SOFIE_common.cxx
    ${PROTO_SRCS}
  LIBRARIES
    ${Protobuf_LIBRARIES}
    ${PYTHON_LIBRARIES_Development_Main}
  DEPENDENCIES
    TMVA
)

target_include_directories(ROOTTMVASofie PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_include_directories(ROOTTMVASofie PUBLIC
  ${Protobuf_INCLUDE_DIRS})
set_target_properties(ROOTTMVASofie PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE)
