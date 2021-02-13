#ifndef ROOT7_RNTuple_Test
#define ROOT7_RNTuple_Test

#include <ROOT/RColumnModel.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RFieldVisitor.hxx>
#include <ROOT/RMiniFile.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <ROOT/RNTupleMerger.hxx>
#include <ROOT/RNTupleMetrics.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RNTupleZip.hxx>
#include <ROOT/RPageAllocator.hxx>
#include <ROOT/RPagePool.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RPageStorageFile.hxx>
#include <ROOT/RRawFile.hxx>
#include <ROOT/RVec.hxx>

#include <RZip.h>
#include <TClass.h>
#include <TFile.h>
#include <TRandom3.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "CustomStruct.hxx"

#include <array>
#include <chrono>
#include <cstdio>
#include <exception>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#if __cplusplus >= 201703L
#include <variant>
#endif
#include <vector>

using DescriptorId_t = ROOT::Experimental::DescriptorId_t;
using EColumnType = ROOT::Experimental::EColumnType;
using ENTupleContainerFormat = ROOT::Experimental::ENTupleContainerFormat;
using ENTupleStructure = ROOT::Experimental::ENTupleStructure;
using NTupleSize_t = ROOT::Experimental::NTupleSize_t;
using RColumnModel = ROOT::Experimental::RColumnModel;
using RDanglingFieldDescriptor = ROOT::Experimental::RDanglingFieldDescriptor;
using RException = ROOT::Experimental::RException;
template <class T>
using RField = ROOT::Experimental::RField<T>;
using RFieldBase = ROOT::Experimental::Detail::RFieldBase;
using RFieldDescriptor = ROOT::Experimental::RFieldDescriptor;
using RFieldFuse = ROOT::Experimental::Detail::RFieldFuse;
using RFieldMerger = ROOT::Experimental::RFieldMerger;
using RFieldValue = ROOT::Experimental::Detail::RFieldValue;
using RMiniFileReader = ROOT::Experimental::Internal::RMiniFileReader;
using RNTuple = ROOT::Experimental::RNTuple;
using RNTupleAtomicCounter = ROOT::Experimental::Detail::RNTupleAtomicCounter;
using RNTupleAtomicTimer = ROOT::Experimental::Detail::RNTupleAtomicTimer;
using RNTupleCalcPerf = ROOT::Experimental::Detail::RNTupleCalcPerf;
using RNTupleCompressor = ROOT::Experimental::Detail::RNTupleCompressor;
using RNTupleDecompressor = ROOT::Experimental::Detail::RNTupleDecompressor;
using RNTupleDescriptor = ROOT::Experimental::RNTupleDescriptor;
using RNTupleDescriptorBuilder = ROOT::Experimental::RNTupleDescriptorBuilder;
using RNTupleFileWriter = ROOT::Experimental::Internal::RNTupleFileWriter;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleReadOptions = ROOT::Experimental::RNTupleReadOptions;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleWriteOptions = ROOT::Experimental::RNTupleWriteOptions;
using RNTupleMetrics = ROOT::Experimental::Detail::RNTupleMetrics;
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTuplePlainCounter = ROOT::Experimental::Detail::RNTuplePlainCounter;
using RNTuplePlainTimer = ROOT::Experimental::Detail::RNTuplePlainTimer;
using RNTupleVersion = ROOT::Experimental::RNTupleVersion;
using RPage = ROOT::Experimental::Detail::RPage;
using RPageAllocatorHeap = ROOT::Experimental::Detail::RPageAllocatorHeap;
using RPageDeleter = ROOT::Experimental::Detail::RPageDeleter;
using RPagePool = ROOT::Experimental::Detail::RPagePool;
using RPageSink = ROOT::Experimental::Detail::RPageSink;
using RPageSinkFile = ROOT::Experimental::Detail::RPageSinkFile;
using RPageSource = ROOT::Experimental::Detail::RPageSource;
using RPageSourceFile = ROOT::Experimental::Detail::RPageSourceFile;
using RPrepareVisitor = ROOT::Experimental::RPrepareVisitor;
using RPrintSchemaVisitor = ROOT::Experimental::RPrintSchemaVisitor;
using RRawFile = ROOT::Internal::RRawFile;
template <class T>
using RResult = ROOT::Experimental::RResult<T>;

/**
 * An RAII wrapper around an open temporary file on disk. It cleans up the guarded file when the wrapper object
 * goes out of scope.
 */
class FileRaii {
private:
   std::string fPath;
public:
   explicit FileRaii(const std::string &path) : fPath(path) { }
   FileRaii(const FileRaii&) = delete;
   FileRaii& operator=(const FileRaii&) = delete;
   ~FileRaii() { std::remove(fPath.c_str()); }
   std::string GetPath() const { return fPath; }
};

#endif
