// Author: Jakob Blomer CERN  07/2018

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// clang-format off
/** \class ROOT::RDF::RSqliteDS
    \ingroup dataframe
    \brief RDataFrame data source class for reading SQlite files.
*/

// clang-format on

#include <ROOT/RSqliteDS.hxx>
#include <ROOT/RConfig.hxx>
#include <ROOT/RDF/Utils.hxx>
#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RRawFile.hxx>

#include "TError.h"
#include "TRandom.h"
#include "TSystem.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring> // for memcpy
#include <ctime>
#include <memory> // for placement new
#include <stdexcept>
#include <utility>

#include <sqlite3.h>

namespace {

// In order to provide direct access to remote sqlite files through HTTP and HTTPS, this datasource provides a custom
// "SQlite VFS module" that uses Davix for data access. The SQlite VFS modules are roughly what TSystem is
// for ROOT -- an abstraction of the operating system interface.
//
// SQlite allows for registering custom VFS modules, which are a set of C callback functions that SQlite invokes when
// it needs to read from a file, write to a file, etc. More information is available under https://sqlite.org/vfs.html
//
// In the context of a data source, SQlite will only ever call reading functions from the VFS module, the sqlite
// files are not modified. Therefore, only a subset of the callback functions provide a non-trivial implementation.
// The custom VFS module uses a RRawFile for the byte access, thereby it can access local and remote files.

////////////////////////////////////////////////////////////////////////////
/// SQlite VFS modules are identified by string names. The name has to be unique for the entire application.
constexpr char const *gSQliteVfsName = "ROOT-Davix-readonly";

////////////////////////////////////////////////////////////////////////////
/// Holds the state of an open sqlite database. Objects of this struct are created in VfsRdOnlyOpen()
/// and then passed by sqlite to the file I/O callbacks (read, close, etc.). This uses C style inheritance
/// where the struct starts with a sqlite3_file member (base class) which is extended by members related to
/// this particular VFS module. Every callback here thus casts the sqlite3_file input parameter to its "derived"
/// type VfsRootFile.
struct VfsRootFile {
   VfsRootFile() = default;

   sqlite3_file pFile;
   std::unique_ptr<ROOT::Internal::RRawFile> fRawFile;
};

// The following callbacks implement the I/O operations of an open database

////////////////////////////////////////////////////////////////////////////
/// Releases the resources associated to a file opened with davix
int VfsRdOnlyClose(sqlite3_file *pFile)
{
   VfsRootFile *p = reinterpret_cast<VfsRootFile *>(pFile);
   // We can't use delete because the storage for p is managed by sqlite
   p->~VfsRootFile();
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// Issues a byte range request for a chunk to the raw file
int VfsRdOnlyRead(sqlite3_file *pFile, void *zBuf, int count, sqlite_int64 offset)
{
   VfsRootFile *p = reinterpret_cast<VfsRootFile *>(pFile);
   auto nbytes = p->fRawFile->ReadAt(zBuf, count, offset);
   return (nbytes != static_cast<unsigned int>(count)) ? SQLITE_IOERR : SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// We do not write to a database in the RDataSource and therefore can simply return an error for this callback
int VfsRdOnlyWrite(sqlite3_file * /*pFile*/, const void * /*zBuf*/, int /*iAmt*/, sqlite_int64 /*iOfst*/)
{
   return SQLITE_OPEN_READONLY;
}

////////////////////////////////////////////////////////////////////////////
/// We do not write to a database in the RDataSource and therefore can simply return an error for this callback
int VfsRdOnlyTruncate(sqlite3_file * /*pFile*/, sqlite_int64 /*size*/)
{
   return SQLITE_OPEN_READONLY;
}

////////////////////////////////////////////////////////////////////////////
/// As the database is read-only, syncing data to disc is a no-op and always succeeds
int VfsRdOnlySync(sqlite3_file * /*pFile*/, int /*flags*/)
{
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// Returns the cached file size
int VfsRdOnlyFileSize(sqlite3_file *pFile, sqlite_int64 *pSize)
{
   VfsRootFile *p = reinterpret_cast<VfsRootFile *>(pFile);
   *pSize = p->fRawFile->GetSize();
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// As the database is read-only, locks for concurrent access are no-ops and always succeeds
int VfsRdOnlyLock(sqlite3_file * /*pFile*/, int /*level*/)
{
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// As the database is read-only, locks for concurrent access are no-ops and always succeeds
int VfsRdOnlyUnlock(sqlite3_file * /*pFile*/, int /*level*/)
{
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// As the database is read-only, locks for concurrent access are no-ops and always succeeds
int VfsRdOnlyCheckReservedLock(sqlite3_file * /*pFile*/, int *pResOut)
{
   *pResOut = 0;
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// As the database is read-only, we know there are no additional control files such as a database journal
int VfsRdOnlyFileControl(sqlite3_file * /*p*/, int /*op*/, void * /*pArg*/)
{
   return SQLITE_NOTFOUND;
}

////////////////////////////////////////////////////////////////////////////
/// The database device's sector size is only needed for writing
int VfsRdOnlySectorSize(sqlite3_file * /*pFile*/)
{
   return SQLITE_OPEN_READONLY;
}

////////////////////////////////////////////////////////////////////////////
/// The database device's properties are only needed for writing
int VfsRdOnlyDeviceCharacteristics(sqlite3_file * /*pFile*/)
{
   return SQLITE_OPEN_READONLY;
}

////////////////////////////////////////////////////////////////////////////
/// Set the function pointers of the custom VFS I/O operations in a
/// forward-compatible way
static sqlite3_io_methods GetSqlite3IoMethods()
{
   // The C style initialization is compatible with version 1 and later versions of the struct.
   // Version 1 was introduced with sqlite 3.6, version 2 with sqlite 3.7.8, version 3 with sqlite 3.7.17
   sqlite3_io_methods io_methods;
   memset(&io_methods, 0, sizeof(io_methods));
   io_methods.iVersion               = 1;
   io_methods.xClose                 = VfsRdOnlyClose;
   io_methods.xRead                  = VfsRdOnlyRead;
   io_methods.xWrite                 = VfsRdOnlyWrite;
   io_methods.xTruncate              = VfsRdOnlyTruncate;
   io_methods.xSync                  = VfsRdOnlySync;
   io_methods.xFileSize              = VfsRdOnlyFileSize;
   io_methods.xLock                  = VfsRdOnlyLock;
   io_methods.xUnlock                = VfsRdOnlyUnlock;
   io_methods.xCheckReservedLock     = VfsRdOnlyCheckReservedLock;
   io_methods.xFileControl           = VfsRdOnlyFileControl;
   io_methods.xSectorSize            = VfsRdOnlySectorSize;
   io_methods.xDeviceCharacteristics = VfsRdOnlyDeviceCharacteristics;
   return io_methods;
}

////////////////////////////////////////////////////////////////////////////
/// Fills a new VfsRootFile struct enclosing a Davix file
int VfsRdOnlyOpen(sqlite3_vfs * /*vfs*/, const char *zName, sqlite3_file *pFile, int flags, int * /*pOutFlags*/)
{
   // Storage for the VfsRootFile structure has been already allocated by sqlite, so we use placement new
   VfsRootFile *p = new (pFile) VfsRootFile();
   p->pFile.pMethods = nullptr;

   // This global struct contains the function pointers to all the callback operations that act on an open database.
   // It is passed via the pFile struct back to sqlite so that it can call back to the functions provided above.
   static const sqlite3_io_methods io_methods = GetSqlite3IoMethods();

   if (flags & (SQLITE_OPEN_READWRITE | SQLITE_OPEN_DELETEONCLOSE | SQLITE_OPEN_EXCLUSIVE))
      return SQLITE_IOERR;

   p->fRawFile = ROOT::Internal::RRawFile::Create(zName);
   if (!p->fRawFile) {
      ::Error("VfsRdOnlyOpen", "Cannot open %s\n", zName);
      return SQLITE_IOERR;
   }

   if (!(p->fRawFile->GetFeatures() & ROOT::Internal::RRawFile::kFeatureHasSize)) {
      ::Error("VfsRdOnlyOpen", "cannot determine file size of %s\n", zName);
      return SQLITE_IOERR;
   }

   p->pFile.pMethods = &io_methods;
   return SQLITE_OK;
}

// The following callbacks implement operating system specific functionality. In contrast to the previous callbacks,
// there is no need to implement any customized logic for the following ones. An implementation has to be
// provided nevertheless to have a fully functional VFS module.

////////////////////////////////////////////////////////////////////////////
/// This VFS module cannot remove files
int VfsRdOnlyDelete(sqlite3_vfs * /*vfs*/, const char * /*zName*/, int /*syncDir*/)
{
   return SQLITE_IOERR_DELETE;
}

////////////////////////////////////////////////////////////////////////////
/// Access control always allows read-only access to databases
int VfsRdOnlyAccess(sqlite3_vfs * /*vfs*/, const char * /*zPath*/, int flags, int *pResOut)
{
   *pResOut = 0;
   if (flags == SQLITE_ACCESS_READWRITE) {
      return SQLITE_OPEN_READONLY;
   }
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// No distinction between relative and full paths for URLs, returns the input path name
int VfsRdOnlyFullPathname(sqlite3_vfs * /*vfs*/, const char *zPath, int nOut, char *zOut)
{
   zOut[nOut - 1] = '\0';
   sqlite3_snprintf(nOut, zOut, "%s", zPath);
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// Let TRandom fill the buffer with random bytes
int VfsRdOnlyRandomness(sqlite3_vfs * /*vfs*/, int nBuf, char *zBuf)
{
   for (int i = 0; i < nBuf; ++i) {
      zBuf[i] = (char)gRandom->Integer(256);
   }
   return nBuf;
}

////////////////////////////////////////////////////////////////////////////
/// Use ROOT's platform independent sleep wrapper
int VfsRdOnlySleep(sqlite3_vfs * /*vfs*/, int microseconds)
{
   // Millisecond precision but sleep at least number of given microseconds as requested
   gSystem->Sleep((microseconds + 1000 - 1) / 1000);
   return microseconds;
}

////////////////////////////////////////////////////////////////////////////
/// Use sqlite default implementation
int VfsRdOnlyGetLastError(sqlite3_vfs * /*vfs*/, int /*not_used1*/, char * /*not_used2*/)
{
   return errno;
}

////////////////////////////////////////////////////////////////////////////
/// Return UTC as being done in the sqlite unix VFS without gettimeofday()
int VfsRdOnlyCurrentTimeInt64(sqlite3_vfs * /*vfs*/, sqlite3_int64 *piNow)
{
   static constexpr sqlite3_int64 unixEpoch = 24405875 * (sqlite3_int64)8640000;
   time_t t;
   time(&t);
   *piNow = ((sqlite3_int64)t) * 1000 + unixEpoch;
   return SQLITE_OK;
}

////////////////////////////////////////////////////////////////////////////
/// Wrapper around VfsRdOnlyCurrentTimeInt64
int VfsRdOnlyCurrentTime(sqlite3_vfs *vfs, double *prNow)
{
   sqlite3_int64 i = 0;
   int rc = VfsRdOnlyCurrentTimeInt64(vfs, &i);
   *prNow = i / 86400000.0;
   return rc;
}

////////////////////////////////////////////////////////////////////////////
/// Set the function pointers of the VFS implementation in a
/// forward-compatible way
static sqlite3_vfs GetSqlite3Vfs()
{
   // The C style initialization is compatible with version 1 and later versions of the struct.
   // Version 1 was introduced with sqlite 3.5, version 2 with sqlite 3.7, version 3 with sqlite 3.7.6
   sqlite3_vfs vfs;
   memset(&vfs, 0, sizeof(vfs));
   vfs.iVersion      = 1;
   vfs.szOsFile      = sizeof(VfsRootFile);
   vfs.mxPathname    = 2000;
   vfs.zName         = gSQliteVfsName;
   vfs.xOpen         = VfsRdOnlyOpen;
   vfs.xDelete       = VfsRdOnlyDelete;
   vfs.xAccess       = VfsRdOnlyAccess;
   vfs.xFullPathname = VfsRdOnlyFullPathname;
   vfs.xRandomness   = VfsRdOnlyRandomness;
   vfs.xSleep        = VfsRdOnlySleep;
   vfs.xCurrentTime  = VfsRdOnlyCurrentTime;
   vfs.xGetLastError = VfsRdOnlyGetLastError;
   return vfs;
}

////////////////////////////////////////////////////////////////////////////
/// A global struct of function pointers and details on the VfsRootFile class that together constitue a VFS module
static struct sqlite3_vfs kSqlite3Vfs = GetSqlite3Vfs();

static bool RegisterSqliteVfs()
{
   int retval;
   retval = sqlite3_vfs_register(&kSqlite3Vfs, false);
   return (retval == SQLITE_OK);
}

} // anonymous namespace

namespace ROOT {

namespace RDF {

namespace Internal {
////////////////////////////////////////////////////////////////////////////
/// The state of an open dataset in terms of the sqlite3 C library.
struct RSqliteDSDataSet {
   sqlite3 *fDb = nullptr;
   sqlite3_stmt *fQuery = nullptr;
};
}

RSqliteDS::Value_t::Value_t(RSqliteDS::ETypes type)
   : fType(type), fIsActive(false), fInteger(0), fReal(0.0), fText(), fBlob(), fNull(nullptr)
{
   switch (type) {
   case ETypes::kInteger: fPtr = &fInteger; break;
   case ETypes::kReal: fPtr = &fReal; break;
   case ETypes::kText: fPtr = &fText; break;
   case ETypes::kBlob: fPtr = &fBlob; break;
   case ETypes::kNull: fPtr = &fNull; break;
   default: throw std::runtime_error("Internal error");
   }
}

constexpr char const *RSqliteDS::fgTypeNames[];

////////////////////////////////////////////////////////////////////////////
/// \brief Build the dataframe
/// \param[in] fileName The path to an sqlite3 file, will be opened read-only
/// \param[in] query A valid sqlite3 SELECT query
///
/// The constructor opens the sqlite file, prepares the query engine and determines the column names and types.
RSqliteDS::RSqliteDS(const std::string &fileName, const std::string &query)
   : fDataSet(std::make_unique<Internal::RSqliteDSDataSet>()), fNSlots(0), fNRow(0)
{
   static bool hasSqliteVfs = RegisterSqliteVfs();
   if (!hasSqliteVfs)
      throw std::runtime_error("Cannot register SQlite VFS in RSqliteDS");

   int retval;

   retval = sqlite3_open_v2(fileName.c_str(), &fDataSet->fDb, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX,
                            gSQliteVfsName);
   if (retval != SQLITE_OK)
      SqliteError(retval);

   retval = sqlite3_prepare_v2(fDataSet->fDb, query.c_str(), -1, &fDataSet->fQuery, nullptr);
   if (retval != SQLITE_OK)
      SqliteError(retval);

   int colCount = sqlite3_column_count(fDataSet->fQuery);
   retval = sqlite3_step(fDataSet->fQuery);
   if ((retval != SQLITE_ROW) && (retval != SQLITE_DONE))
      SqliteError(retval);

   fValues.reserve(colCount);
   for (int i = 0; i < colCount; ++i) {
      fColumnNames.emplace_back(sqlite3_column_name(fDataSet->fQuery, i));
      int type = SQLITE_NULL;
      // Try first with the declared column type and then with the dynamic type
      // for expressions
      const char *declTypeCstr = sqlite3_column_decltype(fDataSet->fQuery, i);
      if (declTypeCstr == nullptr) {
         if (retval == SQLITE_ROW)
            type = sqlite3_column_type(fDataSet->fQuery, i);
      } else {
         std::string declType(declTypeCstr);
         std::transform(declType.begin(), declType.end(), declType.begin(), ::toupper);
         if (declType == "INTEGER")
            type = SQLITE_INTEGER;
         else if (declType == "FLOAT")
            type = SQLITE_FLOAT;
         else if (declType == "TEXT")
            type = SQLITE_TEXT;
         else if (declType == "BLOB")
            type = SQLITE_BLOB;
         else
            throw std::runtime_error("Unexpected column decl type");
      }

      switch (type) {
      case SQLITE_INTEGER:
         fColumnTypes.push_back(ETypes::kInteger);
         fValues.emplace_back(ETypes::kInteger);
         break;
      case SQLITE_FLOAT:
         fColumnTypes.push_back(ETypes::kReal);
         fValues.emplace_back(ETypes::kReal);
         break;
      case SQLITE_TEXT:
         fColumnTypes.push_back(ETypes::kText);
         fValues.emplace_back(ETypes::kText);
         break;
      case SQLITE_BLOB:
         fColumnTypes.push_back(ETypes::kBlob);
         fValues.emplace_back(ETypes::kBlob);
         break;
      case SQLITE_NULL:
         // TODO: Null values in first rows are not well handled
         fColumnTypes.push_back(ETypes::kNull);
         fValues.emplace_back(ETypes::kNull);
         break;
      default: throw std::runtime_error("Unhandled data type");
      }
   }
}

////////////////////////////////////////////////////////////////////////////
/// Frees the sqlite resources and closes the file.
RSqliteDS::~RSqliteDS()
{
   // sqlite3_finalize returns the error code of the most recent operation on fQuery.
   sqlite3_finalize(fDataSet->fQuery);
   // Closing can possibly fail with SQLITE_BUSY, in which case resources are leaked. This should not happen
   // the way it is used in this class because we cleanup the prepared statement before.
   sqlite3_close(fDataSet->fDb);
}

////////////////////////////////////////////////////////////////////////////
/// Returns the SELECT queries names. The column names have been cached in the constructor.
/// For expressions, the column name is the string of the expression unless the query defines a column name with as
/// like in "SELECT 1 + 1 as mycolumn FROM table"
const std::vector<std::string> &RSqliteDS::GetColumnNames() const
{
   return fColumnNames;
}

////////////////////////////////////////////////////////////////////////////
/// Activates the given column's result value.
RDataSource::Record_t RSqliteDS::GetColumnReadersImpl(std::string_view name, const std::type_info &ti)
{
   const auto index = std::distance(fColumnNames.begin(), std::find(fColumnNames.begin(), fColumnNames.end(), name));
   const auto type = fColumnTypes[index];

   if ((type == ETypes::kInteger && typeid(Long64_t) != ti) || (type == ETypes::kReal && typeid(double) != ti) ||
       (type == ETypes::kText && typeid(std::string) != ti) ||
       (type == ETypes::kBlob && typeid(std::vector<unsigned char>) != ti) ||
       (type == ETypes::kNull && typeid(void *) != ti)) {
      std::string errmsg = "The type selected for column \"";
      errmsg += name;
      errmsg += "\" does not correspond to column type, which is ";
      errmsg += GetTypeName(name);
      throw std::runtime_error(errmsg);
   }

   fValues[index].fIsActive = true;
   return std::vector<void *>{fNSlots, &fValues[index].fPtr};
}

////////////////////////////////////////////////////////////////////////////
/// Returns a range of size 1 as long as more rows are available in the SQL result set.
/// This inherently serialized the RDF independent of the number of slots.
std::vector<std::pair<ULong64_t, ULong64_t>> RSqliteDS::GetEntryRanges()
{
   std::vector<std::pair<ULong64_t, ULong64_t>> entryRanges;
   int retval = sqlite3_step(fDataSet->fQuery);
   switch (retval) {
   case SQLITE_DONE: return entryRanges;
   case SQLITE_ROW:
      entryRanges.emplace_back(fNRow, fNRow + 1);
      fNRow++;
      return entryRanges;
   default:
      SqliteError(retval);
      // Never here
      abort();
   }
}

////////////////////////////////////////////////////////////////////////////
/// Returns the C++ type for a given column name, implemented as a linear search through all the columns.
std::string RSqliteDS::GetTypeName(std::string_view colName) const
{
   unsigned N = fColumnNames.size();

   for (unsigned i = 0; i < N; ++i) {
      if (colName == fColumnNames[i]) {
         return fgTypeNames[static_cast<int>(fColumnTypes[i])];
      }
   }
   throw std::runtime_error("Unknown column: " + std::string(colName));
}

////////////////////////////////////////////////////////////////////////////
/// A linear search through the columns for the given name
bool RSqliteDS::HasColumn(std::string_view colName) const
{
   return std::find(fColumnNames.begin(), fColumnNames.end(), colName) != fColumnNames.end();
}

////////////////////////////////////////////////////////////////////////////
/// Resets the SQlite query engine at the beginning of the event loop.
void RSqliteDS::Initialise()
{
   fNRow = 0;
   int retval = sqlite3_reset(fDataSet->fQuery);
   if (retval != SQLITE_OK)
      throw std::runtime_error("SQlite error, reset");
}

std::string RSqliteDS::GetLabel()
{
   return "RSqliteDS";
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Factory method to create a SQlite RDataFrame.
/// \param[in] fileName Path of the sqlite file.
/// \param[in] query SQL query that defines the data set.
RDataFrame MakeSqliteDataFrame(std::string_view fileName, std::string_view query)
{
   ROOT::RDataFrame rdf(std::make_unique<RSqliteDS>(std::string(fileName), std::string(query)));
   return rdf;
}

////////////////////////////////////////////////////////////////////////////
/// Stores the result of the current active sqlite query row as a C++ value.
bool RSqliteDS::SetEntry(unsigned int /* slot */, ULong64_t entry)
{
   R__ASSERT(entry + 1 == fNRow);
   unsigned N = fValues.size();
   for (unsigned i = 0; i < N; ++i) {
      if (!fValues[i].fIsActive)
         continue;

      int nbytes;
      switch (fValues[i].fType) {
      case ETypes::kInteger: fValues[i].fInteger = sqlite3_column_int64(fDataSet->fQuery, i); break;
      case ETypes::kReal: fValues[i].fReal = sqlite3_column_double(fDataSet->fQuery, i); break;
      case ETypes::kText:
         nbytes = sqlite3_column_bytes(fDataSet->fQuery, i);
         if (nbytes == 0) {
            fValues[i].fText = "";
         } else {
            fValues[i].fText = reinterpret_cast<const char *>(sqlite3_column_text(fDataSet->fQuery, i));
         }
         break;
      case ETypes::kBlob:
         nbytes = sqlite3_column_bytes(fDataSet->fQuery, i);
         fValues[i].fBlob.resize(nbytes);
         if (nbytes > 0) {
            std::memcpy(fValues[i].fBlob.data(), sqlite3_column_blob(fDataSet->fQuery, i), nbytes);
         }
         break;
      case ETypes::kNull: break;
      default: throw std::runtime_error("Unhandled column type");
      }
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// Almost a no-op, many slots can in fact reduce the performance due to thread synchronization.
void RSqliteDS::SetNSlots(unsigned int nSlots)
{
   if (nSlots > 1) {
      ::Warning("SetNSlots", "Currently the SQlite data source faces performance degradation in multi-threaded mode. "
                             "Consider turning off IMT.");
   }
   fNSlots = nSlots;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// Helper function to throw an exception if there is a fatal sqlite error, e.g. an I/O error.
void RSqliteDS::SqliteError(int errcode)
{
   std::string errmsg = "SQlite error: ";
#if SQLITE_VERSION_NUMBER < 3007015
   errmsg += std::to_string(errcode);
#else
   errmsg += sqlite3_errstr(errcode);
#endif
   throw std::runtime_error(errmsg);
}

} // namespace RDF

} // namespace ROOT
