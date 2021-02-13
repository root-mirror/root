// Author: Enric Tejedor CERN  10/2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// clang-format off
/** \class ROOT::RDF::RCsvDS
    \ingroup dataframe
    \brief RDataFrame data source class for reading CSV files.

The RCsvDS class implements a CSV file reader for RDataFrame.

A RDataFrame that reads from a CSV file can be constructed using the factory method
ROOT::RDF::MakeCsvDataFrame, which accepts three parameters:
1. Path to the CSV file.
2. Boolean that specifies whether the first row of the CSV file contains headers or
not (optional, default `true`). If `false`, header names will be automatically generated as Col0, Col1, ..., ColN.
3. Delimiter (optional, default ',').

The types of the columns in the CSV file are automatically inferred. The supported
types are:
- Integer: stored as a 64-bit long long int.
- Floating point number: stored with double precision.
- Boolean: matches the literals `true` and `false`.
- String: stored as an std::string, matches anything that does not fall into any of the
previous types.

These are some formatting rules expected by the RCsvDS implementation:
- All records must have the same number of fields, in the same order.
- Any field may be quoted.
~~~
    "1997","Ford","E350"
~~~
- Fields with embedded delimiters (e.g. comma) must be quoted.
~~~
    1997,Ford,E350,"Super, luxurious truck"
~~~
- Fields with double-quote characters must be quoted, and each of the embedded
double-quote characters must be represented by a pair of double-quote characters.
~~~
    1997,Ford,E350,"Super, ""luxurious"" truck"
~~~
- Fields with embedded line breaks are not supported, even when quoted.
~~~
    1997,Ford,E350,"Go get one now
    they are going fast"
~~~
- Spaces are considered part of a field and are not ignored.
~~~
    1997, Ford , E350
    not same as
    1997,Ford,E350
    but same as
    1997, "Ford" , E350
~~~
- If a header row is provided, it must contain column names for each of the fields.
~~~
    Year,Make,Model
    1997,Ford,E350
    2000,Mercury,Cougar
~~~

The current implementation of RCsvDS reads the entire CSV file content into memory before
RDataFrame starts processing it. Therefore, before creating a CSV RDataFrame, it is
important to check both how much memory is available and the size of the CSV file.
*/
// clang-format on

#include <ROOT/RDF/Utils.hxx>
#include <ROOT/TSeq.hxx>
#include <ROOT/RCsvDS.hxx>
#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RRawFile.hxx>
#include <TError.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

namespace ROOT {

namespace RDF {

std::string RCsvDS::AsString()
{
   return "CSV data source";
}

// Regular expressions for type inference
const TRegexp RCsvDS::fgIntRegex("^[-+]?[0-9]+$");
const TRegexp RCsvDS::fgDoubleRegex1("^[-+]?[0-9]+\\.[0-9]*$");
const TRegexp RCsvDS::fgDoubleRegex2("^[-+]?[0-9]*\\.[0-9]+$");
const TRegexp RCsvDS::fgDoubleRegex3("^[-+]?[0-9]*\\.[0-9]+[eEdDqQ][-+]?[0-9]+$");
const TRegexp RCsvDS::fgTrueRegex("^true$");
const TRegexp RCsvDS::fgFalseRegex("^false$");

const std::map<RCsvDS::ColType_t, std::string>
   RCsvDS::fgColTypeMap({{'b', "bool"}, {'d', "double"}, {'l', "Long64_t"}, {'s', "std::string"}});

void RCsvDS::FillHeaders(const std::string &line)
{
   auto columns = ParseColumns(line);
   for (auto &col : columns) {
      fHeaders.emplace_back(col);
   }
}

void RCsvDS::FillRecord(const std::string &line, Record_t &record)
{
   std::istringstream lineStream(line);
   auto i = 0U;

   auto columns = ParseColumns(line);

   for (auto &col : columns) {
      auto colType = fColTypes[fHeaders[i]];

      switch (colType) {
      case 'd': {
         record.emplace_back(new double(std::stod(col)));
         break;
      }
      case 'l': {
         record.emplace_back(new Long64_t(std::stoll(col)));
         break;
      }
      case 'b': {
         auto b = new bool();
         record.emplace_back(b);
         std::istringstream is(col);
         is >> std::boolalpha >> *b;
         break;
      }
      case 's': {
         record.emplace_back(new std::string(col));
         break;
      }
      }
      ++i;
   }
}

void RCsvDS::GenerateHeaders(size_t size)
{
   for (size_t i = 0; i < size; ++i) {
      fHeaders.push_back("Col" + std::to_string(i));
   }
}

std::vector<void *> RCsvDS::GetColumnReadersImpl(std::string_view colName, const std::type_info &ti)
{
   const auto colType = GetType(colName);

   if ((colType == 'd' && typeid(double) != ti) || (colType == 'l' && typeid(Long64_t) != ti) ||
       (colType == 's' && typeid(std::string) != ti) || (colType == 'b' && typeid(bool) != ti)) {
      std::string err = "The type selected for column \"";
      err += colName;
      err += "\" does not correspond to column type, which is ";
      err += fgColTypeMap.at(colType);
      throw std::runtime_error(err);
   }

   const auto &colNames = GetColumnNames();
   const auto index = std::distance(colNames.begin(), std::find(colNames.begin(), colNames.end(), colName));
   std::vector<void *> ret(fNSlots);
   for (auto slot : ROOT::TSeqU(fNSlots)) {
      auto &val = fColAddresses[index][slot];
      if (ti == typeid(double)) {
         val = &fDoubleEvtValues[index][slot];
      } else if (ti == typeid(Long64_t)) {
         val = &fLong64EvtValues[index][slot];
      } else if (ti == typeid(std::string)) {
         val = &fStringEvtValues[index][slot];
      } else {
         val = &fBoolEvtValues[index][slot];
      }
      ret[slot] = &val;
   }
   return ret;
}

void RCsvDS::InferColTypes(std::vector<std::string> &columns)
{
   auto i = 0U;
   for (auto &col : columns) {
      InferType(col, i);
      ++i;
   }
}

void RCsvDS::InferType(const std::string &col, unsigned int idxCol)
{
   ColType_t type;
   int dummy;

   if (fgIntRegex.Index(col, &dummy) != -1) {
      type = 'l'; // Long64_t
   } else if (fgDoubleRegex1.Index(col, &dummy) != -1 ||
              fgDoubleRegex2.Index(col, &dummy) != -1 ||
              fgDoubleRegex3.Index(col, &dummy) != -1) {
      type = 'd'; // double
   } else if (fgTrueRegex.Index(col, &dummy) != -1 || fgFalseRegex.Index(col, &dummy) != -1) {
      type = 'b'; // bool
   } else {       // everything else is a string
      type = 's'; // std::string
   }
   // TODO: Date

   fColTypes[fHeaders[idxCol]] = type;
   fColTypesList.push_back(type);
}

std::vector<std::string> RCsvDS::ParseColumns(const std::string &line)
{
   std::vector<std::string> columns;

   for (size_t i = 0; i < line.size(); ++i) {
      i = ParseValue(line, columns, i);
   }

   return columns;
}

size_t RCsvDS::ParseValue(const std::string &line, std::vector<std::string> &columns, size_t i)
{
   std::stringstream val;
   bool quoted = false;

   for (; i < line.size(); ++i) {
      if (line[i] == fDelimiter && !quoted) {
         break;
      } else if (line[i] == '"') {
         // Keep just one quote for escaped quotes, none for the normal quotes
         if (line[i + 1] != '"') {
            quoted = !quoted;
         } else {
            val << line[++i];
         }
      } else {
         val << line[i];
      }
   }

   columns.emplace_back(val.str());

   return i;
}

////////////////////////////////////////////////////////////////////////
/// Constructor to create a CSV RDataSource for RDataFrame.
/// \param[in] fileName Path or URL of the CSV file.
/// \param[in] readHeaders `true` if the CSV file contains headers as first row, `false` otherwise
///                        (default `true`).
/// \param[in] delimiter Delimiter character (default ',').
RCsvDS::RCsvDS(std::string_view fileName, bool readHeaders, char delimiter, Long64_t linesChunkSize) // TODO: Let users specify types?
   : fReadHeaders(readHeaders),
     fCsvFile(ROOT::Internal::RRawFile::Create(fileName)),
     fDelimiter(delimiter),
     fLinesChunkSize(linesChunkSize)
{
   std::string line;

   // Read the headers if present
   if (fReadHeaders) {
      if (fCsvFile->Readln(line)) {
         FillHeaders(line);
      } else {
         std::string msg = "Error reading headers of CSV file ";
         msg += fileName;
         throw std::runtime_error(msg);
      }
   }

   fDataPos = fCsvFile->GetFilePos();
   bool eof = false;
   do {
      eof = !fCsvFile->Readln(line);
   } while (line.empty() && !eof);
   if (!eof) {
      auto columns = ParseColumns(line);

      // Generate headers if not present
      if (!fReadHeaders) {
         GenerateHeaders(columns.size());
      }

      // Infer types of columns with first record
      InferColTypes(columns);

      // rewind
      fCsvFile->Seek(fDataPos);
   } else {
      std::string msg = "Could not infer column types of CSV file ";
      msg += fileName;
      throw std::runtime_error(msg);
   }
}

void RCsvDS::FreeRecords()
{
   for (auto &record : fRecords) {
      for (size_t i = 0; i < record.size(); ++i) {
         void *p = record[i];
         const auto colType = fColTypes[fHeaders[i]];
         switch (colType) {
         case 'd': {
            delete static_cast<double *>(p);
            break;
         }
         case 'l': {
            delete static_cast<Long64_t *>(p);
            break;
         }
         case 'b': {
            delete static_cast<bool *>(p);
            break;
         }
         case 's': {
            delete static_cast<std::string *>(p);
            break;
         }
         }
      }
   }
   fRecords.clear();
}

////////////////////////////////////////////////////////////////////////
/// Destructor.
RCsvDS::~RCsvDS()
{
   FreeRecords();
}

void RCsvDS::Finalise()
{
   fCsvFile->Seek(fDataPos);
   fProcessedLines = 0ULL;
   fEntryRangesRequested = 0ULL;
   FreeRecords();
}

const std::vector<std::string> &RCsvDS::GetColumnNames() const
{
   return fHeaders;
}

std::vector<std::pair<ULong64_t, ULong64_t>> RCsvDS::GetEntryRanges()
{

   // Read records and store them in memory
   auto linesToRead = fLinesChunkSize;
   FreeRecords();

   std::string line;
   while ((-1LL == fLinesChunkSize || 0 != linesToRead) && fCsvFile->Readln(line)) {
      if (line.empty()) continue; // skip empty lines
      fRecords.emplace_back();
      FillRecord(line, fRecords.back());
      --linesToRead;
   }

   if (gDebug > 0) {
      if (fLinesChunkSize == -1LL) {
         Info("GetEntryRanges", "Attempted to read entire CSV file into memory, %zu lines read", fRecords.size());
      } else {
         Info("GetEntryRanges", "Attempted to read chunk of %lld lines of CSV file into memory, %zu lines read", fLinesChunkSize, fRecords.size());
      }
   }

   std::vector<std::pair<ULong64_t, ULong64_t>> entryRanges;
   const auto nRecords = fRecords.size();
   if (0 == nRecords)
      return entryRanges;

   const auto chunkSize = nRecords / fNSlots;
   const auto remainder = 1U == fNSlots ? 0 : nRecords % fNSlots;
   auto start = 0ULL == fEntryRangesRequested ? 0ULL : fProcessedLines;
   auto end = start;

   for (auto i : ROOT::TSeqU(fNSlots)) {
      start = end;
      end += chunkSize;
      entryRanges.emplace_back(start, end);
      (void)i;
   }
   entryRanges.back().second += remainder;

   fProcessedLines += nRecords;
   fEntryRangesRequested++;

   return entryRanges;
}

RCsvDS::ColType_t RCsvDS::GetType(std::string_view colName) const
{
   if (!HasColumn(colName)) {
      std::string msg = "The dataset does not have column ";
      msg += colName;
      throw std::runtime_error(msg);
   }

   return fColTypes.at(colName.data());
}

std::string RCsvDS::GetTypeName(std::string_view colName) const
{
   return fgColTypeMap.at(GetType(colName));
}

bool RCsvDS::HasColumn(std::string_view colName) const
{
   return fHeaders.end() != std::find(fHeaders.begin(), fHeaders.end(), colName);
}

bool RCsvDS::SetEntry(unsigned int slot, ULong64_t entry)
{
   // Here we need to normalise the entry to the number of lines we already processed.
   const auto offset = (fEntryRangesRequested - 1) * fLinesChunkSize;
   const auto recordPos = entry - offset;
   int colIndex = 0;
   for (auto &colType : fColTypesList) {
      auto dataPtr = fRecords[recordPos][colIndex];
      switch (colType) {
      case 'd': {
         fDoubleEvtValues[colIndex][slot] = *static_cast<double *>(dataPtr);
         break;
      }
      case 'l': {
         fLong64EvtValues[colIndex][slot] = *static_cast<Long64_t *>(dataPtr);
         break;
      }
      case 'b': {
         fBoolEvtValues[colIndex][slot] = *static_cast<bool *>(dataPtr);
         break;
      }
      case 's': {
         fStringEvtValues[colIndex][slot] = *static_cast<std::string *>(dataPtr);
         break;
      }
      }
      colIndex++;
   }
   return true;
}

void RCsvDS::SetNSlots(unsigned int nSlots)
{
   R__ASSERT(0U == fNSlots && "Setting the number of slots even if the number of slots is different from zero.");

   fNSlots = nSlots;

   const auto nColumns = fHeaders.size();
   // Initialise the entire set of addresses
   fColAddresses.resize(nColumns, std::vector<void *>(fNSlots, nullptr));

   // Initialize the per event data holders
   fDoubleEvtValues.resize(nColumns, std::vector<double>(fNSlots));
   fLong64EvtValues.resize(nColumns, std::vector<Long64_t>(fNSlots));
   fStringEvtValues.resize(nColumns, std::vector<std::string>(fNSlots));
   fBoolEvtValues.resize(nColumns, std::deque<bool>(fNSlots));
}

std::string RCsvDS::GetLabel()
{
   return "RCsv";
}

RDataFrame MakeCsvDataFrame(std::string_view fileName, bool readHeaders, char delimiter, Long64_t linesChunkSize)
{
   ROOT::RDataFrame tdf(std::make_unique<RCsvDS>(fileName, readHeaders, delimiter, linesChunkSize));
   return tdf;
}

} // ns RDF

} // ns ROOT
